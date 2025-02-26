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
module BPACK_Utilities
   use MISC_Utilities
   use Bplus_Utilities

contains


   recursive subroutine Hmat_GetBlkLst(blocks, ptree, lstblks, Maxlevel)
      implicit none
      ! type(Hoption)::option
      ! type(Hstat)::stats
      ! type(mesh)::msh
      type(proctree)::ptree
      integer nth, bidx, level_c, Maxlevel
      integer ii, jj, idx, row_group, col_group
      type(list)::lstblks(0:Maxlevel)
      type(nod), pointer::cur
      class(*), pointer::ptr
      type(matrixblock), pointer::blocks, blocks_son
      integer flag
      type(block_ptr)::blk_ptr

      row_group = blocks%row_group
      col_group = blocks%col_group
      ! if (IOwnPgrp(ptree, blocks%pgno)) then
         if (blocks%style == 4) then ! divided blocks
            do ii = 1, 2
            do jj = 1, 2
               blocks_son => blocks%sons(ii, jj)
               call Hmat_GetBlkLst(blocks_son, ptree, lstblks, Maxlevel)
            enddo
            enddo
         else
            blk_ptr%ptr => blocks
            call append(lstblks(blocks%level), blk_ptr)
         endif
      ! endif
   end subroutine Hmat_GetBlkLst


   subroutine HODLR_copy(ho_bf_i, ho_bf_o, ptree)


      implicit none

      type(proctree)::ptree
      type(hobf)::ho_bf_i, ho_bf_o

      integer ii
      integer level_c

! real(kind=8),optional::memory
! real(kind=8)::rtemp

      ho_bf_o%Maxlevel = ho_bf_i%Maxlevel
      ho_bf_o%N = ho_bf_i%N

      allocate (ho_bf_o%levels(ho_bf_o%Maxlevel + 1))
      do level_c = 1, ho_bf_o%Maxlevel + 1
         ho_bf_o%levels(level_c)%level = ho_bf_i%levels(level_c)%level
         ho_bf_o%levels(level_c)%N_block_forward = ho_bf_i%levels(level_c)%N_block_forward
         ho_bf_o%levels(level_c)%N_block_inverse = ho_bf_i%levels(level_c)%N_block_inverse
         ho_bf_o%levels(level_c)%Bidxs = ho_bf_i%levels(level_c)%Bidxs
         ho_bf_o%levels(level_c)%Bidxe = ho_bf_i%levels(level_c)%Bidxe

         allocate (ho_bf_o%levels(level_c)%BP(ho_bf_o%levels(level_c)%N_block_forward))
         ! write(*,*)ho_bf_o%levels(level_c)%N_block_inverse,'g'
         allocate (ho_bf_o%levels(level_c)%BP_inverse(ho_bf_o%levels(level_c)%N_block_inverse))
         ! write(*,*)ho_bf_o%levels(level_c)%N_block_inverse,'g1'
         allocate (ho_bf_o%levels(level_c)%BP_inverse_update(ho_bf_o%levels(level_c)%N_block_forward))
         allocate (ho_bf_o%levels(level_c)%BP_inverse_schur(ho_bf_o%levels(level_c)%N_block_inverse))

         do ii = 1, ho_bf_o%levels(level_c)%N_block_forward
            call Bplus_copy(ho_bf_i%levels(level_c)%BP(ii), ho_bf_o%levels(level_c)%BP(ii))
            call Bplus_copy(ho_bf_i%levels(level_c)%BP_inverse_update(ii), ho_bf_o%levels(level_c)%BP_inverse_update(ii))
         end do
         ! if(level_c/=ho_bf_o%Maxlevel+1)then
         do ii = 1, ho_bf_o%levels(level_c)%N_block_inverse
            ! write(*,*)ii,'6642'
            call Bplus_copy(ho_bf_i%levels(level_c)%BP_inverse(ii), ho_bf_o%levels(level_c)%BP_inverse(ii))
            if (level_c < ho_bf_i%Maxlevel + 1)call Bplus_copy(ho_bf_i%levels(level_c)%BP_inverse_schur(ii), ho_bf_o%levels(level_c)%BP_inverse_schur(ii))
         end do
         ! end if
      end do

   end subroutine HODLR_copy


   subroutine Hmat_copy(h_mat_i, h_mat_o, ptree)

      implicit none
      integer ii,jj,bm,bn,level
      type(Hmat)::h_mat_i, h_mat_o
      type(matrixblock), pointer :: blocks1,blocks2
      type(proctree)::ptree

      h_mat_o%Maxlevel=h_mat_i%Maxlevel
      h_mat_o%N=h_mat_i%N
      h_mat_o%Dist_level=h_mat_i%Dist_level
      h_mat_o%idxs=h_mat_i%idxs
      h_mat_o%idxe=h_mat_i%idxe
      h_mat_o%myArows=h_mat_i%myArows
      h_mat_o%myAcols=h_mat_i%myAcols

      if(associated(h_mat_i%N_p))then
         allocate(h_mat_o%N_p(size(h_mat_i%N_p,1),size(h_mat_i%N_p,2)))
         h_mat_o%N_p=h_mat_i%N_p
      endif

      if(allocated(h_mat_i%basis_group))then
         allocate (h_mat_o%basis_group(size(h_mat_i%basis_group,1)))
         h_mat_o%basis_group = h_mat_i%basis_group
      endif

      if(associated(h_mat_i%Local_blocks))then
         bm = size(h_mat_i%Local_blocks, 1)
         bn = size(h_mat_i%Local_blocks, 2)
         allocate (h_mat_o%Local_blocks(bm, bn))
         do ii = 1, bm
         do jj = 1, bn
            blocks1 => h_mat_i%Local_blocks(ii, jj)
            blocks2 => h_mat_o%Local_blocks(ii, jj)
            call Hmat_block_copy('N',blocks2,blocks1)
         enddo
         enddo
      endif

      if(associated(h_mat_i%Local_blocks_copy))then
         bm = size(h_mat_i%Local_blocks_copy, 1)
         bn = size(h_mat_i%Local_blocks_copy, 2)
         allocate (h_mat_o%Local_blocks_copy(bm, bn))
         do ii = 1, bm
         do jj = 1, bn
            blocks1 => h_mat_i%Local_blocks_copy(ii, jj)
            blocks2 => h_mat_o%Local_blocks_copy(ii, jj)
            call Hmat_block_copy('N',blocks2,blocks1)
         enddo
         enddo
      endif

      if(allocated(h_mat_i%colorsets))then
         allocate(h_mat_o%colorsets(0:h_mat_i%Maxlevel))
         do level = 0, h_mat_i%Maxlevel
            allocate(h_mat_o%colorsets(level)%dat(2**level))
            h_mat_o%colorsets(level)%dat=h_mat_i%colorsets(level)%dat
            h_mat_o%colorsets(level)%idx=h_mat_i%colorsets(level)%idx
         enddo
      endif

      if(allocated(h_mat_i%fullmat))then
         allocate(h_mat_o%fullmat(size(h_mat_i%fullmat,1),size(h_mat_i%fullmat,2)))
         h_mat_o%fullmat=h_mat_i%fullmat
      endif

      if (allocated(h_mat_i%lstblks) .and. associated(h_mat_i%Local_blocks)) then
         allocate (h_mat_o%lstblks(0:h_mat_o%Maxlevel))
         do level = 0, h_mat_o%Maxlevel
            h_mat_o%lstblks(level) = list()
         enddo
         bm = size(h_mat_i%Local_blocks, 1)
         bn = size(h_mat_i%Local_blocks, 2)
         do ii = 1, bm
            do jj = 1, bn
               blocks1 => h_mat_o%Local_blocks(ii, jj)
               call Hmat_GetBlkLst(blocks1, ptree, h_mat_o%lstblks, h_mat_o%Maxlevel)
            enddo
         enddo
         do level = 0, h_mat_o%Maxlevel
            call MergeSort(h_mat_o%lstblks(level)%head, node_score_block_ptr_row)
         enddo
      endif
   end subroutine Hmat_copy



   subroutine Hmat_delete(h_mat)

      implicit none

      type(Hmat)::h_mat
      integer bm, bn, ii, jj, level

      if (associated(h_mat%N_p)) deallocate (h_mat%N_p)
      if (allocated(h_mat%basis_group)) deallocate (h_mat%basis_group)
      if (associated(h_mat%Local_blocks)) then
         bm = size(h_mat%Local_blocks, 1)
         bn = size(h_mat%Local_blocks, 2)
         do ii = 1, bm
         do jj = 1, bn
            call Hmat_block_delete(h_mat%Local_blocks(ii, jj))
         enddo
         enddo
         deallocate (h_mat%Local_blocks)
      endif

      if (associated(h_mat%Local_blocks_copy)) then
         bm = size(h_mat%Local_blocks_copy, 1)
         bn = size(h_mat%Local_blocks_copy, 2)
         do ii = 1, bm
         do jj = 1, bn
            call Hmat_block_delete(h_mat%Local_blocks_copy(ii, jj))
         enddo
         enddo
         deallocate (h_mat%Local_blocks_copy)
      endif

      if (allocated(h_mat%lstblks)) then
         do level = 0, h_mat%Maxlevel
            call list_finalizer(h_mat%lstblks(level))
         enddo
         deallocate (h_mat%lstblks)
      endif

      if(allocated(h_mat%colorsets))then
         do level = 0, h_mat%Maxlevel
            deallocate(h_mat%colorsets(level)%dat)
         enddo
         deallocate (h_mat%colorsets)
      endif
      if(allocated(h_mat%fullmat))then
         deallocate(h_mat%fullmat)
      endif
      if(allocated(h_mat%fullmat2D))then
         deallocate(h_mat%fullmat2D)
      endif
   end subroutine Hmat_delete

   subroutine HODLR_delete(ho_bf_o)

      implicit none

      type(hobf)::ho_bf_o

      integer ii
      integer level_c

      do level_c = 1, ho_bf_o%Maxlevel + 1
         do ii = 1, ho_bf_o%levels(level_c)%N_block_forward
            call Bplus_delete(ho_bf_o%levels(level_c)%BP(ii))
            call Bplus_delete(ho_bf_o%levels(level_c)%BP_inverse_update(ii))
         end do
         do ii = 1, ho_bf_o%levels(level_c)%N_block_inverse
            call Bplus_delete(ho_bf_o%levels(level_c)%BP_inverse(ii))
            call Bplus_delete(ho_bf_o%levels(level_c)%BP_inverse_schur(ii))
         end do
         deallocate (ho_bf_o%levels(level_c)%BP)
         deallocate (ho_bf_o%levels(level_c)%BP_inverse_update)
         deallocate (ho_bf_o%levels(level_c)%BP_inverse)
         deallocate (ho_bf_o%levels(level_c)%BP_inverse_schur)
      end do
      deallocate (ho_bf_o%levels)
      if(allocated(ho_bf_o%fullmat2D))then
         deallocate(ho_bf_o%fullmat2D)
      endif
   end subroutine HODLR_delete


   subroutine HSSBF_delete(hss_bf1)

      implicit none

      type(hssbf)::hss_bf1

      call Bplus_delete(hss_bf1%BP)
      call Bplus_delete(hss_bf1%BP_inverse)

      if(allocated(hss_bf1%fullmat2D))then
         deallocate(hss_bf1%fullmat2D)
      endif
   end subroutine HSSBF_delete



   subroutine HSSBF_MD_delete(hss_bf1_md)
      implicit none
      integer ll,bb
      type(hssbf_md)::hss_bf1_md

      if (associated(hss_bf1_md%BP%LL)) then
         do ll = 1, LplusMax
            if (hss_bf1_md%BP%LL(ll)%Nbound > 0) then
               if (associated(hss_bf1_md%BP%LL(ll)%matrices_block)) then
               do bb = 1, hss_bf1_md%BP%LL(ll)%Nbound
                  ! write(*,*)ll,hss_bf1_md%BP%Lplus,bb,hss_bf1_md%BP%LL(ll)%Nbound,'fff'
                  call BF_MD_delete(hss_bf1_md%Ndim, hss_bf1_md%BP%LL(ll)%matrices_block(bb), 1)
               end do
               deallocate (hss_bf1_md%BP%LL(ll)%matrices_block)
               endif
               if (allocated(hss_bf1_md%BP%LL(ll)%boundary_map)) deallocate (hss_bf1_md%BP%LL(ll)%boundary_map)
            end if
         end do
         deallocate (hss_bf1_md%BP%LL)
      endif

      if(allocated(hss_bf1_md%BP%row_group))then
         deallocate(hss_bf1_md%BP%row_group)
      endif
      if(allocated(hss_bf1_md%BP%col_group))then
         deallocate(hss_bf1_md%BP%col_group)
      endif
      if(allocated(hss_bf1_md%N))then
         deallocate(hss_bf1_md%N)
      endif

   end subroutine HSSBF_MD_delete



   subroutine BPACK_delete(bmat)

      implicit none
      type(Bmatrix)::bmat
      if (associated(bmat%ho_bf)) then
         call HODLR_delete(bmat%ho_bf)
         deallocate (bmat%ho_bf)
         bmat%ho_bf => null()
      endif
      if (associated(bmat%h_mat)) then
         call Hmat_delete(bmat%h_mat)
         deallocate (bmat%h_mat)
         bmat%h_mat => null()
      endif
      if (associated(bmat%hss_bf)) then
         call HSSBF_delete(bmat%hss_bf)
         deallocate (bmat%hss_bf)
         bmat%hss_bf => null()
      endif

      if (associated(bmat%hss_bf_md)) then
         call HSSBF_MD_delete(bmat%hss_bf_md)
         deallocate (bmat%hss_bf_md)
         bmat%hss_bf_md => null()
      endif

   end subroutine BPACK_delete


   subroutine BPACK_copy(bmat_i,bmat_o,ptree)

      implicit none
      type(Bmatrix)::bmat_i,bmat_o
      type(proctree)::ptree
      if (associated(bmat_i%ho_bf)) then
         if(.not. associated(bmat_o%ho_bf))allocate (bmat_o%ho_bf)
         call HODLR_copy(bmat_i%ho_bf, bmat_o%ho_bf, ptree)
      endif
      if (associated(bmat_i%h_mat)) then
         if(.not. associated(bmat_o%h_mat))allocate (bmat_o%h_mat)
         call Hmat_copy(bmat_i%h_mat, bmat_o%h_mat, ptree)
      endif
      if (associated(bmat_i%hss_bf)) then
         write(*,*)'HSS-BF copy not yet implemented'
         stop
      endif
   end subroutine BPACK_copy



   subroutine delete_kernelquant(ker)

      implicit none
      type(kernelquant)::ker
      if (allocated(ker%matZ_glo)) deallocate (ker%matZ_glo)
   end subroutine delete_kernelquant

   subroutine delete_mesh(msh,stats)

      implicit none
      type(mesh)::msh
      type(Hstat),optional::stats
      integer ii

      if (allocated(msh%xyz))then
         if(present(stats))call LogMemory(stats, -SIZEOF(msh%xyz)/1024.0d3)
         deallocate (msh%xyz)
      endif
      if (allocated(msh%nns))then
         if(present(stats))call LogMemory(stats, -SIZEOF(msh%nns)/1024.0d3)
         deallocate (msh%nns)
      endif
      if (allocated(msh%new2old))then
         if(present(stats))call LogMemory(stats, -SIZEOF(msh%new2old)/1024.0d3)
         deallocate (msh%new2old)
      endif
      if (allocated(msh%old2new))then
         if(present(stats))call LogMemory(stats, -SIZEOF(msh%old2new)/1024.0d3)
         deallocate (msh%old2new)
      endif
      if (allocated(msh%pretree)) deallocate (msh%pretree)
      if (allocated(msh%basis_group)) then
      do ii = 1, msh%Maxgroup
         if (allocated(msh%basis_group(ii)%center)) deallocate (msh%basis_group(ii)%center)
         if (allocated(msh%basis_group(ii)%nlist)) deallocate (msh%basis_group(ii)%nlist)
         msh%basis_group(ii)%nn = 0
      enddo
      if(present(stats))call LogMemory(stats, -SIZEOF(msh%basis_group)/1024.0d3)
      deallocate (msh%basis_group)
      endif

   end subroutine delete_mesh

   subroutine delete_proctree(ptree)

      implicit none
      type(proctree)::ptree
      integer ii, Maxgrp
      integer ierr

      if (allocated(ptree%pgrp)) then
         Maxgrp = 2**(ptree%nlevel) - 1
         do ii = 1, Maxgrp
            ! if(associated(ptree%pgrp(ii)%gd))then
            ! call delete_grid(ptree%pgrp(ii)%gd)
            ! deallocate(ptree%pgrp(ii)%gd)
            ! ptree%pgrp(ii)%gd=>null()
            ! endif
            if (ptree%pgrp(ii)%ctxt /= -1) call blacs_gridexit_wrp(ptree%pgrp(ii)%ctxt)
            if (ptree%pgrp(ii)%ctxt1D /= -1) call blacs_gridexit_wrp(ptree%pgrp(ii)%ctxt1D)
            if (ptree%pgrp(ii)%ctxt1DCol /= -1) call blacs_gridexit_wrp(ptree%pgrp(ii)%ctxt1DCol)
            if (ptree%pgrp(ii)%ctxt_head /= -1) call blacs_gridexit_wrp(ptree%pgrp(ii)%ctxt_head)
            if (ptree%pgrp(ii)%Comm /= MPI_COMM_NULL) call MPI_Comm_free(ptree%pgrp(ii)%Comm, ierr)
         enddo
         deallocate (ptree%pgrp)
      endif
      if (ptree%Comm /= MPI_COMM_NULL) call MPI_Comm_free(ptree%Comm, ierr)

   end subroutine delete_proctree

! recursive subroutine delete_grid(gd)

! implicit none
! type(grid)::gd
! integer ierr

! if(.not. associated(gd%gdc))then
   ! if(gd%ctxt/=-1)call blacs_gridexit_wrp(gd%ctxt)
   ! if(gd%Comm/=MPI_COMM_NULL)call MPI_Comm_free(gd%Comm,ierr)
   ! return
! else
   ! call delete_grid(gd%gdc(1))
   ! call delete_grid(gd%gdc(2))
   ! deallocate(gd%gdc)
   ! gd%gdc=>null()
! endif
! end subroutine delete_grid

   subroutine delete_Hstat(stats)

      implicit none
      type(Hstat)::stats

      if (allocated(stats%leafs_of_level)) deallocate (stats%leafs_of_level)
      if (allocated(stats%rankmax_of_level)) deallocate (stats%rankmax_of_level)
      if (allocated(stats%rankmin_of_level)) deallocate (stats%rankmin_of_level)
      if (allocated(stats%rankmax_of_level_global)) deallocate (stats%rankmax_of_level_global)
      if (allocated(stats%rankmax_of_level_global_factor)) deallocate (stats%rankmax_of_level_global_factor)
      if (allocated(stats%Add_random_CNT)) deallocate (stats%Add_random_CNT)
      if (allocated(stats%Mul_random_CNT)) deallocate (stats%Mul_random_CNT)
      if (allocated(stats%XLUM_random_CNT)) deallocate (stats%XLUM_random_CNT)
      if (allocated(stats%Add_random_Time)) deallocate (stats%Add_random_Time)
      if (allocated(stats%Mul_random_Time)) deallocate (stats%Mul_random_Time)
      if (allocated(stats%XLUM_random_Time)) deallocate (stats%XLUM_random_Time)

   end subroutine delete_Hstat

   recursive subroutine copy_basis_group(basis_group1, node1, Maxgroup1, basis_group2, node2, Maxgroup2, offset)
      implicit none
      type(basisgroup):: basis_group1(:), basis_group2(:)
      integer node1, node2, Maxgroup1, Maxgroup2, offset
      integer Dimn

      if (node2 <= Maxgroup2 .and. node1 <= Maxgroup1) then

         basis_group2(node2)%head = basis_group1(node1)%head + offset
         basis_group2(node2)%tail = basis_group1(node1)%tail + offset
         ! basis_group2(node2)%pgno = basis_group1(node1)%pgno
         basis_group2(node2)%radius = basis_group1(node1)%radius
         if(allocated(basis_group1(node1)%center))then
            Dimn = size(basis_group1(node1)%center, 1)
            allocate(basis_group2(node2)%center(Dimn))
            basis_group2(node2)%center = basis_group1(node1)%center
         endif

         call copy_basis_group(basis_group1, node1*2, Maxgroup1, basis_group2, node2*2, Maxgroup2, offset)
         call copy_basis_group(basis_group1, node1*2 + 1, Maxgroup1, basis_group2, node2*2 + 1, Maxgroup2, offset)
      endif

   end subroutine copy_basis_group


   subroutine CopyMesh(msh_i,msh_o)
      implicit none
      type(mesh)msh_i,msh_o
      msh_o%Nunk=msh_i%Nunk
      msh_o%Dist_level=msh_i%Dist_level
      msh_o%Maxgroup=msh_i%Maxgroup
      msh_o%idxs=msh_i%idxs
      msh_o%idxe=msh_i%idxe
      if(allocated(msh_i%xyz))then
         allocate(msh_o%xyz(size(msh_i%xyz,1),size(msh_i%xyz,2)))
         msh_o%xyz=msh_i%xyz
      endif
      if(allocated(msh_i%new2old))then
         allocate(msh_o%new2old(size(msh_i%new2old,1)))
         msh_o%new2old=msh_i%new2old
      endif
      if(allocated(msh_i%old2new))then
         allocate(msh_o%old2new(size(msh_i%old2new,1)))
         msh_o%old2new=msh_i%old2new
      endif
      if(allocated(msh_i%nns))then
         allocate(msh_o%nns(size(msh_i%nns,1),size(msh_i%nns,2)))
         msh_o%nns=msh_i%nns
      endif
      if(allocated(msh_i%basis_group))then
         allocate(msh_o%basis_group(msh_o%Maxgroup))
         call copy_basis_group(msh_i%basis_group, 1, msh_i%Maxgroup, msh_o%basis_group, 1, msh_o%Maxgroup, 0)
      endif
   end subroutine CopyMesh


   subroutine CopyStat(stats_i,stats_o)
      implicit none
      type(Hstat)stats_i,stats_o

      stats_o%Time_random = stats_i%Time_random
      stats_o%Time_Sblock = stats_i%Time_Sblock
      stats_o%Time_Sol = stats_i%Time_Sol
      stats_o%Time_C_Mult = stats_i%Time_C_Mult
      stats_o%Time_BLK_MVP = stats_i%Time_BLK_MVP
      stats_o%Time_C_Extract = stats_i%Time_C_Extract
      stats_o%Time_Inv = stats_i%Time_Inv
      stats_o%Time_RedistB = stats_i%Time_RedistB
      stats_o%Time_RedistV = stats_i%Time_RedistV
      stats_o%Time_SMW = stats_i%Time_SMW
      stats_o%Time_PartialUpdate = stats_i%Time_PartialUpdate
      stats_o%Time_Fill = stats_i%Time_Fill
      stats_o%Time_Entry = stats_i%Time_Entry
      stats_o%Time_Entry_Traverse = stats_i%Time_Entry_Traverse
      stats_o%Time_Entry_BF = stats_i%Time_Entry_BF
      stats_o%Time_Entry_Comm = stats_i%Time_Entry_Comm
      stats_o%Mem_peak = stats_i%Mem_peak
      stats_o%Mem_Current = stats_i%Mem_Current
      stats_o%Mem_Sblock = stats_i%Mem_Sblock
      stats_o%Mem_SMW = stats_i%Mem_SMW
      stats_o%Mem_Direct_for = stats_i%Mem_Direct_for
      stats_o%Mem_Direct_inv = stats_i%Mem_Direct_inv
      stats_o%Mem_int_vec = stats_i%Mem_int_vec
      stats_o%Mem_Comp_for = stats_i%Mem_Comp_for
      stats_o%Mem_Fill = stats_i%Mem_Fill
      stats_o%Mem_Factor = stats_i%Mem_Factor
      stats_o%Flop_Fill = stats_i%Flop_Fill
      stats_o%Flop_Factor = stats_i%Flop_Factor
      stats_o%Flop_Sol = stats_i%Flop_Sol
      stats_o%Flop_C_Mult = stats_i%Flop_C_Mult
      stats_o%Flop_C_Extract = stats_i%Flop_C_Extract

      stats_o%Time_Direct_LU = stats_i%Time_Direct_LU
      stats_o%Time_Add_Multiply = stats_i%Time_Add_Multiply
      stats_o%Time_Multiply = stats_i%Time_Multiply
      stats_o%Time_XLUM = stats_i%Time_XLUM
      stats_o%Time_Split = stats_i%Time_Split
      stats_o%Time_Comm = stats_i%Time_Comm
      stats_o%Time_Idle = stats_i%Time_Idle
      stats_o%Time_Factor = stats_i%Time_Factor

      if(allocated(stats_i%rankmax_of_level))then
         allocate(stats_o%rankmax_of_level(0:size(stats_i%rankmax_of_level,1)))
         stats_o%rankmax_of_level=stats_i%rankmax_of_level
      endif
      if(allocated(stats_i%rankmin_of_level))then
         allocate(stats_o%rankmin_of_level(0:size(stats_i%rankmin_of_level,1)))
         stats_o%rankmin_of_level=stats_i%rankmin_of_level
      endif
      if(allocated(stats_i%rankmax_of_level_global))then
         allocate(stats_o%rankmax_of_level_global(0:size(stats_i%rankmax_of_level_global,1)))
         stats_o%rankmax_of_level_global=stats_i%rankmax_of_level_global
      endif
      if(allocated(stats_i%rankmax_of_level_global_factor))then
         allocate(stats_o%rankmax_of_level_global_factor(0:size(stats_i%rankmax_of_level_global_factor,1)))
         stats_o%rankmax_of_level_global_factor=stats_i%rankmax_of_level_global_factor
      endif
      if(allocated(stats_i%Add_random_CNT))then
         allocate(stats_o%Add_random_CNT(0:size(stats_i%Add_random_CNT,1)))
         stats_o%Add_random_CNT=stats_i%Add_random_CNT
      endif
      if(allocated(stats_i%Mul_random_CNT))then
         allocate(stats_o%Mul_random_CNT(0:size(stats_i%Mul_random_CNT,1)))
         stats_o%Mul_random_CNT=stats_i%Mul_random_CNT
      endif
      if(allocated(stats_i%XLUM_random_CNT))then
         allocate(stats_o%XLUM_random_CNT(0:size(stats_i%XLUM_random_CNT,1)))
         stats_o%XLUM_random_CNT=stats_i%XLUM_random_CNT
      endif
      if(allocated(stats_i%Add_random_Time))then
         allocate(stats_o%Add_random_Time(0:size(stats_i%Add_random_Time,1)))
         stats_o%Add_random_Time=stats_i%Add_random_Time
      endif
      if(allocated(stats_i%Mul_random_Time))then
         allocate(stats_o%Mul_random_Time(0:size(stats_i%Mul_random_Time,1)))
         stats_o%Mul_random_Time=stats_i%Mul_random_Time
      endif
      if(allocated(stats_i%XLUM_random_Time))then
         allocate(stats_o%XLUM_random_Time(0:size(stats_i%XLUM_random_Time,1)))
         stats_o%XLUM_random_Time=stats_i%XLUM_random_Time
      endif
      if(allocated(stats_i%leafs_of_level))then
         allocate(stats_o%leafs_of_level(0:size(stats_i%leafs_of_level,1)))
         stats_o%leafs_of_level=stats_i%leafs_of_level
      endif
   end subroutine CopyStat
   subroutine InitStat(stats)
      implicit none
      type(Hstat)::stats

      stats%Time_random = 0  ! Intialization, MVP, Reconstruction
      stats%Time_Sblock = 0
      stats%Time_Sol = 0
      stats%Time_C_Mult = 0
      stats%Time_BLK_MVP = 0
      stats%Time_C_Extract = 0
      stats%Time_Inv = 0
      stats%Time_RedistB = 0
      stats%Time_RedistV = 0
      stats%Time_SMW = 0
      stats%Time_PartialUpdate = 0
      stats%Time_Fill = 0
      stats%Time_Entry = 0
      stats%Time_Entry_Traverse = 0
      stats%Time_Entry_BF = 0
      stats%Time_Entry_Comm = 0
      stats%Mem_peak = 0
      stats%Mem_Current = 0
      stats%Mem_Sblock = 0
      stats%Mem_SMW = 0
      stats%Mem_Direct_for = 0
      stats%Mem_Direct_inv = 0
      stats%Mem_int_vec = 0
      stats%Mem_Comp_for = 0
      stats%Mem_Fill = 0
      stats%Mem_Factor = 0
      stats%Flop_Fill = 0
      stats%Flop_Factor = 0
      stats%Flop_Sol = 0
      stats%Flop_C_Mult = 0
      stats%Flop_C_Extract = 0

      stats%Time_Direct_LU = 0
      stats%Time_Add_Multiply = 0
      stats%Time_Multiply = 0
      stats%Time_XLUM = 0
      stats%Time_Split = 0
      stats%Time_Comm = 0
      stats%Time_Idle = 0
      stats%Time_Factor = 0

      time_tmp = 0
      time_tmp1 = 0
      time_tmp2 = 0
      time_tmp3 = 0
      time_tmp4 = 0
      time_tmp5 = 0
   end subroutine InitStat

   subroutine PrintStat(stats, ptree)
      implicit none
      type(Hstat)::stats
      type(proctree)::ptree
      real(kind=8)::rtemp, rtemp1, rtemp2
      integer ierr

      ! stats%Time_random=0  ! Intialization, MVP, Reconstruction
      ! stats%Time_Sblock=0
      ! stats%Time_Sol=0
      ! stats%Time_Inv=0
      ! stats%Time_RedistB=0
      ! stats%Time_RedistV=0
      ! stats%Time_SMW=0
      ! stats%Time_Fill=0
      ! stats%Mem_peak=0
      ! stats%Mem_Sblock=0
      ! stats%Mem_SMW=0
      ! stats%Mem_Direct_for=0
      ! stats%Mem_Direct_inv=0
      ! stats%Mem_int_vec=0
      ! stats%Mem_Comp_for=0
      ! stats%Flop_Fill=0
      ! stats%Flop_Factor=0
      ! stats%Flop_Sol=0

      call MPI_ALLREDUCE(stats%Time_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2,A8)') 'Constr time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Time_Entry, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2,A8)') 'EntryEval time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Time_Entry_Traverse, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2,A8)') '   traversal time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Time_Entry_BF, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2,A8)') '  BF compute time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Time_Entry_Comm, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2,A8)') '  communicate time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Mem_Comp_for, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(stats%Mem_Direct_for, rtemp1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2,A3)') 'Tot constr mem:', rtemp + rtemp1, 'MB'
      call MPI_ALLREDUCE(stats%Flop_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2)') 'Constr flops:', rtemp
      if (ptree%MyID == Main_ID .and. allocated(stats%rankmax_of_level_global)) write (*, '(A21,I14)') 'Rank before factor:', maxval(stats%rankmax_of_level_global)

      call MPI_ALLREDUCE(stats%Time_Factor, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2,A8)') 'Factor time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Mem_Factor, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2,A3)') 'Tot factor mem:', rtemp, 'MB'
      call MPI_ALLREDUCE(stats%Flop_Factor, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2)') 'Factor flops:', rtemp
      if(allocated(stats%rankmax_of_level_global_factor))then
      if (ptree%MyID == Main_ID .and. allocated(stats%rankmax_of_level_global_factor)) write (*, '(A21,I14)') 'Rank after factor:', maxval(stats%rankmax_of_level_global_factor)

      call MPI_ALLREDUCE(stats%Time_Sol, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2,A8)') 'Solve time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Flop_Sol, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2)') 'Solve flops:', rtemp

      call MPI_ALLREDUCE(stats%Time_BLK_MVP, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2,A8)') 'BLK_mult time:', rtemp, 'Seconds'

      call MPI_ALLREDUCE(stats%Time_C_Mult, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2,A8)') 'C_mult time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Flop_C_Mult, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2)') 'C_mult flops:', rtemp

      call MPI_ALLREDUCE(stats%Time_RedistV, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if(ptree%MyID==Main_ID)write(*,'(A21,Es14.2,A8)') 'RedistV time: ', rtemp, 'Seconds'

      call MPI_ALLREDUCE(stats%Time_C_Extract, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2,A8)') 'C_extract time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Flop_C_Extract, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2)') 'C_extract flops:', rtemp

      call MPI_ALLREDUCE(stats%Mem_peak, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) write (*, '(A21,Es14.2,A3)') 'Peak mem:', rtemp, 'MB'
      if (ptree%MyID == Main_ID) write (*, *) 'time_tmp', time_tmp

   end subroutine PrintStat

   subroutine SetDefaultOptions(option)
      implicit none
      type(Hoption)::option

      option%Nmin_leaf = 200
      option%tol_comp = 1d-4
      option%tol_LS = 1d-12
      option%tol_itersol = 1d-6
      option%n_iter = 1000
      option%tol_rand = option%tol_comp
      option%tol_Rdetect = option%tol_comp*1d-1
      option%level_check = 10000
      option%precon = DIRECT
      option%xyzsort = TM
      option%lnoBP = 40000
      option%TwoLayerOnly = 1
      option%touch_para = 0d0
      option%schulzorder = 3
      option%schulzhardstart = 0
      option%schulzsplitlevel = 1
      option%schulzlevel = 3000
      option%LRlevel = 0
      option%ErrFillFull = 0
      option%BACA_Batch = 64
      option%RecLR_leaf = BACA
      option%nogeo = 0
      option%per_geo = 0
      option%hextralevel = 0
      option%periods = 0
      option%ErrSol = 0
      option%LR_BLK_NUM = 1
      option%rank0 = 32
      option%rankrate = 1.2d0
      option%itermax = 10
      option%powiter = 0
      option%ILU = 0
      option%Nbundle = 1
      option%near_para = BPACK_SafeEps
      option%knn_near_para = 20d0
      option%format = HODLR
      option%verbosity = 0
      option%scale_factor = 1d0
      option%jitter = 1d-13
      option%rmax = 3000
      option%forwardN15flag = 0
      option%fastsample_tensor = 0
      option%sample_para = 2.0d0
      option%sample_para_outer = 2.0d0
      option%use_zfp = 0
      ! option%sample_heuristic = 1
      option%pat_comp = 3
      option%elem_extract = 0
      option%knn = 0
      option%cpp = 0
      option%bp_cnt_lr = 1
      option%less_adapt = 0

   end subroutine SetDefaultOptions

   subroutine ReadOption(option, ptree, ii)
      implicit none
      type(Hoption)::option
      type(proctree)::ptree
      integer ii
      integer nargs, flag
      character(len=1024)  :: strings, strings1

      nargs = iargc()
      flag = 1
      do while (flag == 1)
         ii = ii + 1
         if (ii <= nargs) then
            call getarg(ii, strings)
            if (strings(1:2) == '--') then
               ii = ii + 1
               call getarg(ii, strings1)

               if (trim(strings) == '--nmin_leaf') then
                  read (strings1, *) option%Nmin_leaf
               else if (trim(strings) == '--tol_comp') then
                  read (strings1, *) option%tol_comp
                  option%tol_rand = option%tol_comp
                  option%tol_Rdetect = option%tol_comp*1d-1
               else if (trim(strings) == '--tol_rand') then
                  read (strings1, *) option%tol_rand
                  option%tol_Rdetect = option%tol_rand*1d-1
               else if (trim(strings) == '--tol_Rdetect') then
                  read (strings1, *) option%tol_Rdetect
               else if (trim(strings) == '--tol_itersol') then
                  read (strings1, *) option%tol_itersol
               else if (trim(strings) == '--n_iter') then
                  read (strings1, *) option%n_iter
               else if (trim(strings) == '--level_check') then
                  read (strings1, *) option%level_check
               else if (trim(strings) == '--precon') then
                  read (strings1, *) option%precon
               else if (trim(strings) == '--xyzsort') then
                  read (strings1, *) option%xyzsort
               else if (trim(strings) == '--schulzorder') then
                  read (strings1, *) option%schulzorder
               else if (trim(strings) == '--schulzhardstart') then
                  read (strings1, *) option%schulzhardstart
               else if (trim(strings) == '--schulzsplitlevel') then
                  read (strings1, *) option%schulzsplitlevel
               else if (trim(strings) == '--schulzlevel') then
                  read (strings1, *) option%schulzlevel
               else if (trim(strings) == '--lrlevel') then
                  read (strings1, *) option%LRlevel
               else if (trim(strings) == '--errfillfull') then
                  read (strings1, *) option%ErrFillFull
               else if (trim(strings) == '--forwardN15flag') then
                  read (strings1, *) option%forwardN15flag
               else if (trim(strings) == '--fastsample_tensor') then
                  read (strings1, *) option%fastsample_tensor
               else if (trim(strings) == '--baca_batch') then
                  read (strings1, *) option%BACA_Batch
               else if (trim(strings) == '--reclr_leaf') then
                  read (strings1, *) option%RecLR_leaf
               else if (trim(strings) == '--nogeo') then
                  read (strings1, *) option%nogeo
               else if (trim(strings) == '--per_geo') then
                  read (strings1, *) option%per_geo
               else if (trim(strings) == '--hextralevel') then
                  read (strings1, *) option%hextralevel
               else if (trim(strings) == '--period1') then
                  read (strings1, *) option%periods(1)
               else if (trim(strings) == '--period2') then
                  read (strings1, *) option%periods(2)
               else if (trim(strings) == '--period3') then
                  read (strings1, *) option%periods(3)
               else if (trim(strings) == '--less_adapt') then
                  read (strings1, *) option%less_adapt
               else if (trim(strings) == '--errsol') then
                  read (strings1, *) option%ErrSol
               else if (trim(strings) == '--lr_blk_num') then
                  read (strings1, *) option%LR_BLK_NUM
               else if (trim(strings) == '--rank0') then
                  read (strings1, *) option%rank0
               else if (trim(strings) == '--rankrate') then
                  read (strings1, *) option%rankrate
               else if (trim(strings) == '--itermax') then
                  read (strings1, *) option%itermax
               else if (trim(strings) == '--powiter') then
                  read (strings1, *) option%powiter
               else if (trim(strings) == '--ilu') then
                  read (strings1, *) option%ILU
               else if (trim(strings) == '--nbundle') then
                  read (strings1, *) option%Nbundle
               else if (trim(strings) == '--near_para') then
                  read (strings1, *) option%near_para
               else if (trim(strings) == '--knn_near_para') then
                  read (strings1, *) option%knn_near_para
               else if (trim(strings) == '--format') then
                  read (strings1, *) option%format
               else if (trim(strings) == '--verbosity') then
                  read (strings1, *) option%verbosity
               else if (trim(strings) == '--rmax') then
                  read (strings1, *) option%rmax
               else if (trim(strings) == '--sample_para') then
                  read (strings1, *) option%sample_para
               else if (trim(strings) == '--sample_para_outer') then
                  read (strings1, *) option%sample_para_outer
               else if (trim(strings) == '--use_zfp') then
                  read (strings1, *) option%use_zfp
               ! else if (trim(strings) == '--sample_heuristic') then
               !    read (strings1, *) option%sample_heuristic
               else if (trim(strings) == '--pat_comp') then
                  read (strings1, *) option%pat_comp
               else if (trim(strings) == '--elem_extract') then
                  read (strings1, *) option%elem_extract
               else if (trim(strings) == '--knn') then
                  read (strings1, *) option%knn
               else if (trim(strings) == '--cpp') then
                  read (strings1, *) option%cpp
               else if (trim(strings) == '--lnobp') then
                  read (strings1, *) option%lnoBP
               else if (trim(strings) == '--bp_cnt_lr') then
                  read (strings1, *) option%bp_cnt_lr
               else if (trim(strings) == '--touch_para') then
                  read (strings1, *) option%touch_para
               else
                  if (ptree%MyID == Main_ID) write (*, *) 'ignoring unknown option: ', trim(strings)
               endif
            else
               flag = 0
            endif
         else
            flag = 0
         endif
      enddo

   end subroutine ReadOption

   subroutine CopyOptions(option, option1)
      implicit none
      type(Hoption)::option, option1

      option1%Nmin_leaf = option%Nmin_leaf
      option1%tol_comp = option%tol_comp
      option1%tol_Rdetect = option%tol_Rdetect
      option1%tol_LS = option%tol_LS
      option1%tol_itersol = option%tol_itersol
      option1%n_iter = option%n_iter
      option1%tol_rand = option%tol_rand
      option1%level_check = option%level_check
      option1%precon = option%precon
      option1%xyzsort = option%xyzsort
      option1%lnoBP = option%lnoBP
      option1%bp_cnt_lr = option%bp_cnt_lr
      option1%TwoLayerOnly = option%TwoLayerOnly
      option1%touch_para = option%touch_para
      option1%schulzorder = option%schulzorder
      option1%schulzhardstart = option%schulzhardstart
      option1%schulzsplitlevel = option%schulzsplitlevel
      option1%schulzlevel = option%schulzlevel
      option1%LRlevel = option%LRlevel
      option1%ErrFillFull = option%ErrFillFull
      option1%BACA_Batch = option%BACA_Batch
      option1%RecLR_leaf = option%RecLR_leaf
      option1%nogeo = option%nogeo
      option1%per_geo = option%per_geo
      option1%hextralevel = option%hextralevel
      option1%periods = option%periods
      option1%ErrSol = option%ErrSol
      option1%LR_BLK_NUM = option%LR_BLK_NUM
      option1%rank0 = option%rank0
      option1%rankrate = option%rankrate
      option1%itermax = option%itermax
      option1%powiter = option%powiter
      option1%ILU = option%ILU
      option1%Nbundle = option%Nbundle
      option1%near_para = option%near_para
      option1%knn_near_para = option%knn_near_para
      option1%format = option%format
      option1%verbosity = option%verbosity
      option1%scale_factor = option%scale_factor
      option1%jitter = option%jitter
      option1%rmax = option%rmax
      option1%forwardN15flag = option%forwardN15flag
option1%fastsample_tensor = option%fastsample_tensor
      option1%sample_para = option%sample_para
      option1%sample_para_outer = option%sample_para_outer
      option1%use_zfp = option%use_zfp
      ! option1%sample_heuristic = option%sample_heuristic
      option1%pat_comp = option%pat_comp
      option1%elem_extract = option%elem_extract
      option1%cpp = option%cpp
      option1%knn = option%knn
      option1%less_adapt = option%less_adapt

   end subroutine CopyOptions

   subroutine PrintOptions(option, ptree)
      implicit none
      type(Hoption)::option
      type(proctree)::ptree

      if (ptree%MyID == Main_ID) then
         write (*, *) ' '
         write (*, *) '***************************'
         write (*, '(A25)') 'Printing Solver Options:'
         write (*, '(A18,I8)') 'format', option%format
         write (*, '(A18,I8)') 'lrlevel', option%LRlevel
         write (*, '(A18,Es14.7)') 'near_para', option%near_para
         if(option%format==HMAT)then
            if(option%LRlevel==0)then
               write (*, '(A18,A10)') 'algorithm', 'H-LR'
            else
               write (*, '(A18,A10)') 'algorithm', 'H-BF'
            endif
         else if(option%format==HODLR)then
            if(option%LRlevel==0)then
               write (*, '(A18,A10)') 'algorithm', 'HODLR'
            else
               write (*, '(A18,A10)') 'algorithm', 'HODBF'
            endif
         else if(option%format==BLR)then
            if(option%LRlevel==0 .or. option%Hextralevel==0)then
               write (*, '(A18,A10)') 'algorithm', 'Block-LR'
            else
               write (*, '(A18,A10)') 'algorithm', 'Block-BF'
            endif
         else if(option%format==HSS)then
            if(option%LRlevel==0)then
               write (*, '(A18,A10)') 'algorithm', 'INVALID'
               stop
            else
               if(option%near_para>1.0)then
                  write (*, '(A18,A10)') 'algorithm', 'SHNBF'
               else
                  write (*, '(A18,A10)') 'algorithm', 'HSSBF'
               endif
            endif
         else if(option%format==HSS_MD)then
            if(option%LRlevel==0)then
               write (*, '(A18,A10)') 'algorithm', 'INVALID'
               stop
            else
               if(option%near_para>1.0)then
                  write (*, '(A18,A10)') 'algorithm', 'SHNBF_MD'
               else
                  write (*, '(A18,A10)') 'algorithm', 'HSSBF_MD'
               endif
            endif
         endif

         write (*, '(A18,I8)') 'nmin_leaf', option%Nmin_leaf
         write (*, '(A18,I8)') 'n_iter', option%n_iter
         write (*, '(A18,I8)') 'level_check', option%level_check
         write (*, '(A18,I8)') 'precon', option%precon
         write (*, '(A18,I8)') 'xyzsort', option%xyzsort
         write (*, '(A18,I8)') 'lnoBP', option%lnoBP
         write (*, '(A18,I8)') 'bp_cnt_lr', option%bp_cnt_lr
         write (*, '(A18,I8)') 'TwoLayerOnly', option%TwoLayerOnly
         write (*, '(A18,I8)') 'schulzorder', option%schulzorder
         write (*, '(A18,I8)') 'schulzhardstart', option%schulzhardstart
         write (*, '(A18,I8)') 'schulzsplitlevel', option%schulzsplitlevel
         write (*, '(A18,I8)') 'schulzlevel', option%schulzlevel
         write (*, '(A18,I8)') 'baca_batch', option%BACA_Batch
         write (*, '(A18,I8)') 'reclr_leaf', option%RecLR_leaf
         write (*, '(A18,I8)') 'nogeo', option%nogeo
         write (*, '(A18,I8)') 'per_geo', option%per_geo
         write (*, '(A18,I8)') 'hextralevel', option%hextralevel
         write (*, '(A18,Es14.7)') 'period1', option%periods(1)
         write (*, '(A18,Es14.7)') 'period2', option%periods(2)
         write (*, '(A18,Es14.7)') 'period3', option%periods(3)
         write (*, '(A18,I8)') 'lr_blk_num', option%LR_BLK_NUM
         write (*, '(A18,I8)') 'rank0', option%rank0
         write (*, '(A18,I8)') 'itermax', option%itermax
         write (*, '(A18,I8)') 'powiter', option%powiter
         write (*, '(A18,I8)') 'ilu', option%ILU
         write (*, '(A18,I8)') 'nbundle', option%Nbundle
         write (*, '(A18,I8)') 'verbosity', option%verbosity
         write (*, '(A18,I8)') 'rmax', option%rmax
         write (*, '(A18,I8)') 'forwardN15flag', option%forwardN15flag
         write (*, '(A18,I8)') 'fastsample_tensor', option%fastsample_tensor
         write (*, '(A18,I8)') 'pat_comp', option%pat_comp
         write (*, '(A18,I8)') 'elem_extract', option%elem_extract
         write (*, '(A18,I8)') 'cpp', option%cpp
         write (*, '(A18,I8)') 'knn', option%knn
         write (*, '(A18,I8)') 'errfillfull', option%ErrFillFull
         write (*, '(A18,I8)') 'errsol', option%ErrSol
         ! write (*, '(A18,I8)') 'sample_heuristic', option%sample_heuristic
         write (*, '(A18,I8)') 'less_adapt', option%less_adapt
         write (*, '(A18,I8)') 'use_zfp', option%use_zfp

         write (*, '(A18,Es14.7)') 'rankrate', option%rankrate
         write (*, '(A18,Es14.7)') 'tol_comp', option%tol_comp
         write (*, '(A18,Es14.7)') 'tol_Rdetect', option%tol_Rdetect
         write (*, '(A18,Es14.7)') 'tol_LS', option%tol_LS
         write (*, '(A18,Es14.7)') 'tol_itersol', option%tol_itersol
         write (*, '(A18,Es14.7)') 'tol_rand', option%tol_rand
         write (*, '(A18,Es14.7)') 'touch_para', option%touch_para
         write (*, '(A18,Es14.7)') 'knn_near_para', option%knn_near_para
         write (*, '(A18,Es14.7)') 'scale_factor', option%scale_factor
         write (*, '(A18,Es14.7)') 'jitter_factor', option%jitter
         write (*, '(A18,Es14.7)') 'sample_para', option%sample_para
         write (*, '(A18,Es14.7)') 'sample_para_outer', option%sample_para_outer
         write (*, *) '***************************'
         write (*, *) ' '
      endif
   end subroutine PrintOptions

   subroutine BPACK_GetVersionNumber(v_major, v_minor, v_bugfix)
      implicit none
      integer v_major, v_minor, v_bugfix
      v_major = BPACK_MAJOR_VERSION
      v_minor = BPACK_MINOR_VERSION
      v_bugfix = BPACK_PATCH_VERSION

   end subroutine BPACK_GetVersionNumber

end module BPACK_Utilities
