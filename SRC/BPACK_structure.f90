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
module BPACK_structure
   use BPACK_Utilities
   use BPACK_DEFS
contains

! far is returning 1, near if returning 0
   integer function near_or_far_user(group_m, group_n, msh, option, ker, para)
      implicit none
      type(mesh)::msh
      type(Hoption)::option
      type(kernelquant)::ker
      integer group_m, group_n
      real*8 para
      procedure(F_Compressibility), POINTER :: proc1
      procedure(C_Compressibility), POINTER :: proc1_c

      if (option%nogeo == 0 .or. option%nogeo == 4) then  ! geometrical information is provided
         near_or_far_user = near_or_far_geo(group_m, group_n, msh, option, ker, para)
      else if (option%nogeo == 1) then ! no geometrical information is provided, ! only return 0 for self elements
         if (group_m == group_n) then
            near_or_far_user = 0
         else
            near_or_far_user = 1
         endif
      else if (option%nogeo == 2) then  ! no geometrical information is provided, but user provides a compressibility function
         if (option%cpp == 1) then
            call c_f_procpointer(ker%C_FuncNearFar, proc1_C)
            call proc1_C(group_m, group_n, near_or_far_user, ker%C_QuantApp)
         else
            proc1 => ker%FuncNearFar
            call proc1(group_m, group_n, near_or_far_user, ker%QuantApp)
         endif
      endif

   end function near_or_far_user

! far is returning 1, near if returning 0
   integer function near_or_far_geo(group_m, group_n, msh, option, ker, para)
      implicit none

      type(mesh)::msh
      type(Hoption)::option
      type(kernelquant)::ker
      integer group_m, group_n, farblock, level
      integer i, j, ii, jj
      real*8 dis, rad1, rad2, para
      real*8, allocatable:: a(:), b(:)
      integer Dimn

      Dimn = size(msh%basis_group(group_m)%center, 1)
      allocate (a(Dimn))
      allocate (b(Dimn))
      do i = 1, Dimn
         a(i) = msh%basis_group(group_m)%center(i)
         b(i) = msh%basis_group(group_n)%center(i)
      enddo

      rad1 = msh%basis_group(group_m)%radius
      rad2 = msh%basis_group(group_n)%radius

      dis = 0d0
      do i = 1, Dimn
         dis = dis + (a(i) - b(i))**2
      enddo
      dis = sqrt(dis)

      ! write(*,*)dis/((rad1+rad2)/2)

      ! if (dis>para*max(rad1,rad2)) then
      if (dis > para*(rad1 + rad2)/2) then
         near_or_far_geo = 1
      else
         near_or_far_geo = 0
      endif
      deallocate (a, b)

   end function near_or_far_geo

   real(kind=8) function distance_user(edgem, edgen, ker, msh, option, ptree, stats)

      use BPACK_DEFS
      implicit none

      integer edgem, edgen
      real(kind=8) dis
      integer i, j
      integer Dimn
      type(mesh)::msh
      type(kernelquant)::ker
      type(Hoption)::option
      type(proctree)::ptree
      type(Hstat)::stats
      procedure(F_Dist), POINTER :: proc1
      procedure(C_Dist), POINTER :: proc1_c

      if (option%nogeo == 0 .or. option%nogeo == 4) then  ! geometrical information is provided
         distance_user = distance_geo(edgem, edgen, ker, msh, option, ptree, stats)
      else if (option%nogeo == 1) then ! no geometrical information is provided
         if (option%xyzsort == TM_GRAM) then ! try gram distance
            distance_user = distance_gram(edgem, edgen, ker, msh, option, ptree, stats)
         else                  ! only return 0 for self elements
            if (edgem == edgen) then
               distance_user = 0d0
            else
               distance_user = 1d0
            endif
         endif
      else if (option%nogeo == 2) then  ! no geometrical information is provided, but user provides a distance function
         if (option%cpp == 1) then
            call c_f_procpointer(ker%C_FuncDistmn, proc1_C)
            call proc1_C(edgem, edgen, distance_user, ker%C_QuantApp)
         else
            proc1 => ker%FuncDistmn
            call proc1(edgem, edgen, distance_user, ker%QuantApp)
         endif
      endif

      return

   end function distance_user

   real(kind=8) function distance_geo(edgem, edgen, ker, msh, option, ptree, stats)

   use BPACK_DEFS
   implicit none

   integer edgem, edgen
   real(kind=8) dis
   integer i, j
   integer Dimn
   type(mesh)::msh
   type(kernelquant)::ker
   type(Hoption)::option
   type(proctree)::ptree
   type(Hstat)::stats

   call assert(allocated(msh%xyz), 'xyz is not allocated in distance_geo')
   Dimn = size(msh%xyz, 1)

   dis = 0d0
   do i = 1, Dimn
      dis = dis + (msh%xyz(i, edgem) - msh%xyz(i, edgen))**2
   enddo

   distance_geo = dis

   return

end function distance_geo

!**** l2 gram distance^2 between element edgem and edgen is
!     defined as: Re{Z_ii+Z_jj-Z_ij-Z_ji} for SPD, HPD, general symmetric real, and hermitian matrices
!     undefined otherwise
!**** angular gram distance^2 is
!     defined as 1-Z_ij^2/(Z_iiZ_jj)
!     undefined otherwise
!     Use with caution !!!
   real(kind=8) function distance_gram(edgem, edgen, ker, msh, option, ptree, stats)

      use BPACK_DEFS
      implicit none

      integer edgem, edgen, passflag, rows(1), cols(1)
      type(mesh)::msh
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(kernelquant)::ker
      DT r1(1, 1), r2(1, 1), r3(1, 1), r4(1, 1)

      rows(1) = edgem
      cols(1) = edgen

! l2 distance
#if 0
      call element_Zmn_block_user(1, 1, rows, rows, r1, msh, option, ker, 0, passflag, ptree, stats)
      call element_Zmn_block_user(1, 1, cols, cols, r2, msh, option, ker, 0, passflag, ptree, stats)
      call element_Zmn_block_user(1, 1, rows, cols, r3, msh, option, ker, 0, passflag, ptree, stats)
      call element_Zmn_block_user(1, 1, cols, rows, r4, msh, option, ker, 0, passflag, ptree, stats)
      distance_gram = dble(r1(1, 1) + r2(1, 1) - r3(1, 1) - r4(1, 1))
! angular distance
#else
      call element_Zmn_block_user(1, 1, rows, rows, r1, msh, option, ker, 0, passflag, ptree, stats)
      call element_Zmn_block_user(1, 1, cols, cols, r2, msh, option, ker, 0, passflag, ptree, stats)
      call element_Zmn_block_user(1, 1, rows, cols, r3, msh, option, ker, 0, passflag, ptree, stats)
      distance_gram = dble(1d0 - r3(1, 1)**2d0/(r1(1, 1)*r2(1, 1)))
#endif
      return

   end function distance_gram

!**** l2 gram distance^2 between element edgem and edgen is
!     defined as: Re{Z_ii+Z_jj-Z_ij-Z_ji} for SPD, HPD, general symmetric real, and hermitian matrices
!     undefined otherwise
!**** angular gram distance^2 is
!     defined as 1-Z_ij^2/(Z_iiZ_jj)
!     undefined otherwise
!     Use with caution !!!
   subroutine distance_gram_block(nrow, ncol, rows, cols, dists, ker, msh, option, ptree, stats)

      use BPACK_DEFS
      implicit none

      integer nrow, ncol, ii, jj
      integer edgem, edgen, passflag, rows(nrow), cols(ncol), onerow(1), onecol(1)
      type(mesh)::msh
      type(Hoption)::option
      type(Hstat)::stats
      type(proctree)::ptree
      type(kernelquant)::ker
      DT r1(nrow, 1), r2(1, ncol)
      DT r3(nrow, ncol), r4(ncol, nrow)
      real*8:: dists(nrow, ncol)

      do ii = 1, nrow
         onerow(1) = rows(ii)
         call element_Zmn_block_user(1, 1, onerow, onerow, r1(ii, 1), msh, option, ker, 0, passflag, ptree, stats)
      enddo
      do ii = 1, ncol
         onecol(1) = cols(ii)
         call element_Zmn_block_user(1, 1, onecol, onecol, r2(1, ii), msh, option, ker, 0, passflag, ptree, stats)
      enddo
      call element_Zmn_block_user(nrow, ncol, rows, cols, r3, msh, option, ker, 0, passflag, ptree, stats)
      call element_Zmn_block_user(ncol, nrow, cols, rows, r4, msh, option, ker, 0, passflag, ptree, stats)

! l2 distance
#if 0
      do ii = 1, nrow
      do jj = 1, ncol
         dists(ii, jj) = dble(r1(ii, 1) + r2(1, jj) - r3(ii, jj) - r4(jj, ii))
      enddo
      enddo
! angular distance
#else
      do ii = 1, nrow
      do jj = 1, ncol
         dists(ii, jj) = dble(1d0 - r3(ii, jj)**2d0/(r1(ii, 1)*r2(1, jj)))
      enddo
      enddo
#endif
      return

   end subroutine distance_gram_block


   recursive subroutine Hmat_GetBlkLst(blocks, option, stats, msh, ptree, h_mat)
      use BPACK_DEFS
      implicit none
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(proctree)::ptree
      integer nth, bidx, level_c
      integer ii, jj, idx, row_group, col_group
      type(Hmat)::h_mat
      type(nod), pointer::cur
      class(*), pointer::ptr
      type(matrixblock), pointer::blocks, blocks_son
      integer flag
      type(block_ptr)::blk_ptr

      row_group = blocks%row_group
      col_group = blocks%col_group
      if (IOwnPgrp(ptree, blocks%pgno)) then
         if (blocks%style == 4) then ! divided blocks
            do ii = 1, 2
            do jj = 1, 2
               blocks_son => blocks%sons(ii, jj)
               call Hmat_GetBlkLst(blocks_son, option, stats, msh, ptree, h_mat)
            enddo
            enddo
         else
            blk_ptr%ptr => blocks
            call append(h_mat%lstblks(blocks%level), blk_ptr)
         endif
      endif
   end subroutine Hmat_GetBlkLst

   recursive subroutine Hmat_construct_local_tree(blocks, option, stats, msh, ker, ptree, Maxlevel)

      implicit none

      type(matrixblock):: blocks
      type(matrixblock), pointer :: blocks_son
      integer group_m, group_n, i, j, k, level
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer Maxlevel

      group_m = blocks%row_group
      group_n = blocks%col_group
      level = blocks%level

      if (level >= msh%Dist_level .and. near_or_far_user(group_m, group_n, msh, option, ker, option%near_para) == 1) then
         stats%leafs_of_level(level) = stats%leafs_of_level(level) + 1
         blocks%style = 2
         ! blocks%prestyle=2
         ! blocks%data_type=1
      else
         if (level == Maxlevel) then
            blocks%style = 1
            ! blocks%prestyle=1
            ! blocks%data_type=1
            stats%leafs_of_level(level) = stats%leafs_of_level(level) + 1
         else
            blocks%style = 4
            ! blocks%prestyle=4
            ! blocks%data_type=1
            allocate (blocks%sons(2, 2))
            call LogMemory(stats, SIZEOF(blocks%sons)/1024.0d3)
            do j = 1, 2
               do i = 1, 2
                  blocks%sons(i, j)%level = blocks%level + 1
                  ! blocks%sons(i,j)%father=>blocks
                  blocks%sons(i, j)%row_group = 2*blocks%row_group + i - 1
                  blocks%sons(i, j)%col_group = 2*blocks%col_group + j - 1
                  ! blocks%sons(i,j)%data_type=1

                  if (blocks%sons(i, j)%level > option%LRlevel) then
                     blocks%sons(i, j)%level_butterfly = 0 ! low rank below LRlevel
                  else
                     blocks%sons(i, j)%level_butterfly = Maxlevel - blocks%sons(i, j)%level   ! butterfly
                  endif

                  group_m = blocks%sons(i, j)%row_group
                  group_n = blocks%sons(i, j)%col_group

                  blocks%sons(i, j)%pgno = blocks%pgno
                  blocks%sons(i, j)%M = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
                  blocks%sons(i, j)%N = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
                  blocks%sons(i, j)%headm = msh%basis_group(group_m)%head
                  blocks%sons(i, j)%headn = msh%basis_group(group_n)%head
                  call ComputeParallelIndices(blocks%sons(i, j), blocks%sons(i, j)%pgno, ptree, msh)

               enddo
            enddo
            blocks_son => blocks%sons(1, 1)
            call Hmat_construct_local_tree(blocks_son, option, stats, msh, ker, ptree, Maxlevel)
            blocks_son => blocks%sons(2, 1)
            call Hmat_construct_local_tree(blocks_son, option, stats, msh, ker, ptree, Maxlevel)
            blocks_son => blocks%sons(1, 2)
            call Hmat_construct_local_tree(blocks_son, option, stats, msh, ker, ptree, Maxlevel)
            blocks_son => blocks%sons(2, 2)
            call Hmat_construct_local_tree(blocks_son, option, stats, msh, ker, ptree, Maxlevel)

         endif
      endif

      return

   end subroutine Hmat_construct_local_tree

   recursive subroutine Hmat_construct_global_tree(blocks, treelevel, stats)
      implicit none
      integer treelevel
      type(global_matricesblock), pointer :: blocks, blocks_son
      integer group_m, group_n, i, j, k, level
      type(Hstat)::stats

      group_m = blocks%row_group
      group_n = blocks%col_group
      level = blocks%level

      if (level < treelevel) then
         allocate (blocks%sons(2, 2))
         call LogMemory(stats, SIZEOF(blocks%sons)/1024.0d3)
         do j = 1, 2
            do i = 1, 2
               blocks%sons(i, j)%level = blocks%level + 1
               blocks%sons(i, j)%father => blocks
               blocks%sons(i, j)%row_group = 2*blocks%row_group + i - 1
               blocks%sons(i, j)%col_group = 2*blocks%col_group + j - 1
            enddo
         enddo
         blocks_son => blocks%sons(1, 1)
         call Hmat_construct_global_tree(blocks_son, treelevel, stats)
         blocks_son => blocks%sons(2, 1)
         call Hmat_construct_global_tree(blocks_son, treelevel, stats)
         blocks_son => blocks%sons(1, 2)
         call Hmat_construct_global_tree(blocks_son, treelevel, stats)
         blocks_son => blocks%sons(2, 2)
         call Hmat_construct_global_tree(blocks_son, treelevel, stats)
      endif

      return

   end subroutine Hmat_construct_global_tree

   subroutine Cluster_partition(bmat, option, msh, ker, stats, ptree)

      use BPACK_DEFS
      use MISC_Utilities
      implicit none

      integer Cflag
      integer i, j, ii, jj, iii, jjj, kk
      integer level, edge, node, patch, group, group_m, group_n, col_group, row_group, fidx
      integer blocks
      integer center_edge

      integer index_temp
      DT r1, r2, r3, r4, rr(2, 2)
      integer rows(2), cols(2)
      real(kind=8) a, b, c, d, para, xmax, xmin, ymax, ymin, zmax, zmin, seperator, r, theta, phi, phi_tmp
      real(kind=8) radius, radiusmax, radius2, radiusmax2
      real(kind=8), allocatable:: xyzrange(:), xyzmin(:), xyzmax(:), auxpoint(:), groupcenter(:)
      real(kind=8), allocatable :: distance(:), array(:, :), dist_gram(:, :)
      integer level_c, sortdirec, mm, phi_end, Ninfo_edge, ind_i, ind_j
      real(kind=8) t1, t2
      integer Maxgroup, nlevel_pre, passflag
      character(len=1024)  :: strings
      integer, allocatable :: order(:), edge_temp(:), map_temp(:)
      integer dimn, groupsize, idxstart, Nsmp
      type(Hoption)::option
      type(Hstat)::stats
      type(Bmatrix)::bmat
      integer Maxlevel, Maxgrp
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer, allocatable:: perms(:), rows_gram(:), cols_gram(:)
      integer Navr, Bidxs, Bidxe, ierr

      !*************Initialize permutation vector ********
      allocate (msh%new2old(msh%Nunk))
      call LogMemory(stats, SIZEOF(msh%new2old)/1024.0d3)
      do ii = 1, msh%Nunk
         msh%new2old(ii) = ii
      end do

      !************Compute Maxlevel of hodlr tree*******************
      nlevel_pre = 0
      if (allocated(msh%pretree)) then
         nlevel_pre = ceiling_safe(log(dble(size(msh%pretree, 1)))/log(2d0))
      endif
      level = 0; i = 1
      do while (int(msh%Nunk/i) > option%Nmin_leaf)
         level = level + 1
         i = 2**level
      enddo
      Maxlevel = level
      if (Maxlevel < nlevel_pre) Maxlevel = nlevel_pre

      !!!!! make refinement to make sure Maxlevel>=ptree%nlevel, i.e., one processor handles at least two leaf boxes. This refinement can have performance penalties

      ! if(Maxlevel<ptree%nlevel-1)then
      ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'too many processes for paralleling leaf boxes, keep refining the tree ...'
      ! Maxlevel = ptree%nlevel-1
      ! endif

      if (Maxlevel < ptree%nlevel) then
         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'warning: too many processes for paralleling leaf boxes, keep refining the tree ...'
         Maxlevel = ptree%nlevel
      endif

      ! the following is needed when bplus is used as bplus only support even number of levels for now.
      if (Maxlevel == ptree%nlevel .and. option%lnoBP < Maxlevel) then
         Maxlevel = ptree%nlevel + 1
      endif

      select case (option%format)
      case (HODLR)
         allocate (bmat%ho_bf)
         bmat%ho_bf%Maxlevel = Maxlevel
      case (HMAT)
         allocate (bmat%h_mat)
         bmat%h_mat%Maxlevel = Maxlevel
      case (HSS)
         allocate (bmat%hss_bf)
         bmat%hss_bf%Maxlevel = Maxlevel
      end select

      !************** check whether the sorting option is valid
      if (Maxlevel > nlevel_pre) then
         if (.not. allocated(msh%xyz)) then
            if (option%xyzsort == CKD .or. option%xyzsort == TM) then
               write (*, *) 'Geometrical information is not provided. Try use NATRUAL or TM_GRAM ordering'
               stop
            endif
         endif

         if (option%xyzsort == TM_GRAM) then
            call random_number(a)
            ind_i = floor_safe(a*(msh%Nunk - 1)) + 1
            rows(1) = ind_i
            cols(1) = ind_i
            ind_j = ind_i
            do while (ind_i == ind_j)
               call random_number(a)
               ind_j = floor_safe(a*(msh%Nunk - 1)) + 1
            enddo
            rows(2) = ind_j
            cols(2) = ind_j
            call element_Zmn_block_user(2, 2, rows, cols, rr, msh, option, ker, 0, passflag, ptree, stats)
            if (abs(aimag(cmplx(rr(1, 1), kind=8))) > BPACK_SafeUnderflow .or. abs(aimag(cmplx(rr(2, 2), kind=8))) > BPACK_SafeUnderflow .or. abs(rr(1, 2) - conjg(cmplx(rr(2, 1), kind=8))) > abs(rr(1, 2))*BPACK_SafeEps) then
               write (*, *) 'Matrix not hermitian. The gram distance is undefined'
            endif
         endif
      endif

      !***************************************************

      Maxgroup = 2**(Maxlevel + 1) - 1
      msh%Maxgroup = Maxgroup
      allocate (msh%basis_group(Maxgroup))
      call LogMemory(stats, SIZEOF(msh%basis_group)/1024.0d3)

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Maxlevel_for_blocks:', Maxlevel
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'N_leaf:', int(msh%Nunk/(2**Maxlevel))
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Constructing basis groups...'

      dimn = 0
      if (allocated(msh%xyz)) Dimn = size(msh%xyz, 1)
      if (dimn > 0) then
         allocate (xyzrange(dimn))
         allocate (xyzmin(dimn))
         allocate (xyzmax(dimn))
         allocate (auxpoint(dimn))
         allocate (groupcenter(dimn))
      endif

      !**** construct the top few levels whose ordering is provided by the user
      msh%basis_group(1)%head = 1; msh%basis_group(1)%tail = msh%Nunk; msh%basis_group(1)%pgno = 1
      do level = nlevel_pre, 0, -1
         idxstart = 1
         do group = 2**level, 2**(level + 1) - 1
            ! msh%basis_group(group)%level=level

            if (level == nlevel_pre) then
               if (nlevel_pre == 0) then
                  groupsize = msh%Nunk
               else
                  groupsize = msh%pretree(group - 2**nlevel_pre + 1)
               endif
               call assert(groupsize > 0, 'zero leafsize may not be handled')
               msh%basis_group(group)%head = idxstart
               msh%basis_group(group)%tail = idxstart + groupsize - 1
               idxstart = idxstart + groupsize
            else
               msh%basis_group(group)%head = msh%basis_group(2*group)%head
               msh%basis_group(group)%tail = msh%basis_group(2*group + 1)%tail
            endif

            !***** the following is needed for the near_or_far function in H matrix, this needs to be improved
            if (allocated(msh%xyz)) then
               Dimn = size(msh%xyz, 1)
               groupcenter(1:dimn) = 0.0d0
               do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                  do ii = 1, dimn
                     groupcenter(ii) = groupcenter(ii) + msh%xyz(ii, msh%new2old(edge))
                  enddo
               enddo
               do ii = 1, dimn
                  groupcenter(ii) = groupcenter(ii)/(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1)
               enddo

               radiusmax = 0.
               do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                  radius = 0
                  do ii = 1, dimn
                     radius = radius + (msh%xyz(ii, msh%new2old(edge)) - groupcenter(ii))**2
                  enddo
                  radius = sqrt(radius)
                  if (radius > radiusmax) then
                     radiusmax = radius
                     center_edge = edge
                  endif
               enddo

               allocate (msh%basis_group(group)%center(dimn))
               msh%basis_group(group)%center = groupcenter
               msh%basis_group(group)%radius = radiusmax
            endif

         enddo
      enddo

      if (ptree%MyID == Main_ID) then

         !**** if necessary, continue ordering the sub-trees using clustering method specified by option%xyzsort
         do level = nlevel_pre, Maxlevel
            do group = 2**level, 2**(level + 1) - 1
               ! msh%basis_group(group)%level=level

               if (allocated(msh%xyz)) then
               if (.not. allocated(msh%basis_group(group)%center)) then
                  groupcenter(1:dimn) = 0.0d0
                  ! !$omp parallel do default(shared) private(edge,ii) reduction(+:groupcenter)
                  do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                     do ii = 1, dimn
                        groupcenter(ii) = groupcenter(ii) + msh%xyz(ii, msh%new2old(edge))
                     enddo
                  enddo
                  ! !$omp end parallel do
                  do ii = 1, dimn
                     groupcenter(ii) = groupcenter(ii)/(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1)
                  enddo

                  radiusmax = 0.
                  do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                     radius = 0
                     do ii = 1, dimn
                        radius = radius + (msh%xyz(ii, msh%new2old(edge)) - groupcenter(ii))**2
                     enddo
                     radius = sqrt(radius)
                     if (radius > radiusmax) then
                        radiusmax = radius
                        center_edge = edge
                     endif
                  enddo

                  allocate (msh%basis_group(group)%center(dimn))
                  msh%basis_group(group)%center = groupcenter
                  msh%basis_group(group)%radius = radiusmax
               endif
               endif

               if (option%xyzsort == NATURAL) then !natural ordering
                  mm = msh%basis_group(group)%tail - msh%basis_group(group)%head + 1
                  allocate (distance(mm))
                  do i = msh%basis_group(group)%head, msh%basis_group(group)%tail
                     distance(i - msh%basis_group(group)%head + 1) = dble(i)
                  enddo

               else if (option%xyzsort == CKD) then !msh%xyz sort
                  xyzmin = 1d300
                  xyzmax = -1d300
                  do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                     do ii = 1, Dimn
                        xyzmax(ii) = max(xyzmax(ii), msh%xyz(ii, msh%new2old(edge)))
                        xyzmin(ii) = min(xyzmin(ii), msh%xyz(ii, msh%new2old(edge)))
                     enddo
                  enddo
                  xyzrange(1:Dimn) = xyzmax(1:Dimn) - xyzmin(1:Dimn)

                  mm = msh%basis_group(group)%tail - msh%basis_group(group)%head + 1
                  allocate (distance(mm))
                  sortdirec = maxloc(xyzrange(1:Dimn), 1)
                  ! write(*,*)'gaw',sortdirec,xyzrange(1:Dimn)

                  ! ! if(ker%Kernel==EMSURF)then
                  ! if(mod(level,2)==1)then           !!!!!!!!!!!!!!!!!!!!!!!!! note: applys only to plates
                  ! sortdirec=1
                  ! else
                  ! sortdirec=2
                  ! end if
                  ! ! endif

                  !$omp parallel do default(shared) private(i)
                  do i = msh%basis_group(group)%head, msh%basis_group(group)%tail
                     distance(i - msh%basis_group(group)%head + 1) = msh%xyz(sortdirec, msh%new2old(i))
                  enddo
                  !$omp end parallel do

               else if (option%xyzsort == TM) then !cobblestone sort

                  mm = msh%basis_group(group)%tail - msh%basis_group(group)%head + 1
                  allocate (distance(mm))

                  distance(1:mm) = BPACK_Bigvalue
                  !$omp parallel do default(shared) private(i)
                  do i = msh%basis_group(group)%head, msh%basis_group(group)%tail
                     distance(i - msh%basis_group(group)%head + 1) = distance_user(msh%new2old(i), msh%new2old(center_edge), ker, msh, option, ptree, stats)
                  enddo
                  !$omp end parallel do

               else if (option%xyzsort == TM_GRAM) then !GRAM-distance-based cobblestone sort

                  Nsmp = min(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1, 50)
                  allocate (perms(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1))
                  call rperm(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1, perms)

                  allocate (dist_gram(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1, Nsmp))
                  allocate (rows_gram(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1))
                  allocate (cols_gram(Nsmp))
                  do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                     rows_gram(edge - msh%basis_group(group)%head + 1) = msh%new2old(edge)
                  enddo
                  do ii = 1, Nsmp
                     cols_gram(ii) = msh%new2old(perms(ii) + msh%basis_group(group)%head - 1)
                  enddo
                  call distance_gram_block(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1, Nsmp, rows_gram, cols_gram, dist_gram, ker, msh, option, ptree, stats)

                  radiusmax2 = 0.
                  do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                     radius2 = 0
                     do ii = 1, Nsmp  ! take average of distance^2 to Nsmp samples as the distance^2 to the group center
                        radius2 = radius2 + dist_gram(edge - msh%basis_group(group)%head + 1, ii)
                     enddo
                     ! call assert(radius2>0,'radius2<0 cannot take square root')
                     ! radius2 = sqrt(radius2)
                     radius2 = radius2/Nsmp
                     if (radius2 > radiusmax2) then
                        radiusmax2 = radius2
                        center_edge = edge
                     endif
                  enddo
                  deallocate (dist_gram)
                  deallocate (rows_gram)
                  deallocate (cols_gram)

                  mm = msh%basis_group(group)%tail - msh%basis_group(group)%head + 1
                  allocate (distance(mm))

                  distance(1:mm) = BPACK_Bigvalue

                  allocate (dist_gram(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1, 1))
                  allocate (rows_gram(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1))
                  allocate (cols_gram(1))
                  do i = msh%basis_group(group)%head, msh%basis_group(group)%tail
                     rows_gram(i - msh%basis_group(group)%head + 1) = msh%new2old(i)
                  enddo
                  cols_gram = msh%new2old(center_edge)
                  call distance_gram_block(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1, 1, rows_gram, cols_gram, dist_gram, ker, msh, option, ptree, stats)
                  distance = dist_gram(:, 1)
                  deallocate (dist_gram)
                  deallocate (rows_gram)
                  deallocate (cols_gram)

                  deallocate (perms)

               end if

               allocate (order(mm))
               allocate (map_temp(mm))

               call quick_sort(distance, order, mm)
               !$omp parallel do default(shared) private(ii)
               do ii = 1, mm
                  map_temp(ii) = msh%new2old(order(ii) + msh%basis_group(group)%head - 1)
               enddo
               !$omp end parallel do

               !$omp parallel do default(shared) private(ii)
               do ii = 1, mm
                  msh%new2old(ii + msh%basis_group(group)%head - 1) = map_temp(ii)
               enddo
               !$omp end parallel do
               deallocate (map_temp)
               deallocate (order)

               deallocate (distance)

               if (level < Maxlevel) then

                  call assert(msh%basis_group(group)%tail /= msh%basis_group(group)%head, 'detected zero-sized group, try larger leafsizes or smaller MPI counts')
                  msh%basis_group(2*group)%head = msh%basis_group(group)%head
                  msh%basis_group(2*group)%tail = int((msh%basis_group(group)%head + msh%basis_group(group)%tail)/2)
                  msh%basis_group(2*group + 1)%head = msh%basis_group(2*group)%tail + 1
                  msh%basis_group(2*group + 1)%tail = msh%basis_group(group)%tail
               endif
            enddo
         enddo
      endif

      call MPI_Bcast(msh%new2old, msh%Nunk, MPI_integer, Main_ID, ptree%Comm, ierr)

      !**** generate tree structures on other processes
      do level = nlevel_pre, Maxlevel
         do group = 2**level, 2**(level + 1) - 1
            ! msh%basis_group(group)%level=level

            if (allocated(msh%xyz)) then
            if (.not. allocated(msh%basis_group(group)%center)) then
               groupcenter(1:dimn) = 0.0d0
               ! !$omp parallel do default(shared) private(edge,ii) reduction(+:groupcenter)
               do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                  do ii = 1, dimn
                     groupcenter(ii) = groupcenter(ii) + msh%xyz(ii, msh%new2old(edge))
                  enddo
               enddo
               ! !$omp end parallel do
               do ii = 1, dimn
                  groupcenter(ii) = groupcenter(ii)/(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1)
               enddo

               radiusmax = 0.
               do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                  radius = 0
                  do ii = 1, dimn
                     radius = radius + (msh%xyz(ii, msh%new2old(edge)) - groupcenter(ii))**2
                  enddo
                  radius = sqrt(radius)
                  if (radius > radiusmax) then
                     radiusmax = radius
                  endif
               enddo

               allocate (msh%basis_group(group)%center(dimn))
               msh%basis_group(group)%center = groupcenter
               msh%basis_group(group)%radius = radiusmax
            endif
            endif

            if (option%xyzsort == CKD) then !msh%xyz sort
               xyzmin = 1d300
               xyzmax = -1d300
               do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                  do ii = 1, Dimn
                     xyzmax(ii) = max(xyzmax(ii), msh%xyz(ii, msh%new2old(edge)))
                     xyzmin(ii) = min(xyzmin(ii), msh%xyz(ii, msh%new2old(edge)))
                  enddo
               enddo
               xyzrange(1:Dimn) = xyzmax(1:Dimn) - xyzmin(1:Dimn)
               sortdirec = maxloc(xyzrange(1:Dimn), 1)
               seperator = msh%xyz(sortdirec, msh%new2old(int((msh%basis_group(group)%head + msh%basis_group(group)%tail)/2)))
               msh%basis_group(group)%boundary(1) = sortdirec
               msh%basis_group(group)%boundary(2) = seperator
            end if

            if (level < Maxlevel) then
               call assert(msh%basis_group(group)%tail /= msh%basis_group(group)%head, 'detected zero-sized group, try larger leafsizes or smaller MPI counts')
               msh%basis_group(2*group)%head = msh%basis_group(group)%head
               msh%basis_group(2*group)%tail = int((msh%basis_group(group)%head + msh%basis_group(group)%tail)/2)
               msh%basis_group(2*group + 1)%head = msh%basis_group(2*group)%tail + 1
               msh%basis_group(2*group + 1)%tail = msh%basis_group(group)%tail
            endif
         enddo
      enddo

      if (dimn > 0) then
         deallocate (xyzrange)
         deallocate (xyzmin)
         deallocate (xyzmax)
         deallocate (auxpoint)
         deallocate (groupcenter)
      endif

      allocate (msh%old2new(msh%Nunk))
      do ii = 1, msh%Nunk
         msh%old2new(msh%new2old(ii)) = ii
      end do

      ! do ii=1,msh%Nunk
      ! write(110,*)msh%old2new(ii)
      ! enddo

      !**********Dump the ordering into a file********************************

#if        0
      write (strings, *) Dimn
      do level = 0, Maxlevel
         do group = 2**level, 2**(level + 1) - 1
            do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
               write (113, '(I5,I8,'//TRIM(strings)//'Es16.8)') level, group, msh%xyz(1:Dimn, msh%new2old(edge))
            enddo
         enddo
      enddo
#endif

      if (option%nogeo == 1 .and. option%knn > 0) then
         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "no geometrical information or distance function provided, force option%knn to be 0"
         option%knn = 0
      endif

!**** construct a list of k-nearest neighbours for each point
      if (option%knn > 0 .and. option%nogeo /= 3 .and. option%nogeo /= 4) then
         call FindKNNs(option, msh, ker, stats, ptree, 1, 1)
      endif

      return

   end subroutine Cluster_partition

   subroutine FindKNNs(option, msh, ker, stats, ptree, groupm_start, groupn_start)
      implicit none
      real(kind=8), allocatable :: distance(:, :)
      integer, allocatable :: order(:, :), edge_temp(:, :)

      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer ii, iii, jjj, kk, jj, Bidxs, Bidxe, Navr, Maxgroup, Maxlevel, ierr
      integer num_threads
      integer, save:: my_tid = 0
      integer groupm_start, groupn_start
      real(kind=8) t1, t2, tmp

      real(kind=8), allocatable :: dist_gram(:, :)
      integer, allocatable:: rows_gram(:), cols_gram(:)

!$omp threadprivate(my_tid)

!$omp parallel default(shared)
!$omp master
      num_threads = omp_get_num_threads()
!$omp end master
      my_tid = omp_get_thread_num()
!$omp end parallel


      tmp=0

      t1 = OMP_get_wtime()
      allocate (msh%nns(msh%Nunk, option%knn))
      call LogMemory(stats, SIZEOF(msh%nns)/1024.0d3)
      msh%nns = 0
      call MPI_barrier(ptree%Comm, ierr)
      ! if(ptree%MyID==Main_ID)write(*,*)'nn0',tmp,'nns:',  SIZEOF(msh%nns)/1024.0d3, msh%Nunk*option%knn*4/1024.0d3


      allocate (distance(msh%Nunk, num_threads))
      allocate (order(msh%Nunk, num_threads))
      allocate (edge_temp(msh%Nunk, num_threads))
      distance = BPACK_Bigvalue
      ! call MPI_barrier(ptree%Comm, ierr)
      ! if(ptree%MyID==Main_ID)write(*,*)'nn0',tmp,'nns:',  SIZEOF(msh%nns)/1024.0d3, msh%Nunk*option%knn*4/1024.0d3

      Maxgroup = size(msh%basis_group, 1)
      Maxlevel = GetTreelevel(Maxgroup) - 1
      Navr = 2**Maxlevel/ptree%nproc
      Bidxs = 2**Maxlevel + ptree%MyID*Navr
      Bidxe = 2**Maxlevel + (ptree%MyID + 1)*Navr - 1
      if (ptree%MyID == ptree%nproc - 1) Bidxe = 2**(Maxlevel + 1) - 1


      ! call MPI_barrier(ptree%Comm, ierr)
      ! if(ptree%MyID==Main_ID)write(*,*)ptree%MyID,'nn1.0',tmp,'nns:',  SIZEOF(msh%nns)/1024.0d3,Bidxe-Bidxs+1,stats%Mem_Peak



      do ii = Bidxs, Bidxe
         msh%basis_group(ii)%nn = 0
      enddo
      call append_nlist(ker, option, stats, msh, ptree, groupm_start, groupn_start, 0,Bidxs, Bidxe)


      do ii = Bidxs, Bidxe
         if (msh%basis_group(ii)%nn > 0) then
            allocate (msh%basis_group(ii)%nlist(msh%basis_group(ii)%nn))
            tmp=tmp+SIZEOF(msh%basis_group(ii)%nlist)/1024.0d3
            msh%basis_group(ii)%nn = 0
         endif
      enddo
      call append_nlist(ker, option, stats, msh, ptree, groupm_start, groupn_start, 1,Bidxs, Bidxe)


      ! if(ptree%MyID==Main_ID)write(*,*)'nn2',tmp,'nns:',  SIZEOF(msh%nns)/1024.0d3,stats%Mem_Peak
      ! call MPI_barrier(ptree%Comm, ierr)

      if (option%xyzsort == 3) then
         do ii = Bidxs, Bidxe

            kk = 0
            do jj = 1, msh%basis_group(ii)%nn
               kk = kk + msh%basis_group(msh%basis_group(ii)%nlist(jj))%tail - msh%basis_group(msh%basis_group(ii)%nlist(jj))%head + 1
            enddo
            allocate (dist_gram(msh%basis_group(ii)%tail - msh%basis_group(ii)%head + 1, kk))
            allocate (rows_gram(msh%basis_group(ii)%tail - msh%basis_group(ii)%head + 1))
            allocate (cols_gram(kk))

            kk = 0
            do jj = 1, msh%basis_group(ii)%nn
               do jjj = msh%basis_group(msh%basis_group(ii)%nlist(jj))%head, msh%basis_group(msh%basis_group(ii)%nlist(jj))%tail
                  kk = kk + 1
                  cols_gram(kk) = msh%new2old(jjj)
               enddo
            enddo

            do iii = msh%basis_group(ii)%head, msh%basis_group(ii)%tail
               rows_gram(iii - msh%basis_group(ii)%head + 1) = msh%new2old(iii)
            enddo

            call distance_gram_block(msh%basis_group(ii)%tail - msh%basis_group(ii)%head + 1, kk, rows_gram, cols_gram, dist_gram, ker, msh, option, ptree, stats)

            !$omp parallel do default(shared) private(iii,kk,jj,jjj)
            do iii = msh%basis_group(ii)%head, msh%basis_group(ii)%tail
               kk = 0
               do jj = 1, msh%basis_group(ii)%nn
                  do jjj = msh%basis_group(msh%basis_group(ii)%nlist(jj))%head, msh%basis_group(msh%basis_group(ii)%nlist(jj))%tail
                     kk = kk + 1
                     distance(kk, my_tid + 1) = dist_gram(iii - msh%basis_group(ii)%head + 1, kk)
                     edge_temp(kk, my_tid + 1) = jjj
                  enddo
               enddo
               call quick_sort(distance(:, my_tid + 1), order(:, my_tid + 1), kk)
               kk = min(kk, option%knn)
               msh%nns(iii, 1:kk) = edge_temp(order(1:kk, my_tid + 1), my_tid + 1)
            enddo
            !$omp end parallel do

            deallocate (dist_gram)
            deallocate (rows_gram)
            deallocate (cols_gram)

         enddo
      else
         !$omp parallel do default(shared) private(ii,iii,kk,jj,jjj)
         do ii = Bidxs, Bidxe
            do iii = msh%basis_group(ii)%head, msh%basis_group(ii)%tail
               kk = 0
               do jj = 1, msh%basis_group(ii)%nn
                  do jjj = msh%basis_group(msh%basis_group(ii)%nlist(jj))%head, msh%basis_group(msh%basis_group(ii)%nlist(jj))%tail
                     if (iii /= jjj) then
                        kk = kk + 1
                        distance(kk, my_tid + 1) = distance_user(msh%new2old(iii), msh%new2old(jjj), ker, msh, option, ptree, stats)
                        edge_temp(kk, my_tid + 1) = jjj
                     endif
                  enddo
               enddo
               call quick_sort(distance(:, my_tid + 1), order(:, my_tid + 1), kk)
               kk = min(kk, option%knn)
               msh%nns(iii, 1:kk) = edge_temp(order(1:kk, my_tid + 1), my_tid + 1)
            enddo
         enddo
         !$omp end parallel do
      endif

      do ii=1,option%knn
         call MPI_ALLREDUCE(MPI_IN_PLACE, msh%nns(:,ii), msh%Nunk, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      enddo

      do ii = Bidxs, Bidxe
         if (msh%basis_group(ii)%nn > 0) deallocate(msh%basis_group(ii)%nlist)
      enddo

      deallocate (distance)
      deallocate (order)
      deallocate (edge_temp)

      t2 = OMP_get_wtime()
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Finding neighbours time: ", t2 - t1

   end subroutine FindKNNs

   recursive subroutine append_nlist(ker, option, stats, msh, ptree, group_m, group_n, flag,Bidxs, Bidxe)
      use BPACK_DEFS
      implicit none
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(proctree)::ptree
      type(kernelquant)::ker
      integer flag, group_m, group_n
      integer ii, jj,Bidxs, Bidxe


      if (group_m*2 + 1 > size(msh%basis_group)) then
         if (group_m == group_n) then
            if(group_m>=Bidxs .and. group_m<=Bidxe)then
            msh%basis_group(group_m)%nn = msh%basis_group(group_m)%nn + 1
            if (flag == 1) then
               msh%basis_group(group_m)%nlist(msh%basis_group(group_m)%nn) = group_n
            endif
            endif
         else if (group_m < group_n) then ! only search in the upper block triangular matrix
            if(group_m>=Bidxs .and. group_m<=Bidxe)msh%basis_group(group_m)%nn = msh%basis_group(group_m)%nn + 1
            if(group_n>=Bidxs .and. group_n<=Bidxe)msh%basis_group(group_n)%nn = msh%basis_group(group_n)%nn + 1
            if (flag == 1) then
               if(group_m>=Bidxs .and. group_m<=Bidxe)msh%basis_group(group_m)%nlist(msh%basis_group(group_m)%nn) = group_n
               if(group_n>=Bidxs .and. group_n<=Bidxe)msh%basis_group(group_n)%nlist(msh%basis_group(group_n)%nn) = group_m
            endif
         endif
      else
         do ii = 1, 2
         do jj = 1, 2
            if (near_or_far_user(group_m*2 + ii - 1, group_n*2 + jj - 1, msh, option, ker, option%knn_near_para) == 0) call append_nlist(ker, option, stats, msh, ptree, group_m*2 + ii - 1, group_n*2 + jj - 1, flag,Bidxs, Bidxe)
         enddo
         enddo
      endif
   end subroutine append_nlist

   subroutine BPACK_structuring(bmat, option, msh, ker, ptree, stats)
      implicit none
      type(Hoption)::option
      type(mesh)::msh
      type(kernelquant)::ker
      type(Hstat)::stats
      type(Bmatrix)::bmat
      type(proctree)::ptree

      select case (option%format)
      case (HODLR)
         call HODLR_structuring(bmat%ho_bf, option, msh, ker, ptree, stats)
      case (HMAT)
         call Hmat_structuring(bmat%h_mat, option, msh, ker, ptree, stats)
      case (HSS)
         call HSS_structuring(bmat%hss_bf, option, msh, ker, ptree, stats)
      end select
   end subroutine BPACK_structuring

   subroutine HSS_structuring(hss_bf1, option, msh, ker, ptree, stats)
      use BPACK_DEFS
      use MISC_Utilities
      implicit none

      integer i, j, ii, jj, kk, iii, jjj, ll, bb, sortdirec, ii_sch, pgno_bplus
      integer level, edge, patch, node, group, group_touch
      integer rank, index_near, m, n, length, flag, itemp, cnt, detection
      real T0
      real(kind=8):: tolerance, rtemp, rel_error, seperator, dist
      real(kind=8) Memory_direct_forward, Memory_butterfly_forward
      integer mm, nn, header_m, header_n, edge_m, edge_n, group_m, group_n, group_m1, group_n1, group_m2, group_n2, levelm, groupm_start, index_i_m, index_j_m
      integer level_c, iter, level_cc, level_BP, Nboundall, level_butterfly
      type(matrixblock), pointer::blocks, block_f, block_sch, block_inv
      real(kind=8)::minbound, theta, phi, r, rmax, phi_tmp, measure
      real(kind=8), allocatable::Centroid_M(:, :), Centroid_N(:, :)
      integer, allocatable::Isboundary_M(:), Isboundary_N(:)
      integer Dimn, col_group, row_group, Maxgrp
      type(Hoption)::option
      type(mesh)::msh
      type(kernelquant)::ker
      type(Hstat)::stats
      type(hssbf)::hss_bf1
      character(len=1024)  :: strings
      type(proctree)::ptree

      Maxgrp = 2**(ptree%nlevel) - 1

      msh%basis_group(1)%pgno = 1
      do level = 0, hss_bf1%Maxlevel
         do group = 2**level, 2**(level + 1) - 1
            if (level < hss_bf1%Maxlevel) then
            if (msh%basis_group(group)%pgno*2 <= Maxgrp) then
               msh%basis_group(2*group)%pgno = msh%basis_group(group)%pgno*2
            else
               msh%basis_group(2*group)%pgno = msh%basis_group(group)%pgno
            endif
            if (msh%basis_group(group)%pgno*2 + 1 <= Maxgrp) then
               msh%basis_group(2*group + 1)%pgno = msh%basis_group(group)%pgno*2 + 1
            else
               msh%basis_group(2*group + 1)%pgno = msh%basis_group(group)%pgno
            endif
            endif
         enddo
      enddo

      hss_bf1%N = msh%Nunk
      hss_bf1%BP%level = 0
      hss_bf1%BP%col_group = 1
      hss_bf1%BP%row_group = 1
      hss_bf1%BP%pgno = 1

      allocate (hss_bf1%BP%LL(LplusMax))
      do ll = 1, LplusMax
         hss_bf1%BP%LL(ll)%Nbound = 0
      end do

      hss_bf1%BP%LL(1)%Nbound = 1
      allocate (hss_bf1%BP%LL(1)%matrices_block(1))
      block_f => hss_bf1%BP%LL(1)%matrices_block(1)
      block_f%level = hss_bf1%BP%level
      block_f%level_butterfly = int((hss_bf1%Maxlevel - block_f%level)/2)*2 ! butterfly plus needs even number of levels

      block_f%col_group = hss_bf1%BP%col_group
      block_f%row_group = hss_bf1%BP%row_group
      block_f%pgno = hss_bf1%BP%pgno
      ! pgno_bplus=block_f%pgno

      block_f%M = msh%basis_group(block_f%row_group)%tail - msh%basis_group(block_f%row_group)%head + 1
      block_f%N = msh%basis_group(block_f%col_group)%tail - msh%basis_group(block_f%col_group)%head + 1
      block_f%headm = msh%basis_group(block_f%row_group)%head
      block_f%headn = msh%basis_group(block_f%col_group)%head

      call ComputeParallelIndices(block_f, block_f%pgno, ptree, msh)

      block_f%style = 2
      allocate (hss_bf1%BP%LL(1)%boundary_map(1))
      hss_bf1%BP%LL(1)%boundary_map(1) = block_f%col_group
      hss_bf1%BP%Lplus = 0

      do ll = 1, LplusMax - 1
         if (hss_bf1%BP%LL(ll)%Nbound > 0) then
            hss_bf1%BP%Lplus = hss_bf1%BP%Lplus + 1
            call assert(hss_bf1%BP%Lplus <= LplusMax, 'increase LplusMax')

            block_f => hss_bf1%BP%LL(ll)%matrices_block(1)

            if (ll == LplusMax - 1 .or. block_f%level_butterfly == 0) then
               hss_bf1%BP%LL(ll + 1)%Nbound = 0
            else
               level_butterfly = block_f%level_butterfly
               level_BP = hss_bf1%BP%level
               levelm = ceiling_safe(dble(level_butterfly)/2d0)
               groupm_start = block_f%row_group*2**levelm
               Nboundall = 2**(block_f%level + levelm - level_BP)
               allocate (hss_bf1%BP%LL(ll + 1)%boundary_map(Nboundall))
               do bb = 1, Nboundall
                  hss_bf1%BP%LL(ll + 1)%boundary_map(bb) = bb + groupm_start - 1
               enddo
               hss_bf1%BP%LL(ll + 1)%Nbound = Nboundall

               allocate (hss_bf1%BP%LL(ll + 1)%matrices_block(hss_bf1%BP%LL(ll + 1)%Nbound))
               cnt = 0
               do bb = 1, Nboundall
                  if (hss_bf1%BP%LL(ll + 1)%boundary_map(bb) /= -1) then
                     cnt = cnt + 1
                     group_m = bb + groupm_start - 1
                     group_n = hss_bf1%BP%LL(ll + 1)%boundary_map(bb)
                     blocks => hss_bf1%BP%LL(ll + 1)%matrices_block(cnt)
                     blocks%row_group = group_m
                     blocks%col_group = group_n
                     blocks%level = GetTreelevel(group_m) - 1
                     blocks%level_butterfly = int((hss_bf1%Maxlevel - blocks%level)/2)*2
                     blocks%pgno = msh%basis_group(group_m)%pgno

                     blocks%M = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
                     blocks%N = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
                     blocks%headm = msh%basis_group(group_m)%head
                     blocks%headn = msh%basis_group(group_n)%head

                     call ComputeParallelIndices(blocks, blocks%pgno, ptree, msh)
                     if (blocks%level_butterfly > 0) then
                        blocks%style = 2
                     else
                        blocks%style = 1  ! leaflevel or leaflevel+1 is dense
                        if(ptree%pgrp(blocks%pgno)%nproc>1)then
                           write(*,*)'more than one process sharing a dense block, try to reduce number of processes'
                           stop
                        endif
                     endif

                  end if
               end do
            end if
         else
            exit
         end if
      end do

      call Bplus_copy(hss_bf1%BP, hss_bf1%BP_inverse)
      call LogMemory(stats, SIZEOF(hss_bf1%BP)/1024.0d3 + SIZEOF(hss_bf1%BP_inverse)/1024.0d3)


      msh%idxs = hss_bf1%BP%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1, 1)
      msh%idxe = hss_bf1%BP%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1, 2)

      if (allocated(msh%xyz)) then
         call LogMemory(stats, - SIZEOF(msh%xyz)/1024.0d3)
         ! deallocate (msh%xyz)
      endif

   end subroutine HSS_structuring

   subroutine HODLR_structuring(ho_bf1, option, msh, ker, ptree, stats)
      use BPACK_DEFS
      use MISC_Utilities
      implicit none

      integer i, j, ii, jj, kk, iii, jjj, ll, bb, sortdirec, ii_sch, pgno_bplus
      integer level, edge, patch, node, group, group_touch
      integer rank, index_near, m, n, length, flag, itemp, cnt, detection
      real T0
      real(kind=8):: tolerance, rtemp, rel_error, seperator, dist
      real(kind=8) Memory_direct_forward, Memory_butterfly_forward
      integer mm, nn, header_m, header_n, edge_m, edge_n, group_m, group_n, group_m1, group_n1, group_m2, group_n2, levelm, groupm_start, index_i_m, index_j_m
      integer level_c, iter, level_cc, level_BP, Nboundall, level_butterfly
      type(matrixblock), pointer::blocks, block_f, block_sch, block_inv
      real(kind=8)::minbound, theta, phi, r, rmax, phi_tmp, measure
      real(kind=8), allocatable::Centroid_M(:, :), Centroid_N(:, :)
      integer, allocatable::Isboundary_M(:), Isboundary_N(:)
      integer Dimn, col_group, row_group, Maxgrp
      type(Hoption)::option
      type(mesh)::msh
      type(kernelquant)::ker
      type(Hstat)::stats
      type(hobf)::ho_bf1
      character(len=1024)  :: strings
      type(proctree)::ptree

      Maxgrp = 2**(ptree%nlevel) - 1

      msh%basis_group(1)%pgno = 1
      do level = 0, ho_bf1%Maxlevel
         do group = 2**level, 2**(level + 1) - 1
            if (level < ho_bf1%Maxlevel) then
            if (msh%basis_group(group)%pgno*2 <= Maxgrp) then
               msh%basis_group(2*group)%pgno = msh%basis_group(group)%pgno*2
            else
               msh%basis_group(2*group)%pgno = msh%basis_group(group)%pgno
            endif
            if (msh%basis_group(group)%pgno*2 + 1 <= Maxgrp) then
               msh%basis_group(2*group + 1)%pgno = msh%basis_group(group)%pgno*2 + 1
            else
               msh%basis_group(2*group + 1)%pgno = msh%basis_group(group)%pgno
            endif
            endif
         enddo
      enddo

      ho_bf1%N = msh%Nunk
      allocate (ho_bf1%levels(ho_bf1%Maxlevel + 1))
      call LogMemory(stats, SIZEOF(ho_bf1%levels)/1024.0d3)


      do level_c = 1, ho_bf1%Maxlevel + 1
         ho_bf1%levels(level_c)%level = level_c
         if (level_c == ho_bf1%Maxlevel + 1) then
            ho_bf1%levels(level_c)%N_block_forward = 2**(level_c - 1)
         else
            ho_bf1%levels(level_c)%N_block_forward = 2**level_c
         endif
         ho_bf1%levels(level_c)%N_block_inverse = 2**(level_c - 1)
         ho_bf1%levels(level_c)%Bidxs = 2**(ho_bf1%Maxlevel + 1)
         ho_bf1%levels(level_c)%Bidxe = -2**(ho_bf1%Maxlevel + 1)

         allocate (ho_bf1%levels(level_c)%BP(ho_bf1%levels(level_c)%N_block_forward))
         call LogMemory(stats, SIZEOF(ho_bf1%levels(level_c)%BP)/1024.0d3)
         allocate (ho_bf1%levels(level_c)%BP_inverse(ho_bf1%levels(level_c)%N_block_inverse))
         call LogMemory(stats, SIZEOF(ho_bf1%levels(level_c)%BP_inverse)/1024.0d3)
         allocate (ho_bf1%levels(level_c)%BP_inverse_update(ho_bf1%levels(level_c)%N_block_forward))
         call LogMemory(stats, SIZEOF(ho_bf1%levels(level_c)%BP_inverse_update)/1024.0d3)
         allocate (ho_bf1%levels(level_c)%BP_inverse_schur(ho_bf1%levels(level_c)%N_block_inverse))
         call LogMemory(stats, SIZEOF(ho_bf1%levels(level_c)%BP_inverse_schur)/1024.0d3)
      end do

      ho_bf1%levels(1)%BP_inverse(1)%level = 0
      ho_bf1%levels(1)%BP_inverse(1)%col_group = 1
      ho_bf1%levels(1)%BP_inverse(1)%row_group = 1
      ho_bf1%levels(1)%BP_inverse(1)%pgno = 1
      ! ho_bf1%levels(1)%BP_inverse(1)%style = 2

      ! treat hodlr as a full matrix if Maxlevel=0
      if (ho_bf1%Maxlevel == 0) then
         ho_bf1%levels(1)%BP(1)%level = 0
         ho_bf1%levels(1)%BP(1)%col_group = 1
         ho_bf1%levels(1)%BP(1)%row_group = 1
         ho_bf1%levels(1)%BP(1)%pgno = 1
      endif

      do level_c = 1, ho_bf1%Maxlevel
         do ii = 1, ho_bf1%levels(level_c)%N_block_inverse
            col_group = ho_bf1%levels(level_c)%BP_inverse(ii)%col_group
            row_group = ho_bf1%levels(level_c)%BP_inverse(ii)%row_group

            allocate (ho_bf1%levels(level_c)%BP_inverse(ii)%LL(LplusMax))
            do ll = 1, LplusMax
               ho_bf1%levels(level_c)%BP_inverse(ii)%LL(ll)%Nbound = 0
            end do
            ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%Nbound = 1

            allocate (ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1))
            block_inv => ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)
            block_inv%col_group = col_group
            block_inv%row_group = row_group
            block_inv%level = ho_bf1%levels(level_c)%BP_inverse(ii)%level
            block_inv%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
            block_inv%headm = msh%basis_group(row_group)%head
            block_inv%headn = msh%basis_group(col_group)%head
            block_inv%M = msh%basis_group(row_group)%tail - msh%basis_group(row_group)%head + 1
            block_inv%N = msh%basis_group(col_group)%tail - msh%basis_group(col_group)%head + 1

            block_inv%level_butterfly = ho_bf1%Maxlevel - block_inv%level

            call ComputeParallelIndices(block_inv, block_inv%pgno, ptree, msh)

            if (IOwnPgrp(ptree, block_inv%pgno)) then
               ho_bf1%levels(level_c)%Bidxs = min(ho_bf1%levels(level_c)%Bidxs, ii)
               ho_bf1%levels(level_c)%Bidxe = max(ho_bf1%levels(level_c)%Bidxe, ii)
            endif

            if (GetTreelevel(ho_bf1%levels(level_c)%BP_inverse(ii)%pgno) == ptree%nlevel) then
               ho_bf1%levels(level_c)%BP(ii*2 - 1)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
               ho_bf1%levels(level_c)%BP(ii*2)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
               ho_bf1%levels(level_c + 1)%BP_inverse(ii*2 - 1)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
               ho_bf1%levels(level_c + 1)%BP_inverse(ii*2)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
               ho_bf1%levels(level_c)%BP_inverse_schur(ii)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
            else
               ho_bf1%levels(level_c)%BP(ii*2 - 1)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno*2
               ho_bf1%levels(level_c)%BP(ii*2)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno*2 + 1
               ho_bf1%levels(level_c + 1)%BP_inverse(ii*2 - 1)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno*2
               ho_bf1%levels(level_c + 1)%BP_inverse(ii*2)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno*2 + 1
               ho_bf1%levels(level_c)%BP_inverse_schur(ii)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
            endif

            ! off-diagonal blocks and their updates
            ho_bf1%levels(level_c)%BP(ii*2 - 1)%level = level_c
            ho_bf1%levels(level_c)%BP(ii*2 - 1)%col_group = col_group*2 + 1
            ho_bf1%levels(level_c)%BP(ii*2 - 1)%row_group = row_group*2
            ho_bf1%levels(level_c)%BP(ii*2)%level = level_c
            ho_bf1%levels(level_c)%BP(ii*2)%col_group = col_group*2
            ho_bf1%levels(level_c)%BP(ii*2)%row_group = row_group*2 + 1

            ! schur complement of every two off-diagonal blocks
            ho_bf1%levels(level_c)%BP_inverse_schur(ii)%level = level_c + 1
            ho_bf1%levels(level_c)%BP_inverse_schur(ii)%col_group = col_group*2
            ho_bf1%levels(level_c)%BP_inverse_schur(ii)%row_group = row_group*2
            ho_bf1%levels(level_c)%BP_inverse_schur(ii)%Lplus = 1

            ! diagonal blocks and their inverses at bottom level
            if (level_c == ho_bf1%Maxlevel) then
               ho_bf1%levels(level_c + 1)%BP(ii*2 - 1)%level = level_c + 1
               ho_bf1%levels(level_c + 1)%BP(ii*2 - 1)%col_group = col_group*2
               ho_bf1%levels(level_c + 1)%BP(ii*2 - 1)%row_group = row_group*2
               ho_bf1%levels(level_c + 1)%BP(ii*2)%level = level_c + 1
               ho_bf1%levels(level_c + 1)%BP(ii*2)%col_group = col_group*2 + 1
               ho_bf1%levels(level_c + 1)%BP(ii*2)%row_group = row_group*2 + 1
               if (GetTreelevel(ho_bf1%levels(level_c)%BP_inverse(ii)%pgno) == ptree%nlevel) then
                  ho_bf1%levels(level_c + 1)%BP(ii*2 - 1)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
                  ho_bf1%levels(level_c + 1)%BP(ii*2)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
               else
                  ho_bf1%levels(level_c + 1)%BP(ii*2 - 1)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno*2
                  ho_bf1%levels(level_c + 1)%BP(ii*2)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno*2 + 1
               endif
            end if

            ho_bf1%levels(level_c + 1)%BP_inverse(ii*2 - 1)%level = level_c
            ho_bf1%levels(level_c + 1)%BP_inverse(ii*2 - 1)%col_group = col_group*2
            ho_bf1%levels(level_c + 1)%BP_inverse(ii*2 - 1)%row_group = row_group*2
            ho_bf1%levels(level_c + 1)%BP_inverse(ii*2)%level = level_c
            ho_bf1%levels(level_c + 1)%BP_inverse(ii*2)%col_group = col_group*2 + 1
            ho_bf1%levels(level_c + 1)%BP_inverse(ii*2)%row_group = row_group*2 + 1
         end do
      end do

      ! do level_c = 1,ho_bf1%Maxlevel+1
      ! deallocate(ho_bf1%levels(level_c)%BP_inverse)
      ! enddo

      Dimn = 0
      if (allocated(msh%xyz)) Dimn = size(msh%xyz, 1)

      do level_c = 1, ho_bf1%Maxlevel + 1
         do ii = 1, ho_bf1%levels(level_c)%N_block_forward
            ! if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(ii)%pgno))then
            if (level_c == ho_bf1%Maxlevel + 1) then

               ! bottom level dense blocks
               ho_bf1%levels(level_c)%BP(ii)%Lplus = 1
               allocate (ho_bf1%levels(level_c)%BP(ii)%LL(LplusMax))
               do ll = 1, LplusMax
                  ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound = 0
               end do
               ho_bf1%levels(level_c)%BP(ii)%LL(1)%Nbound = 1
               allocate (ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1))
               block_f => ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)
               block_f%level = ho_bf1%levels(level_c)%BP(ii)%level
               block_f%col_group = ho_bf1%levels(level_c)%BP(ii)%col_group
               block_f%row_group = ho_bf1%levels(level_c)%BP(ii)%row_group
               block_f%style = 1  !!!!! be careful here
               block_f%pgno = msh%basis_group(block_f%row_group)%pgno

               block_f%M = msh%basis_group(block_f%row_group)%tail - msh%basis_group(block_f%row_group)%head + 1
               block_f%N = msh%basis_group(block_f%col_group)%tail - msh%basis_group(block_f%col_group)%head + 1
               block_f%headm = msh%basis_group(block_f%row_group)%head
               block_f%headn = msh%basis_group(block_f%col_group)%head
               block_f%level_butterfly = 0
               ! call ComputeParallelIndices(block_f,block_f%pgno,ptree,msh)  ! block_f%level= Maxlevel+1 causes a bug in ComputeParallelIndices
               block_f%M_loc = block_f%M
               block_f%N_loc = block_f%N
               allocate (block_f%M_p(1, 2))
               block_f%M_p(1, 1) = 1
               block_f%M_p(1, 2) = block_f%M
               allocate (block_f%N_p(1, 2))
               block_f%N_p(1, 1) = 1
               block_f%N_p(1, 2) = block_f%N

               ! bottom level dense blocks' inverse
               ho_bf1%levels(level_c)%BP_inverse(ii)%Lplus = 1
               allocate (ho_bf1%levels(level_c)%BP_inverse(ii)%LL(LplusMax))
               do ll = 1, LplusMax
                  ho_bf1%levels(level_c)%BP_inverse(ii)%LL(ll)%Nbound = 0
               end do
               ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%Nbound = 1
               allocate (ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1))
               block_inv => ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)
               block_inv%level = ho_bf1%levels(level_c)%BP_inverse(ii)%level
               block_inv%col_group = ho_bf1%levels(level_c)%BP_inverse(ii)%col_group
               block_inv%row_group = ho_bf1%levels(level_c)%BP_inverse(ii)%row_group
               block_inv%style = 1  !!!!! be careful here
               block_inv%pgno = msh%basis_group(block_inv%row_group)%pgno

               block_inv%M = msh%basis_group(block_inv%row_group)%tail - msh%basis_group(block_inv%row_group)%head + 1
               block_inv%N = msh%basis_group(block_inv%col_group)%tail - msh%basis_group(block_inv%col_group)%head + 1
               block_inv%headm = msh%basis_group(block_inv%row_group)%head
               block_inv%headn = msh%basis_group(block_inv%col_group)%head
               block_inv%level_butterfly = 0
               call ComputeParallelIndices(block_inv, block_inv%pgno, ptree, msh)
               if (IOwnPgrp(ptree, block_inv%pgno)) then
                  ho_bf1%levels(level_c)%Bidxs = min(ho_bf1%levels(level_c)%Bidxs, ii)
                  ho_bf1%levels(level_c)%Bidxe = max(ho_bf1%levels(level_c)%Bidxe, ii)
               endif
            else
               allocate (ho_bf1%levels(level_c)%BP(ii)%LL(LplusMax))
               do ll = 1, LplusMax
                  ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound = 0
               end do

               ho_bf1%levels(level_c)%BP(ii)%LL(1)%Nbound = 1
               allocate (ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1))
               block_f => ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)
               block_f%level = ho_bf1%levels(level_c)%BP(ii)%level

               if (level_c > option%LRlevel) then
                  block_f%level_butterfly = 0 ! low rank below LRlevel
               else
                  if (ho_bf1%Maxlevel - block_f%level < option%lnoBP) then
                     block_f%level_butterfly = ho_bf1%Maxlevel - ho_bf1%levels(level_c)%BP(ii)%level   ! butterfly
                  else
                     block_f%level_butterfly = int((ho_bf1%Maxlevel - ho_bf1%levels(level_c)%BP(ii)%level)/2)*2 ! butterfly plus needs even number of levels
                  endif
               endif

               block_f%col_group = ho_bf1%levels(level_c)%BP(ii)%col_group
               block_f%row_group = ho_bf1%levels(level_c)%BP(ii)%row_group
               block_f%pgno = msh%basis_group(block_f%row_group)%pgno
               pgno_bplus = block_f%pgno

               ! compute the partial indices when BP is shared by double number of processes
               ii_sch = ceiling_safe(ii/2d0)
               block_inv => ho_bf1%levels(level_c)%BP_inverse(ii_sch)%LL(1)%matrices_block(1)
               block_f%pgno_db = block_inv%pgno

               block_f%M = msh%basis_group(block_f%row_group)%tail - msh%basis_group(block_f%row_group)%head + 1
               block_f%N = msh%basis_group(block_f%col_group)%tail - msh%basis_group(block_f%col_group)%head + 1
               block_f%headm = msh%basis_group(block_f%row_group)%head
               block_f%headn = msh%basis_group(block_f%col_group)%head

               call ComputeParallelIndices(block_f, block_f%pgno, ptree, msh)
               ! call ComputeParallelIndices(block_f,block_f%pgno_db,ptree,msh,1)
               ! if(block_f%M==2500)write(*,*)ptree%myID,block_f%pgno,block_f%pgno_db,block_f%N_loc,block_f%N_loc_db,'eref'

               block_f%style = 2
               allocate (ho_bf1%levels(level_c)%BP(ii)%LL(1)%boundary_map(1))
               ho_bf1%levels(level_c)%BP(ii)%LL(1)%boundary_map(1) = block_f%col_group
               ho_bf1%levels(level_c)%BP(ii)%Lplus = 0

               group = floor((ii - 1 + 2**level_c)/2d0)
               sortdirec = NINT(msh%basis_group(group)%boundary(1))
               seperator = msh%basis_group(group)%boundary(2)

               do ll = 1, LplusMax - 1
                  if (ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound > 0) then
                     ho_bf1%levels(level_c)%BP(ii)%Lplus = ho_bf1%levels(level_c)%BP(ii)%Lplus + 1
                     call assert(ho_bf1%levels(level_c)%BP(ii)%Lplus <= LplusMax, 'increase LplusMax')
                     ! write(*,*)'nini',level_c,ho_bf1%Maxlevel - ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(1)%level,option%lnoBP,ll

                     block_f => ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(1)

                     if (ho_bf1%Maxlevel - block_f%level < option%lnoBP .or. ll == LplusMax - 1 .or. block_f%level_butterfly == 0) then
                        ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%Nbound = 0
                     else
                        ! write(*,*)'gggggg'
                        ! level_butterfly = int((ho_bf1%Maxlevel - ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(1)%level)/2)*2
                        level_butterfly = block_f%level_butterfly
                        level_BP = ho_bf1%levels(level_c)%BP(ii)%level
                        levelm = ceiling_safe(dble(level_butterfly)/2d0)
                        groupm_start = block_f%row_group*2**levelm
                        Nboundall = 2**(block_f%level + levelm - level_BP)
                        allocate (ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%boundary_map(Nboundall))
                        ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%boundary_map = -1
                        ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%Nbound = 0

                        do bb = 1, ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound
                           blocks => ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(bb)

                           allocate (Centroid_M(2**levelm, Dimn))
                           allocate (Isboundary_M(2**levelm))
                           Isboundary_M = 0
                           Centroid_M = 0

                           do index_i_m = 1, 2**levelm
                              group_m = blocks%row_group   ! Note: row_group and col_group interchanged here
                              group_m = group_m*2**levelm - 1 + index_i_m

                              CNT = 0
                              if (option%xyzsort == CKD) then
                                 do nn = msh%basis_group(group_m)%head, msh%basis_group(group_m)%tail
                                    measure = abs(msh%xyz(sortdirec, msh%new2old(nn)) - seperator)
                                    if (measure < option%touch_para) then
                                       Isboundary_M(index_i_m) = 1
                                       CNT = CNT + 1
                                       Centroid_M(index_i_m, 1:Dimn) = Centroid_M(index_i_m, 1:Dimn) + msh%xyz(1:Dimn, msh%new2old(nn))
                                    end if
                                 end do
                                 if (Isboundary_M(index_i_m) == 1) Centroid_M(index_i_m, :) = Centroid_M(index_i_m, :)/CNT

                                 ! if(blocks%col_group==8 .or. blocks%col_group==9)then
                                 ! write(*,*)'wocaoo',group_m,Isboundary_M(index_i_m),CNT,sortdirec,seperator
                                 ! endif

                              end if
                           end do

                           allocate (Centroid_N(2**(level_butterfly - levelm), Dimn))
                           allocate (Isboundary_N(2**(level_butterfly - levelm)))
                           Isboundary_N = 0
                           Centroid_N = 0

                           do index_j_m = 1, 2**(level_butterfly - levelm)
                              group_n = blocks%col_group
                              group_n = group_n*2**(level_butterfly - levelm) - 1 + index_j_m

                              CNT = 0
                              if (option%xyzsort == CKD) then
                                 do nn = msh%basis_group(group_n)%head, msh%basis_group(group_n)%tail
                                    measure = abs(msh%xyz(sortdirec, msh%new2old(nn)) - seperator)
                                    if (measure < option%touch_para) then
                                       Isboundary_N(index_j_m) = 1
                                       CNT = CNT + 1

                                       Centroid_N(index_j_m, 1:Dimn) = Centroid_N(index_j_m, 1:Dimn) + msh%xyz(1:Dimn, msh%new2old(nn))
                                    end if
                                 end do
                                 if (Isboundary_N(index_j_m) == 1) Centroid_N(index_j_m, :) = Centroid_N(index_j_m, :)/CNT
                              end if
                           end do

                           ! if(level_c==1)then
                           ! ! write(*,*)Isboundary_N,Isboundary_M
                           ! do kk=1,2**levelm
                           ! if(Isboundary_M(kk)==1)then
                           ! ! write(*,*)Centroid_M(kk,1),Centroid_M(kk,2),Centroid_M(kk,3)
                           ! write(777,*)Centroid_M(kk,1),Centroid_M(kk,2),Centroid_M(kk,3)
                           ! end if
                           ! end do

                           ! do kk=1,2**(level_butterfly-levelm)
                           ! if(Isboundary_N(kk)==1)then
                           ! write(777,*)Centroid_N(kk,1),Centroid_N(kk,2),Centroid_N(kk,3)
                           ! end if
                           ! end do
                           ! end if

                           do index_i_m = 1, 2**levelm
                              group_m = blocks%row_group   ! Note: row_group and col_group interchanged here
                              group_m = group_m*2**levelm - 1 + index_i_m

                              if (Isboundary_M(index_i_m) == 1) then
                                 ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%Nbound = ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%Nbound + 1
                                 dist = 100000000d0
                                 do index_j_m = 1, 2**(level_butterfly - levelm)
                                    group_n = blocks%col_group
                                    group_n = group_n*2**(level_butterfly - levelm) - 1 + index_j_m

                                    ! if(blocks%col_group==8 .or. blocks%col_group==9)then
                                    ! write(*,*)group_m,group_n,sqrt(sum((Centroid_N(index_j_m,:)-Centroid_M(index_i_m,:))**2d0)),'nima'
                                    ! end        if

                                    if (Isboundary_N(index_j_m) == 1) then
                                       if (dist > sqrt(sum((Centroid_N(index_j_m, :) - Centroid_M(index_i_m, :))**2d0))) then
                                          ! if(level_c==1)write(*,*)index_i_m,index_j_m
                                          dist = sqrt(sum((Centroid_N(index_j_m, :) - Centroid_M(index_i_m, :))**2d0))
                                          ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%boundary_map(group_m - groupm_start + 1) = group_n
                                       end if
                                    end if
                                 end do
                              end if
                           enddo
                           deallocate (Isboundary_M)
                           deallocate (Isboundary_N)
                           deallocate (Centroid_M)
                           deallocate (Centroid_N)

                        end do

                        if (ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%Nbound > 1) then
                           ! write(*,*)level_c,ii,ll,ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%col_group,ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%row_group,ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%Nbound,'niamaa'
                        endif

                        call assert(ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%Nbound > 0, 'why is no boundary group detected')

                        allocate (ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%matrices_block(ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%Nbound))

                        cnt = 0
                        do bb = 1, Nboundall
                           if (ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%boundary_map(bb) /= -1) then
                              cnt = cnt + 1
                              group_m = bb + groupm_start - 1
                              group_n = ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%boundary_map(bb)
                              blocks => ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%matrices_block(cnt)
                              blocks%row_group = group_m
                              blocks%col_group = group_n
                              blocks%level = GetTreelevel(group_m) - 1
                              ! blocks%level_butterfly = int((ho_bf1%Maxlevel - blocks%level)/2)*2
                              blocks%level_butterfly = 0 ! only two layer butterfly plus here

                              blocks%pgno = msh%basis_group(group_m)%pgno
                              do while (blocks%pgno > pgno_bplus)
                                 if (level_butterfly < ptree%nlevel - GetTreelevel(blocks%pgno)) exit
                                 blocks%pgno = blocks%pgno/2
                              enddo

                              blocks%pgno_db = blocks%pgno
                              blocks%M = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
                              blocks%N = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
                              blocks%headm = msh%basis_group(group_m)%head
                              blocks%headn = msh%basis_group(group_n)%head

                              blocks%style = 2
                              call ComputeParallelIndices(blocks, blocks%pgno, ptree, msh)
                              ! call ComputeParallelIndices(blocks,blocks%pgno_db,ptree,msh,1)
                           end if
                        end do
                     end if
                  else
                     exit
                  end if
               end do

               ! write(*,*)level_c,ii,ho_bf1%levels(level_c)%BP(ii)%Lplus,'gaogao '

               if (mod(ii, 2) == 1) then  ! in the beginning only even block hold information about the schur complement
                  ii_sch = ceiling_safe(ii/2d0)

                  allocate (ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(LplusMax))
                  ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%Lplus = ho_bf1%levels(level_c)%BP(ii)%Lplus
                  do ll = 1, LplusMax
                     ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%Nbound = 0
                  end do

                  ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(1)%Nbound = 1

                  allocate (ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(1)%boundary_map(1))
                  ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(1)%boundary_map(1) = ho_bf1%levels(level_c)%BP(ii)%row_group

                  do ll = 1, LplusMax - 1
                     if (ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound > 0) then

                        ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%rankmax = 0
                        ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%Nbound = ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound

                        allocate (ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%Nbound))

                        do bb = 1, ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%Nbound
                           block_f => ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(bb)
                           block_sch => ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(bb)

                           row_group = block_f%row_group

                           block_sch%row_group = row_group
                           block_sch%col_group = row_group

                           if (msh%basis_group(row_group)%pgno /= msh%basis_group(INT(row_group/2d0))%pgno) then
                              block_sch%pgno = msh%basis_group(INT(row_group/2d0))%pgno
                           else
                              block_sch%pgno = msh%basis_group(row_group)%pgno
                           end if

                           block_sch%style = block_f%style
                           block_sch%level = block_f%level
                           block_sch%level_butterfly = block_f%level_butterfly

                           block_sch%M = msh%basis_group(row_group)%tail - msh%basis_group(row_group)%head + 1
                           block_sch%N = msh%basis_group(row_group)%tail - msh%basis_group(row_group)%head + 1
                           block_sch%headm = msh%basis_group(row_group)%head
                           block_sch%headn = msh%basis_group(row_group)%head
                           call ComputeParallelIndices(block_sch, block_sch%pgno, ptree, msh)
                        end do

                        if (ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%Nbound == 0) then
                           ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll + 1)%Nbound = 0
                        else
                           level_butterfly = ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(1)%level_butterfly
                           level_BP = ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%level
                           levelm = ceiling_safe(dble(level_butterfly)/2d0)
                           groupm_start = ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(1)%row_group*2**levelm
                           Nboundall = 2**(ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(1)%level + levelm - level_BP)

                           allocate (ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll + 1)%boundary_map(Nboundall))

                           ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll + 1)%boundary_map = ho_bf1%levels(level_c)%BP(ii)%LL(ll + 1)%boundary_map
                           do bb = 1, Nboundall
                              if (ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll + 1)%boundary_map(bb) /= -1) then
                                 ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll + 1)%boundary_map(bb) = bb + groupm_start - 1
                              end if
                           end do
                        end if
                     else
                        exit
                     end if
                  end do

               end if

               ! ! if(level_c==1 .and. ii==1)then

               ! write(strings , *) 2*dimn
               ! ! write(177,*)'Bplus:', level_c,ii
               ! do ll=1,ho_bf1%levels(level_c)%BP(ii)%Lplus
               ! ! write(*,*)ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound,'ddd'
               ! do bb = 1,ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound
               ! write(177,'(I3,I7,I3,I3,'//TRIM(strings)//'Es16.7)')level_c,ii,ll,ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%level,msh%basis_group(ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%row_group)%center(1:dimn),msh%basis_group(ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%col_group)%center(1:dimn)
               ! end do
               ! end do
               ! ! end if
            end if
            ! end if

            call Bplus_copy(ho_bf1%levels(level_c)%BP(ii), ho_bf1%levels(level_c)%BP_inverse_update(ii))

         end do
      end do

      msh%idxs = ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1, 1)
      msh%idxe = ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1, 2)

      if (allocated(msh%xyz)) then
         call LogMemory(stats, - SIZEOF(msh%xyz)/1024.0d3)
         ! deallocate (msh%xyz)
      endif

   end subroutine HODLR_structuring

   subroutine Hmat_structuring(h_mat, option, msh, ker, ptree, stats)
      use BPACK_DEFS
      use MISC_Utilities
      implicit none

      type(Hmat)::h_mat
      type(Hoption)::option
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(Hstat)::stats
      integer i, j, ii, jj, iii, jjj, k, kk, kkk
      integer level, edge, node, patch, group, group_m, group_n
      integer mm, nn, num_blocks,group_start,gg
      type(matrixblock), pointer :: blocks
      type(matrixblock) :: blocks_dummy
      integer Maxgrp, ierr, row_group, col_group
      type(global_matricesblock), pointer::global_block
      integer nprow,npcol,myrow,mycol,myArows,myAcols,mypgno

      msh%idxs = 1000000000
      msh%idxe = -1000000000
      h_mat%myArows=0
      h_mat%myAcols=0

      nprow = ptree%pgrp(1)%nprow
      npcol = ptree%pgrp(1)%npcol
      ii = max(nprow,npcol)
      level = 0
      do while (2**level<ii)
         level = level + 1
      enddo
      msh%Dist_level = level
      h_mat%Dist_level = level

      Maxgrp = 2**(ptree%nlevel) - 1
      msh%basis_group(1)%pgno = 1
      do level = 0, h_mat%Maxlevel
         do group = 2**level, 2**(level + 1) - 1
            if (level < h_mat%Maxlevel) then
               if (msh%basis_group(group)%pgno*2 <= Maxgrp) then
                  msh%basis_group(2*group)%pgno = msh%basis_group(group)%pgno*2
               else
                  msh%basis_group(2*group)%pgno = msh%basis_group(group)%pgno
               endif
               if (msh%basis_group(group)%pgno*2 + 1 <= Maxgrp) then
                  msh%basis_group(2*group + 1)%pgno = msh%basis_group(group)%pgno*2 + 1
               else
                  msh%basis_group(2*group + 1)%pgno = msh%basis_group(group)%pgno
               endif
            else
               if(ptree%pgrp(msh%basis_group(group)%pgno)%head==ptree%MyID)then
                  msh%idxs = min(msh%idxs,msh%basis_group(group)%head)
                  msh%idxe = max(msh%idxe,msh%basis_group(group)%tail)
                  mypgno = msh%basis_group(group)%pgno
               endif
            endif
         enddo
      enddo

      h_mat%N = msh%Nunk
      h_mat%idxs = msh%idxs
      h_mat%idxe = msh%idxe

      allocate (h_mat%blocks_root)
      h_mat%blocks_root%level = 0
      h_mat%blocks_root%row_group = 1
      h_mat%blocks_root%col_group = 1
      nullify (h_mat%blocks_root%father)

      blocks_dummy%pgno=1
      blocks_dummy%M=h_mat%N
      blocks_dummy%N=h_mat%N
      blocks_dummy%row_group=1
      blocks_dummy%col_group=1
      blocks_dummy%level=0
      call ComputeParallelIndices(blocks_dummy, blocks_dummy%pgno, ptree, msh)
      allocate(h_mat%N_p(size(blocks_dummy%N_p,1),size(blocks_dummy%N_p,2)))
      h_mat%N_p = blocks_dummy%N_p
      deallocate(blocks_dummy%N_p)
      deallocate(blocks_dummy%M_p)

      num_blocks = 2**msh%Dist_level
      allocate (h_mat%basis_group(num_blocks))
      group_start = num_blocks - 1
      do gg=1,num_blocks
         h_mat%basis_group(gg)%head = msh%basis_group(gg+group_start)%head
         h_mat%basis_group(gg)%tail = msh%basis_group(gg+group_start)%tail
      enddo

      allocate (stats%leafs_of_level(0:h_mat%Maxlevel))
      stats%leafs_of_level = 0

      allocate (stats%Add_random_CNT(0:h_mat%Maxlevel))
      stats%Add_random_CNT = 0
      allocate (stats%Add_random_Time(0:h_mat%Maxlevel))
      stats%Add_random_Time = 0
      allocate (stats%Mul_random_CNT(0:h_mat%Maxlevel))
      stats%Mul_random_CNT = 0
      allocate (stats%Mul_random_Time(0:h_mat%Maxlevel))
      stats%Mul_random_Time = 0
      allocate (stats%XLUM_random_CNT(0:h_mat%Maxlevel))
      stats%XLUM_random_CNT = 0
      allocate (stats%XLUM_random_Time(0:h_mat%Maxlevel))
      stats%XLUM_random_Time = 0

      call Hmat_construct_global_tree(h_mat%blocks_root, msh%Dist_level, stats)

      allocate (h_mat%First_block_eachlevel(0:h_mat%Maxlevel))
      h_mat%First_block_eachlevel(0)%father => h_mat%blocks_root
      global_block => h_mat%blocks_root
      do level = 1, msh%Dist_level
         global_block => global_block%sons(1, 1)
         h_mat%First_block_eachlevel(level)%father => global_block
      enddo


      call blacs_gridinfo(ptree%pgrp(1)%ctxt, nprow, npcol, myrow, mycol)
      if (nprow /= -1 .and. npcol /= -1) then

         num_blocks = 2**msh%Dist_level
         myArows = numroc_wp(num_blocks, 1, myrow, 0, nprow)
         myAcols = numroc_wp(num_blocks, 1, mycol, 0, npcol)
         h_mat%myArows = myArows
         h_mat%myAcols = myAcols

         allocate (h_mat%Local_blocks(myAcols, myArows))
         allocate (h_mat%Local_blocks_copy(myAcols, myArows))
         do i = 1, myArows
            call l2g(i, myrow, num_blocks, nprow, 1, ii)
            do j = 1, myAcols
               call l2g(j, mycol, num_blocks, npcol, 1, jj)
               blocks => h_mat%Local_blocks(j, i)
               blocks%level = msh%Dist_level
               blocks%row_group = num_blocks + ii - 1
               blocks%col_group = num_blocks + jj - 1

               if (blocks%level > option%LRlevel) then
                  blocks%level_butterfly = 0 ! low rank below LRlevel
               else
                  blocks%level_butterfly = h_mat%Maxlevel - blocks%level   ! butterfly
               endif

               nullify (blocks%father)
               row_group = blocks%row_group
               col_group = blocks%col_group
               blocks%pgno = mypgno
               blocks%M = msh%basis_group(row_group)%tail - msh%basis_group(row_group)%head + 1
               blocks%N = msh%basis_group(col_group)%tail - msh%basis_group(col_group)%head + 1
               blocks%headm = msh%basis_group(row_group)%head
               blocks%headn = msh%basis_group(col_group)%head
               call ComputeParallelIndices(blocks, blocks%pgno, ptree, msh)
               call Hmat_construct_local_tree(blocks, option, stats, msh, ker, ptree, h_mat%Maxlevel)
            enddo
         enddo
      endif

      call MPI_allreduce(MPI_IN_PLACE, stats%leafs_of_level(0:h_mat%Maxlevel), h_mat%Maxlevel + 1, MPI_integer, MPI_sum, ptree%Comm, ierr)


      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
         do level = 1, h_mat%Maxlevel
            write (*, *) "Level:", level, stats%leafs_of_level(level)
         enddo
      endif

      call Hmat_assign_admissible(h_mat, 1, 1, 0, option, stats, msh, ker, ptree)
      call Hmat_compute_colorset(h_mat, option, stats, msh, ker, ptree)



      allocate (h_mat%lstblks(0:h_mat%Maxlevel))
      do level = 0, h_mat%Maxlevel
         h_mat%lstblks(level) = list()
      enddo

      do i = 1, h_mat%myArows
         do j = 1, h_mat%myAcols
            blocks => h_mat%Local_blocks(j, i)
            call Hmat_GetBlkLst(blocks, option, stats, msh, ptree, h_mat)
         enddo
      enddo

      do level = 0, h_mat%Maxlevel
         call MergeSort(h_mat%lstblks(level)%head, node_score_block_ptr_row)
      enddo


      !***************************************************************************************

      if (allocated(msh%xyz)) then
         call LogMemory(stats, - SIZEOF(msh%xyz)/1024.0d3)
         ! deallocate (msh%xyz)
      endif

      return

   end subroutine Hmat_structuring




   recursive subroutine Hmat_assign_admissible(h_mat, group_m, group_n, level, option, stats, msh, ker, ptree)

      implicit none

      type(Hmat)::h_mat
      integer group_m, group_n, group_m_c, group_n_c, i, j, k, level,ll, group_start
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(ipair)::p

      if(level==0)then
         allocate(h_mat%admissibles(msh%Maxgroup))
         call LogMemory(stats, SIZEOF(h_mat%admissibles)/1024.0d3)
      endif

      if (level >= msh%Dist_level .and. near_or_far_user(group_m, group_n, msh, option, ker, option%near_para) == 1) then
         p%i=group_n
         p%j=2
         call append(h_mat%admissibles(group_m), p)
         call LogMemory(stats, SIZEOF(p)/1024.0d3)
      else
         p%i=group_n
         p%j=1
         call append(h_mat%admissibles(group_m), p)
         call LogMemory(stats, SIZEOF(p)/1024.0d3)
         if (level == h_mat%Maxlevel) then
         else
            do j = 1, 2
               do i = 1, 2
                  group_m_c = group_m*2+i-1
                  group_n_c = group_n*2+j-1
                  call Hmat_assign_admissible(h_mat, group_m_c, group_n_c, level+1, option, stats, msh, ker, ptree)
               enddo
            enddo
         endif
      endif

      if(level==0)then
         do ll=0,h_mat%Maxlevel
            group_start = 2**ll - 1
            do i = 1, 2**ll
               call MergeSort(h_mat%admissibles(group_start+i)%head, nod_score_ipair)
               ! if(ptree%MyID==Main_ID)write(*,*)'admissibles: ',group_start+i,h_mat%admissibles(group_start+i)%num_nods
            enddo

         enddo

      endif

      return

   end subroutine Hmat_assign_admissible



   subroutine Hmat_compute_colorset(h_mat, option, stats, msh, ker, ptree)

      implicit none

      type(Hmat)::h_mat
      integer group_m, group_n, group_m_c, group_n_c, i, j, k, ii, jj, level,ll, group_start,mm,nn,nnz,num_colors, max_admissibles, ierr
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(list),allocatable::graph_color(:)

      type(nod), pointer::curr, curc
      class(*), pointer::ptrr, ptrc
      integer,allocatable::row_ptr(:),col_ind(:)
      type(ipair)::p

      if(ptree%MyID==Main_ID)then
         allocate(h_mat%colorsets(0:h_mat%Maxlevel))
         call LogMemory(stats, SIZEOF(h_mat%colorsets)/1024.0d3)
         do level=0,h_mat%Maxlevel
            allocate(graph_color(2**level))
            group_start = 2**level - 1
            h_mat%colorsets(level)%idx = 0
            do i = 1, 2**level
               group_m = group_start + i
               curr => h_mat%admissibles(group_m)%head
               do mm = 1, h_mat%admissibles(group_m)%num_nods
                  ptrr=>curr%item
                  select type (ptrr)
                  type is (ipair)
                     ii = ptrr%i - group_start
                     curc => h_mat%admissibles(group_m)%head
                     do nn = 1, h_mat%admissibles(group_m)%num_nods
                        ptrc=>curc%item
                        select type (ptrc)
                        type is (ipair)
                           jj = ptrc%i - group_start
                           call append(graph_color(ii),jj)
                           h_mat%colorsets(level)%idx = max(h_mat%colorsets(level)%idx,ptrc%j) ! record the maximum admissible values at each level
                        end select
                        curc => curc%next
                     enddo
                  end select
                  curr => curr%next
               enddo
            enddo

            allocate(row_ptr(2**level+1))
            row_ptr(1)=1
            nnz=0
            do i = 1, 2**level
               call MergeSortUnique(graph_color(i), nod_score_integer)
               nnz = nnz + graph_color(i)%num_nods
               row_ptr(i+1) = row_ptr(i) + graph_color(i)%num_nods
            enddo

            allocate(col_ind(nnz))
            nnz=0
            do i = 1, 2**level
               curc => graph_color(i)%head
               do nn = 1, graph_color(i)%num_nods
                  ptrc=>curc%item
                  select type (ptrc)
                  type is (integer)
                     nnz = nnz +1
                     col_ind(nnz) = ptrc
                  end select
                  curc => curc%next
               enddo
            enddo

            h_mat%colorsets(level)%num_nods= 2**level
            allocate(h_mat%colorsets(level)%dat(2**level))
            call LogMemory(stats, SIZEOF(h_mat%colorsets(level)%dat)/1024.0d3)
            call get_graph_colors_JP(2**level,row_ptr,col_ind,h_mat%colorsets(level)%dat)
            call MPI_Bcast(h_mat%colorsets(level)%dat, 2**level, MPI_INTEGER, Main_ID, ptree%pgrp(1)%Comm, ierr) ! this broadcast is needed as the JP algorithm is randomized.
            call MPI_Bcast(h_mat%colorsets(level)%idx, 1, MPI_INTEGER, Main_ID, ptree%pgrp(1)%Comm, ierr)

            do i = 1, 2**level
               call list_finalizer(graph_color(i))
            enddo
            deallocate(graph_color)
            deallocate(row_ptr)
            deallocate(col_ind)
         enddo
      else
         allocate(h_mat%colorsets(0:h_mat%Maxlevel))
         call LogMemory(stats, SIZEOF(h_mat%colorsets)/1024.0d3)
         do level=0,h_mat%Maxlevel
            h_mat%colorsets(level)%num_nods= 2**level
            allocate(h_mat%colorsets(level)%dat(2**level))
            call MPI_Bcast(h_mat%colorsets(level)%dat, 2**level, MPI_INTEGER, Main_ID, ptree%pgrp(1)%Comm, ierr) ! this broadcast is needed as the JP algorithm is randomized.
            call MPI_Bcast(h_mat%colorsets(level)%idx, 1, MPI_INTEGER, Main_ID, ptree%pgrp(1)%Comm, ierr)
         enddo
      endif

      do group_m=1,msh%Maxgroup
         call list_finalizer(h_mat%admissibles(group_m))
         call LogMemory(stats, -SIZEOF(p)*h_mat%admissibles(group_m)%num_nods/1024.0d3)
      enddo
      call LogMemory(stats, -SIZEOF(h_mat%admissibles)/1024.0d3)
      deallocate(h_mat%admissibles)


   end subroutine Hmat_compute_colorset


end module BPACK_structure
