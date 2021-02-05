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
   use Bplus_randomizedop

contains

!**** user-defined subroutine to sample a list of intersections from the bmat of Z
   subroutine Zelem_block_Extraction(Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps, quant)
      use BPACK_DEFS
      implicit none

      class(*), pointer :: quant
      integer:: Ninter
      integer:: allrows(:), allcols(:)
      DT::alldat_loc(:)
      integer::colidx(Ninter), rowidx(Ninter), pgidx(Ninter)
      integer::Npmap, pmaps(Npmap, 3)

      select TYPE (quant)
      type is (quant_bmat)
         call BPACK_ExtractElement(quant%bmat, quant%option, quant%msh, quant%stats, quant%ptree, Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps)
      class default
         write (*, *) "unexpected type"
         stop
      end select

      return

   end subroutine Zelem_block_Extraction

!**** Initialization of the construction phase
   ! N is matrix dimension
   ! P is the permutation vector returned
   ! N_loc is the local number of rows/columns
   ! bmat is the meta-data storing the compressed matrix
   ! Coordinates(optional) of dimension dim*N is the array of Cartesian coordinates corresponding to each row or column
   ! clustertree(optional) is an array of leafsizes in a user-provided cluster tree. clustertree has length 2*nl with nl denoting level of the clustertree.
   ! If clustertree is incomplete with 0 element, ButterflyPACK will adjust it to a complete tree and return a modified clustertree.
   ! If the hierarchical matrix has more levels than clustertree, the code will generate more levels according to option%xyzsort, option%nogeo, and option%Nmin_leaf
   ! nns(optional) of dimension option%knn*N is the array of user provided nearest neighbours
   subroutine BPACK_construction_Init(Nunk, Permutation, Nunk_loc, bmat, option, stats, msh, ker, ptree, Coordinates, tree, nns)
      implicit none
      integer Nunk, Ndim
      real(kind=8), optional:: Coordinates(:, :)
      integer, optional:: nns(:, :)

      real(kind=8) para
      real(kind=8) tolerance
      integer nn, mm, Maxlevel, give, need
      integer i, j, k, ii, edge, Dimn, kk
      integer nlevel, level
      integer Permutation(Nunk)
      integer, optional:: tree(:)
      integer Nunk_loc
      integer groupm

      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(Bmatrix)::bmat
      type(proctree)::ptree

      real(kind=8) t1, t2
      character(len=1024)  :: strings
      integer threads_num

      call assert(associated(ker%QuantApp), 'ker%QuantApp is not assigned')
      call assert(associated(ker%FuncZmn) .or. associated(ker%FuncHMatVec), 'neither ker%FuncZmn nor ker%FuncHMatVec is assigned')

      stats%Flop_Fill = 0
      stats%Time_Fill = 0
      stats%Time_Entry = 0
      stats%Time_Entry_Traverse = 0
      stats%Time_Entry_BF = 0
      stats%Time_Entry_Comm = 0

      !**** set thread number here
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'NUMBER_MPI=', ptree%nproc
      threads_num = 1
      CALL getenv("OMP_NUM_THREADS", strings)
      strings = TRIM(strings)
      if (LEN_TRIM(strings) > 0) then
         read (strings, *) threads_num
      endif
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'OMP_NUM_THREADS=', threads_num
      call OMP_set_num_threads(threads_num)

      msh%Nunk = Nunk

      t1 = OMP_get_wtime()
      nlevel = 0
      if (present(tree)) then
         nlevel = ceiling_safe(log(dble(size(tree, 1)))/log(2d0))
         Maxlevel = nlevel
         allocate (msh%pretree(2**Maxlevel))
         msh%pretree(1:2**Maxlevel) = tree(1:2**Maxlevel)

         !**** make 0-element node a 1-element node

         ! write(*,*)'before adjustment:',msh%pretree
         need = 0
         do ii = 1, 2**Maxlevel
            if (msh%pretree(ii) == 0) need = need + 1
         enddo
         do while (need > 0)
            give = ceiling_safe(need/dble(2**Maxlevel - need))
            do ii = 1, 2**Maxlevel
               nn = msh%pretree(ii)
               if (nn > 1) then
                  msh%pretree(ii) = msh%pretree(ii) - min(min(nn - 1, give), need)
                  need = need - min(min(nn - 1, give), need)
               endif
            enddo
         enddo
         do ii = 1, 2**Maxlevel
            if (msh%pretree(ii) == 0) msh%pretree(ii) = 1
         enddo
         ! write(*,*)'after adjustment:',msh%pretree
         tree(1:2**Maxlevel) = msh%pretree(1:2**Maxlevel)
      endif

      !**** copy geometry points if present
      if (option%nogeo == 0) then
         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "User-supplied kernel requiring reorder"
         call assert(present(Coordinates), 'geometry points should be provided if option%nogeo==0')
         Ndim = size(Coordinates, 1)
         Dimn = Ndim
         allocate (msh%xyz(Dimn, 1:msh%Nunk))
         stats%Mem_Peak = stats%Mem_Peak + SIZEOF(msh%xyz)/1024.0d3
         msh%xyz = Coordinates
      endif

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "    "
      t2 = OMP_get_wtime()

      t1 = OMP_get_wtime()
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Hierarchical format......"
      call Cluster_partition(bmat, option, msh, ker, stats, ptree)
      call BPACK_structuring(bmat, option, msh, ker, ptree, stats)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Hierarchical format finished"
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "    "
      t2 = OMP_get_wtime()

      if (option%nogeo == 3 .and. option%knn > 0) then
         call assert(present(nns), 'nearest neighbours should be provided if option%nogeo==3')
         allocate (msh%nns(msh%Nunk, option%knn))
         do ii = 1, msh%Nunk
         do kk = 1, option%knn
            if (nns(kk, msh%new2old(ii)) /= 0) then
               msh%nns(ii, kk) = msh%old2new(nns(kk, msh%new2old(ii)))
            else
               msh%nns(ii, kk) = 0
            endif
         enddo
         enddo
      endif

      !**** return the permutation vector
      select case (option%format)
      case (HODLR)
         msh%idxs = bmat%ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1, 1)
         msh%idxe = bmat%ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1, 2)
      case (HMAT)
         msh%idxs = bmat%h_mat%Local_blocks(1, 1)%headm
         msh%idxe = bmat%h_mat%Local_blocks(1, 1)%headm + bmat%h_mat%Local_blocks(1, 1)%M - 1
      end select

      Nunk_loc = msh%idxe - msh%idxs + 1
      if (ptree%MyID == Main_ID) then
         do edge = 1, Nunk
            Permutation(edge) = msh%new2old(edge)
         enddo
      endif

   end subroutine BPACK_construction_Init

!**** Computation of the full matrix with matrix entry evaluation
   subroutine FULLMAT_Element(option, stats, msh, ker, ptree)

      implicit none

      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      DT, allocatable::fullmat(:, :), fullmat_tmp(:, :)
      integer N_unk_loc, N_unk_loc1, N_unk_locmax, ii, jj, pp
      integer, allocatable::mrange(:), nrange(:)
      integer ierr, passflag, myflag
      type(intersect)::submats(1)

      N_unk_loc = msh%idxe - msh%idxs + 1
      allocate (fullmat(msh%Nunk, N_unk_loc))
      allocate (mrange(msh%Nunk))
      do ii = 1, msh%Nunk
         mrange(ii) = ii
      enddo
      allocate (nrange(N_unk_loc))
      do jj = 1, N_unk_loc
         nrange(jj) = jj - msh%idxs + 1
      enddo

      submats(1)%nr = msh%Nunk
      submats(1)%nc = N_unk_loc
      allocate(submats(1)%rows(submats(1)%nr))
      submats(1)%rows = mrange(1:submats(1)%nr)
      allocate(submats(1)%cols(submats(1)%nc))
      submats(1)%cols = nrange(1:submats(1)%nc)
      allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
      call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
      fullmat = submats(1)%dat
      deallocate(submats(1)%rows)
      deallocate(submats(1)%cols)
      deallocate(submats(1)%dat)


      ! call element_Zmn_block_user(msh%Nunk, N_unk_loc, mrange, nrange, fullmat, msh, option, ker, 0, passflag, ptree, stats)



      call MPI_ALLREDUCE(N_unk_loc, N_unk_locmax, 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      allocate (fullmat_tmp(msh%Nunk, N_unk_locmax))

      do pp = 1, ptree%nproc
         if (ptree%MyID == pp - 1) then
            N_unk_loc1 = N_unk_loc
            fullmat_tmp = fullmat
         endif
         call MPI_Bcast(N_unk_loc1, 1, MPI_integer, pp - 1, ptree%Comm, ierr)
         call MPI_Bcast(fullmat_tmp, N_unk_loc1*msh%Nunk, MPI_DT, pp - 1, ptree%Comm, ierr)

         if (ptree%MyID == Main_ID) then
            do jj = 1, N_unk_loc1
            do ii = 1, msh%Nunk
               write (665, *) dble(cmplx(fullmat_tmp(ii, jj), kind=8)), aimag(cmplx(fullmat_tmp(ii, jj), kind=8)), ''
            enddo
            enddo
         endif
      enddo

      deallocate (fullmat)
      deallocate (mrange)
      deallocate (nrange)
      deallocate (fullmat_tmp)

   end subroutine FULLMAT_Element

!**** Computation of the construction phase with matrix entry evaluation
   subroutine BPACK_construction_Element(bmat, option, stats, msh, ker, ptree)

      implicit none

      type(Hoption)::option
      type(Hstat)::stats
      type(Bmatrix)::bmat
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer Maxlevel
      if (allocated(msh%xyz)) deallocate (msh%xyz)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Matrix construction......"

      select case (option%format)
      case (HODLR)
         call HODLR_construction(bmat%ho_bf, option, stats, msh, ker, ptree)
      case (HMAT)
         call Hmat_construction(bmat%h_mat, option, stats, msh, ker, ptree)
      case (HSS)
         call HSS_construction(bmat%hss_bf, option, stats, msh, ker, ptree)
      end select

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Matrix construction finished"

      select case (option%format)
      case (HODLR)
         Maxlevel = bmat%ho_bf%Maxlevel
      case (HMAT)
         Maxlevel = bmat%h_mat%Maxlevel
      case (HSS)
         Maxlevel = bmat%hss_bf%Maxlevel
      end select
      if (option%lnoBP > Maxlevel .and. option%verbosity >= 0) call BPACK_CheckError(bmat, option, msh, ker, stats, ptree)

   end subroutine BPACK_construction_Element

   subroutine Hmat_construction(h_mat, option, stats, msh, ker, ptree)

      implicit none

      type(Hoption)::option
      type(Hstat)::stats
      type(Hmat)::h_mat
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer ierr
      integer i, j, ii, jj, iii, jjj, k, kk
      integer num_blocks, level
      real*8 T0, T1, T3, T4, rtemp1, rtemp2
      real*8 rtemp
      type(matrixblock), pointer :: blocks, blocks_copy
      integer Maxtmp
      integer::mrange(1), nrange(1)
      DT::mat(1, 1)
      DT:: ctemp
      real(kind=8):: scale_factor
      integer::passflag = 0
      integer::mrange_dummy(1), nrange_dummy(1)
      DT::mat_dummy(1, 1)
      type(nod), pointer::cur
      class(*), pointer::ptr
      type(intersect)::submats(1)

      call MPI_barrier(ptree%Comm, ierr)

      T0 = OMP_get_wtime()

      allocate (stats%rankmax_of_level(0:h_mat%Maxlevel))
      stats%rankmax_of_level = 0
      allocate (stats%rankmax_of_level_global(0:h_mat%Maxlevel))
      stats%rankmax_of_level_global = 0

      num_blocks = 2**msh%Dist_level
      if (ptree%MyID == ptree%nproc - 1 .and. option%verbosity >= 0) then
         write (*, *) "   "
         write (*, *) num_blocks*Rows_per_processor, 'total blocks'
      endif

      scale_factor = 0
      ! compute the largest diagonal entry as the scaling factor
      do i = 1, Rows_per_processor
         blocks => h_mat%Local_blocks(ptree%MyID + 1, i)
         do ii = blocks%headm, blocks%headm + blocks%M - 1
            mrange(1) = ii
            nrange(1) = ii


            submats(1)%nr = 1
            submats(1)%nc = 1
            allocate(submats(1)%rows(submats(1)%nr))
            submats(1)%rows = mrange(1:submats(1)%nr)
            allocate(submats(1)%cols(submats(1)%nc))
            submats(1)%cols = nrange(1:submats(1)%nc)
            allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
            call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
            mat = submats(1)%dat
            deallocate(submats(1)%rows)
            deallocate(submats(1)%cols)
            deallocate(submats(1)%dat)


            ! call element_Zmn_block_user(1, 1, mrange, nrange, mat, msh, option, ker, 0, passflag, ptree, stats)

            scale_factor = max(scale_factor, abs(mat(1, 1)/option%scale_factor))
            ! write(*,*)ii,abs(ctemp)
         enddo
      enddo
      if (scale_factor < SafeUnderflow) scale_factor = 1d0

      passflag = 0
      do while (passflag == 0)
         call element_Zmn_blocklist_user(submats, 0, msh, option, ker, 2, passflag, ptree, stats)
         ! call element_Zmn_block_user(0, 0, mrange, nrange, mat, msh, option, ker, 2, passflag, ptree, stats)
      enddo

      option%scale_factor = 1d0/scale_factor
      call MPI_ALLREDUCE(MPI_IN_PLACE, option%scale_factor, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ptree%Comm, ierr)

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
         write (*, *) 'element_Zmn is scaled by a factor of:', option%scale_factor
      endif

      allocate (h_mat%lstblks(0:h_mat%Maxlevel))
      do level = 0, h_mat%Maxlevel
         h_mat%lstblks(level) = list()
      enddo
      do i = 1, Rows_per_processor
         do j = 1, num_blocks
            T3 = OMP_get_wtime()
            blocks => h_mat%Local_blocks(j, i)
            call Hmat_GetBlkLst(blocks, option, stats, msh, ptree, h_mat)
         enddo
      enddo
      do level = 0, h_mat%Maxlevel
         call MergeSort(h_mat%lstblks(level)%head, node_score_block_ptr_row)
      enddo

      do level = 0, h_mat%Maxlevel
         ! write(*,*)h_mat%lstblks%num_nods,'niam'
         T3 = OMP_get_wtime()
         cur => h_mat%lstblks(level)%head
         rtemp1 = 0.; rtemp2 = 0.
         do ii = 1, h_mat%lstblks(level)%num_nods
            select type (ptr=>cur%item)
            type is (block_ptr)
               call Hmat_block_construction(ptr%ptr, rtemp1, rtemp2, option, stats, msh, ker, ptree)
            end select
            cur => cur%next
         enddo
         stats%Mem_Comp_for = stats%Mem_Comp_for + rtemp1
         stats%Mem_Direct_for = stats%Mem_Direct_for + rtemp2
         stats%Mem_Peak = stats%Mem_Peak + rtemp1 + rtemp2 + rtemp1 + rtemp2

         T4 = OMP_get_wtime()
         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
            write (*, '(A10,I6,A10,Es14.7,A8,Es14.7,A8,Es14.7)') 'Level:', level, 'finished', T4 - T3, 'secnds', rtemp1 + rtemp2, 'Mbytes'
         endif

         passflag = 0
         do while (passflag < 2)
            call element_Zmn_blocklist_user(submats, 0, msh, option, ker, 2, passflag, ptree, stats)
            ! call element_Zmn_block_user(0, 0, mrange_dummy, nrange_dummy, mat_dummy, msh, option, ker, 2, passflag, ptree, stats)
         enddo
      enddo

      do i = 1, Rows_per_processor
         do j = 1, num_blocks
            blocks => h_mat%Local_blocks(j, i)
            blocks_copy => h_mat%Local_blocks_copy(j, i)
            call Hmat_block_copy('N', blocks_copy, blocks)
         enddo
      enddo

      T1 = OMP_get_wtime()
      stats%Time_Fill = stats%Time_Fill + T1 - T0

      call MPI_ALLREDUCE(stats%rankmax_of_level(0:h_mat%Maxlevel), stats%rankmax_of_level_global(0:h_mat%Maxlevel), h_mat%Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)

      stats%Mem_Fill = stats%Mem_Comp_for + stats%Mem_Direct_for

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'rankmax_of_level:', stats%rankmax_of_level_global
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      call MPI_ALLREDUCE(stats%Time_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Total construction time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Time_Entry, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Total entry eval time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Flop_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A26Es14.2)') 'Total construction flops:', rtemp

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      call MPI_ALLREDUCE(stats%Mem_Comp_for, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) rtemp, 'MB costed for butterfly forward blocks'
      call MPI_ALLREDUCE(stats%Mem_Direct_for, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) rtemp, 'MB costed for direct forward blocks'
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      ! stop
      return

   end subroutine Hmat_construction

   subroutine HODLR_construction(ho_bf1, option, stats, msh, ker, ptree)

      use BPACK_DEFS
      implicit none
      real(kind=8) n1, n2, n3, n4, n5
      integer i, j, ii, ii_inv, jj, kk, iii, jjj, ll
      integer level, blocks, edge, patch, node, group
      integer rank, index_near, m, n, length, flag, itemp, rank0_inner, rank0_outter, ierr
      real T0
      real(kind=8):: rtemp, rel_error, error, t1, t2, tim_tmp, rankrate_inner, rankrate_outter
      integer mm, nn, header_m, header_n, edge_m, edge_n, group_m, group_n, group_m1, group_n1, group_m2, group_n2
      type(matrixblock)::block_tmp, block_tmp1
      DT, allocatable::fullmat(:, :)
      integer level_c, iter, level_cc, level_butterfly, Bidxs, Bidxe
      type(Hoption)::option
      type(Hstat)::stats
      type(hobf)::ho_bf1
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer::passflag = 0
      integer::mrange_dummy(1), nrange_dummy(1)
      DT::mat_dummy(1, 1)
      type(intersect)::submats(1)

      ! Memory_direct_forward=0
      ! Memory_butterfly_forward=0
      tim_tmp = 0
      !tolerance=0.001
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      ! write (*,*) 'ACA error threshold',tolerance
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'SVD error threshold', option%tol_comp
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "constructing Leaf-blocks......"

      n1 = OMP_get_wtime()
      level = 0
      flag = 0
      ! ForwardSymmetricFlag = 0
      allocate (stats%rankmax_of_level(0:ho_bf1%Maxlevel))
      stats%rankmax_of_level = 0
      allocate (stats%rankmax_of_level_global(0:ho_bf1%Maxlevel))
      stats%rankmax_of_level_global = 0

      do level_c = 1, ho_bf1%Maxlevel + 1
         ! do level_c = ho_bf1%Maxlevel+1,ho_bf1%Maxlevel+1
         if (level_c /= ho_bf1%Maxlevel + 1) then
            Bidxs = ho_bf1%levels(level_c)%Bidxs*2 - 1
            Bidxe = ho_bf1%levels(level_c)%Bidxe*2
         else
            Bidxs = ho_bf1%levels(level_c)%Bidxs
            Bidxe = ho_bf1%levels(level_c)%Bidxe
         endif
         n3 = OMP_get_wtime()
         do ii = Bidxs, Bidxe
            ! do ii =Bidxs,Bidxs
            if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(ii)%pgno)) then
               if (level_c /= ho_bf1%Maxlevel + 1) then
                  if (ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level /= level) then
                     level = ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level
                     if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'constructing level', level
                  endif

                  ! if(mod(ii,2)==1)then
                  call Bplus_compress_entry(ho_bf1%levels(level_c)%BP(ii), option, rtemp, stats, msh, ker, ptree)
                  ! else
                  ! call BF_delete(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1),1)
                  ! call BF_copy('T',ho_bf1%levels(level_c)%BP(ii-1)%LL(1)%matrices_block(1),ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1))
                  ! endif

                  if (level == option%level_check) then
                     ! call Bplus_randomized_Exact_test(ho_bf1%levels(level_c)%BP(ii))

                     rank0_inner = ho_bf1%levels(level_c)%BP(ii)%LL(2)%rankmax
                     rankrate_inner = 1.2d0
                     rank0_outter = ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%rankmax
                     rankrate_outter = 1.2d0
                     level_butterfly = ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level_butterfly

                     t1 = OMP_GET_WTIME()
                     call BF_randomized(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%pgno, level_butterfly, rank0_outter, rankrate_outter, ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1), ho_bf1%levels(level_c)%BP(ii), Bplus_block_MVP_Exact_dat, error, 'Exact', option, stats, ptree, msh)

                     call BF_ComputeMemory(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1), rtemp)

                     t2 = OMP_GET_WTIME()
                     tim_tmp = tim_tmp + t2 - t1

                     ! call Bplus_randomized_constr(level_butterfly,ho_bf1%levels(level_c)%BP(ii),ho_bf1%levels(level_c)%BP(ii),rank0_inner,rankrate_inner,Bplus_block_MVP_Exact_dat,rank0_outter,rankrate_outter,Bplus_block_MVP_Outter_Exact_dat,error,'Exact',option,stats,ptree,msh)

                     if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
                        write (*, *) 'time_tmp', time_tmp, 'randomized_bf time,', tim_tmp, 'stats%Time_random,', stats%Time_random, 'mem', rtemp
                        stop
                     endif
                  end if

                  stats%Mem_Comp_for = stats%Mem_Comp_for + rtemp
               else

                  if (ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level /= level) then
                     level = ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level
                     if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'constructing level', level
                  endif
                  call Full_construction(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1), msh, ker, stats, option, ptree)
                  stats%Mem_Direct_for = stats%Mem_Direct_for + SIZEOF(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%fullmat)/1024.0d3
               endif
               ! ! write(*,*)level_c,ii,ho_bf1%levels(level_c)%N_block_forward
               ! if(level>=option%level_check .and. level_c/=ho_bf1%Maxlevel+1)then
               ! call BF_compress_test(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1),msh,ker,element_Zmn,ptree,option,stats)
               ! !stop
               ! end if
            endif
         end do

         ! call MPI_barrier(ptree%Comm, ierr)
         ! do while(option%elem_extract == 1 .and. level==2)

         ! enddo
         passflag = 0
         do while (passflag < 2)
            ! call element_Zmn_block_user(0, 0, mrange_dummy, nrange_dummy, mat_dummy, msh, option, ker, 2, passflag, ptree, stats)
            ! write(*,*)ptree%MyID,level, 'level done'
            call element_Zmn_blocklist_user(submats, 0, msh, option, ker, 2, passflag, ptree, stats)
         enddo

         n4 = OMP_get_wtime()
         n5 = n4 - n3

         call MPI_ALLREDUCE(MPI_IN_PLACE, n5, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time', n5, 'rankmax_of_level so far:', stats%rankmax_of_level
      end do
      n2 = OMP_get_wtime()
      stats%Time_Fill = stats%Time_Fill + n2 - n1

      stats%Mem_Fill = stats%Mem_Comp_for + stats%Mem_Direct_for
      stats%Mem_Peak = stats%Mem_Peak + stats%Mem_Fill

      call MPI_ALLREDUCE(stats%rankmax_of_level(0:ho_bf1%Maxlevel), stats%rankmax_of_level_global(0:ho_bf1%Maxlevel), ho_bf1%Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'rankmax_of_level:', stats%rankmax_of_level_global
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      call MPI_ALLREDUCE(stats%Time_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Total construction time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Time_Entry, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Total entry eval time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Flop_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A26Es14.2)') 'Total construction flops:', rtemp
      call MPI_ALLREDUCE(time_tmp, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp', time_tmp

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      call MPI_ALLREDUCE(stats%Mem_Comp_for, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) rtemp, 'MB costed for butterfly forward blocks'
      call MPI_ALLREDUCE(stats%Mem_Direct_for, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) rtemp, 'MB costed for direct forward blocks'
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      ! stop

      return

   end subroutine HODLR_construction

   subroutine HSS_construction(hss_bf1, option, stats, msh, ker, ptree)

      use BPACK_DEFS
      implicit none
      real(kind=8) n1, n2, n3, n4, n5
      integer i, j, ii, ii_inv, jj, kk, iii, jjj, ll
      integer level, blocks, edge, patch, node, group
      integer rank, index_near, m, n, length, flag, itemp, rank0_inner, rank0_outter, ierr
      real T0
      real(kind=8):: rtemp, rel_error, error, t1, t2, tim_tmp, rankrate_inner, rankrate_outter
      integer mm, nn, header_m, header_n, edge_m, edge_n, group_m, group_n, group_m1, group_n1, group_m2, group_n2
      type(matrixblock)::block_tmp, block_tmp1
      DT, allocatable::fullmat(:, :)
      integer level_c, iter, level_cc, level_butterfly, Bidxs, Bidxe
      type(Hoption)::option
      type(Hstat)::stats
      type(hssbf)::hss_bf1
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer::passflag = 0
      type(intersect)::submats(1)
      integer::mrange_dummy(1), nrange_dummy(1)
      DT::mat_dummy(1, 1)

      ! Memory_direct_forward=0
      ! Memory_butterfly_forward=0
      tim_tmp = 0
      !tolerance=0.001
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      ! write (*,*) 'ACA error threshold',tolerance
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'SVD error threshold', option%tol_comp
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "constructing HSS-BF......"

      n1 = OMP_get_wtime()

      allocate (stats%rankmax_of_level(0:hss_bf1%Maxlevel))
      stats%rankmax_of_level = 0
      allocate (stats%rankmax_of_level_global(0:hss_bf1%Maxlevel))
      stats%rankmax_of_level_global = 0

      call Bplus_compress_entry(hss_bf1%BP, option, rtemp, stats, msh, ker, ptree)
      stats%Mem_Comp_for = stats%Mem_Comp_for + rtemp

      passflag = 0
      do while (passflag < 2)
         ! call element_Zmn_block_user(0, 0, mrange_dummy, nrange_dummy, mat_dummy, msh, option, ker, 2, passflag, ptree, stats)
         call element_Zmn_blocklist_user(submats, 0, msh, option, ker, 2, passflag, ptree, stats)
      enddo

      do ll = 1, hss_bf1%BP%Lplus
         call MPI_ALLREDUCE(MPI_IN_PLACE, hss_bf1%BP%LL(ll)%rankmax, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(hss_bf1%BP%LL(1)%matrices_block(1)%pgno)%Comm, ierr)
      enddo


      n2 = OMP_get_wtime()
      stats%Time_Fill = stats%Time_Fill + n2 - n1

      stats%Mem_Fill = stats%Mem_Comp_for + stats%Mem_Direct_for
      stats%Mem_Peak = stats%Mem_Peak + stats%Mem_Fill

      call MPI_ALLREDUCE(stats%rankmax_of_level(0:hss_bf1%Maxlevel), stats%rankmax_of_level_global(0:hss_bf1%Maxlevel), hss_bf1%Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'rankmax_of_level:', stats%rankmax_of_level_global
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      call MPI_ALLREDUCE(stats%Time_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Total construction time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Time_Entry, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Total entry eval time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Flop_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A26Es14.2)') 'Total construction flops:', rtemp
      call MPI_ALLREDUCE(time_tmp, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp', time_tmp

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      call MPI_ALLREDUCE(stats%Mem_Comp_for, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) rtemp, 'MB costed for butterfly forward blocks'
      call MPI_ALLREDUCE(stats%Mem_Direct_for, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) rtemp, 'MB costed for direct forward blocks'
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      ! stop

      return

   end subroutine HSS_construction

   recursive subroutine Hmat_block_construction(blocks, Memory_far, Memory_near, option, stats, msh, ker, ptree)

      implicit none

      type(matrixblock), pointer :: blocks_son
      type(matrixblock) :: blocks
      real*8 Memory_far, Memory_near, rtemp, Memory_tmp, t1, t2
      integer i, j, k, flag, conv, m, n, ii, jj

      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer, allocatable:: boundary_map(:)
      integer level_butterfly, levelm, groupm_start, Nboundall

      ! t1=OMP_GET_WTIME()
      if (blocks%style == 2) then

         groupm_start = 0
         Nboundall = 0

         if (option%forwardN15flag == 1) then
            call BF_compress_N15(blocks, boundary_map, Nboundall, groupm_start, option, Memory_tmp, stats, msh, ker, ptree, 1)
            call BF_sym2asym(blocks)
         else
            call BF_compress_NlogN(blocks, boundary_map, Nboundall, groupm_start, option, Memory_tmp, stats, msh, ker, ptree, 1)
         end if
         Memory_far = Memory_far + Memory_tmp

      elseif (blocks%style == 1) then
         call Full_construction(blocks, msh, ker, stats, option, ptree)
         Memory_near = Memory_near + SIZEOF(blocks%fullmat)/1024.0d3
      elseif (blocks%style == 4) then
         do ii = 1, 2
         do jj = 1, 2
            blocks_son => blocks%sons(ii, jj)
            call Hmat_block_construction(blocks_son, Memory_far, Memory_near, option, stats, msh, ker, ptree)
         enddo
         enddo
      endif
! t2=OMP_GET_WTIME()
! if(blocks%level==1)write(*,*)blocks%row_group,blocks%col_group,ptree%MyID,t2-t1
      return

   end subroutine Hmat_block_construction

!!!!!!! check error of BF construction using parallel element extraction
   subroutine BF_checkError(blocks, option, msh, ker, stats, ptree)
      use BPACK_DEFS
      implicit none

      type(Hoption)::option
      type(Hstat)::stats
      ! type(Bmatrix)::bmat
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(intersect), allocatable::inters(:)
      real(kind=8)::n1, n2, n3, n4
      integer Ntest, passflag
      integer nsproc1, nsproc2, nprow, npcol, nprow1D, npcol1D, myrow, mycol, nprow1, npcol1, myrow1, mycol1, nprow2, npcol2, myrow2, mycol2, myArows, myAcols, rank1, rank2, ierr, MyID
      integer:: cridx, info
      integer, allocatable::rows(:), cols(:)
      DT, allocatable:: Vin(:, :), Vout1(:, :), Vout2(:, :), Vinter(:, :), Mat(:, :)
      integer::descVin(9), descVout(9), descVinter(9), descFull(9), descButterflyU(9), descButterflyV(9)
      integer N, M, i, j, ii, jj, nn, pp, myi, myj, iproc, jproc, rmax
      integer edge_n, edge_m, rank
      real(kind=8):: fnorm1, fnorm0, rtemp1 = 0, rtemp0 = 0
      real(kind=8):: a, v1, v2, v3
      DT:: value1, value2, value3
      type(list)::lstr, lstc, lst, lstblk
      type(nod), pointer::cur, curr, curc, curri, curci
      class(*), pointer::ptr, ptrr, ptrc, ptrri, ptrci
      integer::head, tail, idx, pgno, ctxt, nr_loc, nc_loc
      type(matrixblock), pointer::blocks
      integer num_blocks, idx_row, idx_col, idx_dat
      integer, allocatable:: allrows(:), allcols(:), pmaps(:, :)
      integer, allocatable::datidx(:), colidx(:), rowidx(:), pgidx(:)
      DT, allocatable::alldat_loc(:)
      integer:: Ninter, nr, nc, ntot_loc, level, Npmap, nproc, npavr, np
      type(intersect)::submats(1)

      ! Ninter=2**level
      ! nr=2500
      ! nc=2500

      Ninter = 4

      nr = 100
      nc = 100

      allocate (colidx(Ninter))
      allocate (rowidx(Ninter))
      allocate (pgidx(Ninter))
      ! allocate(datidx(Ninter))

      allocate (allrows(Ninter*nr))
      allocate (allcols(Ninter*nc))

      ! pgno=1
      ! ctxt = ptree%pgrp(pgno)%ctxt
      ! call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
      ! nprow = ptree%pgrp(pgno)%nprow
      ! npcol = ptree%pgrp(pgno)%npcol

      nproc = ptree%nproc
      Npmap = min(Ninter, nproc)
      npavr = nproc/Npmap
      allocate (pmaps(Npmap, 3))
      do nn = 1, Npmap
         nprow = floor_safe(sqrt(dble(npavr)))
         npcol = floor_safe(npavr/dble(nprow))
         pmaps(nn, 1) = nprow
         pmaps(nn, 2) = npcol
         pmaps(nn, 3) = (nn - 1)*npavr
      enddo

      idx_row = 0
      idx_col = 0
      idx_dat = 0
      ! ntot_loc=0
      pp = 0
      do nn = 1, Ninter
         rowidx(nn) = nr
         colidx(nn) = nc
         pp = pp + 1
         pp = mod(pp - 1, Npmap) + 1
         pgidx(nn) = pp
         nprow = pmaps(pgidx(nn), 1)
         npcol = pmaps(pgidx(nn), 2)
         call Gridinfo_2D(pmaps(pgidx(nn), :), ptree%MyID, myrow, mycol)
         ! datidx(nn)=ntot_loc
         if (myrow /= -1 .and. mycol /= -1) then
            myArows = numroc_wp(nr, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(nc, nbslpk, mycol, 0, npcol)
            idx_dat = idx_dat + myArows*myAcols
         endif

         do ii = 1, nr
            call random_number(a)
            call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, Main_ID, ptree%Comm, ierr)

            allrows(idx_row + 1) = max(floor_safe(blocks%M*a), 1)
            ! allrows(idx_row+1)=msh%basis_group(2**level+nn-1)%head+ii-1
            idx_row = idx_row + 1
         enddo

         do ii = 1, nc
            call random_number(a)
            call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, Main_ID, ptree%Comm, ierr)
            allcols(idx_col + 1) = max(floor_safe(blocks%N*a), 1) + blocks%M
            ! allcols(idx_col+1)=msh%basis_group(2**level+1-(nn-1))%head+ii-1
            idx_col = idx_col + 1
         enddo
      enddo

      allocate (alldat_loc(idx_dat))
      if (idx_dat > 0) alldat_loc = 0

      n1 = OMP_get_wtime()
      call BF_ExtractElement(blocks, option, msh, stats, ptree, Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps)
      n2 = OMP_get_wtime()

      ! compare extracted values with element_Zmn
      v1 = 0
      v2 = 0
      v3 = 0
      idx_row = 0
      idx_col = 0
      idx_dat = 0
      do nn = 1, Ninter
         nr = rowidx(nn)
         nc = colidx(nn)
         nprow = pmaps(pgidx(nn), 1)
         npcol = pmaps(pgidx(nn), 2)
         call Gridinfo_2D(pmaps(pgidx(nn), :), ptree%MyID, myrow, mycol)
         if (myrow /= -1 .and. mycol /= -1) then
            nr_loc = numroc_wp(nr, nbslpk, myrow, 0, nprow)
            nc_loc = numroc_wp(nc, nbslpk, mycol, 0, npcol)
            allocate (rows(nr_loc))
            allocate (cols(nc_loc))
            allocate (Mat(nr_loc, nc_loc))
            do myi = 1, nr_loc
               call l2g(myi, myrow, nr, nprow, nbslpk, ii)
               rows(myi) = allrows(ii + idx_row)
            enddo
            do myj = 1, nc_loc
               call l2g(myj, mycol, nc, npcol, nbslpk, jj)
               cols(myj) = allcols(jj + idx_col)
            enddo


            submats(1)%nr = nr_loc
            submats(1)%nc = nc_loc
            allocate(submats(1)%rows(submats(1)%nr))
            submats(1)%rows = rows(1:submats(1)%nr)
            allocate(submats(1)%cols(submats(1)%nc))
            submats(1)%cols = cols(1:submats(1)%nc)
            allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
            call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
            Mat = submats(1)%dat
            deallocate(submats(1)%rows)
            deallocate(submats(1)%cols)
            deallocate(submats(1)%dat)

            ! call element_Zmn_block_user(nr_loc, nc_loc, rows, cols, Mat, msh, option, ker, 0, passflag, ptree, stats)

            do myi = 1, nr_loc
            do myj = 1, nc_loc
               value2 = alldat_loc(idx_dat + myi + (myj - 1)*nr_loc)
               value1 = Mat(myi, myj)
               v1 = v1 + abs(value1)**2d0
               v2 = v2 + abs(value2)**2d0
               v3 = v3 + abs(value2 - value1)**2d0
               ! if(abs(value2-value1)**2d0/abs(value1)**2d0>1D-3)write(*,*)ptree%MyID,nn,myi,myj,abs(value2-value1)**2d0/abs(value1)**2d0,value2,value1
            enddo
            enddo
            idx_dat = idx_dat + nr_loc*nc_loc
         else
            nr_loc = 0
            nc_loc = 0
            allocate (rows(nr_loc))
            allocate (cols(nc_loc))
            allocate (Mat(nr_loc, nc_loc))
            call element_Zmn_blocklist_user(submats, 0, msh, option, ker, 2, passflag, ptree, stats)
            ! call element_Zmn_block_user(nr_loc, nc_loc, rows, cols, Mat, msh, option, ker, 2, passflag, ptree, stats)
         endif
         deallocate (rows)
         deallocate (cols)
         deallocate (Mat)
         idx_row = idx_row + nr
         idx_col = idx_col + nc
      enddo

      deallocate (rowidx)
      deallocate (colidx)
      deallocate (pgidx)
      ! deallocate(datidx)
      deallocate (allrows)
      deallocate (allcols)
      deallocate (alldat_loc)
      deallocate (pmaps)

      call MPI_ALLREDUCE(MPI_IN_PLACE, v1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, v2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, v3, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A25,Es14.7,Es14.7,A6,Es9.2,A7,Es9.2)') 'BF_CheckError: fnorm:', sqrt(v1), sqrt(v2), ' acc: ', sqrt(v3/v1), ' time: ', n2 - n1

      !stop

   end subroutine BF_checkError

!!!!!!! extract a list of intersections from a block
   subroutine BF_ExtractElement(blocks_o, option, msh, stats, ptree, Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps)
      use BPACK_DEFS
      implicit none

      type(Hoption)::option
      type(Hstat)::stats
      type(Bmatrix)::bmat
      type(mesh)::msh
      type(proctree)::ptree
      type(intersect), allocatable::inters(:)
      real(kind=8)::n0, n1, n2, n3, n4, n5
      integer Ntest, passflag
      integer nsproc1, nsproc2, nprow, npcol, nprow1D, npcol1D, myrow, mycol, nprow1, npcol1, myrow1, mycol1, nprow2, npcol2, myrow2, mycol2, myArows, myAcols, rank1, rank2, ierr, MyID
      integer:: cridx, info
      integer, allocatable::rows(:), cols(:)
      DT, allocatable:: Vin(:, :), Vout1(:, :), Vout2(:, :), Vinter(:, :), Mat(:, :)
      integer::descVin(9), descVout(9), descVinter(9), descFull(9), descButterflyU(9), descButterflyV(9)
      integer N, M, i, j, ii, jj, nn, myi, myj, iproc, jproc, rmax
      integer edge_n, edge_m, rank
      real(kind=8):: fnorm1, fnorm0, rtemp1 = 0, rtemp0 = 0
      real(kind=8):: a, v1, v2, v3
      DT:: value1, value2, value3
      type(list)::lstr, lstc, lstblk
      type(nod), pointer::cur, curr, curc, curri, curci
      class(*), pointer::ptr, ptrr, ptrc, ptrri, ptrci
      integer::head, tail, idx, pp, pgno, nr_loc, nc_loc
      type(matrixblock), pointer::blocks
      integer num_blocks, idx_row, idx_col, idx_dat
      integer:: allrows(:), allcols(:)
      integer, allocatable::datidx(:)
      integer:: Ninter, nr, nc, ntot_loc
      DT::alldat_loc(:)
      integer::colidx(Ninter), rowidx(Ninter), pgidx(Ninter)
      integer::Npmap, pmaps(Npmap, 3), flag2D
      type(iarray)::lst
      type(matrixblock), pointer::blocks_o

      stats%Flop_Tmp = 0d0

      n0 = OMP_get_wtime()
      flag2D = 0
      allocate (inters(Ninter))
      lstr = list()
      lstc = list()
      ! lst=list()
      lstblk = list()
      idx_row = 0
      idx_col = 0
      do nn = 1, Ninter
         nr = rowidx(nn)
         nc = colidx(nn)
         inters(nn)%nr = nr
         inters(nn)%nc = nc
         inters(nn)%pg = pgidx(nn)
         nprow = pmaps(pgidx(nn), 1)
         npcol = pmaps(pgidx(nn), 2)
         call Gridinfo_2D(pmaps(pgidx(nn), :), ptree%MyID, myrow, mycol)
         if (myrow /= -1 .and. mycol /= -1) then
            myArows = numroc_wp(nr, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(nc, nbslpk, mycol, 0, npcol)
            allocate (inters(nn)%dat_loc(myArows, myAcols))
            if (myArows > 0 .and. myAcols > 0) inters(nn)%dat_loc = 0
         endif
         if (nprow*npcol > 1) flag2D = 1
         allocate (inters(nn)%rows(nr))
         lst%num_nods = nr
         lst%idx = nn
         allocate (lst%dat(lst%num_nods))
         do ii = 1, nr
            idx_row = idx_row + 1
            inters(nn)%rows(ii) = allrows(idx_row)
            lst%dat(ii) = ii
         enddo
         call append(lstr, lst)
         call iarray_finalizer(lst)

         allocate (inters(nn)%cols(nc))
         lst%num_nods = nc
         lst%idx = nn
         allocate (lst%dat(lst%num_nods))
         do ii = 1, nc
            idx_col = idx_col + 1
            inters(nn)%cols(ii) = allcols(idx_col)
            lst%dat(ii) = ii
         enddo
         call append(lstc, lst)
         call iarray_finalizer(lst)
      enddo

      n1 = OMP_get_wtime()

      curr => lstr%head
      curc => lstc%head
      do nn = 1, Ninter
         ptrr=>curr%item
         select type (ptrr)
         type is (iarray)
            ptrc=>curc%item
            select type (ptrc)
            type is (iarray)
               call Hmat_MapIntersec2Block_Loc(blocks_o, option, stats, msh, ptree, inters, nn, ptrr, ptrc, lstblk)
            end select
         end select
         curr => curr%next
         curc => curc%next
      enddo

      call MergeSort(lstblk%head, node_score_block_ptr_row)

      n2 = OMP_get_wtime()
      stats%Time_Entry_Traverse = stats%Time_Entry_Traverse + n2-n1

#if 0
      write (*, *) lstblk%num_nods
      cur => lstblk%head
      do ii = 1, lstblk%num_nods
         select type (ptr=>cur%item)
         type is (block_ptr)
            curr => ptr%ptr%lstr%head
            curc => ptr%ptr%lstc%head
            do nn = 1, ptr%ptr%lstr%num_nods
               ptrr=>curr%item
               select type (ptrr)
               type is (iarray)
                  ptrc=>curc%item
                  select type (ptrc)
                  type is (iarray)
                     write (*, *) ptree%MyID, ptr%ptr%row_group, ptr%ptr%col_group, nn, ptrr%num_nods, ptrc%num_nods
                  end select
               end select
               curr => curr%next
               curc => curc%next
            enddo
         end select
         cur => cur%next
      enddo
#endif

#ifdef HAVE_TASKLOOP
      !$omp parallel
      !$omp single
#endif
      ! construct intersections for each block from the block's lists
      cur => lstblk%head
      do ii = 1, lstblk%num_nods ! loop all blocks
         select type (ptr=>cur%item)
         type is (block_ptr)
            blocks => ptr%ptr
            pp = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1
            head = blocks%M_p(pp, 1) + blocks%headm - 1
            tail = blocks%M_p(pp, 2) + blocks%headm - 1

            curr => blocks%lstr%head
            curc => blocks%lstc%head
            allocate (blocks%inters(blocks%lstr%num_nods))
            do nn = 1, blocks%lstr%num_nods ! loop all lists of list of rows and columns
               ptrr=>curr%item
               select type (ptrr)
               type is (iarray)
                  ptrc=>curc%item
                  select type (ptrc)
                  type is (iarray)
                     blocks%inters(nn)%nr = ptrr%num_nods
                     allocate (blocks%inters(nn)%rows(ptrr%num_nods))
                     blocks%inters(nn)%rows = ptrr%dat
                     blocks%inters(nn)%nc = ptrc%num_nods
                     allocate (blocks%inters(nn)%cols(ptrc%num_nods))
                     blocks%inters(nn)%cols = ptrc%dat
                     blocks%inters(nn)%idx = ptrr%idx
                  end select
               end select
               curr => curr%next
               curc => curc%next

               blocks%inters(nn)%nr_loc = 0
               do jj = 1, blocks%inters(nn)%nr
                  idx = blocks%inters(nn)%rows(jj)
                  if (inters(blocks%inters(nn)%idx)%rows(idx) >= head .and. inters(blocks%inters(nn)%idx)%rows(idx) <= tail) then
                     blocks%inters(nn)%nr_loc = blocks%inters(nn)%nr_loc + 1
                  endif
               enddo
               allocate (blocks%inters(nn)%rows_loc(blocks%inters(nn)%nr_loc))
               allocate (blocks%inters(nn)%glo2loc(blocks%inters(nn)%nr))
               blocks%inters(nn)%glo2loc = -1
               blocks%inters(nn)%nr_loc = 0
               do jj = 1, blocks%inters(nn)%nr
                  idx = blocks%inters(nn)%rows(jj)
                  ! write(*,*)inters(blocks%inters(nn)%idx)%rows(idx),head,tail
                  if (inters(blocks%inters(nn)%idx)%rows(idx) >= head .and. inters(blocks%inters(nn)%idx)%rows(idx) <= tail) then
                     blocks%inters(nn)%nr_loc = blocks%inters(nn)%nr_loc + 1
                     blocks%inters(nn)%rows_loc(blocks%inters(nn)%nr_loc) = jj ! rows_loc stores indices in rows
                     blocks%inters(nn)%glo2loc(jj) = blocks%inters(nn)%nr_loc
                  endif
               enddo
               allocate (blocks%inters(nn)%dat_loc(blocks%inters(nn)%nr_loc, blocks%inters(nn)%nc))
            enddo

            ! extract entries on an array of intersections for each block

            if (blocks%style == 1) then
               call Full_block_extraction(blocks, inters, ptree, msh, stats)
            else
               if (blocks%level_butterfly == 0) then
                  call LR_block_extraction(blocks, inters, ptree, msh, stats)
               else
                  call BF_block_extraction(blocks, inters, ptree, msh, stats)
               endif
            endif

            ! finalize the lists of lists of rows and columns for each block because they have been transferred to intersections
            call list_finalizer(blocks%lstr)
            call list_finalizer(blocks%lstc)
         end select
         cur => cur%next
      enddo
#ifdef HAVE_TASKLOOP
      !$omp end single
      !$omp end parallel
#endif
      n3 = OMP_get_wtime()
      stats%Time_Entry_BF = stats%Time_Entry_BF + n3-n2

      call MPI_barrier(ptree%Comm, ierr)
      n3 = OMP_get_wtime()


      ! redistribute from blocks' intersections to the global intersecions inters
      ! if (flag2D == 1) then ! if each intersection is only needed by one processor, the communication can be optimized
         call BPACK_all2all_inters(inters, lstblk, stats, ptree, ptree%nproc, Npmap, pmaps)
      ! else
      !    call BPACK_all2all_inters_optimized(inters, lstblk, stats, ptree, ptree%nproc, Npmap, pmaps)
      ! endif

      n4 = OMP_get_wtime()
      stats%Time_Entry_Comm = stats%Time_Entry_Comm + n4-n3


      ntot_loc = 0
      do nn = 1, Ninter
         if (allocated(inters(nn)%dat_loc)) then
            nr_loc = size(inters(nn)%dat_loc, 1)
            nc_loc = size(inters(nn)%dat_loc, 2)
            do jj = 1, nc_loc
            do ii = 1, nr_loc
               alldat_loc(ntot_loc + ii + (jj - 1)*nr_loc) = inters(nn)%dat_loc(ii, jj)
            enddo
            enddo
            ntot_loc = ntot_loc + nr_loc*nc_loc
         endif
      enddo

      ! deallocate intersections at each block
      cur => lstblk%head
      do ii = 1, lstblk%num_nods
         select type (ptr=>cur%item)
         type is (block_ptr)
            blocks => ptr%ptr
            do nn = 1, size(blocks%inters, 1)
               if (allocated(blocks%inters(nn)%dat)) deallocate (blocks%inters(nn)%dat)
               if (allocated(blocks%inters(nn)%dat_loc)) deallocate (blocks%inters(nn)%dat_loc)
               if (allocated(blocks%inters(nn)%rows)) deallocate (blocks%inters(nn)%rows)
               if (allocated(blocks%inters(nn)%cols)) deallocate (blocks%inters(nn)%cols)
               if (allocated(blocks%inters(nn)%rows_loc)) deallocate (blocks%inters(nn)%rows_loc)
               blocks%inters(nn)%nr = 0
               blocks%inters(nn)%nr_loc = 0
               blocks%inters(nn)%nc = 0
               blocks%inters(nn)%idx = 0
            enddo
            deallocate (blocks%inters)
         end select
         cur => cur%next
      enddo

      ! finalize the list of block_ptr
      call list_finalizer(lstblk)

      ! deallocate global intersections
      do nn = 1, Ninter
         if (allocated(inters(nn)%dat)) deallocate (inters(nn)%dat)
         if (allocated(inters(nn)%dat_loc)) deallocate (inters(nn)%dat_loc)
         if (allocated(inters(nn)%rows)) deallocate (inters(nn)%rows)
         if (allocated(inters(nn)%cols)) deallocate (inters(nn)%cols)
         if (allocated(inters(nn)%rows_loc)) deallocate (inters(nn)%rows_loc)
      enddo
      deallocate (inters)

      call list_finalizer(lstr)
      call list_finalizer(lstc)

      n5 = OMP_get_wtime()

      ! if(ptree%MyID==Main_ID)then
      ! write(*,*)n1-n0,n2-n1,n3-n2,n4-n3,n5-n4
      ! endif

   end subroutine BF_ExtractElement

!!!!!!! extract a list of intersections from a bmat
   subroutine BPACK_ExtractElement(bmat, option, msh, stats, ptree, Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps)
      use BPACK_DEFS
      implicit none

      type(Hoption)::option
      type(Hstat)::stats
      type(Bmatrix)::bmat
      type(mesh)::msh
      type(proctree)::ptree
      type(intersect), allocatable::inters(:)
      real(kind=8)::n0, n1, n2, n3, n4, n5
      integer Ntest, passflag
      integer nsproc1, nsproc2, nprow, npcol, nprow1D, npcol1D, myrow, mycol, nprow1, npcol1, myrow1, mycol1, nprow2, npcol2, myrow2, mycol2, myArows, myAcols, rank1, rank2, ierr, MyID
      integer:: cridx, info
      integer, allocatable::rows(:), cols(:)
      DT, allocatable:: Vin(:, :), Vout1(:, :), Vout2(:, :), Vinter(:, :), Mat(:, :)
      integer::descVin(9), descVout(9), descVinter(9), descFull(9), descButterflyU(9), descButterflyV(9)
      integer N, M, i, j, ii, jj, nn, myi, myj, iproc, jproc, rmax
      integer edge_n, edge_m, rank
      real(kind=8):: fnorm1, fnorm0, rtemp1 = 0, rtemp0 = 0
      real(kind=8):: a, v1, v2, v3
      DT:: value1, value2, value3
      type(list)::lstr, lstc, lstblk
      type(nod), pointer::cur, curr, curc, curri, curci
      class(*), pointer::ptr, ptrr, ptrc, ptrri, ptrci
      integer::head, tail, idx, pp, pgno, nr_loc, nc_loc
      type(matrixblock), pointer::blocks
      integer num_blocks, idx_row, idx_col, idx_dat
      integer:: allrows(:), allcols(:)
      integer, allocatable::datidx(:)
      integer:: Ninter, nr, nc, ntot_loc
      DT::alldat_loc(:)
      integer::colidx(Ninter), rowidx(Ninter), pgidx(Ninter)
      integer::Npmap, pmaps(Npmap, 3), flag2D
      type(iarray)::lst

      stats%Flop_Tmp = 0d0

      n0 = OMP_get_wtime()
      flag2D = 0
      allocate (inters(Ninter))
      lstr = list()
      lstc = list()
      ! lst=list()
      lstblk = list()
      idx_row = 0
      idx_col = 0
      do nn = 1, Ninter
         nr = rowidx(nn)
         nc = colidx(nn)
         inters(nn)%nr = nr
         inters(nn)%nc = nc
         inters(nn)%pg = pgidx(nn)
         nprow = pmaps(pgidx(nn), 1)
         npcol = pmaps(pgidx(nn), 2)
         call Gridinfo_2D(pmaps(pgidx(nn), :), ptree%MyID, myrow, mycol)
         if (myrow /= -1 .and. mycol /= -1) then
            myArows = numroc_wp(nr, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(nc, nbslpk, mycol, 0, npcol)
            allocate (inters(nn)%dat_loc(myArows, myAcols))
            if (myArows > 0 .and. myAcols > 0) inters(nn)%dat_loc = 0
         endif
         if (nprow*npcol > 1) flag2D = 1
         allocate (inters(nn)%rows(nr))
         lst%num_nods = nr
         lst%idx = nn
         allocate (lst%dat(lst%num_nods))
         do ii = 1, nr
            idx_row = idx_row + 1
            inters(nn)%rows(ii) = allrows(idx_row)
            lst%dat(ii) = ii
         enddo
         call append(lstr, lst)
         call iarray_finalizer(lst)

         allocate (inters(nn)%cols(nc))
         lst%num_nods = nc
         lst%idx = nn
         allocate (lst%dat(lst%num_nods))
         do ii = 1, nc
            idx_col = idx_col + 1
            inters(nn)%cols(ii) = allcols(idx_col)
            lst%dat(ii) = ii
         enddo
         call append(lstc, lst)
         call iarray_finalizer(lst)
      enddo

      n1 = OMP_get_wtime()

      curr => lstr%head
      curc => lstc%head
      do nn = 1, Ninter
         ptrr=>curr%item
         select type (ptrr)
         type is (iarray)
            ptrc=>curc%item
            select type (ptrc)
            type is (iarray)
               select case (option%format)
               case (HODLR)
                  call HODLR_MapIntersec2Block(bmat%ho_bf, option, stats, msh, ptree, inters, nn, ptrr, ptrc, lstblk, 1, 1, 0)
               case (HMAT)
                  num_blocks = 2**msh%Dist_level
                  call Hmat_MapIntersec2Block(bmat%h_mat, option, stats, msh, ptree, inters, nn, ptrr, ptrc, lstblk, num_blocks)
               case (HSS)
                  call HSS_MapIntersec2Block(bmat%hss_bf, option, stats, msh, ptree, inters, nn, ptrr, ptrc, lstblk, 1, bmat%hss_bf%BP%LL(1)%Nbound)
               end select
            end select
         end select
         curr => curr%next
         curc => curc%next
      enddo

      call MergeSort(lstblk%head, node_score_block_ptr_row)

      n2 = OMP_get_wtime()
      stats%Time_Entry_Traverse = stats%Time_Entry_Traverse + n2-n1

#if 0
      write (*, *) lstblk%num_nods
      cur => lstblk%head
      do ii = 1, lstblk%num_nods
         select type (ptr=>cur%item)
         type is (block_ptr)
            curr => ptr%ptr%lstr%head
            curc => ptr%ptr%lstc%head
            do nn = 1, ptr%ptr%lstr%num_nods
               ptrr=>curr%item
               select type (ptrr)
               type is (iarray)
                  ptrc=>curc%item
                  select type (ptrc)
                  type is (iarray)
                     write (*, *) ptree%MyID, ptr%ptr%row_group, ptr%ptr%col_group, nn, ptrr%num_nods*ptrc%num_nods
                  end select
               end select
               curr => curr%next
               curc => curc%next
            enddo
         end select
         cur => cur%next
      enddo
#endif

#ifdef HAVE_TASKLOOP
      !$omp parallel
      !$omp single
#endif
      ! construct intersections for each block from the block's lists
      cur => lstblk%head
      do ii = 1, lstblk%num_nods ! loop all blocks
         select type (ptr=>cur%item)
         type is (block_ptr)
            blocks => ptr%ptr
            pp = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1
            head = blocks%M_p(pp, 1) + blocks%headm - 1
            tail = blocks%M_p(pp, 2) + blocks%headm - 1

            curr => blocks%lstr%head
            curc => blocks%lstc%head
            allocate (blocks%inters(blocks%lstr%num_nods))
            do nn = 1, blocks%lstr%num_nods ! loop all lists of list of rows and columns
               ptrr=>curr%item
               select type (ptrr)
               type is (iarray)
                  ptrc=>curc%item
                  select type (ptrc)
                  type is (iarray)
                     blocks%inters(nn)%nr = ptrr%num_nods
                     allocate (blocks%inters(nn)%rows(ptrr%num_nods))
                     blocks%inters(nn)%rows = ptrr%dat
                     blocks%inters(nn)%nc = ptrc%num_nods
                     allocate (blocks%inters(nn)%cols(ptrc%num_nods))
                     blocks%inters(nn)%cols = ptrc%dat
                     blocks%inters(nn)%idx = ptrr%idx
                  end select
               end select
               curr => curr%next
               curc => curc%next

               blocks%inters(nn)%nr_loc = 0
               do jj = 1, blocks%inters(nn)%nr
                  idx = blocks%inters(nn)%rows(jj)
                  if (inters(blocks%inters(nn)%idx)%rows(idx) >= head .and. inters(blocks%inters(nn)%idx)%rows(idx) <= tail) then
                     blocks%inters(nn)%nr_loc = blocks%inters(nn)%nr_loc + 1
                  endif
               enddo
               allocate (blocks%inters(nn)%rows_loc(blocks%inters(nn)%nr_loc))
               allocate (blocks%inters(nn)%glo2loc(blocks%inters(nn)%nr))
               blocks%inters(nn)%glo2loc = -1
               blocks%inters(nn)%nr_loc = 0
               do jj = 1, blocks%inters(nn)%nr
                  idx = blocks%inters(nn)%rows(jj)
                  ! write(*,*)inters(blocks%inters(nn)%idx)%rows(idx),head,tail
                  if (inters(blocks%inters(nn)%idx)%rows(idx) >= head .and. inters(blocks%inters(nn)%idx)%rows(idx) <= tail) then
                     blocks%inters(nn)%nr_loc = blocks%inters(nn)%nr_loc + 1
                     blocks%inters(nn)%rows_loc(blocks%inters(nn)%nr_loc) = jj ! rows_loc stores indices in rows
                     blocks%inters(nn)%glo2loc(jj) = blocks%inters(nn)%nr_loc
                  endif
               enddo
               allocate (blocks%inters(nn)%dat_loc(blocks%inters(nn)%nr_loc, blocks%inters(nn)%nc))
            enddo

            ! extract entries on an array of intersections for each block

            if (blocks%style == 1) then
               call Full_block_extraction(blocks, inters, ptree, msh, stats)
            else
               if (blocks%level_butterfly == 0) then
                  call LR_block_extraction(blocks, inters, ptree, msh, stats)
               else
                  call BF_block_extraction(blocks, inters, ptree, msh, stats)
               endif
            endif

            ! finalize the lists of lists of rows and columns for each block because they have been transferred to intersections
            call list_finalizer(blocks%lstr)
            call list_finalizer(blocks%lstc)
         end select
         cur => cur%next
      enddo
#ifdef HAVE_TASKLOOP
      !$omp end single
      !$omp end parallel
#endif
      n3 = OMP_get_wtime()
      stats%Time_Entry_BF = stats%Time_Entry_BF + n3-n2

      call MPI_barrier(ptree%Comm, ierr)
      n3 = OMP_get_wtime()
      ! redistribute from blocks' intersections to the global intersecions inters
      ! if (flag2D == 1) then ! if each intersection is only needed by one processor, the communication can be optimized
         call BPACK_all2all_inters(inters, lstblk, stats, ptree, ptree%nproc, Npmap, pmaps)
      ! else
      !    call BPACK_all2all_inters_optimized(inters, lstblk, stats, ptree, ptree%nproc, Npmap, pmaps)
      ! endif

      n4 = OMP_get_wtime()
      stats%Time_Entry_Comm = stats%Time_Entry_Comm + n4-n3

      ntot_loc = 0
      do nn = 1, Ninter
         if (allocated(inters(nn)%dat_loc)) then
            nr_loc = size(inters(nn)%dat_loc, 1)
            nc_loc = size(inters(nn)%dat_loc, 2)
            do jj = 1, nc_loc
            do ii = 1, nr_loc
               alldat_loc(ntot_loc + ii + (jj - 1)*nr_loc) = inters(nn)%dat_loc(ii, jj)
            enddo
            enddo
            ntot_loc = ntot_loc + nr_loc*nc_loc
         endif
      enddo

      ! deallocate intersections at each block
      cur => lstblk%head
      do ii = 1, lstblk%num_nods
         select type (ptr=>cur%item)
         type is (block_ptr)
            blocks => ptr%ptr
            do nn = 1, size(blocks%inters, 1)
               if (allocated(blocks%inters(nn)%dat)) deallocate (blocks%inters(nn)%dat)
               if (allocated(blocks%inters(nn)%dat_loc)) deallocate (blocks%inters(nn)%dat_loc)
               if (allocated(blocks%inters(nn)%rows)) deallocate (blocks%inters(nn)%rows)
               if (allocated(blocks%inters(nn)%cols)) deallocate (blocks%inters(nn)%cols)
               if (allocated(blocks%inters(nn)%rows_loc)) deallocate (blocks%inters(nn)%rows_loc)
               blocks%inters(nn)%nr = 0
               blocks%inters(nn)%nr_loc = 0
               blocks%inters(nn)%nc = 0
               blocks%inters(nn)%idx = 0
            enddo
            deallocate (blocks%inters)
         end select
         cur => cur%next
      enddo

      ! finalize the list of block_ptr
      call list_finalizer(lstblk)

      ! deallocate global intersections
      do nn = 1, Ninter
         if (allocated(inters(nn)%dat)) deallocate (inters(nn)%dat)
         if (allocated(inters(nn)%dat_loc)) deallocate (inters(nn)%dat_loc)
         if (allocated(inters(nn)%rows)) deallocate (inters(nn)%rows)
         if (allocated(inters(nn)%cols)) deallocate (inters(nn)%cols)
         if (allocated(inters(nn)%rows_loc)) deallocate (inters(nn)%rows_loc)
      enddo
      deallocate (inters)

      n5 = OMP_get_wtime()

      call list_finalizer(lstr)
      call list_finalizer(lstc)
      ! time_tmp = time_tmp + n2- n1
      ! if(ptree%MyID==Main_ID)then
      ! write(*,*)n1-n0,n2-n1,n3-n2,n4-n3,n5-n4
      ! endif

   end subroutine BPACK_ExtractElement

!!!!!!! check error of BPACK construction using parallel element extraction
   subroutine BPACK_CheckError(bmat, option, msh, ker, stats, ptree)
      use BPACK_DEFS
      implicit none

      type(Hoption)::option
      type(Hstat)::stats
      type(Bmatrix)::bmat
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(intersect), allocatable::inters(:)
      real(kind=8)::n1, n2, n3, n4
      integer Ntest, passflag
      integer nsproc1, nsproc2, nprow, npcol, nprow1D, npcol1D, myrow, mycol, nprow1, npcol1, myrow1, mycol1, nprow2, npcol2, myrow2, mycol2, myArows, myAcols, rank1, rank2, ierr, MyID
      integer:: cridx, info
      integer, allocatable::rows(:), cols(:)
      DT, allocatable:: Vin(:, :), Vout1(:, :), Vout2(:, :), Vinter(:, :), Mat(:, :)
      integer::descVin(9), descVout(9), descVinter(9), descFull(9), descButterflyU(9), descButterflyV(9)
      integer N, M, i, j, ii, jj, nn, pp, myi, myj, iproc, jproc, rmax
      integer edge_n, edge_m, rank
      real(kind=8):: fnorm1, fnorm0, rtemp1 = 0, rtemp0 = 0
      real(kind=8):: a, v1, v2, v3
      DT:: value1, value2, value3
      type(list)::lstr, lstc, lst, lstblk
      type(nod), pointer::cur, curr, curc, curri, curci
      class(*), pointer::ptr, ptrr, ptrc, ptrri, ptrci
      integer::head, tail, idx, pgno, ctxt, nr_loc, nc_loc
      type(matrixblock), pointer::blocks
      integer num_blocks, idx_row, idx_col, idx_dat
      integer, allocatable:: allrows(:), allcols(:), pmaps(:, :)
      integer, allocatable::datidx(:), colidx(:), rowidx(:), pgidx(:)
      DT, allocatable::alldat_loc(:)
      integer:: Ninter, nr, nc, ntot_loc, level, Npmap, nproc, npavr, np
      type(intersect)::submats(1)

      ! select case(option%format)
      ! case(HODLR)
      ! level=bmat%ho_bf%Maxlevel
      ! case(HMAT)
      ! level=bmat%h_mat%Maxlevel
      ! end select

      ! level=1

      ! Ninter=2**level
      ! nr=2500
      ! nc=2500

      Ninter = 4
      ! nr=msh%Nunk
      ! nc=msh%Nunk

      nr = 100
      nc = 100

      allocate (colidx(Ninter))
      allocate (rowidx(Ninter))
      allocate (pgidx(Ninter))
      ! allocate(datidx(Ninter))

      allocate (allrows(Ninter*nr))
      allocate (allcols(Ninter*nc))

      ! pgno=1
      ! ctxt = ptree%pgrp(pgno)%ctxt
      ! call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
      ! nprow = ptree%pgrp(pgno)%nprow
      ! npcol = ptree%pgrp(pgno)%npcol

      nproc = ptree%nproc
      Npmap = min(Ninter, nproc)
      npavr = nproc/Npmap
      allocate (pmaps(Npmap, 3))
      do nn = 1, Npmap
         nprow = floor_safe(sqrt(dble(npavr)))
         npcol = floor_safe(npavr/dble(nprow))
         pmaps(nn, 1) = nprow
         pmaps(nn, 2) = npcol
         pmaps(nn, 3) = (nn - 1)*npavr
      enddo

      idx_row = 0
      idx_col = 0
      idx_dat = 0
      ! ntot_loc=0
      pp = 0
      do nn = 1, Ninter
         rowidx(nn) = nr
         colidx(nn) = nc
         pp = pp + 1
         pp = mod(pp - 1, Npmap) + 1
         pgidx(nn) = pp
         nprow = pmaps(pgidx(nn), 1)
         npcol = pmaps(pgidx(nn), 2)
         call Gridinfo_2D(pmaps(pgidx(nn), :), ptree%MyID, myrow, mycol)
         ! datidx(nn)=ntot_loc
         if (myrow /= -1 .and. mycol /= -1) then
            myArows = numroc_wp(nr, nbslpk, myrow, 0, nprow)
            myAcols = numroc_wp(nc, nbslpk, mycol, 0, npcol)
            idx_dat = idx_dat + myArows*myAcols
         endif

         do ii = 1, nr
            call random_number(a)
            call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, Main_ID, ptree%Comm, ierr)
            allrows(idx_row + 1) = max(floor_safe(msh%Nunk*a), 1)
            ! allrows(idx_row + 1) = max(floor_safe(3125*a), 1)+3125*0
            ! allrows(idx_row + 1) = max(floor_safe(7812*a), 1)+7812*0
            ! allrows(idx_row + 1) = max(floor_safe(19531*a), 1)+19531*0
            ! allrows(idx_row+1)=msh%basis_group(2**level+nn-1)%head+ii-1
            idx_row = idx_row + 1
         enddo

         do ii = 1, nc
            call random_number(a)
            call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, Main_ID, ptree%Comm, ierr)
            allcols(idx_col + 1) = max(floor_safe(msh%Nunk*a), 1)
            ! allcols(idx_col + 1) = max(floor_safe(3125*a), 1)+3125*1
            ! allcols(idx_col + 1) = max(floor_safe(7812*a), 1)+7812*1
            ! allcols(idx_col + 1) = max(floor_safe(19531*a), 1)+19531*1
            ! allcols(idx_col+1)=msh%basis_group(2**level+1-(nn-1))%head+ii-1
            idx_col = idx_col + 1
         enddo
      enddo

      allocate (alldat_loc(idx_dat))
      if (idx_dat > 0) alldat_loc = 0

      n1 = OMP_get_wtime()
      call BPACK_ExtractElement(bmat, option, msh, stats, ptree, Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps)
      n2 = OMP_get_wtime()

      ! compare extracted values with element_Zmn
      v1 = 0
      v2 = 0
      v3 = 0
      idx_row = 0
      idx_col = 0
      idx_dat = 0
      do nn = 1, Ninter
         nr = rowidx(nn)
         nc = colidx(nn)
         nprow = pmaps(pgidx(nn), 1)
         npcol = pmaps(pgidx(nn), 2)
         call Gridinfo_2D(pmaps(pgidx(nn), :), ptree%MyID, myrow, mycol)
         if (myrow /= -1 .and. mycol /= -1) then
            nr_loc = numroc_wp(nr, nbslpk, myrow, 0, nprow)
            nc_loc = numroc_wp(nc, nbslpk, mycol, 0, npcol)
            allocate (rows(nr_loc))
            allocate (cols(nc_loc))
            allocate (Mat(nr_loc, nc_loc))
            do myi = 1, nr_loc
               call l2g(myi, myrow, nr, nprow, nbslpk, ii)
               rows(myi) = allrows(ii + idx_row)
            enddo
            do myj = 1, nc_loc
               call l2g(myj, mycol, nc, npcol, nbslpk, jj)
               cols(myj) = allcols(jj + idx_col)
            enddo

            submats(1)%nr = nr_loc
            submats(1)%nc = nc_loc
            allocate(submats(1)%rows(submats(1)%nr))
            submats(1)%rows = rows(1:submats(1)%nr)
            allocate(submats(1)%cols(submats(1)%nc))
            submats(1)%cols = cols(1:submats(1)%nc)
            allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
            call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
            Mat = submats(1)%dat
            deallocate(submats(1)%rows)
            deallocate(submats(1)%cols)
            deallocate(submats(1)%dat)

            ! call element_Zmn_block_user(nr_loc, nc_loc, rows, cols, Mat, msh, option, ker, 0, passflag, ptree, stats)

            do myi = 1, nr_loc
            do myj = 1, nc_loc
               value2 = alldat_loc(idx_dat + myi + (myj - 1)*nr_loc)
               value1 = Mat(myi, myj)
               v1 = v1 + abs(value1)**2d0
               v2 = v2 + abs(value2)**2d0
               v3 = v3 + abs(value2 - value1)**2d0
               ! if(abs(value2-value1)**2d0/abs(value1)**2d0>1D-3)write(*,*)ptree%MyID,nn,myi,myj,abs(value2-value1)**2d0/abs(value1)**2d0,value2,value1
            enddo
            enddo
            idx_dat = idx_dat + nr_loc*nc_loc
         else
            nr_loc = 0
            nc_loc = 0
            allocate (rows(nr_loc))
            allocate (cols(nc_loc))
            allocate (Mat(nr_loc, nc_loc))
            call element_Zmn_blocklist_user(submats, 0, msh, option, ker, 2, passflag, ptree, stats)
            ! call element_Zmn_block_user(nr_loc, nc_loc, rows, cols, Mat, msh, option, ker, 2, passflag, ptree, stats)
         endif
         deallocate (rows)
         deallocate (cols)
         deallocate (Mat)
         idx_row = idx_row + nr
         idx_col = idx_col + nc
      enddo

      deallocate (rowidx)
      deallocate (colidx)
      deallocate (pgidx)
      ! deallocate(datidx)
      deallocate (allrows)
      deallocate (allcols)
      deallocate (alldat_loc)
      deallocate (pmaps)

      call MPI_ALLREDUCE(MPI_IN_PLACE, v1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, v2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, v3, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A25,Es14.7,Es14.7,A6,Es9.2,A7,Es9.2)') 'BPACK_CheckError: fnorm:', sqrt(v1), sqrt(v2), ' acc: ', sqrt(v3/v1), ' time: ', n2 - n1

   end subroutine BPACK_CheckError

!*********** all to all communication of element extraction results from local layout to 2D block-cyclic layout of each intersection (each process knows where to send, but doesn't know where to receive without communication)
   subroutine BPACK_all2all_inters(inters, lstblk, stats, ptree, nproc, Npmap, pmaps)

      use BPACK_DEFS
      implicit none
      integer i, j, k
      integer mm, nn, index_i, index_j, bb, ii, jj, ij, pp, tt, idx, idxs
      real(kind=8) flop
      type(Hstat)::stats
      type(proctree)::ptree
      integer ierr, nsendrecv, pid, tag, nproc, Nreqr, Nreqs, recvid, sendid
      type(commquant1D)::sendquant(nproc), recvquant(nproc)
      integer::sendactive(nproc), recvactive(nproc)
      integer::S_req(nproc),R_req(nproc)
      integer:: statuss(MPI_status_size, nproc), statusr(MPI_status_size,nproc)
      real(kind=8)::n1, n2, n3, n4, n5, n6, n7
      integer::sendIDactive(nproc), recvIDactive(nproc)
      integer Nsendactive, Nrecvactive
      integer::dist, pgno
      type(intersect)::inters(:)
      type(list)::lstblk
      type(nod), pointer::cur
      class(*), pointer::ptr
      type(matrixblock), pointer::blocks
      integer :: nprow, npcol, idstart, myi, myj, iproc, jproc, myrow, mycol, ri, ci
      integer,allocatable:: tmpidx(:,:)
      DT::val
      integer::Npmap, pmaps(Npmap, 3)
      integer*8:: cnt
      type(dat_pack)::foo(2)

      integer :: blocklen(2), type(2), newtype,newarrtype
      integer(KIND=MPI_ADDRESS_KIND) :: disp(2), base, lb, extent

      ! call MPI_GET_ADDRESS(foo(1)%idx(1), disp(1), ierr)
      ! call MPI_GET_ADDRESS(foo(1)%dat(1), disp(2), ierr)
      ! base = disp(1)
      ! disp(1) = disp(1) - base
      ! disp(2) = disp(2) - base
      ! blocklen(1) = 3
      ! blocklen(2) = 1
      ! type(1) = MPI_INTEGER
      ! type(2) = MPI_DT
      ! call MPI_TYPE_CREATE_STRUCT(2, blocklen, disp, type, newtype, ierr)
      ! call MPI_TYPE_COMMIT(newtype, ierr)
      ! call MPI_GET_ADDRESS(foo(1), disp(1), ierr)
      ! call MPI_GET_ADDRESS(foo(2), disp(2), ierr)
      ! extent = disp(2) - disp(1)
      ! lb = 0
      ! call MPI_TYPE_CREATE_RESIZED(newtype, lb, extent, newarrtype, ierr)
      ! call MPI_TYPE_COMMIT(newarrtype, ierr)


#ifdef HAVE_TASKLOOP
      !$omp parallel
      !$omp single
#endif
      n1 = OMP_get_wtime()
      pgno = 1
      ! nproc = ptree%pgrp(pgno)%nproc
      tag = pgno

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
      sendactive = 0
      recvactive = 0

      ! calculate send buffer sizes in the first pass
      cur => lstblk%head
      do bb = 1, lstblk%num_nods
         select type (ptr=>cur%item)
         type is (block_ptr)
            blocks => ptr%ptr
            do nn = 1, size(blocks%inters, 1)
               idx = blocks%inters(nn)%idx

               nprow = pmaps(inters(idx)%pg, 1)
               npcol = pmaps(inters(idx)%pg, 2)
               idstart = pmaps(inters(idx)%pg, 3)
               allocate(tmpidx(blocks%inters(nn)%nc,2))
               do jj=1,blocks%inters(nn)%nc
                  ci = blocks%inters(nn)%cols(jj)
                  call g2l(ci, inters(idx)%nc, npcol, nbslpk, tmpidx(jj,1), tmpidx(jj,2))
               enddo
               do ii = 1, blocks%inters(nn)%nr_loc
                  ri = blocks%inters(nn)%rows(blocks%inters(nn)%rows_loc(ii))
                  call g2l(ri, inters(idx)%nr, nprow, nbslpk, iproc, myi)
                  do jj = 1, blocks%inters(nn)%nc
                     ! ci = blocks%inters(nn)%cols(jj)
                     ! call g2l(ci, inters(idx)%nc, npcol, nbslpk, jproc, myj)
                     jproc = tmpidx(jj,1)
                     myj = tmpidx(jj,2)
                     pp = jproc*nprow + iproc + idstart + 1
                     if (sendquant(pp)%active == 0) then
                        sendquant(pp)%active = 1
                        sendactive(pp) = 1
                        Nsendactive = Nsendactive + 1
                        sendIDactive(Nsendactive) = pp
                     endif
                     sendquant(pp)%size = sendquant(pp)%size + 4                ! ri,ci,idx,value
                  enddo
               enddo
               deallocate(tmpidx)
            enddo
         end select
         cur => cur%next
      enddo

      ! compute recvquant(pp)%active by doing alltoall since receivers don't know where the data come from
      call MPI_ALLTOALL(sendactive, 1, MPI_INTEGER, recvactive, 1, MPI_INTEGER, ptree%pgrp(pgno)%Comm, ierr)
      Nreqr = 0
      do pp = 1, nproc
         if (recvactive(pp) == 1) then
            recvquant(pp)%active = 1
            Nrecvactive = Nrecvactive + 1
            recvIDactive(Nrecvactive) = pp
            sendid = pp - 1 + ptree%pgrp(pgno)%head
            if (sendid /= ptree%MyID) then
               Nreqr = Nreqr + 1
            endif
         endif
      enddo

      ! ! communicate receive buffer sizes
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
      enddo

      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         sendquant(pp)%size = 0
      enddo

      n2 = OMP_get_wtime()

      ! pack the send buffer in the second pass
      cur => lstblk%head
      do bb = 1, lstblk%num_nods
         select type (ptr=>cur%item)
         type is (block_ptr)
            blocks => ptr%ptr
#ifdef HAVE_TASKLOOP
            !$omp taskloop default(shared) private(nn,idx,nprow,npcol,idstart,cnt,ii,ri,iproc,myi,jj,ci,jproc,myj,pp,idxs)
#endif
            do nn = 1, size(blocks%inters, 1)
               idx = blocks%inters(nn)%idx
               nprow = pmaps(inters(idx)%pg, 1)
               npcol = pmaps(inters(idx)%pg, 2)
               idstart = pmaps(inters(idx)%pg, 3)
               ! allocate(tmpidx(blocks%inters(nn)%nc,2))
               ! do jj=1,blocks%inters(nn)%nc
               !    ci = blocks%inters(nn)%cols(jj)
               !    call g2l(ci, inters(idx)%nc, npcol, nbslpk, tmpidx(jj,1), tmpidx(jj,2))
               ! enddo
               !!$omp taskloop default(shared) private(cnt,ii,ri,iproc,myi,jj,ci,jproc,myj,pp,idxs)
               do cnt = 1, blocks%inters(nn)%nr_loc*blocks%inters(nn)%nc
                  jj = mod(cnt - 1, blocks%inters(nn)%nc) + 1
                  ii = (cnt - 1)/blocks%inters(nn)%nc + 1
                  ri = blocks%inters(nn)%rows(blocks%inters(nn)%rows_loc(ii))
                  call g2l(ri, inters(idx)%nr, nprow, nbslpk, iproc, myi)
                  ! jproc = tmpidx(jj,1)
                  ! myj = tmpidx(jj,2)
                  ci = blocks%inters(nn)%cols(jj)
                  call g2l(ci, inters(idx)%nc, npcol, nbslpk, jproc, myj)
                  pp = jproc*nprow + iproc + idstart + 1
                  !$omp atomic capture
                  idxs = sendquant(pp)%size
                  sendquant(pp)%size = sendquant(pp)%size + 4
                  !$omp end atomic
                  sendquant(pp)%dat(idxs + 1, 1) = myi
                  sendquant(pp)%dat(idxs + 2, 1) = myj
                  sendquant(pp)%dat(idxs + 3, 1) = idx
                  sendquant(pp)%dat(idxs + 4, 1) = blocks%inters(nn)%dat_loc(ii, jj)
               enddo
               !!$omp end taskloop
               ! deallocate(tmpidx)
            enddo
#ifdef HAVE_TASKLOOP
            !$omp end taskloop
#endif
         end select
         cur => cur%next
      enddo


      n3 = OMP_get_wtime()

      Nreqs = 0
      do tt = 1, Nsendactive
         ! n6 = OMP_get_wtime()
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
         endif
         ! n7 = OMP_get_wtime()
         ! write(*,*)ptree%MyID,'to',pp-1,n7-n6,sendquant(pp)%size
      enddo

      n4 = OMP_get_wtime()

      ! copy data from buffer to target
      cnt = 0
      do tt = 1, Nrecvactive
         if(tt==1 .and. Nreqr+1== Nrecvactive)then
            pp = ptree%MyID + 1 - ptree%pgrp(pgno)%head
            recvquant(pp)%size=sendquant(pp)%size
            if(recvquant(pp)%size>0)then
            allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
            recvquant(pp)%dat = sendquant(pp)%dat
            endif
         else
            call MPI_Probe(MPI_ANY_SOURCE, tag+1, ptree%pgrp(pgno)%Comm, statusr,ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
            call MPI_Get_count(statusr, MPI_DT, recvquant(pp)%size,ierr)
            allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
            cnt = cnt + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, R_req(cnt), ierr)
         endif
      enddo

      if (Nreqr > 0) then
         call MPI_waitall(Nreqr, R_req, statusr, ierr)
      endif
#ifdef HAVE_TASKLOOP
      !$omp taskloop default(shared) private(tt,pp,i,myi,myj,idx,val)
#endif
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         do i=1,recvquant(pp)%size/4
            myi = NINT(dble(recvquant(pp)%dat((i-1)*4+1, 1)))
            myj = NINT(dble(recvquant(pp)%dat((i-1)*4+2, 1)))
            idx = NINT(dble(recvquant(pp)%dat((i-1)*4+3, 1)))
            val = recvquant(pp)%dat((i-1)*4+4, 1)
            ! !$omp atomic
            inters(idx)%dat_loc(myi, myj) = inters(idx)%dat_loc(myi, myj) + val
            ! !$omp end atomic
         enddo
      enddo
#ifdef HAVE_TASKLOOP
      !$omp end taskloop
#endif



      n5 = OMP_get_wtime()
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

      ! n2 = OMP_get_wtime()
      ! time_tmp = time_tmp + n2 - n1
#ifdef HAVE_TASKLOOP
      !$omp end single
      !$omp end parallel
#endif
      ! write(*,*)n2-n1,n3-n2,n4-n3,n5-n4,'wordi',ptree%MyID

   end subroutine BPACK_all2all_inters

!*********** all to all communication of element extraction results from local layout to each entire intersection (each process knows where to send, but doesn't know where to receive without communication)
! YL: This subroutine seems to be slower than BPACK_all2all_inters
   subroutine BPACK_all2all_inters_optimized(inters, lstblk, stats, ptree, nproc, Npmap, pmaps)

      use BPACK_DEFS
      implicit none
      integer i, j, k
      integer mm, nn, index_i, index_j, bb, ii, jj, ij, pp, tt, idx, nr, nc, nr_max, nc_max
      real(kind=8) flop
      type(Hstat)::stats
      type(proctree)::ptree
      integer ierr, nsendrecv, pid, tag, nproc, Nreqr, Nreqs, recvid, sendid
      type(commquant1D)::sendquant(nproc), recvquant(nproc)
      integer::sendactive(nproc), recvactive(nproc)
      integer::S_req(nproc), R_req(nproc), Req
      integer:: statuss(MPI_status_size, nproc), statusr(MPI_status_size, nproc), status(MPI_status_size)
      real(kind=8)::n1, n2, n3, n4, n5
      integer::sendIDactive(nproc), recvIDactive(nproc)
      integer Nsendactive, Nrecvactive, Nsendactive_min, Nrecvactive_min
      logical all2all
      integer::sdispls(nproc), sendcounts(nproc), rdispls(nproc), recvcounts(nproc)
      DT, allocatable::sendbufall2all(:), recvbufall2all(:)
      integer::dist, pgno
      type(intersect)::inters(:)
      type(list)::lstblk
      type(nod), pointer::cur
      class(*), pointer::ptr
      type(matrixblock), pointer::blocks
      integer :: nprow, npcol, idstart, myi, myj, iproc, jproc, myrow, mycol, ri, ci, idxs
      DT::val
      integer::Npmap, pmaps(Npmap, 3), num_threads
      integer, save:: my_tid = 0
      integer, allocatable::ridx(:,:), cidx(:,:)
      integer*8:: cnt
#ifdef HAVE_TASKLOOP
      !$omp threadprivate(my_tid)
      !$omp parallel
      !$omp single
#endif
      num_threads = omp_get_num_threads()

      nr_max = 0
      nc_max = 0
      do nn = 1, size(inters, 1)
         nr_max = max(nr_max, inters(nn)%nr)
         nc_max = max(nc_max, inters(nn)%nc)
      enddo
      allocate (ridx(nr_max,num_threads))
      allocate (cidx(nc_max,num_threads))

      n1 = OMP_get_wtime()
      pgno = 1
      ! nproc = ptree%pgrp(pgno)%nproc
      tag = pgno

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
      sendactive = 0
      recvactive = 0

      ! calculate send buffer sizes in the first pass
      cur => lstblk%head
      do bb = 1, lstblk%num_nods
         select type (ptr=>cur%item)
         type is (block_ptr)
            blocks => ptr%ptr
            do nn = 1, size(blocks%inters, 1)
               idx = blocks%inters(nn)%idx
               pp = pmaps(inters(idx)%pg, 3) + 1
               nr = blocks%inters(nn)%nr_loc
               nc = blocks%inters(nn)%nc
               if (sendquant(pp)%active == 0) then
                  sendquant(pp)%active = 1
                  sendactive(pp) = 1
                  Nsendactive = Nsendactive + 1
                  sendIDactive(Nsendactive) = pp
               endif
               sendquant(pp)%size = sendquant(pp)%size + 3 + nr + nc + nr*nc                ! idx,nr,nc,indices(nr+nc),values(nr*nc)
            enddo
         end select
         cur => cur%next
      enddo

      ! compute recvquant(pp)%active by doing alltoall since receivers don't know where the data come from
#ifdef HAVE_MPI3
      call MPI_IALLTOALL(sendactive, 1, MPI_INTEGER, recvactive, 1, MPI_INTEGER, ptree%pgrp(pgno)%Comm, Req, ierr)
      call MPI_wait(Req, status, ierr)
#else
      call MPI_ALLTOALL(sendactive, 1, MPI_INTEGER, recvactive, 1, MPI_INTEGER, ptree%pgrp(pgno)%Comm, ierr)
#endif

      do pp = 1, nproc
         if (recvactive(pp) == 1) then
            recvquant(pp)%active = 1
            Nrecvactive = Nrecvactive + 1
            recvIDactive(Nrecvactive) = pp
         endif
      enddo

      ! communicate receive buffer sizes
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         allocate (sendquant(pp)%dat(sendquant(pp)%size, 1))
      enddo

      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
      enddo

      n2 = OMP_get_wtime()

      ! pack the send buffer in the second pass
      cur => lstblk%head
      do bb = 1, lstblk%num_nods
         select type (ptr=>cur%item)
         type is (block_ptr)
            blocks => ptr%ptr
#ifdef HAVE_TASKLOOP
            !$omp taskloop default(shared) private(nn,idx,idxs,pp,nr,nc,ii,myi,jj,myj,cnt)
#endif
            do nn = 1, size(blocks%inters, 1)
               idx = blocks%inters(nn)%idx
               pp = pmaps(inters(idx)%pg, 3) + 1
               nr = blocks%inters(nn)%nr_loc
               nc = blocks%inters(nn)%nc
#ifdef HAVE_TASKLOOP
               !$omp atomic capture
#endif
               idxs = sendquant(pp)%size
               sendquant(pp)%size = sendquant(pp)%size + 3 + nr + nc + nr*nc
#ifdef HAVE_TASKLOOP
               !$omp end atomic
#endif
               sendquant(pp)%dat(idxs + 1, 1) = idx
               sendquant(pp)%dat(idxs + 2, 1) = nr
               sendquant(pp)%dat(idxs + 3, 1) = nc
               do ii = 1, nr
                  myi = blocks%inters(nn)%rows(blocks%inters(nn)%rows_loc(ii))
                  sendquant(pp)%dat(idxs+3 + ii, 1) = myi
               enddo
               do jj = 1, nc
                  myj = blocks%inters(nn)%cols(jj)
                  sendquant(pp)%dat(idxs+3 +nr + jj, 1) = myj
               enddo
               do cnt = 1, nr*nc
                  jj = mod(cnt - 1, nc) + 1
                  ii = (cnt - 1)/nc + 1
                  sendquant(pp)%dat(idxs+3 +nr +nc + (ii - 1)*nc + jj, 1) = blocks%inters(nn)%dat_loc(ii, jj)
               enddo
            enddo
#ifdef HAVE_TASKLOOP
            !$omp end taskloop
#endif
         end select
         cur => cur%next
      enddo

      n3 = OMP_get_wtime()

      Nreqs = 0
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
         endif
      enddo

      Nreqr = 0
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         sendid = pp - 1 + ptree%pgrp(pgno)%head
         if (sendid /= ptree%MyID) then
            Nreqr = Nreqr + 1
            ! call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, R_req(Nreqr), ierr)
         endif
      enddo

      cnt = 0
      do tt = 1, Nrecvactive
         if(tt==1 .and. Nreqr+1== Nrecvactive)then
            pp = ptree%MyID + 1 - ptree%pgrp(pgno)%head
            recvquant(pp)%size=sendquant(pp)%size
            if(recvquant(pp)%size>0)then
            allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
            recvquant(pp)%dat = sendquant(pp)%dat
            endif
         else
            call MPI_Probe(MPI_ANY_SOURCE, tag+1, ptree%pgrp(pgno)%Comm, statusr,ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
            call MPI_Get_count(statusr, MPI_DT, recvquant(pp)%size,ierr)
            allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
            cnt = cnt + 1
            call MPI_Irecv(recvquant(pp)%dat, recvquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, R_req(cnt), ierr)
         endif
      enddo

      if (Nreqr > 0) then
         call MPI_waitall(Nreqr, R_req, statusr, ierr)
      endif
#ifdef HAVE_TASKLOOP
      !$omp taskloop default(shared) private(tt,pp,i,myi,myj,idx,nr,nc,ii,jj,cnt)
#endif
      do tt = 1, Nrecvactive
         my_tid = omp_get_thread_num()
         pp = recvIDactive(tt)
         i = 0
         do while (i < recvquant(pp)%size)
            i = i + 1
            idx = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            nr = NINT(dble(recvquant(pp)%dat(i, 1)))
            i = i + 1
            nc = NINT(dble(recvquant(pp)%dat(i, 1)))

            do ii = 1, nr
               ridx(ii,my_tid+1) = NINT(dble(recvquant(pp)%dat(i + ii, 1)))
            enddo
            i = i + nr
            do jj = 1, nc
               cidx(jj,my_tid+1) = NINT(dble(recvquant(pp)%dat(i + jj, 1)))
            enddo
            i = i + nc

            !!$omp taskloop default(shared) private(cnt,ii,jj)
            do cnt = 1, nr*nc
               jj = mod(cnt - 1, nc) + 1
               ii = (cnt - 1)/nc + 1
               !!$omp atomic
               inters(idx)%dat_loc(ridx(ii,my_tid+1), cidx(jj,my_tid+1)) = inters(idx)%dat_loc(ridx(ii,my_tid+1), cidx(jj,my_tid+1)) + recvquant(pp)%dat(i + (ii - 1)*nc + jj, 1)
               !!$omp end atomic
            enddo
            !!$omp end taskloop
            i = i + nr*nc

         enddo
      enddo
#ifdef HAVE_TASKLOOP
      !$omp end taskloop
#endif

      n5 = OMP_get_wtime()
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

      deallocate (ridx)
      deallocate (cidx)

#ifdef HAVE_TASKLOOP
      !$omp end single
      !$omp end parallel
#endif
      ! n2 = OMP_get_wtime()
      ! time_tmp = time_tmp + n2 - n1

      ! write(*,*)n2-n1,n3-n2,n4-n3,n5-n4,'wori',ptree%MyID

   end subroutine BPACK_all2all_inters_optimized

   recursive subroutine HODLR_MapIntersec2Block(ho_bf1, option, stats, msh, ptree, inters, nth, lstr, lstc, lstblk, level_c, bidx, flag)
      use BPACK_DEFS
      implicit none
      type(Hoption)::option
      type(Hstat)::stats
      type(hobf)::ho_bf1
      type(mesh)::msh
      type(proctree)::ptree
      type(intersect)::inters(:)
      integer nth, bidx, level_c
      integer ii, ll, idx, row_group, col_group
      type(list)::lstblk
      type(iarray)::lstr, lstc, clstr(2), clstc(2)
      type(nod), pointer::cur
      class(*), pointer::ptr
      type(matrixblock), pointer::blocks
      integer flag, num_nods
      type(block_ptr)::blk_ptr
      real(kind=8)::n1, n0

      if (flag == 0) then ! inverse blocks
         blocks => ho_bf1%levels(level_c)%BP_inverse(bidx)%LL(1)%matrices_block(1)
         row_group = blocks%row_group
         col_group = blocks%col_group
         if (IOwnPgrp(ptree, blocks%pgno)) then
         if (level_c == ho_bf1%Maxlevel + 1) then
            call HODLR_MapIntersec2Block(ho_bf1, option, stats, msh, ptree, inters, nth, lstr, lstc, lstblk, level_c, bidx, 1)
         else
            clstr(1)%idx = nth
            clstr(2)%idx = nth
            allocate (clstr(1)%dat(lstr%num_nods))
            allocate (clstr(2)%dat(lstr%num_nods))
            clstr(1)%num_nods = 0
            clstr(2)%num_nods = 0
            ! n0 = OMP_get_wtime()

            num_nods = msh%basis_group(row_group*2)%tail - msh%basis_group(row_group*2)%head + 1
            ! if(lstr%num_nods==msh%basis_group(row_group)%tail-msh%basis_group(row_group)%head+1)then
            ! clstr(1)%num_nods=num_nods
            ! clstr(1)%dat(1:clstr(1)%num_nods)=lstr%dat(1:clstr(1)%num_nods)
            ! clstr(2)%num_nods=lstr%num_nods-num_nods
            ! clstr(2)%dat(1:clstr(2)%num_nods)=lstr%dat(1+clstr(1)%num_nods:lstr%num_nods)
            ! else
            do ii = 1, lstr%num_nods
               ll = 2
               if (inters(nth)%rows(lstr%dat(ii)) <= msh%basis_group(row_group*2)%tail) ll = 1
               clstr(ll)%num_nods = clstr(ll)%num_nods + 1
               clstr(ll)%dat(clstr(ll)%num_nods) = lstr%dat(ii)
            enddo
            ! endif
            ! n1 = OMP_get_wtime()
            ! time_tmp = time_tmp + n1 - n0

            clstc(1)%idx = nth
            clstc(2)%idx = nth
            allocate (clstc(1)%dat(lstc%num_nods))
            allocate (clstc(2)%dat(lstc%num_nods))
            clstc(1)%num_nods = 0
            clstc(2)%num_nods = 0
            ! n0 = OMP_get_wtime()

            num_nods = msh%basis_group(col_group*2)%tail - msh%basis_group(col_group*2)%head + 1
            ! if(lstc%num_nods==msh%basis_group(col_group)%tail-msh%basis_group(col_group)%head+1)then
            ! clstc(1)%num_nods=num_nods
            ! clstc(1)%dat(1:clstc(1)%num_nods)=lstc%dat(1:clstc(1)%num_nods)
            ! clstc(2)%num_nods=lstc%num_nods-num_nods
            ! clstc(2)%dat(1:clstc(2)%num_nods)=lstc%dat(1+clstc(1)%num_nods:lstc%num_nods)
            ! else
            do ii = 1, lstc%num_nods
               ll = 2
               if (inters(nth)%cols(lstc%dat(ii)) <= msh%basis_group(col_group*2)%tail) ll = 1
               clstc(ll)%num_nods = clstc(ll)%num_nods + 1
               clstc(ll)%dat(clstc(ll)%num_nods) = lstc%dat(ii)
            enddo
            ! endif

            ! n1 = OMP_get_wtime()
            ! time_tmp = time_tmp + n1 - n0

            if (clstr(1)%num_nods > 0 .and. clstc(2)%num_nods > 0) call HODLR_MapIntersec2Block(ho_bf1, option, stats, msh, ptree, inters, nth, clstr(1), clstc(2), lstblk, level_c, 2*bidx - 1, 1)
            if (clstr(2)%num_nods > 0 .and. clstc(1)%num_nods > 0) call HODLR_MapIntersec2Block(ho_bf1, option, stats, msh, ptree, inters, nth, clstr(2), clstc(1), lstblk, level_c, 2*bidx, 1)
            if (clstr(1)%num_nods > 0 .and. clstc(1)%num_nods > 0) call HODLR_MapIntersec2Block(ho_bf1, option, stats, msh, ptree, inters, nth, clstr(1), clstc(1), lstblk, level_c + 1, 2*bidx - 1, 0)
            if (clstr(2)%num_nods > 0 .and. clstc(2)%num_nods > 0) call HODLR_MapIntersec2Block(ho_bf1, option, stats, msh, ptree, inters, nth, clstr(2), clstc(2), lstblk, level_c + 1, 2*bidx, 0)

            deallocate (clstc(1)%dat)
            deallocate (clstc(2)%dat)
            deallocate (clstr(1)%dat)
            deallocate (clstr(2)%dat)

         endif
         endif
      else ! forward blocks
         blocks => ho_bf1%levels(level_c)%BP(bidx)%LL(1)%matrices_block(1)
         row_group = blocks%row_group
         col_group = blocks%col_group

         if (IOwnPgrp(ptree, blocks%pgno)) then
            if (blocks%lstr%num_nods == 0) then
               blk_ptr%ptr => blocks
               call append(lstblk, blk_ptr)
            endif
            call append(blocks%lstr, lstr)
            call append(blocks%lstc, lstc)
         endif
      endif

   end subroutine HODLR_MapIntersec2Block

   recursive subroutine HSS_MapIntersec2Block(hss_bf1, option, stats, msh, ptree, inters, nth, lstr, lstc, lstblk, ll, Nbound)
      use BPACK_DEFS
      implicit none
      type(Hoption)::option
      type(Hstat)::stats
      type(hssbf)::hss_bf1
      type(mesh)::msh
      type(proctree)::ptree
      type(intersect)::inters(:)
      integer nth, bidx, level_c
      integer ii, ll, bb, row_group, col_group
      type(list)::lstblk
      type(iarray)::lstr, lstc
      type(nod), pointer::cur
      class(*), pointer::ptr
      type(matrixblock), pointer::blocks
      integer flag, num_nods
      type(block_ptr)::blk_ptr
      real(kind=8)::n1, n0
      integer, allocatable::rowblocks(:), colblocks(:)
      integer row0, col0, sort
      integer Nbound, level1, level0
      type(iarray)::clstr(Nbound), clstc(Nbound)

      do bb = 1, Nbound
         allocate (clstr(bb)%dat(lstr%num_nods))
         clstr(bb)%num_nods = 0
         clstr(bb)%idx = nth
         allocate (clstc(bb)%dat(lstc%num_nods))
         clstc(bb)%num_nods = 0
         clstc(bb)%idx = nth
      enddo

      allocate (rowblocks(hss_bf1%BP%LL(ll)%Nbound))
      allocate (colblocks(hss_bf1%BP%LL(ll)%Nbound))
      row0 = 0
      col0 = 0
      sort = 0
      do bb = 1, hss_bf1%BP%LL(ll)%Nbound
         rowblocks(bb) = hss_bf1%BP%LL(ll)%matrices_block(bb)%row_group
         colblocks(bb) = hss_bf1%BP%LL(ll)%matrices_block(bb)%col_group
         if (rowblocks(bb) < row0 .or. colblocks(bb) < col0) then
            sort = 1
            exit
         endif
         row0 = rowblocks(bb)
         col0 = colblocks(bb)
      enddo

      call assert(sort == 0, 'the rowblocks and colblocks need sorting first')

      level0 = hss_bf1%BP%LL(1)%matrices_block(1)%level
      level1 = hss_bf1%BP%LL(ll)%matrices_block(1)%level

      do ii = 1, lstr%num_nods
         row0 = findgroup(inters(nth)%rows(lstr%dat(ii)), msh, level1 - level0, hss_bf1%BP%LL(1)%matrices_block(1)%row_group)
         call binary_search(Nbound, rowblocks, row0, bb)
         if (bb /= -1) then
            clstr(bb)%num_nods = clstr(bb)%num_nods + 1
            clstr(bb)%dat(clstr(bb)%num_nods) = lstr%dat(ii)
         endif
      enddo

      do ii = 1, lstc%num_nods
         col0 = findgroup(inters(nth)%cols(lstc%dat(ii)), msh, level1 - level0, hss_bf1%BP%LL(1)%matrices_block(1)%col_group)
         call binary_search(Nbound, colblocks, col0, bb)
         if (bb /= -1) then
            clstc(bb)%num_nods = clstc(bb)%num_nods + 1
            clstc(bb)%dat(clstc(bb)%num_nods) = lstc%dat(ii)
         endif
      enddo

      do bb = 1, Nbound
         blocks => hss_bf1%BP%LL(ll)%matrices_block(bb)
         if (clstr(bb)%num_nods > 0 .and. clstc(bb)%num_nods > 0 .and. IOwnPgrp(ptree, blocks%pgno)) then
            if (blocks%lstr%num_nods == 0) then
               blk_ptr%ptr => blocks
               call append(lstblk, blk_ptr)
            endif
            call append(blocks%lstr, clstr(bb))
            call append(blocks%lstc, clstc(bb))
         endif
      enddo

      do bb = 1, Nbound
         deallocate (clstr(bb)%dat)
         deallocate (clstc(bb)%dat)
      enddo

      deallocate (rowblocks)
      deallocate (colblocks)

      if (ll < hss_bf1%BP%Lplus) then
         if (hss_bf1%BP%LL(ll + 1)%Nbound > 0) then
            call HSS_MapIntersec2Block(hss_bf1, option, stats, msh, ptree, inters, nth, lstr, lstc, lstblk, ll + 1, hss_bf1%BP%LL(ll + 1)%Nbound)
         endif
      endif

   end subroutine HSS_MapIntersec2Block

   subroutine Hmat_MapIntersec2Block(h_mat, option, stats, msh, ptree, inters, nth, lstr, lstc, lstblk, num_blocks)
      use BPACK_DEFS
      implicit none
      type(Hoption)::option
      type(Hstat)::stats
      type(Hmat)::h_mat
      type(mesh)::msh
      type(proctree)::ptree
      type(intersect)::inters(:)
      integer nth, num_blocks
      integer ii, jj, idx, row_group, col_group
      type(list)::lstblk
      type(iarray)::lstr, lstc, clstr_g, clstc_g(num_blocks)

      type(nod), pointer::cur
      class(*), pointer::ptr
      type(matrixblock), pointer::blocks
      type(block_ptr)::blk_ptr

      clstr_g%idx = nth
      clstr_g%num_nods = 0
      do jj = 1, num_blocks
         clstc_g(jj)%idx = nth
         clstc_g(jj)%num_nods = 0
      enddo

      allocate (clstr_g%dat(lstr%num_nods))
      clstr_g%num_nods = 0
      do ii = 1, lstr%num_nods
         row_group = h_mat%Local_blocks(1, 1)%row_group
         if (inters(nth)%rows(lstr%dat(ii)) >= msh%basis_group(row_group)%head .and. inters(nth)%rows(lstr%dat(ii)) <= msh%basis_group(row_group)%tail) then
            clstr_g%num_nods = clstr_g%num_nods + 1
            clstr_g%dat(clstr_g%num_nods) = lstr%dat(ii)
         endif
      enddo

      do jj = 1, num_blocks
         allocate (clstc_g(jj)%dat(lstc%num_nods))
         clstc_g(jj)%num_nods = 0
      enddo
      do ii = 1, lstc%num_nods
         jj = findgroup(inters(nth)%cols(lstc%dat(ii)), msh, msh%Dist_level, 1) - 2**msh%Dist_level + 1
         clstc_g(jj)%num_nods = clstc_g(jj)%num_nods + 1
         clstc_g(jj)%dat(clstc_g(jj)%num_nods) = lstc%dat(ii)
      enddo

      if (clstr_g%num_nods > 0) then
      do jj = 1, num_blocks
         if (clstc_g(jj)%num_nods > 0) then
            blocks => h_mat%Local_blocks(jj, 1)
            call Hmat_MapIntersec2Block_Loc(blocks, option, stats, msh, ptree, inters, nth, clstr_g, clstc_g(jj), lstblk)
         endif
      enddo
      endif

      deallocate (clstr_g%dat)
      do jj = 1, num_blocks
         deallocate (clstc_g(jj)%dat)
      enddo
   end subroutine Hmat_MapIntersec2Block

   recursive subroutine Hmat_MapIntersec2Block_Loc(blocks, option, stats, msh, ptree, inters, nth, lstr, lstc, lstblk)
      use BPACK_DEFS
      implicit none
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(proctree)::ptree
      type(intersect)::inters(:)
      integer nth, bidx, level_c
      integer ii, jj, ll, idx, row_group, col_group
      type(list)::lstblk
      type(iarray)::lstr, lstc, clstr(2), clstc(2)
      type(nod), pointer::cur
      class(*), pointer::ptr
      type(matrixblock), pointer::blocks, blocks_son
      integer flag
      type(block_ptr)::blk_ptr

      row_group = blocks%row_group
      col_group = blocks%col_group
      if (IOwnPgrp(ptree, blocks%pgno)) then
         if (blocks%style == 4) then ! divided blocks
            clstr(1)%idx = nth
            clstr(2)%idx = nth
            clstr(1)%num_nods = 0
            clstr(2)%num_nods = 0

            do ii = 1, lstr%num_nods
               ll = 2
               if (inters(nth)%rows(lstr%dat(ii)) <= msh%basis_group(row_group*2)%tail) ll = 1
               clstr(ll)%num_nods = clstr(ll)%num_nods + 1
            enddo
            allocate (clstr(1)%dat(clstr(1)%num_nods))
            allocate (clstr(2)%dat(clstr(2)%num_nods))
            clstr(1)%num_nods = 0
            clstr(2)%num_nods = 0
            do ii = 1, lstr%num_nods
               ll = 2
               if (inters(nth)%rows(lstr%dat(ii)) <= msh%basis_group(row_group*2)%tail) ll = 1
               clstr(ll)%num_nods = clstr(ll)%num_nods + 1
               clstr(ll)%dat(clstr(ll)%num_nods) = lstr%dat(ii)
            enddo

            clstc(1)%idx = nth
            clstc(2)%idx = nth
            clstc(1)%num_nods = 0
            clstc(2)%num_nods = 0
            do ii = 1, lstc%num_nods
               ll = 2
               if (inters(nth)%cols(lstc%dat(ii)) <= msh%basis_group(col_group*2)%tail) ll = 1
               clstc(ll)%num_nods = clstc(ll)%num_nods + 1
            enddo
            allocate (clstc(1)%dat(clstc(1)%num_nods))
            allocate (clstc(2)%dat(clstc(2)%num_nods))
            clstc(1)%num_nods = 0
            clstc(2)%num_nods = 0
            do ii = 1, lstc%num_nods
               ll = 2
               if (inters(nth)%cols(lstc%dat(ii)) <= msh%basis_group(col_group*2)%tail) ll = 1
               clstc(ll)%num_nods = clstc(ll)%num_nods + 1
               clstc(ll)%dat(clstc(ll)%num_nods) = lstc%dat(ii)
            enddo

            do ii = 1, 2
            do jj = 1, 2
               blocks_son => blocks%sons(ii, jj)
               if (clstr(ii)%num_nods > 0 .and. clstc(jj)%num_nods > 0) call Hmat_MapIntersec2Block_Loc(blocks_son, option, stats, msh, ptree, inters, nth, clstr(ii), clstc(jj), lstblk)
            enddo
            enddo
         else
            if (blocks%lstr%num_nods == 0) then
               blk_ptr%ptr => blocks
               call append(lstblk, blk_ptr)
            endif
            call append(blocks%lstr, lstr)
            call append(blocks%lstc, lstc)
         endif
      endif
   end subroutine Hmat_MapIntersec2Block_Loc

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

end module BPACK_constr
