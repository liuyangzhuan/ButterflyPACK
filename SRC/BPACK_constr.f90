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

!> @file BPACK_constr.f90
!> @brief Top-level subroutines for constructing a BPACK (H/HODBF/HODLR/HSS-BF) matrix or Butterfly block


#include "ButterflyPACK_config.fi"
module BPACK_constr

! use Butterfly_exact
   use Bplus_compress
   use Bplus_randomizedop
   use BPACK_Solve_Mul


   !> a list of blocks for the GPU construction (HODLR_gpu_baca_blocks)
   type hodlr_gpu_blkptr
      type(matrixblock), pointer :: p => null()
   end type hodlr_gpu_blkptr

   !> seconds of the merges of a shared block: moving the decomposed side to
   !> the even layout, the TSQR, the other side's products
   real(kind=8) :: hodlr_gpu_merge_time(3) = 0

contains

!>**** user-defined subroutine to sample a list of intersections from the bmat of Z
   subroutine Zelem_block_Extraction(Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps, quant)
      implicit none

      class(*), pointer :: quant
      integer:: Ninter
      integer:: allrows(:), allcols(:)
      DT,target::alldat_loc(:)
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

!>**** Initialization of the construction phase
   ! Nunk is matrix dimension
   ! Permutation is the permutation vector returned
   ! Nunk_loc is the local number of rows/columns
   ! bmat is the meta-data storing the compressed matrix
   ! Coordinates(optional) of dimension dim*N is the array of Cartesian coordinates corresponding to each row or column
   ! tree(optional) is an array of leafsizes in a user-provided cluster tree. clustertree has length 2*nl with nl denoting level of the clustertree.
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
      integer Permutation(:)
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

      call init_random_seed()

      call assert(associated(ker%QuantApp), 'ker%QuantApp is not assigned')
      if(option%cpp==0)call assert(associated(ker%FuncZmnBlock) .or. associated(ker%FuncZmn) .or. associated(ker%FuncHMatVec), 'one of the following should be assigned: ker%FuncZmn, ker%FuncZmnBlock, ker%FuncHMatVec')

      stats%Flop_Fill = 0
      stats%Time_Fill = 0
      stats%Time_Entry = 0
      stats%Time_Entry_Traverse = 0
      stats%Time_Entry_BF = 0
      stats%Time_Entry_Comm = 0

      !>**** set thread number here
#ifdef HAVE_MPI
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'NUMBER_MPI=', ptree%nproc
#endif
      threads_num = 1
      CALL getenv("OMP_NUM_THREADS", strings)
      strings = TRIM(strings)
      if (LEN_TRIM(strings) > 0) then
         read (strings, *) threads_num
      endif
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'OMP_NUM_THREADS=', threads_num
#ifdef HAVE_OPENMP
      call OMP_set_num_threads(threads_num)
#endif

      msh%Nunk = Nunk

      t1 = MPI_Wtime()
      nlevel = 0
      if (present(tree)) then
         nlevel = ceiling_safe(log(dble(size(tree, 1)))/log(2d0))
         Maxlevel = nlevel
         allocate (msh%pretree(2**Maxlevel))
         msh%pretree(1:2**Maxlevel) = tree(1:2**Maxlevel)

         !>**** make 0-element node a 1-element node

         ! write(*,*)'before adjustment:',msh%pretree
         call assert(Nunk>=2**Maxlevel,'The incomplete tree cannot be made complete. Try decreasing tree levels')
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

      !>**** copy geometry points if present
      if (option%nogeo == 0 .or. option%nogeo == 4) then
         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "User-supplied kernel requiring reorder"
         call assert(present(Coordinates), 'geometry points should be provided if option%nogeo==0 or 4')
         call LogMemory(stats, SIZEOF(Coordinates)/1024.0d3) ! this assumes that the user will deallocate Coordinates at a much later time
         Ndim = size(Coordinates, 1)
         Dimn = Ndim
         allocate (msh%xyz(Dimn, 1:msh%Nunk))
         call LogMemory(stats, SIZEOF(msh%xyz)/1024.0d3)
         msh%xyz = Coordinates
      endif

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "    "
      t2 = MPI_Wtime()

      t1 = MPI_Wtime()
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Hierarchical format......"
      call Cluster_partition(bmat, option, msh, ker, stats, ptree)
      call BPACK_structuring(bmat, option, msh, ker, ptree, stats)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Hierarchical format finished"
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "    "
      t2 = MPI_Wtime()

      if ((option%nogeo == 3 .or. option%nogeo == 4) .and. option%knn > 0) then
         call assert(present(nns), 'nearest neighbours should be provided if option%nogeo==3 or 4')
         allocate (msh%nns(msh%Nunk, option%knn))
         call LogMemory(stats, SIZEOF(msh%nns)/1024.0d3)
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

      !>**** return the permutation vector
      Nunk_loc = msh%idxe - msh%idxs + 1
      if (ptree%MyID == Main_ID) then
         do edge = 1, Nunk
            Permutation(edge) = msh%new2old(edge)
         enddo
         call LogMemory(stats, SIZEOF(Permutation)/1024.0d3) ! this assumes that the user will deallocate Permutation at a much later time
      endif

   end subroutine BPACK_construction_Init






!>**** Initialization of the construction phase for tensor butterfly-based hierarchical tensors
   ! Nunk(Ndim) is tensor size in each dimension
   ! Ndim is dimensionility
   ! Permutation(maxval(Nunk),Ndim) is the permutation vector (per dimension) returned
   ! Nunk_loc(Ndim) is the local number of indices per dimension
   ! bmat is the meta-data storing the compressed matrix
   ! Coordinates(optional) of dimension dim*max(Nunk) is the array of Cartesian coordinates corresponding to each dimension
   ! tree(optional) is an array of leafsizes in a user-provided cluster tree (per dimension). tree(:,dim_i) has length 2*nl with nl denoting level of the tree of dimension dim_i.
   ! If tree is incomplete with 0 element, ButterflyPACK will adjust it to a complete tree and return a modified tree.
   ! If the hierarchical matrix has more levels than tree, the code will generate more levels according to option%xyzsort, option%nogeo, and option%Nmin_leaf
   ! nns(optional) of dimension option%knn*max(Nunk)*Ndim is the array of user provided nearest neighbours
   subroutine BPACK_MD_construction_Init(Nunk, Ndim, Permutation, Nunk_loc, bmat, option, stats, msh, ker, ptree, Coordinates, tree, nns)
      implicit none
      integer Ndim,Nunk(Ndim)
      real(kind=8), optional:: Coordinates(:, :)
      integer, optional:: nns(:, :, :)

      real(kind=8) para
      real(kind=8) tolerance
      integer nn, mm, Maxlevel, give, need
      integer i, j, k, ii, edge, Dimn, kk, dim_i, ierr
      integer nlevel, level
      integer Permutation(:,:)
      integer, optional:: tree(:,:)
      integer Nunk_loc(Ndim)

      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh(Ndim)
      type(kernelquant)::ker
      type(Bmatrix)::bmat
      type(proctree)::ptree

      real(kind=8) t1, t2, structure_time, structure_time_max
      character(len=1024)  :: strings
      integer threads_num

      call init_random_seed()

      call assert(associated(ker%QuantApp), 'ker%QuantApp is not assigned')
      if(option%cpp==0)call assert(associated(ker%FuncZmn_MD) , 'one of the following should be assigned: ker%FuncZmn_MD, ker%FuncZmn_MDBlock')

      stats%Flop_Fill = 0
      stats%Time_Fill = 0
      stats%Time_Entry = 0
      stats%Time_Entry_Traverse = 0
      stats%Time_Entry_BF = 0
      stats%Time_Entry_Comm = 0

      !>**** set thread number here
#ifdef HAVE_MPI
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'NUMBER_MPI=', ptree%nproc
#endif
      threads_num = 1
      CALL getenv("OMP_NUM_THREADS", strings)
      strings = TRIM(strings)
      if (LEN_TRIM(strings) > 0) then
         read (strings, *) threads_num
      endif
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'OMP_NUM_THREADS=', threads_num
#ifdef HAVE_OPENMP
      call OMP_set_num_threads(threads_num)
#endif
      do dim_i=1,Ndim
      msh(dim_i)%Nunk = Nunk(dim_i)
      enddo

      t1 = MPI_Wtime()
      nlevel = 0
      if (present(tree)) then
         do dim_i=1,Ndim
            nlevel = ceiling_safe(log(dble(size(tree, 1)))/log(2d0))
            Maxlevel = nlevel
            allocate (msh(dim_i)%pretree(2**Maxlevel))
            msh(dim_i)%pretree(1:2**Maxlevel) = tree(1:2**Maxlevel,dim_i)

            !>**** make 0-element node a 1-element node

            ! write(*,*)'before adjustment:',msh%pretree
            call assert(Nunk(dim_i)>=2**Maxlevel,'The incomplete tree cannot be made complete. Try decreasing tree levels')
            need = 0
            do ii = 1, 2**Maxlevel
               if (msh(dim_i)%pretree(ii) == 0) need = need + 1
            enddo
            do while (need > 0)
               give = ceiling_safe(need/dble(2**Maxlevel - need))
               do ii = 1, 2**Maxlevel
                  nn = msh(dim_i)%pretree(ii)
                  if (nn > 1) then
                     msh(dim_i)%pretree(ii) = msh(dim_i)%pretree(ii) - min(min(nn - 1, give), need)
                     need = need - min(min(nn - 1, give), need)
                  endif
               enddo
            enddo
            do ii = 1, 2**Maxlevel
               if (msh(dim_i)%pretree(ii) == 0) msh(dim_i)%pretree(ii) = 1
            enddo
            ! write(*,*)'after adjustment:',msh%pretree
            tree(1:2**Maxlevel,dim_i) = msh(dim_i)%pretree(1:2**Maxlevel)
         enddo
      endif

      !>**** copy geometry points if present
      if (option%nogeo == 0 .or. option%nogeo == 4) then
         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "User-supplied kernel requiring reorder"
         call assert(present(Coordinates), 'geometry points should be provided if option%nogeo==0 or 4')
         call assert(Ndim == size(Coordinates, 1),'Ndim should be equal to the first dimension of Coordinates')
         do dim_i=1,Ndim
            allocate (msh(dim_i)%xyz(1, 1:msh(dim_i)%Nunk))
            call LogMemory(stats, SIZEOF(msh(dim_i)%xyz)/1024.0d3)
            msh(dim_i)%xyz(1,:) = Coordinates(dim_i,1:msh(dim_i)%Nunk)
         enddo
      endif

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "    "
      t2 = MPI_Wtime()

      t1 = MPI_Wtime()
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Hierarchical format......"
      call Cluster_partition_MD(Ndim, bmat, option, msh, ker, stats, ptree)
      call BPACK_structuring_MD(Ndim, bmat, option, msh, ker, ptree, stats)
      t2 = MPI_Wtime()
      structure_time = t2 - t1
      call MPI_Allreduce(structure_time, structure_time_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Hierarchical format finished"
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) &
         "Hierarchical format time:", structure_time_max, "Seconds"
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "    "

      if ((option%nogeo == 3 .or. option%nogeo == 4) .and. option%knn > 0) then
         call assert(present(nns), 'nearest neighbours should be provided if option%nogeo==3 or 4')
         do dim_i=1,Ndim
         allocate (msh(dim_i)%nns(msh(dim_i)%Nunk, option%knn))
         call LogMemory(stats, SIZEOF(msh(dim_i)%nns)/1024.0d3)
         do ii = 1, msh(dim_i)%Nunk
         do kk = 1, option%knn
            if (nns(kk, msh(dim_i)%new2old(ii),dim_i) /= 0) then
               msh(dim_i)%nns(ii, kk) = msh(dim_i)%old2new(nns(kk, msh(dim_i)%new2old(ii),dim_i))
            else
               msh(dim_i)%nns(ii, kk) = 0
            endif
         enddo
         enddo
         enddo
      endif

      !>**** return the permutation vector
      do dim_i=1,Ndim
      Nunk_loc(dim_i) = msh(dim_i)%idxe - msh(dim_i)%idxs + 1
      if (ptree%MyID == Main_ID) then
         do edge = 1, Nunk(dim_i)
            Permutation(edge,dim_i) = msh(dim_i)%new2old(edge)
         enddo
      endif
      enddo

   end subroutine BPACK_MD_construction_Init



!>**** Interface of BF construction via blackbox matvec or entry extraction with mshr and mshc already generated
   !> @param M,N: matrix size (in)
   !> @param M_loc,N_loc: number of local row/column indices (out)
   !> @param blocks: the structure containing the block (out)
   !> @param option: the structure containing option (in)
   !> @param stats: the structure containing statistics (in)
   !> @param msh: the structure containing points and ordering information combined from mshr and mshc (out)
   !> @param mshr: the structure containing points and ordering information for the row dimension (in)
   !> @param mshc: the structure containing points and ordering information for the column dimension (in)
   !> @param ker: the structure containing kernel quantities (in)
   !> @param ptree: the structure containing process tree (in)
   !> @param nns_m: (optional) (DIM knn*M) nearest neighbours(indexed from 1 to N) for each row (from 1 to M) provided by user (referenced if nogeo=3 or 4)
   !> @param nns_n: (optional) (DIM knn*N) nearest neighbours(indexed from 1 to M) for each column (from 1 to N) provided by user (referenced if nogeo=3 or 4)
   subroutine BF_Construct_Init_from_mshrc(M, N, M_loc, N_loc, mshr, mshc, blocks, option, stats, msh, ker, ptree, nns_m, nns_n)
      implicit none
      integer M, N

      integer Maxlevel
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam
      integer, allocatable:: groupmembers(:)
      integer level
      integer M_loc, N_loc
      integer, optional:: nns_m(:, :),nns_n(:, :)


      type(Hoption) ::option
      type(Hstat) ::stats
      type(mesh) ::msh, mshr, mshc
      type(kernelquant) ::ker
      type(matrixblock) ::blocks
      type(proctree) ::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2
      character(len=1024)  :: strings
      integer Maxgroup_rc
      integer(kind=8)::idx,kk,knn


      call assert(mshr%Nunk == M, 'mshr%Nunk\=M')
      call assert(mshc%Nunk == N, 'mshc%Nunk\=N')
      Maxgroup_rc = min(mshc%Maxgroup, mshr%Maxgroup)
      ! call assert(mshc%Maxgroup==mshr%Maxgroup,'mshc%Maxgroup\=mshr%Maxgroup')

      msh%Nunk = N + M
      msh%Maxgroup = Maxgroup_rc*2 + 1
      allocate (msh%basis_group(msh%Maxgroup))
      msh%basis_group(1)%head = 1
      msh%basis_group(1)%tail = N + M
      call copy_basis_group(mshr%basis_group, 1, Maxgroup_rc, msh%basis_group, 2, msh%Maxgroup, 0)
      call copy_basis_group(mshc%basis_group, 1, Maxgroup_rc, msh%basis_group, 3, msh%Maxgroup, mshr%Nunk)
      ! msh%new2old maps the indices in a larger (M+N)x(M+N) matrix into the MxM and NxN matrices
      allocate (msh%new2old(msh%Nunk))
      do ii = 1, M
         msh%new2old(ii) = ii
      enddo
      do ii = 1 + M, N + M
         msh%new2old(ii) = -(ii - M)
      enddo
      !>**** generate msh%xyz(1:Dimn,-N:M), needed in KNN
      if (option%nogeo ==0 .or. option%nogeo ==4) then
         Dimn = size(mshr%xyz,1)
         allocate(msh%xyz(1:Dimn,-N:M))
         msh%xyz=0
         do ii=1,M
            msh%xyz(1:Dimn,ii) = mshr%xyz(1:Dimn,mshr%new2old(ii))
         enddo
         do ii=1,N
            msh%xyz(:,-ii) = mshc%xyz(:,mshc%new2old(ii))
         enddo
      endif

      !>**** construct a list of k-nearest neighbours for each point
      if (option%nogeo /= 3 .and. option%nogeo /= 4 .and. option%knn > 0) then
         call FindKNNs(option, msh, ker, stats, ptree, 2, 3)
      endif

      if ((option%nogeo == 3 .or. option%nogeo == 4) .and. option%knn > 0) then
         allocate (msh%nns(msh%Nunk, option%knn))
         do ii = 1, M
         do kk = 1, option%knn
            knn = option%knn
            if (nns_m(kk,ii) /= 0) then
               msh%nns(ii, kk) = nns_m(kk,ii) + M
            else
               msh%nns(ii, kk) = 0
            endif
         enddo
         enddo
         do ii = 1, N
         do kk = 1, option%knn
            knn = option%knn
            if (nns_m(kk,ii) /= 0) then
               msh%nns(ii + M, kk) = nns_n(kk,ii)
            else
               msh%nns(ii + M, kk) = 0
            endif
         enddo
         enddo
      endif

      blocks%level = 1
      blocks%col_group = 3
      blocks%row_group = 2
      blocks%pgno = 1
      blocks%headm = msh%basis_group(blocks%row_group)%head
      blocks%headn = msh%basis_group(blocks%col_group)%head
      blocks%M = M
      blocks%N = N
      blocks%style = 2

      Maxlevel = GetTreelevel(msh%Maxgroup) - 1

      if (blocks%level > option%LRlevel) then
         blocks%level_butterfly = 0 ! low rank below LRlevel
      else
         blocks%level_butterfly = Maxlevel - blocks%level   ! butterfly
      endif

      call ComputeParallelIndices(blocks, blocks%pgno, ptree, msh)

      M_loc = blocks%M_loc
      N_loc = blocks%N_loc

   end subroutine BF_Construct_Init_from_mshrc




!>**** Interface of BF construction via blackbox matvec or entry extraction
   !> @param M,N: matrix size (in)
   !> @param M_loc,N_loc: number of local row/column indices (out)
   !> @param blocks: the structure containing the block (out)
   !> @param option: the structure containing option (in)
   !> @param stats: the structure containing statistics (in)
   !> @param msh: the structure containing points and ordering information combined from mshr and mshc (out)
   !> @param ker: the structure containing kernel quantities (in)
   !> @param ptree: the structure containing process tree (in)
   !> @param Permutation_m,Permutation_n: the permutation vectors on the row and column dimensions (out)
   !> @param tree_m, tree_n: (optional) is an array of leafsizes in a user-provided cluster tree for the row and column dimensions (in)
   !> @param Coordinates_m, Coordinates_n: (optional) of dimension dim*N is the array of Cartesian coordinates corresponding to each row or column
   !> @param nns_m: (optional) (DIM knn*M) nearest neighbours(indexed from 1 to N) for each row (from 1 to M) provided by user (referenced if nogeo=3 or 4)
   !> @param nns_n: (optional) (DIM knn*N) nearest neighbours(indexed from 1 to M) for each column (from 1 to N) provided by user (referenced if nogeo=3 or 4)
   subroutine BF_Construct_Init(M, N, M_loc, N_loc, Permutation_m, Permutation_n, blocks, option, stats, msh, ker, ptree, Coordinates_m, Coordinates_n, tree_m, tree_n, nns_m, nns_n)
      implicit none
      integer M, N

      integer Permutation_m(M),Permutation_n(N)
      real(kind=8), optional:: Coordinates_m(:, :), Coordinates_n(:, :)
      integer, optional:: tree_m(:),tree_n(:)
      integer, optional:: nns_m(:, :),nns_n(:, :)

      integer Maxlevel
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam
      integer, allocatable:: groupmembers(:)
      integer level
      integer M_loc, N_loc

      type(Hoption) ::option
      type(Hstat) ::stats
      type(mesh) ::msh, mshr, mshc
      type(kernelquant) ::ker
      type(matrixblock) ::blocks
      type(proctree) ::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2
      character(len=1024)  :: strings
      integer Maxgroup_rc,verbosity_save
      integer(kind=8)::idx,kk,knn
      type(Bmatrix)::bmat_m,bmat_n ! dummy

#ifdef HAVE_MPI
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'NUMBER_MPI=', ptree%nproc
#endif
      threads_num = 1
      CALL getenv("OMP_NUM_THREADS", strings)
      strings = TRIM(strings)
      if (LEN_TRIM(strings) > 0) then
         read (strings, *) threads_num
      endif
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'OMP_NUM_THREADS=', threads_num
#ifdef HAVE_OPENMP
      call OMP_set_num_threads(threads_num)
#endif

      !>**** create a random seed
      ! call DATE_AND_TIME(values=times)     ! Get the current time
      ! seed_myid(1) = times(4)*(360000*times(5) + 6000*times(6) + 100*times(7) + times(8))
      ! seed_myid(1) = myid*1000
      ! call RANDOM_SEED(PUT=seed_myid)
      call init_random_seed()
      verbosity_save=option%verbosity
      option%verbosity=-1
      call BPACK_construction_Init(M, Permutation_m, M_loc, bmat_m, option, stats, mshr, ker, ptree, Coordinates_m, tree_m)
      call BPACK_construction_Init(N, Permutation_n, N_loc, bmat_n, option, stats, mshc, ker, ptree, Coordinates_n, tree_n)
      option%verbosity=verbosity_save

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Maxlevel_for_blocks_m:', bmat_m%Maxlevel
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'N_leaf_m:', int(mshr%Nunk/(2**bmat_m%Maxlevel))
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Maxlevel_for_blocks_n:', bmat_n%Maxlevel
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'N_leaf_n:', int(mshc%Nunk/(2**bmat_n%Maxlevel))
      call BPACK_delete(bmat_m)
      call BPACK_delete(bmat_n)

      call BF_Construct_Init_from_mshrc(M, N, M_loc, N_loc, mshr, mshc, blocks, option, stats, msh, ker, ptree, nns_m, nns_n)


   end subroutine BF_Construct_Init



!>**** Interface of BF (tensor) construction via blackbox matvec or entry extraction with mshr and mshc already generated
   !> @param Ndim: number of dimensions (in)
   !> @param M(Ndim),N(Ndim): tensor sizes per dimension (in)
   !> @param M_loc(Ndim),N_loc(Ndim): local tensor sizes per dimension (out)
   !> @param blocks: the structure containing the block (out)
   !> @param option: the structure containing option (in)
   !> @param stats: the structure containing statistics (in)
   !> @param msh(Ndim): the structure containing points and ordering information combined from mshr and mshc (out)
   !> @param mshr(Ndim): the structure containing points and ordering information for the row dimension (in)
   !> @param mshc(Ndim): the structure containing points and ordering information for the column dimension (in)
   !> @param ker: the structure containing kernel quantities (in)
   !> @param ptree: the structure containing process tree (in)
   !> @param nns_m: (optional) (DIM knn*max(M)*Ndim) nearest neighbours(indexed from 1 to N) for each row (from 1 to M) provided by user (referenced if nogeo=3 or 4)
   !> @param nns_n: (optional) (DIM knn*max(N)*Ndim) nearest neighbours(indexed from 1 to M) for each column (from 1 to N) provided by user (referenced if nogeo=3 or 4)
   subroutine BF_MD_Construct_Init_from_mshrc(Ndim, M, N, M_loc, N_loc, mshr, mshc, blocks, option, stats, msh, ker, ptree, nns_m, nns_n)
      implicit none
      integer Ndim
      integer M(Ndim), N(Ndim)

      integer Maxlevel
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam
      integer, allocatable:: groupmembers(:)
      integer level, dim_i
      integer M_loc(Ndim), N_loc(Ndim)
      integer, optional:: nns_m(:, :, :),nns_n(:, :, :)


      type(Hoption) ::option
      type(Hstat) ::stats
      type(mesh) ::msh(Ndim), mshr(Ndim), mshc(Ndim)
      type(kernelquant) ::ker
      type(matrixblock_MD) ::blocks
      type(proctree) ::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2
      character(len=1024)  :: strings
      integer Maxgroup_rc
      integer(kind=8)::idx,kk,knn

      do dim_i=1,Ndim
         call assert(mshr(dim_i)%Nunk == M(dim_i), 'mshr%Nunk\=M')
         call assert(mshc(dim_i)%Nunk == N(dim_i), 'mshc%Nunk\=N')
         Maxgroup_rc = min(mshc(dim_i)%Maxgroup, mshr(dim_i)%Maxgroup)
         ! call assert(mshc%Maxgroup==mshr%Maxgroup,'mshc%Maxgroup\=mshr%Maxgroup')

         msh(dim_i)%Nunk = N(dim_i) + M(dim_i)
         msh(dim_i)%Maxgroup = Maxgroup_rc*2 + 1
         allocate (msh(dim_i)%basis_group(msh(dim_i)%Maxgroup))
         msh(dim_i)%basis_group(1)%head = 1
         msh(dim_i)%basis_group(1)%tail = N(dim_i) + M(dim_i)
         call copy_basis_group(mshr(dim_i)%basis_group, 1, Maxgroup_rc, msh(dim_i)%basis_group, 2, msh(dim_i)%Maxgroup, 0)
         call copy_basis_group(mshc(dim_i)%basis_group, 1, Maxgroup_rc, msh(dim_i)%basis_group, 3, msh(dim_i)%Maxgroup, mshr(dim_i)%Nunk)
         ! msh%new2old maps the indices in a larger (M+N)x(M+N) matrix into the MxM and NxN matrices
         allocate (msh(dim_i)%new2old(msh(dim_i)%Nunk))
         do ii = 1, M(dim_i)
            msh(dim_i)%new2old(ii) = ii
         enddo
         do ii = 1 + M(dim_i), N(dim_i) + M(dim_i)
            msh(dim_i)%new2old(ii) = -(ii - M(dim_i))
         enddo
         !>**** generate msh%xyz(1:Dimn,-N:M), needed in KNN
         if (option%nogeo ==0 .or. option%nogeo ==4) then
            allocate(msh(dim_i)%xyz(1,-N(dim_i):M(dim_i)))
            msh(dim_i)%xyz=0
            do ii=1,M(dim_i)
               msh(dim_i)%xyz(1,ii) = mshr(dim_i)%xyz(1,mshr(dim_i)%new2old(ii))
            enddo
            do ii=1,N(dim_i)
               msh(dim_i)%xyz(:,-ii) = mshc(dim_i)%xyz(:,mshc(dim_i)%new2old(ii))
            enddo
         endif

         !>**** construct a list of k-nearest neighbours for each point
         if (option%nogeo /= 3 .and. option%nogeo /= 4 .and. option%knn > 0) then
            call FindKNNs(option, msh(dim_i), ker, stats, ptree, 2, 3)
         endif

         if ((option%nogeo == 3 .or. option%nogeo == 4) .and. option%knn > 0) then
            allocate (msh(dim_i)%nns(msh(dim_i)%Nunk, option%knn))
            do ii = 1, M(dim_i)
            do kk = 1, option%knn
               knn = option%knn
               if (nns_m(kk,ii,dim_i) /= 0) then
                  msh(dim_i)%nns(ii, kk) = nns_m(kk,ii, dim_i) + M(dim_i)
               else
                  msh(dim_i)%nns(ii, kk) = 0
               endif
            enddo
            enddo
            do ii = 1, N(dim_i)
            do kk = 1, option%knn
               knn = option%knn
               if (nns_m(kk,ii,dim_i) /= 0) then
                  msh(dim_i)%nns(ii + M(dim_i), kk) = nns_n(kk,ii,dim_i)
               else
                  msh(dim_i)%nns(ii + M(dim_i), kk) = 0
               endif
            enddo
            enddo
         endif
      enddo



      blocks%Ndim = Ndim
      allocate(blocks%row_group(Ndim))
      allocate(blocks%col_group(Ndim))
      blocks%col_group = 3
      blocks%row_group = 2
      blocks%level = 1
      blocks%pgno = 1
      blocks%style = 2

      allocate(blocks%M(Ndim))
      allocate(blocks%N(Ndim))
      allocate(blocks%headm(Ndim))
      allocate(blocks%headn(Ndim))

      do dim_i=1,Ndim
         blocks%headm(dim_i) = msh(dim_i)%basis_group(blocks%row_group(dim_i))%head
         blocks%headn(dim_i) = msh(dim_i)%basis_group(blocks%col_group(dim_i))%head
         blocks%M(dim_i) = M(dim_i)
         blocks%N(dim_i) = N(dim_i)
      enddo


      Maxlevel = GetTreelevel(msh(1)%Maxgroup) - 1

      if (blocks%level > option%LRlevel) then
         blocks%level_butterfly = 0 ! low rank below LRlevel
      else
         blocks%level_butterfly = Maxlevel - blocks%level   ! butterfly
      endif

      call ComputeParallelIndices_MD(blocks, blocks%pgno, Ndim, ptree, msh)

      M_loc = blocks%M_loc
      N_loc = blocks%N_loc

   end subroutine BF_MD_Construct_Init_from_mshrc



!>**** Interface of BF (tensor) construction via blackbox matvec or entry extraction
   !> @param Ndim: dimensionality (in)
   !> @param M(Ndim),N(Ndim): tensor sizes of each dimension (in)
   !> @param M_loc,N_loc: local tensor sizes of each dimension (out)
   !> @param blocks: the structure containing the block (out)
   !> @param option: the structure containing option (in)
   !> @param stats: the structure containing statistics (in)
   !> @param msh(Ndim): the structure containing points and ordering information combined from mshr(Ndim) and mshc(Ndim) (out)
   !> @param ker: the structure containing kernel quantities (in)
   !> @param ptree: the structure containing process tree (in)
   !> @param Permutation_m(max(M),Ndim),Permutation_n(max(N),Ndim): the permutation vectors in each dimension (out)
   !> @param tree_m, tree_n: (optional) (2^L,Ndim) is an array of leafsizes in a user-provided cluster tree for the row and column dimensions (in)
   !> @param Coordinates_m, Coordinates_n: (optional) of dimension Ndim*N is the array of Cartesian coordinates corresponding to each row or column
   !> @param nns_m: (optional) (DIM knn*max(M)*Ndim) nearest neighbours(indexed from 1 to N) for each row (from 1 to M) provided by user (referenced if nogeo=3 or 4)
   !> @param nns_n: (optional) (DIM knn*max(N)*Ndim) nearest neighbours(indexed from 1 to M) for each column (from 1 to N) provided by user (referenced if nogeo=3 or 4)
   subroutine BF_MD_Construct_Init(Ndim, M, N, M_loc, N_loc, Permutation_m, Permutation_n, blocks, option, stats, msh, ker, ptree, Coordinates_m, Coordinates_n, tree_m, tree_n, nns_m, nns_n)
      implicit none
      integer Ndim
      integer M(Ndim), N(Ndim), N_leaf_m(Ndim), N_leaf_n(Ndim)

      integer Permutation_m(:,:),Permutation_n(:,:)
      real(kind=8), optional:: Coordinates_m(:, :), Coordinates_n(:, :)
      integer, optional:: tree_m(:,:),tree_n(:,:)
      integer, optional:: nns_m(:, :, :),nns_n(:, :, :)

      integer Maxlevel
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam, dim_i
      integer, allocatable:: groupmembers(:)
      integer level
      integer M_loc(Ndim), N_loc(Ndim)

      type(Hoption) ::option
      type(Hstat) ::stats
      type(mesh) ::msh(Ndim), mshr(Ndim), mshc(Ndim)
      type(kernelquant) ::ker
      type(matrixblock_MD) ::blocks
      type(proctree) ::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2
      character(len=1024)  :: strings
      integer Maxgroup_rc,verbosity_save, Maxlevel_m, Maxlevel_n
      integer(kind=8)::idx,kk,knn
      type(Bmatrix)::bmat_m,bmat_n ! dummy

#ifdef HAVE_MPI
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'NUMBER_MPI=', ptree%nproc
#endif
      threads_num = 1
      CALL getenv("OMP_NUM_THREADS", strings)
      strings = TRIM(strings)
      if (LEN_TRIM(strings) > 0) then
         read (strings, *) threads_num
      endif
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'OMP_NUM_THREADS=', threads_num
#ifdef HAVE_OPENMP
      call OMP_set_num_threads(threads_num)
#endif

      !>**** create a random seed
      ! call DATE_AND_TIME(values=times)     ! Get the current time
      ! seed_myid(1) = times(4)*(360000*times(5) + 6000*times(6) + 100*times(7) + times(8))
      ! seed_myid(1) = myid*1000
      ! call RANDOM_SEED(PUT=seed_myid)
      call init_random_seed()
      verbosity_save=option%verbosity
      option%verbosity=-1
      call BPACK_MD_construction_Init(M, Ndim, Permutation_m, M_loc, bmat_m, option, stats, mshr, ker, ptree, Coordinates_m, tree_m)
      Maxlevel_m=bmat_m%Maxlevel
      call BPACK_delete(bmat_m)

      call BPACK_MD_construction_Init(N, Ndim, Permutation_n, N_loc, bmat_n, option, stats, mshc, ker, ptree, Coordinates_n, tree_n)
      call BPACK_delete(bmat_n)
      Maxlevel_n=bmat_n%Maxlevel
      option%verbosity=verbosity_save

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Maxlevel_for_blocks_m:', Maxlevel_m
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Maxlevel_for_blocks_n:', Maxlevel_n
      do dim_i=1,Ndim
         N_leaf_m(dim_i)=int(mshr(dim_i)%Nunk/(2**Maxlevel_m))
         N_leaf_n(dim_i)=int(mshc(dim_i)%Nunk/(2**Maxlevel_n))
      enddo
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'N_leaf_m:', N_leaf_m
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'N_leaf_n:', N_leaf_n

      call BF_MD_Construct_Init_from_mshrc(Ndim, M, N, M_loc, N_loc, mshr, mshc, blocks, option, stats, msh, ker, ptree, nns_m, nns_n)


   end subroutine BF_MD_Construct_Init



!>**** Interface of BP construction via entry extraction with mshr and mshc already generated
   !> @param M,N: matrix size (in)
   !> @param M_loc,N_loc: number of local row/column indices (out)
   !> @param BP: the structure containing the BP (out)
   !> @param option: the structure containing option (in)
   !> @param stats: the structure containing statistics (in)
   !> @param msh: the structure containing points and ordering information combined from mshr and mshc (out)
   !> @param mshr: the structure containing points and ordering information for the row dimension (in)
   !> @param mshc: the structure containing points and ordering information for the column dimension (in)
   !> @param ker: the structure containing kernel quantities (in)
   !> @param ptree: the structure containing process tree (in)
   !> @param nns_m: (optional) (DIM knn*M) nearest neighbours(indexed from 1 to N) for each row (from 1 to M) provided by user (referenced if nogeo=3 or 4)
   !> @param nns_n: (optional) (DIM knn*N) nearest neighbours(indexed from 1 to M) for each column (from 1 to N) provided by user (referenced if nogeo=3 or 4)
   subroutine BP_Construct_Init_from_mshrc(M, N, M_loc, N_loc, mshr, mshc, BP, option, stats, msh, ker, ptree, nns_m, nns_n)
      implicit none
      integer M, N

      integer Maxlevel
      integer i, j, k, ii, jj, ll, edge, threads_num, nth, Dimn, nmpi, ninc, acam, groupm_ll, level_butterfly_ll,level_ll,level_BP,levelm,groupm_start,Nboundall,bb,group_m,group_n,Ninadmissible_max,Ninadmissible_tot,gg,group_m1,nn,group_n1,Nboundall1,Nboundall2,groupm_start1,groupn_start1,cc,cnt
      integer, allocatable:: groupmembers(:)
      integer level,group_tmp(1)
      integer M_loc, N_loc
      integer, optional:: nns_m(:, :),nns_n(:, :)


      type(Hoption) ::option
      type(Hstat) ::stats
      type(mesh) ::msh, mshr, mshc
      type(kernelquant) ::ker
      type(blockplus) ::BP
      type(matrixblock),pointer ::blocks,block_f
      type(proctree) ::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2
      character(len=1024)  :: strings
      integer Maxgroup_rc,ierr
      integer(kind=8)::idx,kk,knn


      call assert(mshr%Nunk == M, 'mshr%Nunk\=M')
      call assert(mshc%Nunk == N, 'mshc%Nunk\=N')
      Maxgroup_rc = min(mshc%Maxgroup, mshr%Maxgroup)
      ! call assert(mshc%Maxgroup==mshr%Maxgroup,'mshc%Maxgroup\=mshr%Maxgroup')

      msh%Nunk = N + M
      msh%Maxgroup = Maxgroup_rc*2 + 1
      allocate (msh%basis_group(msh%Maxgroup))
      msh%basis_group(1)%head = 1
      msh%basis_group(1)%tail = N + M
      call copy_basis_group(mshr%basis_group, 1, Maxgroup_rc, msh%basis_group, 2, msh%Maxgroup, 0)
      call copy_basis_group(mshc%basis_group, 1, Maxgroup_rc, msh%basis_group, 3, msh%Maxgroup, mshr%Nunk)
      ! msh%new2old maps the indices in a larger (M+N)x(M+N) matrix into the MxM and NxN matrices
      allocate (msh%new2old(msh%Nunk))
      do ii = 1, M
         msh%new2old(ii) = ii
      enddo
      do ii = 1 + M, N + M
         msh%new2old(ii) = -(ii - M)
      enddo
      !>**** generate msh%xyz(1:Dimn,-N:M), needed in KNN
      if (option%nogeo ==0 .or. option%nogeo ==4) then
         Dimn = size(mshr%xyz,1)
         allocate(msh%xyz(1:Dimn,-N:M))
         msh%xyz=0
         do ii=1,M
            msh%xyz(1:Dimn,ii) = mshr%xyz(1:Dimn,mshr%new2old(ii))
         enddo
         do ii=1,N
            msh%xyz(:,-ii) = mshc%xyz(:,mshc%new2old(ii))
         enddo
      endif

      !>**** construct a list of k-nearest neighbours for each point
      if (option%nogeo /= 3 .and. option%nogeo /= 4 .and. option%knn > 0) then
         call FindKNNs(option, msh, ker, stats, ptree, 2, 3)
      endif

      if ((option%nogeo == 3 .or. option%nogeo == 4) .and. option%knn > 0) then
         allocate (msh%nns(msh%Nunk, option%knn))
         do ii = 1, M
         do kk = 1, option%knn
            knn = option%knn
            if (nns_m(kk,ii) /= 0) then
               msh%nns(ii, kk) = nns_m(kk,ii) + M
            else
               msh%nns(ii, kk) = 0
            endif
         enddo
         enddo
         do ii = 1, N
         do kk = 1, option%knn
            knn = option%knn
            if (nns_m(kk,ii) /= 0) then
               msh%nns(ii + M, kk) = nns_n(kk,ii)
            else
               msh%nns(ii + M, kk) = 0
            endif
         enddo
         enddo
      endif

      select case (option%format)
      case (HODLR)
         BP%Lplus=1
         BP%pgno=1
         allocate (BP%LL(LplusMax))
         do ll = 1, LplusMax
            BP%LL(ll)%Nbound = 0
         end do
         BP%LL(1)%Nbound = 1
         allocate (BP%LL(1)%matrices_block(1))
         blocks => BP%LL(1)%matrices_block(1)

         blocks%pgno=1
         blocks%level = 1
         blocks%col_group = 3
         blocks%row_group = 2
         blocks%pgno = 1
         blocks%headm = msh%basis_group(blocks%row_group)%head
         blocks%headn = msh%basis_group(blocks%col_group)%head
         blocks%M = M
         blocks%N = N
         blocks%style = 2

         Maxlevel = GetTreelevel(msh%Maxgroup) - 1

         if (blocks%level > option%LRlevel) then
            blocks%level_butterfly = 0 ! low rank below LRlevel
         else
            blocks%level_butterfly = Maxlevel - blocks%level   ! butterfly
         endif

         call ComputeParallelIndices(blocks, blocks%pgno, ptree, msh)

         M_loc = blocks%M_loc
         N_loc = blocks%N_loc


      case (HMAT,BLR)
         msh%Dist_level=0
         Maxlevel = GetTreelevel(msh%Maxgroup) - 1
         allocate (stats%leafs_of_level(0:Maxlevel))
         stats%leafs_of_level = 0

         BP%Lplus=1
         BP%pgno=1
         allocate (BP%LL(1))
         BP%LL(1)%Nbound = 1
         allocate (BP%LL(1)%matrices_block(1))
         blocks => BP%LL(1)%matrices_block(1)

         blocks%pgno=1
         blocks%level = 1
         blocks%col_group = 3
         blocks%row_group = 2
         blocks%headm = msh%basis_group(blocks%row_group)%head
         blocks%headn = msh%basis_group(blocks%col_group)%head
         blocks%M = M
         blocks%N = N
         if (blocks%level > option%LRlevel) then
            blocks%level_butterfly = 0 ! low rank below LRlevel
         else
            blocks%level_butterfly = Maxlevel - blocks%level   ! butterfly
         endif
         call ComputeParallelIndices(blocks, blocks%pgno, ptree, msh)
         M_loc = blocks%M_loc
         N_loc = blocks%N_loc
         call Hmat_construct_local_tree(blocks, option, stats, msh, ker, ptree, Maxlevel)
         call MPI_allreduce(MPI_IN_PLACE, stats%leafs_of_level(0:Maxlevel), Maxlevel + 1, MPI_integer, MPI_sum, ptree%Comm, ierr)

         allocate(blocks%lstblks(0:Maxlevel))
         do level = 0, Maxlevel
            blocks%lstblks(level) = list()
         enddo
         call Hmat_GetBlkLst(blocks, ptree, blocks%lstblks,Maxlevel)
         do level = 0, Maxlevel
            call MergeSort(blocks%lstblks(level)%head, node_score_block_ptr_row)
         enddo

      case (HSS)


         BP%level = 1
         BP%col_group = 3
         BP%row_group = 2
         BP%pgno = 1

         allocate (BP%LL(LplusMax))
         do ll = 1, LplusMax
            BP%LL(ll)%Nbound = 0
         end do

         Maxlevel = GetTreelevel(msh%Maxgroup) - 1

         BP%LL(1)%Nbound = 1
         allocate (BP%LL(1)%matrices_block(1))
         block_f => BP%LL(1)%matrices_block(1)
         block_f%level = BP%level
         block_f%level_butterfly = int((Maxlevel - block_f%level)/2)*2 ! butterfly plus needs even number of levels

         block_f%col_group = BP%col_group
         block_f%row_group = BP%row_group
         block_f%pgno = BP%pgno
         ! pgno_bplus=block_f%pgno

         block_f%M = msh%basis_group(block_f%row_group)%tail - msh%basis_group(block_f%row_group)%head + 1
         block_f%N = msh%basis_group(block_f%col_group)%tail - msh%basis_group(block_f%col_group)%head + 1
         block_f%headm = msh%basis_group(block_f%row_group)%head
         block_f%headn = msh%basis_group(block_f%col_group)%head

         call ComputeParallelIndices(block_f, block_f%pgno, ptree, msh)

         block_f%style = 2
         allocate (BP%LL(1)%boundary_map(1,1))
         BP%LL(1)%boundary_map(1,1) = block_f%col_group
         BP%Lplus = 0
         groupm_ll = block_f%row_group
         level_butterfly_ll = block_f%level_butterfly
         level_ll = GetTreelevel(groupm_ll) - 1
         msh%basis_group(groupm_ll)%nn = 1
         allocate(msh%basis_group(groupm_ll)%nlist(1))
         msh%basis_group(groupm_ll)%nlist(1) = block_f%col_group

         do ll = 1, LplusMax - 1
            if (BP%LL(ll)%Nbound > 0) then
               BP%Lplus = BP%Lplus + 1
               call assert(BP%Lplus <= LplusMax, 'increase LplusMax')

               if (ll == LplusMax - 1 .or. level_butterfly_ll == 0) then
                  BP%LL(ll + 1)%Nbound = 0
               else
                  level_BP = BP%level
                  levelm = ceiling_safe(dble(level_butterfly_ll)/2d0)


                  groupm_start = groupm_ll*2**levelm
                  Nboundall = 2**(level_ll + levelm - level_BP)
                  do bb = 1, Nboundall
                     group_m = bb + groupm_start - 1
                     msh%basis_group(group_m)%nn=0
                  enddo

                  Ninadmissible_max=0
                  Ninadmissible_tot=0
                  do gg=groupm_ll,groupm_ll+2**(level_ll) -1
                     group_m1 = gg
                     do nn=1,msh%basis_group(gg)%nn
                        group_n1 = msh%basis_group(gg)%nlist(nn)
                        Nboundall1 = 2**(levelm)
                        Nboundall2 = 2**(level_butterfly_ll-levelm)
                        groupm_start1 = group_m1*2**levelm
                        groupn_start1 = group_n1*2**(level_butterfly_ll-levelm)
                        do bb = 1, Nboundall1
                        do cc = 1, Nboundall2
                           group_m = bb + groupm_start1 - 1
                           group_n = cc + groupn_start1 - 1
                           if (near_or_far_user(group_m, group_n, msh, option, ker, option%near_para) == 0)then
                              msh%basis_group(group_m)%nn = msh%basis_group(group_m)%nn + 1
                              Ninadmissible_max = max(Ninadmissible_max,msh%basis_group(group_m)%nn)
                              Ninadmissible_tot = Ninadmissible_tot +1
                           endif
                        enddo
                        enddo
                     enddo
                  enddo

                  do bb = 1, Nboundall
                     group_m = bb + groupm_start - 1
                     allocate(msh%basis_group(group_m)%nlist(max(1,msh%basis_group(group_m)%nn)))
                     msh%basis_group(group_m)%nn=0
                  enddo

                  do gg=groupm_ll,groupm_ll+2**(level_ll) -1
                     group_m1 = gg
                     do nn=1,msh%basis_group(gg)%nn
                        group_n1 = msh%basis_group(gg)%nlist(nn)
                        Nboundall1 = 2**(levelm)
                        Nboundall2 = 2**(level_butterfly_ll-levelm)
                        groupm_start1 = group_m1*2**levelm
                        groupn_start1 = group_n1*2**(level_butterfly_ll-levelm)
                        do bb = 1, Nboundall1
                        do cc = 1, Nboundall2
                           group_m = bb + groupm_start1 - 1
                           group_n = cc + groupn_start1 - 1
                           if (near_or_far_user(group_m, group_n, msh, option, ker, option%near_para) == 0)then
                              msh%basis_group(group_m)%nn = msh%basis_group(group_m)%nn + 1
                              msh%basis_group(group_m)%nlist(msh%basis_group(group_m)%nn)=group_n
                           endif
                        enddo
                        enddo
                     enddo
                  enddo

                  allocate (BP%LL(ll + 1)%boundary_map(Nboundall,Ninadmissible_max))
                  BP%LL(ll + 1)%boundary_map=-1
                  do bb = 1, Nboundall
                     group_m = bb + groupm_start - 1
                     do nn=1,msh%basis_group(group_m)%nn
                        group_n = msh%basis_group(group_m)%nlist(nn)
                        BP%LL(ll + 1)%boundary_map(bb,nn)=group_n
                     enddo
                  enddo
                  BP%LL(ll + 1)%Nbound = Ninadmissible_tot


                  allocate (BP%LL(ll + 1)%matrices_block(BP%LL(ll + 1)%Nbound))
                  cnt = 0
                  do bb = 1, Nboundall
                     do jj=1,Ninadmissible_max
                        if (BP%LL(ll + 1)%boundary_map(bb,jj) /= -1) then
                           cnt = cnt + 1
                           group_m = bb + groupm_start - 1
                           group_n = BP%LL(ll + 1)%boundary_map(bb,jj)
                           blocks => BP%LL(ll + 1)%matrices_block(cnt)
                           blocks%row_group = group_m
                           blocks%col_group = group_n
                           blocks%level = GetTreelevel(group_m) - 1
                           blocks%level_butterfly = int((Maxlevel - blocks%level)/2)*2

                           group_tmp(1) = INT(group_m/2)  ! this lines ensures rowgroup 2 has pgno=1 instead pgno=2
                           blocks%pgno = GetMshGroup_Pgno(ptree,1,group_tmp)

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
                     enddo
                  end do
                  groupm_ll = groupm_ll*2**levelm
                  level_ll = GetTreelevel(groupm_ll) - 1
                  level_butterfly_ll = int((Maxlevel - level_ll)/2)*2
               end if
            else
               exit
            end if
         end do

         do gg=1,msh%Maxgroup
            if (allocated(msh%basis_group(gg)%nlist)) then
               deallocate(msh%basis_group(gg)%nlist)
               msh%basis_group(gg)%nn=0
            endif
         enddo

         call LogMemory(stats, SIZEOF(BP)/1024.0d3)

         M_loc = BP%LL(1)%matrices_block(1)%M_loc
         N_loc = BP%LL(1)%matrices_block(1)%N_loc


      end select

   end subroutine BP_Construct_Init_from_mshrc





!>**** Interface of BP construction via entry extraction
   !> @param M,N: matrix size (in)
   !> @param M_loc,N_loc: number of local row/column indices (out)
   !> @param BP: the structure containing the BP (out)
   !> @param option: the structure containing option (in)
   !> @param stats: the structure containing statistics (in)
   !> @param msh: the structure containing points and ordering information combined from mshr and mshc (out)
   !> @param ker: the structure containing kernel quantities (in)
   !> @param ptree: the structure containing process tree (in)
   !> @param Permutation_m,Permutation_n: the permutation vectors on the row and column dimensions (out)
   !> @param tree_m, tree_n: (optional) is an array of leafsizes in a user-provided cluster tree for the row and column dimensions (in)
   !> @param Coordinates_m, Coordinates_n: (optional) of dimension dim*N is the array of Cartesian coordinates corresponding to each row or column
   !> @param nns_m: (optional) (DIM knn*M) nearest neighbours(indexed from 1 to N) for each row (from 1 to M) provided by user (referenced if nogeo=3 or 4)
   !> @param nns_n: (optional) (DIM knn*N) nearest neighbours(indexed from 1 to M) for each column (from 1 to N) provided by user (referenced if nogeo=3 or 4)
   subroutine BP_Construct_Init(M, N, M_loc, N_loc, Permutation_m, Permutation_n, BP, option, stats, msh, ker, ptree, Coordinates_m, Coordinates_n, tree_m, tree_n, nns_m, nns_n)
      implicit none
      integer M, N

      integer Permutation_m(M),Permutation_n(N)
      real(kind=8), optional:: Coordinates_m(:, :), Coordinates_n(:, :)
      integer, optional:: tree_m(:),tree_n(:)
      integer, optional:: nns_m(:, :),nns_n(:, :)

      integer Maxlevel
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam
      integer, allocatable:: groupmembers(:)
      integer level
      integer M_loc, N_loc

      type(Hoption) ::option
      type(Hstat) ::stats
      type(mesh) ::msh, mshr, mshc
      type(kernelquant) ::ker
      type(blockplus) ::BP
      type(proctree) ::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2
      character(len=1024)  :: strings
      integer Maxgroup_rc,verbosity_save,format_save
      integer(kind=8)::idx,kk,knn
      type(Bmatrix)::bmat_m,bmat_n ! dummy

#ifdef HAVE_MPI
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'NUMBER_MPI=', ptree%nproc
#endif
      threads_num = 1
      CALL getenv("OMP_NUM_THREADS", strings)
      strings = TRIM(strings)
      if (LEN_TRIM(strings) > 0) then
         read (strings, *) threads_num
      endif
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'OMP_NUM_THREADS=', threads_num
#ifdef HAVE_OPENMP
      call OMP_set_num_threads(threads_num)
#endif

      !>**** create a random seed
      ! call DATE_AND_TIME(values=times)     ! Get the current time
      ! seed_myid(1) = times(4)*(360000*times(5) + 6000*times(6) + 100*times(7) + times(8))
      ! seed_myid(1) = myid*1000
      ! call RANDOM_SEED(PUT=seed_myid)
      call init_random_seed()
      verbosity_save=option%verbosity
      format_save=option%format
      option%verbosity=-1
      option%format=HODLR
      call BPACK_construction_Init(M, Permutation_m, M_loc, bmat_m, option, stats, mshr, ker, ptree, Coordinates_m, tree_m)
      call BPACK_construction_Init(N, Permutation_n, N_loc, bmat_n, option, stats, mshc, ker, ptree, Coordinates_n, tree_n)
      option%verbosity=verbosity_save
      option%format=format_save

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Maxlevel_for_blocks_m:', bmat_m%Maxlevel
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'N_leaf_m:', int(mshr%Nunk/(2**bmat_m%Maxlevel))
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Maxlevel_for_blocks_n:', bmat_n%Maxlevel
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'N_leaf_n:', int(mshc%Nunk/(2**bmat_n%Maxlevel))

      call BPACK_delete(bmat_m)
      call BPACK_delete(bmat_n)
      call BP_Construct_Init_from_mshrc(M, N, M_loc, N_loc, mshr, mshc, BP, option, stats, msh, ker, ptree, nns_m, nns_n)


   end subroutine BP_Construct_Init


!>**** Interface of BF construction via entry extraction
   !> @param blocks: the structure containing the block (inout)
   !> @param option: the structure containing option (in)
   !> @param stats: the structure containing statistics (inout)
   !> @param msh: the structure containing points and ordering information (in)
   !> @param ker: the structure containing kernel quantities (in)
   !> @param ptree: the structure containing process tree (in)
   subroutine BF_Construct_Element_Compute(blocks, option, stats, msh, ker, ptree)
      implicit none

      integer Maxlevel
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam
      integer, allocatable:: groupmembers(:)
      integer level

      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(matrixblock),target::blocks
      type(matrixblock),pointer::blocks_1
      type(proctree)::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2, error, Memory, tol_comp_tmp
      integer ierr,pp,knn_tmp
      integer:: boundary_map(1,1)
      integer groupm_start, Nboundall,Ninadmissible
      blocks_1 => blocks
      if (allocated(msh%xyz)) deallocate (msh%xyz)

      t1 = MPI_Wtime()
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) " "
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "EntryExtraction-based BF construction......"

      groupm_start = 0
      Nboundall = 0
      Ninadmissible = 0
      if (option%forwardN15flag == 1) then
         knn_tmp = option%knn
         option%knn=0
         call BF_compress_N15(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree, 1)
         option%knn=knn_tmp
         call BF_sym2asym(blocks)
      elseif (option%forwardN15flag == 2) then
         call BF_compress_NlogN(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree, 1)
         call BF_checkError(blocks_1, option, msh, ker, stats, ptree, 0, -1, error)
         if(error>50*option%tol_comp)then
            pp = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1
            if(option%verbosity>=0 .and. pp==1)write(*,*)'warning: error ',error,',  with the N15 algorithm'
            call BF_delete(blocks,0)
            knn_tmp = option%knn
            option%tol_comp=option%tol_comp*5
            option%knn=0
            call BF_compress_N15(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree, 1)
            option%knn=knn_tmp
            option%tol_comp=option%tol_comp/5
            call BF_sym2asym(blocks)
         endif
      else
         call BF_compress_NlogN(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree, 1)
      end if

      if (option%verbosity >= 0) call BF_checkError(blocks_1, option, msh, ker, stats, ptree, 0, option%verbosity)
      call BF_ComputeMemory(blocks, stats%Mem_Comp_for)

      if (option%verbosity >= 2)call BF_print_size(blocks)

      !>**** delete neighours in msh
      if(allocated(msh%nns))then
         call LogMemory(stats, -SIZEOF(msh%nns)/1024.0d3)
         deallocate(msh%nns)
      endif

      t2 = MPI_Wtime()

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "EntryExtraction-based BF construction finished in", t2-t1, 'Seconds with', Memory,'MB Memory'

      if (.not. allocated(stats%rankmax_of_level)) allocate (stats%rankmax_of_level(0:0))
      if (.not. allocated(stats%rankmax_of_level_global)) allocate (stats%rankmax_of_level_global(0:0))
      stats%rankmax_of_level(0) = blocks%rankmax
      stats%rankmax_of_level_global(0) = stats%rankmax_of_level(0)
      stats%Mem_Fill = stats%Mem_Comp_for + stats%Mem_Direct_for
      call LogMemory(stats, stats%Mem_Fill)
      stats%Time_Fill = stats%Time_Fill + t2 - t1
      ! stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp ! Flop_Fill already counted in BF_compress_NlogN

   end subroutine BF_Construct_Element_Compute



!>**** Interface of BP construction via entry extraction
   !> @param BP: the structure containing the BP (inout)
   !> @param option: the structure containing option (in)
   !> @param stats: the structure containing statistics (inout)
   !> @param msh: the structure containing points and ordering information (in)
   !> @param ker: the structure containing kernel quantities (in)
   !> @param ptree: the structure containing process tree (in)
   subroutine BP_Construct_Element_Compute(BP, option, stats, msh, ker, ptree)
      implicit none

      integer Maxlevel
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam
      integer, allocatable:: groupmembers(:)
      integer level

      type(intersect)::submats_dummy(1)
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(blockplus)::BP
      type(matrixblock),pointer::blocks_1,blocks
      type(proctree)::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t0, t1, t2, t3, t4, tt1, tt2, error, Memory, tol_comp_tmp,rtemp1,rtemp2
      integer ierr,pp,knn_tmp,passflag
      integer:: boundary_map(1,1)
      integer groupm_start, Nboundall,Ninadmissible,rankmax
      type(nod), pointer::cur
      class(*), pointer::ptr

      rankmax=0

      if (allocated(msh%xyz)) deallocate (msh%xyz)

      tt1 = MPI_Wtime()
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) " "
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "EntryExtraction-based BP construction......"



      select case (option%format)
      case (HODLR)
         blocks => BP%ll(1)%matrices_block(1)
         groupm_start = 0
         Nboundall = 0
         Ninadmissible = 0
         if (option%forwardN15flag == 1) then
            knn_tmp = option%knn
            option%knn=0
            call BF_compress_N15(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree, 1)
            option%knn=knn_tmp
            call BF_sym2asym(blocks)
         elseif (option%forwardN15flag == 2) then
            call BF_compress_NlogN(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree, 1)
            call BF_checkError(blocks_1, option, msh, ker, stats, ptree, 0, -1, error)
            if(error>50*option%tol_comp)then
               pp = ptree%myid - ptree%pgrp(blocks%pgno)%head + 1
               if(option%verbosity>=0 .and. pp==1)write(*,*)'warning: error ',error,',  with the N15 algorithm'
               call BF_delete(blocks,0)
               knn_tmp = option%knn
               option%tol_comp=option%tol_comp*5
               option%knn=0
               call BF_compress_N15(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree, 1)
               option%knn=knn_tmp
               option%tol_comp=option%tol_comp/5
               call BF_sym2asym(blocks)
            endif
         else
            call BF_compress_NlogN(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree, 1)
         end if
         call BF_ComputeMemory(blocks, stats%Mem_Comp_for)
         if (option%verbosity >= 2)call BF_print_size(blocks)
         rankmax = blocks%rankmax
      case (HMAT,BLR)
         blocks => BP%ll(1)%matrices_block(1)
         Maxlevel = GetTreelevel(msh%Maxgroup) - 1
         allocate (stats%rankmax_of_level(0:Maxlevel))
         stats%rankmax_of_level = 0
         allocate (stats%rankmax_of_level_global(0:Maxlevel))
         stats%rankmax_of_level_global = 0
         T0 = MPI_Wtime()
         do level = 0, Maxlevel
            ! write(*,*)h_mat%lstblks%num_nods,'niam'
            T3 = MPI_Wtime()
            cur => blocks%lstblks(level)%head
            rtemp1 = 0.; rtemp2 = 0.
            do ii = 1, blocks%lstblks(level)%num_nods
               select type (ptr=>cur%item)
               type is (block_ptr)
                  call Hmat_block_construction(ptr%ptr, rtemp1, rtemp2, option, stats, msh, ker, ptree)
               end select
               cur => cur%next
            enddo
            stats%Mem_Comp_for = stats%Mem_Comp_for + rtemp1
            stats%Mem_Direct_for = stats%Mem_Direct_for + rtemp2
            call LogMemory(stats, rtemp1 + rtemp2)

            T4 = MPI_Wtime()
            if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
               write (*, '(A10,I6,A10,Es14.7,A8,Es14.7,A8,Es14.7)') 'Level:', level, 'finished', T4 - T3, 'secnds', rtemp1 + rtemp2, 'Mbytes'
            endif

            passflag = 0
            do while (passflag < 2)
               call element_Zmn_blocklist_user(submats_dummy, 0, msh, option, ker, 2, passflag, ptree, stats)
               ! call element_Zmn_block_user(0, 0, mrange_dummy, nrange_dummy, mat_dummy, msh, option, ker, 2, passflag, ptree, stats)
            enddo
         enddo

         ! !!!!!!! check full error
         if (option%ErrFillFull == 1) call Bplus_CheckError_Full(BP, option, msh, ker, stats, ptree)
         ! !!!!!!! check full error

         T1 = MPI_Wtime()
         stats%Time_Fill = stats%Time_Fill + T1 - T0
         call MPI_ALLREDUCE(stats%rankmax_of_level(0:Maxlevel), stats%rankmax_of_level_global(0:Maxlevel), Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
         rankmax = stats%rankmax_of_level(0)
         Memory = stats%Mem_Comp_for + stats%Mem_Direct_for
      case (HSS)
         call BP_compress_entry(BP, option, Memory, stats, msh, ker, ptree)
         call Bplus_ComputeMemory(BP, stats%Mem_Comp_for, rankmax)
         stats%Mem_Direct_for=0 ! already counted by Bplus_ComputeMemory
      end select


      if (option%verbosity >= 0) call BP_checkError(BP, option, msh, ker, stats, ptree, 0, option%verbosity)


      !>**** delete neighours in msh
      if(allocated(msh%nns))then
         call LogMemory(stats, -SIZEOF(msh%nns)/1024.0d3)
         deallocate(msh%nns)
      endif

      tt2 = MPI_Wtime()

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "EntryExtraction-based BP construction finished in", tt2-tt1, 'Seconds with', Memory,'MB Memory'
      stats%Mem_Fill = stats%Mem_Comp_for + stats%Mem_Direct_for
      if (.not. allocated(stats%rankmax_of_level)) allocate (stats%rankmax_of_level(0:0))
      if (.not. allocated(stats%rankmax_of_level_global)) allocate (stats%rankmax_of_level_global(0:0))
      stats%rankmax_of_level(0) = rankmax
      stats%rankmax_of_level_global(0) = stats%rankmax_of_level(0)
      call LogMemory(stats, stats%Mem_Fill)
      stats%Time_Fill = tt2 - tt1
      ! stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp ! Flop_Fill already counted in BF_compress_NlogN

   end subroutine BP_Construct_Element_Compute





!>**** Interface of BP (generally a rectangular matrix) deletion
   !> @param BP: the structure containing the BP (inout)
   subroutine BP_Delete(BP)
      implicit none
      type(blockplus)::BP
      type(matrixblock),pointer::blocks

      if(associated(BP%ll))then
      if(associated(BP%ll(1)%matrices_block))then
         blocks => BP%ll(1)%matrices_block(1)
         if(blocks%style==4)then !!! H matrix
            call Hmat_block_delete(blocks)
            deallocate(BP%ll(1)%matrices_block)
            deallocate(BP%ll)
         else !!! HSS or HODLR(a single BF)
            call Bplus_delete(BP)
         endif
      endif
      endif
   end subroutine BP_Delete



!>**** Interface of BF construction via entry extraction
   !> @param Ndim: dimensionality (in)
   !> @param blocks: the structure containing the block (inout)
   !> @param option: the structure containing option (in)
   !> @param stats: the structure containing statistics (inout)
   !> @param msh: the structure containing points and ordering information (in)
   !> @param ker: the structure containing kernel quantities (in)
   !> @param ptree: the structure containing process tree (in)
   subroutine BF_MD_Construct_Element_Compute(Ndim, blocks, option, stats, msh, ker, ptree)
      implicit none

      integer Maxlevel,Ndim,dim_i
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam
      integer, allocatable:: groupmembers(:)
      integer level,rmax_start,ntrial

      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh(Ndim)
      type(kernelquant)::ker
      type(matrixblock_MD),target::blocks
      type(matrixblock_MD),pointer::blocks_1
      type(proctree)::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2, error, Memory, tol_comp_tmp
      integer ierr,pp,knn_tmp, level_butterfly_tmp
      integer:: boundary_map(1,1,Ndim)
      integer groupm_start(Ndim), Nboundall,Ninadmissible, format_tmp

      do dim_i =1,Ndim
         if (allocated(msh(dim_i)%xyz)) deallocate (msh(dim_i)%xyz)
      enddo

      t1 = MPI_Wtime()
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) " "
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "EntryExtraction-based BF construction......"
      format_tmp = option%format
      option%format = HTENSOR
      groupm_start = 0
      Nboundall = 0
      Ninadmissible = 0
      if (option%forwardN15flag == 0) then

         call BF_MD_compress_N(Ndim,blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree,1, 0)
         call BF_MD_get_rank(Ndim, blocks, ptree)

         if(blocks%level_butterfly==0)then ! this makes LR compression more robust when rmax is too small
            rmax_start = option%rmax
            ntrial = option%itermax
            do while (ntrial > 1)
               call BF_MD_get_rank(Ndim, blocks, ptree)
               if(blocks%rankmin>=option%rmax)then
                  option%rmax = option%rmax*option%rankrate
                  call BF_MD_delete(Ndim, blocks, 0)
                  call BF_MD_compress_N(Ndim,blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory, stats, msh, ker, ptree,1, 0)
               else
                  exit
               endif
               ntrial = ntrial - 1
            enddo
            option%rmax = rmax_start
         endif

      else
         write(*,*)'forwardN15flag=',option%forwardN15flag,'is not supported in BF_MD_Construct_Element_Compute'
      end if

      ! !!!! the following functions can be expensive as it extracts complete middle-level blocks, currently enabled with option%verbosity >= 2
      blocks_1 => blocks
      if (option%verbosity >= 2)call BF_MD_checkError(Ndim, blocks_1, option, msh, ker, stats, ptree, 0, option%verbosity)
      stats%Mem_Comp_for=Memory
      ! call BF_ComputeMemory(blocks, stats%Mem_Comp_for)


      !>**** delete neighours in msh
      do dim_i =1,Ndim
         if(allocated(msh(dim_i)%nns))then
            call LogMemory(stats, -SIZEOF(msh(dim_i)%nns)/1024.0d3)
            deallocate(msh(dim_i)%nns)
         endif
      enddo

      t2 = MPI_Wtime()

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "EntryExtraction-based BF construction finished in", t2-t1, 'Seconds with', Memory,'MB Memory'

      option%format = format_tmp

      if (.not. allocated(stats%rankmax_of_level)) allocate (stats%rankmax_of_level(0:0))
      if (.not. allocated(stats%rankmax_of_level_global)) allocate (stats%rankmax_of_level_global(0:0))
      stats%rankmax_of_level(0) = blocks%rankmax
      stats%rankmax_of_level_global(0) = stats%rankmax_of_level(0)
      stats%Mem_Fill = stats%Mem_Comp_for + stats%Mem_Direct_for
      call LogMemory(stats, stats%Mem_Fill)
      stats%Time_Fill = stats%Time_Fill + t2 - t1
      ! stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp ! Flop_Fill already counted in BF_compress_NlogN

   end subroutine BF_MD_Construct_Element_Compute




!>**** Computation of the full matrix with matrix entry evaluation
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

!>**** Computation of the construction phase with entry evaluation for the hierarchical tensros
   subroutine BPACK_MD_construction_Element(Ndim, bmat, option, stats, msh, ker, ptree)

      implicit none

      integer Ndim
      type(Hoption)::option
      type(Hstat)::stats
      type(Bmatrix)::bmat
      type(mesh)::msh(Ndim)
      type(kernelquant)::ker
      type(proctree)::ptree
      integer Maxlevel,ii

      do ii=1,Ndim
         if (allocated(msh(ii)%xyz)) deallocate (msh(ii)%xyz)
      enddo
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Hierarchical tensor construction......"

      select case (option%format)
      case (HSS_MD)
         call HSS_MD_construction(Ndim, bmat%hss_bf_md, option, stats, msh, ker, ptree)
      case (HTENSOR)
         call HMAT_MD_construction(Ndim, bmat%h_mat_md, option, stats, msh, ker, ptree)
      case default
         write(*,*)'not supported format in BPACK_MD_construction_Element:', option%format
         stop
      end select

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Hierarchical tensor construction finished"

      if (option%verbosity >= 0) call BPACK_MD_CheckError_SMVP(Ndim, bmat, option, msh, ker, stats, ptree)


   end subroutine BPACK_MD_construction_Element


!>**** Computation of the construction phase with matrix entry evaluation
   subroutine BPACK_construction_Element(bmat, option, stats, msh, ker, ptree)

      implicit none

      type(Hoption)::option
      type(Hstat)::stats
      type(Bmatrix)::bmat
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer Maxlevel

      !write(*,*)stats%Mem_peak,'peak before BPACK_construction_Element'

      call BPACK_gpu_bind(bmat, option, ptree)
      call HODLR_gpu_set_points(bmat, option, msh)  ! (before msh%xyz goes)
      if (allocated(msh%xyz)) deallocate (msh%xyz)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Matrix construction......"

      select case (option%format)
      case (HODLR)
         call HODLR_construction(bmat%ho_bf, option, stats, msh, ker, ptree)
         if (option%HODLR_use_gpu > 0) call HODLR_gpu_upload_forward(bmat%ho_bf, option, stats, ptree)
      case (HMAT,BLR)
         call Hmat_construction(bmat%h_mat, option, stats, msh, ker, ptree)
      case (HSS)
         call HSS_construction(bmat%hss_bf, option, stats, msh, ker, ptree)
      end select

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Matrix construction finished"

      Maxlevel = bmat%Maxlevel
      if (option%lnoBP > Maxlevel .and. option%verbosity >= 0) call BPACK_CheckError_entry(bmat, option, msh, ker, stats, ptree)
      if (option%lnoBP > Maxlevel .and. option%verbosity >= 0) call BPACK_CheckError_SMVP(bmat, option, msh, ker, stats, ptree)

      !write(*,*)stats%Mem_peak,'peak after BPACK_construction_Element'

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
      integer nprow,npcol,myrow,mycol,myArows,myAcols

      call MPI_barrier(ptree%Comm, ierr)

      T0 = MPI_Wtime()

      allocate (stats%rankmax_of_level(0:h_mat%Maxlevel))
      stats%rankmax_of_level = 0
      allocate (stats%rankmax_of_level_global(0:h_mat%Maxlevel))
      stats%rankmax_of_level_global = 0

      scale_factor = 0
      ! compute the largest diagonal entry as the scaling factor
      level=h_mat%Maxlevel
      cur => h_mat%lstblks(level)%head
      do i = 1, h_mat%lstblks(level)%num_nods
         select type (ptr=>cur%item)
         type is (block_ptr)
            if(ptr%ptr%row_group==ptr%ptr%col_group)then
               do ii = ptr%ptr%headm, ptr%ptr%headm + ptr%ptr%M - 1
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
            endif
         end select
         cur => cur%next
      enddo
      if (scale_factor < BPACK_SafeUnderflow) scale_factor = 1d0

      passflag = 0
      do while (passflag == 0)
         call element_Zmn_blocklist_user(submats, 0, msh, option, ker, 2, passflag, ptree, stats)
         ! call element_Zmn_block_user(0, 0, mrange, nrange, mat, msh, option, ker, 2, passflag, ptree, stats)
      enddo


      !!!! It's safer to not modify option%scale_factor inside butterflypack
#if 0
      option%scale_factor = 1d0/scale_factor
      call MPI_ALLREDUCE(MPI_IN_PLACE, option%scale_factor, 1, MPI_DOUBLE_PRECISION, MPI_MIN, ptree%Comm, ierr)

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
         write (*, *) 'element_Zmn is scaled by a factor of:', option%scale_factor
      endif
#endif

      do level = 0, h_mat%Maxlevel
         ! write(*,*)h_mat%lstblks%num_nods,'niam'
         T3 = MPI_Wtime()
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
         call LogMemory(stats, rtemp1 + rtemp2 + rtemp1 + rtemp2)

         T4 = MPI_Wtime()
         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
            write (*, '(A10,I6,A10,Es14.7,A8,Es14.7,A8,Es14.7)') 'Level:', level, 'finished', T4 - T3, 'secnds', rtemp1 + rtemp2, 'Mbytes'
         endif

         passflag = 0
         do while (passflag < 2)
            call element_Zmn_blocklist_user(submats, 0, msh, option, ker, 2, passflag, ptree, stats)
            ! call element_Zmn_block_user(0, 0, mrange_dummy, nrange_dummy, mat_dummy, msh, option, ker, 2, passflag, ptree, stats)
         enddo
      enddo

      if (option%ErrSol == 1 .or. option%precon > 1) then ! no need to save the forward operator if ErrSol=0 and precon=1
         do i = 1, h_mat%myArows
            do j = 1, h_mat%myAcols
               blocks => h_mat%Local_blocks(j, i)
               blocks_copy => h_mat%Local_blocks_copy(j, i)
               call Hmat_block_copy('N', blocks_copy, blocks)
            enddo
         enddo
         call LogMemory(stats, stats%Mem_Comp_for + stats%Mem_Direct_for)
      endif

      T1 = MPI_Wtime()
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
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A26,Es14.2)') 'Total construction flops:', rtemp

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

      implicit none
      real(kind=8) n1, n2, n3, n4, n5,n5_tmp
      integer i, j, ii, ii_inv, jj, kk, iii, jjj, ll
      integer level, blocks, edge, patch, node, group
      integer rank, index_near, m, n, length, flag, itemp, rank0_inner, rank0_outter, ierr
      real T0
      real(kind=8):: rtemp, rel_error, error, t1, t2, tim_tmp, rankrate_inner, rankrate_outter, memory
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
      DT::sym_value
      type(intersect)::submats(1)

      ! Memory_direct_forward=0
      ! Memory_butterfly_forward=0
      call InitStat(stats)
      tim_tmp = 0
      ho_bf1%gpu_host_stale = .false.
      !tolerance=0.001
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      ! write (*,*) 'ACA error threshold',tolerance
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'SVD error threshold', option%tol_comp
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "constructing Leaf-blocks......"

      n1 = MPI_Wtime()
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
         n3 = MPI_Wtime()
         if (HODLR_gpu_construct_ok(ho_bf1, level_c, option, ptree)) then
            ! the whole level on the GPU (HODLR_use_gpu > 0 with a device kernel)
            if (ho_bf1%levels(level_c)%BP(Bidxe)%LL(1)%matrices_block(1)%level /= level) then
               level = ho_bf1%levels(level_c)%BP(Bidxe)%LL(1)%matrices_block(1)%level
               if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'constructing level', level, '(GPU)'
            endif
            if (level_c /= ho_bf1%Maxlevel + 1) then
               call HODLR_gpu_construct_level_checked(ho_bf1, level_c, option, stats, msh, ker, ptree)
            else
               call HODLR_gpu_construct_leaves(ho_bf1, option, stats, msh, ker, ptree)
            endif
         else
         do ii = Bidxs, Bidxe
            if (option%sym > 0 .and. level_c /= ho_bf1%Maxlevel + 1 .and. mod(ii, 2) == 1) cycle
            ! do ii =Bidxs,Bidxs
            if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(ii)%pgno)) then
               if (level_c /= ho_bf1%Maxlevel + 1) then
                  if (ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level /= level) then
                     level = ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level
                     if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'constructing level', level
                  endif

                  ! if(mod(ii,2)==1)then
                  call BP_compress_entry(ho_bf1%levels(level_c)%BP(ii), option, rtemp, stats, msh, ker, ptree)
                  ! else
                  ! call BF_delete(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1),1)
                  ! call BF_copy('T',ho_bf1%levels(level_c)%BP(ii-1)%LL(1)%matrices_block(1),ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1))
                  ! endif

                  if (option%verbosity >= 2)call BF_print_size(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1))

                  if (level == option%level_check) then
                     ! call Bplus_randomized_Exact_test(ho_bf1%levels(level_c)%BP(ii))

                     ! rank0_inner = ho_bf1%levels(level_c)%BP(ii)%LL(2)%rankmax
                     ! rankrate_inner = 1.2d0
                     rank0_outter = ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%rankmax
                     rankrate_outter = option%rankrate !1.2d0
                     level_butterfly = ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%level_butterfly

                     t1 = MPI_Wtime()
                     call BF_randomized(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%pgno, level_butterfly, rank0_outter, rankrate_outter, ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1), ho_bf1%levels(level_c)%BP(ii), Bplus_block_MVP_Exact_dat, error, 'Exact', option, stats, ptree, msh)

                     call BF_ComputeMemory(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1), rtemp)

                     t2 = MPI_Wtime()
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
                  call Full_construction(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1), msh, ker, stats, option, ptree, memory)
                  if (option%sym > 0) then
                     associate (leaf => ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1))
                        call assert(leaf%M == leaf%N, 'a symmetric HODLR dense leaf must be square')
                        do jj = 1, leaf%N
                           do iii = jj + 1, leaf%M
                              sym_value = (leaf%fullmat(iii, jj) + leaf%fullmat(jj, iii))/2d0
                              leaf%fullmat(iii, jj) = sym_value
                              leaf%fullmat(jj, iii) = sym_value
                           enddo
                        enddo
                     end associate
                  endif
                  stats%Mem_Direct_for = stats%Mem_Direct_for + memory
               endif
               ! ! write(*,*)level_c,ii,ho_bf1%levels(level_c)%N_block_forward
               ! if(level>=option%level_check .and. level_c/=ho_bf1%Maxlevel+1)then
               ! call BF_compress_test(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1),msh,ker,element_Zmn,ptree,option,stats)
               ! !stop
               ! end if
            endif
         end do
         endif

         if (option%sym > 0 .and. level_c /= ho_bf1%Maxlevel + 1) then
            do ii = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
               if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP_inverse(ii)%pgno)) then
                  call HODLR_sym_link_pair(ho_bf1%levels(level_c)%BP(2*ii - 1), &
                     ho_bf1%levels(level_c)%BP(2*ii), &
                     ho_bf1%levels(level_c)%BP_inverse(ii)%pgno, ptree, rtemp)
                  stats%Mem_Comp_for = stats%Mem_Comp_for + rtemp
               endif
            enddo
         endif

         ! call MPI_barrier(ptree%Comm, ierr)
         ! do while(option%elem_extract == 1 .and. level==2)

         ! enddo
         passflag = 0
         do while (passflag < 2)
            ! call element_Zmn_block_user(0, 0, mrange_dummy, nrange_dummy, mat_dummy, msh, option, ker, 2, passflag, ptree, stats)
            ! write(*,*)ptree%MyID,level, 'level done'
            call element_Zmn_blocklist_user(submats, 0, msh, option, ker, 2, passflag, ptree, stats)
         enddo

         n4 = MPI_Wtime()
         n5 = n4 - n3
         n5_tmp = n5
         call MPI_ALLREDUCE(n5_tmp, n5, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
         call MPI_ALLREDUCE(stats%rankmax_of_level(0:ho_bf1%Maxlevel), &
            stats%rankmax_of_level_global(0:ho_bf1%Maxlevel), &
            ho_bf1%Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) &
            'time', n5, 'rankmax_of_level so far:', stats%rankmax_of_level_global
      end do
      n2 = MPI_Wtime()
      stats%Time_Fill = stats%Time_Fill + n2 - n1

      stats%Mem_Fill = stats%Mem_Comp_for + stats%Mem_Direct_for
      call LogMemory(stats, stats%Mem_Fill)

      call MPI_ALLREDUCE(stats%rankmax_of_level(0:ho_bf1%Maxlevel), stats%rankmax_of_level_global(0:ho_bf1%Maxlevel), ho_bf1%Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'rankmax_of_level:', stats%rankmax_of_level_global
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      call MPI_ALLREDUCE(stats%Time_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Total construction time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Time_Entry, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Total entry eval time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Flop_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A26,Es14.2)') 'Total construction flops:', rtemp
      call MPI_ALLREDUCE(time_tmp, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp', rtemp
      call MPI_ALLREDUCE(time_tmp1, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp1', rtemp
      call MPI_ALLREDUCE(time_tmp2, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp2', rtemp
      call MPI_ALLREDUCE(time_tmp3, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'time_tmp3', rtemp

      ! if (option%verbosity >= 0) write (*, *) 'time_tmp', time_tmp, 'time_tmp1', time_tmp1, 'time_tmp2', time_tmp2, 'time_tmp3', time_tmp3

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      call MPI_ALLREDUCE(stats%Mem_Comp_for, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) rtemp, 'MB costed for butterfly forward blocks'
      call MPI_ALLREDUCE(stats%Mem_Direct_for, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) rtemp, 'MB costed for direct forward blocks'
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      ! stop

      return

   end subroutine HODLR_construction

   !> Whether level_c of the HODLR construction (Maxlevel+1: the dense
   !> leaves) runs on the GPU: HODLR_use_gpu > 0, a device kernel and the
   !> points registered, and for the low-rank levels BACA
   !> (RecLR_leaf=4) or BACA without overlap (RecLR_leaf=5) on whole blocks
   !> (LR_BLK_NUM=1, forwardN15flag=0)
   logical function HODLR_gpu_construct_ok(ho_bf1, level_c, option, ptree)
      implicit none
      type(hobf)::ho_bf1
      integer level_c
      type(Hoption)::option
      type(proctree)::ptree
      integer(c_int) :: ready


      HODLR_gpu_construct_ok = .false.
      if (option%HODLR_use_gpu <= 0 .or. option%format /= HODLR) return
      if (.not. c_associated(ho_bf1%gpu)) return
      call c_bpack_hodlr_gpu_construct_ready(ho_bf1%gpu, ready)
      if (ready /= 1) return
      if (level_c <= ho_bf1%Maxlevel) then
         if (option%RecLR_leaf /= BACANOVER .and. option%RecLR_leaf /= BACA) return
         if (option%LR_BLK_NUM /= 1 .or. option%forwardN15flag /= 0) return
         if (option%level_check == level_c) return
      else
         if (option%use_zfp /= 0) return
      endif
      HODLR_gpu_construct_ok = .true.
   end function HODLR_gpu_construct_ok

   !> The low-rank blocks of level level_c on the GPU, with BPACK_CHECK=hodlr
   !> compared with the CPU construction of the same blocks from the same
   !> random numbers
   subroutine HODLR_gpu_construct_level_checked(ho_bf1, level_c, option, stats, msh, ker, ptree)
      implicit none
      type(hobf)::ho_bf1
      integer level_c
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer nseed, ierr
      logical differ
      integer, allocatable :: seed0(:), seed1(:), seed2(:)
      real(kind=8) :: flop_gpu

      if (HODLR_gpu_check_level() <= 0) then
         call HODLR_gpu_construct_level(ho_bf1, level_c, option, stats, msh, ker, ptree)
         return
      endif
      call random_seed(size=nseed)
      allocate (seed0(nseed), seed1(nseed), seed2(nseed))
      call random_seed(get=seed0)
      flop_gpu = stats%Flop_Fill
      call HODLR_gpu_construct_level(ho_bf1, level_c, option, stats, msh, ker, ptree)
      flop_gpu = stats%Flop_Fill - flop_gpu
      call random_seed(get=seed1)
      call random_seed(put=seed0)
      call HODLR_gpu_check_construct_level(ho_bf1, level_c, option, msh, ker, ptree, flop_gpu)
      call random_seed(get=seed2)
      differ = any(seed1 /= seed2)
      call MPI_ALLREDUCE(MPI_IN_PLACE, differ, 1, MPI_LOGICAL, MPI_LOR, ptree%Comm, ierr)
      if (differ .and. ptree%MyID == Main_ID) then
         write (*, *) 'HODLR GPU check (construction): the GPU and CPU runs drew different random numbers'
      endif
      call random_seed(put=seed1)
      deallocate (seed0, seed1, seed2)
   end subroutine HODLR_gpu_construct_level_checked

   !> The low-rank blocks of level level_c (1 .. Maxlevel) by BACA without
   !> overlap (RecLR_leaf=5) or BACA (4) on the GPU
   !> (hodlr_gpu/hodlr_construct.hpp).  This is the host side of
   !> LR_BACA_noOverlap / LR_BACA for all blocks of the level at once: the
   !> blocks iterate in lockstep, each through the steps of the CPU routine
   !> in the same order, with the random first columns drawn (rperm) in the
   !> CPU's block order.  The device evaluates the panels and keeps U and
   !> V^T; the host factors the r x r cores, tests convergence from the Gram
   !> matrices of the update, and truncates the SVD of the recompression.
   subroutine HODLR_gpu_construct_level(ho_bf1, level_c, option, stats, msh, ker, ptree)
      implicit none
      type(hobf), target::ho_bf1
      integer level_c
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(matrixblock), pointer::blk
      type(hodlr_gpu_blkptr), allocatable :: blks(:)
      integer, allocatable :: bidx(:), rnk(:)
      integer nb, b, ii, k, npass, nblkit, shared_ii, shared_pg, pgno, Maxgrp
      real(c_double) :: gout(6)
      real(kind=8) :: rtemp, thost, tcore, trec, tbook

      ! the blocks of the level (the symmetric HODLR builds A21 only), each
      ! built by the ranks of its construction group (HODLR_gpu_cgroup); a
      ! block split into pieces goes its own way (one per rank and level): a
      ! block of a group of several ranks, and one of a single rank above the
      ! bottom of the process tree (nproc not a power of 2), whose rank
      ! compresses and merges all the pieces
      ! (the host copies of the factors left unfilled: HODLR_gpu_defer_host, on every rank)
      if (HODLR_gpu_defer_host(option)) ho_bf1%gpu_host_stale = .true.
      Maxgrp = 2**(ptree%nlevel) - 1
      nb = 0
      shared_ii = 0
      shared_pg = 0
      do ii = ho_bf1%levels(level_c)%Bidxs*2 - 1, ho_bf1%levels(level_c)%Bidxe*2
         if (option%sym > 0 .and. mod(ii, 2) == 1) cycle
         pgno = HODLR_gpu_cgroup(ho_bf1, level_c, ii, option)
         if (.not. IOwnPgrp(ptree, pgno)) cycle
         if (ptree%pgrp(pgno)%nproc > 1 .or. pgno*2 <= Maxgrp) then
            call assert(shared_ii == 0, 'HODLR GPU construction: two shared blocks of a level on one rank')
            shared_ii = ii
            shared_pg = pgno
            cycle
         endif
         call assert(IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(ii)%pgno), &
                     'HODLR GPU construction: a local block of another rank')
         nb = nb + 1
      enddo
      if (shared_ii > 0) then
         call assert(nb == 0, 'HODLR GPU construction: shared and local blocks on one level of a rank')
         call HODLR_gpu_construct_shared(ho_bf1, level_c, shared_ii, shared_pg, option, stats, msh, ker, ptree)
         return
      endif
      if (nb == 0) return
      allocate (bidx(nb), blks(nb), rnk(nb))
      b = 0
      do ii = ho_bf1%levels(level_c)%Bidxs*2 - 1, ho_bf1%levels(level_c)%Bidxe*2
         if (option%sym > 0 .and. mod(ii, 2) == 1) cycle
         pgno = HODLR_gpu_cgroup(ho_bf1, level_c, ii, option)
         if (.not. IOwnPgrp(ptree, pgno)) cycle
         if (ptree%pgrp(pgno)%nproc > 1 .or. pgno*2 <= Maxgrp) cycle
         b = b + 1
         bidx(b) = ii
         call assert(ho_bf1%levels(level_c)%BP(ii)%Lplus == 1, 'HODLR GPU construction requires LRlevel=0')
         blks(b)%p => ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)
         call assert(blks(b)%p%level_butterfly == 0, 'HODLR GPU construction: a block is not low rank')
      enddo
      call HODLR_gpu_baca_blocks(ho_bf1%gpu, blks, nb, option, stats, msh, ptree, rnk, gout, thost, tcore, trec, &
                                 npass, nblkit, keep=.true.)

      !>**** the bookkeeping of BF_compress_NlogN and BP_compress_entry
      tbook = MPI_Wtime()
      do b = 1, nb
         blk => blks(b)%p
         k = rnk(b)
         call HODLR_gpu_block_bookkeeping(ho_bf1%levels(level_c)%BP(bidx(b)), blk, k, option, stats, msh, ker, ptree)
      enddo
      tbook = MPI_Wtime() - tbook
      if (ptree%MyID == Main_ID .and. option%verbosity >= 1) then
         write (*, '(A,I3,A,I6,A,4F9.3,A,3F8.3,A,I5,A,I8,A)') ' HODLR GPU construction level', level_c, ':', nb, &
            ' blocks; panels, append, recompression, host', gout(2), gout(3), gout(4), thost, &
            ' s (host: cores, recompression, bookkeeping', tcore, trec, tbook, ' s); iterations', npass, &
            ' (', nblkit, ' block iterations)'
      endif
      deallocate (bidx, blks, rnk)
   end subroutine HODLR_gpu_construct_level

   !> The bookkeeping of BF_compress_NlogN and BP_compress_entry for a block
   !> of rank k built on the GPU (its factors in place)
   subroutine HODLR_gpu_block_bookkeeping(bplus, blk, k, option, stats, msh, ker, ptree)
      implicit none
      type(blockplus)::bplus
      type(matrixblock)::blk
      integer k
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      real(kind=8) :: rtemp

      blk%rankmax = k
      blk%rankmin = k
      blk%ButterflyU%idx = 1
      blk%ButterflyU%inc = 1
      blk%ButterflyU%nblk_loc = 1
      blk%ButterflyV%idx = 1
      blk%ButterflyV%inc = 1
      blk%ButterflyV%nblk_loc = 1
      blk%level_half = BF_Switchlevel(blk%level_butterfly, option%pat_comp)
      bplus%LL(1)%rankmax = k
      if (allocated(stats%rankmax_of_level)) stats%rankmax_of_level(blk%level) = max(k, stats%rankmax_of_level(blk%level))
      call BF_ComputeMemory(blk, rtemp)
      call BF_get_rank(blk, ptree)
      stats%Mem_Comp_for = stats%Mem_Comp_for + rtemp
      if (option%ErrFillFull == 1) call Bplus_CheckError_Full(bplus, option, msh, ker, stats, ptree)
   end subroutine HODLR_gpu_block_bookkeeping

   !> The low-rank block BP(ii) of level level_c, shared by the ranks of its
   !> process group, on the GPU: split among them as LR_HBACA_Leaflevel
   !> does, the pieces of this rank compressed by BACA on the GPU (in the
   !> CPU's order, so they draw the same random columns), merged up the split
   !> tree as LR_HMerge does with the truncated SVDs by TSQR on the GPUs, and
   !> moved into the block's 1D layout.
   !> The group of ranks that builds block ii of level level_c on the GPU:
   !> its own group, or with HODLR_gpu_split_mode 1 in the symmetric HODLR (which
   !> builds A21 only, on half of its node's ranks) its node's group
   integer function HODLR_gpu_cgroup(ho_bf1, level_c, ii, option)
      implicit none
      type(hobf)::ho_bf1
      integer level_c, ii
      type(Hoption)::option

      HODLR_gpu_cgroup = ho_bf1%levels(level_c)%BP(ii)%pgno
      if (option%sym > 0 .and. HODLR_gpu_split_mode() > 0) HODLR_gpu_cgroup = ho_bf1%levels(level_c)%BP_inverse((ii + 1)/2)%pgno
   end function HODLR_gpu_cgroup

   !> Block ii of level level_c, built by the ranks of group cgno (its own
   !> group or its node's, HODLR_gpu_cgroup): split into pieces as in
   !> LR_HBACA_Leaflevel (HODLR_gpu_split), the pieces compressed by BACA on
   !> the GPUs (with HODLR_gpu_split_mode 1 spread over the group by their
   !> estimated cost, HODLR_gpu_balanced_pieces), merged up the split tree
   !> (HODLR_gpu_hmerge), and moved into the block's layout over its group
   subroutine HODLR_gpu_construct_shared(ho_bf1, level_c, ii, cgno, option, stats, msh, ker, ptree)
      implicit none
      type(hobf), target::ho_bf1
      integer level_c, ii
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer cgno
      type(matrixblock), pointer::blk
      type(hodlr_gpu_blkptr), allocatable :: leaves(:)
      integer, allocatable :: up(:, :), vp(:, :), rnk(:), mpo(:, :), npo(:, :), pid(:, :)
      integer nleaf, rank, npass, nblkit, ierr, bpg, nextra, nbp
      integer(c_int) :: uid, vid, bu, bv, rows_c, cols_c, stale_c
      real(c_double) :: fmerge
      logical inblk, defer
      real(c_double) :: gout(6)
      real(kind=8) :: thost, tcore, trec, t0, t1, tb, tmg, tp(2)

      t0 = MPI_Wtime()
      hodlr_gpu_merge_time = 0
      hodlr_gpu_tsqr_time = 0
      blk => ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)
      call assert(ho_bf1%levels(level_c)%BP(ii)%Lplus == 1, 'HODLR GPU construction requires LRlevel=0')
      bpg = blk%pgno
      inblk = IOwnPgrp(ptree, bpg)
      nextra = -1
      ! (the option HODLR_gpu_pieces, rounded down to a power of 2: that many more cuts)
      if (HODLR_gpu_split_mode() > 0 .and. option%HODLR_gpu_pieces > 1) &
         nextra = floor(log(dble(option%HODLR_gpu_pieces))/log(2d0) + 1d-9)
      allocate (leaves(256))
      nleaf = 0
      call HODLR_gpu_split(blk, cgno, 0, nextra, ptree, leaves, nleaf)
      blk%pgno = bpg
      allocate (rnk(max(1, nleaf)), pid(3, max(1, nleaf)))
      pid = 0
      if (nextra >= 0 .and. ptree%pgrp(cgno)%nproc > 1) then
         call HODLR_gpu_balanced_pieces(ho_bf1%gpu, blk, cgno, nextra, leaves, nleaf, option, stats, msh, ptree, gout, pid)
      else
         call HODLR_gpu_baca_blocks(ho_bf1%gpu, leaves, nleaf, option, stats, msh, ptree, rnk, gout, thost, tcore, trec, &
                                    npass, nblkit, dids=pid)
      endif
      tb = MPI_Wtime()
      tp(1) = tb - t0
      tp(2) = -(tb - t0)
      call MPI_ALLREDUCE(MPI_IN_PLACE, tp, 2, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%pgrp(cgno)%Comm, ierr)
      tb = MPI_Wtime()
      call HODLR_gpu_hmerge(ho_bf1%gpu, blk, cgno, 0, option, ptree, leaves, nleaf, pid, uid, vid, up, vp, rank)
      tmg = MPI_Wtime()
      ! into the block's layout (M_p, N_p over its group, known on the ranks of cgno)
      nbp = ptree%pgrp(bpg)%nproc
      allocate (mpo(nbp, 2), npo(nbp, 2))
      mpo = 0
      npo = 0
      if (ptree%MyID == ptree%pgrp(bpg)%head) then
         mpo = blk%M_p
         npo = blk%N_p
      endif
      if (cgno /= bpg) then
         call MPI_ALLREDUCE(MPI_IN_PLACE, mpo, 2*nbp, MPI_INTEGER, MPI_MAX, ptree%pgrp(cgno)%Comm, ierr)
         call MPI_ALLREDUCE(MPI_IN_PLACE, npo, 2*nbp, MPI_INTEGER, MPI_MAX, ptree%pgrp(cgno)%Comm, ierr)
      else
         mpo = blk%M_p
         npo = blk%N_p
      endif
      ! (on the device, then this rank's rows to the host, unless the host copies are deferred)
      defer = HODLR_gpu_defer_host(option)
      stale_c = 0
      if (defer) stale_c = 1
      bu = 0
      bv = 0
      if (inblk) then
         allocate (blk%ButterflyU%blocks(1), blk%ButterflyV%blocks(1))
         allocate (blk%ButterflyU%blocks(1)%matrix(max(1, blk%M_loc), rank), &
                   blk%ButterflyV%blocks(1)%matrix(max(1, blk%N_loc), rank))
         if (.not. defer .or. blk%M_loc == 0) blk%ButterflyU%blocks(1)%matrix = 0
         if (.not. defer .or. blk%N_loc == 0) blk%ButterflyV%blocks(1)%matrix = 0
         rows_c = blk%M_loc
         cols_c = rank
         call c_bpack_gpu_dm_alloc(ho_bf1%gpu, rows_c, cols_c, c_null_ptr, bu)
         rows_c = blk%N_loc
         call c_bpack_gpu_dm_alloc(ho_bf1%gpu, rows_c, cols_c, c_null_ptr, bv)
      endif
      call HODLR_gpu_dm_redist(ho_bf1%gpu, uid, up, cgno, bu, 0, mpo, bpg, rank, ptree)
      call HODLR_gpu_dm_redist(ho_bf1%gpu, vid, vp, cgno, bv, 0, npo, bpg, rank, ptree)
      call c_bpack_gpu_dm_free(ho_bf1%gpu, uid)
      call c_bpack_gpu_dm_free(ho_bf1%gpu, vid)
      if (inblk) then
         if (.not. defer) then
            if (blk%M_loc > 0) call c_bpack_gpu_dm_download(ho_bf1%gpu, bu, c_loc(blk%ButterflyU%blocks(1)%matrix(1, 1)))
            if (blk%N_loc > 0) call c_bpack_gpu_dm_download(ho_bf1%gpu, bv, c_loc(blk%ButterflyV%blocks(1)%matrix(1, 1)))
         endif
         ! (U and V stay on the device for HODLR_gpu_upload_forward, which moves V to its column child's layout)
         if (blk%M_loc > 0) then
            call c_bpack_gpu_dm_keep(ho_bf1%gpu, bu, c_loc(blk%ButterflyU%blocks(1)%matrix(1, 1)), stale_c)
         else
            call c_bpack_gpu_dm_free(ho_bf1%gpu, bu)
         endif
         if (blk%N_loc > 0) then
            call c_bpack_gpu_dm_keep(ho_bf1%gpu, bv, c_loc(blk%ButterflyV%blocks(1)%matrix(1, 1)), stale_c)
         else
            call c_bpack_gpu_dm_free(ho_bf1%gpu, bv)
         endif
      endif
      deallocate (up, vp, leaves, rnk, pid, mpo, npo)
      ! (the flops of the merges: TSQR and the products of the other side)
      call c_bpack_gpu_dm_flops(ho_bf1%gpu, fmerge)
      stats%Flop_Fill = stats%Flop_Fill + fmerge
      if (inblk) call HODLR_gpu_block_bookkeeping(ho_bf1%levels(level_c)%BP(ii), blk, rank, option, stats, msh, ker, ptree)
      t1 = MPI_Wtime()
      if (ptree%MyID == ptree%pgrp(cgno)%head .and. option%verbosity >= 1) then
         write (*, '(A,I3,A,I5,A,I3,A,I6,A,F9.3,A,2F8.3,A,F8.3,A,F8.3,A)') ' HODLR GPU construction level', level_c, &
            ': block', ii, ' over', ptree%pgrp(cgno)%nproc, ' ranks, rank', rank, ',', t1 - t0, &
            ' s (pieces min, max', -tp(2), tp(1), ', merges', tmg - tb, ', layout', t1 - tmg, ')'
         write (*, '(A,3F8.3,A,5F8.3)') '    merges: redistribute, TSQR, other side', hodlr_gpu_merge_time, &
            '; TSQR: local QR, up, root SVD, down, apply', hodlr_gpu_tsqr_time
      endif
   end subroutine HODLR_gpu_construct_shared

   !> The split of a block among the ranks of group pgno as in
   !> LR_HBACA_Leaflevel (LR_BLK_NUM=1), alternating row and column cuts in
   !> proportion to the ranks of the two child groups: the sons of the nodes
   !> of this rank, and its pieces in depth-first order.  nextra < 0: the
   !> CPU's pieces (at even cridx on single ranks); nextra >= 0: each rank's
   !> part cut nextra more times (HODLR_gpu_split_leaf)
   recursive subroutine HODLR_gpu_split(node, pgno, cridx, nextra, ptree, leaves, nleaf)
      implicit none
      type(matrixblock), target::node
      integer pgno, cridx, nextra, nleaf
      type(proctree)::ptree
      type(hodlr_gpu_blkptr) :: leaves(:)
      integer pgno1, pgno2, cut, next
      logical col

      node%pgno = pgno
      if (HODLR_gpu_split_leaf(node%M, node%N, pgno, cridx, nextra, ptree)) then
         nleaf = nleaf + 1
         call assert(nleaf <= size(leaves), 'HODLR_gpu_split: too many pieces on one rank')
         leaves(nleaf)%p => node
         return
      endif
      call HODLR_gpu_split_cut(node%M, node%N, pgno, cridx, nextra, ptree, pgno1, pgno2, col, cut, next)
      allocate (node%sons(2, 1))
      node%sons(1, 1)%headm = node%headm
      node%sons(1, 1)%headn = node%headn
      node%sons(2, 1)%headm = node%headm
      node%sons(2, 1)%headn = node%headn
      if (col) then  ! split along the column dimension
         node%sons(1, 1)%M = node%M
         node%sons(1, 1)%N = cut
         node%sons(2, 1)%M = node%M
         node%sons(2, 1)%headn = node%headn + cut
         node%sons(2, 1)%N = node%N - cut
      else  ! split along the row dimension
         node%sons(1, 1)%N = node%N
         node%sons(1, 1)%M = cut
         node%sons(2, 1)%N = node%N
         node%sons(2, 1)%headm = node%headm + cut
         node%sons(2, 1)%M = node%M - cut
      endif
      if (IOwnPgrp(ptree, pgno1)) call HODLR_gpu_split(node%sons(1, 1), pgno1, cridx + 1, next, ptree, leaves, nleaf)
      if (IOwnPgrp(ptree, pgno2)) call HODLR_gpu_split(node%sons(2, 1), pgno2, cridx + 1, next, ptree, leaves, nleaf)
   end subroutine HODLR_gpu_split

   !> Whether a node (m x n, group pgno, depth cridx) of HODLR_gpu_split is a
   !> piece: on a single rank at the bottom of the process tree, at even
   !> cridx, and (nextra >= 0) after nextra more cuts
   logical function HODLR_gpu_split_leaf(m, n, pgno, cridx, nextra, ptree)
      implicit none
      integer m, n, pgno, cridx, nextra
      type(proctree)::ptree
      integer Maxgrp

      ! (a piece at even cridx: BACA leaves the singular values in U, and the
      ! merges alternate from a column merge, which decomposes [U1 U2],
      ! above it; a row merge above a piece would truncate its orthonormal V)
      Maxgrp = 2**(ptree%nlevel) - 1
      if (pgno*2 <= Maxgrp .or. mod(cridx, 2) /= 0) then
         HODLR_gpu_split_leaf = .false.
      elseif (nextra <= 0) then
         HODLR_gpu_split_leaf = .true.
      else
         HODLR_gpu_split_leaf = m < 2048 .or. n < 2048  ! (small pieces: not worth balancing)
      endif
   end function HODLR_gpu_split_leaf

   !> The cut of a node (m x n, group pgno, depth cridx) of HODLR_gpu_split:
   !> the groups of its sons, along the columns (col) or the rows, the size
   !> of son 1 there, and the sons' nextra
   subroutine HODLR_gpu_split_cut(m, n, pgno, cridx, nextra, ptree, pgno1, pgno2, col, cut, next)
      implicit none
      integer m, n, pgno, cridx, nextra, pgno1, pgno2, cut, next
      type(proctree)::ptree
      logical col
      integer Maxgrp, nsproc1, nsproc2

      Maxgrp = 2**(ptree%nlevel) - 1
      next = nextra
      if (pgno*2 > Maxgrp) then
         pgno1 = pgno
         pgno2 = pgno
         if (nextra > 0) next = nextra - 1
      else
         pgno1 = pgno*2
         pgno2 = pgno*2 + 1
      endif
      nsproc1 = ptree%pgrp(pgno1)%nproc
      nsproc2 = ptree%pgrp(pgno2)%nproc
      col = mod(cridx + 1, 2) == 0
      if (col) then
         cut = INT(n*dble(nsproc1)/(dble(nsproc1 + nsproc2)))
         call assert(cut > 0 .and. cut < n, 'column of blocks%sons(1,1) or blocks%sons(2,1) cannot be empty')
      else
         cut = INT(m*dble(nsproc1)/(dble(nsproc1 + nsproc2)))
         call assert(cut > 0 .and. cut < m, 'row of blocks%sons(1,1) or blocks%sons(2,1) cannot be empty')
      endif
   end subroutine HODLR_gpu_split_cut

   !> All pieces of HODLR_gpu_split (on every rank, in the same order): pl(:,
   !> k) = the first row, rows, first column, columns and owner (MyID) of piece k
   recursive subroutine HODLR_gpu_piece_list(headm, m, headn, n, pgno, cridx, nextra, ptree, pl, np)
      implicit none
      integer headm, m, headn, n, pgno, cridx, nextra, np
      type(proctree)::ptree
      integer pl(:, :)
      integer pgno1, pgno2, cut, next
      logical col

      if (HODLR_gpu_split_leaf(m, n, pgno, cridx, nextra, ptree)) then
         np = np + 1
         call assert(np <= size(pl, 2), 'HODLR_gpu_piece_list: too many pieces')
         pl(:, np) = (/headm, m, headn, n, ptree%pgrp(pgno)%head/)
         return
      endif
      call HODLR_gpu_split_cut(m, n, pgno, cridx, nextra, ptree, pgno1, pgno2, col, cut, next)
      if (col) then
         call HODLR_gpu_piece_list(headm, m, headn, cut, pgno1, cridx + 1, next, ptree, pl, np)
         call HODLR_gpu_piece_list(headm, m, headn + cut, n - cut, pgno2, cridx + 1, next, ptree, pl, np)
      else
         call HODLR_gpu_piece_list(headm, cut, headn, n, pgno1, cridx + 1, next, ptree, pl, np)
         call HODLR_gpu_piece_list(headm + cut, m - cut, headn, n, pgno2, cridx + 1, next, ptree, pl, np)
      endif
   end subroutine HODLR_gpu_piece_list

   !> The estimated cost of compressing a piece (rows headm .. headm+m-1,
   !> columns headn .. headn+n-1): (m + n) r^2 with its rank r taken to fall
   !> with the gap between the bounding boxes of its row and column points
   !> relative to their size (m + n without the points)
   real(kind=8) function HODLR_gpu_piece_cost(msh, headm, m, headn, n)
      implicit none
      type(mesh)::msh
      integer headm, m, headn, n
      integer i, d, nd
      real(kind=8) :: lo(3, 2), hi(3, 2), gap, diam, x

      HODLR_gpu_piece_cost = dble(m + n)
      if (.not. allocated(msh%xyz) .or. .not. allocated(msh%new2old)) return
      nd = size(msh%xyz, 1)
      if (nd < 1 .or. nd > 3) return
      lo = huge(1d0)
      hi = -huge(1d0)
      do i = headm, headm + m - 1
         do d = 1, nd
            x = msh%xyz(d, msh%new2old(i))
            lo(d, 1) = min(lo(d, 1), x)
            hi(d, 1) = max(hi(d, 1), x)
         enddo
      enddo
      do i = headn, headn + n - 1
         do d = 1, nd
            x = msh%xyz(d, msh%new2old(i))
            lo(d, 2) = min(lo(d, 2), x)
            hi(d, 2) = max(hi(d, 2), x)
         enddo
      enddo
      gap = 0
      do d = 1, nd
         gap = gap + max(0d0, lo(d, 1) - hi(d, 2), lo(d, 2) - hi(d, 1))**2
      enddo
      gap = sqrt(gap)
      diam = max(sqrt(sum((hi(1:nd, 1) - lo(1:nd, 1))**2)), sqrt(sum((hi(1:nd, 2) - lo(1:nd, 2))**2)), 1d-300)
      HODLR_gpu_piece_cost = dble(m + n)/(1d0 + 4d0*gap/diam)**2
   end function HODLR_gpu_piece_cost

   !> The pieces of a block split over the ranks of group cgno
   !> (HODLR_gpu_split with nextra >= 0), compressed by BACA on the GPUs of
   !> the group in a static schedule: the pieces, largest estimated cost
   !> first (HODLR_gpu_piece_cost), each to the rank of least estimated load
   !> (its owner if nearly as idle), the random first columns seeded by the
   !> piece (the same whichever rank computes it); the pieces computed off
   !> their owner then move to it.  leaves: this rank's pieces; out: pid(:, i)
   !> the ids of the device factors U, V of leaves(i) and its rank
   subroutine HODLR_gpu_balanced_pieces(gpu, blk, cgno, nextra, leaves, nleaf, option, stats, msh, ptree, gout, pid)
      implicit none
      type(c_ptr) :: gpu
      type(matrixblock)::blk
      integer cgno, nextra, nleaf
      type(hodlr_gpu_blkptr) :: leaves(:)
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(proctree)::ptree
      real(c_double) :: gout(6)
      integer :: pid(:, :)
      integer, allocatable :: pl(:, :), asg(:), order(:), rk(:), seeds(:), mine(:), rnk(:), sreq(:), rreq(:), leafof(:)
      integer, allocatable :: tid(:, :)
      integer(c_int) :: rows_c, cols_c, idu, idv
      real(kind=8), allocatable :: cost(:), load(:)
      type(matrixblock), allocatable, target :: tmp(:)
      type(hodlr_gpu_blkptr), allocatable :: list(:)
      integer np, nproc, head, k, q, i, j, b, nmine, ierr, sb, comm, crank, nsr, nrr, own, npass, nblkit
      integer, parameter :: tag0 = 7400
      real(kind=8) :: x, thost, tcore, trec

      comm = ptree%pgrp(cgno)%Comm
      nproc = ptree%pgrp(cgno)%nproc
      head = ptree%pgrp(cgno)%head
      crank = ptree%MyID - head
      ! all pieces, the same list and schedule on every rank
      allocate (pl(5, 4096))
      np = 0
      call HODLR_gpu_piece_list(blk%headm, blk%M, blk%headn, blk%N, cgno, 0, nextra, ptree, pl, np)
      allocate (cost(np), order(np), asg(np), load(0:nproc - 1), leafof(np), rk(np))
      do k = 1, np
         cost(k) = HODLR_gpu_piece_cost(msh, pl(1, k), pl(2, k), pl(3, k), pl(4, k))
         order(k) = k
      enddo
      do i = 2, np  ! by decreasing cost (insertion sort, stable)
         k = order(i)
         j = i - 1
         do while (j >= 1)
            if (cost(order(j)) >= cost(k)) exit
            order(j + 1) = order(j)
            j = j - 1
         enddo
         order(j + 1) = k
      enddo
      load = 0
      do i = 1, np
         k = order(i)
         own = pl(5, k) - head
         q = minloc(load, 1) - 1
         if (load(own) <= load(q) + 0.1d0*cost(k)) q = own
         asg(k) = q
         load(q) = load(q) + cost(k)
      enddo
      if (crank == 0) then
         call random_number(x)
         sb = int(x*50000d0)
      endif
      call MPI_BCAST(sb, 1, MPI_INTEGER, 0, comm, ierr)
      leafof = 0
      do b = 1, nleaf
         do k = 1, np
            if (pl(1, k) == leaves(b)%p%headm .and. pl(2, k) == leaves(b)%p%M .and. pl(3, k) == leaves(b)%p%headn &
                .and. pl(4, k) == leaves(b)%p%N .and. pl(5, k) == ptree%MyID) leafof(k) = b
         enddo
      enddo
      call assert(count(leafof > 0) == nleaf, 'HODLR_gpu_balanced_pieces: pieces and leaves differ')
      ! this rank's pieces, largest first: its leaves, or temporaries for the pieces of other ranks
      nmine = count(asg == crank)
      allocate (list(max(1, nmine)), seeds(max(1, nmine)), mine(max(1, nmine)), rnk(max(1, nmine)), tmp(max(1, nmine)))
      nmine = 0
      do i = 1, np
         k = order(i)
         if (asg(k) /= crank) cycle
         nmine = nmine + 1
         mine(nmine) = k
         seeds(nmine) = sb + k
         if (leafof(k) > 0) then
            list(nmine)%p => leaves(leafof(k))%p
         else
            tmp(nmine)%headm = pl(1, k)
            tmp(nmine)%M = pl(2, k)
            tmp(nmine)%headn = pl(3, k)
            tmp(nmine)%N = pl(4, k)
            tmp(nmine)%level_butterfly = 0
            list(nmine)%p => tmp(nmine)
         endif
      enddo
      gout = 0
      allocate (tid(3, max(1, nmine)))
      tid = 0
      if (nmine > 0) call HODLR_gpu_baca_blocks(gpu, list, nmine, option, stats, msh, ptree, rnk, gout, thost, tcore, &
                                                trec, npass, nblkit, seeds, dids=tid)
      ! the ranks of all pieces; this rank's own pieces stay on the device, the others move to their
      ! owner (through the host)
      rk = 0
      do b = 1, nmine
         rk(mine(b)) = rnk(b)
         if (leafof(mine(b)) > 0) pid(:, leafof(mine(b))) = tid(:, b)
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, rk, np, MPI_INTEGER, MPI_MAX, comm, ierr)
      allocate (sreq(2*np), rreq(2*np))
      nsr = 0
      nrr = 0
      do b = 1, nmine
         k = mine(b)
         own = pl(5, k) - head
         if (own == crank) cycle
         allocate (tmp(b)%ButterflyU%blocks(1), tmp(b)%ButterflyV%blocks(1))
         allocate (tmp(b)%ButterflyU%blocks(1)%matrix(pl(2, k), rk(k)), tmp(b)%ButterflyV%blocks(1)%matrix(pl(4, k), rk(k)))
         idu = tid(1, b)
         idv = tid(2, b)
         call c_bpack_gpu_dm_download(gpu, idu, c_loc(tmp(b)%ButterflyU%blocks(1)%matrix(1, 1)))
         call c_bpack_gpu_dm_download(gpu, idv, c_loc(tmp(b)%ButterflyV%blocks(1)%matrix(1, 1)))
         call c_bpack_gpu_dm_free(gpu, idu)
         call c_bpack_gpu_dm_free(gpu, idv)
         nsr = nsr + 1
         call MPI_ISEND(tmp(b)%ButterflyU%blocks(1)%matrix, pl(2, k)*rk(k), MPI_DT, own, tag0 + 2*k, comm, sreq(nsr), ierr)
         nsr = nsr + 1
         call MPI_ISEND(tmp(b)%ButterflyV%blocks(1)%matrix, pl(4, k)*rk(k), MPI_DT, own, tag0 + 2*k + 1, comm, sreq(nsr), ierr)
      enddo
      do k = 1, np
         b = leafof(k)
         if (b == 0 .or. asg(k) == crank) cycle
         allocate (leaves(b)%p%ButterflyU%blocks(1), leaves(b)%p%ButterflyV%blocks(1))
         allocate (leaves(b)%p%ButterflyU%blocks(1)%matrix(pl(2, k), rk(k)), leaves(b)%p%ButterflyV%blocks(1)%matrix(pl(4, k), rk(k)))
         nrr = nrr + 1
         call MPI_IRECV(leaves(b)%p%ButterflyU%blocks(1)%matrix, pl(2, k)*rk(k), MPI_DT, asg(k), tag0 + 2*k, comm, &
                        rreq(nrr), ierr)
         nrr = nrr + 1
         call MPI_IRECV(leaves(b)%p%ButterflyV%blocks(1)%matrix, pl(4, k)*rk(k), MPI_DT, asg(k), tag0 + 2*k + 1, comm, &
                        rreq(nrr), ierr)
      enddo
      if (nrr > 0) call MPI_WAITALL(nrr, rreq, MPI_STATUSES_IGNORE, ierr)
      if (nsr > 0) call MPI_WAITALL(nsr, sreq, MPI_STATUSES_IGNORE, ierr)
      do k = 1, np  ! the pieces received: to the device
         b = leafof(k)
         if (b == 0 .or. asg(k) == crank) cycle
         cols_c = rk(k)
         rows_c = pl(2, k)
         call c_bpack_gpu_dm_alloc(gpu, rows_c, cols_c, c_loc(leaves(b)%p%ButterflyU%blocks(1)%matrix(1, 1)), idu)
         rows_c = pl(4, k)
         call c_bpack_gpu_dm_alloc(gpu, rows_c, cols_c, c_loc(leaves(b)%p%ButterflyV%blocks(1)%matrix(1, 1)), idv)
         pid(:, b) = (/idu, idv, rk(k)/)
         deallocate (leaves(b)%p%ButterflyU%blocks(1)%matrix, leaves(b)%p%ButterflyV%blocks(1)%matrix)
         deallocate (leaves(b)%p%ButterflyU%blocks, leaves(b)%p%ButterflyV%blocks)
      enddo
      do b = 1, nmine
         if (leafof(mine(b)) > 0) cycle
         deallocate (tmp(b)%ButterflyU%blocks(1)%matrix, tmp(b)%ButterflyV%blocks(1)%matrix)
         deallocate (tmp(b)%ButterflyU%blocks, tmp(b)%ButterflyV%blocks)
      enddo
      if (option%verbosity >= 2 .and. crank == 0) then
         write (*, '(A,I4,A,I3,A,2ES10.3)') '    balanced pieces:', np, ' over', nproc, &
            ' ranks; estimated load min, max', minval(load), maxval(load)
      endif
      deallocate (pl, cost, order, asg, load, leafof, rk, list, seeds, mine, rnk, tmp, sreq, rreq, tid)
   end subroutine HODLR_gpu_balanced_pieces


   !> The merge of the pieces of a shared block up the split tree
   !> (LR_HMerge): for a column split, [U1 U2] = UU S VV by a truncated SVD
   !> and V = [V1 VV1^T; V2 VV2^T] with S in VV; for a row split the same
   !> with the roles of U and V swapped.  The factors stay in device memory
   !> (ids of DistQr::dm_*, hodlr_gpu/hodlr_distqr.hpp): the pieces' factors
   !> go to the device once, the decomposed side moves between the GPUs into
   !> an even layout over the node's ranks (HODLR_gpu_dm_redist), the TSQR
   !> runs there (c_bpack_gpu_dm_tsqr), and the other side's products too.
   !> The pieces (leaves, nleaf) are on the device already: pid(:, i) the ids of
   !> U, V of leaves(i) and its rank (HODLR_gpu_baca_blocks).  Out: the ids
   !> uid, vid of this rank's rows of the node's factors, their layouts up, vp
   !> over the node's group (rows of the node, from 1), the rank
   recursive subroutine HODLR_gpu_hmerge(gpu, node, pgno, cridx, option, ptree, leaves, nleaf, pid, uid, vid, up, vp, rank)
      implicit none
      type(c_ptr) :: gpu
      type(matrixblock), target::node
      integer pgno, cridx, rank
      type(Hoption)::option
      type(proctree)::ptree
      integer(c_int) :: uid, vid
      integer, allocatable :: up(:, :), vp(:, :)
      type(hodlr_gpu_blkptr) :: leaves(:)
      integer nleaf, pid(:, :)
      DT, allocatable, target :: sw(:, :), swp(:, :)
      integer, allocatable :: up1(:, :), vp1(:, :), up2(:, :), vp2(:, :), dg(:, :), op(:, :), lay(:)
      integer, allocatable :: xp1(:, :), xp2(:, :), op1(:, :), op2(:, :)
      integer pgno1, pgno2, nproc, nproc1, nproc2, Maxgrp, dims(6), p, pp, gid, p1, p2, ierr, l, l_loc, r1, r2
      integer m1, n1, m2, n2, o1, o2, o1_loc, o2_loc, lo, hi, comm, i
      logical col, in1, in2
      integer(c_int) :: u1, v1, u2, v2, xid, xnid, oid, rows_c, cols_c, m_c, k_c, n_c, c0_c, rn_c, a_c
      real(c_double) :: tol_c, uf_c, tsec(5)
      real(kind=8) :: tm0, tm1

      if (.not. associated(node%sons)) then  ! a piece compressed by this rank alone (on the device)
         p = 0
         do i = 1, nleaf
            if (associated(leaves(i)%p, node)) p = i
         enddo
         call assert(p > 0, 'HODLR_gpu_hmerge: a piece not in the list')
         call assert(pid(1, p) /= 0 .and. pid(2, p) /= 0, 'HODLR_gpu_hmerge: a piece without its device factors')
         uid = pid(1, p)
         vid = pid(2, p)
         rank = pid(3, p)
         allocate (up(1, 2), vp(1, 2))
         up(1, :) = (/1, node%M/)
         vp(1, :) = (/1, node%N/)
         return
      endif
      Maxgrp = 2**(ptree%nlevel) - 1
      if (pgno*2 > Maxgrp) then
         pgno1 = pgno
         pgno2 = pgno
      else
         pgno1 = pgno*2
         pgno2 = pgno*2 + 1
      endif
      nproc = ptree%pgrp(pgno)%nproc
      nproc1 = ptree%pgrp(pgno1)%nproc
      nproc2 = ptree%pgrp(pgno2)%nproc
      in1 = IOwnPgrp(ptree, pgno1)
      in2 = IOwnPgrp(ptree, pgno2)
      r1 = 0
      r2 = 0
      u1 = 0
      v1 = 0
      u2 = 0
      v2 = 0
      if (in1) call HODLR_gpu_hmerge(gpu, node%sons(1, 1), pgno1, cridx + 1, option, ptree, leaves, nleaf, pid, u1, v1, up1, &
                                     vp1, r1)
      if (in2) call HODLR_gpu_hmerge(gpu, node%sons(2, 1), pgno2, cridx + 1, option, ptree, leaves, nleaf, pid, u2, v2, up2, &
                                     vp2, r2)
      ! sizes, ranks and the sons' layouts on every rank of the node
      dims = 0
      if (in1) dims(1:3) = (/node%sons(1, 1)%M, node%sons(1, 1)%N, r1/)
      if (in2) dims(4:6) = (/node%sons(2, 1)%M, node%sons(2, 1)%N, r2/)
      call MPI_ALLREDUCE(MPI_IN_PLACE, dims, 6, MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      m1 = dims(1)
      n1 = dims(2)
      r1 = dims(3)
      m2 = dims(4)
      n2 = dims(5)
      r2 = dims(6)
      allocate (lay(4*nproc1 + 4*nproc2))
      lay = 0
      if (in1 .and. ptree%MyID == ptree%pgrp(pgno1)%head) then
         lay(1:2*nproc1) = reshape(up1, (/2*nproc1/))
         lay(2*nproc1 + 1:4*nproc1) = reshape(vp1, (/2*nproc1/))
      endif
      if (in2 .and. ptree%MyID == ptree%pgrp(pgno2)%head) then
         lay(4*nproc1 + 1:4*nproc1 + 2*nproc2) = reshape(up2, (/2*nproc2/))
         lay(4*nproc1 + 2*nproc2 + 1:4*nproc1 + 4*nproc2) = reshape(vp2, (/2*nproc2/))
      endif
      call MPI_ALLREDUCE(MPI_IN_PLACE, lay, size(lay), MPI_INTEGER, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)
      col = mod(cridx + 1, 2) == 0
      allocate (xp1(nproc1, 2), op1(nproc1, 2), xp2(nproc2, 2), op2(nproc2, 2))
      if (col) then  ! the sons share the rows: [U1 U2] is decomposed, V1 and V2 updated
         xp1 = reshape(lay(1:2*nproc1), (/nproc1, 2/))
         op1 = reshape(lay(2*nproc1 + 1:4*nproc1), (/nproc1, 2/))
         xp2 = reshape(lay(4*nproc1 + 1:4*nproc1 + 2*nproc2), (/nproc2, 2/))
         op2 = reshape(lay(4*nproc1 + 2*nproc2 + 1:4*nproc1 + 4*nproc2), (/nproc2, 2/))
         call assert(m1 == m2, 'M1/=M2 in column merge (GPU)')
         l = m1
         o1 = n1
         o2 = n2
      else  ! the sons share the columns: [V1 V2] is decomposed, U1 and U2 updated
         op1 = reshape(lay(1:2*nproc1), (/nproc1, 2/))
         xp1 = reshape(lay(2*nproc1 + 1:4*nproc1), (/nproc1, 2/))
         op2 = reshape(lay(4*nproc1 + 1:4*nproc1 + 2*nproc2), (/nproc2, 2/))
         xp2 = reshape(lay(4*nproc1 + 2*nproc2 + 1:4*nproc1 + 4*nproc2), (/nproc2, 2/))
         call assert(n1 == n2, 'N1/=N2 in row merge (GPU)')
         l = n1
         o1 = m1
         o2 = m2
      endif
      ! the decomposed side in an even split of its l rows over the node's ranks
      allocate (dg(nproc, 2))
      do p = 1, nproc
         dg(p, 1) = int((int(p - 1, 8)*l)/nproc) + 1
         dg(p, 2) = int((int(p, 8)*l)/nproc)
      enddo
      pp = ptree%MyID - ptree%pgrp(pgno)%head + 1
      l_loc = dg(pp, 2) - dg(pp, 1) + 1
      tm0 = MPI_Wtime()
      rows_c = l_loc
      cols_c = r1 + r2
      call c_bpack_gpu_dm_alloc(gpu, rows_c, cols_c, c_null_ptr, xid)
      call HODLR_gpu_dm_redist(gpu, merge(u1, v1, col), xp1, pgno1, xid, 0, dg, pgno, r1, ptree)
      call HODLR_gpu_dm_redist(gpu, merge(u2, v2, col), xp2, pgno2, xid, r1, dg, pgno, r2, ptree)
      tm1 = MPI_Wtime()
      hodlr_gpu_merge_time(1) = hodlr_gpu_merge_time(1) + tm1 - tm0
      comm = ptree%pgrp(pgno)%Comm
      tol_c = option%tol_comp
      uf_c = BPACK_SafeUnderflow
      tsec = 0
      call c_bpack_gpu_dm_tsqr(gpu, comm, xid, tol_c, uf_c, rn_c, tsec, xnid)
      rank = rn_c
      hodlr_gpu_tsqr_time = hodlr_gpu_tsqr_time + tsec
      call c_bpack_gpu_dm_free(gpu, xid)
      allocate (sw(rank, r1 + r2))
      call c_bpack_gpu_dm_fetch_sw(gpu, c_loc(sw(1, 1)))
      tm0 = MPI_Wtime()
      hodlr_gpu_merge_time(2) = hodlr_gpu_merge_time(2) + tm0 - tm1
      ! the other side: son 1's rows (0 .. o1) then son 2's (o1 .. o1+o2), each times its part of sw
      o1_loc = 0
      o2_loc = 0
      if (in1) then
         p1 = ptree%MyID - ptree%pgrp(pgno1)%head + 1
         o1_loc = op1(p1, 2) - op1(p1, 1) + 1
      endif
      if (in2) then
         p2 = ptree%MyID - ptree%pgrp(pgno2)%head + 1
         o2_loc = op2(p2, 2) - op2(p2, 1) + 1
      endif
      rows_c = o1_loc + o2_loc
      cols_c = rank
      call c_bpack_gpu_dm_alloc(gpu, rows_c, cols_c, c_null_ptr, oid)
      n_c = rank
      if (in1 .and. o1_loc > 0 .and. r1 > 0) then
         allocate (swp(rank, r1))
         swp = sw(:, 1:r1)
         a_c = merge(v1, u1, col)
         m_c = o1_loc
         k_c = r1
         c0_c = 0
         call c_bpack_gpu_dm_gemm_nt(gpu, a_c, m_c, k_c, c_loc(swp(1, 1)), n_c, oid, c0_c)
         deallocate (swp)
      endif
      if (in2 .and. o2_loc > 0 .and. r2 > 0) then
         allocate (swp(rank, r2))
         swp = sw(:, r1 + 1:r1 + r2)
         a_c = merge(v2, u2, col)
         m_c = o2_loc
         k_c = r2
         c0_c = o1_loc
         call c_bpack_gpu_dm_gemm_nt(gpu, a_c, m_c, k_c, c_loc(swp(1, 1)), n_c, oid, c0_c)
         deallocate (swp)
      endif
      if (u1 /= 0) call c_bpack_gpu_dm_free(gpu, u1)
      if (v1 /= 0) call c_bpack_gpu_dm_free(gpu, v1)
      if (u2 /= 0) call c_bpack_gpu_dm_free(gpu, u2)
      if (v2 /= 0) call c_bpack_gpu_dm_free(gpu, v2)
      hodlr_gpu_merge_time(3) = hodlr_gpu_merge_time(3) + MPI_Wtime() - tm0
      ! the node's layouts: the decomposed side even, the other son 1's rows then son 2's after o1
      allocate (op(nproc, 2))
      do p = 1, nproc
         gid = ptree%pgrp(pgno)%head + p - 1
         lo = huge(1)
         hi = 0
         if (gid >= ptree%pgrp(pgno1)%head .and. gid < ptree%pgrp(pgno1)%head + nproc1) then
            p1 = gid - ptree%pgrp(pgno1)%head + 1
            lo = op1(p1, 1)
            hi = op1(p1, 2)
         endif
         if (gid >= ptree%pgrp(pgno2)%head .and. gid < ptree%pgrp(pgno2)%head + nproc2) then
            p2 = gid - ptree%pgrp(pgno2)%head + 1
            lo = min(lo, o1 + op2(p2, 1))
            hi = max(hi, o1 + op2(p2, 2))
         endif
         op(p, :) = (/lo, hi/)
      enddo
      if (col) then
         uid = xnid
         vid = oid
         call move_alloc(dg, up)
         call move_alloc(op, vp)
      else
         vid = xnid
         uid = oid
         call move_alloc(dg, vp)
         call move_alloc(op, up)
      endif
      deallocate (lay, sw, xp1, op1, xp2, op2)
      if (allocated(up1)) deallocate (up1)
      if (allocated(vp1)) deallocate (vp1)
      if (allocated(up2)) deallocate (up2)
      if (allocated(vp2)) deallocate (vp2)
      deallocate (node%sons)
   end subroutine HODLR_gpu_hmerge

   !> |U1 V1^T - U2 V2^T| / |U2 V2^T| of a block whose factors (this rank's
   !> rows, in the block's layout M_p, N_p) may be spread over the ranks of
   !> its group, from a few fixed probe vectors (no random numbers drawn)
   real(kind=8) function HODLR_gpu_lr_diff_dist(u1, v1, u2, v2, blk, k1, k2, ptree)
      implicit none
      type(matrixblock)::blk
      integer k1, k2
      DT::u1(:, :), v1(:, :), u2(:, :), v2(:, :)
      type(proctree)::ptree
      DT, allocatable :: x(:, :), t1(:, :), t2(:, :), y1(:, :), y2(:, :)
      integer, parameter :: nv = 4
      integer m, n, pp, j0, i, j, ierr
      real(kind=8) :: nd(2)

      pp = ptree%MyID - ptree%pgrp(blk%pgno)%head + 1
      m = blk%M_p(pp, 2) - blk%M_p(pp, 1) + 1
      n = blk%N_p(pp, 2) - blk%N_p(pp, 1) + 1
      j0 = blk%N_p(pp, 1)
      allocate (x(max(1, n), nv), t1(k1, nv), t2(k2, nv), y1(max(1, m), nv), y2(max(1, m), nv))
      do j = 1, nv
         do i = 1, n
#if DAT==0 || DAT==2
            x(i, j) = cmplx(sin(0.37d0*(j0 + i - 1) + 1.91d0*j), cos(0.73d0*(j0 + i - 1) - 0.53d0*j), kind=8)
#else
            x(i, j) = sin(0.37d0*(j0 + i - 1) + 1.91d0*j)
#endif
         enddo
      enddo
      t1 = 0
      t2 = 0
      y1 = 0
      y2 = 0
      if (n > 0) then
         call gemmf90(v1, size(v1, 1), x, max(1, n), t1, k1, 'T', 'N', k1, nv, n, BPACK_cone, BPACK_czero)
         call gemmf90(v2, size(v2, 1), x, max(1, n), t2, k2, 'T', 'N', k2, nv, n, BPACK_cone, BPACK_czero)
      endif
      call MPI_ALLREDUCE(MPI_IN_PLACE, t1, k1*nv, MPI_DT, MPI_SUM, ptree%pgrp(blk%pgno)%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, t2, k2*nv, MPI_DT, MPI_SUM, ptree%pgrp(blk%pgno)%Comm, ierr)
      if (m > 0) then
         call gemmf90(u1, size(u1, 1), t1, k1, y1, max(1, m), 'N', 'N', m, nv, k1, BPACK_cone, BPACK_czero)
         call gemmf90(u2, size(u2, 1), t2, k2, y2, max(1, m), 'N', 'N', m, nv, k2, BPACK_cone, BPACK_czero)
      endif
      nd(1) = fnorm(y1 - y2, max(1, m), nv)**2
      nd(2) = fnorm(y2, max(1, m), nv)**2
      call MPI_ALLREDUCE(MPI_IN_PLACE, nd, 2, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(blk%pgno)%Comm, ierr)
      HODLR_gpu_lr_diff_dist = sqrt(nd(1)/max(nd(2), tiny(1d0)))
      deallocate (x, t1, t2, y1, y2)
   end function HODLR_gpu_lr_diff_dist


   subroutine HODLR_gpu_baca_blocks(gpu, blks, nb, option, stats, msh, ptree, rnk, gout, thost, tcore, trec, npass, nblkit, &
                                    seeds, keep, dids)
      implicit none
      type(c_ptr) :: gpu
      integer nb
      type(hodlr_gpu_blkptr) :: blks(nb)
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(proctree)::ptree
      integer rnk(nb)
      integer, optional :: seeds(nb)  ! (the random first column of block b from a stream seeded by seeds(b))
      logical, optional :: keep       ! (the blocks' final factors: their device copies kept for HODLR_gpu_upload_forward)
      integer, optional :: dids(3, nb)  ! (the factors left on the device instead: ids of U, V (DistQr::dm_*) and the rank)

      type ilist
         integer, allocatable :: v(:)
      end type ilist
      type(matrixblock), pointer::blk
      type(ilist), allocatable :: kcols(:), krows(:)
      integer b, a, ii, i, j, k, r, u, nold, rn, cr, cc, iii, tmpn, edge, maxmn, nq
      integer total, off, woff, joff, goff, npass, nblkit
      integer, allocatable :: rank(:), rank0(:), itr(:), itrmax(:), rmin(:), perms(:), sel(:)
      integer, allocatable :: coff(:), jof(:), ruw(:), jpw(:)
      integer(c_int), allocatable :: mm(:), nn(:), rest(:), bl(:), rul(:), jpl(:), rskl(:), rnl(:), ncl(:), nrl(:)
      integer(c_int), allocatable :: colsl(:), rowsl(:)
      integer(c_int64_t), allocatable :: r0(:), c0(:)
      integer(c_int) :: nac, nap, bc, mode, knnflag, variant
      real(kind=8), allocatable :: normA(:), normUV(:)
      logical, allocatable :: active(:)
      DT, allocatable, target :: cores(:), ww(:), wl(:), grams(:)
      type(c_ptr), allocatable :: uptr(:), vptr(:)
      real(c_double) :: scale, gout(6)
      real(kind=8) :: tol, flops, fl, s, inner, t0, thost, tcore, trec
      real(c_double) :: tol_c, uf_c
      integer(c_int) :: keep_c, idu, idv, rows_c, cols_c
      integer(c_int), allocatable :: idl(:)
      integer nseed
      integer, allocatable :: seed_save(:), seed_b(:)

      if (nb == 0) return
      if (present(seeds)) then
         call random_seed(size=nseed)
         allocate (seed_save(nseed), seed_b(nseed))
         call random_seed(get=seed_save)
      endif
      maxmn = 1
      allocate (r0(nb), c0(nb), mm(nb), nn(nb), rest(nb), rank(nb), rank0(nb), itr(nb), itrmax(nb), rmin(nb))
      allocate (normA(nb), normUV(nb), active(nb), bl(nb), rul(nb), rskl(nb), rnl(nb), coff(nb + 1), jof(nb + 1))
      allocate (ruw(nb))
      do b = 1, nb
         blk => blks(b)%p
         r0(b) = blk%headm - 1
         c0(b) = blk%headn - 1
         mm(b) = blk%M
         nn(b) = blk%N
         rest(b) = min(option%BACA_Batch, min(blk%M, blk%N))
         rmin(b) = min(blk%M, blk%N)
         itrmax(b) = floor_safe(min(blk%M, blk%N)/dble(rest(b)))*2
         maxmn = max(maxmn, blk%M, blk%N)
      enddo

      scale = option%scale_factor
      mode = option%HODLR_use_gpu
      variant = option%RecLR_leaf
      nac = nb
      call c_bpack_hodlr_gpu_baca_begin(gpu, scale, mode, variant, nac, r0, mm, c0, nn, rest)
      tol = option%tol_comp
      rank = 0
      rank0 = 0
      itr = 0
      normA = 0
      normUV = 0
      active = .true.
      thost = 0
      tcore = 0
      trec = 0
      flops = 0
      allocate (perms(maxmn))

      !>**** the nearest neighbours first (option%knn > 0, BACA without overlap)
      if (option%knn > 0 .and. option%RecLR_leaf == BACANOVER) then
         t0 = MPI_Wtime()
         allocate (kcols(nb), krows(nb), ncl(nb), nrl(nb))
         nac = 0
         do b = 1, nb
            blk => blks(b)%p
            allocate (sel(max(blk%M, blk%N)*option%knn))
            cc = 0
            do i = 1, blk%M
               edge = blk%headm + i - 1
               do iii = 1, option%knn
                  if (msh%nns(edge, iii) >= blk%headn .and. msh%nns(edge, iii) <= blk%headn + blk%N - 1) then
                     cc = cc + 1
                     sel(cc) = msh%nns(edge, iii) + 1 - blk%headn
                  endif
               enddo
            enddo
            tmpn = cc
            if (cc > 0) call remove_dup_int(sel, tmpn, cc)
            if (cc > 0) kcols(b)%v = sel(1:cc)
            cr = 0
            do j = 1, blk%N
               edge = blk%headn + j - 1
               do iii = 1, option%knn
                  if (msh%nns(edge, iii) >= blk%headm .and. msh%nns(edge, iii) <= blk%headm + blk%M - 1) then
                     cr = cr + 1
                     sel(cr) = msh%nns(edge, iii) + 1 - blk%headm
                  endif
               enddo
            enddo
            tmpn = cr
            if (cr > 0) call remove_dup_int(sel, tmpn, cr)
            if (cr > 0) krows(b)%v = sel(1:cr)
            deallocate (sel)
            if (cr > 0 .and. cc > 0) then
               nac = nac + 1
               bl(nac) = b
               ncl(nac) = cc
               nrl(nac) = cr
            endif
         enddo
         if (nac > 0) then
            allocate (colsl(max(1, sum(ncl(1:nac)))), rowsl(max(1, sum(nrl(1:nac)))))
            coff(1) = 0
            jof(1) = 0
            off = 0
            joff = 0
            do a = 1, nac
               b = bl(a)
               colsl(off + 1:off + ncl(a)) = kcols(b)%v
               rowsl(joff + 1:joff + nrl(a)) = krows(b)%v
               off = off + ncl(a)
               joff = joff + nrl(a)
               coff(a + 1) = coff(a) + ncl(a)*nrl(a)
               jof(a + 1) = jof(a) + ncl(a)
            enddo
            allocate (cores(max(1, coff(nac + 1))), ww(max(1, coff(nac + 1))), jpw(max(1, jof(nac + 1))))
            thost = thost + MPI_Wtime() - t0
            call c_bpack_hodlr_gpu_baca_knn_panels(gpu, nac, bl, ncl, colsl, nrl, rowsl, c_loc(cores))
            t0 = MPI_Wtime()
#ifdef HAVE_OPENMP
            !$omp parallel do default(shared) private(a, b, fl) reduction(+:flops) schedule(dynamic) &
            !$omp num_threads(HODLR_gpu_host_threads())
#endif
            do a = 1, nac
               b = bl(a)
               call HODLR_gpu_baca_core(nrl(a), ncl(a), cores(coff(a) + 1), tol, rmin(b), huge(1), ruw(a), &
                                        jpw(jof(a) + 1), ww(coff(a) + 1), fl)
               flops = flops + fl
            enddo
#ifdef HAVE_OPENMP
            !$omp end parallel do
#endif
            allocate (wl(max(1, coff(nac + 1))), jpl(max(1, jof(nac + 1))), grams(1))
            nap = 0
            woff = 0
            joff = 0
            do a = 1, nac
               b = bl(a)
               u = ruw(a)
               if (u <= 0) cycle
               nap = nap + 1
               rul(nap) = u
               rskl(nap) = 0
               jpl(joff + 1:joff + u) = jpw(jof(a) + 1:jof(a) + u)
               joff = joff + u
               wl(woff + 1:woff + u*nrl(a)) = ww(coff(a) + 1:coff(a) + u*nrl(a))
               woff = woff + u*nrl(a)
               rank(b) = u
               rank0(b) = u
               bl(nap) = b
            enddo
            thost = thost + MPI_Wtime() - t0
            knnflag = 1
            if (nap > 0) call c_bpack_hodlr_gpu_baca_append(gpu, knnflag, nap, bl, rul, jpl, c_loc(wl), rskl, &
                                                            c_loc(grams))
            deallocate (colsl, rowsl, cores, ww, jpw, wl, jpl, grams)
         endif
         deallocate (kcols, krows, ncl, nrl)
      endif

      !>**** BACA iterations
      knnflag = 0
      npass = 0
      nblkit = 0
      do
         ! the blocks that go on (the loop test of the CPU), with their first columns
         nac = 0
         do b = 1, nb
            if (.not. active(b)) cycle
            if (.not. (normUV(b) >= tol*normA(b) .and. itr(b) < itrmax(b))) then
               active(b) = .false.
               cycle
            endif
            if (rank(b) == 0) then
               if (present(seeds)) then
                  do i = 1, nseed
                     seed_b(i) = mod(seeds(b), 100000)*7919 + i*104729 + 1
                  enddo
                  call random_seed(put=seed_b)
               endif
               call rperm(nn(b), perms)
               bc = b
               call c_bpack_hodlr_gpu_baca_set_columns(gpu, bc, perms)
            endif
            nac = nac + 1
            bl(nac) = b
         enddo
         if (nac == 0) exit
         npass = npass + 1
         nblkit = nblkit + nac

         coff(1) = 0
         jof(1) = 0
         do a = 1, nac
            coff(a + 1) = coff(a) + rest(bl(a))**2
            jof(a + 1) = jof(a) + rest(bl(a))
         enddo
         allocate (cores(max(1, coff(nac + 1))), ww(max(1, coff(nac + 1))), jpw(max(1, jof(nac + 1))))
         call c_bpack_hodlr_gpu_baca_panels(gpu, nac, bl, c_loc(cores))

         ! the cores: rank-revealing QR, W = R11^{-1} Q1^H, the rank update
         t0 = MPI_Wtime()
#ifdef HAVE_OPENMP
         !$omp parallel do default(shared) private(a, b, r, fl) reduction(+:flops) schedule(dynamic) &
         !$omp num_threads(HODLR_gpu_host_threads())
#endif
         do a = 1, nac
            b = bl(a)
            r = rest(b)
            call HODLR_gpu_baca_core(r, r, cores(coff(a) + 1), tol, huge(1), rmin(b) - rank(b), ruw(a), &
                                     jpw(jof(a) + 1), ww(coff(a) + 1), fl)
            flops = flops + fl
         enddo
#ifdef HAVE_OPENMP
         !$omp end parallel do
#endif
         allocate (wl(max(1, coff(nac + 1))), jpl(max(1, jof(nac + 1))))
         nap = 0
         woff = 0
         joff = 0
         goff = 0
         do a = 1, nac
            b = bl(a)
            r = rest(b)
            u = ruw(a)
            if (u <= 0) then  ! (no rank found, or the rank is full)
               active(b) = .false.
               cycle
            endif
            nap = nap + 1
            bl(nap) = b
            rul(nap) = u
            rskl(nap) = rank0(b)
            jpl(joff + 1:joff + u) = jpw(jof(a) + 1:jof(a) + u)
            joff = joff + u
            wl(woff + 1:woff + u*r) = ww(coff(a) + 1:coff(a) + u*r)
            woff = woff + u*r
            goff = goff + 2*u*u + 2*(rank(b) - rank0(b))*u
         enddo
         thost = thost + MPI_Wtime() - t0
         tcore = tcore + MPI_Wtime() - t0

         if (nap > 0) then
            allocate (grams(max(1, goff)))
            call c_bpack_hodlr_gpu_baca_append(gpu, knnflag, nap, bl, rul, jpl, c_loc(wl), rskl, c_loc(grams))
            ! norm of the update (LR_Fnorm) and of the approximation (LR_FnormUp)
            goff = 0
            do a = 1, nap
               b = bl(a)
               u = rul(a)
               nold = rank(b) - rank0(b)
               s = 0
               do i = 1, u*u
                  s = s + dble(grams(goff + i)*grams(goff + u*u + i))
               enddo
               normUV(b) = sqrt(max(s, 0d0))
               goff = goff + 2*u*u
               inner = 0
               do i = 1, nold*u
                  inner = inner + 2*dble(grams(goff + i)*grams(goff + nold*u + i))
               enddo
               goff = goff + 2*nold*u
               rank(b) = rank(b) + u
               normA(b) = sqrt(normA(b)**2d0 + normUV(b)**2d0 + inner)
               if (normA(b) > BPACK_SafeUnderflow) then
                  itr(b) = itr(b) + 1
               else
                  active(b) = .false.
               endif
            enddo
            deallocate (grams)
         endif
         deallocate (cores, ww, jpw, wl, jpl)
      enddo

      !>**** recompression (LR_ReCompression) of the blocks of positive rank
      nq = 0
      do b = 1, nb
         if (rank(b) > 0) then
            nq = nq + 1
            bl(nq) = b
         endif
      enddo
      if (nq > 0) then
         ! all on the device (HodlrConstruct::recompress), then the new factors to the host, or
         ! left on the device as merge matrices (dids, for HODLR_gpu_hmerge)
         nac = nq
         tol_c = option%tol_comp
         uf_c = BPACK_SafeUnderflow
         call c_bpack_hodlr_gpu_baca_recompress(gpu, nac, bl, tol_c, uf_c, rnl)
         if (present(dids)) then
            allocate (idl(2*nq))
            call c_bpack_hodlr_gpu_baca_to_dm(gpu, nac, bl, rnl, idl)
            do a = 1, nq
               dids(:, bl(a)) = (/idl(2*a - 1), idl(2*a), rnl(a)/)
            enddo
            deallocate (idl)
         else
            allocate (uptr(nq), vptr(nq))
            t0 = MPI_Wtime()
            do a = 1, nq
               b = bl(a)
               rn = rnl(a)
               blk => blks(b)%p
               allocate (blk%ButterflyU%blocks(1), blk%ButterflyV%blocks(1))
               allocate (blk%ButterflyU%blocks(1)%matrix(blk%M, rn), blk%ButterflyV%blocks(1)%matrix(blk%N, rn))
               uptr(a) = c_loc(blk%ButterflyU%blocks(1)%matrix(1, 1))
               vptr(a) = c_loc(blk%ButterflyV%blocks(1)%matrix(1, 1))
            enddo
            thost = thost + MPI_Wtime() - t0
            trec = trec + MPI_Wtime() - t0
            keep_c = 0
            if (present(keep)) then
               if (keep) keep_c = 1
               if (keep .and. HODLR_gpu_defer_host(option)) keep_c = 2  ! (kept, the host arrays left unfilled)
            endif
            call c_bpack_hodlr_gpu_baca_download(gpu, nac, bl, rnl, uptr, vptr, keep_c)
            deallocate (uptr, vptr)
         endif
      endif
      call c_bpack_hodlr_gpu_baca_end(gpu, gout)
      stats%Flop_Fill = stats%Flop_Fill + gout(5) + flops
      stats%Time_Entry = stats%Time_Entry + gout(6)  ! (the device time of the entry evaluations)

      ! the ranks; the blocks of rank 0 get zero factors of rank 1 (as LR_BACA_noOverlap)
      nq = 0
      do b = 1, nb
         blk => blks(b)%p
         if (rank(b) > 0) then
            nq = nq + 1
            rnk(b) = rnl(nq)
         elseif (present(dids)) then
            rnk(b) = 1
            rows_c = blk%M
            cols_c = 1
            call c_bpack_gpu_dm_alloc(gpu, rows_c, cols_c, c_null_ptr, idu)
            rows_c = blk%N
            call c_bpack_gpu_dm_alloc(gpu, rows_c, cols_c, c_null_ptr, idv)
            dids(:, b) = (/idu, idv, 1/)
         else
            rnk(b) = 1
            allocate (blk%ButterflyU%blocks(1), blk%ButterflyV%blocks(1))
            allocate (blk%ButterflyU%blocks(1)%matrix(blk%M, 1), blk%ButterflyV%blocks(1)%matrix(blk%N, 1))
            blk%ButterflyU%blocks(1)%matrix = 0
            blk%ButterflyV%blocks(1)%matrix = 0
         endif
      enddo
      deallocate (r0, c0, mm, nn, rest, rank, rank0, itr, itrmax, rmin, normA, normUV, active, bl, rul, rskl, rnl)
      deallocate (coff, jof, ruw, perms)
      if (present(seeds)) then
         call random_seed(put=seed_save)
         deallocate (seed_save, seed_b)
      endif
   end subroutine HODLR_gpu_baca_blocks

   !> The host step of BACA without overlap on one core C(I, J) (cr x cc),
   !> as LR_BACA_noOverlap takes it: the rank of the update by the truncated
   !> column-pivoted QR (geqp3modf90, tolerance tol), capped at pre_cap
   !> before W and at post_cap after it (the CPU caps the knn step before and
   !> the iterations after), jp(1:rankup) the pivots and w(1:rankup*cr) the
   !> first rankup rows of W = R11^{-1} Q1^H (rankup x cr); rankup <= 0: no
   !> update
   subroutine HODLR_gpu_baca_core(cr, cc, core, tol, pre_cap, post_cap, rankup, jp, w, flops)
      implicit none
      integer cr, cc, pre_cap, post_cap, rankup, jp(*), i, ranknew
      DT::core(cr, cc), w(*)
      real(kind=8) :: tol, flops, flop
      DT, allocatable :: a(:, :), tau(:), wt(:, :)
      integer, allocatable :: jpvt(:)

      flops = 0
      allocate (a(cr, cc), tau(max(cr, cc)), jpvt(max(cr, cc)))
      a = core
      jpvt = 0
      call geqp3modf90(a, jpvt, tau, tol, BPACK_SafeUnderflow, ranknew, flop=flop)
      flops = flops + flop
      rankup = min(ranknew, pre_cap)
      if (rankup > 0) then
         allocate (wt(cr, cr))
         wt = 0
         do i = 1, cr
            wt(i, i) = 1
         enddo
         call un_or_mqrf90(a, tau, wt, 'L', 'C', cr, cr, rankup, flop=flop)
         flops = flops + flop
         call trsmf90(a, wt, 'L', 'U', 'N', 'N', rankup, cr, flop=flop)
         flops = flops + flop
         rankup = min(rankup, post_cap)
         if (rankup > 0) then
            w(1:rankup*cr) = reshape(wt(1:rankup, 1:cr), (/rankup*cr/))
            jp(1:rankup) = jpvt(1:rankup)
         endif
         deallocate (wt)
      endif
      deallocate (a, tau, jpvt)
   end subroutine HODLR_gpu_baca_core

   !> The dense leaves (level Maxlevel+1) evaluated on the GPU, as
   !> Full_construction does on the CPU (symmetrized when option%sym > 0, as
   !> HODLR_construction does); with BPACK_CHECK=hodlr compared with the CPU
   subroutine HODLR_gpu_construct_leaves(ho_bf1, option, stats, msh, ker, ptree)
      implicit none
      real(kind=8) :: tev
      type(hobf), target::ho_bf1
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(matrixblock), pointer::blk
      type(matrixblock)::ref
      type(Hstat)::stats_chk
      integer level_c, ii, nb, b, i, j
      integer(c_int), allocatable :: mm(:), nn(:)
      integer(c_int64_t), allocatable :: r0(:), c0(:)
      type(c_ptr), allocatable :: ptrs(:)
      integer(c_int) :: count
      real(c_double) :: scale
      real(kind=8) :: memory, err, errmax, t0, t1
      DT::sym_value

      level_c = ho_bf1%Maxlevel + 1
      nb = 0
      do ii = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
         if (IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(ii)%pgno)) nb = nb + 1
      enddo
      if (nb == 0) return
      t0 = MPI_Wtime()
      allocate (mm(nb), nn(nb), r0(nb), c0(nb), ptrs(nb))
      b = 0
      do ii = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
         if (.not. IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(ii)%pgno)) cycle
         b = b + 1
         blk => ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)
         allocate (blk%fullmat(blk%M, blk%N))
         if (blk%row_group == blk%col_group) allocate (blk%ipiv(blk%M))
         mm(b) = blk%M
         nn(b) = blk%N
         r0(b) = blk%headm - 1
         c0(b) = blk%headn - 1
         ptrs(b) = c_loc(blk%fullmat(1, 1))
      enddo
      scale = option%scale_factor
      count = nb
      tev = MPI_Wtime()
      call c_bpack_hodlr_gpu_eval_dense(ho_bf1%gpu, scale, count, r0, mm, c0, nn, ptrs)
      t1 = MPI_Wtime()
      stats%Time_Entry = stats%Time_Entry + t1 - tev  ! (the entries, evaluated on the device, on the host)

      errmax = 0
      if (HODLR_gpu_check_level() > 0) call InitStat(stats_chk)
      do ii = ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
         if (.not. IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(ii)%pgno)) cycle
         blk => ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)
         if (HODLR_gpu_check_level() > 0) then
            ref%M = blk%M
            ref%N = blk%N
            ref%headm = blk%headm
            ref%headn = blk%headn
            ref%row_group = blk%row_group
            ref%col_group = blk%col_group
            call Full_construction(ref, msh, ker, stats_chk, option, ptree, memory)
            err = fnorm(ref%fullmat - blk%fullmat, blk%M, blk%N)/max(fnorm(ref%fullmat, blk%M, blk%N), tiny(err))
            errmax = max(errmax, err)
            deallocate (ref%fullmat)
            if (allocated(ref%ipiv)) deallocate (ref%ipiv)
         endif
         if (option%sym > 0) then
            call assert(blk%M == blk%N, 'a symmetric HODLR dense leaf must be square')
            do j = 1, blk%N
               do i = j + 1, blk%M
                  sym_value = (blk%fullmat(i, j) + blk%fullmat(j, i))/2d0
                  blk%fullmat(i, j) = sym_value
                  blk%fullmat(j, i) = sym_value
               enddo
            enddo
         endif
         memory = SIZEOF(blk%fullmat)/1024.0d3
         stats%Mem_Direct_for = stats%Mem_Direct_for + memory
      enddo
      if (HODLR_gpu_check_level() > 0 .and. ptree%MyID == Main_ID) then
         write (*, '(A,I8,A,Es11.3)') ' HODLR GPU check (dense leaves):', nb, &
            ' leaves, largest relative difference to the CPU', errmax
      endif
      if (ptree%MyID == Main_ID .and. option%verbosity >= 1) then
         write (*, '(A,I8,A,F9.3,A)') ' HODLR GPU construction: ', nb, ' dense leaves in', t1 - t0, ' s'
      endif
      deallocate (mm, nn, r0, c0, ptrs)
   end subroutine HODLR_gpu_construct_leaves

   !> BPACK_CHECK=hodlr: rebuild the blocks of level level_c on the CPU
   !> (BP_compress_entry, from the random numbers the GPU run drew) and print
   !> the largest relative difference of U V^T and the rank differences;
   !> the GPU blocks stay.
   subroutine HODLR_gpu_check_construct_level(ho_bf1, level_c, option, msh, ker, ptree, flop_gpu)
      implicit none
      type(hobf), target::ho_bf1
      integer level_c
      type(Hoption)::option
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      real(kind=8) :: flop_gpu  ! (this rank's flops of the GPU construction of the level)
      real(kind=8) :: flop_cpu
      type(matrixblock), pointer::blk
      type(Hstat)::stats_chk
      DT, pointer :: ug(:, :), vg(:, :)
      integer ii, kg, kc, nb, ndiff, maxdk, sumg, sumc, ierr
      real(kind=8) :: rtemp, err, errmax

      call InitStat(stats_chk)
      allocate (stats_chk%rankmax_of_level(0:ho_bf1%Maxlevel))
      stats_chk%rankmax_of_level = 0
      nb = 0
      ndiff = 0
      maxdk = 0
      sumg = 0
      sumc = 0
      errmax = 0
      do ii = ho_bf1%levels(level_c)%Bidxs*2 - 1, ho_bf1%levels(level_c)%Bidxe*2
         if (option%sym > 0 .and. mod(ii, 2) == 1) cycle
         if (.not. IOwnPgrp(ptree, ho_bf1%levels(level_c)%BP(ii)%pgno)) cycle
         nb = nb + 1
         blk => ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)
         ug => blk%ButterflyU%blocks(1)%matrix
         vg => blk%ButterflyV%blocks(1)%matrix
         kg = blk%rankmax
         nullify (blk%ButterflyU%blocks(1)%matrix, blk%ButterflyV%blocks(1)%matrix)
         deallocate (blk%ButterflyU%blocks, blk%ButterflyV%blocks)
         call BP_compress_entry(ho_bf1%levels(level_c)%BP(ii), option, rtemp, stats_chk, msh, ker, ptree)
         kc = blk%rankmax
         err = HODLR_gpu_lr_diff_dist(ug, vg, blk%ButterflyU%blocks(1)%matrix, blk%ButterflyV%blocks(1)%matrix, &
                                      blk, kg, kc, ptree)
         errmax = max(errmax, err)
         if (ptree%MyID == ptree%pgrp(blk%pgno)%head) then  ! (a shared block counts once)
            if (kc /= kg) ndiff = ndiff + 1
            maxdk = max(maxdk, abs(kc - kg))
            sumg = sumg + kg
            sumc = sumc + kc
         else
            nb = nb - 1
         endif
         deallocate (blk%ButterflyU%blocks(1)%matrix, blk%ButterflyV%blocks(1)%matrix)
         blk%ButterflyU%blocks(1)%matrix => ug
         blk%ButterflyV%blocks(1)%matrix => vg
         blk%rankmax = kg
         blk%rankmin = kg
         ho_bf1%levels(level_c)%BP(ii)%LL(1)%rankmax = kg
      enddo
      call MPI_ALLREDUCE(MPI_IN_PLACE, errmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, maxdk, 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, nb, 1, MPI_INTEGER, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, ndiff, 1, MPI_INTEGER, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, sumg, 1, MPI_INTEGER, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, sumc, 1, MPI_INTEGER, MPI_SUM, ptree%Comm, ierr)
      flop_cpu = stats_chk%Flop_Fill
      call MPI_ALLREDUCE(MPI_IN_PLACE, flop_cpu, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE, flop_gpu, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID) then
         write (*, '(A,I3,A,I6,A,Es11.3,A,I6,A,I4,A,I8,A,I8,A,Es11.3,A,Es11.3)') ' HODLR GPU check (construction) level', &
            level_c, ':', nb, ' blocks, largest relative difference of U V^T to the CPU', errmax, '; ranks differ in', &
            ndiff, ' blocks, by up to', maxdk, '; rank sums GPU', sumg, ' CPU', sumc, '; flops GPU', flop_gpu, ' CPU', flop_cpu
      endif
      deallocate (stats_chk%rankmax_of_level)
   end subroutine HODLR_gpu_check_construct_level


   ! Expose A12=A21^T without duplicating A21=U1*V0^T.  The transpose view
   ! stays on the A21 block row and swaps both the factors and their 1D
   ! layouts; callers redistribute its input/output through the parent row.
   subroutine HODLR_sym_link_pair(b12, b21, parent_pgno, ptree, memory)
      implicit none
      type(blockplus)::b12, b21
      type(proctree)::ptree
      type(matrixblock), pointer::block12, block21
      integer parent_pgno, rank, rank_loc, ierr
      real(kind=8)::memory

      call assert(b12%Lplus == 1 .and. b21%Lplus == 1, &
         'symmetric HODLR path requires option%LRlevel=0')
      block12 => b12%LL(1)%matrices_block(1)
      block21 => b21%LL(1)%matrices_block(1)

      rank_loc = 0
      if (IOwnPgrp(ptree, block21%pgno)) rank_loc = block21%rankmax
      call MPI_ALLREDUCE(rank_loc, rank, 1, MPI_INTEGER, MPI_MAX, &
         ptree%pgrp(parent_pgno)%Comm, ierr)

      block12%rankmax = rank
      block12%rankmin = rank
      block21%rankmax = rank
      block21%rankmin = rank
      b12%LL(1)%rankmax = rank
      b12%pgno = b21%pgno
      block12%pgno = block21%pgno
      block12%pgno_db = block21%pgno_db
      block12%M_loc = block21%N_loc
      block12%N_loc = block21%M_loc
      call assert(size(block12%M_p, 1) == size(block21%N_p, 1), &
         'symmetric HODLR transpose-view row layout mismatch')
      call assert(size(block12%N_p, 1) == size(block21%M_p, 1), &
         'symmetric HODLR transpose-view column layout mismatch')
      block12%M_p = block21%N_p
      block12%N_p = block21%M_p
      block12%is_transpose_view = 1

      block12%ButterflyU%num_blk = block21%ButterflyV%num_blk
      block12%ButterflyV%num_blk = block21%ButterflyU%num_blk
      block12%ButterflyU%idx = block21%ButterflyV%idx
      block12%ButterflyV%idx = block21%ButterflyU%idx
      block12%ButterflyU%inc = block21%ButterflyV%inc
      block12%ButterflyV%inc = block21%ButterflyU%inc
      block12%ButterflyU%nblk_loc = 0
      block12%ButterflyV%nblk_loc = 0

      memory = 0
      if (IOwnPgrp(ptree, block21%pgno)) then
         allocate(block12%ButterflyU%blocks(1))
         allocate(block12%ButterflyV%blocks(1))
         block12%ButterflyU%blocks(1)%matrix => block21%ButterflyV%blocks(1)%matrix
         block12%ButterflyV%blocks(1)%matrix => block21%ButterflyU%blocks(1)%matrix
         block12%ButterflyU%nblk_loc = 1
         block12%ButterflyV%nblk_loc = 1
         memory = (SIZEOF(block12%ButterflyU%blocks) + &
            SIZEOF(block12%ButterflyV%blocks))/1024.0d3
      endif
   end subroutine HODLR_sym_link_pair

   subroutine HSS_construction(hss_bf1, option, stats, msh, ker, ptree)

      implicit none
      real(kind=8) n1, n2, n3, n4, n5
      integer i, j, ii, ii_inv, jj, kk, iii, jjj, ll, rankmax
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
      call InitStat(stats)
      tim_tmp = 0
      !tolerance=0.001
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      ! write (*,*) 'ACA error threshold',tolerance
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'SVD error threshold', option%tol_comp
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "constructing HSS-BF......"

      n1 = MPI_Wtime()

      allocate (stats%rankmax_of_level(0:hss_bf1%Maxlevel))
      stats%rankmax_of_level = 0
      allocate (stats%rankmax_of_level_global(0:hss_bf1%Maxlevel))
      stats%rankmax_of_level_global = 0

      call BP_compress_entry(hss_bf1%BP, option, rtemp, stats, msh, ker, ptree)
      stats%Mem_Comp_for = stats%Mem_Comp_for + rtemp

      passflag = 0
      do while (passflag < 2)
         ! call element_Zmn_block_user(0, 0, mrange_dummy, nrange_dummy, mat_dummy, msh, option, ker, 2, passflag, ptree, stats)
         call element_Zmn_blocklist_user(submats, 0, msh, option, ker, 2, passflag, ptree, stats)
      enddo

      do ll = 1, hss_bf1%BP%Lplus
         rankmax = hss_bf1%BP%LL(ll)%rankmax
         call MPI_ALLREDUCE(rankmax, hss_bf1%BP%LL(ll)%rankmax, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(hss_bf1%BP%LL(1)%matrices_block(1)%pgno)%Comm, ierr)
      enddo


      n2 = MPI_Wtime()
      stats%Time_Fill = stats%Time_Fill + n2 - n1

      stats%Mem_Fill = stats%Mem_Comp_for + stats%Mem_Direct_for
      call LogMemory(stats, stats%Mem_Fill)

      call MPI_ALLREDUCE(stats%rankmax_of_level(0:hss_bf1%Maxlevel), stats%rankmax_of_level_global(0:hss_bf1%Maxlevel), hss_bf1%Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'rankmax_of_level:', stats%rankmax_of_level_global
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      call MPI_ALLREDUCE(stats%Time_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Total construction time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Time_Entry, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Total entry eval time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Flop_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A26,Es14.2)') 'Total construction flops:', rtemp
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


   subroutine HSS_MD_construction(Ndim, hss_bf_md1, option, stats, msh, ker, ptree)

      implicit none
      integer Ndim
      real(kind=8) n1, n2, n3, n4, n5
      integer i, j, ii, ii_inv, jj, kk, iii, jjj, ll, rankmax
      integer level, blocks, edge, patch, node, group
      integer rank, index_near, m, n, length, flag, itemp, rank0_inner, rank0_outter, ierr
      real T0
      real(kind=8):: rtemp, rtemp1, rel_error, error, t1, t2, tim_tmp, rankrate_inner, rankrate_outter
      integer mm, nn, header_m, header_n, edge_m, edge_n, group_m, group_n, group_m1, group_n1, group_m2, group_n2
      type(matrixblock)::block_tmp, block_tmp1
      DT, allocatable::fullmat(:, :)
      integer level_c, iter, level_cc, level_butterfly, Bidxs, Bidxe
      type(Hoption)::option
      type(Hstat)::stats
      type(hssbf_md)::hss_bf_md1
      type(mesh)::msh(Ndim)
      type(kernelquant)::ker
      type(proctree)::ptree
      integer::passflag = 0
      type(intersect_MD) :: subtensors_dummy(1)
      integer::mrange_dummy(1), nrange_dummy(1)
      DT::mat_dummy(1, 1)

      ! Memory_direct_forward=0
      ! Memory_butterfly_forward=0
      call InitStat(stats)
      tim_tmp = 0
      !tolerance=0.001
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      ! write (*,*) 'ACA error threshold',tolerance
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'SVD error threshold', option%tol_comp
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "constructing HSS-BF-MD......"

      n1 = MPI_Wtime()

      allocate (stats%rankmax_of_level(0:hss_bf_md1%Maxlevel))
      stats%rankmax_of_level = 0
      allocate (stats%rankmax_of_level_global(0:hss_bf_md1%Maxlevel))
      stats%rankmax_of_level_global = 0

      call BP_MD_compress_entry(Ndim, hss_bf_md1%BP, option, rtemp, rtemp1,stats, msh, ker, ptree, logn_level_flag=0)
      stats%Mem_Comp_for = stats%Mem_Comp_for + rtemp
      stats%Mem_Direct_for = stats%Mem_Direct_for + rtemp1

      passflag = 0
      do while (passflag < 2)
         call element_Zmn_tensorlist_user(Ndim, subtensors_dummy, 0, msh, option, ker, 2, passflag, ptree, stats)
      enddo

      do ll = 1, hss_bf_md1%BP%Lplus
         rankmax = hss_bf_md1%BP%LL(ll)%rankmax
         call MPI_ALLREDUCE(rankmax, hss_bf_md1%BP%LL(ll)%rankmax, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(hss_bf_md1%BP%LL(1)%matrices_block(1)%pgno)%Comm, ierr)
      enddo


      n2 = MPI_Wtime()
      stats%Time_Fill = stats%Time_Fill + n2 - n1

      stats%Mem_Fill = stats%Mem_Comp_for + stats%Mem_Direct_for
      ! if (ptree%MyID == Main_ID)write(*,*)stats%Mem_Fill,stats%Mem_Comp_for,stats%Mem_Direct_for,stats%Mem_Peak,'wocao'
      call LogMemory(stats, stats%Mem_Fill)
      ! if (ptree%MyID == Main_ID)write(*,*)stats%Mem_Fill,stats%Mem_Comp_for,stats%Mem_Direct_for,stats%Mem_Peak,'wocao after logmemory'

      call MPI_ALLREDUCE(stats%rankmax_of_level(0:hss_bf_md1%Maxlevel), stats%rankmax_of_level_global(0:hss_bf_md1%Maxlevel), hss_bf_md1%Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'rankmax_of_level:', stats%rankmax_of_level_global
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      call MPI_ALLREDUCE(stats%Time_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Total construction time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Time_Entry, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Total entry eval time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Flop_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A26,Es14.2)') 'Total construction flops:', rtemp
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

   end subroutine HSS_MD_construction





   subroutine HMAT_MD_construction(Ndim, h_mat_md, option, stats, msh, ker, ptree)

      implicit none
      integer Ndim
      real(kind=8) n1, n2, n3, n4, n5
      integer i, j, ii, ii_inv, jj, kk, iii, jjj, ll, rankmax
      integer level, blocks, edge, patch, node, group
      integer rank, index_near, m, n, length, flag, itemp, rank0_inner, rank0_outter, ierr
      real T0
      real(kind=8):: rtemp,rtemp1, rel_error, error, t1, t2, tim_tmp, rankrate_inner, rankrate_outter
      integer mm, nn, header_m, header_n, edge_m, edge_n, group_m, group_n, group_m1, group_n1, group_m2, group_n2
      type(matrixblock)::block_tmp, block_tmp1
      DT, allocatable::fullmat(:, :)
      integer level_c, iter, level_cc, level_butterfly, Bidxs, Bidxe
      type(Hoption)::option
      type(Hstat)::stats
      type(Hmat_md)::h_mat_md
      type(mesh)::msh(Ndim)
      type(kernelquant)::ker
      type(proctree)::ptree
      integer::passflag = 0, bp_cnt_lr_tmp
      type(intersect_MD) :: subtensors_dummy(1)
      integer::mrange_dummy(1), nrange_dummy(1)
      DT::mat_dummy(1, 1)

      ! Memory_direct_forward=0
      ! Memory_butterfly_forward=0
      call InitStat(stats)
      tim_tmp = 0
      !tolerance=0.001
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      ! write (*,*) 'ACA error threshold',tolerance
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'SVD error threshold', option%tol_comp
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "constructing HTENSOR......"

      n1 = MPI_Wtime()

      allocate (stats%rankmax_of_level(0:h_mat_md%Maxlevel))
      stats%rankmax_of_level = 0
      allocate (stats%rankmax_of_level_global(0:h_mat_md%Maxlevel))
      stats%rankmax_of_level_global = 0
      bp_cnt_lr_tmp=option%bp_cnt_lr
      option%bp_cnt_lr=1
      call BP_MD_compress_entry(Ndim, h_mat_md%BP, option, rtemp, rtemp1, stats, msh, ker, ptree, logn_level_flag=1)
      option%bp_cnt_lr=bp_cnt_lr_tmp
      stats%Mem_Comp_for = stats%Mem_Comp_for + rtemp
      stats%Mem_Direct_for = stats%Mem_Direct_for + rtemp1

      passflag = 0
      do while (passflag < 2)
         call element_Zmn_tensorlist_user(Ndim, subtensors_dummy, 0, msh, option, ker, 2, passflag, ptree, stats)
      enddo

      do ll = 1, h_mat_md%BP%Lplus
         rankmax = h_mat_md%BP%LL(ll)%rankmax
         call MPI_ALLREDUCE(rankmax, h_mat_md%BP%LL(ll)%rankmax, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(h_mat_md%BP%LL(1)%matrices_block(1)%pgno)%Comm, ierr)
      enddo


      n2 = MPI_Wtime()
      stats%Time_Fill = stats%Time_Fill + n2 - n1

      stats%Mem_Fill = stats%Mem_Comp_for + stats%Mem_Direct_for
      call LogMemory(stats, stats%Mem_Fill)

      call MPI_ALLREDUCE(stats%rankmax_of_level(0:h_mat_md%Maxlevel), stats%rankmax_of_level_global(0:h_mat_md%Maxlevel), h_mat_md%Maxlevel + 1, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'rankmax_of_level:', stats%rankmax_of_level_global
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) ''
      call MPI_ALLREDUCE(stats%Time_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Total construction time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Time_Entry, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Total entry eval time:', rtemp, 'Seconds'
      call MPI_ALLREDUCE(stats%Flop_Fill, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A26,Es14.2)') 'Total construction flops:', rtemp
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

   end subroutine HMAT_MD_construction









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
      integer:: boundary_map(1,1)
      integer level_butterfly, levelm, groupm_start, Nboundall, Ninadmissible
      if (IOwnPgrp(ptree, blocks%pgno))then
      ! t1=MPI_Wtime()
      if (blocks%style == 2) then

         groupm_start = 0
         Nboundall = 0
         Ninadmissible = 0

         if (option%forwardN15flag == 1) then
            call BF_compress_N15(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory_tmp, stats, msh, ker, ptree, 1)
            call BF_sym2asym(blocks)
         else
            call BF_compress_NlogN(blocks, boundary_map, Nboundall, Ninadmissible, groupm_start, option, Memory_tmp, stats, msh, ker, ptree, 1)
         end if
         Memory_far = Memory_far + Memory_tmp

      elseif (blocks%style == 1) then
         call Full_construction(blocks, msh, ker, stats, option, ptree, Memory_tmp)
         Memory_near = Memory_near + Memory_tmp
      elseif (blocks%style == 4) then
         do ii = 1, 2
         do jj = 1, 2
            blocks_son => blocks%sons(ii, jj)
            if (IOwnPgrp(ptree, blocks_son%pgno))call Hmat_block_construction(blocks_son, Memory_far, Memory_near, option, stats, msh, ker, ptree)
         enddo
         enddo
      endif
      endif
! t2=MPI_Wtime()
! if(blocks%level==1)write(*,*)blocks%row_group,blocks%col_group,ptree%MyID,t2-t1
      return

   end subroutine Hmat_block_construction

!!!!!!! check error of BF construction using parallel element extraction
   subroutine BF_checkError(blocks, option, msh, ker, stats, ptree, bpackflag, verbosity,error)
      implicit none

      integer bpackflag ! whether blocks is treated as one offdiagonal block in a hierarhical matrix or one standalone block
      type(Hoption)::option
      type(Hstat)::stats
      ! type(Bmatrix)::bmat
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(intersect), allocatable::inters(:)
      real(kind=8)::n1, n2, n3, n4
      integer Ntest, passflag, verbosity
      integer nsproc1, nsproc2, nprow, npcol, nprow1D, npcol1D, myrow, mycol, nprow1, npcol1, myrow1, mycol1, nprow2, npcol2, myrow2, mycol2, myArows, myAcols, rank1, rank2, ierr, MyID
      integer:: cridx, info
      integer, allocatable::rows(:), cols(:)
      DT, allocatable:: Vin(:, :), Vout1(:, :), Vout2(:, :), Vinter(:, :), Mat(:, :)
      integer::descVin(9), descVout(9), descVinter(9), descFull(9), descButterflyU(9), descButterflyV(9)
      integer N, M, i, j, ii, jj, nn, pp, myi, myj, iproc, jproc, rmax
      integer edge_n, edge_m, rank
      real(kind=8):: fnorm1, fnorm0, rtemp1 = 0, rtemp0 = 0
      real(kind=8):: a, v1, v2, v3, vtmp
      real(kind=8),optional:: error
      DT:: value1, value2, value3
      type(list)::lstr, lstc, lst, lstblk
      type(nod), pointer::cur, curr, curc, curri, curci
      class(*), pointer::ptr, ptrr, ptrc, ptrri, ptrci
      integer::head, tail, idx, pgno, ctxt, nr_loc, nc_loc
      type(matrixblock), pointer::blocks
      integer num_blocks, idx_row, idx_col, idx_dat
      integer, allocatable:: allrows(:), allcols(:), pmaps(:, :)
      integer, allocatable::datidx(:), colidx(:), rowidx(:), pgidx(:)
      DT, target, allocatable::alldat_loc(:)
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
      ! call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
      ! nprow = ptree%pgrp(pgno)%nprow
      ! npcol = ptree%pgrp(pgno)%npcol

      nproc = ptree%pgrp(blocks%pgno)%nproc
      Npmap = min(Ninter, nproc)
      npavr = nproc/Npmap
      allocate (pmaps(Npmap, 3))
      do nn = 1, Npmap
         nprow = floor_safe(sqrt(dble(npavr)))
         npcol = floor_safe(npavr/dble(nprow))
         pmaps(nn, 1) = 1   ! nprow   ! this makes sure the intersection is on 1 processor, this makes it easier for cpp user-defined extraction function
         pmaps(nn, 2) = 1   ! npcol
         pmaps(nn, 3) = (nn - 1)*npavr + ptree%pgrp(blocks%pgno)%head
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
            call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, Main_ID, ptree%pgrp(blocks%pgno)%Comm, ierr)
            if(bpackflag==1)then
               allrows(idx_row + 1) = max(floor_safe(blocks%M*a), 1) + blocks%headm -1
            else
               allrows(idx_row + 1) = max(floor_safe(blocks%M*a), 1)
            endif
            ! allrows(idx_row+1)=msh%basis_group(2**level+nn-1)%head+ii-1
            idx_row = idx_row + 1
         enddo

         do ii = 1, nc
            call random_number(a)
            call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, Main_ID, ptree%pgrp(blocks%pgno)%Comm, ierr)
            if(bpackflag==1)then
               allcols(idx_col + 1) = max(floor_safe(blocks%N*a), 1) + blocks%headn -1
            else
               allcols(idx_col + 1) = max(floor_safe(blocks%N*a), 1) + blocks%M
            endif
            ! allcols(idx_col+1)=msh%basis_group(2**level+1-(nn-1))%head+ii-1
            idx_col = idx_col + 1
         enddo
      enddo

      allocate (alldat_loc(idx_dat))
      if (idx_dat > 0) alldat_loc = 0

      n1 = MPI_Wtime()
      call BF_ExtractElement(blocks, option, msh, stats, ptree, Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps)


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

      n2 = MPI_Wtime()
      vtmp = v1
      call MPI_ALLREDUCE(vtmp, v1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(blocks%pgno)%Comm  , ierr)
      vtmp = v2
      call MPI_ALLREDUCE(vtmp, v2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(blocks%pgno)%Comm  , ierr)
      vtmp = v3
      call MPI_ALLREDUCE(vtmp, v3, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(blocks%pgno)%Comm  , ierr)

      if (ptree%MyID - ptree%pgrp(blocks%pgno)%head == Main_ID .and. verbosity >= 0) write (*, '(A25,Es14.7,Es14.7,A6,Es9.2,A7,Es9.2)') 'BF_CheckError: fnorm:', sqrt(v1), sqrt(v2), ' acc: ', sqrt(v3/v1), ' time: ', n2 - n1
      if(present(error)) error=sqrt(v3/v1)

      !stop

   end subroutine BF_checkError



!!!!!!! check error of BP construction using parallel element extraction
   subroutine BP_checkError(BP, option, msh, ker, stats, ptree, bpackflag, verbosity,error)

      implicit none

      integer bpackflag ! whether blocks is treated as one offdiagonal block in a hierarhical matrix or one standalone block
      type(Hoption)::option
      type(Hstat)::stats
      ! type(Bmatrix)::bmat
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(intersect), allocatable::inters(:)
      real(kind=8)::n1, n2, n3, n4
      integer Ntest, passflag, verbosity
      integer nsproc1, nsproc2, nprow, npcol, nprow1D, npcol1D, myrow, mycol, nprow1, npcol1, myrow1, mycol1, nprow2, npcol2, myrow2, mycol2, myArows, myAcols, rank1, rank2, ierr, MyID
      integer:: cridx, info
      integer, allocatable::rows(:), cols(:)
      DT, allocatable:: Vin(:, :), Vout1(:, :), Vout2(:, :), Vinter(:, :), Mat(:, :)
      integer::descVin(9), descVout(9), descVinter(9), descFull(9), descButterflyU(9), descButterflyV(9)
      integer N, M, i, j, ii, jj, nn, pp, myi, myj, iproc, jproc, rmax
      integer edge_n, edge_m, rank
      real(kind=8):: fnorm1, fnorm0, rtemp1 = 0, rtemp0 = 0
      real(kind=8):: a, v1, v2, v3, vtmp
      real(kind=8),optional:: error
      DT:: value1, value2, value3
      type(list)::lstr, lstc, lst, lstblk
      type(nod), pointer::cur, curr, curc, curri, curci
      class(*), pointer::ptr, ptrr, ptrc, ptrri, ptrci
      integer::head, tail, idx, pgno, ctxt, nr_loc, nc_loc
      type(blockplus)::BP
      type(matrixblock), pointer::blocks
      integer num_blocks, idx_row, idx_col, idx_dat
      integer, allocatable:: allrows(:), allcols(:), pmaps(:, :)
      integer, allocatable::datidx(:), colidx(:), rowidx(:), pgidx(:)
      DT, target, allocatable::alldat_loc(:)
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
      ! call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
      ! nprow = ptree%pgrp(pgno)%nprow
      ! npcol = ptree%pgrp(pgno)%npcol


      blocks => BP%ll(1)%matrices_block(1)
      nproc = ptree%pgrp(blocks%pgno)%nproc
      Npmap = min(Ninter, nproc)
      npavr = nproc/Npmap
      allocate (pmaps(Npmap, 3))
      do nn = 1, Npmap
         nprow = floor_safe(sqrt(dble(npavr)))
         npcol = floor_safe(npavr/dble(nprow))
         pmaps(nn, 1) = 1   ! nprow   ! this makes sure the intersection is on 1 processor, this makes it easier for cpp user-defined extraction function
         pmaps(nn, 2) = 1   ! npcol
         pmaps(nn, 3) = (nn - 1)*npavr + ptree%pgrp(blocks%pgno)%head
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
            call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, Main_ID, ptree%pgrp(blocks%pgno)%Comm, ierr)
            if(bpackflag==1)then
               allrows(idx_row + 1) = max(floor_safe(blocks%M*a), 1) + blocks%headm -1
            else
               allrows(idx_row + 1) = max(floor_safe(blocks%M*a), 1)
            endif
            ! allrows(idx_row+1)=msh%basis_group(2**level+nn-1)%head+ii-1
            idx_row = idx_row + 1
         enddo

         do ii = 1, nc
            call random_number(a)
            call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, Main_ID, ptree%pgrp(blocks%pgno)%Comm, ierr)
            if(bpackflag==1)then
               allcols(idx_col + 1) = max(floor_safe(blocks%N*a), 1) + blocks%headn -1
            else
               allcols(idx_col + 1) = max(floor_safe(blocks%N*a), 1) + blocks%M
            endif
            ! allcols(idx_col+1)=msh%basis_group(2**level+1-(nn-1))%head+ii-1
            idx_col = idx_col + 1
         enddo
      enddo

      allocate (alldat_loc(idx_dat))
      if (idx_dat > 0) alldat_loc = 0

      n1 = MPI_Wtime()
      select case (option%format)
      case (HODLR)
         call BF_ExtractElement(blocks, option, msh, stats, ptree, Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps)
      case (HMAT,BLR)
         call BP_ExtractElement(BP, option, msh, stats, ptree, Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps)
      case (HSS)
         call BP_ExtractElement(BP, option, msh, stats, ptree, Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps)
      end select
      n2 = MPI_Wtime()

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
      vtmp=v1
      call MPI_ALLREDUCE(vtmp, v1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(BP%pgno)%Comm  , ierr)
      vtmp=v2
      call MPI_ALLREDUCE(vtmp, v2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(BP%pgno)%Comm  , ierr)
      vtmp=v3
      call MPI_ALLREDUCE(vtmp, v3, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(BP%pgno)%Comm  , ierr)

      if (ptree%MyID - ptree%pgrp(BP%pgno)%head == Main_ID .and. verbosity >= 0) write (*, '(A25,Es14.7,Es14.7,A6,Es9.2,A7,Es9.2)') 'BP_CheckError: fnorm:', sqrt(v1), sqrt(v2), ' acc: ', sqrt(v3/v1), ' time: ', n2 - n1
      if(present(error)) error=sqrt(v3/v1)

      !stop

   end subroutine BP_checkError



   subroutine BP_compress_entry(bplus, option, Memory, stats, msh, ker, ptree)

      implicit none

      type(blockplus)::bplus
      integer:: ii, ll, bb, ierr, pp
      real(kind=8) Memory, rtemp, error
      integer:: level_butterfly, level_BP, levelm, groupm_start, Nboundall, Ninadmissible, statflag, knn_tmp
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(matrixblock), pointer::blocks

      Memory = 0
      do ll = 1, bplus%Lplus
         bplus%LL(ll)%rankmax = 0
         statflag = 0
         if (ll == 1 .or. option%bp_cnt_lr == 1) statflag = 1  !!! only record the rank of the top-layer butterfly in a bplus
         do bb = 1, bplus%LL(ll)%Nbound
            if (IOwnPgrp(ptree, bplus%LL(ll)%matrices_block(bb)%pgno)) then
               if (bplus%LL(ll)%matrices_block(bb)%style == 1) then
                  call Full_construction(bplus%LL(ll)%matrices_block(bb), msh, ker, stats, option, ptree, rtemp)
                  Memory = Memory + rtemp
               else

                  level_butterfly = bplus%LL(ll)%matrices_block(bb)%level_butterfly
                  level_BP = bplus%level

                  bplus%LL(ll)%matrices_block(bb)%level_half = BF_Switchlevel(bplus%LL(ll)%matrices_block(bb)%level_butterfly, option%pat_comp)
                  levelm = bplus%LL(ll)%matrices_block(bb)%level_half
                  groupm_start = bplus%LL(ll)%matrices_block(1)%row_group*2**levelm
                  Nboundall = 0
                  Ninadmissible = 0
                  if (allocated(bplus%LL(ll + 1)%boundary_map)) then
                     Nboundall = size(bplus%LL(ll + 1)%boundary_map, 1)
                     Ninadmissible = size(bplus%LL(ll + 1)%boundary_map, 2)
                  endif
                  if (option%forwardN15flag == 1) then
                     knn_tmp = option%knn
                     option%knn=0
                     call BF_compress_N15(bplus%LL(ll)%matrices_block(bb), bplus%LL(ll + 1)%boundary_map, Nboundall, Ninadmissible, groupm_start, option, rtemp, stats, msh, ker, ptree, statflag)
                     option%knn=knn_tmp
                     call BF_sym2asym(bplus%LL(ll)%matrices_block(bb))

                     ! Move singular values to leftmost factor
                     if(bplus%LL(ll)%matrices_block(bb)%level_butterfly>0)then
                        call BF_ChangePattern(bplus%LL(ll)%matrices_block(bb), 3, 1, stats, ptree)
                        call BF_MoveSingular_Ker(bplus%LL(ll)%matrices_block(bb), 'N', 1, bplus%LL(ll)%matrices_block(bb)%level_butterfly, ptree, stats, option%tol_comp)
                        call BF_ChangePattern(bplus%LL(ll)%matrices_block(bb), 1, 3, stats, ptree)
                     endif
                  elseif (option%forwardN15flag == 2) then
                        call BF_compress_NlogN(bplus%LL(ll)%matrices_block(bb), bplus%LL(ll + 1)%boundary_map, Nboundall, Ninadmissible, groupm_start, option, rtemp, stats, msh, ker, ptree, statflag)
                        allocate(blocks)
                        call BF_copy('N', bplus%LL(ll)%matrices_block(bb), blocks)
                        call BF_checkError(blocks, option, msh, ker, stats, ptree, 1, -1, error)
                        deallocate(blocks)
                        if(error>5*option%tol_comp)then
                           pp = ptree%myid - ptree%pgrp(bplus%LL(ll)%matrices_block(bb)%pgno)%head + 1
                           if(option%verbosity>=0 .and. pp==1)write(*,*)'warning: BF',bplus%LL(ll)%matrices_block(bb)%row_group,bplus%LL(ll)%matrices_block(bb)%col_group,' error ',error,', redo the compression with the N15 algorithm'
                           call BF_delete(bplus%LL(ll)%matrices_block(bb),0)
                           knn_tmp = option%knn
                           option%knn=0
                           call BF_compress_N15(bplus%LL(ll)%matrices_block(bb), bplus%LL(ll + 1)%boundary_map, Nboundall, Ninadmissible, groupm_start, option, rtemp, stats, msh, ker, ptree, statflag)
                           option%knn=knn_tmp
                           call BF_sym2asym(bplus%LL(ll)%matrices_block(bb))
                           ! Move singular values to leftmost factor
                           if(bplus%LL(ll)%matrices_block(bb)%level_butterfly>0)then
                              call BF_ChangePattern(bplus%LL(ll)%matrices_block(bb), 3, 1, stats, ptree)
                              call BF_MoveSingular_Ker(bplus%LL(ll)%matrices_block(bb), 'N', 1, bplus%LL(ll)%matrices_block(bb)%level_butterfly, ptree, stats, option%tol_comp)
                              call BF_ChangePattern(bplus%LL(ll)%matrices_block(bb), 1, 3, stats, ptree)
                           endif
                        endif
                  else
                     call BF_compress_NlogN(bplus%LL(ll)%matrices_block(bb), bplus%LL(ll + 1)%boundary_map, Nboundall, Ninadmissible, groupm_start, option, rtemp, stats, msh, ker, ptree, statflag)
                  end if

                  Memory = Memory + rtemp
                  bplus%LL(ll)%rankmax = max(bplus%LL(ll)%rankmax, bplus%LL(ll)%matrices_block(bb)%rankmax)
               endif
            endif
         end do
      end do

      ! !!!!!!! check error
      if (option%ErrFillFull == 1) call Bplus_CheckError_Full(bplus, option, msh, ker, stats, ptree)
      ! !!!!!!! check error

      ! if(bplus%LL(1)%matrices_block(1)%level==1)write(*,*)bplus%LL(1)%rankmax,'dfdfdfdf'

      return

   end subroutine BP_compress_entry





   subroutine BP_MD_compress_entry(Ndim, bplus, option, Memory_comp, Memory_direct, stats, msh, ker, ptree, logn_level_flag)

      implicit none
      integer Ndim
      type(blockplus_MD)::bplus
      integer:: ii, ll, bb, ierr, pp, rep_bb, trans_reuse_level, trans_nrep_level
      real(kind=8) Memory_comp, Memory_direct, Memory_perlevel, rtemp,rtemp1, error
      integer:: level_butterfly, level_BP, levelm, statflag, knn_tmp
      integer, allocatable:: trans_rep_candidates(:)
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh(Ndim)
      type(kernelquant)::ker
      type(proctree)::ptree
      type(matrixblock_MD), pointer::blocks
      integer groupm_start(Ndim), groupm_start0(Ndim),Nboundall,Ninadmissible
      integer logn_level_flag
      real(kind=8):: n1,n2
      integer level_butterfly_tmp, ntrial, rmax_start
      integer(kind=8) dist_bb

      Memory_comp = 0
      Memory_direct = 0
      groupm_start0=bplus%LL(1)%matrices_block(1)%row_group
      do ll = 1, bplus%Lplus
         Memory_perlevel = 0
         n1 = MPI_Wtime()
         bplus%LL(ll)%rankmax = 0
         trans_reuse_level = 0
         trans_nrep_level = 0
         if (option%trans_invariant /= 0 .and. logn_level_flag == 1) then
            allocate(trans_rep_candidates(max(1, bplus%LL(ll)%Nbound_loc)))
         endif
         statflag = 0
         if (ll == 1 .or. option%bp_cnt_lr == 1) statflag = 1  !!! only record the rank of the top-layer butterfly in a bplus
         do bb = 1, bplus%LL(ll)%Nbound_loc
            if (IOwnPgrp(ptree, bplus%LL(ll)%matrices_block(bb)%pgno)) then
               if(option%verbosity>=2 .and. ptree%MyID==ptree%pgrp(bplus%LL(ll)%matrices_block(bb)%pgno)%head)then
                  write(*,*)'Start compressing BF_MD',bb,' out of #local ', bplus%LL(ll)%Nbound_loc,'on rank',ptree%MyID, 'at level',ll
               endif
               if (bplus%LL(ll)%matrices_block(bb)%style == 1) then
                  rep_bb = 0
                  if (option%trans_invariant /= 0 .and. logn_level_flag == 1) then
                     ! call assert(option%use_zfp == 0 .and. option%use_qtt == 0, "trans_invariant currently requires use_zfp=0 and use_qtt=0")
                     rep_bb = BP_MD_find_trans_rep(Ndim, bplus, ll, bb, ptree, option%trans_invariant, &
                        trans_rep_candidates, trans_nrep_level)
                  endif
                  if (rep_bb > 0) then
                     call BP_MD_alias_trans_data(Ndim, bplus%LL(ll)%matrices_block(rep_bb), bplus%LL(ll)%matrices_block(bb), rep_bb)
                     trans_reuse_level = trans_reuse_level + 1
                  else
                     call Full_construction_MD(Ndim, bplus%LL(ll)%matrices_block(bb), msh, ker, stats, option, ptree)
                     call BF_MD_ComputeMemory(Ndim, bplus%LL(ll)%matrices_block(bb), rtemp,rtemp1)
                     Memory_direct = Memory_direct + rtemp
                     Memory_perlevel = Memory_perlevel + rtemp
                     if (option%trans_invariant /= 0 .and. logn_level_flag == 1 .and. &
                        ptree%pgrp(bplus%LL(ll)%matrices_block(bb)%pgno)%nproc == 1) then
                        trans_nrep_level = trans_nrep_level + 1
                        trans_rep_candidates(trans_nrep_level) = bb
                     endif
                  endif
               elseif (bplus%LL(ll)%matrices_block(bb)%style == 2) then

                  level_butterfly = bplus%LL(ll)%matrices_block(bb)%level_butterfly
                  bplus%LL(ll)%matrices_block(bb)%level_half = BF_Switchlevel(bplus%LL(ll)%matrices_block(bb)%level_butterfly, option%pat_comp)
                  if(logn_level_flag==0)then
                     levelm = bplus%LL(ll)%matrices_block(bb)%level_half
                  else
                     levelm = 1
                  endif
                  groupm_start = groupm_start0*2**levelm
                  Nboundall = 0
                  Ninadmissible = 0
                  if (allocated(bplus%LL(ll + 1)%boundary_map) .and. logn_level_flag==0) then
                     Nboundall = size(bplus%LL(ll + 1)%boundary_map, 1)
                     Nboundall = NINT(exp(log(dble(Nboundall)) / dble(Ndim)))
                     Ninadmissible = size(bplus%LL(ll + 1)%boundary_map, 2)
                  endif
                  call assert(option%forwardN15flag == 0, "only forwardN15flag == 0 can be used for tensor butterfly")

                  rep_bb = 0
                  if (option%trans_invariant /= 0 .and. logn_level_flag == 1) then
                     ! call assert(option%use_zfp == 0 .and. option%use_qtt == 0, "trans_invariant currently requires use_zfp=0 and use_qtt=0")
                     rep_bb = BP_MD_find_trans_rep(Ndim, bplus, ll, bb, ptree, option%trans_invariant, &
                        trans_rep_candidates, trans_nrep_level)
                  endif
                  if (rep_bb > 0) then
                     call BP_MD_alias_trans_data(Ndim, bplus%LL(ll)%matrices_block(rep_bb), bplus%LL(ll)%matrices_block(bb), rep_bb)
                     trans_reuse_level = trans_reuse_level + 1
                  else
                     blocks=>bplus%LL(ll)%matrices_block(bb)
                     call BF_MD_compress_N(Ndim,blocks, bplus%LL(ll + 1)%boundary_map, Nboundall, Ninadmissible, groupm_start, option, rtemp, stats, msh, ker, ptree,statflag, 1)
                     call BF_MD_get_rank(Ndim, blocks, ptree)

                     if(blocks%level_butterfly==0)then ! this makes LR compression more robust when rmax is too small
                        rmax_start = option%rmax
                        ntrial = option%itermax
                        do while (ntrial > 1)
                           call BF_MD_get_rank(Ndim, blocks, ptree)
                           if(blocks%rankmin>=option%rmax)then
                              option%rmax = option%rmax*option%rankrate
                              call BF_MD_delete(Ndim, blocks, 0)
                              call BF_MD_compress_N(Ndim,blocks, bplus%LL(ll + 1)%boundary_map, Nboundall, Ninadmissible, groupm_start, option, rtemp, stats, msh, ker, ptree,statflag, 1)
                           else
                              exit
                           endif
                           ntrial = ntrial - 1
                        enddo
                        option%rmax = rmax_start
                     endif

                     dist_bb = sum(int(bplus%LL(ll)%matrices_block(bb)%row_group - bplus%LL(ll)%matrices_block(bb)%col_group, kind=8)**2)
                     call MPI_ALLREDUCE(rtemp, rtemp1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(bplus%LL(ll)%matrices_block(bb)%pgno)%Comm, ierr)
                     if(option%verbosity>=2 .and. ptree%MyID==ptree%pgrp(bplus%LL(ll)%matrices_block(bb)%pgno)%head)write(*,*)'Compressed block:', bb, 'Memory used:', rtemp1,'dist',dist_bb,'lb',blocks%level_butterfly,'min_rank',blocks%rankmin,'max_rank',blocks%rankmax

                     Memory_comp = Memory_comp + rtemp
                     Memory_perlevel = Memory_perlevel + rtemp
                     if (option%trans_invariant /= 0 .and. logn_level_flag == 1 .and. &
                        ptree%pgrp(bplus%LL(ll)%matrices_block(bb)%pgno)%nproc == 1) then
                        trans_nrep_level = trans_nrep_level + 1
                        trans_rep_candidates(trans_nrep_level) = bb
                     endif
                  endif
                  bplus%LL(ll)%rankmax = max(bplus%LL(ll)%rankmax, bplus%LL(ll)%matrices_block(bb)%rankmax)
               endif
            endif
         end do
         if(logn_level_flag==0)then
            levelm = bplus%LL(ll)%matrices_block(1)%level_half
         else
            levelm = 1
         endif
         n2 = MPI_Wtime()
         groupm_start0 = groupm_start0*2**levelm
         if(option%verbosity>=1 .and. ptree%MyID==ptree%pgrp(bplus%LL(1)%matrices_block(1)%pgno)%head)then
            write(*,*)'Finishing level ', ll, 'in BP_MD_compress_entry, level_butterfly:', &
               bplus%LL(ll)%level_butterfly, 'rankmax at this level:', bplus%LL(ll)%rankmax, 'time:', n2-n1, 'mem(rank0):', Memory_perlevel
            if (option%trans_invariant /= 0 .and. logn_level_flag == 1) then
               write(*,*)'trans_invariant unique/compact local blocks at this level:', &
                  trans_nrep_level, bplus%LL(ll)%trans_ndup
            endif
         endif
         if (allocated(trans_rep_candidates)) deallocate(trans_rep_candidates)
      end do

      ! ! !!!!!!! check error
      ! if (option%ErrFillFull == 1) call Bplus_CheckError_Full(bplus, option, msh, ker, stats, ptree)
      ! ! !!!!!!! check error

      ! if(bplus%LL(1)%matrices_block(1)%level==1)write(*,*)bplus%LL(1)%rankmax,'dfdfdfdf'

      return

   end subroutine BP_MD_compress_entry


   integer function BP_MD_find_trans_rep(Ndim, bplus, ll, bb, ptree, trans_invariant, rep_candidates, nrep)

      implicit none
      integer Ndim, ll, bb, rr, trans_invariant, kk, nsearch
      integer, optional, intent(in)::rep_candidates(:), nrep
      integer(kind=8) dist_rr, dist_bb
      type(blockplus_MD)::bplus
      type(proctree)::ptree

      BP_MD_find_trans_rep = 0
      if (ptree%pgrp(bplus%LL(ll)%matrices_block(bb)%pgno)%nproc /= 1) return
      if (.not. allocated(bplus%LL(ll)%matrices_block(bb)%row_group)) return
      if (.not. allocated(bplus%LL(ll)%matrices_block(bb)%col_group)) return

      nsearch = bb - 1
      if (present(rep_candidates) .and. present(nrep)) nsearch = nrep
      do kk = 1, nsearch
         if (present(rep_candidates) .and. present(nrep)) then
            rr = rep_candidates(kk)
         else
            rr = kk
         endif
         if (.not. IOwnPgrp(ptree, bplus%LL(ll)%matrices_block(rr)%pgno)) cycle
         if (ptree%pgrp(bplus%LL(ll)%matrices_block(rr)%pgno)%nproc /= 1) cycle
         if (bplus%LL(ll)%matrices_block(rr)%trans_invariant_dup == 1) cycle
         if (bplus%LL(ll)%matrices_block(rr)%style /= bplus%LL(ll)%matrices_block(bb)%style) cycle
         if (bplus%LL(ll)%matrices_block(rr)%level /= bplus%LL(ll)%matrices_block(bb)%level) cycle
         if (bplus%LL(ll)%matrices_block(rr)%level_butterfly /= bplus%LL(ll)%matrices_block(bb)%level_butterfly) cycle
         if (.not. allocated(bplus%LL(ll)%matrices_block(rr)%row_group)) cycle
         if (.not. allocated(bplus%LL(ll)%matrices_block(rr)%col_group)) cycle
         if (trans_invariant == 1) then
            if (.not. all((bplus%LL(ll)%matrices_block(rr)%row_group - bplus%LL(ll)%matrices_block(rr)%col_group) == &
               (bplus%LL(ll)%matrices_block(bb)%row_group - bplus%LL(ll)%matrices_block(bb)%col_group))) cycle
         elseif (trans_invariant == 2) then
            dist_rr = sum(int(bplus%LL(ll)%matrices_block(rr)%row_group - bplus%LL(ll)%matrices_block(rr)%col_group, kind=8)**2)
            dist_bb = sum(int(bplus%LL(ll)%matrices_block(bb)%row_group - bplus%LL(ll)%matrices_block(bb)%col_group, kind=8)**2)
            if (dist_rr /= dist_bb) cycle
            if (.not. BP_MD_trans2_compatible(Ndim, bplus%LL(ll)%matrices_block(rr), bplus%LL(ll)%matrices_block(bb))) cycle
         else
            cycle
         endif
         if (trans_invariant /= 2) then
            if (allocated(bplus%LL(ll)%matrices_block(rr)%M) .and. allocated(bplus%LL(ll)%matrices_block(bb)%M)) then
               if (.not. all(bplus%LL(ll)%matrices_block(rr)%M == bplus%LL(ll)%matrices_block(bb)%M)) cycle
            endif
            if (allocated(bplus%LL(ll)%matrices_block(rr)%N) .and. allocated(bplus%LL(ll)%matrices_block(bb)%N)) then
               if (.not. all(bplus%LL(ll)%matrices_block(rr)%N == bplus%LL(ll)%matrices_block(bb)%N)) cycle
            endif
         endif
         if (bplus%LL(ll)%matrices_block(bb)%style == 1) then
#if HAVE_FFTW
            if (.not. associated(bplus%LL(ll)%matrices_block(rr)%fullmat) .and. &
                .not. associated(bplus%LL(ll)%matrices_block(rr)%FullmatFFT)) cycle
#else
            if (.not. associated(bplus%LL(ll)%matrices_block(rr)%fullmat)) cycle
#endif
         elseif (bplus%LL(ll)%matrices_block(bb)%style == 2) then
            if (bplus%LL(ll)%matrices_block(rr)%level_half /= bplus%LL(ll)%matrices_block(bb)%level_half) cycle
            if (.not. allocated(bplus%LL(ll)%matrices_block(rr)%nr_m)) cycle
            if (.not. allocated(bplus%LL(ll)%matrices_block(rr)%nc_m)) cycle
         endif
         BP_MD_find_trans_rep = rr
         return
      enddo

   end function BP_MD_find_trans_rep


   logical function BP_MD_trans2_compatible(Ndim, rep, blk)

      implicit none
      integer Ndim
      type(matrixblock_MD)::rep, blk
      integer used_rep(Ndim)
      integer dim_blk, dim_rep, delta_blk, delta_rep
      logical found

      BP_MD_trans2_compatible = .false.
      if (.not. allocated(rep%row_group) .or. .not. allocated(rep%col_group)) return
      if (.not. allocated(blk%row_group) .or. .not. allocated(blk%col_group)) return

      used_rep = 0
      do dim_blk = 1, Ndim
         delta_blk = blk%row_group(dim_blk) - blk%col_group(dim_blk)
         found = .false.
         do dim_rep = 1, Ndim
            if (used_rep(dim_rep) == 1) cycle
            delta_rep = rep%row_group(dim_rep) - rep%col_group(dim_rep)
            if (abs(delta_rep) /= abs(delta_blk)) cycle
            if (allocated(rep%M) .and. allocated(blk%M)) then
               if (rep%M(dim_rep) /= blk%M(dim_blk)) cycle
            endif
            if (allocated(rep%N) .and. allocated(blk%N)) then
               if (rep%N(dim_rep) /= blk%N(dim_blk)) cycle
            endif
            if (allocated(rep%M_loc) .and. allocated(blk%M_loc)) then
               if (rep%M_loc(dim_rep) /= blk%M_loc(dim_blk)) cycle
            endif
            if (allocated(rep%N_loc) .and. allocated(blk%N_loc)) then
               if (rep%N_loc(dim_rep) /= blk%N_loc(dim_blk)) cycle
            endif
            used_rep(dim_rep) = 1
            found = .true.
            exit
         enddo
         if (.not. found) return
      enddo

      BP_MD_trans2_compatible = .true.

   end function BP_MD_trans2_compatible


   subroutine BP_MD_alias_trans_data(Ndim, rep, blk, rep_idx)

      implicit none
      integer Ndim
      integer, optional:: rep_idx
      type(matrixblock_MD), target::rep
      type(matrixblock_MD), target::blk

      blk%trans_invariant_dup = 1
      if (present(rep_idx)) blk%trans_rep_idx = rep_idx
      blk%trans_rep => rep
      blk%rankmax = rep%rankmax
      blk%rankmin = rep%rankmin
      blk%level_half = rep%level_half

      ! Duplicate blocks retain only their block-specific geometry/layout.
      ! Trans-invariant MVP groups apply rep's numerical factors directly.
      if (blk%style == 1) then
         blk%fullmat => rep%fullmat
         if (associated(rep%FullmatFFT)) blk%FullmatFFT => rep%FullmatFFT
      endif

   end subroutine BP_MD_alias_trans_data


   subroutine BP_MD_copy_iarray(src, dst)

      implicit none
      integer, allocatable, intent(in)::src(:)
      integer, allocatable, intent(inout)::dst(:)

      if (allocated(dst)) deallocate(dst)
      if (allocated(src)) then
         allocate(dst(size(src, 1)))
         dst = src
      endif

   end subroutine BP_MD_copy_iarray


   subroutine BP_MD_alias_uv(src, dst, Ndim)

      implicit none
      integer Ndim, ii, dim_i
      type(butterfly_UV_MD), intent(in)::src
      type(butterfly_UV_MD), intent(inout)::dst

      dst%num_blk = src%num_blk
      dst%nblk_loc = src%nblk_loc
      dst%idx = src%idx
      dst%inc = src%inc
      if (allocated(src%blocks)) then
         allocate(dst%blocks(size(src%blocks, 1), size(src%blocks, 2)))
         do dim_i = 1, Ndim
            do ii = 1, size(src%blocks, 1)
               call BP_MD_alias_butterfly_matrix(src%blocks(ii, dim_i), dst%blocks(ii, dim_i))
            enddo
         enddo
      endif

   end subroutine BP_MD_alias_uv


   subroutine BP_MD_alias_kerl(src, dst)

      implicit none
      integer ii, jj, kk
      type(butterfly_kerl_MD), intent(in)::src
      type(butterfly_kerl_MD), intent(inout)::dst

      call BP_MD_copy_iarray(src%num_col, dst%num_col)
      call BP_MD_copy_iarray(src%num_row, dst%num_row)
      call BP_MD_copy_iarray(src%nc, dst%nc)
      call BP_MD_copy_iarray(src%nr, dst%nr)
      call BP_MD_copy_iarray(src%idx_c, dst%idx_c)
      call BP_MD_copy_iarray(src%idx_r, dst%idx_r)
      call BP_MD_copy_iarray(src%inc_c, dst%inc_c)
      call BP_MD_copy_iarray(src%inc_r, dst%inc_r)
      if (allocated(src%blocks)) then
         allocate(dst%blocks(size(src%blocks, 1), size(src%blocks, 2), size(src%blocks, 3)))
         do kk = 1, size(src%blocks, 3)
            do jj = 1, size(src%blocks, 2)
               do ii = 1, size(src%blocks, 1)
                  call BP_MD_alias_butterfly_matrix(src%blocks(ii, jj, kk), dst%blocks(ii, jj, kk))
               enddo
            enddo
         enddo
      endif

   end subroutine BP_MD_alias_kerl


   subroutine BP_MD_alias_butterfly_matrix(src, dst)

      implicit none
      type(butterflymatrix), intent(in)::src
      type(butterflymatrix), intent(inout)::dst

      dst%ndim = src%ndim
      call BP_MD_copy_iarray(src%dims_m, dst%dims_m)
      call BP_MD_copy_iarray(src%dims_n, dst%dims_n)
      if (associated(src%matrix)) dst%matrix => src%matrix

   end subroutine BP_MD_alias_butterfly_matrix





!!!!!!! check error of BF construction using parallel element extraction
   subroutine BF_MD_checkError(Ndim, blocks, option, msh, ker, stats, ptree, bpackflag, verbosity,error)

      implicit none
      integer Ndim
      integer bpackflag ! whether blocks is treated as one offdiagonal block in a hierarhical matrix or one standalone block
      type(Hoption)::option
      type(Hstat)::stats
      ! type(Bmatrix)::bmat
      type(mesh)::msh(Ndim)
      type(kernelquant)::ker
      type(proctree)::ptree
      type(intersect), allocatable::inters(:)
      real(kind=8)::n1, n2, n3, n4
      integer Ntest, passflag, verbosity
      integer nsproc1, nsproc2, nprow, npcol, nprow1D, npcol1D, myrow, mycol, nprow1, npcol1, myrow1, mycol1, nprow2, npcol2, myrow2, mycol2, myArows, myAcols, rank1, rank2, ierr, MyID
      integer:: cridx, info, dim_i
      integer, allocatable::rows(:), cols(:)
      DT, allocatable:: Vin(:, :), Vout1(:, :), Vout2(:, :), Vinter(:, :), Mat(:, :)
      integer::descVin(9), descVout(9), descVinter(9), descFull(9), descButterflyU(9), descButterflyV(9)
      integer N, M, i, j, ii, jj, nn, pp, myi, myj, iproc, jproc, rmax
      integer edge_n, edge_m, rank
      real(kind=8):: fnorm1, fnorm0, rtemp1 = 0, rtemp0 = 0
      real(kind=8):: a, v1, v2, v3, vtmp
      real(kind=8),optional:: error
      DT:: value1, value2, value3
      type(list)::lstr, lstc, lst, lstblk
      type(nod), pointer::cur, curr, curc, curri, curci
      class(*), pointer::ptr, ptrr, ptrc, ptrri, ptrci
      integer::head, tail, idx, pgno, ctxt, nr_loc, nc_loc
      type(matrixblock_MD), pointer::blocks
      integer num_blocks, idx_row, idx_col, idx_dat
      integer, allocatable:: allrows(:), allcols(:), pmaps(:, :)
      integer, allocatable::datidx(:), colidx(:), rowidx(:), pgidx(:)
      DT, target, allocatable::alldat_loc(:)
      integer:: Ninter, dims_r, dims_c, dim_MD(Ndim*2), idx_MD(Ndim*2), idx_m(Ndim), idx_r_m(Ndim), idx_c_m(Ndim), ntot_loc, level, Npmap, nproc, npavr, np
      type(intersect_MD)::subtensor(1)
      type(intersect_MD),allocatable::inter_MD(:)
      integer level_butterfly,level_half,levelm,receiver, sender, bbm, vtmpi
      integer,allocatable::order(:),order_m(:)
      integer:: dims_row(Ndim),dims_col(Ndim), dims_one(Ndim),group_m(Ndim),group_n(Ndim),dims_1D(1)
      real(kind=8)::dis
      real(kind=8),allocatable::distances(:)

      level_butterfly=blocks%level_butterfly
      level_half=blocks%level_half
      levelm=level_half
      nproc = ptree%pgrp(blocks%pgno)%nproc
      pp = ptree%MyID - ptree%pgrp(blocks%pgno)%head + 1

      dims_r=2**level_half
      dims_c=2**(level_butterfly-level_half)


      dim_MD(1:Ndim) = dims_r
      dim_MD(1+Ndim:Ndim*2) = dims_c
      allocate(distances(product(dim_MD)))
      do bbm=1,product(dim_MD)
         call SingleIndexToMultiIndex(Ndim*2, dim_MD, bbm, idx_MD)
         group_m = blocks%row_group
         group_m = group_m*2**levelm - 1 + idx_MD(1:Ndim)
         group_n = blocks%col_group
         group_n = group_n*2**(level_butterfly-levelm) - 1 + idx_MD(1+Ndim:Ndim*2)
         if (option%nogeo == 1)then
            if(ALL(group_m==group_n))then
               distances(bbm)=0
            else
               distances(bbm)=1
            endif
         else
            dis = 0d0
            do i = 1, Ndim
               dis = dis + (msh(i)%basis_group(group_m(i))%center(1) - msh(i)%basis_group(group_n(i))%center(1))**2
            enddo
            dis = dis**(1d0/2d0)
            distances(bbm)= dis
         endif
      enddo
      allocate(order_m(product(dim_MD)))
      call quick_sort(distances, order_m, product(dim_MD))
      call MPI_Bcast(order_m, product(dim_MD), MPI_INTEGER, 0, ptree%Comm, ierr)


      Ninter = min(4, min(dims_r**Ndim,dims_c**Ndim))
      allocate(inter_MD(Ninter))
      do nn=1,Ninter
         allocate(inter_MD(nn)%nr(Ndim))
         allocate(inter_MD(nn)%rows(Ndim))
         bbm=order_m(nn)
         call SingleIndexToMultiIndex(Ndim*2, dim_MD, bbm, idx_MD)

         do dim_i=1,Ndim
            ! call random_number(a)
            ! call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, Main_ID, ptree%pgrp(blocks%pgno)%Comm, ierr)
            ! idx_r_m(dim_i) = max(floor_safe(dims_r*a), 1)
            idx_r_m(dim_i) = idx_MD(dim_i)
            group_m = blocks%row_group
            group_m(dim_i) = group_m(dim_i)*2**levelm - 1 + idx_r_m(dim_i)
            inter_MD(nn)%nr(dim_i) = msh(dim_i)%basis_group(group_m(dim_i))%tail - msh(dim_i)%basis_group(group_m(dim_i))%head + 1
            allocate(inter_MD(nn)%rows(dim_i)%dat(inter_MD(nn)%nr(dim_i)))
            do ii = 1, inter_MD(nn)%nr(dim_i)
               inter_MD(nn)%rows(dim_i)%dat(ii) = ii + msh(dim_i)%basis_group(group_m(dim_i))%head -1
            enddo
            allocate(order(inter_MD(nn)%nr(dim_i)))
            call quick_sort_int(inter_MD(nn)%rows(dim_i)%dat,order,inter_MD(nn)%nr(dim_i))
            inter_MD(nn)%rows(dim_i)%dat=inter_MD(nn)%rows(dim_i)%dat(order)
            deallocate(order)
         enddo

         allocate(inter_MD(nn)%nc(Ndim))
         allocate(inter_MD(nn)%cols(Ndim))
         do dim_i=1,Ndim
            ! call random_number(a)
            ! call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, Main_ID, ptree%pgrp(blocks%pgno)%Comm, ierr)
            ! idx_c_m(dim_i) = max(floor_safe(dims_c*a), 1)
            idx_c_m(dim_i) = idx_MD(dim_i+Ndim)
            group_n = blocks%col_group
            group_n(dim_i) = group_n(dim_i)*2**(level_butterfly-levelm) - 1 + idx_c_m(dim_i)
            inter_MD(nn)%nc(dim_i) = msh(dim_i)%basis_group(group_n(dim_i))%tail - msh(dim_i)%basis_group(group_n(dim_i))%head + 1
            allocate(inter_MD(nn)%cols(dim_i)%dat(inter_MD(nn)%nc(dim_i)))
            do ii = 1, inter_MD(nn)%nc(dim_i)
               inter_MD(nn)%cols(dim_i)%dat(ii) = ii + msh(dim_i)%basis_group(group_n(dim_i))%head -1
            enddo
            allocate(order(inter_MD(nn)%nc(dim_i)))
            call quick_sort_int(inter_MD(nn)%cols(dim_i)%dat,order,inter_MD(nn)%nc(dim_i))
            inter_MD(nn)%cols(dim_i)%dat=inter_MD(nn)%cols(dim_i)%dat(order)
            deallocate(order)
         enddo

         receiver=0
         sender=0
         idx_m=idx_r_m-blocks%idx_r_m+1
         if (ALL(idx_m >0) .and. ALL(idx_m <=blocks%nr_m))receiver=pp
         idx_m=idx_c_m-blocks%idx_c_m+1
         if (ALL(idx_m >0) .and. ALL(idx_m <=blocks%nc_m))sender=pp
         vtmpi = receiver
         call MPI_ALLREDUCE(vtmpi, receiver, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(blocks%pgno)%Comm , ierr)
         vtmpi = sender
         call MPI_ALLREDUCE(vtmpi, sender, 1, MPI_INTEGER, MPI_MAX, ptree%pgrp(blocks%pgno)%Comm , ierr)

         inter_MD(nn)%receiver=receiver
         inter_MD(nn)%sender=sender
         allocate(inter_MD(nn)%idx_r_m(Ndim))
         inter_MD(nn)%idx_r_m=idx_r_m
         allocate(inter_MD(nn)%idx_c_m(Ndim))
         inter_MD(nn)%idx_c_m=idx_c_m
         ! write(*,*)"inter nn",nn,idx_r_m,idx_c_m
         if(receiver==pp)then
            allocate(inter_MD(nn)%dat(product(inter_MD(nn)%nr),product(inter_MD(nn)%nc)))
         endif
      enddo
      deallocate(order_m)

      n1 = MPI_Wtime()
      call BF_MD_block_extraction(blocks, Ndim, Ninter, inter_MD, ptree, msh, stats, option)


      ! compare extracted values with element_Zmn
      v1 = 0
      v2 = 0
      v3 = 0
      do nn = 1, Ninter
         receiver = inter_MD(nn)%receiver
         if(receiver==pp)then
            if(option%verbosity>=2)write(*,*)'generating subtensor from entry evaluation',nn,'of',Ninter
            allocate(subtensor(1)%nr(Ndim))
            subtensor(1)%nr=inter_MD(nn)%nr
            allocate(subtensor(1)%nc(Ndim))
            subtensor(1)%nc=inter_MD(nn)%nc
            allocate(subtensor(1)%rows(Ndim))
            allocate(subtensor(1)%cols(Ndim))
            do dim_i=1,Ndim
               allocate (subtensor(1)%rows(dim_i)%dat(subtensor(1)%nr(dim_i)))
               subtensor(1)%rows(dim_i)%dat=inter_MD(nn)%rows(dim_i)%dat
               allocate (subtensor(1)%cols(dim_i)%dat(subtensor(1)%nc(dim_i)))
               subtensor(1)%cols(dim_i)%dat=inter_MD(nn)%cols(dim_i)%dat
            enddo
            allocate (subtensor(1)%dat(product(subtensor(1)%nr),product(subtensor(1)%nc)))
            subtensor(1)%dat=0
            call element_Zmn_tensorlist_user(Ndim, subtensor, 1, msh, option, ker, 0, passflag, ptree, stats)
            allocate (Mat(product(subtensor(1)%nr),product(subtensor(1)%nc)))
            Mat = subtensor(1)%dat

            do myi = 1, product(subtensor(1)%nr)
            do myj = 1, product(subtensor(1)%nc)
               value2 = inter_MD(nn)%dat(myi,myj)
               value1 = Mat(myi, myj)
               v1 = v1 + abs(value1)**2d0
               v2 = v2 + abs(value2)**2d0
               v3 = v3 + abs(value2 - value1)**2d0
               ! if(abs(value2-value1)**2d0/abs(value1)**2d0>1D-3)write(*,*)ptree%MyID,nn,myi,myj,abs(value2-value1)**2d0/abs(value1)**2d0,value2,value1
            enddo
            enddo
            dims_one=1
            call BF_MD_delete_subtensors(Ndim, dims_one, subtensor, stats)
            deallocate(Mat)
         else
            call element_Zmn_tensorlist_user(Ndim, subtensor, 0, msh, option, ker, 2, passflag, ptree, stats)
         endif
      enddo
      dims_1D=Ninter
      call BF_MD_delete_subtensors(Ndim, dims_1D, inter_MD, stats)

      n2 = MPI_Wtime()
      vtmp=v1
      call MPI_ALLREDUCE(vtmp, v1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(blocks%pgno)%Comm  , ierr)
      vtmp=v2
      call MPI_ALLREDUCE(vtmp, v2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(blocks%pgno)%Comm  , ierr)
      vtmp=v3
      call MPI_ALLREDUCE(vtmp, v3, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%pgrp(blocks%pgno)%Comm  , ierr)

      if (ptree%MyID - ptree%pgrp(blocks%pgno)%head == Main_ID .and. verbosity >= 0) write (*, '(A25,Es14.7,Es14.7,A6,Es9.2,A7,Es9.2)') 'BF_CheckError: fnorm:', sqrt(v1), sqrt(v2), ' acc: ', sqrt(v3/v1), ' time: ', n2 - n1
      if(present(error)) error=sqrt(v3/v1)

      !stop

   end subroutine BF_MD_checkError



   subroutine BF_compress_entry(bplus, option, Memory, stats, msh, ker, ptree)

      implicit none

      type(blockplus)::bplus
      integer:: ii, ll, bb, ierr, pp
      real(kind=8) Memory, rtemp, error
      integer:: level_butterfly, level_BP, levelm, groupm_start, Nboundall, Ninadmissible, statflag, knn_tmp
      type(Hoption)::option
      type(Hstat)::stats
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      type(matrixblock), pointer::blocks

      Memory = 0
      do ll = 1, bplus%Lplus
         bplus%LL(ll)%rankmax = 0
         statflag = 0
         if (ll == 1 .or. option%bp_cnt_lr == 1) statflag = 1  !!! only record the rank of the top-layer butterfly in a bplus
         do bb = 1, bplus%LL(ll)%Nbound
            if (IOwnPgrp(ptree, bplus%LL(ll)%matrices_block(bb)%pgno)) then
               if (bplus%LL(ll)%matrices_block(bb)%style == 1) then
                  call Full_construction(bplus%LL(ll)%matrices_block(bb), msh, ker, stats, option, ptree, rtemp)
                  Memory = Memory + rtemp
               else

                  level_butterfly = bplus%LL(ll)%matrices_block(bb)%level_butterfly
                  level_BP = bplus%level

                  bplus%LL(ll)%matrices_block(bb)%level_half = BF_Switchlevel(bplus%LL(ll)%matrices_block(bb)%level_butterfly, option%pat_comp)
                  levelm = bplus%LL(ll)%matrices_block(bb)%level_half
                  groupm_start = bplus%LL(ll)%matrices_block(1)%row_group*2**levelm
                  Nboundall = 0
                  Ninadmissible = 0
                  if (allocated(bplus%LL(ll + 1)%boundary_map)) then
                     Nboundall = size(bplus%LL(ll + 1)%boundary_map, 1)
                     Ninadmissible = size(bplus%LL(ll + 1)%boundary_map, 2)
                  endif
                  if (option%forwardN15flag == 1) then
                     knn_tmp = option%knn
                     option%knn=0
                     call BF_compress_N15(bplus%LL(ll)%matrices_block(bb), bplus%LL(ll + 1)%boundary_map, Nboundall, Ninadmissible, groupm_start, option, rtemp, stats, msh, ker, ptree, statflag)
                     option%knn=knn_tmp
                     call BF_sym2asym(bplus%LL(ll)%matrices_block(bb))

                     ! Move singular values to leftmost factor
                     if(bplus%LL(ll)%matrices_block(bb)%level_butterfly>0)then
                        call BF_ChangePattern(bplus%LL(ll)%matrices_block(bb), 3, 1, stats, ptree)
                        call BF_MoveSingular_Ker(bplus%LL(ll)%matrices_block(bb), 'N', 1, bplus%LL(ll)%matrices_block(bb)%level_butterfly, ptree, stats, option%tol_comp)
                        call BF_ChangePattern(bplus%LL(ll)%matrices_block(bb), 1, 3, stats, ptree)
                     endif
                  elseif (option%forwardN15flag == 2) then
                        call BF_compress_NlogN(bplus%LL(ll)%matrices_block(bb), bplus%LL(ll + 1)%boundary_map, Nboundall, Ninadmissible, groupm_start, option, rtemp, stats, msh, ker, ptree, statflag)
                        allocate(blocks)
                        call BF_copy('N', bplus%LL(ll)%matrices_block(bb), blocks)
                        call BF_checkError(blocks, option, msh, ker, stats, ptree, 1, -1, error)
                        deallocate(blocks)
                        if(error>5*option%tol_comp)then
                           pp = ptree%myid - ptree%pgrp(bplus%LL(ll)%matrices_block(bb)%pgno)%head + 1
                           if(option%verbosity>=0 .and. pp==1)write(*,*)'warning: BF',bplus%LL(ll)%matrices_block(bb)%row_group,bplus%LL(ll)%matrices_block(bb)%col_group,' error ',error,', redo the compression with the N15 algorithm'
                           call BF_delete(bplus%LL(ll)%matrices_block(bb),0)
                           knn_tmp = option%knn
                           option%knn=0
                           call BF_compress_N15(bplus%LL(ll)%matrices_block(bb), bplus%LL(ll + 1)%boundary_map, Nboundall, Ninadmissible, groupm_start, option, rtemp, stats, msh, ker, ptree, statflag)
                           option%knn=knn_tmp
                           call BF_sym2asym(bplus%LL(ll)%matrices_block(bb))
                           ! Move singular values to leftmost factor
                           if(bplus%LL(ll)%matrices_block(bb)%level_butterfly>0)then
                              call BF_ChangePattern(bplus%LL(ll)%matrices_block(bb), 3, 1, stats, ptree)
                              call BF_MoveSingular_Ker(bplus%LL(ll)%matrices_block(bb), 'N', 1, bplus%LL(ll)%matrices_block(bb)%level_butterfly, ptree, stats, option%tol_comp)
                              call BF_ChangePattern(bplus%LL(ll)%matrices_block(bb), 1, 3, stats, ptree)
                           endif
                        endif
                  else
                     call BF_compress_NlogN(bplus%LL(ll)%matrices_block(bb), bplus%LL(ll + 1)%boundary_map, Nboundall, Ninadmissible, groupm_start, option, rtemp, stats, msh, ker, ptree, statflag)
                  end if

                  Memory = Memory + rtemp
                  bplus%LL(ll)%rankmax = max(bplus%LL(ll)%rankmax, bplus%LL(ll)%matrices_block(bb)%rankmax)
               endif
            endif
         end do
      end do

      ! !!!!!!! check error
      if (option%ErrFillFull == 1) call Bplus_CheckError_Full(bplus, option, msh, ker, stats, ptree)
      ! !!!!!!! check error

      ! if(bplus%LL(1)%matrices_block(1)%level==1)write(*,*)bplus%LL(1)%rankmax,'dfdfdfdf'

      return

   end subroutine BF_compress_entry



!!!!!!! extract a list of intersections from a block
   subroutine BF_ExtractElement(blocks_o, option, msh, stats, ptree, Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps)

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
      DT,target::alldat_loc(:)
      integer::colidx(Ninter), rowidx(Ninter), pgidx(Ninter)
      integer::Npmap, pmaps(Npmap, 3), flag2D
      type(iarray)::lst
      type(matrixblock), pointer::blocks_o

      stats%Flop_Tmp = 0d0

      n0 = MPI_Wtime()
      flag2D = 0
      allocate (inters(Ninter))
      lstr = list()
      lstc = list()
      ! lst=list()
      lstblk = list()
      idx_row = 0
      idx_col = 0
      ntot_loc = 0
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
            nr_loc = myArows
            nc_loc = myAcols
            if (nr_loc > 0 .and. nc_loc > 0)then
               call Array1DtoPointer2D(alldat_loc(ntot_loc+1:ntot_loc + nr_loc*nc_loc), inters(nn)%dat_loc, nr_loc, nc_loc)
               ntot_loc = ntot_loc + nr_loc*nc_loc
            endif
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
      if(ntot_loc>0)alldat_loc(1:ntot_loc)=0

      n1 = MPI_Wtime()

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

      n2 = MPI_Wtime()
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
               call LogMemory(stats, SIZEOF(blocks%inters(nn)%dat_loc)/1024.0d3)
            enddo

            ! extract entries on an array of intersections for each block

            if (blocks%style == 1) then
               call Full_block_extraction(blocks, inters, ptree, msh, stats, option)
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
      n3 = MPI_Wtime()
      stats%Time_Entry_BF = stats%Time_Entry_BF + n3-n2

      call MPI_barrier(ptree%pgrp(blocks_o%pgno)%Comm, ierr)
      n3 = MPI_Wtime()


      ! redistribute from blocks' intersections to the global intersecions inters
      if (flag2D == 1) then ! if each intersection is only needed by one processor, the communication can be optimized
         call BPACK_all2all_inters(Ninter, inters, lstblk, stats, ptree, blocks_o%pgno, ptree%pgrp(blocks_o%pgno)%nproc, Npmap, pmaps)
      else
         call BPACK_all2all_inters_optimized(Ninter, inters, lstblk, stats, ptree, blocks_o%pgno, ptree%pgrp(blocks_o%pgno)%nproc, Npmap, pmaps)
      endif

      n4 = MPI_Wtime()
      stats%Time_Entry_Comm = stats%Time_Entry_Comm + n4-n3


      ! deallocate intersections at each block
      cur => lstblk%head
      do ii = 1, lstblk%num_nods
         select type (ptr=>cur%item)
         type is (block_ptr)
            blocks => ptr%ptr
            do nn = 1, size(blocks%inters, 1)
               if (associated(blocks%inters(nn)%dat)) then
                  call LogMemory(stats, -SIZEOF(blocks%inters(nn)%dat)/1024.0d3)
                  deallocate (blocks%inters(nn)%dat)
               endif
               if (associated(blocks%inters(nn)%dat_loc))then
                  call LogMemory(stats, -SIZEOF(blocks%inters(nn)%dat_loc)/1024.0d3)
                  deallocate (blocks%inters(nn)%dat_loc)
               endif
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
         if (associated(inters(nn)%dat))then
            call LogMemory(stats, -SIZEOF(inters(nn)%dat)/1024.0d3)
            deallocate (inters(nn)%dat)
         endif
         ! if (associated(inters(nn)%dat_loc)) deallocate (inters(nn)%dat_loc)
         if (allocated(inters(nn)%rows)) deallocate (inters(nn)%rows)
         if (allocated(inters(nn)%cols)) deallocate (inters(nn)%cols)
         if (allocated(inters(nn)%rows_loc)) deallocate (inters(nn)%rows_loc)
      enddo
      deallocate (inters)

      call list_finalizer(lstr)
      call list_finalizer(lstc)

      n5 = MPI_Wtime()

      ! if(ptree%MyID==Main_ID)then
      ! write(*,*)n1-n0,n2-n1,n3-n2,n4-n3,n5-n4
      ! endif

   end subroutine BF_ExtractElement


!!!!!!! extract a list of intersections from a blockplus
   subroutine BP_ExtractElement(BP_o, option, msh, stats, ptree, Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps)

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
      DT,target::alldat_loc(:)
      integer::colidx(Ninter), rowidx(Ninter), pgidx(Ninter)
      integer::Npmap, pmaps(Npmap, 3), flag2D
      type(iarray)::lst
      type(matrixblock), pointer::blocks_o
      type(blockplus)::BP_o
      blocks_o => BP_o%ll(1)%matrices_block(1)

      stats%Flop_Tmp = 0d0

      n0 = MPI_Wtime()
      flag2D = 0
      allocate (inters(Ninter))
      lstr = list()
      lstc = list()
      ! lst=list()
      lstblk = list()
      idx_row = 0
      idx_col = 0
      ntot_loc = 0
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
            nr_loc = myArows
            nc_loc = myAcols
            if (nr_loc > 0 .and. nc_loc > 0)then
               call Array1DtoPointer2D(alldat_loc(ntot_loc+1:ntot_loc + nr_loc*nc_loc), inters(nn)%dat_loc, nr_loc, nc_loc)
               ntot_loc = ntot_loc + nr_loc*nc_loc
            endif
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
      if(ntot_loc>0)alldat_loc(1:ntot_loc)=0

      n1 = MPI_Wtime()

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
                  call Hmat_MapIntersec2Block_Loc(blocks_o, option, stats, msh, ptree, inters, nn, ptrr, ptrc, lstblk)
               case (HMAT,BLR)
                  call Hmat_MapIntersec2Block_Loc(blocks_o, option, stats, msh, ptree, inters, nn, ptrr, ptrc, lstblk)
               case (HSS)
                  call BP_MapIntersec2Block(BP_o, option, stats, msh, ptree, inters, nn, ptrr, ptrc, lstblk, 1, BP_o%LL(1)%Nbound)
               end select

            end select
         end select
         curr => curr%next
         curc => curc%next
      enddo

      call MergeSort(lstblk%head, node_score_block_ptr_row)

      n2 = MPI_Wtime()
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
               call LogMemory(stats, SIZEOF(blocks%inters(nn)%dat_loc)/1024.0d3)
            enddo

            ! extract entries on an array of intersections for each block

            if (blocks%style == 1) then
               call Full_block_extraction(blocks, inters, ptree, msh, stats, option)
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
      n3 = MPI_Wtime()
      stats%Time_Entry_BF = stats%Time_Entry_BF + n3-n2

      call MPI_barrier(ptree%pgrp(BP_o%pgno)%Comm, ierr)
      n3 = MPI_Wtime()


      ! redistribute from blocks' intersections to the global intersecions inters
      if (flag2D == 1) then ! if each intersection is only needed by one processor, the communication can be optimized
         call BPACK_all2all_inters(Ninter, inters, lstblk, stats, ptree, BP_o%pgno, ptree%pgrp(BP_o%pgno)%nproc, Npmap, pmaps)
      else
         call BPACK_all2all_inters_optimized(Ninter, inters, lstblk, stats, ptree, BP_o%pgno, ptree%pgrp(BP_o%pgno)%nproc, Npmap, pmaps)
      endif

      n4 = MPI_Wtime()
      stats%Time_Entry_Comm = stats%Time_Entry_Comm + n4-n3


      ! deallocate intersections at each block
      cur => lstblk%head
      do ii = 1, lstblk%num_nods
         select type (ptr=>cur%item)
         type is (block_ptr)
            blocks => ptr%ptr
            do nn = 1, size(blocks%inters, 1)
               if (associated(blocks%inters(nn)%dat)) then
                  call LogMemory(stats, -SIZEOF(blocks%inters(nn)%dat)/1024.0d3)
                  deallocate (blocks%inters(nn)%dat)
               endif
               if (associated(blocks%inters(nn)%dat_loc))then
                  call LogMemory(stats, -SIZEOF(blocks%inters(nn)%dat_loc)/1024.0d3)
                  deallocate (blocks%inters(nn)%dat_loc)
               endif
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
         if (associated(inters(nn)%dat))then
            call LogMemory(stats, -SIZEOF(inters(nn)%dat)/1024.0d3)
            deallocate (inters(nn)%dat)
         endif
         ! if (associated(inters(nn)%dat_loc)) deallocate (inters(nn)%dat_loc)
         if (allocated(inters(nn)%rows)) deallocate (inters(nn)%rows)
         if (allocated(inters(nn)%cols)) deallocate (inters(nn)%cols)
         if (allocated(inters(nn)%rows_loc)) deallocate (inters(nn)%rows_loc)
      enddo
      deallocate (inters)

      call list_finalizer(lstr)
      call list_finalizer(lstc)

      n5 = MPI_Wtime()

      ! if(ptree%MyID==Main_ID)then
      ! write(*,*)n1-n0,n2-n1,n3-n2,n4-n3,n5-n4
      ! endif

   end subroutine BP_ExtractElement


!!!!!!! extract a list of intersections from a bmat
   subroutine BPACK_ExtractElement(bmat, option, msh, stats, ptree, Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps)

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
      DT,target::alldat_loc(:)
      integer::colidx(Ninter), rowidx(Ninter), pgidx(Ninter)
      integer::Npmap, pmaps(Npmap, 3), flag2D
      type(iarray)::lst

      stats%Flop_Tmp = 0d0

      n0 = MPI_Wtime()
      flag2D = 0
      allocate (inters(Ninter))
      lstr = list()
      lstc = list()
      ! lst=list()
      lstblk = list()
      idx_row = 0
      idx_col = 0
      ntot_loc = 0
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
            nr_loc = myArows
            nc_loc = myAcols
            if (nr_loc > 0 .and. nc_loc > 0)then
               call Array1DtoPointer2D(alldat_loc(ntot_loc+1:ntot_loc + nr_loc*nc_loc), inters(nn)%dat_loc, nr_loc, nc_loc)
               ntot_loc = ntot_loc + nr_loc*nc_loc
            endif
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

      if(ntot_loc>0)alldat_loc(1:ntot_loc)=0

      n1 = MPI_Wtime()

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
               case (HMAT,BLR)
                  num_blocks = 2**msh%Dist_level
                  call Hmat_MapIntersec2Block(bmat%h_mat, option, stats, msh, ptree, inters, nn, ptrr, ptrc, lstblk, num_blocks)
               case (HSS)
                  call BP_MapIntersec2Block(bmat%hss_bf%BP, option, stats, msh, ptree, inters, nn, ptrr, ptrc, lstblk, 1, bmat%hss_bf%BP%LL(1)%Nbound)
               end select
            end select
         end select
         curr => curr%next
         curc => curc%next
      enddo

      call MergeSort(lstblk%head, node_score_block_ptr_row)
      ! (a HODLR built on the GPU whose host copies of the factors are unfilled, HODLR_gpu_defer_host: the rows
      ! read below, from the GPUs, into each block's gpu_u and gpu_v)
      if (option%format == HODLR) then
         if (bmat%ho_bf%gpu_host_stale) call HODLR_gpu_fill_rows(bmat%ho_bf, lstblk, inters, option, ptree)
      endif

      n2 = MPI_Wtime()
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
               call LogMemory(stats, SIZEOF(blocks%inters(nn)%dat_loc)/1024.0d3)
            enddo

            ! extract entries on an array of intersections for each block

            if (blocks%style == 1) then
               call Full_block_extraction(blocks, inters, ptree, msh, stats, option)
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
      n3 = MPI_Wtime()
      stats%Time_Entry_BF = stats%Time_Entry_BF + n3-n2

      call MPI_barrier(ptree%Comm, ierr)
      n3 = MPI_Wtime()
      ! redistribute from blocks' intersections to the global intersecions inters
      if (flag2D == 1) then ! if each intersection is only needed by one processor, the communication can be optimized
         call BPACK_all2all_inters(Ninter, inters, lstblk, stats, ptree, 1, ptree%nproc, Npmap, pmaps)
      else
         call BPACK_all2all_inters_optimized(Ninter, inters, lstblk, stats, ptree, 1, ptree%nproc, Npmap, pmaps)
      endif

      n4 = MPI_Wtime()
      stats%Time_Entry_Comm = stats%Time_Entry_Comm + n4-n3

      ! deallocate intersections at each block
      cur => lstblk%head
      do ii = 1, lstblk%num_nods
         select type (ptr=>cur%item)
         type is (block_ptr)
            blocks => ptr%ptr
            do nn = 1, size(blocks%inters, 1)
               if (associated(blocks%inters(nn)%dat)) then
                  call LogMemory(stats, -SIZEOF(blocks%inters(nn)%dat)/1024.0d3)
                  deallocate (blocks%inters(nn)%dat)
               endif
               if (associated(blocks%inters(nn)%dat_loc))then
                  call LogMemory(stats, -SIZEOF(blocks%inters(nn)%dat_loc)/1024.0d3)
                  deallocate (blocks%inters(nn)%dat_loc)
               endif
               if (allocated(blocks%inters(nn)%rows)) deallocate (blocks%inters(nn)%rows)
               if (allocated(blocks%inters(nn)%cols)) deallocate (blocks%inters(nn)%cols)
               if (allocated(blocks%inters(nn)%rows_loc)) deallocate (blocks%inters(nn)%rows_loc)
               blocks%inters(nn)%nr = 0
               blocks%inters(nn)%nr_loc = 0
               blocks%inters(nn)%nc = 0
               blocks%inters(nn)%idx = 0
            enddo
            deallocate (blocks%inters)
            ! (the rows taken from the GPUs, HODLR_gpu_fill_rows)
            if (allocated(blocks%gpu_urow)) deallocate (blocks%gpu_urow)
            if (allocated(blocks%gpu_vrow)) deallocate (blocks%gpu_vrow)
            if (allocated(blocks%gpu_u)) deallocate (blocks%gpu_u)
            if (allocated(blocks%gpu_v)) deallocate (blocks%gpu_v)
         end select
         cur => cur%next
      enddo

      ! finalize the list of block_ptr
      call list_finalizer(lstblk)

      ! deallocate global intersections
      do nn = 1, Ninter
         if (associated(inters(nn)%dat))then
            call LogMemory(stats, -SIZEOF(inters(nn)%dat)/1024.0d3)
            deallocate (inters(nn)%dat)
         endif
         ! if (associated(inters(nn)%dat_loc)) deallocate (inters(nn)%dat_loc)
         if (allocated(inters(nn)%rows)) deallocate (inters(nn)%rows)
         if (allocated(inters(nn)%cols)) deallocate (inters(nn)%cols)
         if (allocated(inters(nn)%rows_loc)) deallocate (inters(nn)%rows_loc)
      enddo
      deallocate (inters)

      n5 = MPI_Wtime()

      call list_finalizer(lstr)
      call list_finalizer(lstc)
      ! time_tmp = time_tmp + n2- n1
      ! if(ptree%MyID==Main_ID)then
      ! write(*,*)n1-n0,n2-n1,n3-n2,n4-n3,n5-n4
      ! endif

   end subroutine BPACK_ExtractElement

!!!!!!! check error of BPACK construction using parallel element extraction
   subroutine BPACK_CheckError_entry(bmat, option, msh, ker, stats, ptree)

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
      real(kind=8):: a, v1, v2, v3, vtmp
      DT:: value1, value2, value3
      type(list)::lstr, lstc, lst, lstblk
      type(nod), pointer::cur, curr, curc, curri, curci
      class(*), pointer::ptr, ptrr, ptrc, ptrri, ptrci
      integer::head, tail, idx, pgno, ctxt, nr_loc, nc_loc
      type(matrixblock), pointer::blocks
      integer num_blocks, idx_row, idx_col, idx_dat
      integer, allocatable:: allrows(:), allcols(:), pmaps(:, :)
      integer, allocatable::datidx(:), colidx(:), rowidx(:), pgidx(:)
      DT, target, allocatable::alldat_loc(:)
      DT, allocatable::dat_gpu(:)
      integer:: Ninter, nr, nrmax, nc, ncmax, ntot_loc, level, Npmap, nproc, npavr, np
      type(intersect)::submats(1)

      ! select case(option%format)
      ! case(HODLR)
      ! level=bmat%ho_bf%Maxlevel
      ! case(HMAT,BLR)
      ! level=bmat%h_mat%Maxlevel
      ! end select

      ! level=1

      ! Ninter=2**level
      ! nr=2500
      ! nc=2500

      Ninter = 4
      ! nr=msh%Nunk
      ! nc=msh%Nunk

      nrmax = 100
      ncmax = 100

      allocate (colidx(Ninter))
      allocate (rowidx(Ninter))
      allocate (pgidx(Ninter))
      ! allocate(datidx(Ninter))

      allocate (allrows(Ninter*nrmax))
      allocate (allcols(Ninter*ncmax))

      ! pgno=1
      ! ctxt = ptree%pgrp(pgno)%ctxt
      ! call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
      ! nprow = ptree%pgrp(pgno)%nprow
      ! npcol = ptree%pgrp(pgno)%npcol

      nproc = ptree%nproc
      Npmap = min(Ninter, nproc)
      npavr = nproc/Npmap
      allocate (pmaps(Npmap, 3))
      do nn = 1, Npmap
         nprow = floor_safe(sqrt(dble(npavr)))
         npcol = floor_safe(npavr/dble(nprow))
         pmaps(nn, 1) = 1   ! nprow   ! this makes sure the intersection is on 1 processor, this makes it easier for cpp user-defined extraction function
         pmaps(nn, 2) = 1   ! npcol
         pmaps(nn, 3) = (nn - 1)*npavr
      enddo

      idx_row = 0
      idx_col = 0
      idx_dat = 0
      ! ntot_loc=0
      pp = 0
      do nn = 1, Ninter

         do ii = 1, nrmax
            call random_number(a)
            call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, Main_ID, ptree%Comm, ierr)
            allrows(idx_row + 1) = max(floor_safe(msh%Nunk*a), 1)
            ! allrows(idx_row + 1) = max(floor_safe(3125*a), 1)+3125*0
            ! allrows(idx_row + 1) = max(floor_safe(7812*a), 1)+7812*0
            ! allrows(idx_row + 1) = max(floor_safe(19531*a), 1)+19531*0
            ! allrows(idx_row+1)=msh%basis_group(2**level+nn-1)%head+ii-1
            idx_row = idx_row + 1
         enddo
         call remove_dup_int(allrows(idx_row-nrmax+1:idx_row), nrmax, nr)
         idx_row = idx_row - (nrmax-nr)


         do ii = 1, ncmax
            call random_number(a)
            call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, Main_ID, ptree%Comm, ierr)
            allcols(idx_col + 1) = max(floor_safe(msh%Nunk*a), 1)
            ! allcols(idx_col + 1) = max(floor_safe(3125*a), 1)+3125*1
            ! allcols(idx_col + 1) = max(floor_safe(7812*a), 1)+7812*1
            ! allcols(idx_col + 1) = max(floor_safe(19531*a), 1)+19531*1
            ! allcols(idx_col+1)=msh%basis_group(2**level+1-(nn-1))%head+ii-1
            idx_col = idx_col + 1
         enddo
         call remove_dup_int(allcols(idx_col-ncmax+1:idx_col), ncmax, nc)
         idx_col = idx_col - (ncmax-nc)


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


      enddo

      allocate (alldat_loc(idx_dat))
      if (idx_dat > 0) alldat_loc = 0

      n1 = MPI_Wtime()
      call BPACK_ExtractElement(bmat, option, msh, stats, ptree, Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps)
      n2 = MPI_Wtime()
      ! (BPACK_CHECK=hodlr: the same entries with the host factors overwritten by a sentinel and marked unfilled, so
      ! that the extraction takes the rows it reads from the GPUs; then the host factors restored from the GPUs)
      if (option%format == HODLR .and. option%HODLR_use_gpu > 0 .and. HODLR_gpu_check_level() > 0) then
         if (c_associated(bmat%ho_bf%gpu) .and. .not. bmat%ho_bf%gpu_host_stale) then
            allocate (dat_gpu(max(1, idx_dat)))
            dat_gpu = 0
            call HODLR_gpu_poison_host(bmat%ho_bf, option, ptree)
            call BPACK_ExtractElement(bmat, option, msh, stats, ptree, Ninter, allrows, allcols, dat_gpu, rowidx, &
                                      colidx, pgidx, Npmap, pmaps)
            call HODLR_gpu_fetch_host(bmat%ho_bf, option, ptree)
            vtmp = 0
            if (idx_dat > 0) vtmp = maxval(abs(dat_gpu(1:idx_dat) - alldat_loc(1:idx_dat)))
            call MPI_ALLREDUCE(vtmp, v3, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
            vtmp = 0
            if (idx_dat > 0) vtmp = maxval(abs(alldat_loc(1:idx_dat)))
            call MPI_ALLREDUCE(vtmp, v1, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
            if (ptree%MyID == Main_ID) write (*, '(A,Es11.3,A,Es11.3)') ' HODLR GPU check (entries): largest difference'// &
               ' of the entries from the rows on the GPUs to those from the host factors', v3, ', largest entry', v1
            deallocate (dat_gpu)
         endif
      endif

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
      vtmp=v1
      call MPI_ALLREDUCE(vtmp, v1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      vtmp=v2
      call MPI_ALLREDUCE(vtmp, v2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      vtmp=v3
      call MPI_ALLREDUCE(vtmp, v3, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A32,Es14.7,Es14.7,A6,Es9.2,A7,Es9.2)') 'BPACK_CheckError(entry): fnorm:', sqrt(v1), sqrt(v2), ' acc: ', sqrt(v3/v1), ' time: ', n2 - n1

   end subroutine BPACK_CheckError_entry




!!!!!!! check error of BPACK construction using BPACK matvec with a sparse vector
   subroutine BPACK_CheckError_SMVP(bmat, option, msh, ker, stats, ptree)

      implicit none

      type(Hoption)::option
      type(Hstat)::stats
      type(Bmatrix)::bmat
      type(mesh)::msh
      type(kernelquant)::ker
      type(proctree)::ptree
      integer nvec
      DT, allocatable:: x_loc(:, :), rhs_loc(:,:),rhs_loc_ref(:,:)
      DT::tmp
      integer dim_i,ij,ii,ii1,ij1
      integer:: Nunk_n_loc,idxs,idxe,idx_1,idx_2
      integer,allocatable:: idx_src(:)
      integer:: Npt_src, N_glo, Npt_src_tmp
      real(kind=8):: a, v1, v2, v3, vtmp
      integer ierr
      type(intersect) :: submats(1)
      integer passflag
      real(kind=8)::n1, n2, n3, n4
      integer:: dims_one

      select case (option%format)
      case (HODLR)
         N_glo=bmat%ho_bf%N
      case (HMAT,BLR)
         N_glo=bmat%h_mat%N
      case (HSS)
         N_glo=bmat%hss_bf%N
      case default
         write(*,*)'not supported format in BPACK_CheckError:', option%format
         stop
      end select
      Nunk_n_loc = msh%idxe - msh%idxs + 1
      idxs= msh%idxs

      if (allocated(bmat%xtrue)) then
         call LogMemory(stats, -SIZEOF(bmat%xtrue)/1024.0d3)
         deallocate(bmat%xtrue)
      endif
      if (allocated(bmat%b_true)) then
         call LogMemory(stats, -SIZEOF(bmat%b_true)/1024.0d3)
         deallocate(bmat%b_true)
      endif

      nvec=1 !! currently this can only be 1
      allocate(x_loc(Nunk_n_loc,nvec))
      x_loc=0
      call LogMemory(stats, SIZEOF(x_loc)/1024.0d3)

      Npt_src = min(20,N_glo)
      allocate(idx_src(Npt_src))
      do ij=1,Npt_src
         call random_number(a)
         idx_src(ij) = max(floor_safe(N_glo*a), 1)
      enddo
      call MPI_Bcast(idx_src, Npt_src, MPI_INTEGER, Main_ID, ptree%Comm, ierr)

      Npt_src_tmp = Npt_src
      call remove_dup_int(idx_src, Npt_src_tmp, Npt_src)

      ! Npt_src=1
      ! allocate(idx_src(Npt_src))
      ! idx_src(1) =1

      idxe = idxs + Nunk_n_loc -1
      do ij=1,Npt_src
         idx_1=idx_src(ij)
         if(idx_1>=idxs .and. idx_1<=idxe)then
            idx_1 = idx_1 - idxs + 1
            x_loc(idx_1,1) = x_loc(idx_1,1) + 1
         endif
      enddo

      n1 = MPI_Wtime()

      !! Generate rhs_loc by using BPACK_MD_Mult
      allocate(rhs_loc(Nunk_n_loc,nvec))
      call LogMemory(stats, SIZEOF(rhs_loc)/1024.0d3)
      rhs_loc=0

      call BPACK_Mult('N', Nunk_n_loc, nvec, x_loc, rhs_loc, bmat, ptree, option, stats,0)


      !! Generate the reference rhs_loc_ref by using element_Zmn_tensorlist_user
      allocate(rhs_loc_ref(Nunk_n_loc,nvec))
      call LogMemory(stats, SIZEOF(rhs_loc_ref)/1024.0d3)
      rhs_loc_ref=0


      submats(1)%nr = Nunk_n_loc
      submats(1)%nc = Npt_src
      allocate(submats(1)%rows(submats(1)%nr))
      allocate(submats(1)%cols(submats(1)%nc))
      allocate(submats(1)%dat(submats(1)%nr,submats(1)%nc))
      call LogMemory(stats, SIZEOF(submats(1)%dat)/1024.0d3)
      submats(1)%dat = 0
      do ii=1,submats(1)%nr
         submats(1)%rows(ii) = msh%idxs + ii - 1
      enddo
      do ij=1,Npt_src
         idx_1=idx_src(ij)
         submats(1)%cols(ij) = idx_1
      enddo

      call element_Zmn_blocklist_user(submats, 1, msh, option, ker, 0, passflag, ptree, stats)
      do ij=1,Npt_src
         rhs_loc_ref(:,1) = rhs_loc_ref(:,1) + submats(1)%dat(:,ij)
      enddo

      deallocate(submats(1)%rows)
      deallocate(submats(1)%cols)
      call LogMemory(stats, -SIZEOF(submats(1)%dat)/1024.0d3)
      deallocate(submats(1)%dat)

      n2 = MPI_Wtime()

      v1 =(fnorm(rhs_loc,Nunk_n_loc,nvec))**2d0
      v2 =(fnorm(rhs_loc_ref,Nunk_n_loc,nvec))**2d0
      rhs_loc = rhs_loc - rhs_loc_ref
      v3 =(fnorm(rhs_loc,Nunk_n_loc,nvec))**2d0

      vtmp = v1
      call MPI_ALLREDUCE(vtmp, v1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      vtmp = v2
      call MPI_ALLREDUCE(vtmp, v2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      vtmp = v3
      call MPI_ALLREDUCE(vtmp, v3, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A30,Es14.7,Es14.7,A6,Es9.2,A7,Es9.2)') 'BPACK_CheckError(mvp): fnorm:', sqrt(v1), sqrt(v2), ' acc: ', sqrt(v3/v1), ' time: ', n2 - n1

      call LogMemory(stats, -SIZEOF(rhs_loc)/1024.0d3)
      deallocate(rhs_loc)
      call move_alloc(x_loc, bmat%xtrue)
      call move_alloc(rhs_loc_ref, bmat%b_true)
      deallocate(idx_src)

   end subroutine BPACK_CheckError_SMVP




!!!!!!! check error of BPACK_MD construction using BPACK_MD matvec with a sparse vector
   subroutine BPACK_MD_CheckError_SMVP(Ndim, bmat, option, msh, ker, stats, ptree)

      implicit none

      integer Ndim
      type(Hoption)::option
      type(Hstat)::stats
      type(Bmatrix)::bmat
      type(mesh)::msh(Ndim)
      type(kernelquant)::ker
      type(proctree)::ptree
      integer nvec
      DT, allocatable:: x_loc(:, :), rhs_loc(:,:),rhs_loc_ref(:,:)
      DT::tmp
      integer dim_i,ij,ii,ii1,ij1
      integer:: Nunk_n_loc(Ndim),idxs(Ndim),idxe(Ndim),idx_1(Ndim),idx_2(Ndim)
      integer,allocatable:: idx_src(:)
      integer:: Npt_src, N_glo(Ndim)
      real(kind=8):: a, v1, v2, v3, vtmp
      integer ierr
      type(intersect_MD) :: subtensors(1)
      integer passflag
      real(kind=8)::n1, n2, n3, n4
      integer:: dims_one(Ndim)

      select case (option%format)
      case (HSS_MD)
         N_glo=bmat%hss_bf_md%N
         idxs = bmat%hss_bf_md%BP%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1, 1, :)
      case (HTENSOR)
         N_glo=bmat%h_mat_md%N
         idxs = bmat%h_mat_md%BP%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1, 1, :)
      case default
         write(*,*)'not supported format in BPACK_MD_CheckError:', option%format
         stop
      end select
      do dim_i=1,Ndim
      Nunk_n_loc(dim_i) = msh(dim_i)%idxe - msh(dim_i)%idxs + 1
      enddo
      nvec=1 !! currently this can only be 1
      allocate(x_loc(product(Nunk_n_loc),nvec))
      call LogMemory(stats, SIZEOF(x_loc)/1024.0d3)
      x_loc=0

      Npt_src = min(40,product(N_glo))
      allocate(idx_src(Npt_src))
      do ij=1,Npt_src
         call random_number(a)
         idx_src(ij) = max(floor_safe(product(N_glo)*a), 1)
      enddo
      call MPI_Bcast(idx_src, Npt_src, MPI_INTEGER, Main_ID, ptree%Comm, ierr)


      ! Npt_src=1
      ! allocate(idx_src(Npt_src))
      ! idx_src(1) =1

      idxe = idxs + Nunk_n_loc -1
      do ij=1,Npt_src
         call SingleIndexToMultiIndex(Ndim,N_glo, idx_src(ij), idx_1)
         if(ALL(idx_1>=idxs) .and. ALL(idx_1<=idxe))then
            idx_1 = idx_1 - idxs + 1
            call MultiIndexToSingleIndex(Ndim,Nunk_n_loc, ij1, idx_1)
            x_loc(ij1,1) = x_loc(ij1,1) + 1
         endif
      enddo

      n1 = MPI_Wtime()

      !! Generate rhs_loc by using BPACK_MD_Mult
      allocate(rhs_loc(product(Nunk_n_loc),nvec))
      call LogMemory(stats, SIZEOF(rhs_loc)/1024.0d3)
      rhs_loc=0
      call BPACK_MD_Mult(Ndim, 'N', Nunk_n_loc, nvec, x_loc, rhs_loc, bmat, ptree, option, stats, msh)


      !! Generate the reference rhs_loc_ref by using element_Zmn_tensorlist_user
      allocate(rhs_loc_ref(product(Nunk_n_loc),nvec))
      call LogMemory(stats, SIZEOF(rhs_loc_ref)/1024.0d3)
      rhs_loc_ref=0
      allocate(subtensors(1)%nr(Ndim))
      allocate(subtensors(1)%nc(Ndim))
      allocate(subtensors(1)%rows(Ndim))
      allocate(subtensors(1)%cols(Ndim))
      subtensors(1)%nr = Nunk_n_loc
      subtensors(1)%nc = 1
      do dim_i=1,Ndim
         allocate (subtensors(1)%rows(dim_i)%dat(subtensors(1)%nr(dim_i)))
         allocate (subtensors(1)%cols(dim_i)%dat(subtensors(1)%nc(dim_i)))
      enddo
      allocate(subtensors(1)%dat(product(subtensors(1)%nr),product(subtensors(1)%nc)))
      call LogMemory(stats, SIZEOF(subtensors(1)%dat)/1024.0d3)

      do ij=1,Npt_src
         subtensors(1)%dat = 0
         call SingleIndexToMultiIndex(Ndim,N_glo, idx_src(ij), idx_1)
         do dim_i=1,Ndim
            do ii=1,subtensors(1)%nr(dim_i)
               subtensors(1)%rows(dim_i)%dat(ii) = msh(dim_i)%idxs + ii - 1
            enddo
            subtensors(1)%cols(dim_i)%dat(1) = idx_1(dim_i)
         enddo

         call element_Zmn_tensorlist_user(Ndim, subtensors, 1, msh, option, ker, 0, passflag, ptree, stats)
         rhs_loc_ref = rhs_loc_ref + subtensors(1)%dat
      enddo
      dims_one=1
      call BF_MD_delete_subtensors(Ndim, dims_one, subtensors, stats)
      n2 = MPI_Wtime()

      v1 =(fnorm(rhs_loc,product(Nunk_n_loc),nvec))**2d0
      v2 =(fnorm(rhs_loc_ref,product(Nunk_n_loc),nvec))**2d0
      rhs_loc = rhs_loc - rhs_loc_ref
      v3 =(fnorm(rhs_loc,product(Nunk_n_loc),nvec))**2d0

      vtmp=v1
      call MPI_ALLREDUCE(vtmp, v1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      vtmp=v2
      call MPI_ALLREDUCE(vtmp, v2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      vtmp=v3
      call MPI_ALLREDUCE(vtmp, v3, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, '(A28,Es14.7,Es14.7,A6,Es9.2,A7,Es9.2)') 'BPACK_MD_CheckError: fnorm:', sqrt(v1), sqrt(v2), ' acc: ', sqrt(v3/v1), ' time: ', n2 - n1

      call LogMemory(stats, -SIZEOF(x_loc)/1024.0d3)
      deallocate(x_loc)
      call LogMemory(stats, -SIZEOF(rhs_loc)/1024.0d3)
      deallocate(rhs_loc)
      call LogMemory(stats, -SIZEOF(rhs_loc_ref)/1024.0d3)
      deallocate(rhs_loc_ref)
      deallocate(idx_src)

   end subroutine BPACK_MD_CheckError_SMVP







!>*********** all to all communication of element extraction results from local layout to 2D block-cyclic layout of each intersection (each process knows where to send, but doesn't know where to receive without communication)
   subroutine BPACK_all2all_inters(Ninter, inters, lstblk, stats, ptree, pgno, nproc, Npmap, pmaps)


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
      integer Ninter
      type(intersect)::inters(Ninter)
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
      n1 = MPI_Wtime()
      ! pgno = 1
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
               idstart = pmaps(inters(idx)%pg, 3) - ptree%pgrp(pgno)%head
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
         call LogMemory(stats, SIZEOF(sendquant(pp)%dat)/1024.0d3)
      enddo

      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         sendquant(pp)%size = 0
      enddo

      n2 = MPI_Wtime()

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
               idstart = pmaps(inters(idx)%pg, 3) - ptree%pgrp(pgno)%head
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
#ifdef HAVE_TASKLOOP
                  !$omp atomic capture
#endif
                  idxs = sendquant(pp)%size
                  sendquant(pp)%size = sendquant(pp)%size + 4
#ifdef HAVE_TASKLOOP
                  !$omp end atomic
#endif
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


      n3 = MPI_Wtime()

      Nreqs = 0
      do tt = 1, Nsendactive
         ! n6 = MPI_Wtime()
         pp = sendIDactive(tt)
         recvid = pp - 1 + ptree%pgrp(pgno)%head
         if (recvid /= ptree%MyID) then
            Nreqs = Nreqs + 1
            call MPI_Isend(sendquant(pp)%dat, sendquant(pp)%size, MPI_DT, pp - 1, tag + 1, ptree%pgrp(pgno)%Comm, S_req(Nreqs), ierr)
         endif
         ! n7 = MPI_Wtime()
         ! write(*,*)ptree%MyID,'to',pp-1,n7-n6,sendquant(pp)%size
      enddo

      n4 = MPI_Wtime()

      ! copy data from buffer to target
      cnt = 0
      do tt = 1, Nrecvactive
         if(tt==1 .and. Nreqr+1== Nrecvactive)then
            pp = ptree%MyID + 1 - ptree%pgrp(pgno)%head
            recvquant(pp)%size=sendquant(pp)%size
            if(recvquant(pp)%size>0)then
            allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
            call LogMemory(stats, SIZEOF(recvquant(pp)%dat)/1024.0d3)
            recvquant(pp)%dat = sendquant(pp)%dat
            endif
         else
            call MPI_Probe(MPI_ANY_SOURCE, tag+1, ptree%pgrp(pgno)%Comm, statusr(:,1),ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
            call MPI_Get_count(statusr(:,1), MPI_DT, recvquant(pp)%size,ierr)
            allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
            call LogMemory(stats, SIZEOF(recvquant(pp)%dat)/1024.0d3)
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



      n5 = MPI_Wtime()
      if (Nreqs > 0) then
         call MPI_waitall(Nreqs, S_req, statuss, ierr)
      endif
      ! deallocation
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         if (allocated(sendquant(pp)%dat))then
            call LogMemory(stats, -SIZEOF(sendquant(pp)%dat)/1024.0d3)
            deallocate (sendquant(pp)%dat)
         endif
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat))then
            call LogMemory(stats, -SIZEOF(recvquant(pp)%dat)/1024.0d3)
            deallocate (recvquant(pp)%dat)
         endif
      enddo

      ! n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1
#ifdef HAVE_TASKLOOP
      !$omp end single
      !$omp end parallel
#endif
      ! write(*,*)n2-n1,n3-n2,n4-n3,n5-n4,'wordi',ptree%MyID

   end subroutine BPACK_all2all_inters

!>*********** all to all communication of element extraction results from local layout to each entire intersection (each process knows where to send, but doesn't know where to receive without communication)
! YL: This subroutine seems to be slower than BPACK_all2all_inters
   subroutine BPACK_all2all_inters_optimized(Ninter, inters, lstblk, stats, ptree, pgno, nproc, Npmap, pmaps)

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
      integer Ninter
      type(intersect)::inters(Ninter)
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
#ifdef HAVE_OPENMP
      num_threads = omp_get_num_threads()
#else
      num_threads = 1
#endif

      nr_max = 0
      nc_max = 0
      do nn = 1, size(inters, 1)
         nr_max = max(nr_max, inters(nn)%nr)
         nc_max = max(nc_max, inters(nn)%nc)
      enddo
      allocate (ridx(nr_max,num_threads))
      allocate (cidx(nc_max,num_threads))

      n1 = MPI_Wtime()
      ! pgno = 1
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
               pp = pmaps(inters(idx)%pg, 3) + 1 - ptree%pgrp(pgno)%head
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
         call LogMemory(stats, SIZEOF(sendquant(pp)%dat)/1024.0d3)
      enddo

      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         sendquant(pp)%size = 0
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
      enddo

      n2 = MPI_Wtime()

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
               pp = pmaps(inters(idx)%pg, 3) + 1 - ptree%pgrp(pgno)%head
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

      n3 = MPI_Wtime()

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
            call LogMemory(stats, SIZEOF(recvquant(pp)%dat)/1024.0d3)
            recvquant(pp)%dat = sendquant(pp)%dat
            endif
         else
            call MPI_Probe(MPI_ANY_SOURCE, tag+1, ptree%pgrp(pgno)%Comm, statusr(:,1),ierr)
            pp = statusr(MPI_SOURCE, 1) + 1
            call MPI_Get_count(statusr(:,1), MPI_DT, recvquant(pp)%size,ierr)
            allocate (recvquant(pp)%dat(recvquant(pp)%size, 1))
            call LogMemory(stats, SIZEOF(recvquant(pp)%dat)/1024.0d3)
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
#ifdef HAVE_OPENMP
         my_tid = omp_get_thread_num()
#else
         my_tid = 0
#endif
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

      n5 = MPI_Wtime()
      if (Nreqs > 0) then
         call MPI_waitall(Nreqs, S_req, statuss, ierr)
      endif
      ! deallocation
      do tt = 1, Nsendactive
         pp = sendIDactive(tt)
         if (allocated(sendquant(pp)%dat))then
            call LogMemory(stats, -SIZEOF(sendquant(pp)%dat)/1024.0d3)
            deallocate (sendquant(pp)%dat)
         endif
      enddo
      do tt = 1, Nrecvactive
         pp = recvIDactive(tt)
         if (allocated(recvquant(pp)%dat))then
            call LogMemory(stats, -SIZEOF(recvquant(pp)%dat)/1024.0d3)
            deallocate (recvquant(pp)%dat)
         endif
      enddo

      deallocate (ridx)
      deallocate (cidx)

#ifdef HAVE_TASKLOOP
      !$omp end single
      !$omp end parallel
#endif
      ! n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1

      ! write(*,*)n2-n1,n3-n2,n4-n3,n5-n4,'wori',ptree%MyID

   end subroutine BPACK_all2all_inters_optimized

   recursive subroutine HODLR_MapIntersec2Block(ho_bf1, option, stats, msh, ptree, inters, nth, lstr, lstc, lstblk, level_c, bidx, flag)

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
            ! n0 = MPI_Wtime()

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
            ! n1 = MPI_Wtime()
            ! time_tmp = time_tmp + n1 - n0

            clstc(1)%idx = nth
            clstc(2)%idx = nth
            allocate (clstc(1)%dat(lstc%num_nods))
            allocate (clstc(2)%dat(lstc%num_nods))
            clstc(1)%num_nods = 0
            clstc(2)%num_nods = 0
            ! n0 = MPI_Wtime()

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

            ! n1 = MPI_Wtime()
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

   recursive subroutine BP_MapIntersec2Block(BP, option, stats, msh, ptree, inters, nth, lstr, lstc, lstblk, ll, Nbound)

      implicit none
      type(Hoption)::option
      type(Hstat)::stats
      type(blockplus)::BP
      type(mesh)::msh
      type(proctree)::ptree
      type(intersect)::inters(:)
      integer nth, bidx, level_c
      integer ii, ll, bb, cc, row_group, col_group
      type(list)::lstblk
      type(iarray)::lstr, lstc
      type(nod), pointer::cur
      class(*), pointer::ptr
      type(matrixblock), pointer::blocks
      integer flag, num_nods
      type(block_ptr)::blk_ptr
      real(kind=8)::n1, n0
      integer, allocatable::rowblocks(:), colblocks(:), order_r(:), order_c(:)
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

      allocate (rowblocks(BP%LL(ll)%Nbound))
      allocate (order_r(BP%LL(ll)%Nbound))
      allocate (colblocks(BP%LL(ll)%Nbound))
      allocate (order_c(BP%LL(ll)%Nbound))
      ! row0 = 0
      ! col0 = 0
      ! sort = 0
      do bb = 1, BP%LL(ll)%Nbound
         rowblocks(bb) = BP%LL(ll)%matrices_block(bb)%row_group
         colblocks(bb) = BP%LL(ll)%matrices_block(bb)%col_group
         ! if (rowblocks(bb) < row0 .or. colblocks(bb) < col0) then
         !    sort = 1
         !    exit
         ! endif
         ! row0 = rowblocks(bb)
         ! col0 = colblocks(bb)
      enddo
      call quick_sort_int(rowblocks,order_r,BP%LL(ll)%Nbound)
      call quick_sort_int(colblocks,order_c,BP%LL(ll)%Nbound)

      ! call assert(sort == 0, 'the rowblocks and colblocks need sorting first')

      level0 = BP%LL(1)%matrices_block(1)%level
      level1 = BP%LL(ll)%matrices_block(1)%level

      do ii = 1, lstr%num_nods
         row0 = findgroup(inters(nth)%rows(lstr%dat(ii)), msh, level1 - level0, BP%LL(1)%matrices_block(1)%row_group)
         call binary_search(Nbound, rowblocks, row0, bb)
         if (bb /= -1) then
            do cc=bb,1,-1
               if(rowblocks(cc)==row0)then
                  clstr(order_r(cc))%num_nods = clstr(order_r(cc))%num_nods + 1
                  clstr(order_r(cc))%dat(clstr(order_r(cc))%num_nods) = lstr%dat(ii)
               else
                  exit
               endif
            enddo
            do cc=bb+1,BP%LL(ll)%Nbound
               if(rowblocks(cc)==row0)then
                  clstr(order_r(cc))%num_nods = clstr(order_r(cc))%num_nods + 1
                  clstr(order_r(cc))%dat(clstr(order_r(cc))%num_nods) = lstr%dat(ii)
               else
                  exit
               endif
            enddo
         endif
      enddo

      do ii = 1, lstc%num_nods
         col0 = findgroup(inters(nth)%cols(lstc%dat(ii)), msh, level1 - level0, BP%LL(1)%matrices_block(1)%col_group)
         call binary_search(Nbound, colblocks, col0, bb)
         if (bb /= -1) then
            do cc=bb,1,-1
               if(colblocks(cc)==col0)then
                  clstc(order_c(cc))%num_nods = clstc(order_c(cc))%num_nods + 1
                  clstc(order_c(cc))%dat(clstc(order_c(cc))%num_nods) = lstc%dat(ii)
               else
                  exit
               endif
            enddo
            do cc=bb+1,BP%LL(ll)%Nbound
               if(colblocks(cc)==col0)then
                  clstc(order_c(cc))%num_nods = clstc(order_c(cc))%num_nods + 1
                  clstc(order_c(cc))%dat(clstc(order_c(cc))%num_nods) = lstc%dat(ii)
               else
                  exit
               endif
            enddo
         endif
      enddo

      do bb = 1, Nbound
         blocks => BP%LL(ll)%matrices_block(bb)
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
      deallocate (order_r)
      deallocate (order_c)

      if (ll < BP%Lplus) then
         if (BP%LL(ll + 1)%Nbound > 0) then
            call BP_MapIntersec2Block(BP, option, stats, msh, ptree, inters, nth, lstr, lstc, lstblk, ll + 1, BP%LL(ll + 1)%Nbound)
         endif
      endif

   end subroutine BP_MapIntersec2Block

   subroutine Hmat_MapIntersec2Block(h_mat, option, stats, msh, ptree, inters, nth, lstr, lstc, lstblk, num_blocks)

      implicit none
      type(Hoption)::option
      type(Hstat)::stats
      type(Hmat)::h_mat
      type(mesh)::msh
      type(proctree)::ptree
      type(intersect)::inters(:)
      integer nth, num_blocks
      integer ii, jj, idx, row_group, col_group, i, j
      type(list)::lstblk
      type(iarray)::lstr, lstc
      type(iarray),allocatable:: clstr_g(:), clstc_g(:)

      type(nod), pointer::cur
      class(*), pointer::ptr
      type(matrixblock), pointer::blocks
      type(block_ptr)::blk_ptr

      if(h_mat%myArows>0 .and. h_mat%myAcols>0)then

      allocate(clstr_g(h_mat%myArows))
      do i = 1, h_mat%myArows
         clstr_g(i)%idx = nth
         clstr_g(i)%num_nods = 0
         allocate (clstr_g(i)%dat(lstr%num_nods))
      enddo

      allocate(clstc_g(h_mat%myAcols))
      do j = 1, h_mat%myAcols
         clstc_g(j)%idx = nth
         clstc_g(j)%num_nods = 0
         allocate (clstc_g(j)%dat(lstc%num_nods))
      enddo

      do ii = 1, lstr%num_nods
      do i = 1, h_mat%myArows
         blocks => h_mat%Local_blocks(1, i)
         if (inters(nth)%rows(lstr%dat(ii)) >= msh%basis_group(blocks%row_group)%head .and. inters(nth)%rows(lstr%dat(ii)) <= msh%basis_group(blocks%row_group)%tail) then
            clstr_g(i)%num_nods = clstr_g(i)%num_nods + 1
            clstr_g(i)%dat(clstr_g(i)%num_nods) = lstr%dat(ii)
         endif
      enddo
      enddo


      do jj = 1, lstc%num_nods
      do j = 1, h_mat%myAcols
         blocks => h_mat%Local_blocks(j, 1)
         if (inters(nth)%cols(lstc%dat(jj)) >= msh%basis_group(blocks%col_group)%head .and. inters(nth)%cols(lstc%dat(jj)) <= msh%basis_group(blocks%col_group)%tail) then
            clstc_g(j)%num_nods = clstc_g(j)%num_nods + 1
            clstc_g(j)%dat(clstc_g(j)%num_nods) = lstc%dat(jj)
         endif
      enddo
      enddo

      do i = 1, h_mat%myArows
      if (clstr_g(i)%num_nods > 0) then
      do j = 1, h_mat%myAcols
         if (clstc_g(j)%num_nods > 0) then
            blocks => h_mat%Local_blocks(j, i)
            call Hmat_MapIntersec2Block_Loc(blocks, option, stats, msh, ptree, inters, nth, clstr_g(i), clstc_g(j), lstblk)
         endif
      enddo
      endif
      enddo


      do i = 1, h_mat%myArows
         deallocate (clstr_g(i)%dat)
      enddo
      deallocate(clstr_g)
      do j = 1, h_mat%myAcols
         deallocate (clstc_g(j)%dat)
      enddo
      deallocate(clstc_g)

      endif

   end subroutine Hmat_MapIntersec2Block

   recursive subroutine Hmat_MapIntersec2Block_Loc(blocks, option, stats, msh, ptree, inters, nth, lstr, lstc, lstblk)

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

end module BPACK_constr
