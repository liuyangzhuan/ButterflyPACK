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
!> @file BPACK_defs.f90
!> @brief This file defines all data types, variables and constants used in ButterflyPACK


#include "ButterflyPACK_config.fi"

module BPACK_DEFS

#ifdef MPIMODULE
#ifdef HAVE_MPI
    use MPI
#endif
    use iso_c_binding
#ifdef HAVE_ZFP
    use zfp
#endif
#if __GNUC__ < 5
#else
    use ieee_arithmetic
#endif
    use BPACK_linkedlist
    implicit none
#else


    use iso_c_binding
#ifdef HAVE_ZFP
    use zfp
#endif
#if __GNUC__ < 5
#else
    use ieee_arithmetic
#endif
    use BPACK_linkedlist
    implicit none
#ifdef HAVE_MPI
    INCLUDE 'mpif.h'
#endif
#endif

#ifndef HAVE_MPI
    INCLUDE "mpi_dummy.fi"
#endif


    !>**** the version numbers are automatically replaced with those defined in CMakeList.txt
    integer, parameter:: BPACK_MAJOR_VERSION = 3
    integer, parameter:: BPACK_MINOR_VERSION = 0
    integer, parameter:: BPACK_PATCH_VERSION = 0

    !>**** common parameters
#if defined(PGI) || defined(CRAY)
    integer, external :: iargc
#endif
    real(kind=8), parameter :: BPACK_pi = 4d0*atan(1d0)
    complex(kind=8), parameter :: BPACK_junit = (0d0, 1d0)
    real(kind=8), parameter :: BPACK_Bigvalue = 1d300
    integer, parameter :: BPACK_BigINT = 2147483647
    real(kind=8), parameter:: BPACK_SafeUnderflow = 1D-30
    real(kind=8), parameter:: BPACK_SafeEps = 1D-14
    ! real(kind=8), parameter:: BPACK_Jitter = 1D-5
    DT, parameter :: BPACK_cone = 1d0
    DT, parameter :: BPACK_czero = 0d0
    integer, parameter :: Main_ID = 0 !< Head MPI rank
    integer, parameter :: nbslpk = 32 !< blacs/scalapack block size
    integer, parameter :: Rows_per_processor = 1 !< depreciated
    integer, parameter :: MPI_Header = 11 !< number of integers in the MPI header
    integer, parameter :: INDEX_Header = 4 !< number of integers in header of Butterfly_index_MPI
    integer, parameter :: vec_oversample = 5 !< number of extra vectors adding onto estimated rank in the randomized scheme
    integer, parameter:: msg_chunk = 100000 !< used to determine message tag and hence the massage size
    integer, parameter:: Ndim_max = 100 !< maximum dimension of the inputs

    !>**** parameters for CEM
    real(kind=8), parameter :: BPACK_cd = 299792458d0 !< free-space speed of light
    real(kind=8), parameter :: BPACK_eps0 = 1d7/(4d0*BPACK_pi*BPACK_cd**2) !< free-space permittivity
    real(kind=8), parameter :: BPACK_mu0 = BPACK_pi*4d-7 !< free-space permeability
    real(kind=8), parameter :: BPACK_gamma = 1.781072418d0 !< BPACK_gamma constant
    real(kind=8), parameter :: BPACK_impedence0 = sqrt(BPACK_mu0/BPACK_eps0) !< free-space wave impedance

    !>**** solver parameters
    integer, parameter:: DIRECT = 1         !< use factored HODLR as direct solver
    integer, parameter:: NOPRECON = 2  !< use compressed HODLR as fast matvec
    integer, parameter:: HODLRPRECON = 3        !< use factored HODLR as preconditioner

    integer, parameter:: LplusMax = 10

    integer, parameter:: HODLR = 1  !< use hodlr solver
    integer, parameter:: HMAT = 2  !< use H matrix solver
    integer, parameter:: HSS = 3  !< use hss_bf solver
    integer, parameter:: HSS_MD = 4  !< use hss_bf_md solver

    !>**** construction parameters
    integer, parameter:: SVD = 1
    integer, parameter:: RRQR = 2
    integer, parameter:: ACA = 3
    integer, parameter:: BACA = 4
    integer, parameter:: BACANOVER = 5
    integer, parameter:: PS = 6
    integer, parameter:: ACANMERGE = 7
    integer, parameter:: NATURAL = 0  !< natural order
    integer, parameter:: CKD = 1  !< cartesian kd tree
    integer, parameter:: TM = 2   !< cobble-like ordering
    integer, parameter:: TM_GRAM = 3   !< gram distance-based cobble-like ordering, the behaviour is undefined if matrix is not SPD, HPD, general symmetric or general hermitian

    !>**** hierarchical process grid associated with one process node in the process tree (used for parallel recursive LR compression)
    type grid
        integer :: nsprow=0, nspcol=0 !< number of process rows/columns as square as possible, it's designed that nspcol>=nsprow
        integer :: hprow=-1, hpcol=-1 !< head process in the row/column process tree
        integer :: ctxt=-1 !< blacs context
        integer :: Comm = MPI_COMM_NULL !< MPI communicator for this grid
        integer :: gprow=-1, gpcol=-1 !< the group number in the row and column dimension, no longer needed once constructed
        type(grid), pointer::gdc(:) => null() ! pointer to its two children
    end type grid

    !>**** process groups/nodes in the process tree
    type procgroup
        integer :: nprow=0, npcol=0, nproc=0 !< number of processors and 2D grids
        integer :: head=-1, tail=-1 !< start and end process in the Comm of proctree
        integer :: ctxt=-1 !< blacs context
        integer :: ctxt1D=-1 !< blacs context 1D Row noncyclic (used to distribute from 2D grids to customized noncyclic 1D grid)
        integer :: ctxt1Dcol=-1 !< blacs context 1D Col noncyclic (used to pass into pgemr2df90, otherwise pgemr2d will create an internal 1D col grid, see http://www.netlib.org/scalapack/explore-html/dd/dcd/pdgemr_8c_source.html)
        integer :: ctxt_head=-1 !< blacs context only involving the head process (used to gather and broadcast)
        integer :: Comm = MPI_COMM_NULL  !< MPI communicator for all processes in this node
        ! type(grid),pointer::gd=>null() ! the hierarchical process grid structure associated with each process group
    end type procgroup

    !>**** binary process tree
    type proctree
        integer:: nlevel=-1 !< number of tree levels
        integer :: Comm = MPI_COMM_NULL !< MPI communicator for all processes in this tree
        integer :: nproc = 0 !< # of processes in this tree
        integer :: MyID = 0 !< MPI Rank in Comm
        type(procgroup), allocatable::pgrp(:) !< tree nodes

        DT, allocatable:: send_buff_dat(:), recv_buff_dat(:) !< needs to be moved: communication buffers needed by Bplus_pack_unpack_for_MPI
    end type proctree

    !>**** data packet containing 3 integers and 1 DT
    type dat_pack
        integer::idx(3)
        DT::dat(1)
    end type dat_pack

    !>**** communication buffer for all to all communication
    type commquant1D
        integer:: offset=-1 !< offset in my local array
        integer:: size=0 !< size of the message along first dimension
        integer:: size_i=0 !< size of the message along first dimension
        integer:: active=0 !< whether this communication pair is active
        DT, allocatable::dat(:, :) !< communication buffer
        integer, allocatable::dat_i(:, :) !< communication buffer
        type(dat_pack), allocatable::dat_pk(:, :) !< communication buffer
    end type commquant1D

    !>**** cluster of points/indices
    type basisgroup
        integer:: head=-1 !< head index
        integer:: tail=-1 !< tail index
        ! integer level ! level of this cluster group
        real(kind=8):: radius = 0 !< geomerical radius of this group
        real(kind=8):: boundary(2) = 0 !< seperators used to split this group into children group
        real(kind=8), allocatable:: center(:) !< geometrical center of this group
        integer, allocatable:: nlist(:) !< list of nearfield groups
        integer::nn = 0 !< # of nearfield groups
    end type basisgroup

    !>**** near-interaction list
    type nil_MD
        integer, allocatable:: nlist(:,:) !< list of nearfield groups
        integer::nn = 0 !< # of nearfield groups
    end type nil_MD

    !>**** near-interaction list of all groups at one level
    type nil_onelevel_MD
        integer::len = 0 !< # length of the list
        type(nil_MD),allocatable::list(:) !< list of nlist for each group
    end type nil_onelevel_MD


    !>**** input and output vectors for applying a Bplus
    type vectorsblock
        ! integer style
        ! integer head
        ! integer tail
        DT, allocatable :: vector(:, :)
    end type vectorsblock

    !>**** a vector used to extract one element of a butterfly
    type vectorset
        DT, allocatable :: vector(:)
    end type vectorset

    !>**** information used for one ACA iteration
    type acaquant
        integer M,N !< matrix dimensions
        integer:: header_m=0, header_n=0 !< matrix offsets
        DT, allocatable:: matU(:,:),matV(:,:)
        integer,allocatable:: select_column(:),select_row(:)
        DTR,allocatable:: Singular(:)
        logical:: finish=.False. !< finish flag
        real(kind=8):: normA=0,normUV=0
        integer,allocatable:: rows(:),columns(:)
        integer:: itr=0 !< iteration count
        integer:: rank=0 !< rank of the aca
        integer:: rank0=0 !< the first rank0 columns of matU don't count towards normA
        integer:: itrmax=0 !< max iteration count
    end type acaquant


    !>**** one rank*rank butterfly block
    type butterflymatrix
        DT, pointer :: matrix(:, :)=> null() !< entries of the block
        ! integer::mdim,ndim         ! dimensions of the block
        type(list):: lst !< a list of intersection#s
        integer, allocatable::index(:, :) !< an array of intersection#s
        integer::ndim = 0 !< number of skeletons
        integer,allocatable ::dims_m(:), dims_n(:) !< dimensions of the rows and columns (used in BF_MD)
    end type butterflymatrix

    !>**** index set for one butterfly block
    type butterflyindex
        integer:: size !< length of array
        integer, allocatable :: array(:) !< skeleton columns or rows for one butterfly block
    end type butterflyindex

    !>**** keep track of skeleton columns and rows for one butterfly level
    type butterfly_skel
        integer:: nc = 0, nr = 0 !< # local block rows/columns
        integer:: idx_c = 0, idx_r = 0 !< column and row number of the first local block
        integer:: inc_c = 0, inc_r = 0 !< increment of local block row and columns
        type(butterflyindex), allocatable :: inds(:, :)
    end type butterfly_skel

    !>**** keep track of skeleton columns and rows for one butterfly level
    type butterfly_skel_MD
        integer, allocatable:: nc(:) !< # local block columns per dimension. L: 2^l.  L: 2^l. R: 2^(level_half-l)
        integer, allocatable:: nr(:) !< # local block rows per dimension. L: 2^(level_butterfly-level_half-l). R: 2^l
        integer, allocatable:: idx_c(:), idx_r(:) !< column and row number of the first local block
        integer, allocatable :: inc_c(:), inc_r(:) !< increment of local block row and columns
        type(butterflyindex), allocatable :: inds(:, :, :) !< L: sizes nr(1) * product(nc) * Ndim. !< R: sizes product(nr) * nc(1) * Ndim
    end type butterfly_skel_MD

    !>**** one interior factor
    type butterfly_kerl
        integer num_col !< # block columns
        integer num_row !< # block rows
        integer:: nc = 0, nr = 0 !< # local block rows/columns
        integer:: idx_c = 0, idx_r = 0 !< column and row number of the first local block
        integer:: inc_c = 0, inc_r = 0 !< increment of local block row and columns
        type(butterflymatrix), allocatable :: blocks(:, :)
        type(list):: lst!< a list of active blocks
        integer, allocatable::index(:, :) !< an array of id of active blocks
integer, allocatable::index_MD(:, :, :) !< an array of block offsets
    end type butterfly_kerl

    !>**** one interior factor at level l (L: l=1,level_butterflyL. R: l=1:level_butterflyR)
    type butterfly_kerl_MD
        integer, allocatable:: num_col(:) !< # block columns per dimension. L: 2^l.  L: 2^l. R: 2^(level_half-l+1)
        integer, allocatable:: num_row(:) !< # block rows per dimension. L: 2^(level_butterfly-level_half-l+1). R: 2^l
        integer, allocatable:: nc(:), nr(:) !< # local block rows/columns, should be the same as num_col and num_row
        integer, allocatable :: idx_c(:), idx_r(:) !< column and row number of the first local block
        integer, allocatable :: inc_c(:), inc_r(:) !< increment of local block row and columns
        type(butterflymatrix), allocatable :: blocks(:, :, :) !< L: sizes nr(1) * product(nc) * Ndim. !< R: sizes product(nr) * nc(1) * Ndim
        type(list):: lst!< a list of active blocks
        integer, allocatable::index(:, :, :) !< an array of id of active blocks
    end type butterfly_kerl_MD


    !>**** one outter most factor
    type butterfly_UV
        integer num_blk
        integer:: nblk_loc = 0 !< # local block rows/columns
        integer:: idx = 0 !< column or row number of the first local block
        integer:: inc = 0 !< increment of local block row or columns
        type(butterflymatrix), allocatable :: blocks(:)
    end type butterfly_UV

    !>**** one outter most factor
    type butterfly_UV_MD
        integer :: num_blk
        integer :: nblk_loc !< # local block rows/columns
        integer :: idx  !< column or row number of the first local block
        integer :: inc !< increment of local block row or columns
        type(butterflymatrix), allocatable :: blocks(:,:)
    end type butterfly_UV_MD


    !>**** a derived type containing an integer array
    type:: iarray
        integer:: num_nods = 0 !< length of the array
        integer:: idx = 0 !< user-defined index
        integer, allocatable::dat(:)
    contains
#ifdef HAVE_FINAL
        final :: iarray_finalizer
#endif
    end type iarray

    !>**** intersections of a block row and column
    type:: intersect
        integer::pg !< the index in the process group
        integer::idx
        integer::nc, nr
        integer::nr_loc
        integer, allocatable:: mmap(:) !< map based on masks from the nonzero rows to original rows
        integer, allocatable:: nmap(:) !< map based on masks from the nonzero columns to original columns
        integer, allocatable::masks(:,:) !< store mask of each entry that determine whether to skip the evaluation
        integer, allocatable::rows(:), cols(:) !< store indices in bmat or global list of intersections
        integer, allocatable::rows_loc(:) !< store indices in rows
        integer, allocatable::glo2loc(:) !< store index mapping from rows to rows_loc
        DT, pointer::dat(:, :) => null()
        DT, pointer::dat_loc(:, :) => null()
    end type intersect


    !>**** subtensor of a full tensor
    type:: intersect_MD
        integer::pg !< the index in the process group
        integer::idx
        integer,allocatable::idx_r_m(:),idx_c_m(:)
        integer::sender, receiver
        integer,allocatable::nc(:), nr(:)
        ! integer::nr_loc
        type(iarray), allocatable:: idxmap(:) !< map based on masks from the nonzero indices to original indices
        integer, allocatable::masks(:,:) !< store mask of each entry that determine whether to skip the evaluation
        type(iarray), allocatable::rows(:)  !< store indices in bmat or global list of intersections
        type(iarray), allocatable::cols(:)  !< store indices in bmat or global list of intersections
        ! integer, allocatable::rows_loc(:) !< store indices in rows
        ! integer, allocatable::glo2loc(:) !< store index mapping from rows to rows_loc
        DT, pointer::dat(:, :) => null() !< shape product(nr), product(nc)
        ! DT, pointer::dat_loc(:, :) => null()
    end type intersect_MD

    !>**** ZFP quantity (used for arrays of ZFP compressed data)
    type zfpquant
#ifdef HAVE_ZFP
        type(zFORp_stream) :: stream_r !< ZFP stream for the real part compression
        type(zFORp_stream) :: stream_i !< ZFP stream for the imaginary part compression
#endif
        character, allocatable :: buffer_r(:) !< ZFP buffer for the real part
        character, allocatable :: buffer_i(:) !< ZFP buffer for the imaginary part
    end type zfpquant

    !>**** butterfly or LR structure
    type matrixblock
        integer pgno !< process group
        integer pgno_db !< process group when MPI count is doubled
        integer level !< level in HODLR
        integer col_group !< column group number
        integer row_group !< row group number
        integer style !< 1: full block 2: compressed block 4: hierarchical block
        integer level_butterfly !< butterfly levels
        integer:: level_half = 0 !< the butterfly level where the row-wise and column-wise orderings meet
        integer:: rankmax=0 !< maximum butterfly ranks
        integer:: rankmin=BPACK_BigINT !< minimum butterfly ranks
        integer dimension_rank !< estimated maximum rank
        integer M, N !< size of the block
        integer M_loc, N_loc !< local size of the block
        integer headm, headn !< header indices in row and column dimension
        integer, pointer:: M_p(:, :) => null() !< row sizes of all processes sharing this block
        integer, pointer:: N_p(:, :) => null() !< column sizes of all processes sharing this block
        integer, pointer:: ms(:) => null() !< sizes of accummulated local leaf row blocks
        integer, pointer:: ns(:) => null() !< sizes of accummulated local leaf column blocks
        DT, pointer :: fullmat(:, :) => null() !< full matrix entries
#ifdef HAVE_ZFP
        type(zfpquant):: FullmatZFP !< ZFP quantity for compressing fullmat
#endif
        type(butterfly_UV) :: ButterflyU !< leftmost factor
        type(butterfly_UV) :: ButterflyV !< rightmost factor
        type(butterflymatrix), allocatable :: ButterflyMiddle(:, :) !< middle factor
        type(butterfly_kerl), allocatable :: ButterflyKerl(:) !< interior factors
        type(butterfly_skel), allocatable :: ButterflySkel(:) !< keep track of skeleton columns or rows of each level

        ! the following is for blocks in H matrix solver
        type(matrixblock), pointer :: father => null() !< pointer to its fater
        type(matrixblock), pointer :: sons(:, :) => null() !< pointer to its children
        type(list), allocatable::lstblks(:) !< lstblks(level) is the list of blocks at that level
        ! integer prestyle   !< the block style before the split operation 1: full block 2: compressed block 4: hierarchical block
        ! integer data_type  !< the block data_type, need better documentation later
        ! integer nested_num !< depreciated
        integer, allocatable :: ipiv(:)        !< permutation of the LU of the dense diagonal blocks
        integer blockinfo_MPI(MPI_Header) !< high-level data extracted from the index message: 1. level 2. row_group 3. col_group 4. nested_num(depreciated) 5. style 6. prestyle(depreciated) 7. data_type(depreciated) 8. level_butterfly 9. length_Butterfly_index_MPI 10. length_Butterfly_data_MPI 11. memory (depreciated)
        integer length_Butterfly_index_MPI !< length of the index message, the first INDEX_Header integers are 1. decpreciated 2. rankmax 3. level_butterfly. 4. num_blocks
        integer length_Butterfly_data_MPI !< length of the value message
        DT, allocatable :: fullmat_MPI(:) !< massage for the dense blocks
        integer, allocatable :: Butterfly_index_MPI(:) !< index message the first 4 entries are: 1. depreciated 2. depreciated 3. level_butterfly 4. num_blocks
        DT, allocatable :: Butterfly_data_MPI(:) !< value message
        type(list):: lstr, lstc !< a list of intersections
        type(intersect), allocatable::inters(:) !< an array of intersections
        DT, allocatable:: R(:,:), Rc(:,:), MVP(:,:),MVPc(:,:) !< temporary results for non-transposed and conjugate transposed MVP results and input
    end type matrixblock


    !>**** butterfly or LR structure in tensor format
    type matrixblock_MD
        integer Ndim !< dimensionality
        integer pgno !< process group
        integer pgno_db !< process group when MPI count is doubled
        integer level !< level in HODLR
        integer, allocatable:: col_group(:) !< column group number per dimension
        integer, allocatable:: row_group(:) !< row group number per dimension
        integer style !< 1: full block 2: compressed block 4: hierarchical block
        integer level_butterfly !< butterfly levels
        integer:: level_half = 0 !< the butterfly level where the row-wise and column-wise orderings meet
        integer:: rankmax=0 !< maximum butterfly ranks
        integer:: rankmin=BPACK_BigINT !< minimum butterfly ranks
        integer dimension_rank !< estimated maximum rank
        integer, allocatable:: M(:), N(:) !< size of the block
        integer, allocatable:: M_loc(:), N_loc(:) !< local size of the block
        integer, allocatable:: headm(:), headn(:) !< header indices in row and column dimension
        integer, pointer:: M_p(:, :, :) => null() !< row sizes of all processes sharing this block
        integer, pointer:: N_p(:, :, :) => null() !< column sizes of all processes sharing this block
        integer, pointer:: ms(:,:) => null() !< sizes of accummulated local leaf row blocks
        integer, pointer:: ns(:,:) => null() !< sizes of accummulated local leaf column blocks
        DT, pointer :: fullmat(:, :) => null() !< full matrix entries
#ifdef HAVE_ZFP
        type(zfpquant):: FullmatZFP !< ZFP quantity for compressing fullmat
        type(zfpquant), allocatable :: MiddleZFP(:) ! ZFP quantity array for compressing ButterflyMiddle
#endif
        integer, allocatable::nr_m(:),nc_m(:) !< local number of middle-level row and column groups per dimension. The global number will be nr_m(dim_i)=2^level_half and nc_m(dim_i)=2^(level_butterfly-level_half)
        integer, allocatable:: idx_r_m(:), idx_c_m(:) !< starting index for the middle-level groups per dimension
        type(butterfly_UV_MD), allocatable :: ButterflyU(:) !< leftmost factor of length product(nr_m)
        type(butterfly_UV_MD), allocatable :: ButterflyV(:) !< rightmost factor of length product(nc_m)
        type(butterflymatrix), allocatable :: ButterflyMiddle(:) !< middle factor of sizes product(nr_m) * 2^(level_butterfly-level_half)
        type(butterfly_kerl_MD), allocatable :: ButterflyKerl_L(:,:) !< left half interior factors of sizes product(nr_m) * (level_butterfly-level_half)
        type(butterfly_kerl_MD), allocatable :: ButterflyKerl_R(:,:) !< right half interior factors of sizes product(nr_c) * level_half
        type(butterfly_skel_MD), allocatable :: ButterflySkel_L(:,:) !< keep track of skeleton rows of each level (left half), sizes product(nr_m) * (level_butterfly-level_half+1)
        type(butterfly_skel_MD), allocatable :: ButterflySkel_R(:,:) !< keep track of skeleton columns of each level (right half), sizes product(nr_c) * (level_half+1)

        ! the following is for blocks in H matrix solver
        type(matrixblock_MD), pointer :: father => null() !< pointer to its fater
        type(matrixblock_MD), pointer :: sons(:, :) => null() !< pointer to its children
        ! integer prestyle   !< the block style before the split operation 1: full block 2: compressed block 4: hierarchical block
        ! integer data_type  !< the block data_type, need better documentation later
        ! integer nested_num !< depreciated
        integer, allocatable :: ipiv(:)        !< permutation of the LU of the dense diagonal blocks
        integer blockinfo_MPI(MPI_Header) !< high-level data extracted from the index message: 1. level 2. row_group 3. col_group 4. nested_num(depreciated) 5. style 6. prestyle(depreciated) 7. data_type(depreciated) 8. level_butterfly 9. length_Butterfly_index_MPI 10. length_Butterfly_data_MPI 11. memory (depreciated)
        integer length_Butterfly_index_MPI !< length of the index message, the first INDEX_Header integers are 1. decpreciated 2. rankmax 3. level_butterfly. 4. num_blocks
        integer length_Butterfly_data_MPI !< length of the value message
        DT, allocatable :: fullmat_MPI(:) !< massage for the dense blocks
        integer, allocatable :: Butterfly_index_MPI(:) !< index message the first 4 entries are: 1. depreciated 2. depreciated 3. level_butterfly 4. num_blocks
        DT, allocatable :: Butterfly_data_MPI(:) !< value message
        type(list):: lstr, lstc !< a list of intersections
        type(intersect), allocatable::inters(:) !< an array of intersections
    end type matrixblock_MD




    !>**** one layer in a Bplus
    type onelplus
        integer Nbound !< # of corrected blocks that are further decomposed into deeper layers
        integer rankmax !< maximum butterfly rank on this layer
        type(matrixblock), pointer:: matrices_block(:) => null()
        integer, allocatable::boundary_map(:,:) !< inadmisible subgroups for each subgroup
    end type onelplus


    !>**** one layer in a Bplus
    type onelplus_MD
        integer Nbound !< # of corrected blocks that are further decomposed into deeper layers
        integer rankmax !< maximum butterfly rank on this layer
        type(matrixblock_MD), pointer:: matrices_block(:) => null()
        integer, allocatable::boundary_map(:,:,:) !< inadmisible subgroups for each subgroup
    end type onelplus_MD


    !>**** Bplus structure
    type blockplus
        integer level !< block level in HODLR
        integer col_group !< column group number
        integer row_group !< row group number
        integer pgno   !< process group number
        integer Lplus  !< Number of Bplus layers
        integer ind_ll, ind_bk !< iterator of level and block number in a blockplus
        real(kind=8):: boundary(2) = 0 !< A analytical seperator defined by one out of three coordinates, boundary(1): direction, boundary(2): value
        type(onelplus), pointer:: LL(:) => null() !
    end type blockplus


    !>**** Bplus structure
    type blockplus_MD
        integer level !< block level in HODLR
        integer, allocatable:: col_group(:) !< column group number
        integer, allocatable:: row_group(:) !< row group number
        integer pgno   !< process group number
        integer Lplus  !< Number of Bplus layers
        integer ind_ll, ind_bk !< iterator of level and block number in a blockplus
        type(onelplus_MD), pointer:: LL(:) => null() !
    end type blockplus_MD


    !>**** Structure holding block-diagonal approximations as the Schultz initial guess
    type bdiag
        integer splitlevel  !< level number
        integer N_block !< # of diagonal blocks
        integer Bidxs, Bidxe   !< indice range of my local groups
        type(matrixblock), pointer:: BF_inverse(:) => null() !< inverse blocks
    end type bdiag

    !>**** Structure holding operand in Schulz iteration
    type schulz_operand
        type(bdiag)::bdiags !< block diagonal preconditioner as the initial guess
        type(matrixblock):: matrices_block !< the original butterfly B in I+B
        real(kind=8)::A2norm !< largest singular value of B in I+B
        real(kind=8)::scale !< scalar carried on Identities
        real(kind=8), allocatable::diags(:)
        integer order !< order of schulz iteration
        integer hardstart !< 1: use X0=alphaA^* as the initial guess 0: use block-diagonal approximation of A with recursive inversion as the intial guess
    end type schulz_operand

    !>**** One level in HODLR
    type cascadingfactors
        integer level  !< level number
        integer N_block_forward !< # of forward blocks
        integer N_block_inverse !< # of inverse blocks
        integer Bidxs, Bidxe   !< indice range of my local groups
        type(blockplus), pointer:: BP(:) => null()  !< forward blocks
        type(blockplus), pointer:: BP_inverse(:) => null() !< inverse blocks
        type(blockplus), pointer:: BP_inverse_update(:) => null() !< updated blocks dimension-wise matching forward blocks
        type(blockplus), pointer:: BP_inverse_schur(:) => null() !< schur complement blocks
    end type cascadingfactors

    !>**** HODLR structure
    type hobf
        integer Maxlevel, N !< HODLR levels and sizes
        integer ind_lv, ind_bk !< iterator of level and block number in a HODLR
        type(cascadingfactors), allocatable::levels(:) !
        DT,allocatable::fullmat2D(:,:) !< store the full matrix in 2D block-cyclic fashions
    end type hobf


    !>**** Hmatrix structure
    type Hmat
        integer Maxlevel, N !< H matrix levels and sizes
        integer Dist_level !< used in Hmatrix solver, the level at which parallelization is performed
        integer idxs, idxe !< same as msh%idxs and msh%idxe
        integer, pointer:: N_p(:, :) => null() !< row sizes of all processes sharing this Hmat
        type(basisgroup), allocatable:: basis_group(:) !< basis_group at the Dist_level level, this is needed as msh is not currently passed to matvec interface
        integer:: myArows=0,myAcols=0 !< local number of row and column blocks
        type(matrixblock), pointer :: Local_blocks(:, :) => null()
        type(matrixblock), pointer :: Local_blocks_copy(:, :) => null() !< copy of the forward matrix
      type(matrixblock), pointer :: Computing_matricesblock_m(:, :) => null(), Computing_matricesblock_l(:, :) => null(), Computing_matricesblock_u(:, :) => null()
        type(matrixblock), pointer:: blocks_1 => null(), blocks_2 => null()
        type(list), allocatable::lstblks(:) !< lstblks(level) is the list of blocks at that level
        type(list),allocatable::admissibles(:) !< a list of admissible and inadmissible groups per each group
        type(iarray),allocatable::colorsets(:) !< the colorset (an integer array) of each level
        DT,allocatable::fullmat(:,:) !< store the full matrix for debugging purpose
        DT,allocatable::fullmat2D(:,:) !< store the full matrix in 2D block-cyclic fashions
    end type Hmat

    !>**** HSS structure
    type hssbf
        integer Maxlevel, N !< HSS levels and sizes
        ! ! integer ind_lv,ind_bk ! iterator of level and block number in a HODLR
        type(blockplus)::BP, BP_inverse !< a single butterfly plus for the entire matrix
        DT,allocatable::fullmat2D(:,:) !< store the full matrix in 2D block-cyclic fashions
    end type hssbf

    !>**** HSS_MD structure
    type hssbf_md
        integer Maxlevel, Ndim !< HSS_MD levels and dimensionility
        integer, allocatable:: N(:) !< size per dimension
        ! ! integer ind_lv,ind_bk ! iterator of level and block number in a HODLR
        type(blockplus_MD)::BP !< a single butterfly plus for the entire matrix
    end type hssbf_md


    type SVD_quant
        DT, allocatable:: matU(:,:),matV(:,:)
        DTR,allocatable:: Singular(:)
    end type SVD_quant

    type Bmatrix
        integer Maxlevel
        type(hobf), pointer::ho_bf => null()
        type(Hmat), pointer::h_mat => null()
        type(hssbf), pointer::hss_bf => null()
        type(hssbf_md), pointer::hss_bf_md => null()
    end type Bmatrix

    !>**** intermidate vectors for applying a butterfly
    type RandomBlock
        integer level_butterfly
        type(butterfly_kerl), allocatable :: RandomVectorRR(:)
        type(butterfly_kerl), allocatable :: RandomVectorLL(:)
    end type RandomBlock

    !>**** intermidate vectors for applying a butterfly
    type butterfly_vec
        type(butterfly_kerl), allocatable :: vec(:)
    end type butterfly_vec

    !>**** HODLR solver options
    type Hoption

        integer::format !< HODLR or HMAT or format
        integer::verbosity !< printlevel -1: no printing except error and warning. 0: default printing. 1: print info for each high-level operation 2: print information for each low-level operation

        ! options for Bplus, Butterfly or LR
        integer::LRlevel  !< The top LRlevel level blocks are butterfly or Bplus
        integer:: lnoBP !< the bottom lnoBP levels are either Butterfly or LR, but not Bplus
        integer:: bp_cnt_lr !< only print the rank in the top-layer butterfly of a Bplus
        integer:: TwoLayerOnly  !< restrict Bplus as Butterfly + LR
        real(kind=8) touch_para   !< parameters used to determine whether one patch is closer to seperator
        real(kind=8) sample_para   !< parameters used for linear-complexity ID-butterfly, # of row/columns samples is sample_para*2*butterfly_rank
        real(kind=8) sample_para_outer   !< parameters used for linear-complexity ID-butterfly, # of row/columns samples is sample_para*2*butterfly_rank
        ! integer sample_heuristic   !< 1: use skeleton rows/columns from the previous block during BF compression assuming they should share similar skeletons
        integer:: pat_comp !< pattern of entry-evaluation-based butterfly compression: 1 from right to left, 2 from left to right, 3 from outer to inner
        integer:: use_zfp  !< 1: use zfp for the dense blocks (zfp must be used to install ButterflyPACK) 0: do not use zfp

        ! options for matrix construction
        integer forwardN15flag !< 1 use N^1.5 algorithm. 0: use NlogN pseudo skeleton algorithm. 2: use NlogN first, if not accurate enough, switch to N^1.5
        real(kind=8) tol_comp      !< matrix construction tolerance
        integer::Nmin_leaf !< leaf sizes of HODLR tree
        integer nogeo !< 1: no geometrical information available to hodlr, use NATUTAL or TM_GRAM clustering        0: geometrical points are available for TM or CKD clustering 2: no geometrical information available, but a user-defined distance function and compressibility function is provided. 3: no geometrical information available, but an array of knn*N indicating the knn neighbours of each element is provided. 4: geometrical information available for TM or CKD clustering, and an array of knn*N indicating the knn neighbours of each element is provided
        integer per_geo !< 1: the geomerical points are periodical. 0: the points are not periodical
        real(kind=8):: periods(Ndim_max) !< the periods in each dimension (currently only supports maximum of 3 dimensions) of the geometry points when per_geo=1
        integer xyzsort !< clustering methods given geometrical points: CKD: cartesian kd tree SKD: spherical kd tree (only for 3D points) TM: (2 mins no recursive)
        integer::RecLR_leaf !< bottom level operations in a recursive merge-based LR compression: SVD, RRQR, ACA, BACA
        real(kind=8):: near_para !< parameters used to determine whether two groups are nearfield or farfield pair
        real(kind=8):: knn_near_para !< parameters used to determine whether two groups are nearfield or farfield pair, used for knn search
        real(kind=8):: scale_factor !< parameters used to scale matrix entries
        integer::rmax !< maximum rank truncation
        integer:: elem_extract !< 1: use user-defined element extraction 0: use user-defined formula
        integer:: cpp !< 1: use user-defined c/cpp functions 0: use user-defined fortran functions
        integer:: knn !< number of nearest neighbour points for each point
integer:: fastsample_tensor !< 0: uniformly sample each dimension. 1: uniformly sample the rows of the unfolded matrices on top of 0. 2: use translation invariance

        ! options for inversion
        real(kind=8) tol_LS       !< tolerance in pseudo inverse
        real(kind=8) tol_Rdetect  !< tolerance to detect numerical ranks
        real(kind=8) tol_rand     !< tolerance for randomized contruction
        real(kind=8) jitter     !< jittering for dense diagonal blocks
        integer powiter     !< order of power iteration in randomized LR
        integer less_adapt     !< 0 for rank adaptation for all BF levels, 1 for rank adaptation for the outtermost BF levels
        integer::schulzorder !< order (>=2) of schultz iteration
        integer::schulzhardstart !< 1: use X0=alphaA^* as the initial guess 0: use block-diagonal approximation of A with recursive inversion as the intial guess
        integer::schulzsplitlevel !< number of levels to split A for block-diagonal approximation
        integer::schulzlevel !< (I+B)^-1 is computed by schultz iteration for butterfly with more than schulzlevel levels
        integer::rank0 !< intial guess of ranks
        real(kind=8)::rankrate !< increasing rate of rank estimates per iteration
        integer::itermax !< max number of iteration in randomized schemes
        integer::ILU !< only perform LU on dense diagonal blocks, used only in the context of H-LU
        integer::Nbundle !< multiply Nbundle sets of vectors together in randomized schemes

        ! options for solve
        real(kind=8) tol_itersol  !< tolerance for iterative solvers
        integer n_iter  !< maximum number of iterations for iterative solver
        integer precon  !< DIRECT: use factored HODLR as direct solver, HODLRPRECON: use factored HODLR as preconditioner, NOPRECON: use forward HODLR as fast matvec,

        ! options for error checking
        integer::level_check !< check compression quality by picking random entries at level_check (only work for nmpi=1 now)
        integer::ErrFillFull !< check compression quality by computing all block elements
        integer::ErrSol !< check solution quality by using artificially generated true solution vector
        integer::BACA_Batch !< batch size in batch ACA
        integer::LR_BLK_NUM !< sqrt of number of bottom-level subblocks in blocked LR

    end type Hoption

    !>**** statistics
    type Hstat
        real(kind=8) Time_random(5)  !< Intialization, MVP, Reconstruction, Reconstruction of one subblock
        real(kind=8) Time_Sblock, Time_Inv, Time_SMW, Time_PartialUpdate, Time_Fill, Time_RedistB, Time_RedistV, Time_Sol, Time_BLK_MVP, Time_C_Mult, Time_C_Extract, Time_Entry, Time_Entry_Traverse, Time_Entry_BF, Time_Entry_Comm
        real(kind=8) Time_Direct_LU, Time_Add_Multiply, Time_Multiply, Time_XLUM, Time_Split, Time_Comm, Time_Idle, Time_Factor
        real(kind=8) Mem_Current, Mem_peak, Mem_Sblock, Mem_SMW, Mem_Direct_inv, Mem_Direct_for, Mem_int_vec, Mem_Comp_for, Mem_Fill, Mem_Factor
        real(kind=8) Flop_Fill, Flop_Factor, Flop_Sol, Flop_C_Mult, Flop_C_Extract, Flop_Tmp
        integer, allocatable:: rankmax_of_level(:), rankmin_of_level(:) !< maximum and minimum ranks observed at each level of HODLR during matrix fill
        integer, allocatable:: rankmax_of_level_global(:) !< maximum ranks among all processes observed at each level of HODLR during matrix fill
        integer, allocatable:: rankmax_of_level_global_factor(:) !< maximum ranks among all processes observed at each level of HODLR during matrix factorization
        integer, allocatable:: Add_random_CNT(:), Mul_random_CNT(:), XLUM_random_CNT(:) !< record number of randomized operations
        real(kind=8), allocatable:: Add_random_Time(:), Mul_random_Time(:), XLUM_random_Time(:) !< record number of randomized operations
        integer, allocatable:: leafs_of_level(:) !< number of compressed blocks at each level
    end type Hstat

    !>**** quantities related to geometries, meshes, unknowns and points
    type mesh
        integer Nunk !< size of the matrix
        integer Dist_level !< used in Hmatrix solver, the level at which parallelization is performed
        integer Maxgroup !< number of nodes in the partition tree
        integer idxs, idxe  !< range of local row/column indices after reordering
        real(kind=8), allocatable:: xyz(:, :)   !< coordinates of the points
        integer, allocatable:: new2old(:) !< index mapping from new ordering to old ordering
        integer, allocatable:: old2new(:) !< index mapping from old ordering to new ordering
        integer, allocatable::pretree(:) !< dimension 2**Maxlevel containing box size of each leaf node
        integer, allocatable::nns(:, :) !< index (after permutation) of k-nearest neighbours for each point
        type(basisgroup), allocatable:: basis_group(:)
    end type mesh

    !>**** quantities related to specific matrix kernels
    type kernelquant

        DT, allocatable :: matZ_glo(:, :) !< Full Matrix: the full matrix to sample its entries

        class(*), pointer :: QuantApp => null() !< Kernels Defined in Fortran: pointer to the user-supplied derived type for computing one element of Z
        procedure(F_Zelem), nopass, pointer :: FuncZmn => null() !< Kernels Defined in Fortran: procedure pointer to the user-supplied derived type for computing one element of Z
        procedure(F_Zelem_MD), nopass, pointer :: FuncZmn_MD => null() !< Kernels Defined in Fortran: procedure pointer to the user-supplied derived type for computing one element of Z (as tensors)
        procedure(F_Dist), nopass, pointer :: FuncDistmn => null() !< Kernels Defined in Fortran: procedure pointer to the user-supplied derived type for computing distance between one row and one column
        procedure(F_Compressibility), nopass, pointer :: FuncNearFar => null() !< Kernels Defined in Fortran: procedure pointer to the user-supplied derived type for determining whether a block in Z (after permutation) is compressible or not
        procedure(F_Zelem_block), nopass, pointer :: FuncZmnBlock => null() !< Kernels Defined in Fortran: procedure pointer to the user-supplied derived type for computing a list of intersection of indices from Z (data layout needs to be provided)
        procedure(F_HMatVec), nopass, pointer :: FuncHMatVec => null() !< Kernels Defined in Fortran: procedure pointer to the user-supplied derived type for computing matvec of Z

        type(c_ptr), pointer :: C_QuantApp => null() !< Kernels Defined in C: c_pointer to the user-supplied object for computing one element of Z
        type(c_funptr), pointer :: C_FuncZmn => null() !< Kernels Defined in C: c_function_pointer to the user-supplied function for computing one element of Z
        type(c_funptr), pointer :: C_FuncDistmn => null() !< Kernels Defined in C: c_function_pointer to the user-supplied function for computing distance between one row and one column
        type(c_funptr), pointer :: C_FuncNearFar => null() !< Kernels Defined in C: c_function_pointer to the user-supplied function for determine whether a block in Z is compressible or not
        type(c_funptr), pointer :: C_FuncZmnBlock => null() !< Kernels Defined in C: c_function_pointer to the user-supplied function for computing a list of intersection of indices from Z (data layout needs to be provided)
        type(c_funptr), pointer :: C_FuncHMatVec => null() !< Kernels Defined in C: procedure pointer to the user-supplied derived type for computing matvec of Z
        type(c_funptr), pointer :: C_FuncBMatVec => null() !< Kernels Defined in C: procedure pointer to the user-supplied derived type for computing matvec of a block
    end type kernelquant

    type quant_bmat
        type(Bmatrix), pointer::bmat !< Use this metadata in blocked element extraction
        type(mesh), pointer::msh   !< Use this metadata in blocked element extraction
        type(proctree), pointer::ptree !< Use this metadata in blocked element extraction
        type(Hstat), pointer::stats !< Use this metadata in blocked element extraction
        type(Hoption), pointer::option !< Use this metadata in blocked element extraction
        type(kernelquant), pointer::ker !< Use this metadata in blocked element extraction
    end type quant_bmat

    !>**** a derived type containing a pointer to a block
    type:: block_ptr
        type(matrixblock), pointer::ptr
    end type block_ptr

    !>*** a derived type for a pair of integers
    type:: ipair
        integer i, j
    end type ipair

    abstract interface
        subroutine BMatVec(operand, block_o, trans, M, N, num_vect_sub, Vin, ldi, Vout, ldo, a, b, ptree, stats, operand1)
            import :: matrixblock, proctree, Hstat
            implicit none
            class(*)::operand
            class(*), optional::operand1
            type(matrixblock)::block_o
            character trans
            integer M, N, num_vect_sub
            type(proctree)::ptree
            type(Hstat)::stats
            integer ldi, ldo
            DT :: Vin(ldi, *), Vout(ldo, *), a, b
        end subroutine BMatVec

        subroutine Zelem(edge_m, edge_n, value, msh, option, ker)
            import::mesh, Hoption, kernelquant
            implicit none
            integer edge_m, edge_n
            DT value
            type(mesh)::msh
            type(Hoption)::option
            type(kernelquant)::ker
        end subroutine Zelem

        subroutine Zelem_block(nrow, ncol, mrange, nrange, values, msh, option, ker, myflag, passflag, ptree, stats)
            import::mesh, Hoption, kernelquant, proctree, Hstat
            implicit none
            integer nrow, ncol, myflag, passflag
            integer mrange(nrow)
            integer nrange(ncol)
            DT:: values(nrow, ncol)
            type(mesh)::msh
            type(Hoption)::option
            type(Hstat)::stats
            type(kernelquant)::ker
            type(proctree)::ptree
        end subroutine Zelem_block

        !> interface of user-defined element evaluation routine in Fortran. m,n represents indices in natural order
        subroutine F_Zelem(m, n, val, quant)
            import::mesh, kernelquant
            class(*), pointer :: quant
            integer, INTENT(IN):: m, n
            DT::val
        end subroutine F_Zelem


        !> interface of user-defined element evaluation routine in Fortran. m,n represents (multi-dimensional) indices in natural order
        subroutine F_Zelem_MD(Ndim, m, n, val, quant)
            import::mesh, kernelquant
            class(*), pointer :: quant
            integer, INTENT(IN):: Ndim
            integer, INTENT(IN):: m(Ndim), n(Ndim)
            DT::val
        end subroutine F_Zelem_MD


        !> interface of user-defined distance computation routine in Fortran. m,n represents indices in natural order
        subroutine F_Dist(m, n, val, quant)
            import::mesh, kernelquant
            class(*), pointer :: quant
            integer, INTENT(IN):: m, n
            real(kind=8)::val
        end subroutine F_Dist

        !> interface of user-defined compressibility routine in Fortran. groupm,groupn represents groups in the permuted order
        subroutine F_Compressibility(groupm, groupn, val, quant)
            import::mesh, kernelquant
            class(*), pointer :: quant
            integer, INTENT(IN):: groupm, groupn
            integer::val
        end subroutine F_Compressibility

        !> interface of user-defined element extraction routine in Fortran. allrows,allcols represents indices in natural order
        subroutine F_Zelem_block(Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps, quant)
            class(*), pointer :: quant
            integer:: Ninter
            integer:: allrows(:), allcols(:)
            DT,target::alldat_loc(:)
            integer::colidx(Ninter), rowidx(Ninter), pgidx(Ninter)
            integer::Npmap, pmaps(Npmap, 3)
        end subroutine F_Zelem_block

        !> interface of user-defined element evaluation routine in C. m,n represents indices in natural order
        subroutine C_Zelem(m, n, val, quant)
            USE, INTRINSIC :: ISO_C_BINDING
            type(c_ptr) :: quant
            integer(kind=C_INT), INTENT(IN):: m, n
            CBIND_DT::val
        end subroutine C_Zelem

        !> interface of user-defined distance computation routine in C. m,n represents indices in natural order
        subroutine C_Dist(m, n, val, quant)
            USE, INTRINSIC :: ISO_C_BINDING
            type(c_ptr) :: quant
            integer(kind=C_INT), INTENT(IN):: m, n
            real(kind=C_DOUBLE)::val
        end subroutine C_Dist

        !> interface of user-defined distance compressibility routine in C. groupm,groupn represents groups in the permuted order
        subroutine C_Compressibility(groupm, groupn, val, quant)
            USE, INTRINSIC :: ISO_C_BINDING
            type(c_ptr) :: quant
            integer(kind=C_INT), INTENT(IN):: groupm, groupn
            integer(kind=C_INT)::val
        end subroutine C_Compressibility

        !> interface of user-defined element extraction routine in C. allrows,allcols represents indices in natural order
        subroutine C_Zelem_block(Ninter, Nallrows, Nallcols, Nalldat_loc, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps, quant)
            USE, INTRINSIC :: ISO_C_BINDING
            type(c_ptr) :: quant
            integer(kind=C_INT64_T):: Nalldat_loc
            integer(kind=C_INT):: Ninter, Nallrows, Nallcols
            integer(kind=C_INT):: allrows(Nallrows), allcols(Nallcols)
            CBIND_DT,target::alldat_loc(Nalldat_loc)
            integer(kind=C_INT)::colidx(Ninter), rowidx(Ninter), pgidx(Ninter)
            integer(kind=C_INT)::Npmap, pmaps(Npmap, 3)
        end subroutine C_Zelem_block

        subroutine HMatVec(trans, M, N, num_vect, Vin, Vout, ker)
            import::mesh, kernelquant, proctree, Hstat
            implicit none
            integer, INTENT(IN):: M, N, num_vect
            DT::Vin(:, :), Vout(:, :)
            type(kernelquant)::ker
            character trans
        end subroutine HMatVec

        !> interface of user-defined HODLR/H MatVec routine in Fortran.
        subroutine F_HMatVec(trans, M, N, num_vect, Vin, Vout, quant)
            import::mesh, proctree, Hstat
            class(*), pointer :: quant
            integer, INTENT(IN):: M, N, num_vect
            DT::Vin(:, :), Vout(:, :)
            character trans
        end subroutine F_HMatVec

        !> interface of user-defined HODLR/H MatVec routine in C.
        subroutine C_HMatVec(trans, Nin, Nout, num_vect, Vin, Vout, quant)
            USE, INTRINSIC :: ISO_C_BINDING
            import::mesh, proctree, Hstat
            type(c_ptr) :: quant
            integer(kind=C_INT), INTENT(IN):: Nin, Nout, num_vect
            CBIND_DT::Vin(Nin, num_vect), Vout(Nout, num_vect)
            type(mesh)::msh
            type(proctree)::ptree
            type(Hstat)::stats
            character(kind=c_char, len=1) :: trans(*)
        end subroutine C_HMatVec

        !> interface of user-defined Block MatVec routine in C.
        subroutine C_BMatVec(trans, Nin, Nout, num_vect, Vin, Vout, quant, a, b)
            USE, INTRINSIC :: ISO_C_BINDING
            type(c_ptr) :: quant
            character(kind=c_char, len=1) :: trans(*)
            integer(kind=C_INT), INTENT(IN):: Nin, Nout, num_vect
            CBIND_DT :: Vin(Nin, num_vect), Vout(Nout, num_vect), a, b
        end subroutine C_BMatVec

    end interface

    !>*** need to further work on the following:
    real(kind=8)::time_tmp,time_tmp1,time_tmp2,time_tmp3,time_tmp4,time_tmp5

    ! integer,allocatable:: basis_group_pre(:,:)
contains

    !> deallocate all entries of an iarray
    subroutine iarray_finalizer(lst)
        implicit none
        type(iarray):: lst
        if (allocated(lst%dat)) deallocate (lst%dat)
        lst%idx = 0
        lst%num_nods = 0
    end subroutine iarray_finalizer
end module BPACK_DEFS
