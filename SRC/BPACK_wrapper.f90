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
!> @file BPACK_wrapper.f90
!> @brief C++ interfaces for the high-level Fortran subroutines using iso_c_binding

#include "ButterflyPACK_config.fi"
module BPACK_wrapper
   use BPACK_DEFS
   use BPACK_structure
   use BPACK_factor
   use BPACK_constr
   use BPACK_randomMVP
#ifdef HAVE_OPENMP
   use omp_lib
#endif
   use MISC_Utilities
   use BPACK_Solve_Mul
   use iso_c_binding

contains

!>****** Fortran interface for the matvec function required by BPACK_construction_Matvec, inside which a c++ function pointer ker%C_FuncHMatVec is called \n
   !>******! It is assumed the ker%C_FuncHMatVec interfaces with local input and output vectors, which assume an already ordered hierarchical matrix
   !> @param ker: the structure containing kernel quantities
   !> @param Vin: input vector
   !> @param Vout: output vector
   !> @param M: (local) row dimension of the matrix
   !> @param N: (local) column dimension of the matrix
   !> @param trans: 'N', 'C' or 'T'
   !> @param num_vect: number of vectors
   subroutine matvec_user_C(trans, M, N, num_vect, Vin, Vout, ker)

      integer, INTENT(IN):: M, N, num_vect
      DT::Vin(:, :), Vout(:, :)
      type(kernelquant)::ker
      character trans

      procedure(C_HMatVec), POINTER :: proc

      call c_f_procpointer(ker%C_FuncHMatVec, proc)
      if (trans == 'N') then  ! note that C_HMatVec takes Nin Nout, instead of takes M,N, so if statement is needed here
         call proc(trans//c_null_char, N, M, num_vect, Vin, Vout, ker%C_QuantApp)
      else
         call proc(trans//c_null_char, M, N, num_vect, Vin, Vout, ker%C_QuantApp)
      endif
      return

   end subroutine matvec_user_C



!>****** Fortran interface for the matvec function required by BF_randomized, inside which a c++ function pointer ker%C_FuncBMatVec is called \n
   !>******! It is assumed the ker%C_FuncBMatVec does not need ldi and ldo! \n
   !>******! It is assumed the ker%C_FuncBMatVec interfaces with local input and output vectors, which assume already ordered rows/columns
   !> @param ker: the structure containing kernel quantities
   !> @param block_o: not referenced
   !> @param trans: 'N', 'C' or 'T'
   !> @param M: (local) row dimension of the block
   !> @param N: (local) column dimension of the block
   !> @param num_vect: number of vectors
   !> @param Vin: input vector
   !> @param ldi: leading dimension of Vin, needs to be M or N depending on trans
   !> @param Vout: output vector
   !> @param ldo: leading dimension of Vout, needs to be M or N depending on trans
   !> @param a: Vout = a*A*Vin + b*Vout
   !> @param b: Vout = a*A*Vin + b*Vout
   !> @param ptree: the structure containing process tree
   !> @param stats: the structure containing statistics
   !> @param operand1: not referenced
   subroutine Bmatvec_user_C(ker, block_o, trans, M, N, num_vect, Vin, ldi, Vout, ldo, a, b, ptree, stats, operand1)
      implicit none
      class(*)::ker
      class(*), optional::operand1
      type(matrixblock)::block_o
      character trans
      integer M, N, num_vect
      type(proctree)::ptree
      type(Hstat)::stats
      integer ldi, ldo
      DT :: Vin(ldi, *), Vout(ldo, *), a, b
      procedure(C_BMatVec), POINTER :: proc

      select TYPE (ker)
      type is (kernelquant)
         call c_f_procpointer(ker%C_FuncBMatVec, proc)
         if (trans == 'N') then  ! note that C_HMatVec takes Nin Nout, instead of takes M,N, so if statement is needed here
            call proc(trans//c_null_char, N, M, num_vect, Vin, Vout, ker%C_QuantApp, a, b)
         else
            call proc(trans//c_null_char, M, N, num_vect, Vin, Vout, ker%C_QuantApp, a, b)
         endif
      end select
   end subroutine Bmatvec_user_C

!>**** C interface of process tree construction
   !> @param nmpi: number of MPIs for one hodlr
   !> @param MPIcomm: MPI communicator from C caller
   !> @param groupmembers: MPI ranks in MPIcomm for one hodlr
   !> @param ptree_Cptr: the structure containing process tree
   subroutine C_BPACK_Createptree(nmpi, groupmembers, MPIcomm, ptree_Cptr) bind(c, name="c_bpack_createptree")
      implicit none
      integer nmpi
      integer MPIcomm
      integer:: groupmembers(*)
      type(c_ptr):: ptree_Cptr
      type(proctree), pointer::ptree

      allocate (ptree)
      call CreatePtree(nmpi, groupmembers, MPIcomm, ptree)
      ptree_Cptr = c_loc(ptree)
   end subroutine C_BPACK_Createptree

!>**** C interface of initializing statistics
   !> @param stats_Cptr: the structure containing statistics
   subroutine C_BPACK_Createstats(stats_Cptr) bind(c, name="c_bpack_createstats")
      implicit none
      type(c_ptr), intent(out) :: stats_Cptr
      type(Hstat), pointer::stats

      allocate (stats)
      !>**** initialize statistics variables
      call InitStat(stats)
      stats_Cptr = c_loc(stats)

   end subroutine C_BPACK_Createstats

!>**** C interface of getting one entry in stats
   !> @param stats_Cptr: the structure containing stats
   !> @param nam: name of the stats
   !> @param val_d: value of the stats
   subroutine C_BPACK_Getstats(stats_Cptr, nam, val_d) bind(c, name="c_bpack_getstats")
      implicit none
      real(kind=8)::val_d
      character(kind=c_char, len=1) :: nam(*)
      character(kind=c_char) :: tmpc
      type(c_ptr) :: stats_Cptr
      type(Hstat), pointer::stats
      ! character::nam(:)
      ! type(c_ptr),value :: val_Cptr
      ! integer,pointer::val_i
      ! real(kind=8),pointer::val_d
      integer strlen
      character(len=:), allocatable :: str
      integer valid_opt

      valid_opt = 0
      strlen = 1
      tmpc=nam(strlen)
      do while (tmpc /= c_null_char)
         strlen = strlen + 1
         tmpc=nam(strlen)
      enddo
      strlen = strlen - 1
      allocate (character(len=strlen) :: str)
      str = transfer(nam(1:strlen), str)

      call c_f_pointer(stats_Cptr, stats)

      if (trim(str) == 'Time_Fill') then
         val_d = stats%Time_Fill
         valid_opt = 1
      endif
      if (trim(str) == 'Time_Entry') then
         val_d = stats%Time_Entry
         valid_opt = 1
      endif
      if (trim(str) == 'Time_Factor') then
         val_d = stats%Time_Factor
         valid_opt = 1
      endif
      if (trim(str) == 'Time_Solve') then
         val_d = stats%Time_Sol
         valid_opt = 1
      endif
      if (trim(str) == 'Time_Sblock') then
         val_d = stats%Time_Sblock
         valid_opt = 1
      endif
      if (trim(str) == 'Time_Inv') then
         val_d = stats%Time_Inv
         valid_opt = 1
      endif
      if (trim(str) == 'Time_SMW') then
         val_d = stats%Time_SMW
         valid_opt = 1
      endif
      if (trim(str) == 'Time_PartialUpdate') then
         val_d = stats%Time_PartialUpdate
         valid_opt = 1
      endif
      if (trim(str) == 'Time_RedistB') then
         val_d = stats%Time_RedistB
         valid_opt = 1
      endif
      if (trim(str) == 'Time_RedistV') then
         val_d = stats%Time_RedistV
         valid_opt = 1
      endif
      if (trim(str) == 'Time_BLK_MVP') then
         val_d = stats%Time_BLK_MVP
         valid_opt = 1
      endif
      if (trim(str) == 'Time_C_Mult') then
         val_d = stats%Time_C_Mult
         valid_opt = 1
      endif
      if (trim(str) == 'Time_C_Extract') then
         val_d = stats%Time_C_Extract
         valid_opt = 1
      endif
      if (trim(str) == 'Time_Entry_Traverse') then
         val_d = stats%Time_Entry_Traverse
         valid_opt = 1
      endif
      if (trim(str) == 'Time_Entry_BF') then
         val_d = stats%Time_Entry_BF
         valid_opt = 1
      endif
      if (trim(str) == 'Time_Entry_Comm') then
         val_d = stats%Time_Entry_Comm
         valid_opt = 1
      endif
      if (trim(str) == 'Time_Direct_LU') then
         val_d = stats%Time_Direct_LU
         valid_opt = 1
      endif
      if (trim(str) == 'Time_Add_Multiply') then
         val_d = stats%Time_Add_Multiply
         valid_opt = 1
      endif
      if (trim(str) == 'Time_Multiply') then
         val_d = stats%Time_Multiply
         valid_opt = 1
      endif
      if (trim(str) == 'Time_XLUM') then
         val_d = stats%Time_XLUM
         valid_opt = 1
      endif
      if (trim(str) == 'Time_Split') then
         val_d = stats%Time_Split
         valid_opt = 1
      endif
      if (trim(str) == 'Time_Comm') then
         val_d = stats%Time_Comm
         valid_opt = 1
      endif
      if (trim(str) == 'Time_Idle') then
         val_d = stats%Time_Idle
         valid_opt = 1
      endif

      if (trim(str) == 'Flop_Fill') then
         val_d = stats%Flop_Fill
         valid_opt = 1
      endif
      if (trim(str) == 'Flop_Factor') then
         val_d = stats%Flop_Factor
         valid_opt = 1
      endif
      if (trim(str) == 'Flop_Solve') then
         val_d = stats%Flop_Sol
         valid_opt = 1
      endif
      if (trim(str) == 'Flop_C_Mult') then
         val_d = stats%Flop_C_Mult
         valid_opt = 1
      endif
      if (trim(str) == 'Flop_C_Extract') then
         val_d = stats%Flop_C_Extract
         valid_opt = 1
      endif

      if (trim(str) == 'Mem_Factor') then
         val_d = stats%Mem_Factor
         valid_opt = 1
      endif
      if (trim(str) == 'Mem_Fill') then
         val_d = stats%Mem_Fill
         valid_opt = 1
      endif
      if (trim(str) == 'Mem_Sblock') then
         val_d = stats%Mem_Sblock
         valid_opt = 1
      endif
      if (trim(str) == 'Mem_SMW') then
         val_d = stats%Mem_SMW
         valid_opt = 1
      endif
      if (trim(str) == 'Mem_Direct_inv') then
         val_d = stats%Mem_Direct_inv
         valid_opt = 1
      endif
      if (trim(str) == 'Mem_Direct_for') then
         val_d = stats%Mem_Direct_for
         valid_opt = 1
      endif
      if (trim(str) == 'Mem_int_vec') then
         val_d = stats%Mem_int_vec
         valid_opt = 1
      endif
      if (trim(str) == 'Mem_Comp_for') then
         val_d = stats%Mem_Comp_for
         valid_opt = 1
      endif
      if (trim(str) == 'Mem_Peak') then
         val_d = stats%Mem_Peak
         valid_opt = 1
      endif
      if (trim(str) == 'Rank_max') then
         val_d = dble(maxval(stats%rankmax_of_level_global))
         if(allocated(stats%rankmax_of_level_global_factor))val_d = NINT(max(dble(val_d),dble(maxval(stats%rankmax_of_level_global_factor))))
         valid_opt = 1
      endif

      if (valid_opt == 0) write (*, *) 'invalid BPACK stats: '//trim(str)

      deallocate (str)

   end subroutine C_BPACK_Getstats

!>**** C interface of printing statistics
   !> @param stats_Cptr: the structure containing statistics
   !> @param ptree_Cptr: the structure containing process tree
   subroutine C_BPACK_Printstats(stats_Cptr, ptree_Cptr) bind(c, name="c_bpack_printstats")
      implicit none
      type(c_ptr) :: stats_Cptr
      type(c_ptr) :: ptree_Cptr
      type(Hstat), pointer::stats
      type(proctree), pointer::ptree

      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)
      !>**** print statistics variables
      call PrintStat(stats, ptree)

   end subroutine C_BPACK_Printstats

!>**** C interface of initializing option
   !> @param option_Cptr: the structure containing option
   subroutine C_BPACK_Createoption(option_Cptr) bind(c, name="c_bpack_createoption")
      implicit none
      type(c_ptr) :: option_Cptr
      type(Hoption), pointer::option

      allocate (option)
      !>**** set default hodlr options
      call SetDefaultOptions(option)

      option_Cptr = c_loc(option)

   end subroutine C_BPACK_Createoption

!>**** C interface of copy option
   !> @param option_Cptr: the structure containing option
   !> @param option_Cptr1: the structure containing option
   subroutine C_BPACK_Copyoption(option_Cptr, option_Cptr1) bind(c, name="c_bpack_copyoption")
      implicit none
      type(c_ptr) :: option_Cptr, option_Cptr1
      type(Hoption), pointer::option, option1

      call c_f_pointer(option_Cptr, option)

      !>****copy hodlr options
      allocate (option1)
      call CopyOptions(option, option1)

      option_Cptr1 = c_loc(option1)

   end subroutine C_BPACK_Copyoption

!>**** C interface of printing option
   !> @param option_Cptr: the structure containing option
   !> @param ptree_Cptr: the structure containing process tree
   subroutine C_BPACK_Printoption(option_Cptr, ptree_Cptr) bind(c, name="c_bpack_printoption")
      implicit none
      type(c_ptr) :: option_Cptr
      type(c_ptr):: ptree_Cptr
      type(Hoption), pointer::option
      type(proctree), pointer::ptree

      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(ptree_Cptr, ptree)
      call PrintOptions(option, ptree)

   end subroutine C_BPACK_Printoption


!>**** C interface of getting one entry in option, always returning double
   !> @param option_Cptr: the structure containing option
   !> @param nam: name of the option
   !> @param val_d: value of the option
   subroutine C_BPACK_Getoption(option_Cptr, nam, val_d) bind(c, name="c_bpack_getoption")
      implicit none
      real(kind=8)::val_d
      character(kind=c_char, len=1) :: nam(*)
      character(kind=c_char) :: tmpc
      type(c_ptr) :: option_Cptr
      type(Hoption), pointer::option
      ! character::nam(:)
      ! type(c_ptr),value :: val_Cptr
      ! integer,pointer::val_i
      ! real(kind=8),pointer::val_d
      integer strlen
      character(len=:), allocatable :: str
      integer valid_opt

      valid_opt = 0
      strlen = 1
      tmpc=nam(strlen)
      do while (tmpc /= c_null_char)
         strlen = strlen + 1
         tmpc=nam(strlen)
      enddo
      strlen = strlen - 1
      allocate (character(len=strlen) :: str)
      str = transfer(nam(1:strlen), str)

      call c_f_pointer(option_Cptr, option)

      if (trim(str) == 'n_iter') then
         val_d = option%n_iter
         valid_opt = 1
      endif
      if (trim(str) == 'precon') then
         val_d = option%precon
         valid_opt = 1
      endif
      if (trim(str) == 'xyzsort') then
         val_d = option%xyzsort
         valid_opt = 1
      endif
      if (trim(str) == 'lnoBP') then
         val_d = option%lnoBP
         valid_opt = 1
      endif
      if (trim(str) == 'bp_cnt_lr') then
         val_d = option%bp_cnt_lr
         valid_opt = 1
      endif
      if (trim(str) == 'TwoLayerOnly') then
         val_d = option%TwoLayerOnly
         valid_opt = 1
      endif
      if (trim(str) == 'schulzorder') then
         val_d = option%schulzorder
         valid_opt = 1
      endif
      if (trim(str) == 'schulzhardstart') then
         val_d = option%schulzhardstart
         valid_opt = 1
      endif
      if (trim(str) == 'schulzsplitlevel') then
         val_d = option%schulzsplitlevel
         valid_opt = 1
      endif
      if (trim(str) == 'schulzsplitlevel') then
         val_d = option%schulzsplitlevel
         valid_opt = 1
      endif
      if (trim(str) == 'schulzlevel') then
         val_d = option%schulzlevel
         valid_opt = 1
      endif
      if (trim(str) == 'LRlevel') then
         val_d = option%LRlevel
         valid_opt = 1
      endif
      if (trim(str) == 'ErrFillFull') then
         val_d = option%ErrFillFull
         valid_opt = 1
      endif
      if (trim(str) == 'BACA_Batch') then
         val_d = option%BACA_Batch
         valid_opt = 1
      endif
      if (trim(str) == 'ErrSol') then
         val_d = option%ErrSol
         valid_opt = 1
      endif
      if (trim(str) == 'nogeo') then
         val_d = option%nogeo
         valid_opt = 1
      endif
      if (trim(str) == 'per_geo') then
         val_d = option%per_geo
         valid_opt = 1
      endif
      if (trim(str) == 'less_adapt') then
         val_d = option%less_adapt
         valid_opt = 1
      endif
      if (trim(str) == 'RecLR_leaf') then
         val_d = option%RecLR_leaf
         valid_opt = 1
      endif
      if (trim(str) == 'Nmin_leaf') then
         val_d = option%Nmin_leaf
         valid_opt = 1
      endif
      if (trim(str) == 'LR_BLK_NUM') then
         val_d = option%LR_BLK_NUM
         valid_opt = 1
      endif
      if (trim(str) == 'rank0') then
         val_d = option%rank0
         valid_opt = 1
      endif
      if (trim(str) == 'period1') then
         val_d = option%periods(1)
         valid_opt = 1
      endif
      if (trim(str) == 'period2') then
         val_d = option%periods(2)
         valid_opt = 1
      endif
      if (trim(str) == 'period3') then
         val_d = option%periods(3)
         valid_opt = 1
      endif
      if (trim(str) == 'itermax') then
         val_d = option%itermax
         valid_opt = 1
      endif
      if (trim(str) == 'powiter') then
         val_d = option%powiter
         valid_opt = 1
      endif
      if (trim(str) == 'ILU') then
         val_d = option%ILU
         valid_opt = 1
      endif
      if (trim(str) == 'Nbundle') then
         val_d = option%Nbundle
         valid_opt = 1
      endif
      if (trim(str) == 'format') then
         val_d = option%format
         valid_opt = 1
      endif
      if (trim(str) == 'verbosity') then
         val_d = option%verbosity
         valid_opt = 1
      endif
      if (trim(str) == 'rmax') then
         val_d = option%rmax
         valid_opt = 1
      endif
      if (trim(str) == 'forwardN15flag') then
         val_d = option%forwardN15flag
         valid_opt = 1
      endif
      if (trim(str) == 'pat_comp') then
         val_d = option%pat_comp
         valid_opt = 1
      endif
      if (trim(str) == 'elem_extract') then
         val_d = option%elem_extract
         valid_opt = 1
      endif
      if (trim(str) == 'knn') then
         val_d = option%knn
         valid_opt = 1
      endif
      if (trim(str) == 'use_zfp') then
         val_d = option%use_zfp
         valid_opt = 1
      endif
      if (trim(str) == 'cpp') then
         val_d = option%cpp
         valid_opt = 1
      endif
      if (trim(str) == 'tol_comp') then
         val_d = option%tol_comp
         valid_opt = 1
      endif
      if (trim(str) == 'tol_Rdetect') then
         val_d = option%tol_Rdetect
         valid_opt = 1
      endif
      if (trim(str) == 'tol_LS') then
         val_d = option%tol_LS
         valid_opt = 1
      endif
      if (trim(str) == 'jitter') then
         val_d = option%jitter
         valid_opt = 1
      endif
      if (trim(str) == 'tol_itersol') then
         val_d = option%tol_itersol
         valid_opt = 1
      endif
      if (trim(str) == 'tol_rand') then
         val_d = option%tol_rand
         valid_opt = 1
      endif
      if (trim(str) == 'touch_para') then
         val_d = option%touch_para
         valid_opt = 1
      endif
      if (trim(str) == 'rankrate') then
         val_d = option%rankrate
         valid_opt = 1
      endif
      if (trim(str) == 'near_para') then
         val_d = option%near_para
         valid_opt = 1
      endif
      if (trim(str) == 'knn_near_para') then
         val_d = option%knn_near_para
         valid_opt = 1
      endif
      if (trim(str) == 'sample_para') then
         val_d = option%sample_para
         valid_opt = 1
      endif
      if (trim(str) == 'sample_para_outer') then
         val_d = option%sample_para_outer
         valid_opt = 1
      endif

      if (valid_opt == 0) write (*, *) 'invalid BPACK option: '//trim(str)
      deallocate (str)

   end subroutine C_BPACK_Getoption


!>**** C interface of set one entry in option
   !> @param option_Cptr: the structure containing option
   !> @param nam: name of the option
   !> @param val_Cptr: value of the option
   subroutine C_BPACK_Setoption(option_Cptr, nam, val_Cptr) bind(c, name="c_bpack_setoption")
      implicit none
      type(c_ptr) :: option_Cptr
      character(kind=c_char, len=1) :: nam(*)
      type(Hoption), pointer::option
      ! character::nam(:)
      type(c_ptr), value :: val_Cptr
      integer, pointer::val_i
      real(kind=8), pointer::val_d
      integer strlen
      character(len=:), allocatable :: str
      integer valid_opt

      valid_opt = 0
      strlen = 1
      do while (nam(strlen) /= c_null_char)
         strlen = strlen + 1
      enddo
      strlen = strlen - 1
      allocate (character(len=strlen) :: str)
      str = transfer(nam(1:strlen), str)

      call c_f_pointer(option_Cptr, option)

!>**** integer parameters
      if (trim(str) == 'n_iter') then
         call c_f_pointer(val_Cptr, val_i)
         option%n_iter = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'precon') then
         call c_f_pointer(val_Cptr, val_i)
         option%precon = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'xyzsort') then
         call c_f_pointer(val_Cptr, val_i)
         option%xyzsort = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'lnoBP') then
         call c_f_pointer(val_Cptr, val_i)
         option%lnoBP = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'bp_cnt_lr') then
         call c_f_pointer(val_Cptr, val_i)
         option%bp_cnt_lr = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'TwoLayerOnly') then
         call c_f_pointer(val_Cptr, val_i)
         option%TwoLayerOnly = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'schulzorder') then
         call c_f_pointer(val_Cptr, val_i)
         option%schulzorder = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'schulzhardstart') then
         call c_f_pointer(val_Cptr, val_i)
         option%schulzhardstart = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'schulzsplitlevel') then
         call c_f_pointer(val_Cptr, val_i)
         option%schulzsplitlevel = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'schulzlevel') then
         call c_f_pointer(val_Cptr, val_i)
         option%schulzlevel = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'LRlevel') then
         call c_f_pointer(val_Cptr, val_i)
         option%LRlevel = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'ErrFillFull') then
         call c_f_pointer(val_Cptr, val_i)
         option%ErrFillFull = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'BACA_Batch') then
         call c_f_pointer(val_Cptr, val_i)
         option%BACA_Batch = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'ErrSol') then
         call c_f_pointer(val_Cptr, val_i)
         option%ErrSol = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'nogeo') then
         call c_f_pointer(val_Cptr, val_i)
         option%nogeo = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'per_geo') then
         call c_f_pointer(val_Cptr, val_i)
         option%per_geo = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'less_adapt') then
         call c_f_pointer(val_Cptr, val_i)
         option%less_adapt = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'RecLR_leaf') then
         call c_f_pointer(val_Cptr, val_i)
         option%RecLR_leaf = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'Nmin_leaf') then
         call c_f_pointer(val_Cptr, val_i)
         option%Nmin_leaf = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'LR_BLK_NUM') then
         call c_f_pointer(val_Cptr, val_i)
         option%LR_BLK_NUM = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'rank0') then
         call c_f_pointer(val_Cptr, val_i)
         option%rank0 = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'itermax') then
         call c_f_pointer(val_Cptr, val_i)
         option%itermax = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'powiter') then
         call c_f_pointer(val_Cptr, val_i)
         option%powiter = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'ILU') then
         call c_f_pointer(val_Cptr, val_i)
         option%ILU = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'Nbundle') then
         call c_f_pointer(val_Cptr, val_i)
         option%Nbundle = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'format') then
         call c_f_pointer(val_Cptr, val_i)
         option%format = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'verbosity') then
         call c_f_pointer(val_Cptr, val_i)
         option%verbosity = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'rmax') then
         call c_f_pointer(val_Cptr, val_i)
         option%rmax = val_i
         valid_opt = 1
      endif
      if (trim(str) == 'forwardN15flag') then
         call c_f_pointer(val_Cptr, val_i)
         option%forwardN15flag = val_i
         valid_opt = 1
      endif

      if (trim(str) == 'pat_comp') then
         call c_f_pointer(val_Cptr, val_i)
         option%pat_comp = val_i
         valid_opt = 1
      endif

      if (trim(str) == 'elem_extract') then
         call c_f_pointer(val_Cptr, val_i)
         option%elem_extract = val_i
         valid_opt = 1
      endif

      if (trim(str) == 'knn') then
         call c_f_pointer(val_Cptr, val_i)
         option%knn = val_i
         valid_opt = 1
      endif

      if (trim(str) == 'cpp') then
         call c_f_pointer(val_Cptr, val_i)
         option%cpp = val_i
         valid_opt = 1
      endif

      if (trim(str) == 'use_zfp') then
         call c_f_pointer(val_Cptr, val_i)
         option%use_zfp = val_i
         valid_opt = 1
      endif

      ! if (trim(str) == 'sample_heuristic') then
      !    call c_f_pointer(val_Cptr, val_i)
      !    option%sample_heuristic = val_i
      !    valid_opt = 1
      ! endif

!>**** double parameters
      if (trim(str) == 'tol_comp') then
         call c_f_pointer(val_Cptr, val_d)
         option%tol_comp = val_d
         valid_opt = 1
      endif
      if (trim(str) == 'tol_Rdetect') then
         call c_f_pointer(val_Cptr, val_d)
         option%tol_Rdetect = val_d
         valid_opt = 1
      endif
      if (trim(str) == 'tol_LS') then
         call c_f_pointer(val_Cptr, val_d)
         option%tol_LS = val_d
         valid_opt = 1
      endif
      if (trim(str) == 'jitter') then
         call c_f_pointer(val_Cptr, val_d)
         option%jitter = val_d
         valid_opt = 1
      endif
      if (trim(str) == 'tol_itersol') then
         call c_f_pointer(val_Cptr, val_d)
         option%tol_itersol = val_d
         valid_opt = 1
      endif
      if (trim(str) == 'tol_rand') then
         call c_f_pointer(val_Cptr, val_d)
         option%tol_rand = val_d
         valid_opt = 1
      endif
      if (trim(str) == 'touch_para') then
         call c_f_pointer(val_Cptr, val_d)
         option%touch_para = val_d
         valid_opt = 1
      endif
      if (trim(str) == 'rankrate') then
         call c_f_pointer(val_Cptr, val_d)
         option%rankrate = val_d
         valid_opt = 1
      endif
      if (trim(str) == 'near_para') then
         call c_f_pointer(val_Cptr, val_d)
         option%near_para = val_d
         valid_opt = 1
      endif
      if (trim(str) == 'knn_near_para') then
         call c_f_pointer(val_Cptr, val_d)
         option%knn_near_para = val_d
         valid_opt = 1
      endif

      if (trim(str) == 'sample_para') then
         call c_f_pointer(val_Cptr, val_d)
         option%sample_para = val_d
         valid_opt = 1
      endif

      if (trim(str) == 'sample_para_outer') then
         call c_f_pointer(val_Cptr, val_d)
         option%sample_para_outer = val_d
         valid_opt = 1
      endif

      if (trim(str) == 'period1') then
         call c_f_pointer(val_Cptr, val_d)
         option%periods(1) = val_d
         valid_opt = 1
      endif

      if (trim(str) == 'period2') then
         call c_f_pointer(val_Cptr, val_d)
         option%periods(2) = val_d
         valid_opt = 1
      endif
      if (trim(str) == 'period3') then
         call c_f_pointer(val_Cptr, val_d)
         option%periods(3) = val_d
         valid_opt = 1
      endif


      if (valid_opt == 0) write (*, *) 'invalid BPACK option: '//trim(str)

      deallocate (str)
      option_Cptr = c_loc(option)

   end subroutine C_BPACK_Setoption












!>**** C interface of matrix construction
   !> @param bmat_Cptr: the structure containing HODLR
   !> @param option_Cptr: the structure containing option
   !> @param stats_Cptr: the structure containing statistics
   !> @param msh_Cptr: the structure containing points and ordering information
   !> @param ker_Cptr: the structure containing kernel quantities
   !> @param ptree_Cptr: the structure containing process tree
   !> @param C_FuncZmn: the C_pointer to user-provided function to sample mn^th entry of the matrix
   !> @param C_FuncZmnBlock: the C_pointer to user-provided function to sample a list of intersections of entries of the matrix
   !> @param C_QuantApp: the C_pointer to user-defined quantities required to for entry evaluation,sampling,distance and compressibility test (in)
   subroutine C_BPACK_Construct_Element_Compute(bmat_Cptr, option_Cptr, stats_Cptr, msh_Cptr, ker_Cptr, ptree_Cptr, C_FuncZmn, C_FuncZmnBlock, C_QuantApp) bind(c, name="c_bpack_construct_element_compute")
      implicit none

      real(kind=8) para
      real(kind=8) tolerance, h, lam
      integer Primary_block, nn, mm, MyID_old, Maxlevel, give, need
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam
      integer level
      integer groupm
      type(c_ptr) :: bmat_Cptr
      type(c_ptr) :: option_Cptr
      type(c_ptr) :: stats_Cptr
      type(c_ptr) :: msh_Cptr
      type(c_ptr) :: ker_Cptr
      type(c_ptr) :: ptree_Cptr
      type(c_ptr), intent(in), target :: C_QuantApp
      type(c_funptr), intent(in), value, target :: C_FuncZmn
      type(c_funptr), intent(in), value, target :: C_FuncZmnBlock

      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(mesh), pointer::msh
      type(kernelquant), pointer::ker
      type(Bmatrix), pointer::bmat
      type(proctree), pointer::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2, x, y, z, r, theta, phi
      character(len=1024)  :: strings

      !>**** allocate HODLR solver structures
      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)
      call c_f_pointer(bmat_Cptr, bmat)
      call c_f_pointer(msh_Cptr, msh)
      call c_f_pointer(ker_Cptr, ker)

      !>**** register the user-defined function and type in ker
      ker%C_QuantApp => C_QuantApp
      ker%C_FuncZmn => C_FuncZmn
      ker%C_FuncZmnBlock => C_FuncZmnBlock

      t1 = MPI_Wtime()
      !>**** computation of the construction phase
      call BPACK_construction_Element(bmat, option, stats, msh, ker, ptree)
      ! call BPACK_CheckError(bmat,option,msh,ker,stats,ptree)
      t2 = MPI_Wtime()

      !>**** delete neighours in msh
      if(allocated(msh%nns))then
         call LogMemory(stats, -SIZEOF(msh%nns)/1024.0d3)
         deallocate(msh%nns)
      endif

      !>**** return the C address of hodlr structures to C caller
      bmat_Cptr = c_loc(bmat)
      option_Cptr = c_loc(option)
      stats_Cptr = c_loc(stats)
      msh_Cptr = c_loc(msh)
      ker_Cptr = c_loc(ker)
      ptree_Cptr = c_loc(ptree)

   end subroutine C_BPACK_Construct_Element_Compute

!>**** C interface of matrix construction via entry evaluation
   !> @param N: matrix size (in)
   !> @param Ndim: data set dimensionality (not used if nogeo=1)
   !> @param Locations: coordinates used for clustering (not used if nogeo=1)
   !> @param nns: nearest neighbours provided by user (referenced if nogeo=3 or 4)
   !> @param nlevel: the number of top levels that have been ordered (in)
   !> @param tree: the order tree provided by the caller, if incomplete, the init routine will make it complete (inout)
   !> @param Permutation: return the permutation vector new2old (indexed from 1) (out)
   !> @param N_loc: number of local row/column indices (out)
   !> @param bmat_Cptr: the structure containing HODLR (out)
   !> @param option_Cptr: the structure containing option (in)
   !> @param stats_Cptr: the structure containing statistics (inout)
   !> @param msh_Cptr: the structure containing points and ordering information (out)
   !> @param ker_Cptr: the structure containing kernel quantities (out)
   !> @param ptree_Cptr: the structure containing process tree (in)
   !> @param C_FuncDistmn: the C_pointer to user-provided function to compute distance between any row and column of the matrix
   !> @param C_FuncNearFar: the C_pointer to user-provided function to determine whether a block (in permuted order) is compressible or not
   !> @param C_QuantApp: the C_pointer to user-defined quantities required to for entry evaluation,sampling,distance and compressibility test (in)
   subroutine C_BPACK_Construct_Init(N, Ndim, Locations, nns, nlevel, tree, Permutation, N_loc, bmat_Cptr, option_Cptr, stats_Cptr, msh_Cptr, ker_Cptr, ptree_Cptr, C_FuncDistmn, C_FuncNearFar, C_QuantApp) bind(c, name="c_bpack_construct_init")
      implicit none
      integer N, Ndim
      real(kind=8) Locations(*)
      integer nns(*)
      real(kind=8) para
      real(kind=8) tolerance, h, lam
      integer Primary_block, nn, mm, MyID_old, Maxlevel, give, need
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam
      real(kind=8), parameter :: BPACK_cd = 299792458d0
      integer, allocatable:: groupmembers(:)
      integer nlevel, level
      integer Permutation(N), tree(*)
      integer N_loc
      ! type(matricesblock), pointer :: blocks_i
      integer groupm
      type(c_ptr) :: bmat_Cptr
      type(c_ptr) :: option_Cptr
      type(c_ptr) :: stats_Cptr
      type(c_ptr) :: msh_Cptr
      type(c_ptr) :: ker_Cptr
      type(c_ptr) :: ptree_Cptr
      type(c_ptr), intent(in), target :: C_QuantApp
      type(c_funptr), intent(in), value, target :: C_FuncDistmn
      type(c_funptr), intent(in), value, target :: C_FuncNearFar

      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(mesh), pointer::msh
      type(kernelquant), pointer::ker
      type(Bmatrix), pointer::bmat
      type(proctree), pointer::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2, x, y, z, r, theta, phi
      real(kind=8):: Memory = 0d0, error
      character(len=1024)  :: strings
      integer(kind=8) idx,kk,knn

      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)

      call assert(option%xyzsort /= TM_GRAM, 'gram distance based clustering is not supported in this interface')

      !>**** allocate HODLR solver structures
      allocate (bmat)
      ! allocate(option)
      ! allocate(stats)
      allocate (msh)
      allocate (ker)

      stats%Flop_Fill = 0
      stats%Time_Fill = 0
      stats%Time_Entry = 0
      stats%Time_Entry_Traverse = 0
      stats%Time_Entry_BF = 0
      stats%Time_Entry_Comm = 0

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'NUMBER_MPI=', ptree%nproc

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

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
         write (*, *) "HODLR_BUTTERFLY_SOLVER"
         write (*, *) "   "
      endif

      !>**** register the user-defined function and type in ker
      if (option%nogeo == 2) then
         ker%C_QuantApp => C_QuantApp
         ker%C_FuncDistmn => C_FuncDistmn
         ker%C_FuncNearFar => C_FuncNearFar
      endif

      msh%Nunk = N

      t1 = MPI_Wtime()

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "User-supplied kernel:"
      Maxlevel = nlevel
      allocate (msh%pretree(2**Maxlevel))

      msh%pretree(1:2**Maxlevel) = tree(1:2**Maxlevel)

      !>**** make 0-element node a 1-element node

      ! write(*,*)'before adjustment:',msh%pretree
      call assert(N>=2**Maxlevel,'The incomplete tree cannot be made complete. Try decreasing tree levels')
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

      !>**** the geometry points are provided by user
      if (option%nogeo == 0 .or. option%nogeo == 4) then
         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "User-supplied kernel requiring reorder:"
         Dimn = Ndim
         allocate (msh%xyz(Dimn, 1:msh%Nunk))
         ii = 0
         do edge = 1, msh%Nunk
            msh%xyz(1:Dimn, edge) = Locations(ii + 1:ii + Dimn)
            ii = ii + Dimn
         enddo
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
         allocate (msh%nns(msh%Nunk, option%knn))
         do ii = 1, msh%Nunk
         do kk = 1, option%knn
            knn = option%knn
            idx=kk + (msh%new2old(ii) - 1)*knn
            if (nns(idx) /= 0) then
               msh%nns(ii, kk) = msh%old2new(nns(idx))
            else
               msh%nns(ii, kk) = 0
            endif
         enddo
         enddo
      endif

      !>**** return the permutation vector
      N_loc = msh%idxe - msh%idxs + 1
      ! if (ptree%MyID == Main_ID) then
         do edge = 1, N
            Permutation(edge) = msh%new2old(edge)
         enddo
      ! endif

      !>**** return the C address of hodlr structures to C caller
      bmat_Cptr = c_loc(bmat)
      option_Cptr = c_loc(option)
      stats_Cptr = c_loc(stats)
      msh_Cptr = c_loc(msh)
      ker_Cptr = c_loc(ker)
      ptree_Cptr = c_loc(ptree)

   end subroutine C_BPACK_Construct_Init





!>**** C interface of multi-dimensional BF construction via entry evaluation
   !> @param Ns: size for each dimension (in)
   !> @param Nmax: maximum size among all dimensions (in)
   !> @param Ndim: data set dimensionality (in)
   !> @param Locations: coordinates (1D coordinates per dimension) used for clustering (in)
   !> @param Permutation: return the permutation vector (per dimension) new2old (indexed from 1) (out)
   !> @param N_loc: number of local row/column indices (per dimension) (out)
   !> @param bmat_Cptr: the structure containing the compressed operator (out)
   !> @param option_Cptr: the structure containing option (in)
   !> @param stats_Cptr: the structure containing statistics (inout)
   !> @param msh_Cptr: the structure containing points and ordering information per dimension (out)
   !> @param ker_Cptr: the structure containing kernel quantities (out)
   !> @param ptree_Cptr: the structure containing process tree (in)
   !> @param C_FuncNearFar: the C_pointer to user-provided function to determine whether a block (in permuted order) is compressible or not
   !> @param C_QuantApp: the C_pointer to user-defined quantities required to for entry evaluation,sampling,distance and compressibility test (in)
   subroutine C_BPACK_MD_Construct_Init(Ns, Nmax, Ndim, Locations, Permutation, N_loc, bmat_Cptr, option_Cptr, stats_Cptr, msh_Cptr, ker_Cptr, ptree_Cptr, C_FuncNearFar, C_QuantApp) bind(c, name="c_bpack_md_construct_init")
      implicit none
      integer Ndim,Nmax
      integer Ns(Ndim)
      real(kind=8) Locations(*)
      real(kind=8) para
      real(kind=8) tolerance, h, lam
      integer Primary_block, nn, mm, MyID_old, Maxlevel, give, need
      integer i, j, k, ii, edge, threads_num, nth, nmpi, ninc, acam
      real(kind=8), parameter :: BPACK_cd = 299792458d0
      integer, allocatable:: groupmembers(:)
      integer nlevel, level, dim_i
      integer Permutation(Nmax*Ndim)
      integer N_loc(Ndim)
      ! type(matricesblock), pointer :: blocks_i
      integer groupm
      type(c_ptr) :: bmat_Cptr
      type(c_ptr) :: option_Cptr
      type(c_ptr) :: stats_Cptr
      type(c_ptr) :: msh_Cptr
      type(c_ptr) :: ker_Cptr
      type(c_ptr) :: ptree_Cptr
      type(c_ptr), intent(in), target :: C_QuantApp
      type(c_funptr), intent(in), value, target :: C_FuncNearFar

      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(mesh), pointer::msh(:)
      type(kernelquant), pointer::ker
      type(Bmatrix), pointer::bmat
      type(proctree), pointer::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2, x, y, z, r, theta, phi
      real(kind=8):: Memory = 0d0, error
      character(len=1024)  :: strings
      integer(kind=8) idx,kk,knn

      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)

      call assert(option%xyzsort /= TM_GRAM, 'gram distance based clustering is not supported in this interface')

      !>**** allocate HODLR solver structures
      allocate (bmat)
      ! allocate(option)
      ! allocate(stats)
      allocate (msh(Ndim))
      allocate (ker)

      stats%Flop_Fill = 0
      stats%Time_Fill = 0
      stats%Time_Entry = 0
      stats%Time_Entry_Traverse = 0
      stats%Time_Entry_BF = 0
      stats%Time_Entry_Comm = 0

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'NUMBER_MPI=', ptree%nproc

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

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
         write (*, *) "HODLR_BUTTERFLY_SOLVER"
         write (*, *) "   "
      endif

      !>**** register the user-defined function and type in ker
      if (option%nogeo == 2) then
         write(*,*)"option%nogeo not yet implemented for BPACK_MD_Construct_Init"
         ker%C_QuantApp => C_QuantApp
         ker%C_FuncNearFar => C_FuncNearFar
      endif

      t1 = MPI_Wtime()

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "User-supplied kernel:"
      Maxlevel = nlevel



      !>**** the geometry points are provided by user
      if (option%nogeo == 0 .or. option%nogeo == 4) then
         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "User-supplied kernel requiring reorder:"
         do dim_i=1,Ndim
            msh(dim_i)%Nunk = Ns(dim_i)
            allocate (msh(dim_i)%xyz(1, 1:msh(dim_i)%Nunk))
            ii = dim_i-1
            do edge = 1, msh(dim_i)%Nunk
               msh(dim_i)%xyz(1, edge) = Locations(ii + 1)
               ii = ii + Ndim
            enddo
         enddo
      endif

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "    "
      t2 = MPI_Wtime()

      t1 = MPI_Wtime()
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Hierarchical format......"
      call Cluster_partition_MD(Ndim, bmat, option, msh, ker, stats, ptree)
      call BPACK_structuring_MD(Ndim, bmat, option, msh, ker, ptree, stats)
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Hierarchical format finished"
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "    "
      t2 = MPI_Wtime()

      !>**** return the permutation vector
      Permutation=-1
      do dim_i=1,Ndim
         N_loc(dim_i) = msh(dim_i)%idxe - msh(dim_i)%idxs + 1
         ! if (ptree%MyID == Main_ID) then
            do edge = 1, Ns(dim_i)
               Permutation(edge+(dim_i-1)*Nmax) = msh(dim_i)%new2old(edge)
            enddo
         ! endif
      enddo

      !>**** return the C address of hodlr structures to C caller
      bmat_Cptr = c_loc(bmat)
      option_Cptr = c_loc(option)
      stats_Cptr = c_loc(stats)
      msh_Cptr = c_loc(msh)
      ker_Cptr = c_loc(ker)
      ptree_Cptr = c_loc(ptree)

   end subroutine C_BPACK_MD_Construct_Init



!>**** C interface of converting from new,local index to old, global index, the indexs start from 1
   !> @param newidx_loc: new, local index, from 1 to Nloc
   !> @param oldidx: old, global index, from 1 to N (out)
   !> @param msh_Cptr: the structure containing points and ordering information
   subroutine C_BPACK_New2Old(msh_Cptr, newidx_loc, oldidx) bind(c, name="c_bpack_new2old")
      implicit none
      integer newidx_loc,oldidx
      type(c_ptr) :: msh_Cptr
      type(mesh), pointer::msh

      call c_f_pointer(msh_Cptr, msh)
      oldidx = msh%new2old(msh%idxs + newidx_loc - 1)

   end subroutine C_BPACK_New2Old




!>**** C interface of converting from new,local index to old, global index, the indexs start from 1
   !> @param newidx_loc: new, local index, from 1 to Nloc (in)
   !> @param oldidx: old, global index, from 1 to N (out)
   !> @param mshr_Cptr: the structure containing points and ordering information for the row dimension (in)
   subroutine C_BF_New2Old_Row(mshr_Cptr, newidx_loc, oldidx) bind(c, name="c_bf_new2old_row")
      implicit none
      integer newidx_loc,oldidx
      type(c_ptr) :: msh_Cptr
      type(c_ptr) :: mshr_Cptr
      type(mesh), pointer::msh,mshr

      call c_f_pointer(mshr_Cptr, mshr)
      oldidx = mshr%new2old(mshr%idxs + newidx_loc - 1)
   end subroutine C_BF_New2Old_Row

!>**** C interface of converting from new,local index to old, global index, the indexs start from 1
   !> @param newidx_loc: new, local index, from 1 to Nloc (in)
   !> @param oldidx: old, global index, from 1 to N (out)
   !> @param mshr_Cptr: the structure containing points and ordering information for the column dimension (in)
   subroutine C_BF_New2Old_Col(mshc_Cptr, newidx_loc, oldidx) bind(c, name="c_bf_new2old_col")
      implicit none
      integer newidx_loc,oldidx
      type(c_ptr) :: msh_Cptr
      type(c_ptr) :: mshc_Cptr
      type(mesh), pointer::msh,mshc

      call c_f_pointer(mshc_Cptr, mshc)
      oldidx = mshc%new2old(mshc%idxs + newidx_loc - 1)
   end subroutine C_BF_New2Old_Col

!>**** C interface of matrix construction via entry evaluation and using it for gram distance
   !> @param N: matrix size (in)
   !> @param Ndim: data set dimensionality (not used if nogeo=1)
   !> @param Locations: coordinates used for clustering (not used if nogeo=1)
   !> @param nns: nearest neighbours provided by user (referenced if nogeo=3 or 4)
   !> @param nlevel: the number of top levels that have been ordered (in)
   !> @param tree: the order tree provided by the caller, if incomplete, the init routine will make it complete (inout)
   !> @param Permutation: return the permutation vector new2old (indexed from 1) (out)
   !> @param N_loc: number of local row/column indices (out)
   !> @param bmat_Cptr: the structure containing HODLR (out)
   !> @param option_Cptr: the structure containing option (in)
   !> @param stats_Cptr: the structure containing statistics (inout)
   !> @param msh_Cptr: the structure containing points and ordering information (out)
   !> @param ker_Cptr: the structure containing kernel quantities (out)
   !> @param ptree_Cptr: the structure containing process tree (in)
   !> @param C_FuncZmn: the C_pointer to user-provided function to sample mn^th entry of the matrix
   !> @param C_FuncZmnBlock: the C_pointer to user-provided function to sample a list of intersections of entries of the matrix
   !> @param C_QuantApp: the C_pointer to user-defined quantities required to for entry evaluation,sampling,distance and compressibility test (in)
   subroutine C_BPACK_Construct_Init_Gram(N, Ndim, Locations, nns, nlevel, tree, Permutation, N_loc, bmat_Cptr, option_Cptr, stats_Cptr, msh_Cptr, ker_Cptr, ptree_Cptr, C_FuncZmn, C_FuncZmnBlock, C_QuantApp) bind(c, name="c_bpack_construct_init_gram")
      implicit none
      integer N, Ndim
      real(kind=8) Locations(*)
      integer nns(*)
      real(kind=8) para
      real(kind=8) tolerance, h, lam
      integer Primary_block, nn, mm, MyID_old, Maxlevel, give, need
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam
      real(kind=8), parameter :: BPACK_cd = 299792458d0
      integer, allocatable:: groupmembers(:)
      integer nlevel, level
      integer Permutation(N), tree(*)
      integer N_loc
      ! type(matricesblock), pointer :: blocks_i
      integer groupm
      type(c_ptr) :: bmat_Cptr
      type(c_ptr) :: option_Cptr
      type(c_ptr) :: stats_Cptr
      type(c_ptr) :: msh_Cptr
      type(c_ptr) :: ker_Cptr
      type(c_ptr) :: ptree_Cptr
      type(c_ptr), intent(in), target :: C_QuantApp
      type(c_funptr), intent(in), value, target :: C_FuncZmn
      type(c_funptr), intent(in), value, target :: C_FuncZmnBlock

      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(mesh), pointer::msh
      type(kernelquant), pointer::ker
      type(Bmatrix), pointer::bmat
      type(proctree), pointer::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2, x, y, z, r, theta, phi
      real(kind=8):: Memory = 0d0, error
      character(len=1024)  :: strings
      integer(kind=8)::idx,kk,knn

      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)

      call assert(option%xyzsort == TM_GRAM, 'only gram distance based clustering is supported in this interface')

      !>**** allocate HODLR solver structures
      allocate (bmat)
      ! allocate(option)
      ! allocate(stats)
      allocate (msh)
      allocate (ker)

      stats%Flop_Fill = 0
      stats%Time_Fill = 0
      stats%Time_Entry = 0
      stats%Time_Entry_Traverse = 0
      stats%Time_Entry_BF = 0
      stats%Time_Entry_Comm = 0

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'NUMBER_MPI=', ptree%nproc

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

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) then
         write (*, *) "HODLR_BUTTERFLY_SOLVER"
         write (*, *) "   "
      endif

      !>**** register the user-defined function and type in ker
      ker%C_QuantApp => C_QuantApp
      ker%C_FuncZmn => C_FuncZmn
      ker%C_FuncZmnBlock => C_FuncZmnBlock

      msh%Nunk = N

      t1 = MPI_Wtime()

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "User-supplied kernel:"
      Maxlevel = nlevel
      allocate (msh%pretree(2**Maxlevel))

      msh%pretree(1:2**Maxlevel) = tree(1:2**Maxlevel)

      !>**** make 0-element node a 1-element node

      ! write(*,*)'before adjustment:',msh%pretree
      call assert(N>=2**Maxlevel,'The incomplete tree cannot be made complete. Try decreasing tree levels')
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

      !>**** the geometry points are provided by user
      if (option%nogeo == 0 .or. option%nogeo == 4) then
         if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "User-supplied kernel requiring reorder:"
         Dimn = Ndim
         allocate (msh%xyz(Dimn, 1:msh%Nunk))
         ii = 0
         do edge = 1, msh%Nunk
            msh%xyz(1:Dimn, edge) = Locations(ii + 1:ii + Dimn)
            ii = ii + Dimn
         enddo
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
         allocate (msh%nns(msh%Nunk, option%knn))
         do ii = 1, msh%Nunk
         do kk = 1, option%knn
            knn = option%knn
            idx=kk + (msh%new2old(ii) - 1)*knn

            if (nns(idx) /= 0) then
               msh%nns(ii, kk) = msh%old2new(nns(idx))
            else
               msh%nns(ii, kk) = 0
            endif
         enddo
         enddo
      endif

      !>**** return the permutation vector
      N_loc = msh%idxe - msh%idxs + 1
      if (ptree%MyID == Main_ID) then
         do edge = 1, N
            Permutation(edge) = msh%new2old(edge)
         enddo
      endif

      !>**** return the C address of hodlr structures to C caller
      bmat_Cptr = c_loc(bmat)
      option_Cptr = c_loc(option)
      stats_Cptr = c_loc(stats)
      msh_Cptr = c_loc(msh)
      ker_Cptr = c_loc(ker)
      ptree_Cptr = c_loc(ptree)

   end subroutine C_BPACK_Construct_Init_Gram

!>**** C interface of matrix construction via blackbox matvec
   !> @param bmat_Cptr: the structure containing HODLR (inout)
   !> @param option_Cptr: the structure containing option (in)
   !> @param stats_Cptr: the structure containing statistics (inout)
   !> @param msh_Cptr: the structure containing points and ordering information (in)
   !> @param ker_Cptr: the structure containing kernel quantities (inout)
   !> @param ptree_Cptr: the structure containing process tree (in)
   !> @param C_FuncHMatVec: the C_pointer to user-provided function to multiply A and A* with vectors (in)
   !> @param C_QuantApp: the C_pointer to user-defined quantities required to for entry evaluation,sampling,distance and compressibility test (in)
   subroutine C_BPACK_Construct_Matvec_Compute(bmat_Cptr, option_Cptr, stats_Cptr, msh_Cptr, ker_Cptr, ptree_Cptr, C_FuncHMatVec, C_QuantApp) bind(c, name="c_bpack_construct_matvec_compute")
      implicit none

      real(kind=8) para
      real(kind=8) tolerance, h, lam
      integer Primary_block, nn, mm, MyID_old, Maxlevel, give, need
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam
      real(kind=8), parameter :: BPACK_cd = 299792458d0
      integer, allocatable:: groupmembers(:)
      integer level
      ! type(matricesblock), pointer :: blocks_i
      integer groupm
      type(c_ptr) :: bmat_Cptr
      type(c_ptr) :: option_Cptr
      type(c_ptr) :: stats_Cptr
      type(c_ptr) :: msh_Cptr
      type(c_ptr) :: ker_Cptr
      type(c_ptr) :: ptree_Cptr
      type(c_ptr), intent(in), target :: C_QuantApp
      type(c_funptr), intent(in), value, target :: C_FuncHMatVec

      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(mesh), pointer::msh
      type(kernelquant), pointer::ker
      type(Bmatrix), pointer::bmat
      type(proctree), pointer::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2, x, y, z, r, theta, phi
      real(kind=8):: Memory = 0d0, error
      character(len=1024)  :: strings

      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)
      call c_f_pointer(bmat_Cptr, bmat)
      call c_f_pointer(msh_Cptr, msh)
      call c_f_pointer(ker_Cptr, ker)

      !>**** register the user-defined function and type in ker
      ker%C_QuantApp => C_QuantApp
      ker%C_FuncHMatVec => C_FuncHMatVec

      !>**** computation of the construction phase
      option%less_adapt=0
      call BPACK_construction_Matvec(bmat, matvec_user_C, Memory, error, option, stats, ker, ptree, msh)
      option%less_adapt=1
      !>**** return the C address of hodlr structures to C caller
      bmat_Cptr = c_loc(bmat)
      option_Cptr = c_loc(option)
      stats_Cptr = c_loc(stats)
      msh_Cptr = c_loc(msh)
      ker_Cptr = c_loc(ker)
      ptree_Cptr = c_loc(ptree)

   end subroutine C_BPACK_Construct_Matvec_Compute

!>**** C interface of BF construction via blackbox matvec or entry extraction
   !> @param M,N: matrix size (in)
   !> @param M_loc,N_loc: number of local row/column indices (out)
   !> @param nnsr: (DIM knn*M) nearest neighbours(indexed from 1 to N) for each row (from 1 to M) provided by user (referenced if nogeo=3 or 4)
   !> @param nnsc: (DIM knn*N) nearest neighbours(indexed from 1 to M) for each column (from 1 to N) provided by user (referenced if nogeo=3 or 4)
   !> @param bf_Cptr: the structure containing the block (out)
   !> @param option_Cptr: the structure containing option (in)
   !> @param stats_Cptr: the structure containing statistics (inout)
   !> @param msh_Cptr: the structure containing points and ordering information combined from mshr_Cptr and mshc_Cptr (out)
   !> @param mshr_Cptr: the structure containing points and ordering information for the row dimension (in)
   !> @param mshc_Cptr: the structure containing points and ordering information for the column dimension (in)
   !> @param ker_Cptr: the structure containing kernel quantities (out)
   !> @param ptree_Cptr: the structure containing process tree (in)
   !> @param C_FuncDistmn: the C_pointer to user-provided function to compute distance between any row and column of the matrix
   !> @param C_FuncNearFar: the C_pointer to user-provided function to determine whether a block (in permuted order) is compressible or not
   !> @param C_QuantApp: the C_pointer to user-defined quantities required to for entry evaluation,sampling,distance and compressibility test (in)
   subroutine C_BF_Construct_Init(M, N, M_loc, N_loc, nnsr, nnsc, mshr_Cptr, mshc_Cptr, bf_Cptr, option_Cptr, stats_Cptr, msh_Cptr, ker_Cptr, ptree_Cptr, C_FuncDistmn, C_FuncNearFar, C_QuantApp) bind(c, name="c_bf_construct_init")
      implicit none
      integer M, N

      integer Maxlevel
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam
      integer, allocatable:: groupmembers(:)
      integer level
      integer M_loc, N_loc
      integer nnsr(*), nnsc(*)
      type(c_ptr) :: bf_Cptr
      type(c_ptr) :: option_Cptr
      type(c_ptr) :: stats_Cptr
      type(c_ptr) :: msh_Cptr, mshr_Cptr, mshc_Cptr
      type(c_ptr) :: ker_Cptr
      type(c_ptr) :: ptree_Cptr

      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(mesh), pointer::msh, mshr, mshc
      type(kernelquant), pointer::ker
      type(matrixblock), pointer::blocks
      type(proctree), pointer::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2
      character(len=1024)  :: strings
      integer Maxgroup_rc
      integer(kind=8)::idx,kk,knn

      type(c_ptr), intent(in), target :: C_QuantApp
      type(c_funptr), intent(in), value, target :: C_FuncDistmn
      type(c_funptr), intent(in), value, target :: C_FuncNearFar

      !>**** allocate HODLR solver structures
      allocate (blocks)
      allocate (msh)
      allocate (ker)

      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)
      call c_f_pointer(mshr_Cptr, mshr)
      call c_f_pointer(mshc_Cptr, mshc)

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'NUMBER_MPI=', ptree%nproc

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

      ! if(ptree%MyID==Main_ID)then
      ! write(*,*) "HODLR_BUTTERFLY_SOLVER"
      ! write(*,*) "   "
      ! endif

      !>**** register the user-defined function and type in ker
      if (option%nogeo == 2) then
         ker%C_QuantApp => C_QuantApp
         ker%C_FuncDistmn => C_FuncDistmn
         ker%C_FuncNearFar => C_FuncNearFar
      endif

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
            idx=kk + (ii - 1)*knn
            if (nnsr(idx) /= 0) then
               msh%nns(ii, kk) = nnsr(idx) + M
            else
               msh%nns(ii, kk) = 0
            endif
         enddo
         enddo
         do ii = 1, N
         do kk = 1, option%knn
            knn = option%knn
            idx=kk + (ii - 1)*knn
            msh%nns(ii + M, kk) = nnsc(idx)
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

      !>**** return the C address of hodlr structures to C caller
      bf_Cptr = c_loc(blocks)
      option_Cptr = c_loc(option)
      stats_Cptr = c_loc(stats)
      msh_Cptr = c_loc(msh)
      ker_Cptr = c_loc(ker)
      ptree_Cptr = c_loc(ptree)

   end subroutine C_BF_Construct_Init

!>**** C interface of BF construction via blackbox matvec
   !> @param bf_Cptr: the structure containing the block (inout)
   !> @param option_Cptr: the structure containing option (in)
   !> @param stats_Cptr: the structure containing statistics (inout)
   !> @param msh_Cptr: the structure containing points and ordering information (in)
   !> @param ker_Cptr: the structure containing kernel quantities (inout)
   !> @param ptree_Cptr: the structure containing process tree (in)
   !> @param C_FuncBMatVec: the C_pointer to user-provided function to multiply A and A* with vectors (in)
   !> @param C_QuantApp: the C_pointer to user-defined quantities required to for entry evaluation,sampling,distance and compressibility test (in)
   subroutine C_BF_Construct_Matvec_Compute(bf_Cptr, option_Cptr, stats_Cptr, msh_Cptr, ker_Cptr, ptree_Cptr, C_FuncBMatVec, C_QuantApp) bind(c, name="c_bf_construct_matvec_compute")
      implicit none

      integer Maxlevel
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam
      integer, allocatable:: groupmembers(:)
      integer level
      type(c_ptr) :: bf_Cptr
      type(c_ptr) :: option_Cptr
      type(c_ptr) :: stats_Cptr
      type(c_ptr) :: msh_Cptr
      type(c_ptr) :: ker_Cptr
      type(c_ptr) :: ptree_Cptr
      type(c_ptr), intent(in), target :: C_QuantApp
      type(c_funptr), intent(in), value, target :: C_FuncBMatVec

      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(mesh), pointer::msh
      type(kernelquant), pointer::ker
      type(matrixblock), pointer::blocks
      type(proctree), pointer::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2, error
      integer ierr

      !>**** allocate HODLR solver structures

      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)
      call c_f_pointer(msh_Cptr, msh)
      call c_f_pointer(bf_Cptr, blocks)
      call c_f_pointer(ker_Cptr, ker)
      if (allocated(msh%xyz)) deallocate (msh%xyz)

      if (.not. allocated(stats%rankmax_of_level)) allocate (stats%rankmax_of_level(0:0))
      stats%rankmax_of_level(0) = 0
      if (.not. allocated(stats%rankmax_of_level_global)) allocate (stats%rankmax_of_level_global(0:0))
      stats%rankmax_of_level_global(0) = 0

      ! !>**** register the user-defined function and type in ker
      ker%C_QuantApp => C_QuantApp
      ker%C_FuncBMatVec => C_FuncBMatVec

      t1 = MPI_Wtime()
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) " "
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "FastMATVEC-based BF construction......"

      call BF_randomized(blocks%pgno, blocks%level_butterfly, option%rank0, option%rankrate, blocks, ker, Bmatvec_user_C, error, 'CMatVec', option, stats, ptree, msh)

      call BF_ComputeMemory(blocks, stats%Mem_Comp_for)

      t2 = MPI_Wtime()

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "FastMATVEC-based BF construction finished in", t2-t1, 'Seconds with', stats%Mem_Comp_for,'MB Memory'




      stats%rankmax_of_level(0) = blocks%rankmax
      stats%rankmax_of_level_global(0) = stats%rankmax_of_level(0)
      stats%Mem_Fill = stats%Mem_Comp_for + stats%Mem_Direct_for
      call LogMemory(stats, stats%Mem_Fill)
      stats%Time_Fill = stats%Time_Fill + t2 - t1
      stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp
      !>**** return the C address of hodlr structures to C caller
      bf_Cptr = c_loc(blocks)
      option_Cptr = c_loc(option)
      stats_Cptr = c_loc(stats)
      msh_Cptr = c_loc(msh)
      ker_Cptr = c_loc(ker)
      ptree_Cptr = c_loc(ptree)


      ! if (ptree%MyID == Main_ID ) then
      !    write (*, *) 'time_tmp', time_tmp, 'randomized_bf time,', t2-t1, 'stats%Time_random,', stats%Time_random, 'mem', stats%Mem_Comp_for
      ! endif


   end subroutine C_BF_Construct_Matvec_Compute

!>**** C interface of BF construction via entry extraction
   !> @param bf_Cptr: the structure containing the block (inout)
   !> @param option_Cptr: the structure containing option (in)
   !> @param stats_Cptr: the structure containing statistics (inout)
   !> @param msh_Cptr: the structure containing points and ordering information (in)
   !> @param ker_Cptr: the structure containing kernel quantities (inout)
   !> @param ptree_Cptr: the structure containing process tree (in)
   !> @param C_FuncZmn: the C_pointer to user-provided function to sample mn^th entry of a block (in)
   !> @param C_FuncZmnBlock: the C_pointer to user-provided function to extract a list of intersections from a block (in)
   !> @param C_QuantApp: the C_pointer to user-defined quantities required to for entry evaluation,sampling,distance and compressibility test (in)
   subroutine C_BF_Construct_Element_Compute(bf_Cptr, option_Cptr, stats_Cptr, msh_Cptr, ker_Cptr, ptree_Cptr, C_FuncZmn, C_FuncZmnBlock, C_QuantApp) bind(c, name="c_bf_construct_element_compute")
      implicit none

      integer Maxlevel
      integer i, j, k, ii, edge, threads_num, nth, Dimn, nmpi, ninc, acam
      integer, allocatable:: groupmembers(:)
      integer level
      type(c_ptr) :: bf_Cptr
      type(c_ptr) :: option_Cptr
      type(c_ptr) :: stats_Cptr
      type(c_ptr) :: msh_Cptr
      type(c_ptr) :: ker_Cptr
      type(c_ptr) :: ptree_Cptr
      type(c_ptr), intent(in), target :: C_QuantApp
      type(c_funptr), intent(in), value, target :: C_FuncZmn
      type(c_funptr), intent(in), value, target :: C_FuncZmnBlock

      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(mesh), pointer::msh
      type(kernelquant), pointer::ker
      type(matrixblock), pointer::blocks
      type(proctree), pointer::ptree
      integer seed_myid(50)
      integer times(8)
      real(kind=8) t1, t2, error, Memory, tol_comp_tmp
      integer ierr,pp,knn_tmp
      integer:: boundary_map(1,1)
      integer groupm_start, Nboundall,Ninadmissible

      !>**** allocate HODLR solver structures

      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)
      call c_f_pointer(msh_Cptr, msh)
      call c_f_pointer(bf_Cptr, blocks)
      call c_f_pointer(ker_Cptr, ker)
      if (allocated(msh%xyz)) deallocate (msh%xyz)
      ! !>**** register the user-defined function and type in ker
      ker%C_QuantApp => C_QuantApp
      ker%C_FuncZmnBlock => C_FuncZmnBlock
      ker%C_FuncZmn => C_FuncZmn

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
         call BF_checkError(blocks, option, msh, ker, stats, ptree, 0, -1, error)
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

      if (option%verbosity >= 0) call BF_checkError(blocks, option, msh, ker, stats, ptree, 0, option%verbosity)
      call BF_ComputeMemory(blocks, stats%Mem_Comp_for)

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
      !>**** return the C address of hodlr structures to C caller
      bf_Cptr = c_loc(blocks)
      option_Cptr = c_loc(option)
      stats_Cptr = c_loc(stats)
      msh_Cptr = c_loc(msh)
      ker_Cptr = c_loc(ker)
      ptree_Cptr = c_loc(ptree)

   end subroutine C_BF_Construct_Element_Compute

!>**** C interface of HODLR factorization
   !> @param bmat_Cptr: the structure containing HODLR
   !> @param option_Cptr: the structure containing option
   !> @param stats_Cptr: the structure containing statistics
   !> @param ptree_Cptr: the structure containing process tree
   !> @param msh_Cptr: the structure containing points and ordering information (in)
   subroutine C_BPACK_Factor(bmat_Cptr, option_Cptr, stats_Cptr, ptree_Cptr, msh_Cptr) bind(c, name="c_bpack_factor")
      implicit none

      type(c_ptr), intent(inout) :: bmat_Cptr
      type(c_ptr), intent(in) :: ptree_Cptr
      type(c_ptr), intent(in) :: msh_Cptr
      type(c_ptr) :: option_Cptr
      type(c_ptr), intent(inout) :: stats_Cptr

      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(Bmatrix), pointer::bmat
      type(proctree), pointer::ptree
      type(mesh), pointer::msh

      ! real(kind=8):: tol_fact

      call c_f_pointer(bmat_Cptr, bmat)
      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)
      call c_f_pointer(msh_Cptr, msh)

      stats%Flop_Factor = 0
      stats%Time_Factor = 0

      call BPACK_Factorization(bmat, option, stats, ptree, msh)

      ! return the C address of hodlr structures to C caller
      bmat_Cptr = c_loc(bmat)

   end subroutine C_BPACK_Factor

!>**** C interface of HODLR solve
   !> @param x: local solution vector
   !> @param b: local RHS
   !> @param Nloc: size of local RHS
   !> @param Nrhs: number of RHSs
   !> @param bmat_Cptr: the structure containing HODLR
   !> @param option_Cptr: the structure containing option
   !> @param stats_Cptr: the structure containing statistics
   !> @param ptree_Cptr: the structure containing process tree
   subroutine C_BPACK_Solve(x, b, Nloc, Nrhs, bmat_Cptr, option_Cptr, stats_Cptr, ptree_Cptr) bind(c, name="c_bpack_solve")
      implicit none

      integer Nloc, Nrhs
      DT::x(Nloc, Nrhs), b(Nloc, Nrhs)

      type(c_ptr), intent(in) :: bmat_Cptr
      type(c_ptr), intent(in) :: ptree_Cptr
      type(c_ptr), intent(in) :: option_Cptr
      type(c_ptr), intent(in) :: stats_Cptr

      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(Bmatrix), pointer::bmat

      type(proctree), pointer::ptree

      call c_f_pointer(bmat_Cptr, bmat)

      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)

      stats%Flop_Sol = 0
      stats%Time_Sol = 0

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Solve ......"

      if (option%ErrSol == 1) then
         call BPACK_Test_Solve_error(bmat, Nloc, option, ptree, stats)
      endif

      call BPACK_Solution(bmat, x, b, Nloc, Nrhs, option, ptree, stats)

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "Solve finished"
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "    "

   end subroutine C_BPACK_Solve


!>**** C interface of a blackbox tfqmr without preconditioner, or assuming preconditioner is applied in the blackbox matvec
   !> @param x: local solution vector
   !> @param b: local RHS
   !> @param Nloc: size of local RHS
   !> @param Nrhs: number of RHSs
   !> @param option_Cptr: the structure containing option
   !> @param stats_Cptr: the structure containing statistics
   !> @param ptree_Cptr: the structure containing process tree
   !> @param C_QuantApp: the C_pointer to user-defined quantities required to for entry evaluation,sampling,distance and compressibility test (in)
   !> @param ker_Cptr: the structure containing kernel quantities
   !> @param C_FuncHMatVec: the C_pointer to user-provided function to multiply A and A* with vectors (in)
   subroutine C_BPACK_TFQMR_Noprecon(x, b, Nloc, Nrhs, option_Cptr, stats_Cptr, ptree_Cptr, ker_Cptr, C_FuncHMatVec, C_QuantApp) bind(c, name="c_bpack_tfqmr_noprecon")
      implicit none

      integer Nloc, Nrhs
      DT::x(Nloc, Nrhs), b(Nloc, Nrhs), r0_initial(Nloc,1)

      type(c_ptr), intent(in) :: ptree_Cptr
      type(c_ptr), intent(in) :: option_Cptr
      type(c_ptr), intent(in) :: stats_Cptr
      type(c_ptr) :: ker_Cptr
      type(c_funptr), intent(in), value, target :: C_FuncHMatVec
      type(c_ptr), intent(in), target :: C_QuantApp

      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(proctree), pointer::ptree
      type(kernelquant), pointer::ker

      real(kind=8):: rel_error, error, n1, n2
      integer ii, iter

      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)
      allocate (ker)

      !>**** register the user-defined function and type in ker
      ker%C_QuantApp => C_QuantApp
      ker%C_FuncHMatVec => C_FuncHMatVec

      stats%Flop_Sol = 0
      stats%Time_Sol = 0

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "TFQMR Solve ......"

      n1 = MPI_Wtime()
      do ii = 1, Nloc
         call random_dp_number(r0_initial(ii,1))
      end do

      do ii = 1, Nrhs
         iter = 0
         rel_error = option%tol_itersol
         call BPACK_Ztfqmr_usermatvec_noprecon(option%n_iter, Nloc, b(:, ii), x(:, ii), rel_error, iter, r0_initial, matvec_user_C, ptree, option, stats, ker)
      end do
      n2 = MPI_Wtime()
      stats%Time_Sol = stats%Time_Sol + n2 - n1

      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "TFQMR Solve finished"
      if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "    "

      ker_Cptr = c_loc(ker)

   end subroutine C_BPACK_TFQMR_Noprecon




!>**** C interface of butterfly-vector multiplication
   !> @param xin: input vector
   !> @param Ninloc: size of local input vectors
   !> @param xout: output vector
   !> @param Noutloc: size of local output vectors
   !> @param Ncol: number of vectors
   !> @param bf_for_Cptr: the structure containing butterfly
   !> @param option_Cptr: the structure containing options
   !> @param stats_Cptr: the structure containing statistics
   !> @param ptree_Cptr: the structure containing process tree
   !> @param trans: 'N', 'C' or 'T'
   subroutine C_BF_Mult(trans, xin, xout, Ninloc, Noutloc, Ncol, bf_for_Cptr, option_Cptr, stats_Cptr, ptree_Cptr) bind(c, name="c_bf_mult")
      implicit none
      real(kind=8) t1, t2
      integer Ninloc, Noutloc, Ncol
      DT::xin(Ninloc, Ncol), xout(Noutloc, Ncol)

      character(kind=c_char, len=1) :: trans(*)
      type(c_ptr), intent(in) :: bf_for_Cptr
      type(c_ptr), intent(in) :: ptree_Cptr
      type(c_ptr), intent(in) :: stats_Cptr
      type(c_ptr), intent(in) :: option_Cptr

      integer strlen
      character(len=:), allocatable :: str
      type(Hstat), pointer::stats
      type(Hoption), pointer::option
      type(matrixblock), pointer::blocks

      type(proctree), pointer::ptree

      t1 = MPI_Wtime()

      strlen = 1
      ! do while(trans(strlen) /= c_null_char)
      ! strlen = strlen + 1
      ! enddo
      ! strlen = strlen -1
      allocate (character(len=strlen) :: str)
      str = transfer(trans(1:strlen), str)

      call c_f_pointer(bf_for_Cptr, blocks)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(ptree_Cptr, ptree)
      stats%Flop_Tmp = 0
      stats%Flop_C_Mult = 0
      stats%Time_C_Mult = 0
      xout = 0
      if (trim(str) == 'N') then
         call BF_block_MVP_dat(blocks, trim(str), Noutloc, Ninloc, Ncol, xin, Ninloc, xout,Noutloc, BPACK_cone, BPACK_czero, ptree, stats)
      else
         call BF_block_MVP_dat(blocks, trim(str), Ninloc, Noutloc, Ncol, xin, Ninloc, xout, Noutloc, BPACK_cone, BPACK_czero, ptree, stats)
      endif

      t2 = MPI_Wtime()

      xout = xout/option%scale_factor

      stats%Time_C_Mult = stats%Time_C_Mult + t2 - t1
      stats%Flop_C_Mult = stats%Flop_C_Mult + stats%Flop_Tmp

      ! write(*,*)t2-t1
      deallocate (str)
   end subroutine C_BF_Mult

!>**** C interface of parallel extraction of a list of intersections from a block
   !> @param block_Cptr: the structure containing the block
   !> @param option_Cptr: the structure containing option
   !> @param stats_Cptr: the structure containing statistics
   !> @param ptree_Cptr: the structure containing process tree
   !> @param msh_Cptr: the structure containing points and ordering information
   !> @param pgidx: 1D array containing the process group number of each intersection, the number starts from 0
   !> @param Npmap: number of process groups
   !> @param pmaps: 2D array (Npmapx3) containing number of process rows, number of process columns, and starting process ID in each intersection
   !> @param Ninter: number of intersections
   !> @param allrows: 1D array containing the global row indices (in original order starting from 1 to M) stacked together
   !> @param allcols: 1D array containing the global column indices (in original order starting from 1 to N) stacked together
   !> @param alldat_loc: 1D array containing the local entry values defined by pmaps (in column major) stacked together
   !> @param rowidx: 1D array containing sizes of rows of each intersection
   !> @param colidx: 1D array containing sizes of columns of each intersection
   !> @param Nallrows: total number of rows
   !> @param Nallcols: total number of columns
   !> @param Nalldat_loc: total number of local entries
   subroutine C_BF_ExtractElement(block_Cptr, option_Cptr, msh_Cptr, stats_Cptr, ptree_Cptr, Ninter, Nallrows, Nallcols, Nalldat_loc, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps) bind(c, name="c_bf_extractelement")
      use BPACK_DEFS
      implicit none

      type(intersect), allocatable::inters(:)
      real(kind=8)::n0, n1, n2, n3, n4, n5
      integer Ntest, passflag
      integer Nallrows, Nallcols, Nalldat_loc
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
      integer:: allrows(Nallrows), allcols(Nallcols)
      integer, allocatable::datidx(:)
      integer:: Ninter, nr, nc, ntot_loc
      DT,target::alldat_loc(Nalldat_loc)
      integer::colidx(Ninter), rowidx(Ninter), pgidx(Ninter)
      integer::Npmap, pmaps(Npmap, 3), flag2D
      type(iarray)::lst
      type(matrixblock), pointer::blocks_o
      integer, allocatable:: allrows1(:), allcols1(:), pgidx1(:)

      type(c_ptr), intent(in) :: block_Cptr
      type(c_ptr), intent(in) :: ptree_Cptr
      type(c_ptr), intent(in) :: option_Cptr
      type(c_ptr), intent(in) :: msh_Cptr
      type(c_ptr), intent(in) :: stats_Cptr

      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(Bmatrix), pointer::bmat
      type(proctree), pointer::ptree
      type(mesh), pointer::msh
      real(kind=8) t1, t2

      t1 = MPI_Wtime()

      call c_f_pointer(block_Cptr, blocks_o)
      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)
      call c_f_pointer(msh_Cptr, msh)

      stats%Flop_C_Extract = 0
      stats%Time_C_Extract = 0

      idx_row = sum(rowidx)
      allocate (allrows1(idx_row))
      allrows1 = allrows
      idx_col = sum(colidx)
      allocate (allcols1(idx_col))
      allcols1 = allcols + blocks_o%M
      allocate (pgidx1(Ninter))
      pgidx1 = pgidx + 1

      ! do ii=1,idx_row
      ! write(*,*)'row',ii,allrows1(ii)
      ! enddo

      ! do jj=1,idx_col
      ! write(*,*)'col',jj,allcols1(jj)
      ! enddo

      call BF_ExtractElement(blocks_o, option, msh, stats, ptree, Ninter, allrows1, allcols1, alldat_loc, rowidx, colidx, pgidx1, Npmap, pmaps)

      deallocate (allrows1)
      deallocate (allcols1)
      deallocate (pgidx1)

      t2 = MPI_Wtime()

      stats%Time_C_Extract = stats%Time_C_Extract + t2 - t1
      stats%Flop_C_Extract = stats%Flop_C_Extract + stats%Flop_Tmp

   end subroutine C_BF_ExtractElement

!>**** C interface of parallel extraction of a list of intersections from a Bmat matrix
   !> @param Ninter: number of intersections
   !> @param allrows: 1D array containing the global row indices (in original order starting from 0) stacked together
   !> @param allcols: 1D array containing the global column indices (in original order starting from 0) stacked together
   !> @param alldat_loc: 1D array containing the local entry values defined by pmaps (in column major) stacked together
   !> @param rowidx: 1D array containing sizes of rows of each intersection
   !> @param colidx: 1D array containing sizes of columns of each intersection
   !> @param pgidx: 1D array containing the process group number of each intersection, the number starts from 0
   !> @param Npmap: number of process groups
   !> @param pmaps: 2D array (Npmapx3) containing number of process rows, number of process columns, and starting process ID in each intersection
   !> @param bmat_Cptr: the structure containing HODLR
   !> @param option_Cptr: the structure containing option
   !> @param stats_Cptr: the structure containing statistics
   !> @param ptree_Cptr: the structure containing process tree
   !> @param msh_Cptr: the structure containing points and ordering information
   !> @param Nallrows: total number of rows
   !> @param Nallcols: total number of columns
   !> @param Nalldat_loc: total number of local entries
   subroutine C_BPACK_ExtractElement(bmat_Cptr, option_Cptr, msh_Cptr, stats_Cptr, ptree_Cptr, Ninter, Nallrows, Nallcols, Nalldat_loc, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps) bind(c, name="c_bpack_extractelement")
      use BPACK_DEFS
      implicit none
      type(c_ptr), intent(in) :: bmat_Cptr
      type(c_ptr), intent(in) :: ptree_Cptr
      type(c_ptr), intent(in) :: option_Cptr
      type(c_ptr), intent(in) :: msh_Cptr
      type(c_ptr), intent(in) :: stats_Cptr

      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(Bmatrix), pointer::bmat
      type(proctree), pointer::ptree
      type(mesh), pointer::msh

      integer Nallrows, Nallcols, Nalldat_loc
      integer:: allrows(Nallrows), allcols(Nallcols)
      integer, allocatable:: allrows1(:), allcols1(:), pgidx1(:)
      integer:: Ninter, idx_row, idx_col
      DT,target::alldat_loc(Nalldat_loc)
      integer::colidx(Ninter), rowidx(Ninter), pgidx(Ninter)
      integer::Npmap, pmaps(Npmap, 3)
      real(kind=8) t1, t2

      t1 = MPI_Wtime()

      call c_f_pointer(bmat_Cptr, bmat)
      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)
      call c_f_pointer(msh_Cptr, msh)

      stats%Flop_C_Extract = 0
      stats%Time_C_Extract = 0

      idx_row = sum(rowidx)
      allocate (allrows1(idx_row))
      allrows1 = allrows
      idx_col = sum(colidx)
      allocate (allcols1(idx_col))
      allcols1 = allcols
      allocate (pgidx1(Ninter))
      pgidx1 = pgidx + 1

      call BPACK_ExtractElement(bmat, option, msh, stats, ptree, Ninter, allrows1, allcols1, alldat_loc, rowidx, colidx, pgidx1, Npmap, pmaps)

      deallocate (allrows1)
      deallocate (allcols1)
      deallocate (pgidx1)

      t2 = MPI_Wtime()

      stats%Time_C_Extract = stats%Time_C_Extract + t2 - t1
      stats%Flop_C_Extract = stats%Flop_C_Extract + stats%Flop_Tmp

   end subroutine C_BPACK_ExtractElement

!>**** C interface of HODLR-vector multiplication
   !> @param xin: input vector
   !> @param Ninloc: size of local input vectors
   !> @param xout: output vector
   !> @param Noutloc: size of local output vectors
   !> @param Ncol: number of vectors
   !> @param bmat_Cptr: the structure containing HODLR
   !> @param option_Cptr: the structure containing option
   !> @param stats_Cptr: the structure containing statistics
   !> @param ptree_Cptr: the structure containing process tree
   !> @param trans: 'N', 'C' or 'T'
   subroutine C_BPACK_Mult(trans, xin, xout, Ninloc, Noutloc, Ncol, bmat_Cptr, option_Cptr, stats_Cptr, ptree_Cptr) bind(c, name="c_bpack_mult")
      implicit none
      real(kind=8) t1, t2
      integer Ninloc, Noutloc, Ncol
      DT::xin(Ninloc, Ncol), xout(Noutloc, Ncol)

      character(kind=c_char, len=1) :: trans(*)
      type(c_ptr), intent(in) :: bmat_Cptr
      type(c_ptr), intent(in) :: ptree_Cptr
      type(c_ptr), intent(in) :: option_Cptr
      type(c_ptr), intent(in) :: stats_Cptr

      integer strlen
      character(len=:), allocatable :: str
      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(Bmatrix), pointer::bmat
      type(proctree), pointer::ptree

      t1 = MPI_Wtime()

      strlen = 1
      ! do while(trans(strlen) /= c_null_char)
      ! strlen = strlen + 1
      ! enddo
      ! strlen = strlen -1
      allocate (character(len=strlen) :: str)
      str = transfer(trans(1:strlen), str)

      call c_f_pointer(bmat_Cptr, bmat)

      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)

      stats%Flop_C_Mult = 0
      stats%Time_C_Mult = 0

      ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Multiply ......"

      call assert(Noutloc == Ninloc, "not square Z")
      call BPACK_Mult(trim(str), Noutloc, Ncol, xin, xout, bmat, ptree, option, stats)
      ! need to use another Flop counter for this operation in future

      ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Multiply finished"
      ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "

      t2 = MPI_Wtime()

      stats%Time_C_Mult = stats%Time_C_Mult + t2 - t1
      stats%Flop_C_Mult = stats%Flop_C_Mult + stats%Flop_Tmp

      ! write(*,*)t2-t1
      deallocate (str)
   end subroutine C_BPACK_Mult

!>**** C interface of HODLR(inverse)-vector multiplication
   !> @param xin: input vector
   !> @param Ninloc: size of local input vectors
   !> @param xout: output vector
   !> @param Noutloc: size of local output vectors
   !> @param Ncol: number of vectors
   !> @param bmat_Cptr: the structure containing HODLR
   !> @param option_Cptr: the structure containing option
   !> @param stats_Cptr: the structure containing statistics
   !> @param ptree_Cptr: the structure containing process tree
   !> @param trans: 'C', 'T' or 'N'
   subroutine C_BPACK_Inv_Mult(trans, xin, xout, Ninloc, Noutloc, Ncol, bmat_Cptr, option_Cptr, stats_Cptr, ptree_Cptr) bind(c, name="c_bpack_inv_mult")
      implicit none
      real(kind=8) t1, t2
      integer Ninloc, Noutloc, Ncol
      DT::xin(Ninloc, Ncol), xout(Noutloc, Ncol)

      character(kind=c_char, len=1) :: trans(*)
      type(c_ptr), intent(in) :: bmat_Cptr
      type(c_ptr), intent(in) :: ptree_Cptr
      type(c_ptr), intent(in) :: option_Cptr
      type(c_ptr), intent(in) :: stats_Cptr

      integer strlen
      character(len=:), allocatable :: str
      type(Hoption), pointer::option
      type(Hstat), pointer::stats
      type(Bmatrix), pointer::bmat

      type(proctree), pointer::ptree

      t1 = MPI_Wtime()

      strlen = 1
      ! do while(trans(strlen) /= c_null_char)
      ! strlen = strlen + 1
      ! enddo
      ! strlen = strlen -1
      allocate (character(len=strlen) :: str)
      str = transfer(trans(1:strlen), str)

      call c_f_pointer(bmat_Cptr, bmat)

      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)

      stats%Flop_C_Mult = 0
      stats%Time_C_Mult = 0

      ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Multiply ......"

      call assert(Noutloc == Ninloc, "not square Z")
      call BPACK_Inv_Mult(trim(str), Noutloc, Ncol, xin, xout, bmat, ptree, option, stats)
      ! need to use another Flop counter for this operation in future

      ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Multiply finished"
      ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "

      t2 = MPI_Wtime()

      stats%Time_C_Mult = stats%Time_C_Mult + t2 - t1
      stats%Flop_C_Mult = stats%Flop_C_Mult + stats%Flop_Tmp

      ! write(*,*)t2-t1
      deallocate (str)
   end subroutine C_BPACK_Inv_Mult

!>**** C interface of deleting statistics
   !> @param stats_Cptr: the structure containing statistics
   subroutine C_BPACK_Deletestats(stats_Cptr) bind(c, name="c_bpack_deletestats")
      implicit none
      type(c_ptr), intent(inout) :: stats_Cptr
      type(Hstat), pointer::stats

      call c_f_pointer(stats_Cptr, stats)
      call delete_Hstat(stats)
      deallocate (stats)
      stats_Cptr = c_null_ptr

   end subroutine C_BPACK_Deletestats

!>**** C interface of deleting process tree
   !> @param ptree_Cptr: the structure containing process tree
   subroutine C_BPACK_Deleteproctree(ptree_Cptr) bind(c, name="c_bpack_deleteproctree")
      implicit none
      type(c_ptr), intent(inout) :: ptree_Cptr
      type(proctree), pointer::ptree

      call c_f_pointer(ptree_Cptr, ptree)
      call delete_proctree(ptree)
      deallocate (ptree)
      ptree_Cptr = c_null_ptr

   end subroutine C_BPACK_Deleteproctree

!>**** C interface of deleting mesh
   !> @param msh_Cptr: the structure containing mesh
   subroutine C_BPACK_Deletemesh(msh_Cptr) bind(c, name="c_bpack_deletemesh")
      implicit none
      type(c_ptr), intent(inout) :: msh_Cptr
      type(mesh), pointer::msh

      call c_f_pointer(msh_Cptr, msh)
      call delete_mesh(msh)
      deallocate (msh)
      msh_Cptr = c_null_ptr

   end subroutine C_BPACK_Deletemesh

!>**** C interface of deleting kernelquant
   !> @param ker_Cptr: the structure containing kernelquant
   subroutine C_BPACK_Deletekernelquant(ker_Cptr) bind(c, name="c_bpack_deletekernelquant")
      implicit none
      type(c_ptr), intent(inout) :: ker_Cptr
      type(kernelquant), pointer::ker

      call c_f_pointer(ker_Cptr, ker)
      call delete_kernelquant(ker)
      deallocate (ker)
      ker_Cptr = c_null_ptr

   end subroutine C_BPACK_Deletekernelquant

!>**** C interface of deleting HOBF
   !> @param bmat_Cptr: the structure containing HOBF
   subroutine C_BPACK_Delete(bmat_Cptr) bind(c, name="c_bpack_delete")
      implicit none
      type(c_ptr), intent(inout) :: bmat_Cptr
      type(Bmatrix), pointer::bmat

      call c_f_pointer(bmat_Cptr, bmat)
      call BPACK_delete(bmat)
      deallocate (bmat)
      bmat_Cptr = c_null_ptr

   end subroutine C_BPACK_Delete

!>**** C interface of deleting a BF
   !> @param bf_Cptr: the structure containing BF
   subroutine C_BF_DeleteBF(bf_Cptr) bind(c, name="c_bf_deletebf")
      implicit none
      type(c_ptr), intent(inout) :: bf_Cptr
      type(matrixblock), pointer::blocks

      call c_f_pointer(bf_Cptr, blocks)
      call BF_delete(blocks, 1)
      deallocate (blocks)
      bf_Cptr = c_null_ptr
   end subroutine C_BF_DeleteBF

!>**** C interface of deleting Hoption
   !> @param option_Cptr: the structure containing Hoption
   subroutine C_BPACK_Deleteoption(option_Cptr) bind(c, name="c_bpack_deleteoption")
      implicit none
      type(c_ptr), intent(inout) :: option_Cptr
      type(Hoption), pointer::option

      call c_f_pointer(option_Cptr, option)
      deallocate (option)
      option_Cptr = c_null_ptr

   end subroutine C_BPACK_Deleteoption

!>**** C interface of getting the version number of ButterflyPACK
!> @param v_major: major version number
!> @param v_minor: minor version number
!> @param v_bugfix: bugfix version number
   subroutine C_BPACK_GetVersionNumber(v_major, v_minor, v_bugfix) bind(c, name="c_bpack_getversionnumber")
      implicit none
      integer::v_major, v_minor, v_bugfix

      call BPACK_GetVersionNumber(v_major, v_minor, v_bugfix)

   end subroutine C_BPACK_GetVersionNumber

!>**** C interface of converting the tree index to the index in one of its two child trees
!> @param idx_merge: node number in the parent tree
!> @param idx_child: node number in one of the two child trees
   subroutine C_BPACK_TreeIndex_Merged2Child(idx_merge, idx_child) bind(c, name="c_bpack_treeindex_merged2child")
      implicit none
      integer idx_merge, idx_child ! node index in the merged and child tree
      integer level, nth  ! level and the order on that level in the merged tree
      integer level_c, nth_c ! level and the order on that level in the child tree
      integer l_or_r ! 1: left 0: right
      level = GetTreelevel(idx_merge) - 1
      call assert(level >= 1, 'cannot convert the root node')

      level_c = level - 1
      nth = idx_merge - 2**level + 1
      nth_c = nth
      l_or_r = 1
      if (nth_c > 2**(level - 1)) then
         nth_c = nth_c - 2**(level - 1) ! this is an index in the right child tree
         l_or_r = 0
      endif

      idx_child = 2**level_c + nth_c - 1
      if (l_or_r == 0) idx_child = -idx_child

   end subroutine C_BPACK_TreeIndex_Merged2Child

end module BPACK_wrapper

