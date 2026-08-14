! “ButterflyPACK” Copyright (c) 2018, The Regents of the University of California, through
! Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.

module m_BPACK_utilities
   use c_BPACK_DEFS, only: c_Bmatrix, c_Hmat, c_Hoption, c_Hstat, c_block_ptr, c_blockplus, &
      c_basisgroup, c_butterfly_kerl, c_hobf, c_kernelquant, c_matrixblock, c_mesh, c_proctree, &
      LplusMax, DIRECT, HODLR, HMAT, BLR, HSS, Main_ID, MPI_Wtime
   use z_BPACK_DEFS, only: z_Bmatrix, z_Hmat, z_Hoption, z_Hstat, z_blockplus, &
      z_basisgroup, z_butterfly_kerl, z_hobf, z_kernelquant, z_matrixblock, z_mesh, z_proctree
   use c_BPACK_linkedlist, only: c_append => append, c_list, c_MergeSort
   use c_BPACK_utilities, only: c_BPACK_delete, c_delete_mesh
   use c_Bplus_utilities, only: c_Bplus_delete, c_node_score_block_ptr_row
   use c_BPACK_Solve_Mul, only: c_BPACK_Inv_Mult
   use z_BPACK_Utilities, only: z_delete_kernelquant
   use z_BPACK_Solve_Mul, only: z_BPACK_Inv_Mult, z_BPACK_Mult, z_BPACK_Z_iter_usermatvec_precon
#if HAVE_ZFP   
   use c_MISC_Utilities, only: c_ZFP_Compress
   use z_MISC_Utilities, only: z_ZFP_Compress, z_ZFP_Decompress, z_copymatT
#else 
   use z_MISC_Utilities, only: z_copymatT
#endif

   implicit none

   private
   public :: z2c_BPACK_copy, z2c_BPACK_Solution, z2c_mesh_copy, z2c_CopyOptions
   public :: z2c_Bplus_copy, z2c_BF_copy
   public :: z_blackbox_MVP, c_blackbox_precon_MVP

   type z2c_BPACK_quant
      type(z_Bmatrix), pointer :: z_bmat => null()
      type(c_Bmatrix), pointer :: c_bmat => null()
      type(z_Hoption), pointer :: z_option => null()
      type(z_Hstat), pointer :: z_stats => null()
      type(z_proctree), pointer :: z_ptree => null()
      type(c_Hoption), pointer :: c_option => null()
      type(c_Hstat), pointer :: c_stats => null()
      type(c_proctree), pointer :: c_ptree => null()
   end type z2c_BPACK_quant

contains

   subroutine z2c_CopyOptions(option_i, option_o)
      implicit none
      type(z_Hoption)::option_i
      type(c_Hoption)::option_o

      option_o%format = option_i%format
      option_o%verbosity = option_i%verbosity

      option_o%LRlevel = option_i%LRlevel
      option_o%lnoBP = option_i%lnoBP
      option_o%bp_cnt_lr = option_i%bp_cnt_lr
      option_o%TwoLayerOnly = option_i%TwoLayerOnly
      option_o%touch_para = option_i%touch_para
      option_o%sample_para = option_i%sample_para
      option_o%sample_para_outer = option_i%sample_para_outer
      option_o%pat_comp = option_i%pat_comp
      option_o%use_zfp = option_i%use_zfp
      option_o%use_parsec = option_i%use_parsec
      option_o%use_qtt = option_i%use_qtt

      option_o%Hextralevel = option_i%Hextralevel
      option_o%forwardN15flag = option_i%forwardN15flag
      option_o%tol_comp = option_i%tol_comp
      option_o%Nmin_leaf = option_i%Nmin_leaf
      option_o%nogeo = option_i%nogeo
      option_o%per_geo = option_i%per_geo
      option_o%periods = option_i%periods
      option_o%xyzsort = option_i%xyzsort
      option_o%RecLR_leaf = option_i%RecLR_leaf
      option_o%near_para = option_i%near_para
      option_o%knn_near_para = option_i%knn_near_para
      option_o%scale_factor = option_i%scale_factor
      option_o%rmax = option_i%rmax
      option_o%elem_extract = option_i%elem_extract
      option_o%cpp = option_i%cpp
      option_o%knn = option_i%knn
      option_o%fastsample_tensor = option_i%fastsample_tensor

      option_o%tol_LS = option_i%tol_LS
      option_o%tol_Rdetect = option_i%tol_Rdetect
      option_o%tol_rand = option_i%tol_rand
      option_o%jitter = option_i%jitter
      option_o%iter_solver = option_i%iter_solver
      option_o%powiter = option_i%powiter
      option_o%less_adapt = option_i%less_adapt
      option_o%schulzorder = option_i%schulzorder
      option_o%schulzhardstart = option_i%schulzhardstart
      option_o%schulzsplitlevel = option_i%schulzsplitlevel
      option_o%schulzlevel = option_i%schulzlevel
      option_o%rank0 = option_i%rank0
      option_o%rankrate = option_i%rankrate
      option_o%itermax = option_i%itermax
      option_o%ILU = option_i%ILU
      option_o%Nbundle = option_i%Nbundle

      option_o%tol_itersol = option_i%tol_itersol
      option_o%n_iter = option_i%n_iter
      option_o%precon = option_i%precon

      option_o%level_check = option_i%level_check
      option_o%ErrFillFull = option_i%ErrFillFull
      option_o%ErrSol = option_i%ErrSol
      option_o%BACA_Batch = option_i%BACA_Batch
      option_o%LR_BLK_NUM = option_i%LR_BLK_NUM
   end subroutine z2c_CopyOptions

   subroutine z2c_BPACK_copy(bmat_i, bmat_o, ptree)
      implicit none
      type(z_Bmatrix)::bmat_i
      type(c_Bmatrix)::bmat_o
      type(z_proctree)::ptree

      call c_BPACK_delete(bmat_o)

      bmat_o%Maxlevel = bmat_i%Maxlevel
      if (associated(bmat_i%ho_bf)) then
         if (.not. associated(bmat_o%ho_bf)) allocate(bmat_o%ho_bf)
         call z2c_HODLR_copy(bmat_i%ho_bf, bmat_o%ho_bf, ptree)
      endif
      if (associated(bmat_i%h_mat)) then
         if (.not. associated(bmat_o%h_mat)) allocate(bmat_o%h_mat)
         call z2c_Hmat_copy(bmat_i%h_mat, bmat_o%h_mat, ptree)
      endif
      if (associated(bmat_i%hss_bf)) then
         write(*,*) 'z2c_BPACK_copy: HSS-BF copy not yet implemented'
         stop
      endif
      if (associated(bmat_i%hss_bf_md) .or. associated(bmat_i%h_mat_md)) then
         write(*,*) 'z2c_BPACK_copy: tensor copy not yet implemented'
         stop
      endif
   end subroutine z2c_BPACK_copy

   subroutine z2c_HODLR_copy(ho_bf_i, ho_bf_o, ptree)
      implicit none
      type(z_proctree)::ptree
      type(z_hobf)::ho_bf_i
      type(c_hobf)::ho_bf_o
      integer ii, level_c

      ho_bf_o%Maxlevel = ho_bf_i%Maxlevel
      ho_bf_o%N = ho_bf_i%N
      ho_bf_o%logabsdet = real(ho_bf_i%logabsdet, kind=4)
      ho_bf_o%phase = cmplx(real(ho_bf_i%phase, kind=4), real(aimag(ho_bf_i%phase), kind=4), kind=4)

      allocate(ho_bf_o%levels(ho_bf_o%Maxlevel + 1))
      do level_c = 1, ho_bf_o%Maxlevel + 1
         ho_bf_o%levels(level_c)%level = ho_bf_i%levels(level_c)%level
         ho_bf_o%levels(level_c)%N_block_forward = ho_bf_i%levels(level_c)%N_block_forward
         ho_bf_o%levels(level_c)%N_block_inverse = ho_bf_i%levels(level_c)%N_block_inverse
         ho_bf_o%levels(level_c)%Bidxs = ho_bf_i%levels(level_c)%Bidxs
         ho_bf_o%levels(level_c)%Bidxe = ho_bf_i%levels(level_c)%Bidxe

         allocate(ho_bf_o%levels(level_c)%BP(ho_bf_o%levels(level_c)%N_block_forward))
         allocate(ho_bf_o%levels(level_c)%BP_inverse(ho_bf_o%levels(level_c)%N_block_inverse))
         allocate(ho_bf_o%levels(level_c)%BP_inverse_update(ho_bf_o%levels(level_c)%N_block_forward))
         allocate(ho_bf_o%levels(level_c)%BP_inverse_schur(ho_bf_o%levels(level_c)%N_block_inverse))

         do ii = 1, ho_bf_o%levels(level_c)%N_block_forward
            call z2c_Bplus_copy(ho_bf_i%levels(level_c)%BP(ii), ho_bf_o%levels(level_c)%BP(ii))
            call z2c_Bplus_copy(ho_bf_i%levels(level_c)%BP_inverse_update(ii), &
               ho_bf_o%levels(level_c)%BP_inverse_update(ii))
         enddo
         do ii = 1, ho_bf_o%levels(level_c)%N_block_inverse
            call z2c_Bplus_copy(ho_bf_i%levels(level_c)%BP_inverse(ii), &
               ho_bf_o%levels(level_c)%BP_inverse(ii))
            if (level_c < ho_bf_i%Maxlevel + 1) call z2c_Bplus_copy( &
               ho_bf_i%levels(level_c)%BP_inverse_schur(ii), &
               ho_bf_o%levels(level_c)%BP_inverse_schur(ii))
         enddo
      enddo
   end subroutine z2c_HODLR_copy

   subroutine z2c_Hmat_copy(h_mat_i, h_mat_o, ptree)
      implicit none
      integer ii, jj, bm, bn, level
      type(z_Hmat)::h_mat_i
      type(c_Hmat)::h_mat_o
      type(z_matrixblock), pointer :: blocks_i
      type(c_matrixblock), pointer :: blocks_o
      type(z_proctree)::ptree

      h_mat_o%Maxlevel = h_mat_i%Maxlevel
      h_mat_o%N = h_mat_i%N
      h_mat_o%Dist_level = h_mat_i%Dist_level
      h_mat_o%idxs = h_mat_i%idxs
      h_mat_o%idxe = h_mat_i%idxe
      h_mat_o%myArows = h_mat_i%myArows
      h_mat_o%myAcols = h_mat_i%myAcols
      h_mat_o%logabsdet = real(h_mat_i%logabsdet, kind=4)
      h_mat_o%phase = cmplx(real(h_mat_i%phase, kind=4), real(aimag(h_mat_i%phase), kind=4), kind=4)

      if (associated(h_mat_i%N_p)) then
         allocate(h_mat_o%N_p(size(h_mat_i%N_p, 1), size(h_mat_i%N_p, 2)))
         h_mat_o%N_p = h_mat_i%N_p
      endif

      if (allocated(h_mat_i%basis_group)) then
         allocate(h_mat_o%basis_group(size(h_mat_i%basis_group, 1)))
         do ii = 1, size(h_mat_i%basis_group, 1)
            call z2c_basisgroup_copy(h_mat_i%basis_group(ii), h_mat_o%basis_group(ii))
         enddo
      endif

      if (associated(h_mat_i%Local_blocks)) then
         bm = size(h_mat_i%Local_blocks, 1)
         bn = size(h_mat_i%Local_blocks, 2)
         allocate(h_mat_o%Local_blocks(bm, bn))
         do ii = 1, bm
            do jj = 1, bn
               blocks_i => h_mat_i%Local_blocks(ii, jj)
               blocks_o => h_mat_o%Local_blocks(ii, jj)
               call z2c_Hmat_block_copy('N', blocks_i, blocks_o)
            enddo
         enddo
      endif

      if (associated(h_mat_i%Local_blocks_copy)) then
         bm = size(h_mat_i%Local_blocks_copy, 1)
         bn = size(h_mat_i%Local_blocks_copy, 2)
         allocate(h_mat_o%Local_blocks_copy(bm, bn))
         do ii = 1, bm
            do jj = 1, bn
               blocks_i => h_mat_i%Local_blocks_copy(ii, jj)
               blocks_o => h_mat_o%Local_blocks_copy(ii, jj)
               call z2c_Hmat_block_copy('N', blocks_i, blocks_o)
            enddo
         enddo
      endif

      if (allocated(h_mat_i%colorsets)) then
         allocate(h_mat_o%colorsets(0:h_mat_i%Maxlevel))
         do level = 0, h_mat_i%Maxlevel
            allocate(h_mat_o%colorsets(level)%dat(size(h_mat_i%colorsets(level)%dat)))
            h_mat_o%colorsets(level)%dat = h_mat_i%colorsets(level)%dat
            h_mat_o%colorsets(level)%idx = h_mat_i%colorsets(level)%idx
            h_mat_o%colorsets(level)%num_nods = h_mat_i%colorsets(level)%num_nods
         enddo
      endif

      if (allocated(h_mat_i%fullmat)) then
         allocate(h_mat_o%fullmat(size(h_mat_i%fullmat, 1), size(h_mat_i%fullmat, 2)))
         h_mat_o%fullmat = cmplx(real(h_mat_i%fullmat, kind=4), real(aimag(h_mat_i%fullmat), kind=4), kind=4)
      endif

      if (allocated(h_mat_i%lstblks) .and. associated(h_mat_i%Local_blocks)) then
         allocate(h_mat_o%lstblks(0:h_mat_o%Maxlevel))
         do level = 0, h_mat_o%Maxlevel
            h_mat_o%lstblks(level) = c_list()
         enddo
         bm = size(h_mat_o%Local_blocks, 1)
         bn = size(h_mat_o%Local_blocks, 2)
         do ii = 1, bm
            do jj = 1, bn
               blocks_o => h_mat_o%Local_blocks(ii, jj)
               call z2c_Hmat_GetBlkLst(blocks_o, h_mat_o%lstblks, h_mat_o%Maxlevel)
            enddo
         enddo
         do level = 0, h_mat_o%Maxlevel
            call c_MergeSort(h_mat_o%lstblks(level)%head, c_node_score_block_ptr_row)
         enddo
      endif

      if (ptree%MyID < -1) return
   end subroutine z2c_Hmat_copy

   subroutine z2c_Bplus_copy(bplus_i, bplus_o, memory)
      implicit none
      type(z_blockplus)::bplus_i
      type(c_blockplus)::bplus_o
      integer ll, bb, Nboundall, Ninadmissible
      real(kind=8), optional::memory
      real(kind=8)::rtemp

      call c_Bplus_delete(bplus_o)
      if (present(memory)) memory = 0

      if (.not. associated(bplus_i%LL)) return

      allocate(bplus_o%LL(LplusMax))
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
            allocate(bplus_o%LL(ll)%matrices_block(bplus_i%LL(ll)%Nbound))
            do bb = 1, bplus_i%LL(ll)%Nbound
               call z2c_BF_copy('N', bplus_i%LL(ll)%matrices_block(bb), &
                  bplus_o%LL(ll)%matrices_block(bb), rtemp)
               if (present(memory)) memory = memory + rtemp
            enddo
            if (allocated(bplus_i%LL(ll)%boundary_map)) then
               Nboundall = size(bplus_i%LL(ll)%boundary_map, 1)
               Ninadmissible = size(bplus_i%LL(ll)%boundary_map, 2)
               allocate(bplus_o%LL(ll)%boundary_map(Nboundall, Ninadmissible))
               bplus_o%LL(ll)%boundary_map = bplus_i%LL(ll)%boundary_map
               if (present(memory)) memory = memory + dble(size(bplus_o%LL(ll)%boundary_map))*4d0/1024d3
            endif
         endif
      enddo
   end subroutine z2c_Bplus_copy

   subroutine z2c_BF_copy(trans, block_i, block_o, memory)
      implicit none
      type(z_matrixblock)::block_i
      type(c_matrixblock)::block_o
      character::trans
      real(kind=8), optional::memory
      integer ii, jj, index_i_m, index_j_m, level, level_butterfly
      integer mm, nn, rank
      real(kind=8)::tol_used
      complex(kind=8), allocatable::tmp_z(:, :)

      if (present(memory)) memory = 0

      if (trans == 'N') then
         call z2c_BF_copy_metadata('N', block_i, block_o, memory)
         level_butterfly = block_i%level_butterfly

         if (block_i%style == 2) then
            if (allocated(block_i%ButterflyU%blocks)) then
               if (level_butterfly /= 0) allocate(block_o%ButterflyKerl(level_butterfly))

               block_o%ButterflyV%num_blk = block_i%ButterflyV%num_blk
               block_o%ButterflyV%nblk_loc = block_i%ButterflyV%nblk_loc
               block_o%ButterflyV%idx = block_i%ButterflyV%idx
               block_o%ButterflyV%inc = block_i%ButterflyV%inc
               allocate(block_o%ButterflyV%blocks(block_o%ButterflyV%nblk_loc))
               do jj = 1, block_o%ButterflyV%nblk_loc
                  nn = size(block_i%ButterflyV%blocks(jj)%matrix, 1)
                  rank = size(block_i%ButterflyV%blocks(jj)%matrix, 2)
                  allocate(block_o%ButterflyV%blocks(jj)%matrix(nn, rank))
                  call z2c_matrix_copy(block_i%ButterflyV%blocks(jj)%matrix, &
                     block_o%ButterflyV%blocks(jj)%matrix)
                  if (present(memory)) memory = memory + z2c_cmat_memory(block_o%ButterflyV%blocks(jj)%matrix)
               enddo

               do level = 1, level_butterfly
                  call z2c_BF_copy_ker_level(block_i%ButterflyKerl(level), block_o%ButterflyKerl(level), memory)
               enddo

               block_o%ButterflyU%num_blk = block_i%ButterflyU%num_blk
               block_o%ButterflyU%nblk_loc = block_i%ButterflyU%nblk_loc
               block_o%ButterflyU%idx = block_i%ButterflyU%idx
               block_o%ButterflyU%inc = block_i%ButterflyU%inc
               allocate(block_o%ButterflyU%blocks(block_o%ButterflyU%nblk_loc))
               do ii = 1, block_o%ButterflyU%nblk_loc
                  nn = size(block_i%ButterflyU%blocks(ii)%matrix, 1)
                  rank = size(block_i%ButterflyU%blocks(ii)%matrix, 2)
                  allocate(block_o%ButterflyU%blocks(ii)%matrix(nn, rank))
                  call z2c_matrix_copy(block_i%ButterflyU%blocks(ii)%matrix, &
                     block_o%ButterflyU%blocks(ii)%matrix)
                  if (present(memory)) memory = memory + z2c_cmat_memory(block_o%ButterflyU%blocks(ii)%matrix)
               enddo
            endif

            if (allocated(block_i%ButterflyMiddle)) then
               allocate(block_o%ButterflyMiddle(size(block_i%ButterflyMiddle, 1), size(block_i%ButterflyMiddle, 2)))
               do index_i_m = 1, size(block_i%ButterflyMiddle, 1)
                  do index_j_m = 1, size(block_i%ButterflyMiddle, 2)
                     if (associated(block_i%ButterflyMiddle(index_i_m, index_j_m)%matrix)) then
                        mm = size(block_i%ButterflyMiddle(index_i_m, index_j_m)%matrix, 1)
                        nn = size(block_i%ButterflyMiddle(index_i_m, index_j_m)%matrix, 2)
                        allocate(block_o%ButterflyMiddle(index_i_m, index_j_m)%matrix(mm, nn))
                        call z2c_matrix_copy(block_i%ButterflyMiddle(index_i_m, index_j_m)%matrix, &
                           block_o%ButterflyMiddle(index_i_m, index_j_m)%matrix)
                     endif
                  enddo
               enddo
            endif
         else if (block_i%style == 1) then
            call z2c_BF_copy_fullmat(block_i, block_o, memory)
         endif
      else if (trans == 'T') then
         call z2c_BF_copy_metadata('T', block_i, block_o, memory)
         if (block_i%style == 1) then
            if (associated(block_i%fullmat)) then
               mm = size(block_i%fullmat, 1)
               nn = size(block_i%fullmat, 2)
               allocate(tmp_z(nn, mm))
               call z_copymatT(block_i%fullmat, tmp_z, mm, nn)
               allocate(block_o%fullmat(nn, mm))
               call z2c_matrix_copy(tmp_z, block_o%fullmat)
               deallocate(tmp_z)
               if (present(memory)) memory = memory + z2c_cmat_memory(block_o%fullmat)
            endif
#if HAVE_ZFP
            if (allocated(block_i%FullmatZFP%buffer_r)) then
               write(*,*) 'z2c_BF_copy: transposed ZFP copy not yet implemented'
               stop
            endif
#endif
         else
            write(*,*) 'z2c_BF_copy: transposed butterfly copy not yet implemented'
            stop
         endif
      endif
   end subroutine z2c_BF_copy

   subroutine z2c_BF_copy_metadata(trans, block_i, block_o, memory)
      implicit none
      character::trans
      type(z_matrixblock)::block_i
      type(c_matrixblock)::block_o
      real(kind=8), optional::memory

      block_o%level = block_i%level
      block_o%style = block_i%style
      block_o%level_butterfly = block_i%level_butterfly
      block_o%level_half = block_i%level_half
      block_o%rankmax = block_i%rankmax
      block_o%rankmin = block_i%rankmin
      block_o%dimension_rank = block_i%dimension_rank
      block_o%logabsdet = real(block_i%logabsdet, kind=4)
      block_o%phase = cmplx(real(block_i%phase, kind=4), real(aimag(block_i%phase), kind=4), kind=4)
      block_o%pgno = block_i%pgno
      block_o%pgno_db = block_i%pgno_db

      if (trans == 'N') then
         block_o%col_group = block_i%col_group
         block_o%row_group = block_i%row_group
         block_o%M = block_i%M
         block_o%N = block_i%N
         block_o%headm = block_i%headm
         block_o%headn = block_i%headn
         block_o%M_loc = block_i%M_loc
         block_o%N_loc = block_i%N_loc
         call z2c_int_ptr_copy(block_i%M_p, block_o%M_p, memory)
         call z2c_int_ptr_copy(block_i%N_p, block_o%N_p, memory)
         call z2c_int_vec_ptr_copy(block_i%ms, block_o%ms, memory)
         call z2c_int_vec_ptr_copy(block_i%ns, block_o%ns, memory)
      else
         block_o%col_group = block_i%row_group
         block_o%row_group = block_i%col_group
         block_o%M = block_i%N
         block_o%N = block_i%M
         block_o%headm = block_i%headn
         block_o%headn = block_i%headm
         block_o%M_loc = block_i%N_loc
         block_o%N_loc = block_i%M_loc
         call z2c_int_ptr_copy(block_i%N_p, block_o%M_p, memory)
         call z2c_int_ptr_copy(block_i%M_p, block_o%N_p, memory)
         call z2c_int_vec_ptr_copy(block_i%ns, block_o%ms, memory)
         call z2c_int_vec_ptr_copy(block_i%ms, block_o%ns, memory)
      endif

      if (allocated(block_i%ipiv)) then
         allocate(block_o%ipiv(size(block_i%ipiv, 1)))
         block_o%ipiv = block_i%ipiv
      endif
   end subroutine z2c_BF_copy_metadata

   subroutine z2c_BF_copy_fullmat(block_i, block_o, memory)
      implicit none
      type(z_matrixblock)::block_i
      type(c_matrixblock)::block_o
      real(kind=8), optional::memory
      integer mm, nn
      real(kind=8)::tol_used

      if (associated(block_i%fullmat)) then
         mm = size(block_i%fullmat, 1)
         nn = size(block_i%fullmat, 2)
         allocate(block_o%fullmat(mm, nn))
         call z2c_matrix_copy(block_i%fullmat, block_o%fullmat)
         if (present(memory)) memory = memory + z2c_cmat_memory(block_o%fullmat)
      endif
#if HAVE_ZFP
      if (allocated(block_i%FullmatZFP%buffer_r)) then
         call z_ZFP_Decompress(block_i%fullmat, block_i%FullmatZFP, block_i%M, block_i%N, tol_used, 1)
         mm = size(block_i%fullmat, 1)
         nn = size(block_i%fullmat, 2)
         if (.not. associated(block_o%fullmat)) allocate(block_o%fullmat(mm, nn))
         call z2c_matrix_copy(block_i%fullmat, block_o%fullmat)
         call z_ZFP_Compress(block_i%fullmat, block_i%FullmatZFP, block_i%M, block_i%N, tol_used, 1)
         call c_ZFP_Compress(block_o%fullmat, block_o%FullmatZFP, block_o%M, block_o%N, tol_used, 0)
         if (present(memory) .and. allocated(block_o%FullmatZFP%buffer_r)) &
            memory = memory + dble(size(block_o%FullmatZFP%buffer_r))/1024d3
         if (present(memory) .and. allocated(block_o%FullmatZFP%buffer_i)) &
            memory = memory + dble(size(block_o%FullmatZFP%buffer_i))/1024d3
      endif
#endif
   end subroutine z2c_BF_copy_fullmat

   subroutine z2c_BF_copy_ker_level(ker_i, ker_o, memory)
      implicit none
      type(z_butterfly_kerl)::ker_i
      type(c_butterfly_kerl)::ker_o
      real(kind=8), optional::memory
      integer ii, jj, mm, nn

      ker_o%num_row = ker_i%num_row
      ker_o%num_col = ker_i%num_col
      ker_o%nc = ker_i%nc
      ker_o%nr = ker_i%nr
      ker_o%idx_c = ker_i%idx_c
      ker_o%idx_r = ker_i%idx_r
      ker_o%inc_c = ker_i%inc_c
      ker_o%inc_r = ker_i%inc_r
      if (ker_o%nr > 0 .and. ker_o%nc > 0) then
         allocate(ker_o%blocks(ker_o%nr, ker_o%nc))
         do ii = 1, ker_o%nr
            do jj = 1, ker_o%nc
               if (associated(ker_i%blocks(ii, jj)%matrix)) then
                  mm = size(ker_i%blocks(ii, jj)%matrix, 1)
                  nn = size(ker_i%blocks(ii, jj)%matrix, 2)
                  allocate(ker_o%blocks(ii, jj)%matrix(mm, nn))
                  call z2c_matrix_copy(ker_i%blocks(ii, jj)%matrix, ker_o%blocks(ii, jj)%matrix)
                  if (present(memory)) memory = memory + z2c_cmat_memory(ker_o%blocks(ii, jj)%matrix)
               endif
            enddo
         enddo
      endif
   end subroutine z2c_BF_copy_ker_level

   recursive subroutine z2c_Hmat_block_copy(trans, block_i, block_o, memory)
      implicit none
      character::trans
      type(z_matrixblock), pointer :: block_i, blocks_son_i
      type(c_matrixblock), pointer :: block_o, blocks_son_o
      real(kind=8), optional::memory
      real(kind=8)::memory_tmp
      integer i, j

      block_o%style = block_i%style
      if (block_i%style == 4) then
         call z2c_BF_copy_metadata(trans, block_i, block_o, memory)
         allocate(block_o%sons(2, 2))
         do j = 1, 2
            do i = 1, 2
               block_o%sons(i, j)%father => block_o
            enddo
         enddo
         do j = 1, 2
            do i = 1, 2
               blocks_son_i => block_i%sons(i, j)
               blocks_son_o => block_o%sons(i, j)
               call z2c_Hmat_block_copy(trans, blocks_son_i, blocks_son_o, memory)
            enddo
         enddo
      else
         call z2c_BF_copy(trans, block_i, block_o, memory_tmp)
         if (present(memory)) memory = memory + memory_tmp
      endif
   end subroutine z2c_Hmat_block_copy

   recursive subroutine z2c_Hmat_GetBlkLst(blocks, lstblks, Maxlevel)
      implicit none
      type(c_matrixblock), pointer::blocks, blocks_son
      type(c_list)::lstblks(0:Maxlevel)
      type(c_block_ptr)::blk_ptr
      integer Maxlevel, i, j

      if (blocks%style == 4) then
         do j = 1, 2
            do i = 1, 2
               blocks_son => blocks%sons(i, j)
               call z2c_Hmat_GetBlkLst(blocks_son, lstblks, Maxlevel)
            enddo
         enddo
      else
         blk_ptr%ptr => blocks
         call c_append(lstblks(blocks%level), blk_ptr)
      endif
   end subroutine z2c_Hmat_GetBlkLst

   subroutine z2c_mesh_copy(msh_i, msh_o)
      implicit none
      type(z_mesh)::msh_i
      type(c_mesh)::msh_o
      integer ii

      call c_delete_mesh(msh_o)

      msh_o%Nunk = msh_i%Nunk
      msh_o%Dist_level = msh_i%Dist_level
      msh_o%Maxgroup = msh_i%Maxgroup
      msh_o%idxs = msh_i%idxs
      msh_o%idxe = msh_i%idxe

      if (allocated(msh_i%xyz)) then
         allocate(msh_o%xyz(size(msh_i%xyz, 1), size(msh_i%xyz, 2)))
         msh_o%xyz = msh_i%xyz
      endif
      if (allocated(msh_i%new2old)) then
         allocate(msh_o%new2old(size(msh_i%new2old, 1)))
         msh_o%new2old = msh_i%new2old
      endif
      if (allocated(msh_i%old2new)) then
         allocate(msh_o%old2new(size(msh_i%old2new, 1)))
         msh_o%old2new = msh_i%old2new
      endif
      if (allocated(msh_i%pretree)) then
         allocate(msh_o%pretree(size(msh_i%pretree, 1)))
         msh_o%pretree = msh_i%pretree
      endif
      if (allocated(msh_i%nns)) then
         allocate(msh_o%nns(size(msh_i%nns, 1), size(msh_i%nns, 2)))
         msh_o%nns = msh_i%nns
      endif
      if (allocated(msh_i%basis_group)) then
         allocate(msh_o%basis_group(size(msh_i%basis_group, 1)))
         do ii = 1, size(msh_i%basis_group, 1)
            call z2c_basisgroup_copy(msh_i%basis_group(ii), msh_o%basis_group(ii))
         enddo
      endif
   end subroutine z2c_mesh_copy

   subroutine z2c_BPACK_Solution(z_bmat, c_bmat, x, b, Ns_loc, num_vectors, option, ptree, stats, &
         c_option, c_ptree, c_stats)
      implicit none
      type(z_Bmatrix), target::z_bmat
      type(c_Bmatrix), target::c_bmat
      integer::Ns_loc, num_vectors
      complex(kind=8)::x(Ns_loc, num_vectors), b(Ns_loc, num_vectors)
      type(z_Hoption), target::option
      type(z_proctree), target::ptree
      type(z_Hstat), target::stats
      type(c_Hoption), optional, target::c_option
      type(c_proctree), optional, target::c_ptree
      type(c_Hstat), optional, target::c_stats
      type(z_kernelquant)::ker
      type(z2c_BPACK_quant), target::mvquant
      complex(kind=8), allocatable::rhs(:,:), sol(:,:)
      real(kind=8)::rel_error, n1, n2
      integer ii, iter

      n1 = MPI_Wtime()
      if (option%precon == DIRECT) then
         call z_BPACK_Inv_Mult('N', Ns_loc, num_vectors, b, x, z_bmat, ptree, option, stats)
         stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
      else
         if (.not. present(c_option) .or. .not. present(c_ptree) .or. .not. present(c_stats)) then
            if (ptree%MyID == Main_ID) write(*,*) &
               'z2c_BPACK_Solution needs c_option/c_ptree/c_stats for mixed-precision preconditioning'
            stop
         endif

         mvquant%z_bmat => z_bmat
         mvquant%c_bmat => c_bmat
         mvquant%z_option => option
         mvquant%z_stats => stats
         mvquant%z_ptree => ptree
         mvquant%c_option => c_option
         mvquant%c_stats => c_stats
         mvquant%c_ptree => c_ptree
         ker%QuantApp => mvquant

         allocate(rhs(Ns_loc, 1), sol(Ns_loc, 1))
         do ii = 1, num_vectors
            rhs(:, 1) = b(:, ii)
            sol(:, 1) = x(:, ii)
            rel_error = option%tol_itersol
            iter = 0
            call z_BPACK_Z_iter_usermatvec_precon(option%n_iter, Ns_loc, rhs, sol, rel_error, iter, &
               z_blackbox_MVP, c_blackbox_precon_MVP, ptree, option, stats, ker)
            x(:, ii) = sol(:, 1)
         enddo
         deallocate(rhs, sol)

         call z_delete_kernelquant(ker)
      endif
      n2 = MPI_Wtime()
      stats%Time_Sol = stats%Time_Sol + n2 - n1
   end subroutine z2c_BPACK_Solution

   subroutine z_blackbox_MVP(trans, M, N, num_vect, Vin, Vout, ker)
      implicit none
      character trans
      integer, intent(in)::M, N, num_vect
      complex(kind=8)::Vin(:, :), Vout(:, :)
      type(z_kernelquant)::ker

      select type (mvquant => ker%QuantApp)
      type is (z2c_BPACK_quant)
         if (M /= N) then
            write(*,*) 'z_blackbox_MVP expects a square local operator:', M, N
            stop
         endif
         call z_BPACK_Mult(trans, M, num_vect, Vin, Vout, mvquant%z_bmat, &
            mvquant%z_ptree, mvquant%z_option, mvquant%z_stats)
      class default
         write(*,*) 'unexpected QuantApp type in z_blackbox_MVP'
         stop
      end select
   end subroutine z_blackbox_MVP

   subroutine c_blackbox_precon_MVP(trans, M, N, num_vect, Vin, Vout, ker)
      implicit none
      character trans
      integer, intent(in)::M, N, num_vect
      complex(kind=8)::Vin(:, :), Vout(:, :)
      type(z_kernelquant)::ker
      complex(kind=4), allocatable::vin_sp(:, :), vout_sp(:, :)

      select type (mvquant => ker%QuantApp)
      type is (z2c_BPACK_quant)
         if (M /= N) then
            write(*,*) 'c_blackbox_precon_MVP expects a square local operator:', M, N
            stop
         endif

         allocate(vin_sp(M, num_vect), vout_sp(M, num_vect))
         vin_sp = cmplx(real(Vin(1:M, 1:num_vect), kind=4), &
            real(aimag(Vin(1:M, 1:num_vect)), kind=4), kind=4)
         call c_BPACK_Inv_Mult(trans, M, num_vect, vin_sp, vout_sp, mvquant%c_bmat, &
            mvquant%c_ptree, mvquant%c_option, mvquant%c_stats)
         mvquant%c_stats%Flop_Sol = mvquant%c_stats%Flop_Sol + mvquant%c_stats%Flop_Tmp
         Vout(1:M, 1:num_vect) = cmplx(real(vout_sp, kind=8), &
            real(aimag(vout_sp), kind=8), kind=8)
         deallocate(vin_sp, vout_sp)
      class default
         write(*,*) 'unexpected QuantApp type in c_blackbox_precon_MVP'
         stop
      end select
   end subroutine c_blackbox_precon_MVP

   subroutine z2c_basisgroup_copy(group_i, group_o)
      implicit none
      type(z_basisgroup)::group_i
      type(c_basisgroup)::group_o

      group_o%head = group_i%head
      group_o%tail = group_i%tail
      group_o%radius = group_i%radius
      group_o%boundary = group_i%boundary
      group_o%nn = group_i%nn
      if (allocated(group_i%center)) then
         allocate(group_o%center(size(group_i%center, 1)))
         group_o%center = group_i%center
      endif
      if (allocated(group_i%nlist)) then
         allocate(group_o%nlist(size(group_i%nlist, 1)))
         group_o%nlist = group_i%nlist
      endif
   end subroutine z2c_basisgroup_copy

   subroutine z2c_matrix_copy(mat_i, mat_o)
      implicit none
      complex(kind=8), intent(in)::mat_i(:, :)
      complex(kind=4), intent(out)::mat_o(:, :)

      mat_o = cmplx(real(mat_i, kind=4), real(aimag(mat_i), kind=4), kind=4)
   end subroutine z2c_matrix_copy

   real(kind=8) function z2c_cmat_memory(mat)
      implicit none
      complex(kind=4), intent(in)::mat(:, :)

      z2c_cmat_memory = dble(size(mat))*8d0/1024d3
   end function z2c_cmat_memory

   subroutine z2c_int_ptr_copy(ptr_i, ptr_o, memory)
      implicit none
      integer, pointer::ptr_i(:, :)
      integer, pointer::ptr_o(:, :)
      real(kind=8), optional::memory

      if (associated(ptr_i)) then
         if (associated(ptr_o)) deallocate(ptr_o)
         allocate(ptr_o(size(ptr_i, 1), size(ptr_i, 2)))
         ptr_o = ptr_i
         if (present(memory)) memory = memory + dble(size(ptr_o))*4d0/1024d3
      endif
   end subroutine z2c_int_ptr_copy

   subroutine z2c_int_vec_ptr_copy(ptr_i, ptr_o, memory)
      implicit none
      integer, pointer::ptr_i(:)
      integer, pointer::ptr_o(:)
      real(kind=8), optional::memory

      if (associated(ptr_i)) then
         if (associated(ptr_o)) deallocate(ptr_o)
         allocate(ptr_o(size(ptr_i, 1)))
         ptr_o = ptr_i
         if (present(memory)) memory = memory + dble(size(ptr_o))*4d0/1024d3
      endif
   end subroutine z2c_int_vec_ptr_copy

end module m_BPACK_utilities
