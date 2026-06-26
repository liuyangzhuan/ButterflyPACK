! Fortran side of the PaRSEC PTG bridge.  This file is compiled when
! enable_parsec is ON.

#include "cButterflyPACK_config.fi"
module c_BPACK_hmat_ptg_bridge
   use iso_c_binding
   use c_BPACK_DEFS
   use c_BPACK_factor
   use c_BPACK_block_sendrecv
   use c_BPACK_Utilities
   use c_MISC_Utilities
   use c_Bplus_Utilities
contains

   subroutine c_c_bpack_hmat_ptg_rank(ptree_Cptr, rank) bind(c, name="c_c_bpack_hmat_ptg_rank")
      implicit none
      type(c_ptr), intent(in) :: ptree_Cptr
      integer(c_int), intent(out) :: rank
      type(proctree), pointer :: ptree

      call c_f_pointer(ptree_Cptr, ptree)
      rank = ptree%MyID
   end subroutine c_c_bpack_hmat_ptg_rank

   subroutine c_c_bpack_hmat_ptg_size(ptree_Cptr, size) bind(c, name="c_c_bpack_hmat_ptg_size")
      implicit none
      type(c_ptr), intent(in) :: ptree_Cptr
      integer(c_int), intent(out) :: size
      type(proctree), pointer :: ptree

      call c_f_pointer(ptree_Cptr, ptree)
      size = ptree%nproc
   end subroutine c_c_bpack_hmat_ptg_size

   subroutine c_c_bpack_hmat_ptg_owner_of(hmat_Cptr, ptree_Cptr, row0, col0, owner) bind(c, name="c_c_bpack_hmat_ptg_owner_of")
      implicit none
      type(c_ptr), intent(in) :: hmat_Cptr
      type(c_ptr), intent(in) :: ptree_Cptr
      integer(c_int), intent(in) :: row0, col0
      integer(c_int), intent(out) :: owner
      type(Hmat), pointer :: h_mat
      type(proctree), pointer :: ptree
      integer :: nprow, npcol, myrow, mycol
      integer :: iproc, jproc, myi, myj

      call c_f_pointer(hmat_Cptr, h_mat)
      call c_f_pointer(ptree_Cptr, ptree)
      call blacs_gridinfo_wrp(ptree%pgrp(1)%ctxt, nprow, npcol, myrow, mycol)
      call g2l(row0 + 1, 2**h_mat%Dist_level, nprow, 1, iproc, myi)
      call g2l(col0 + 1, 2**h_mat%Dist_level, npcol, 1, jproc, myj)
      owner = blacs_pnum_wp(nprow, npcol, iproc, jproc)
   end subroutine c_c_bpack_hmat_ptg_owner_of

   subroutine c_c_bpack_hmat_ptg_diag(hmat_Cptr, option_Cptr, stats_Cptr, ptree_Cptr, msh_Cptr, k0, out_ptr, out_bytes) bind(c, name="c_c_bpack_hmat_ptg_diag")
      implicit none
      type(c_ptr), intent(in) :: hmat_Cptr, option_Cptr, stats_Cptr, ptree_Cptr, msh_Cptr
      integer(c_int), intent(in) :: k0
      type(c_ptr), intent(out) :: out_ptr
      integer(c_int64_t), intent(out) :: out_bytes
      type(Hmat), pointer :: h_mat
      type(Hoption), pointer :: option
      type(Hstat), pointer :: stats
      type(proctree), pointer :: ptree
      type(mesh), pointer :: msh
      type(matrixblock), pointer :: block
      DT, pointer :: buffer(:)
      integer :: k, myi, myj, count_dt

      call c_f_pointer(hmat_Cptr, h_mat)
      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)
      call c_f_pointer(msh_Cptr, msh)

      k = k0 + 1
      if (option%verbosity >= 0) write (*, *) "starting panel ", k
      call bpack_hmat_ptg_local_indices(h_mat, ptree, k, k, myi, myj)
      block => h_mat%Local_blocks(myj, myi)
      call Hmat_LU(block, h_mat, option, stats, ptree, msh)
      h_mat%phase = h_mat%phase * block%phase
      h_mat%logabsdet = h_mat%logabsdet + block%logabsdet
      call pack_all_blocks_one_node(block, msh, option)
      call bpack_hmat_ptg_pack_block(block, msh, buffer, count_dt)
      out_ptr = c_loc(buffer(1))
      out_bytes = int(count_dt, c_int64_t)*int(storage_size(buffer(1))/8, c_int64_t)
   end subroutine c_c_bpack_hmat_ptg_diag

   subroutine c_c_bpack_hmat_ptg_upanel(hmat_Cptr, option_Cptr, stats_Cptr, ptree_Cptr, msh_Cptr, k0, j0, diag_ptr, out_ptr, out_bytes) bind(c, name="c_c_bpack_hmat_ptg_upanel")
      implicit none
      type(c_ptr), intent(in) :: hmat_Cptr, option_Cptr, stats_Cptr, ptree_Cptr, msh_Cptr
      integer(c_int), intent(in) :: k0, j0
      type(c_ptr), value :: diag_ptr
      type(c_ptr), intent(out) :: out_ptr
      integer(c_int64_t), intent(out) :: out_bytes
      type(Hmat), pointer :: h_mat
      type(Hoption), pointer :: option
      type(Hstat), pointer :: stats
      type(proctree), pointer :: ptree
      type(mesh), pointer :: msh
      type(matrixblock), pointer :: diag_block, u_block
      DT, pointer :: buffer(:)
      integer :: k, j, myi, myj, count_dt

      call c_f_pointer(hmat_Cptr, h_mat)
      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)
      call c_f_pointer(msh_Cptr, msh)

      k = k0 + 1
      j = j0 + 1
      allocate(diag_block)
      call bpack_hmat_ptg_unpack_block(diag_block, msh, diag_ptr)
      call unpack_all_blocks_one_node(diag_block, h_mat%Maxlevel, ptree, msh, bpack_hmat_ptg_pgno(h_mat, ptree, k, j), option)

      call bpack_hmat_ptg_local_indices(h_mat, ptree, k, j, myi, myj)
      u_block => h_mat%Local_blocks(myj, myi)
      call Hmat_LXM(diag_block, u_block, h_mat, option, stats, ptree, msh)
      call pack_all_blocks_one_node(u_block, msh, option)
      call bpack_hmat_ptg_pack_block(u_block, msh, buffer, count_dt)

      call Hmat_block_delete(diag_block)
      deallocate(diag_block)
      out_ptr = c_loc(buffer(1))
      out_bytes = int(count_dt, c_int64_t)*int(storage_size(buffer(1))/8, c_int64_t)
   end subroutine c_c_bpack_hmat_ptg_upanel

   subroutine c_c_bpack_hmat_ptg_lpanel(hmat_Cptr, option_Cptr, stats_Cptr, ptree_Cptr, msh_Cptr, k0, i0, diag_ptr, out_ptr, out_bytes) bind(c, name="c_c_bpack_hmat_ptg_lpanel")
      implicit none
      type(c_ptr), intent(in) :: hmat_Cptr, option_Cptr, stats_Cptr, ptree_Cptr, msh_Cptr
      integer(c_int), intent(in) :: k0, i0
      type(c_ptr), value :: diag_ptr
      type(c_ptr), intent(out) :: out_ptr
      integer(c_int64_t), intent(out) :: out_bytes
      type(Hmat), pointer :: h_mat
      type(Hoption), pointer :: option
      type(Hstat), pointer :: stats
      type(proctree), pointer :: ptree
      type(mesh), pointer :: msh
      type(matrixblock), pointer :: diag_block, l_block
      DT, pointer :: buffer(:)
      integer :: k, i, myi, myj, count_dt

      call c_f_pointer(hmat_Cptr, h_mat)
      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)
      call c_f_pointer(msh_Cptr, msh)

      k = k0 + 1
      i = i0 + 1
      allocate(diag_block)
      call bpack_hmat_ptg_unpack_block(diag_block, msh, diag_ptr)
      call unpack_all_blocks_one_node(diag_block, h_mat%Maxlevel, ptree, msh, bpack_hmat_ptg_pgno(h_mat, ptree, i, k), option)

      call bpack_hmat_ptg_local_indices(h_mat, ptree, i, k, myi, myj)
      l_block => h_mat%Local_blocks(myj, myi)
      call Hmat_XUM(diag_block, l_block, h_mat, option, stats, ptree, msh)
      call pack_all_blocks_one_node(l_block, msh, option)
      call bpack_hmat_ptg_pack_block(l_block, msh, buffer, count_dt)

      call Hmat_block_delete(diag_block)
      deallocate(diag_block)
      out_ptr = c_loc(buffer(1))
      out_bytes = int(count_dt, c_int64_t)*int(storage_size(buffer(1))/8, c_int64_t)
   end subroutine c_c_bpack_hmat_ptg_lpanel

   subroutine c_c_bpack_hmat_ptg_update(hmat_Cptr, option_Cptr, stats_Cptr, ptree_Cptr, msh_Cptr, k0, i0, j0, l_ptr, u_ptr) bind(c, name="c_c_bpack_hmat_ptg_update")
      implicit none
      type(c_ptr), intent(in) :: hmat_Cptr, option_Cptr, stats_Cptr, ptree_Cptr, msh_Cptr
      integer(c_int), intent(in) :: k0, i0, j0
      type(c_ptr), value :: l_ptr, u_ptr
      type(Hmat), pointer :: h_mat
      type(Hoption), pointer :: option
      type(Hstat), pointer :: stats
      type(proctree), pointer :: ptree
      type(mesh), pointer :: msh
      type(matrixblock), pointer :: l_block, u_block, a_block
      real(kind=8) :: memory
      integer :: i, j, myi, myj, pgno

      call c_f_pointer(hmat_Cptr, h_mat)
      call c_f_pointer(option_Cptr, option)
      call c_f_pointer(stats_Cptr, stats)
      call c_f_pointer(ptree_Cptr, ptree)
      call c_f_pointer(msh_Cptr, msh)

      i = i0 + 1
      j = j0 + 1
      pgno = bpack_hmat_ptg_pgno(h_mat, ptree, i, j)

      allocate(l_block)
      call bpack_hmat_ptg_unpack_block(l_block, msh, l_ptr)
      call unpack_all_blocks_one_node(l_block, h_mat%Maxlevel, ptree, msh, pgno, option)
      memory = 0d0
      call Hmat_block_ComputeMemory(l_block, memory)
      call LogMemory(stats, memory)

      allocate(u_block)
      call bpack_hmat_ptg_unpack_block(u_block, msh, u_ptr)
      call unpack_all_blocks_one_node(u_block, h_mat%Maxlevel, ptree, msh, pgno, option)
      memory = 0d0
      call Hmat_block_ComputeMemory(u_block, memory)
      call LogMemory(stats, memory)

      call bpack_hmat_ptg_local_indices(h_mat, ptree, i, j, myi, myj)
      a_block => h_mat%Local_blocks(myj, myi)
      call Hmat_add_multiply(a_block, '-', l_block, u_block, h_mat, option, stats, ptree, msh)

      memory = 0d0
      call Hmat_block_ComputeMemory(l_block, memory)
      call LogMemory(stats, -memory)
      memory = 0d0
      call Hmat_block_ComputeMemory(u_block, memory)
      call LogMemory(stats, -memory)
      call Hmat_block_delete(l_block)
      call Hmat_block_delete(u_block)
      deallocate(l_block)
      deallocate(u_block)
   end subroutine c_c_bpack_hmat_ptg_update

   subroutine c_c_bpack_hmat_ptg_release_buffer(ptr) bind(c, name="c_c_bpack_hmat_ptg_release_buffer")
      implicit none
      type(c_ptr), intent(inout) :: ptr
      DT, pointer :: buffer(:)

      if (c_associated(ptr)) then
         call c_f_pointer(ptr, buffer, [1])
         deallocate(buffer)
         ptr = c_null_ptr
      endif
   end subroutine c_c_bpack_hmat_ptg_release_buffer

   subroutine bpack_hmat_ptg_local_indices(h_mat, ptree, grow, gcol, myi, myj)
      implicit none
      type(Hmat) :: h_mat
      type(proctree) :: ptree
      integer, intent(in) :: grow, gcol
      integer, intent(out) :: myi, myj
      integer :: nprow, npcol, myrow, mycol, iproc, jproc

      call blacs_gridinfo_wrp(ptree%pgrp(1)%ctxt, nprow, npcol, myrow, mycol)
      call g2l(grow, 2**h_mat%Dist_level, nprow, 1, iproc, myi)
      call g2l(gcol, 2**h_mat%Dist_level, npcol, 1, jproc, myj)
   end subroutine bpack_hmat_ptg_local_indices

   integer function bpack_hmat_ptg_pgno(h_mat, ptree, grow, gcol)
      implicit none
      type(Hmat) :: h_mat
      type(proctree) :: ptree
      integer, intent(in) :: grow, gcol
      integer :: myi, myj

      call bpack_hmat_ptg_local_indices(h_mat, ptree, grow, gcol, myi, myj)
      bpack_hmat_ptg_pgno = h_mat%Local_blocks(myj, myi)%pgno
   end function bpack_hmat_ptg_pgno

   subroutine bpack_hmat_ptg_pack_block(block, msh, buffer, count_dt)
      implicit none
      type(matrixblock), pointer :: block
      type(mesh) :: msh
      DT, pointer :: buffer(:)
      integer, intent(out) :: count_dt
      integer :: send_count_ind, send_count_dat, success

      send_count_ind = 0
      send_count_dat = 0
      call bpack_hmat_ptg_structure2array(block, send_count_ind, send_count_dat, 0, msh, buffer)
      count_dt = send_count_ind + send_count_dat + 2
      allocate(buffer(count_dt), stat=success)
      call assert(success == 0, 'PTG packed buffer allocation failed')

      buffer(1) = send_count_ind
      buffer(2) = count_dt
      send_count_dat = send_count_ind + 2
      send_count_ind = 2
      call bpack_hmat_ptg_structure2array(block, send_count_ind, send_count_dat, 1, msh, buffer)
      call assert(send_count_dat == count_dt, 'PTG packed buffer count mismatch')
   end subroutine bpack_hmat_ptg_pack_block

   subroutine bpack_hmat_ptg_unpack_block(block, msh, buffer_ptr)
      implicit none
      type(matrixblock), pointer :: block
      type(mesh) :: msh
      type(c_ptr), value :: buffer_ptr
      DT, pointer :: buffer(:), buffer_header(:)
      integer :: recv_count_ind, recv_count_dat, count_dt

      call c_f_pointer(buffer_ptr, buffer_header, [2])
      count_dt = nint(dble(buffer_header(2)))
      call assert(count_dt >= 2, 'PTG packed buffer length invalid')
      call c_f_pointer(buffer_ptr, buffer, [count_dt])
      recv_count_ind = 2
      recv_count_dat = 2 + nint(dble(buffer(1)))
      call bpack_hmat_ptg_array2structure(block, recv_count_ind, recv_count_dat, msh, buffer)
   end subroutine bpack_hmat_ptg_unpack_block

   recursive subroutine bpack_hmat_ptg_structure2array(block, count_ind, count_dat, flag, msh, buffer)
      implicit none
      type(matrixblock) :: block
      integer, intent(inout) :: count_ind, count_dat
      integer, intent(in) :: flag
      type(mesh) :: msh
      DT, pointer :: buffer(:)
      type(matrixblock), pointer :: son
      integer :: style, count1, count2, group_m, group_n, mm, nn

      if (flag == 1) buffer(count_ind + 1:count_ind + MPI_Header) = block%blockinfo_MPI
      count_ind = count_ind + MPI_Header

      style = block%style
      if (style == 4) then
         son => block%sons(1, 1)
         call bpack_hmat_ptg_structure2array(son, count_ind, count_dat, flag, msh, buffer)
         son => block%sons(2, 1)
         call bpack_hmat_ptg_structure2array(son, count_ind, count_dat, flag, msh, buffer)
         son => block%sons(1, 2)
         call bpack_hmat_ptg_structure2array(son, count_ind, count_dat, flag, msh, buffer)
         son => block%sons(2, 2)
         call bpack_hmat_ptg_structure2array(son, count_ind, count_dat, flag, msh, buffer)
      elseif (style == 2) then
         count1 = block%length_Butterfly_index_MPI
         count2 = block%length_Butterfly_data_MPI
         if (flag == 1) buffer(count_ind + 1:count_ind + count1) = block%Butterfly_index_MPI
         count_ind = count_ind + count1
         if (flag == 1) buffer(count_dat + 1:count_dat + count2) = block%Butterfly_data_MPI
         count_dat = count_dat + count2
      elseif (style == 1) then
         group_m = block%row_group
         group_n = block%col_group
         mm = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
         nn = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
         if (flag == 1) buffer(count_dat + 1:count_dat + mm*nn) = block%fullmat_MPI
         count_dat = count_dat + mm*nn
         if (group_m == group_n) then
            if (flag == 1) buffer(count_ind + 1:count_ind + mm) = block%ipiv
            count_ind = count_ind + mm
         endif
      endif
   end subroutine bpack_hmat_ptg_structure2array

   recursive subroutine bpack_hmat_ptg_array2structure(block, count_ind, count_dat, msh, buffer)
      implicit none
      type(matrixblock) :: block
      integer, intent(inout) :: count_ind, count_dat
      type(mesh) :: msh
      DT, pointer :: buffer(:)
      type(matrixblock), pointer :: son
      integer :: style, count1, count2, group_m, group_n, mm, nn, success

      block%blockinfo_MPI = nint(dble(buffer(count_ind + 1:count_ind + MPI_Header)))
      count_ind = count_ind + MPI_Header
      block%level = block%blockinfo_MPI(1)
      block%row_group = block%blockinfo_MPI(2)
      block%col_group = block%blockinfo_MPI(3)
      block%style = block%blockinfo_MPI(5)
      block%level_butterfly = block%blockinfo_MPI(8)
      block%length_Butterfly_index_MPI = block%blockinfo_MPI(9)
      block%length_Butterfly_data_MPI = block%blockinfo_MPI(10)

      style = block%style
      if (style == 4) then
         allocate(block%sons(2, 2), stat=success)
         call assert(success == 0, 'PTG sons allocation failed')
         son => block%sons(1, 1)
         call bpack_hmat_ptg_array2structure(son, count_ind, count_dat, msh, buffer)
         son => block%sons(2, 1)
         call bpack_hmat_ptg_array2structure(son, count_ind, count_dat, msh, buffer)
         son => block%sons(1, 2)
         call bpack_hmat_ptg_array2structure(son, count_ind, count_dat, msh, buffer)
         son => block%sons(2, 2)
         call bpack_hmat_ptg_array2structure(son, count_ind, count_dat, msh, buffer)
      elseif (style == 2) then
         count1 = block%length_Butterfly_index_MPI
         count2 = block%length_Butterfly_data_MPI
         allocate(block%Butterfly_index_MPI(count1), stat=success)
         call assert(success == 0, 'PTG Butterfly_index_MPI allocation failed')
         allocate(block%Butterfly_data_MPI(count2), stat=success)
         call assert(success == 0, 'PTG Butterfly_data_MPI allocation failed')
         block%Butterfly_index_MPI = nint(dble(buffer(count_ind + 1:count_ind + count1)))
         count_ind = count_ind + count1
         block%Butterfly_data_MPI = buffer(count_dat + 1:count_dat + count2)
         count_dat = count_dat + count2
      elseif (style == 1) then
         group_m = block%row_group
         group_n = block%col_group
         mm = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
         nn = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
         allocate(block%fullmat_MPI(mm*nn), stat=success)
         call assert(success == 0, 'PTG fullmat_MPI allocation failed')
         block%fullmat_MPI = buffer(count_dat + 1:count_dat + mm*nn)
         count_dat = count_dat + mm*nn
         if (group_m == group_n) then
            allocate(block%ipiv(mm), stat=success)
            call assert(success == 0, 'PTG ipiv allocation failed')
            block%ipiv = nint(dble(buffer(count_ind + 1:count_ind + mm)))
            count_ind = count_ind + mm
         endif
      endif
   end subroutine bpack_hmat_ptg_array2structure

end module c_BPACK_hmat_ptg_bridge
