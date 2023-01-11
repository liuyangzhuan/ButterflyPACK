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
!> @file Bplus_pack_unpack_for_MPI.f90
!> @brief Low-level routines for packing, unpacking and communicating butterfly blocks. Only used in the H matrix solver.
!

#include "ButterflyPACK_config.fi"
module BPACK_block_sendrecv
   use BPACK_DEFS
   use BPACK_Utilities
   use MISC_Utilities
contains

   subroutine MPI_verbose_barrier(msg, ptree, option)
      implicit none
      type(proctree)::ptree
      type(Hoption)::option
      character(*)::msg
      integer ierr
      ! write(ptree%MyID+10000,*)'waiting ',msg
      call MPI_barrier(ptree%Comm, ierr)
      if(option%verbosity>=0 .and. ptree%MyID==Main_ID)write(*,*)msg
      ! write(ptree%MyID+10000,*)'waiting done ',msg

   end subroutine MPI_verbose_barrier

   subroutine blocks_partial_bcast(block_s, block_r, send, recv, send_ID, msh, ptree, option)
      use BPACK_DEFS
      use MISC_Utilities
      implicit none
      type(matrixblock), pointer :: block_s, block_r
      integer send, recv, send_ID
      integer groupm, groupn
      integer Local_COMM, Local_Myid, Local_Np, Local_send_ID, color
      integer send_count_ind, send_count_dat, S_request1, S_request2, success, recv_count_ind, recv_count_dat, send_count_tot, recv_count_tot
      type(proctree)::ptree
      type(mesh)::msh
      type(Hoption)::option
      integer ierr
      integer, allocatable::Localmyid2myid(:)

      color = 0
      if (send == 1 .or. recv == 1) color = 1

      call MPI_BARRIER(ptree%Comm, ierr)
      call MPI_COMM_SPLIT(ptree%Comm, color, ptree%MyID, Local_COMM, ierr)
      call MPI_Comm_rank(Local_COMM, Local_Myid, ierr)
      call MPI_Comm_size(Local_COMM, Local_Np, ierr)

      if (color == 1) then
         allocate (Localmyid2myid(Local_Np))
         Localmyid2myid = -1
         call MPI_ALLGATHER(ptree%MyID, 1, MPI_INTEGER, Localmyid2myid, 1, MPI_INTEGER, Local_COMM, ierr)

         Local_send_ID = minloc(abs(Localmyid2myid - send_ID), DIM=1) - 1
         call assert(Localmyid2myid(Local_send_ID + 1) == send_ID, 'Local_send_ID not found')
         ! write(*,*)Local_send_ID,send_ID,Localmyid2myid,num_processors
         if (Local_send_ID == Local_Myid) then
            send_count_ind = 0
            send_count_dat = 0
            groupm = block_s%row_group
            groupn = block_s%col_group
            call blocks_structure2buff(block_s, send_count_ind, send_count_dat, 0, msh, ptree)
            call assert(send_count_dat /= 0, 'send_count_dat==0 not handled yet in blocks_partial_bcast')
            send_count_tot = send_count_ind + send_count_dat + 1
            call assert(send_count_tot < 1d9 .and. send_count_tot > 0, 'send_count_tot incorrect')
            recv_count_tot = send_count_tot
         end if

         call MPI_Bcast(recv_count_tot, 1, MPI_integer, Local_send_ID, Local_COMM, ierr)
         allocate (ptree%recv_buff_dat(recv_count_tot), stat=success)
         call assert(success == 0, 'ptree%recv_buff_dat allocation fail')

         if (Local_send_ID == Local_Myid) then
            allocate (ptree%send_buff_dat(send_count_tot), stat=success)
            call assert(success == 0, 'ptree%send_buff_dat allocation fail')
            ptree%send_buff_dat(1) = send_count_ind
            send_count_dat = send_count_ind + 1
            send_count_ind = 1
            call blocks_structure2buff(block_s, send_count_ind, send_count_dat, 1, msh, ptree)
            ptree%recv_buff_dat = ptree%send_buff_dat
            deallocate (ptree%send_buff_dat)
         end if

         if (option%verbosity >= 2) write (ptree%MyID + 10000, *) 'cast: ', send_ID, ptree%MyID, recv_count_tot
         call MPI_Bcast(ptree%recv_buff_dat, recv_count_tot, MPI_DT, Local_send_ID, Local_COMM, ierr)
         if (option%verbosity >= 2) write (ptree%MyID + 10000, *) 'castF: ', send_ID, ptree%MyID, recv_count_tot

         if (recv == 1) then
            recv_count_ind = 1
            recv_count_dat = 1 + NINT(dble(ptree%recv_buff_dat(1)))
            call blocks_buff2structure(block_r, recv_count_ind, recv_count_dat, msh, ptree)
            call assert(recv_count_dat == recv_count_tot, 'recv_count_dat/=recv_count_tot')
         end if

         deallocate (ptree%recv_buff_dat)
         deallocate (Localmyid2myid)
      end if

      call MPI_COMM_FREE(Local_COMM, ierr)

   end subroutine blocks_partial_bcast

   subroutine blocks_send(block, indices, recv_ID, send_count, msh, ptree, option)

      use BPACK_DEFS
      use MISC_Utilities
      implicit none

      type(proctree)::ptree
      type(mesh)::msh
      type(Hoption)::option
      integer ierr
      integer count1, count2, requests, rank, group_m, group_n
      integer i, ii, j, jj, style, recv_ID, mm, nn, indices, send_count, topflag, send_count_ind, send_count_dat, send_count_tot, S_request1, S_request2, success
      character chara

      type(matrixblock) :: block
      type(matrixblock), pointer :: blocks_son
      integer statuss(MPI_STATUS_SIZE)
      integer send_tag

      send_count_ind = 0
      send_count_dat = 0
      call blocks_structure2buff(block, send_count_ind, send_count_dat, 0, msh, ptree)
      send_count_tot = send_count_ind + send_count_dat + 1
      call assert(send_count_dat /= 0, 'send_count_dat==0 not handled yet')
      call assert(send_count_tot < 1d9 .and. send_count_tot > 0, 'send_count_tot incorrect')
      allocate (ptree%send_buff_dat(send_count_tot), stat=success)
      call assert(success == 0, 'ptree%send_buff_dat allocation fail')
      ptree%send_buff_dat(1) = send_count_ind
      send_count_dat = send_count_ind + 1
      send_count_ind = 1
      call blocks_structure2buff(block, send_count_ind, send_count_dat, 1, msh, ptree)

      send_tag = ceiling_safe(send_count_tot/dble(msg_chunk))

      call MPI_Isend(ptree%send_buff_dat, send_count_tot, MPI_DT, recv_ID, send_tag, ptree%Comm, S_request2, ierr)

      if (option%verbosity >= 2) write (ptree%MyID + 10000, *) 'send: ', ptree%MyID, recv_ID, send_count_tot, block%row_group, block%col_group

      call MPI_wait(S_request2, statuss, ierr)
      if (ierr /= MPI_SUCCESS) then
         if (option%verbosity >= 2) write (ptree%MyID + 10000, *) 'send dat fail ', ierr
         stop
      end if

      if (option%verbosity >= 2) write (ptree%MyID + 10000, *) 'sendF: ', ptree%MyID, recv_ID, send_count_tot, block%row_group, block%col_group
      ! deallocate(send_buff_ind)
      deallocate (ptree%send_buff_dat)

      return

   end subroutine blocks_send

   recursive subroutine blocks_structure2buff(block, send_count_ind, send_count_dat, flag, msh, ptree)

      use BPACK_DEFS
      use MISC_Utilities
      implicit none

      integer count1, count2, requests, rank, group_m, group_n
      integer i, ii, j, jj, style, recv_ID, mm, nn, indices, send_count_ind, send_count_dat, flag
      character chara

      type(matrixblock) :: block
      type(matrixblock), pointer :: blocks_son
      type(proctree)::ptree
      type(mesh)::msh
      integer ierr

      if (flag == 1) ptree%send_buff_dat(send_count_ind + 1:send_count_ind + MPI_Header) = block%blockinfo_MPI
      send_count_ind = send_count_ind + MPI_Header
      call assert(send_count_ind < 1d9 .and. send_count_ind > 0, 'send_count_ind incorrect')

      style = block%style
      if (style == 4) then

         blocks_son => block%sons(1, 1)
         call blocks_structure2buff(blocks_son, send_count_ind, send_count_dat, flag, msh, ptree)
         blocks_son => block%sons(2, 1)
         call blocks_structure2buff(blocks_son, send_count_ind, send_count_dat, flag, msh, ptree)
         blocks_son => block%sons(1, 2)
         call blocks_structure2buff(blocks_son, send_count_ind, send_count_dat, flag, msh, ptree)
         blocks_son => block%sons(2, 2)
         call blocks_structure2buff(blocks_son, send_count_ind, send_count_dat, flag, msh, ptree)

      elseif (style == 2) then
         count1 = block%length_Butterfly_index_MPI
         count2 = block%length_Butterfly_data_MPI

         if (flag == 1) ptree%send_buff_dat(send_count_ind + 1:send_count_ind + count1) = block%Butterfly_index_MPI
         send_count_ind = send_count_ind + count1
         call assert(send_count_ind < 1d9 .and. send_count_ind > 0, 'send_count_ind incorrect')

         if (flag == 1) ptree%send_buff_dat(send_count_dat + 1:send_count_dat + count2) = block%Butterfly_data_MPI
         send_count_dat = send_count_dat + count2

         if (.not. (send_count_dat < 1d9 .and. send_count_dat > 0)) write (*, *) ptree%MyID, 'here1', send_count_dat, send_count_ind, count1, count2

         call assert(send_count_dat < 1d9 .and. send_count_dat > 0, 'send_count_dat incorrect')

      elseif (style == 1) then
         group_m = block%row_group
         group_n = block%col_group
         mm = (msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1)
         nn = (msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1)

         if (flag == 1) ptree%send_buff_dat(send_count_dat + 1:send_count_dat + mm*nn) = block%fullmat_MPI
         send_count_dat = send_count_dat + mm*nn

         if (.not. (send_count_dat < 1d9 .and. send_count_dat > 0)) write (*, *) ptree%MyID, 'here2', send_count_dat
         call assert(send_count_dat < 1d9 .and. send_count_dat > 0, 'send_count_dat incorrect')

         if (group_m == group_n) then
            if (flag == 1) ptree%send_buff_dat(send_count_ind + 1:send_count_ind + mm) = block%ipiv
            send_count_ind = send_count_ind + mm
            call assert(send_count_ind < 1d9 .and. send_count_ind > 0, 'send_count_ind incorrect')
         end if

      endif

      return

   end subroutine blocks_structure2buff

   subroutine blocks_recv(block, indices, send_ID, recv_count, msh, ptree, option)

      use BPACK_DEFS
      use MISC_Utilities
      implicit none

      integer blocks, flag_recv, count1, count2, recv_count, recv_count_loc, mm, nn, rank, mcnt
      logical rflag
      integer i, ii, j, jj, style, send_ID, group_m, group_n, indices, requestr, success, tag, recv_count_ind, recv_count_dat, recv_count_tot, R_request1, R_request2
      character chara
      type(proctree)::ptree
      type(mesh)::msh
      type(Hoption)::option
      integer ierr

      type(matrixblock) :: block
      type(matrixblock), pointer :: blocks_son
      integer statuss(MPI_STATUS_SIZE)

      recv_count_tot = 0
      mcnt = 0
      do while (mcnt /= 1)
         call MPI_IPROBE(send_ID, MPI_ANY_TAG, ptree%Comm, rflag, statuss, ierr)
         if (rflag .eqv. .true.) then
            mcnt = mcnt + 1
            tag = statuss(MPI_TAG)

            recv_count_tot = tag*msg_chunk
            allocate (ptree%recv_buff_dat(recv_count_tot), stat=success)
            call assert(success == 0, 'ptree%recv_buff_dat allocation fail')
            call MPI_Irecv(ptree%recv_buff_dat, recv_count_tot, MPI_DT, send_ID, tag, ptree%Comm, R_request2, ierr)

         end if
      enddo

      if (option%verbosity >= 2) write (ptree%MyID + 10000, *) 'recv: ', send_ID, ptree%MyID, recv_count_tot, block%row_group, block%col_group

      call MPI_wait(R_request2, statuss, ierr)
      if (ierr /= MPI_SUCCESS) then
         if (option%verbosity >= 2) write (ptree%MyID + 10000, *) 'recv dat fail ', ierr
         stop
      end if

      if (option%verbosity >= 2) write (ptree%MyID + 10000, *) 'recvF: ', send_ID, ptree%MyID, recv_count_tot, block%row_group, block%col_group

      recv_count_ind = 1
      recv_count_dat = 1 + NINT(dble(ptree%recv_buff_dat(1)))
      call blocks_buff2structure(block, recv_count_ind, recv_count_dat, msh, ptree)
      ! call assert(recv_count_dat==recv_count_tot,'recv_count_dat/=recv_count_tot')
      if (option%verbosity >= 2) write (ptree%MyID + 10000, *) 'copyF: ', send_ID, ptree%MyID, recv_count_dat, block%row_group, block%col_group

      ! deallocate(recv_buff_ind)
      deallocate (ptree%recv_buff_dat)

      return

   end subroutine blocks_recv

   recursive subroutine blocks_buff2structure(block, recv_count_ind, recv_count_dat, msh, ptree)

      use BPACK_DEFS
      use MISC_Utilities
      implicit none
      type(proctree)::ptree
      type(mesh)::msh
      integer ierr
      integer blocks, flag_recv, count1, count2, recv_count_ind, recv_count_dat, mm, nn, rank
      integer i, ii, j, jj, style, send_ID, group_m, group_n, indices, requestr, success
      character chara

      type(matrixblock):: block
      type(matrixblock), pointer :: blocks_son
      integer statuss(MPI_STATUS_SIZE)

      block%blockinfo_MPI = NINT(dble(ptree%recv_buff_dat(recv_count_ind + 1:recv_count_ind + MPI_Header)))
      recv_count_ind = recv_count_ind + MPI_Header

      block%level = block%blockinfo_MPI(1)
      block%row_group = block%blockinfo_MPI(2)
      block%col_group = block%blockinfo_MPI(3)
      ! block%nested_num=block%blockinfo_MPI(4)
      block%style = block%blockinfo_MPI(5)
      ! block%prestyle=block%blockinfo_MPI(6)
      ! block%data_type=block%blockinfo_MPI(7)
      block%level_butterfly = block%blockinfo_MPI(8)
      block%length_Butterfly_index_MPI = block%blockinfo_MPI(9)
      block%length_Butterfly_data_MPI = block%blockinfo_MPI(10)
      ! block%memory=dble(block%blockinfo_MPI(MPI_Header))

      style = block%style
      if (style == 4) then
         allocate (block%sons(2, 2), STAT=success)
         call assert(success == 0, 'butter allocation fail 3')
         !!! Be careful: commenting out this may cause problems

         ! do j=1,2
         ! do i=1,2
         ! block%sons(i,j)%father=>block
         ! enddo
         ! enddo

         blocks_son => block%sons(1, 1)
         call blocks_buff2structure(blocks_son, recv_count_ind, recv_count_dat, msh, ptree)
         blocks_son => block%sons(2, 1)
         call blocks_buff2structure(blocks_son, recv_count_ind, recv_count_dat, msh, ptree)
         blocks_son => block%sons(1, 2)
         call blocks_buff2structure(blocks_son, recv_count_ind, recv_count_dat, msh, ptree)
         blocks_son => block%sons(2, 2)
         call blocks_buff2structure(blocks_son, recv_count_ind, recv_count_dat, msh, ptree)

      elseif (style == 2) then
         count1 = block%length_Butterfly_index_MPI
         count2 = block%length_Butterfly_data_MPI
         call assert(count1 > 0 .and. count1 < 1d9, 'count1 incorrect')
         allocate (block%Butterfly_index_MPI(count1), STAT=success)
         call assert(success == 0, 'butter allocation fail 1')
         call assert(count2 > 0 .and. count2 < 1d9, 'count2 incorrect')
         allocate (block%Butterfly_data_MPI(count2), STAT=success)
         call assert(success == 0, 'butter allocation fail 2')

         block%Butterfly_index_MPI = NINT(dble(ptree%recv_buff_dat(recv_count_ind + 1:recv_count_ind + count1)))
         recv_count_ind = recv_count_ind + count1

         block%Butterfly_data_MPI = ptree%recv_buff_dat(recv_count_dat + 1:recv_count_dat + count2)
         recv_count_dat = recv_count_dat + count2

      elseif (style == 1) then
         group_m = block%row_group
         group_n = block%col_group
         mm = (msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1)
         nn = (msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1)
         !write(*,*)ptree%MyID, 'Irecv3',mm*nn,blocks,send_ID
         call assert(mm*nn > 0 .and. mm*nn < 1d9, 'mm*nn incorrect')
         allocate (block%fullmat_MPI(mm*nn), STAT=success)
         call assert(success == 0, 'fullmat allocation fail')

         block%fullmat_MPI = ptree%recv_buff_dat(recv_count_dat + 1:recv_count_dat + mm*nn)
         recv_count_dat = recv_count_dat + mm*nn

         if (group_m == group_n) then
            allocate (block%ipiv(mm), STAT=success)
            call assert(success == 0, 'fullmat allocation fail')
            block%ipiv = NINT(dble(ptree%recv_buff_dat(recv_count_ind + 1:recv_count_ind + mm)))
            recv_count_ind = recv_count_ind + mm
         end if

      endif

      return

   end subroutine blocks_buff2structure

   recursive subroutine Hmat_block_copy_MPIdata(block2, block1, msh)

      use BPACK_DEFS
      implicit none

      integer blocks, flag_recv, count1, count2, recv_count, mm, nn, length
      integer i, ii, j, jj, style, send_ID, group_m, group_n, indices, requests
      character chara

      type(matrixblock), pointer :: block1, block2, blocks_son1, blocks_son2
      type(mesh)::msh

      block2%blockinfo_MPI = block1%blockinfo_MPI

      block2%level = block2%blockinfo_MPI(1)
      ! write(*,*)block2%level,'ahaaa'
      block2%row_group = block2%blockinfo_MPI(2)
      block2%col_group = block2%blockinfo_MPI(3)
      ! block2%nested_num=block2%blockinfo_MPI(4)
      block2%style = block2%blockinfo_MPI(5)
      ! block2%prestyle=block2%blockinfo_MPI(6)
      ! block2%data_type=block2%blockinfo_MPI(7)
      block2%level_butterfly = block2%blockinfo_MPI(8)
      block2%length_Butterfly_index_MPI = block2%blockinfo_MPI(9)
      block2%length_Butterfly_data_MPI = block2%blockinfo_MPI(10)
      ! block2%memory=dble(block2%blockinfo_MPI(11))

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
         call Hmat_block_copy_MPIdata(blocks_son2, blocks_son1, msh)
         blocks_son1 => block1%sons(2, 1)
         blocks_son2 => block2%sons(2, 1)
         call Hmat_block_copy_MPIdata(blocks_son2, blocks_son1, msh)
         blocks_son1 => block1%sons(1, 2)
         blocks_son2 => block2%sons(1, 2)
         call Hmat_block_copy_MPIdata(blocks_son2, blocks_son1, msh)
         blocks_son1 => block1%sons(2, 2)
         blocks_son2 => block2%sons(2, 2)
         call Hmat_block_copy_MPIdata(blocks_son2, blocks_son1, msh)

      elseif (style == 2) then
         length = block2%length_Butterfly_index_MPI
         allocate (block2%Butterfly_index_MPI(length))
         block2%Butterfly_index_MPI = block1%Butterfly_index_MPI
         length = block2%length_Butterfly_data_MPI
         allocate (block2%Butterfly_data_MPI(length))
         block2%Butterfly_data_MPI = block1%Butterfly_data_MPI
      elseif (style == 1) then
         group_m = block2%row_group
         group_n = block2%col_group
         mm = (msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1)
         nn = (msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1)
         allocate (block2%fullmat_MPI(mm*nn))
         block2%fullmat_MPI = block1%fullmat_MPI
         if (group_m == group_n) then
            allocate (block2%ipiv(mm))
            block2%ipiv = block1%ipiv
         endif
      endif

      return

   end subroutine Hmat_block_copy_MPIdata

   subroutine pack_butterfly_blocks(block, msh)

      use BPACK_DEFS
      implicit none

      integer i, j, k, ii, jj, kk, num_blocks, level_butterfly, level_blocks, level
      integer count1, count2, num_row, num_col, indices, mm, nn
      type(matrixblock), pointer :: block
      type(mesh) :: msh

      level_butterfly = block%level_butterfly
      level_blocks = block%level

      if (level_butterfly == 0) then
         num_blocks = 1
      else
         num_blocks = 2**level_butterfly
      endif
      count1 = 4 + num_blocks*2*2
      count2 = 0
      do i = 1, num_blocks
         mm = size(block%ButterflyU%blocks(i)%matrix, 1)
         nn = size(block%ButterflyU%blocks(i)%matrix, 2)
         count2 = count2 + mm*nn
         mm = size(block%ButterflyV%blocks(i)%matrix, 1)
         nn = size(block%ButterflyV%blocks(i)%matrix, 2)
         count2 = count2 + mm*nn
      enddo
      if (level_butterfly /= 0) then
         do level = 1, level_butterfly
            num_row = block%ButterflyKerl(level)%num_row
            num_col = block%ButterflyKerl(level)%num_col
            count1 = count1 + 2
            do i = 1, num_row
               do j = 1, num_col
                  mm = size(block%ButterflyKerl(level)%blocks(i, j)%matrix, 1)
                  nn = size(block%ButterflyKerl(level)%blocks(i, j)%matrix, 2)
                  count1 = count1 + 2
                  count2 = count2 + mm*nn
               enddo
            enddo
         enddo
      endif

      block%length_Butterfly_index_MPI = count1
      block%length_Butterfly_data_MPI = count2
      block%blockinfo_MPI(9) = count1
      block%blockinfo_MPI(10) = count2
      allocate (block%Butterfly_index_MPI(count1))
      allocate (block%Butterfly_data_MPI(count2))

      !block%Butterfly_index_MPI(1)=blocks
      level_blocks = block%level
      block%Butterfly_index_MPI(2) = block%rankmax
      level_butterfly = block%level_butterfly
      block%Butterfly_index_MPI(3) = level_butterfly
      block%Butterfly_index_MPI(4) = num_blocks
      count1 = INDEX_Header
      count2 = 0

      do i = 1, num_blocks
         count1 = count1 + 1
         mm = size(block%ButterflyU%blocks(i)%matrix, 1)
         block%Butterfly_index_MPI(count1) = mm
         count1 = count1 + 1
         nn = size(block%ButterflyU%blocks(i)%matrix, 2)
         block%Butterfly_index_MPI(count1) = nn
         !$omp parallel do default(shared) private(ii,jj,indices)
         do jj = 1, nn
            do ii = 1, mm
               indices = (jj - 1)*mm + ii
               block%Butterfly_data_MPI(count2 + indices) = block%ButterflyU%blocks(i)%matrix(ii, jj)
            enddo
         enddo
         !$omp end parallel do
         count2 = count2 + mm*nn
      enddo
      do i = 1, num_blocks
         count1 = count1 + 1
         mm = size(block%ButterflyV%blocks(i)%matrix, 1)
         block%Butterfly_index_MPI(count1) = mm
         count1 = count1 + 1
         nn = size(block%ButterflyV%blocks(i)%matrix, 2)
         block%Butterfly_index_MPI(count1) = nn
         !$omp parallel do default(shared) private(ii,jj,indices)
         do jj = 1, nn
            do ii = 1, mm
               indices = (jj - 1)*mm + ii
               block%Butterfly_data_MPI(count2 + indices) = block%ButterflyV%blocks(i)%matrix(ii, jj)
            enddo
         enddo
         !$omp end parallel do
         count2 = count2 + mm*nn
      enddo

      if (level_butterfly /= 0) then
         do level = 1, level_butterfly
            num_row = block%ButterflyKerl(level)%num_row
            count1 = count1 + 1
            block%Butterfly_index_MPI(count1) = num_row
            num_col = block%ButterflyKerl(level)%num_col
            count1 = count1 + 1
            block%Butterfly_index_MPI(count1) = num_col
            do i = 1, num_row
               do j = 1, num_col
                  mm = size(block%ButterflyKerl(level)%blocks(i, j)%matrix, 1)
                  count1 = count1 + 1
                  block%Butterfly_index_MPI(count1) = mm
                  nn = size(block%ButterflyKerl(level)%blocks(i, j)%matrix, 2)
                  count1 = count1 + 1
                  block%Butterfly_index_MPI(count1) = nn
                  !$omp parallel do default(shared) private(ii,jj,indices)
                  do jj = 1, nn
                     do ii = 1, mm
                        indices = (jj - 1)*mm + ii
                        block%Butterfly_data_MPI(count2 + indices) = block%ButterflyKerl(level)%blocks(i, j)%matrix(ii, jj)
                     enddo
                  enddo
                  !$omp end parallel do
                  count2 = count2 + mm*nn
               enddo
            enddo
         enddo
      endif

      !$omp parallel do default(shared) private(i)
      do i = 1, num_blocks
         deallocate (block%ButterflyU%blocks(i)%matrix)
         deallocate (block%ButterflyV%blocks(i)%matrix)
      enddo
      !$omp end parallel do
      deallocate (block%ButterflyU%blocks)
      deallocate (block%ButterflyV%blocks)

      if (level_butterfly /= 0) then
         !$omp parallel do default(shared) private(level,i,j,num_col,num_row)
         do level = 1, level_butterfly
            num_col = block%ButterflyKerl(level)%num_col
            num_row = block%ButterflyKerl(level)%num_row
            do j = 1, num_col
               do i = 1, num_row
                  deallocate (block%ButterflyKerl(level)%blocks(i, j)%matrix)
               enddo
            enddo
            deallocate (block%ButterflyKerl(level)%blocks)
         enddo
         !$omp end parallel do
         deallocate (block%ButterflyKerl)
      endif

      return

   end subroutine pack_butterfly_blocks

   subroutine unpack_butterfly_blocks(block, Maxlevel, ptree, msh, pgno)

      use BPACK_DEFS
      implicit none

      integer i, j, k, ii, jj, kk, num_blocks, level_butterfly, level_blocks, level
      integer count1, count2, num_row, num_col, indices, mm, nn
      type(matrixblock), pointer :: block
      integer Maxlevel
      type(proctree)::ptree
      type(mesh)::msh
      integer pgno

      block%rankmax = block%Butterfly_index_MPI(2)
      level_butterfly = block%Butterfly_index_MPI(3)
      block%level_butterfly = level_butterfly
      num_blocks = block%Butterfly_index_MPI(4)
      count1 = 4
      count2 = 0

      call BF_Init_blocks(level_butterfly, block%row_group, block%col_group, pgno, block, msh, ptree)

      allocate (block%ButterflyU%blocks(num_blocks))
      do i = 1, num_blocks
         count1 = count1 + 1
         mm = block%Butterfly_index_MPI(count1)
         count1 = count1 + 1
         nn = block%Butterfly_index_MPI(count1)
         allocate (block%ButterflyU%blocks(i)%matrix(mm, nn))
         !$omp parallel do default(shared) private(ii,jj,indices)
         do jj = 1, nn
            do ii = 1, mm
               indices = (jj - 1)*mm + ii
               block%ButterflyU%blocks(i)%matrix(ii, jj) = block%Butterfly_data_MPI(count2 + indices)
            enddo
         enddo
         !$omp end parallel do
         count2 = count2 + mm*nn
      enddo
      allocate (block%ButterflyV%blocks(num_blocks))
      do i = 1, num_blocks
         count1 = count1 + 1
         mm = block%Butterfly_index_MPI(count1)
         count1 = count1 + 1
         nn = block%Butterfly_index_MPI(count1)
         allocate (block%ButterflyV%blocks(i)%matrix(mm, nn))
         !$omp parallel do default(shared) private(ii,jj,indices)
         do jj = 1, nn
            do ii = 1, mm
               indices = (jj - 1)*mm + ii
               block%ButterflyV%blocks(i)%matrix(ii, jj) = block%Butterfly_data_MPI(count2 + indices)
            enddo
         enddo
         !$omp end parallel do
         count2 = count2 + mm*nn
      enddo

      if (level_butterfly /= 0) then
         if (.not. allocated(block%ButterflyKerl)) allocate (block%ButterflyKerl(level_butterfly))
         do level = 1, level_butterfly
            count1 = count1 + 1
            num_row = block%Butterfly_index_MPI(count1)
            block%ButterflyKerl(level)%num_row = num_row
            count1 = count1 + 1
            num_col = block%Butterfly_index_MPI(count1)
            block%ButterflyKerl(level)%num_col = num_col
            allocate (block%ButterflyKerl(level)%blocks(num_row, num_col))
            do i = 1, num_row
               do j = 1, num_col
                  count1 = count1 + 1
                  mm = block%Butterfly_index_MPI(count1)
                  count1 = count1 + 1
                  nn = block%Butterfly_index_MPI(count1)
                  allocate (block%ButterflyKerl(level)%blocks(i, j)%matrix(mm, nn))
                  !$omp parallel do default(shared) private(ii,jj,indices)
                  do jj = 1, nn
                     do ii = 1, mm
                        indices = (jj - 1)*mm + ii
                        block%ButterflyKerl(level)%blocks(i, j)%matrix(ii, jj) = block%Butterfly_data_MPI(count2 + indices)
                     enddo
                  enddo
                  !$omp end parallel do
                  count2 = count2 + mm*nn
               enddo
            enddo
         enddo
      endif

      deallocate (block%Butterfly_index_MPI)
      deallocate (block%Butterfly_data_MPI)
      return

   end subroutine unpack_butterfly_blocks

   subroutine pack_full_blocks(block, msh)

      use BPACK_DEFS
      implicit none

      integer blocks, i, j, k, ii, jj, kk, group_m, group_n, mm, nn, flag
      integer count1, count2, num_row, num_col, indices
      type(matrixblock), pointer :: block
      type(mesh)::msh
      real(kind=8)::tol_used

      group_m = block%row_group
      group_n = block%col_group
      mm = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
      nn = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
      allocate (block%fullmat_MPI(mm*nn))

#if HAVE_ZFP
      call ZFP_Decompress(block,tol_used) ! no need to recompress as fullmat will be deleted before exiting
#endif

      !$omp parallel do default(shared) private(i,j,indices)
      do j = 1, nn
         do i = 1, mm
            indices = mm*(j - 1) + i
            block%fullmat_MPI(indices) = block%fullmat(i, j)
         enddo
      enddo
      !$omp end parallel do
      deallocate (block%fullmat)

      return

   end subroutine pack_full_blocks

   subroutine unpack_full_blocks(block, Maxlevel, ptree, msh, pgno, option)

      use BPACK_DEFS
      implicit none

      integer blocks, i, j, k, ii, jj, kk, group_m, group_n, mm, nn, flag
      integer count1, count2, num_row, num_col, indices
      type(matrixblock), pointer :: block
      integer Maxlevel
      type(proctree)::ptree
      type(Hoption)::option
      type(mesh)::msh
      integer pgno

      group_m = block%row_group
      group_n = block%col_group
      mm = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
      nn = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
      allocate (block%fullmat(mm, nn))
      !$omp parallel do default(shared) private(i,j,indices)
      do j = 1, nn
         do i = 1, mm
            indices = mm*(j - 1) + i
            block%fullmat(i, j) = block%fullmat_MPI(indices)
         enddo
      enddo
      !$omp end parallel do
      deallocate (block%fullmat_MPI)
#if HAVE_ZFP
      call ZFP_Compress(block,option%tol_comp)
#endif
      return

   end subroutine unpack_full_blocks

   recursive subroutine pack_all_blocks_one_node(block, msh)

      use BPACK_DEFS
      implicit none

      integer level_blocks
      integer i, ii, j, jj, k, kk, level, style
      type(matrixblock), pointer :: block, blocks_son
      type(mesh)::msh

      block%blockinfo_MPI(1) = block%level
      block%blockinfo_MPI(2) = block%row_group
      block%blockinfo_MPI(3) = block%col_group
      ! block%blockinfo_MPI(4)=block%nested_num
      block%blockinfo_MPI(5) = block%style
      ! block%blockinfo_MPI(6)=block%prestyle
      ! block%blockinfo_MPI(7)=block%data_type
      block%blockinfo_MPI(8) = block%level_butterfly
      block%blockinfo_MPI(9) = 0
      block%blockinfo_MPI(10) = 0
      ! block%blockinfo_MPI(11)=int(block%memory)

      if (block%style == 4) then

         blocks_son => block%sons(1, 1)
         call pack_all_blocks_one_node(blocks_son, msh)
         blocks_son => block%sons(2, 1)
         call pack_all_blocks_one_node(blocks_son, msh)
         blocks_son => block%sons(1, 2)
         call pack_all_blocks_one_node(blocks_son, msh)
         blocks_son => block%sons(2, 2)
         call pack_all_blocks_one_node(blocks_son, msh)

      elseif (block%style == 2) then

         call pack_butterfly_blocks(block, msh)

      elseif (block%style == 1) then

         call pack_full_blocks(block, msh)

      endif

      return

   end subroutine pack_all_blocks_one_node

   recursive subroutine unpack_all_blocks_one_node(block, Maxlevel, ptree, msh, pgno,option)

      use BPACK_DEFS
      implicit none

      integer level_blocks
      integer i, ii, j, jj, k, kk, level, style
      type(matrixblock), pointer :: block, blocks_son
      integer Maxlevel
      type(proctree)::ptree
      type(Hoption)::option
      type(mesh)::msh
      integer pgno, pgno1, pgno2, Maxgrp

      call assert(block%row_group > 0 .and. block%row_group <= msh%Maxgroup .and. block%col_group > 0 .and. block%col_group <= msh%Maxgroup, 'wrong row_group or col_group in unpack_all_blocks_one_node')
      block%pgno = pgno
      block%M = msh%basis_group(block%row_group)%tail - msh%basis_group(block%row_group)%head + 1
      block%N = msh%basis_group(block%col_group)%tail - msh%basis_group(block%col_group)%head + 1
      block%headm = msh%basis_group(block%row_group)%head
      block%headn = msh%basis_group(block%col_group)%head
      call ComputeParallelIndices(block, block%pgno, ptree, msh)

      if (block%style == 4) then
         Maxgrp = 2**(ptree%nlevel) - 1
         if (pgno*2 <= Maxgrp) then
            pgno1 = pgno*2
         else
            pgno1 = pgno
         endif
         if (pgno*2 + 1 <= Maxgrp) then
            pgno2 = pgno*2 + 1
         else
            pgno2 = pgno
         endif

         blocks_son => block%sons(1, 1)
         call unpack_all_blocks_one_node(blocks_son, Maxlevel, ptree, msh, pgno1, option)
         blocks_son => block%sons(2, 1)
         call unpack_all_blocks_one_node(blocks_son, Maxlevel, ptree, msh, pgno2, option)
         blocks_son => block%sons(1, 2)
         call unpack_all_blocks_one_node(blocks_son, Maxlevel, ptree, msh, pgno1, option)
         blocks_son => block%sons(2, 2)
         call unpack_all_blocks_one_node(blocks_son, Maxlevel, ptree, msh, pgno2, option)

      elseif (block%style == 2) then
         call unpack_butterfly_blocks(block, Maxlevel, ptree, msh, pgno)

      elseif (block%style == 1) then

         call unpack_full_blocks(block, Maxlevel, ptree, msh, pgno, option)

      endif

      return

   end subroutine unpack_all_blocks_one_node

end module BPACK_block_sendrecv

