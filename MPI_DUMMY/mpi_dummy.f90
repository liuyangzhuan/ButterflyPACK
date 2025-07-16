subroutine MPI_COMM_GET_PARENT(PARENT, IERROR)
  implicit none
  include "mpi_dummy.fi"
  INTEGER    PARENT, IERROR
  ierror = MPI_SUCCESS
  PARENT = MPI_COMM_NULL+1 ! return a non-MPI_COMM_NULL number
end subroutine MPI_COMM_GET_PARENT

subroutine MPI_Comm_disconnect(COMM, IERROR)
  implicit none
  include "mpi_dummy.fi"
  INTEGER    COMM, IERROR
  ierror = MPI_SUCCESS
end subroutine MPI_Comm_disconnect

subroutine mpi_allgather ( sendbuf, sendcount, sendtype, recvbuf, recvcount, &
  recvtype, comm, ierror )


!*****************************************************************************80
!
!! MPI_ALLGATHER gathers data from all the processes in a communicator.
!
!  Discussion:
!
!    The block of data sent from the J-th process is placed in the
!    J-th block of the receive buffer of every process.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type,
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, SENDTYPE SENDBUF(SENDCOUNT), the data to be sent.
!
!    Inout, integer SENDCOUNT, the number of data items being sent.
!
!    Input, integer SENDTYPE, the datatype of the data being sent.
!
!    Output, RECVTYPE RECVBUF(RECVCOUNT,GROUPSIZE), the data as received.
!
!    Input, integer RECVCOUNT, the number of data items to be received
!    from any one process.
!
!    Input, integer RECVTYPE, the datatype of the data being received.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!

  implicit none

  include "mpi_dummy.fi"

  integer recvcount
  integer sendcount

  integer comm
  integer ierror
  integer recvbuf(recvcount,*)
  integer recvtype
  integer sendbuf(sendcount)
  integer sendtype

  ierror = MPI_SUCCESS

  if ( sendtype == mpi_double_precision ) then
    call mpi_copy_double_precision ( sendbuf, recvbuf(:,1), sendcount, ierror )
  else if ( sendtype == mpi_integer ) then
    call mpi_copy_integer ( sendbuf, recvbuf(:,1), sendcount, ierror )
  else if ( sendtype == mpi_float ) then
    call mpi_copy_real ( sendbuf, recvbuf(:,1), sendcount, ierror )
  else if ( sendtype == mpi_complex ) then
    call mpi_copy_complex ( sendbuf, recvbuf(:,1), sendcount, ierror )
  else if ( sendtype == mpi_double_complex ) then
    call mpi_copy_double_complex ( sendbuf, recvbuf(:,1), sendcount, ierror )
  else
    ierror = MPI_FAILURE
  end if

  return
end subroutine mpi_allgather


subroutine mpi_allgatherv ( sendbuf, sendcount, sendtype, recvbuf, &
  recvcounts, displs, recvtype, comm, ierror )

!*****************************************************************************80
!
!! MPI_ALLGATHERV gathers data from all the processes in a communicator.
!
!  Discussion:
!
!    Copy values from DATA1 to DATA2.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type,
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, SENDTYPE SENDBUF(SENDCOUNT), the data to be sent.
!
!    Inout, integer SENDCOUNT, the number of data items being sent.
!
!    Input, integer SENDTYPE, the datatype of the data being sent.
!
!    Output, RECVTYPE RECVBUF(*), the buffer to store the received data.
!
!    Input, integer RECVCOUNTS(0:GROUP_SIZE-1), the number of data items to be 
!    received from each process.
!
!    Input, integer DISPLS(0:GROUP_SIZE-1), the I-th element is the displacement
!    in RECVBUF at which to place the input from process I.
!
!    Input, integer RECVTYPE, the datatype of the data being received.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  include "mpi_dummy.fi"

  integer sendcount

  integer comm
  integer displs(0:*)
  integer ierror
  integer nrecv
  integer recvbuf(*)
  integer recvcounts(0:*)
  integer recvtype
  integer sendbuf(sendcount)
  integer sendtype

  ierror = MPI_SUCCESS

  if ( sendtype == mpi_double_precision ) then
    call mpi_copy_double_precision ( &
      sendbuf, recvbuf(displs(0)), recvcounts(0), ierror )
  else if ( sendtype == mpi_integer ) then
    call mpi_copy_integer ( &
      sendbuf, recvbuf(displs(0)), recvcounts(0), ierror )
  else if ( sendtype == mpi_float ) then
    call mpi_copy_real ( &
      sendbuf, recvbuf(displs(0)), recvcounts(0), ierror )
  else if ( sendtype == mpi_complex ) then
    call mpi_copy_complex ( &
      sendbuf, recvbuf(displs(0)), recvcounts(0), ierror )
  else if ( sendtype == mpi_double_complex ) then
    call mpi_copy_double_complex ( &
      sendbuf, recvbuf(displs(0)), recvcounts(0), ierror )
  else
    ierror = MPI_FAILURE
  end if

  return
end



subroutine mpi_allreduce ( data1, data2, n, datatype, operation, comm, ierror )

!*****************************************************************************80
!
!! MPI_ALLREDUCE carries out a reduction operation.
!
!  Discussion:
!
!    The reduction operations are MAXIMUM, MINIMUM, PRODUCT and SUM.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type,
!    but using the other types should generally not cause a problem.
!
!    Thanks to Simppa Akaslompolo for correcting this routine!
!    12 January 2012.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, DATATYPE DATA1(N), the data to be processed.
!
!    Output, DATATYPE DATA2(N), the value of the reduction operation.
!
!    Input, integer N, the number of items in DATA1.
!
!    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
!
!    Input, integer OPERATION, should have the value of one of the symbolic
!    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  include "mpi_dummy.fi"

  integer n

  integer comm
  integer data1(*)
  integer data2(*)
  integer datatype
  integer ierror
  integer operation


  ierror = MPI_SUCCESS

  if(MPI_IN_PLACE==data1(1) .or. MPI_IN_PLACE==data2(1) .or. MPI_BOTTOM==data1(1) .or. MPI_BOTTOM==data2(1))return


  if ( datatype == mpi_double_precision ) then

    call mpi_reduce_double_precision ( data1, data2, n, operation, ierror )

  else if ( datatype == mpi_integer ) then

    call mpi_reduce_integer ( data1, data2, n, operation, ierror )

  else if ( datatype == mpi_float ) then

    call mpi_reduce_real ( data1, data2, n, operation, ierror )
  else if ( datatype == mpi_complex ) then

    call mpi_reduce_complex ( data1, data2, n, operation, ierror )
  else if ( datatype == mpi_double_complex ) then

    call mpi_reduce_double_complex ( data1, data2, n, operation, ierror )    
  else

    ierror = MPI_FAILURE

  end if

  return
end


subroutine MPI_ALLTOALL(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT,RECVTYPE, COMM, IERROR)
  implicit none

  include "mpi_dummy.fi"
    
    
    Integer    SENDBUF(*), RECVBUF(*)
    INTEGER    SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE
    INTEGER    COMM, IERROR
    if(RECVTYPE/=SENDTYPE)then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPI_ALLTOALL - Error!'
      write ( *, '(a)' )  '  RECVTYPE needs to match SENDTYPE in mpi_dummy. '
    endif 
    
  if ( sendtype == mpi_double_precision ) then
    call mpi_copy_double_precision ( sendbuf, recvbuf, sendcount, ierror )
  else if ( sendtype == mpi_integer ) then
    call mpi_copy_integer ( sendbuf, recvbuf, sendcount, ierror )
  else if ( sendtype == mpi_float ) then
    call mpi_copy_real ( sendbuf, recvbuf, sendcount, ierror )
  else if ( sendtype == mpi_complex ) then
    call mpi_copy_complex ( sendbuf, recvbuf, sendcount, ierror )
  else if ( sendtype == mpi_double_complex ) then
    call mpi_copy_double_complex ( sendbuf, recvbuf, sendcount, ierror )
  else
    ierror = MPI_FAILURE
  end if


end subroutine MPI_ALLTOALL



subroutine mpi_barrier ( comm, ierror )
  implicit none

  include "mpi_dummy.fi"

  integer comm
  integer ierror

  ierror = MPI_SUCCESS

  return
end subroutine mpi_barrier

subroutine mpi_bcast ( data, n, datatype, node, comm, ierror )

  implicit none

  include "mpi_dummy.fi"

  integer n

  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer node

  ierror = MPI_SUCCESS

  return 
end subroutine mpi_bcast

subroutine MPI_IPROBE(source, tag, comm, flag, statuss, ierror)
  implicit none

  include "mpi_dummy.fi"
  logical flag
  integer comm, source, tag
  integer statuss(MPI_STATUS_SIZE)
  integer ierror
  flag = .true.
  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_IPROBE - Error!'
  write ( *, '(a)' ) '  Should not query message from self.'


  return
end subroutine MPI_IPROBE

subroutine MPI_PROBE(source, tag, comm, statuss, ierror)
  implicit none

  include "mpi_dummy.fi"
  integer comm, source, tag
  integer statuss(MPI_STATUS_SIZE)
  integer ierror
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_PROBE - Error!'
  write ( *, '(a)' ) '  Should not query message from self.'


  return
end subroutine MPI_PROBE


subroutine MPI_Comm_Create(comm_base, group, comm,ierror)
  implicit none

  include "mpi_dummy.fi"

  integer comm, group, comm_base
  integer ierror
  if(comm_base /=MPI_COMM_NULL .and. group /= MPI_GROUP_NULL)then
    comm = MPI_COMM_NULL+ 1
  else 
    comm = MPI_COMM_NULL
  endif
  ierror = MPI_SUCCESS

  return
end subroutine MPI_Comm_Create


subroutine MPI_Group_incl(group, n, ranks, newgroup, ierror)
    implicit none
    include "mpi_dummy.fi"    
    
    integer :: group
    integer :: n, ranks(n)
    integer :: newgroup
    integer :: ierror

    if(group /= MPI_GROUP_NULL)then
      newgroup = MPI_GROUP_NULL+1
    endif
    ierror = MPI_SUCCESS


end subroutine MPI_Group_incl


subroutine MPI_Comm_group(comm, group, ierror)
  implicit none

  include "mpi_dummy.fi"

  integer comm, group
  integer ierror
  if(comm /=MPI_COMM_NULL)then
    group = MPI_GROUP_NULL+1
  else 
    group = MPI_GROUP_NULL
  endif
  ierror = MPI_SUCCESS

  return
end subroutine MPI_Comm_group

subroutine mpi_group_free ( group, ierror )

  implicit none

  include "mpi_dummy.fi"

  integer group
  integer ierror

  ierror = MPI_SUCCESS

  return
end subroutine mpi_group_free


subroutine mpi_comm_free ( comm, ierror )
  implicit none

  include "mpi_dummy.fi"

  integer comm
  integer ierror

  ierror = MPI_SUCCESS

  return
end subroutine mpi_comm_free




subroutine mpi_comm_rank ( comm, me, ierror )

  implicit none

  include "mpi_dummy.fi"

  integer comm
  integer ierror
  integer me

  ierror = MPI_SUCCESS
  me = 0

  return
end subroutine mpi_comm_rank

subroutine mpi_comm_size ( comm, nprocs, ierror )

  implicit none

  include "mpi_dummy.fi"

  integer comm
  integer ierror
  integer nprocs

  ierror = MPI_SUCCESS
  nprocs = 1

  return
end subroutine mpi_comm_size

subroutine mpi_comm_split ( comm, icolor, ikey, comm_new, ierror )
  implicit none

  include "mpi_dummy.fi"

  integer comm
  integer comm_new
  integer icolor
  integer ierror
  integer ikey

  ierror = MPI_SUCCESS

  return
end subroutine mpi_comm_split

subroutine mpi_copy_double_precision ( data1, data2, n, ierror )

!*****************************************************************************80
!
!! MPI_COPY_DOUBLE copies a double precision vector.
!
!  Discussion:
!
!    This routine is not part of the MPI standard.  However, it is
!    needed by other routines which do emulate standard MPI routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, double precision DATA1(N), the data to be copied.
!
!    Output, double precision DATA2(N), the copied data.
!
!    Input, integer N, the number of items of data.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  include "mpi_dummy.fi"

  integer n

  double precision data1(n)
  double precision data2(n)
  integer ierror

  ierror = MPI_SUCCESS

  data2(1:n) = data1(1:n)

  return
end

subroutine mpi_copy_integer ( data1, data2, n, ierror )

!*****************************************************************************80
!
!! MPI_COPY_INTEGER copies an integer vector.
!
!  Discussion:
!
!    This routine is not part of the MPI standard.  However, it is
!    needed by other routines which do emulate standard MPI routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer DATA1(N), the data to be copied.
!
!    Output, integer DATA2(N), the copied data.
!
!    Input, integer N, the number of items of data.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  include "mpi_dummy.fi"

  integer n

  integer data1(n)
  integer data2(n)
  integer ierror

  ierror = MPI_SUCCESS

  data2(1:n) = data1(1:n)

  return
end

subroutine mpi_copy_real ( data1, data2, n, ierror )

!*****************************************************************************80
!
!! MPI_COPY_REAL copies a real vector.
!
!  Discussion:
!
!    This routine is not part of the MPI standard.  However, it is
!    needed by other routines which do emulate standard MPI routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, real DATA1(N), the data to be copied.
!
!    Output, real DATA2(N), the copied data.
!
!    Input, integer N, the number of items of data.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  include "mpi_dummy.fi"

  integer n

  real data1(n)
  real data2(n)
  integer ierror

  ierror = MPI_SUCCESS

  data2(1:n) = data1(1:n)

  return
end


subroutine mpi_copy_complex ( data1, data2, n, ierror )
  implicit none

  include "mpi_dummy.fi"

  integer n

  complex(kind=4) data1(n)
  complex(kind=4) data2(n)
  integer ierror

  ierror = MPI_SUCCESS

  data2(1:n) = data1(1:n)

  return
end subroutine mpi_copy_complex


subroutine mpi_copy_double_complex ( data1, data2, n, ierror )
  implicit none

  include "mpi_dummy.fi"

  integer n

  complex(kind=8) data1(n)
  complex(kind=8) data2(n)
  integer ierror

  ierror = MPI_SUCCESS

  data2(1:n) = data1(1:n)

  return
end subroutine mpi_copy_double_complex



subroutine mpi_finalize ( ierror )
  implicit none

  include "mpi_dummy.fi"

  integer ierror

  ierror = MPI_SUCCESS

  return
end


subroutine mpi_get_count ( status, datatype, icount, ierror )


!*****************************************************************************80
!
!! MPI_GET_COUNT reports the number of items actually transmitted.
!
!  Discussion:
!
!    Warn against querying message from self, since no data copy is done.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer STATUS(MPI_STATUS_SIZE), the MPI message status.
!
!    Input, integer DATATYPE, the datatype of the data.
!
!    Output, integer ICOUNT, the number of data items transmitted.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  include "mpi_dummy.fi"

  integer datatype
  integer icount
  integer ierror
  integer status(mpi_status_size)

  ierror = MPI_FAILURE

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_GET_COUNT - Error!'
  write ( *, '(a)' ) '  Should not query message from self.'

  return
end

subroutine mpi_init ( ierror )

  implicit none

  include "mpi_dummy.fi"

  integer ierror

  ierror = MPI_SUCCESS

  return
end subroutine mpi_init

subroutine mpi_init_thread (required, provided, ierror )

  implicit none

  include "mpi_dummy.fi"

  integer required, provided, ierror

  ierror = MPI_SUCCESS

  return
end subroutine mpi_init_thread


subroutine mpi_irecv ( data, n, datatype, iproc, itag, comm, irequest, ierror )

!*****************************************************************************80
!
!! MPI_IRECV performs an immediate receive of data from another process.
!
!  Discussion:
!
!    For an immediate or nonblocking receive, the call to mpi_irecv begins
!    a receive operation, but program execution may proceed to the next
!    statement without waiting for confirmation that the receive has 
!    been completed.  It is up to the user to issue an appropriate
!    statement later, such as a call to MPI_WAIT, with a copy of the 
!    value of IREQUEST, to verify that the receive has completed.
!
!    Warn against receiving message from self, since no data copy is done.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type,
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, DATATYPE DATA(N), a buffer which will contain the data.
!
!    Input, integer N, the number of items expected.
!
!    Input, integer DATATYPE, the MPI datatype of the data.
!
!    Input, integer IPROC, the MPI process from which the data is to
!    be received.
!
!    Input, integer ITAG, a tag for the message.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IREQUEST, the request handle.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  include "mpi_dummy.fi"

  integer n

  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer irequest
  integer itag

  ierror = MPI_FAILURE

  write ( *, '(/,a)' ) 'MPI_IRECV - Error: Should not "recv" message from self.'

  return
end

subroutine mpi_isend ( data, n, datatype, iproc, itag, comm, request, ierror )

!*****************************************************************************80
!
!! MPI_ISEND sends data to another process using nonblocking transmission.
!
!  Discussion:
!
!    Warn against sending message to self, since no data copy is done.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type,
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, datatype DATA(N), the data to be sent.
!
!    Input, integer N, the number of data items to send.
!
!    Input, integer DATAYTPE, the MPI code for the datatype.
!
!    Input, integer IPROC, the rank of the process within the communicator
!    that is to receive the message.
!
!    Input, integer ITAG, a tag for the message.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer REQUEST, a handle.  To determine if the data has been 
!    received yet, call MPI_Test or MPI_Wait, including the value of REQUEST.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  include "mpi_dummy.fi"

  integer n

  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer itag
  integer request

  request = 0
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MPI_ISEND - Error!'
  write ( *, '(a)' )  '  Should not send message to self.'

  return
end

subroutine mpi_recv ( data, n, datatype, iproc, itag, comm, status, ierror )

!*****************************************************************************80
!
!! MPI_RECV receives data from another process within a communicator.
!
!  Discussion:
!
!    Warn against receiving message from self, since no data copy is done.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type,
!    but using the other types should generally not cause a problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, DATATYPE DATA(N), a buffer which will contain the data.
!
!    Input, integer N, the number of items expected.
!
!    Input, integer DATATYPE, the MPI datatype of the data.
!
!    Input, integer IPROC, the MPI process from which the data is to
!    be received.
!
!    Input, integer ITAG, a tag for the message.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer STATUS(MPI_STATUS_SIZE), the MPI message status.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  include "mpi_dummy.fi"

  integer n

  integer comm
  integer data(n)
  integer datatype
  integer ierror
  integer iproc
  integer itag
  integer status(mpi_status_size)

  ierror = MPI_FAILURE

  write ( *, '(/,a)' ) 'MPI_RECV - Error: Should not "recv" message from self.'

  return
end
subroutine mpi_reduce ( data1, data2, n, datatype, operation, receiver, &
  comm, ierror )

!*****************************************************************************80
!
!! MPI_REDUCE carries out a reduction operation.
!
!  Discussion:
!
!    The reduction operations are sum, maximum, minimum, product.
!
!    The first two arguments must not overlap or share memory in any way.
!
!    The data to be transferred can be integer, real, or double precision.
!    In this routine, it is declared and documented as INTEGER type,
!    but using the other types should generally not cause a problem.
!
!    Thanks to Simppa Akaslompolo for correcting this routine!
!    11 January 2012.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, DATATYPE DATA1(N), the data to be processed.
!
!    Output (to RECEIVER only), DATATYPE DATA2(N), the value of the
!    reduction operation.
!
!    Input, integer N, the number of items in DATA1.
!
!    Input, integer DATATYPE, indicates the datatype of DATA1 and DATA2.
!
!    Input, integer OPERATION, should have the value of one of the symbolic
!    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
!
!    Input, integer RECEIVER, the process that is to receive the
!    result.
!
!    Input, integer COMM, the MPI communicator.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  include "mpi_dummy.fi"

  integer n

  integer comm
  integer data1(*)
  integer data2(*)
  integer datatype
  integer ierror
  integer operation
  integer receiver

  if(MPI_IN_PLACE==data1(1) .or. MPI_IN_PLACE==data2(1) .or. MPI_BOTTOM==data1(1) .or. MPI_BOTTOM==data2(1))then
    ierror = MPI_SUCCESS
    return
  endif

  if ( datatype == mpi_double_precision ) then

    call mpi_reduce_double_precision ( data1, data2, n, operation, ierror )

  else if ( datatype == mpi_integer ) then

    call mpi_reduce_integer ( data1, data2, n, operation, ierror )

  else if ( datatype == mpi_float ) then

    call mpi_reduce_real ( data1, data2, n, operation, ierror )

  else if ( datatype == mpi_complex ) then

    call mpi_reduce_complex ( data1, data2, n, operation, ierror )
  else if ( datatype == mpi_double_complex ) then

    call mpi_reduce_double_complex ( data1, data2, n, operation, ierror )   

  else

    ierror = MPI_FAILURE

  end if

  return
end
subroutine mpi_reduce_double_precision ( data1, data2, n, operation, ierror )

!*****************************************************************************80
!
!! MPI_REDUCE_DOUBLE_PRECISION: reduction operation on double precision values.
!
!  Discussion:
!
!    The reduction operations are sum, maximum, minimum, product.
!
!    Thanks to Simppa Akaslompolo for correcting this routine!
!    11 January 2012.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, double precision DATA1(N), the data to be processed.
!
!    Output, double precision DATA2(N), the value of the reduction operation.
!
!    Input, integer N, the number of items in DATA1.
!
!    Input, integer OPERATION, should have the value of one of the symbolic
!    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  include "mpi_dummy.fi"

  integer n

  double precision data1(n)
  double precision data2(n)
  integer ierror
  integer operation

  ierror = MPI_SUCCESS

  if ( operation == mpi_max ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_min ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_product ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_sum ) then

    data2(1:n) = data1(1:n)

  else

    ierror = MPI_FAILURE

  end if

  return
end


subroutine mpi_reduce_complex ( data1, data2, n, operation, ierror )

  implicit none

  include "mpi_dummy.fi"

  integer n

  complex(kind=4) data1(n)
  complex(kind=4) data2(n)
  integer ierror
  integer operation

  ierror = MPI_SUCCESS

  if ( operation == mpi_max ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_min ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_product ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_sum ) then

    data2(1:n) = data1(1:n)

  else

    ierror = MPI_FAILURE

  end if

  return
end subroutine mpi_reduce_complex



subroutine mpi_reduce_double_complex ( data1, data2, n, operation, ierror )

  implicit none

  include "mpi_dummy.fi"

  integer n

  complex(kind=8) data1(n)
  complex(kind=8) data2(n)
  integer ierror
  integer operation

  ierror = MPI_SUCCESS

  if ( operation == mpi_max ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_min ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_product ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_sum ) then

    data2(1:n) = data1(1:n)

  else

    ierror = MPI_FAILURE

  end if

  return
end subroutine mpi_reduce_double_complex


subroutine mpi_reduce_integer ( data1, data2, n, operation, ierror )

!*****************************************************************************80
!
!! MPI_REDUCE_INTEGER carries out a reduction operation on integers.
!
!  Discussion:
!
!    The reduction operations are sum, maximum, minimum, product.
!
!    Thanks to Simppa Akaslompolo for correcting this routine!
!    11 January 2012.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, integer DATA1(N), the data to be processed.
!
!    Output, integer DATA2(N), the value of the reduction operation.
!
!    Input, integer N, the number of items in DATA1.
!
!    Input, integer OPERATION, should have the value of one of the symbolic
!    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  include "mpi_dummy.fi"

  integer n

  integer data1(n)
  integer data2(n)
  integer ierror
  integer operation

  ierror = MPI_SUCCESS

  if ( operation == mpi_max ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_min ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_product ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_sum ) then

    data2(1:n) = data1(1:n)

  else

    ierror = MPI_FAILURE

  end if

  return
end
subroutine mpi_reduce_real ( data1, data2, n, operation, ierror )

!*****************************************************************************80
!
!! MPI_REDUCE_REAL carries out a reduction operation on reals.
!
!  Discussion:
!
!    The reduction operations are sum, maximum, minimum, product.
!
!    Thanks to Simppa Akaslompolo for correcting this routine!
!    11 January 2012.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Input, real DATA1(N), the data to be processed.
!
!    Output, real DATA2(N), the value of the reduction operation.
!
!    Input, integer N, the number of items in DATA1.
!
!    Input, integer OPERATION, should have the value of one of the symbolic
!    constants MPI_MAX, MPI_MIN, MPI_PRODUCT or MPI_SUM.
!
!    Output, integer IERROR, is nonzero if an error occurred.
!
  implicit none

  include "mpi_dummy.fi"

  integer n

  real data1(n)
  real data2(n)
  integer ierror
  integer operation

  ierror = MPI_SUCCESS

  if ( operation == mpi_max ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_min ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_product ) then

    data2(1:n) = data1(1:n)

  else if ( operation == mpi_sum ) then

    data2(1:n) = data1(1:n)

  else

    ierror = MPI_FAILURE

  end if

  return
end


subroutine mpi_wait ( request, status, ierror )
  implicit none

  include "mpi_dummy.fi"

  integer ierror
  integer request
  integer status(mpi_status_size)

  ierror = MPI_SUCCESS

  return
 end
subroutine mpi_waitall ( icount, irequest, status, ierror )

  implicit none

  include "mpi_dummy.fi"

  integer icount
  integer ierror
  integer irequest
  integer status(mpi_status_size)

  ierror = MPI_SUCCESS
  return
end
subroutine mpi_waitany ( icount, requests, index, status, ierror )

  implicit none

  include "mpi_dummy.fi"

  integer requests(*)
  integer icount
  integer ierror
  integer index
  integer status(mpi_status_size)

  ierror = MPI_SUCCESS

  return
end

function mpi_wtime ( )

!*****************************************************************************80
!
!! MPI_WTIME returns the elapsed wall clock time.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323.
!
!  Parameters:
!
!    Output, real ( kind = 8 ) MPI_WTIME, the elapsed wall clock time.
!
  implicit none

  integer count
  integer count_max
  integer count_rate
  real ( kind = 8 ) mpi_wtime

  call system_clock ( count, count_rate, count_max )

  mpi_wtime = real ( count, kind = 8 ) / real ( count_rate, kind = 8 )

  return
end
