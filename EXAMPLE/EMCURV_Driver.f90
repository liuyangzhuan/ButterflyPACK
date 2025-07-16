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
!> @file
!> @brief This is an example that solves a 2D EFIE system for electromagnetics scattering.
!> @details Note that the use of the following \n
!> #define DAT 0 \n
!> #include "zButterflyPACK_config.fi" \n
!> which will macro replace subroutine, function, type names with those defined in SRC_DOUBLECOMLEX with double-complex precision

! This exmple works with double-complex precision data
PROGRAM ButterflyPACK_IE_2D
    use z_BPACK_DEFS
    use EMCURV_MODULE

	use z_BPACK_structure
	use z_BPACK_Solve_Mul
	use z_BPACK_factor
	use z_BPACK_constr
	use z_BPACK_Utilities
#ifdef HAVE_OPENMP
	use omp_lib
#endif
	use z_MISC_Utilities
	use z_BPACK_wrapper
#ifdef HAVE_ZFP
  use zfp
#endif
  use iso_c_binding

    implicit none

	! include "mkl_vml.fi"

    real(kind=8) para
    real(kind=8) tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seed_myid(50)
	integer times(8)
	integer edge
	real(kind=8) t1,t2,t3, x,y,z,r,theta,phi,retval(1)
	complex(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:)

	character(len=:),allocatable  :: string
	character(len=1024)  :: strings,strings1
	character(len=6)  :: info_env
	integer :: length
	integer :: ierr,provided
	integer*8 oldmode,newmode
	type(z_Hoption),target::option,option1
	type(z_Hstat),target::stats,stats1
	type(z_mesh),target::msh,msh1
	type(z_kernelquant),target::ker,ker1
	type(z_Bmatrix),target::bmat,bmat1
	integer,allocatable:: groupmembers(:)
	integer nmpi
	type(z_proctree),target::ptree,ptree1
	type(quant_EMCURV),target::quant
	type(z_quant_bmat),target::quant1
	CHARACTER (LEN=1000) DATA_DIR
	integer:: randsize=50
	real(kind=8),allocatable::xyz(:,:)
	integer,allocatable::Permutation(:),tree(:)
	integer Nunk_loc,Maxlevel
	integer nargs,flag
	integer v_major,v_minor,v_bugfix
	integer parent



#if 0
  ! input/decompressed arrays
  integer xLen, yLen
  real (kind=8), dimension(:, :), allocatable, target :: input_array
  real (kind=8), dimension(:, :), allocatable, target :: decompressed_array
  type(c_ptr) :: array_c_ptr
  real (kind=8) error, max_abs_error

  ! zfp_field
  type(zFORp_field) :: field

  ! bitstream
  character, dimension(:), allocatable, target :: buffer,buffer1
  type(c_ptr) :: buffer_c_ptr
  integer (kind=8) buffer_size_bytes, bitstream_offset_bytes
  type(zFORp_bitstream) :: bitstream, queried_bitstream

  ! zfp_stream
  type(zFORp_stream) :: stream, dstream
  real (kind=8) :: desired_rate, rate_result, tol_result
  integer :: dims, wra
  integer :: zfp_type, is_success


  ! initialize input and decompressed arrays
  xLen = 512
  yLen = 512
  allocate(input_array(xLen, yLen))
  do i = 1, xLen
    do j = 1, yLen
      CALL RANDOM_NUMBER(r)
      input_array(i, j) = r ! i * i + j * (j + 1)
    enddo
  enddo

  allocate(decompressed_array(xLen, yLen))

  ! setup zfp_field
  array_c_ptr = c_loc(input_array)
  zfp_type = zFORp_type_double
  field = zFORp_field_2d(array_c_ptr, zfp_type, xLen, yLen)

  ! setup bitstream
  buffer_size_bytes = xLen*yLen*8
  allocate(buffer(buffer_size_bytes))
  buffer_c_ptr = c_loc(buffer)
  bitstream = zFORp_bitstream_stream_open(buffer_c_ptr, buffer_size_bytes)

  ! setup zfp_stream
  stream = zFORp_stream_open(bitstream)

  tol_result=zFORp_stream_set_accuracy(stream,1d-3)
!   is_success = zFORp_stream_set_execution(stream,zFORp_exec_omp)
!   is_success = zFORp_stream_set_omp_threads(stream,0)



  ! compress
  t1 = MPI_Wtime()
  bitstream_offset_bytes = zFORp_compress(stream, field)
  t2 = MPI_Wtime()
  write(*, *) "After compression, bitstream offset at "
  write(*, *) bitstream_offset_bytes
  write(*,*)'compression time: ', t2-t1
  call zFORp_field_free(field)
  call zFORp_bitstream_stream_close(bitstream)



  ! decompress
  allocate(buffer1(bitstream_offset_bytes))
  buffer1=buffer(1:bitstream_offset_bytes)
  deallocate(buffer)
  buffer_c_ptr = c_loc(buffer1)
  bitstream = zFORp_bitstream_stream_open(buffer_c_ptr, buffer_size_bytes)

#if 1
  ! this works.
  dstream=stream
  call zFORp_stream_set_bit_stream(dstream, bitstream);
#else
  ! this also works. but the decompression stream needs to use the same option as the compression stream
  call zFORp_stream_close(stream)
  dstream = zFORp_stream_open(bitstream)
  tol_result=zFORp_stream_set_accuracy(dstream,1d-3)
#endif


  array_c_ptr = c_loc(decompressed_array)
  field = zFORp_field_2d(array_c_ptr, zfp_type, xLen, yLen)

  t1 = MPI_Wtime()
  bitstream_offset_bytes = zFORp_decompress(dstream, field)
  t2 = MPI_Wtime()
  write(*, *) "After decompression, bitstream offset at "
  write(*, *) bitstream_offset_bytes
  write(*,*)'decompression time: ', t2-t1

  max_abs_error = 0
  do i = 1, xLen
    do j = 1, yLen
      error = abs(decompressed_array(i, j) - input_array(i, j))
      max_abs_error = max(error, max_abs_error)
    enddo
  enddo
  write(*, *) "Max absolute error: "
  write(*, *) max_abs_error

  ! deallocations
  call zFORp_stream_close(dstream)
  call zFORp_bitstream_stream_close(bitstream)
  call zFORp_field_free(field)

  deallocate(buffer1)
  deallocate(input_array)
  deallocate(decompressed_array)

#endif




	! nmpi and groupmembers should be provided by the user
	! call MPI_Init(ierr)
	call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
	call MPI_COMM_GET_PARENT(parent, ierr) ! YL: this is needed if this function is spawned by a master process
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo

	! generate the process tree
	call z_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	deallocate(groupmembers)


	if(ptree%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "ButterflyPACK_IE_2D"
	call z_BPACK_GetVersionNumber(v_major,v_minor,v_bugfix)
	write(*,'(A23,I1,A1,I1,A1,I1,A1)') " ButterflyPACK Version:",v_major,".",v_minor,".",v_bugfix
    write(*,*) "   "
	endif

	!**** initialize stats and option
	call z_InitStat(stats)
	call z_SetDefaultOptions(option)

	!**** intialize the user-defined derived type quant
	quant%RCS_static=1
    quant%RCS_Nsample=1000
	quant%model2d=10
	quant%wavelength=0.08
	quant%freq=1/quant%wavelength/sqrt(BPACK_mu0*BPACK_eps0)
	quant%Nunk=5000

	option%ErrSol=1
	option%format=  HODLR !HMAT!  HODLR !
	option%near_para=0.01d0
	option%verbosity=2
	option%ILU=0
	option%forwardN15flag=0
        ! option%schulzlevel=0
        ! option%LRlevel=100
       ! option%level_check=1
    option%tol_itersol=1d-5
    ! option%sample_heuristic=1
    ! option%sample_para=4d0

	! do ii=2,15
	! call z_C_BPACK_TreeIndex_Merged2Child(ii,jj)
	! write(*,*)ii,jj
	! enddo

	nargs = iargc()
	ii=1
	do while(ii<=nargs)
		call getarg(ii,strings)
		if(trim(strings)=='-quant')then ! user-defined quantity parameters
			flag=1
			do while(flag==1)
				ii=ii+1
				if(ii<=nargs)then
					call getarg(ii,strings)
					if(strings(1:2)=='--')then
						ii=ii+1
						call getarg(ii,strings1)
						if(trim(strings)=='--model2d')then
							read(strings1,*)quant%model2d
						else if	(trim(strings)=='--nunk')then
							read(strings1,*)quant%Nunk
						else if	(trim(strings)=='--wavelength')then
							read(strings1,*)quant%wavelength
							quant%freq=1/quant%wavelength/sqrt(BPACK_mu0*BPACK_eps0)
						else if (trim(strings)=='--freq')then
							read(strings1,*)quant%freq
							quant%wavelength=1/quant%freq/sqrt(BPACK_mu0*BPACK_eps0)
						else
							if(ptree%MyID==Main_ID)write(*,*)'ignoring unknown quant: ', trim(strings)
						endif
					else
						flag=0
					endif
				else
					flag=0
				endif
			enddo
		else if(trim(strings)=='-option')then ! options of ButterflyPACK
			call z_ReadOption(option,ptree,ii)
		else
			if(ptree%MyID==Main_ID)write(*,*)'ignoring unknown argument: ',trim(strings)
			ii=ii+1
		endif
	enddo




    quant%wavenum=2*BPACK_pi/quant%wavelength


   !***********************************************************************
   if(ptree%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'frequency:',quant%freq
   write (*,*) 'wavelength:',quant%wavelength
   write (*,*) ''
   endif
   !***********************************************************************

	!**** geometry generalization and discretization
    call geo_modeling_CURV(quant,ptree%Comm)

	option%touch_para = 3* quant%minedgelength

	!**** register the user-defined function and type in ker
	ker%QuantApp => quant
	ker%FuncZmn => Zelem_EMCURV

	!**** initialization of the construction phase
	allocate(xyz(2,quant%Nunk))
	do edge=1, quant%Nunk
		xyz(:,edge) = quant%xyz(:,edge*2-1)
	enddo
    allocate(Permutation(quant%Nunk))
	call z_PrintOptions(option,ptree)
	call z_BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree,Coordinates=xyz)
	deallocate(Permutation) ! caller can use this permutation vector if needed
	deallocate(xyz)

	!**** computation of the construction phase
    call z_BPACK_construction_Element(bmat,option,stats,msh,ker,ptree)


	select case(option%format)
	case(HODLR)
		Maxlevel=bmat%ho_bf%Maxlevel
	case(HMAT)
		Maxlevel=bmat%h_mat%Maxlevel
	case(HSS)
		Maxlevel=bmat%hss_bf%Maxlevel
	end select
	if(.false.)then
	! if(option%lnoBP>Maxlevel)then	 ! haven't implement the following for Bplus.

		!*********** Construct the second HODLR by using the first HODLR as parallel element extraction
		call z_CopyOptions(option,option1)
		option1%nogeo=1   ! this indicates the second HOLDR construction requires no geometry information
		option1%xyzsort=NATURAL ! this indicates the second HOLDR construction requires no reordering
		option1%elem_extract=1 ! this indicates the second HOLDR construction uses user-defined parallel element extraction

		!**** register the user-defined function and type in ker
		ker1%FuncZmn=>Zelem_EMCURV
		ker1%FuncZmnBlock=>z_Zelem_block_Extraction
		ker1%QuantApp=>quant1

		quant1%bmat=>bmat
		quant1%msh=>msh
		quant1%ptree=>ptree
		quant1%stats=>stats
		quant1%option=>option

		msh1%Nunk = msh%Nunk

		allocate(xyz(2,quant%Nunk))
		do edge=1, quant%Nunk
			xyz(:,edge) = quant%xyz(:,edge*2-1)
		enddo
		allocate(msh1%xyz(2,quant%Nunk)) ! xyz is used to compute group centers to determine admissibility
		do edge=1,quant%Nunk
			msh1%xyz(:,edge) = xyz(:,msh%new2old(edge))
		enddo
		deallocate(xyz)

		call z_InitStat(stats1)

		! !**** generate the process tree for the second HODLR, can use larger number of MPIs if you want to
		! allocate(groupmembers(nmpi))
		! do ii=1,nmpi
			! groupmembers(ii)=(ii-1)
		! enddo
		! call z_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree1)
		! deallocate(groupmembers)


		!**** use the clustering tree from the first HODLR
		allocate (tree(2**Maxlevel))
		do ii=1,2**Maxlevel
			tree(ii)=msh%basis_group(2**Maxlevel+ii-1)%tail-msh%basis_group(2**Maxlevel+ii-1)%head+1
		enddo

		!**** initialization of the construction phase
		allocate(Permutation(msh1%Nunk))
		call z_BPACK_construction_Init(msh1%Nunk,Permutation,Nunk_loc,bmat1,option1,stats1,msh1,ker1,ptree,tree=tree)
		deallocate(Permutation) ! caller can use this permutation vector if needed
		deallocate(tree)

		!**** computation of the construction phase
		call z_BPACK_construction_Element(bmat1,option1,stats1,msh1,ker1,ptree)


		!**** deletion of quantities
		! call z_delete_proctree(ptree1)
		call z_delete_Hstat(stats1)
		call z_delete_mesh(msh1)
		call z_delete_kernelquant(ker1)
		call z_BPACK_delete(bmat1)
	endif



	!**** factorization phase
    call z_BPACK_Factorization(bmat,option,stats,ptree,msh)

	! call z_BPACK_PrintStructure(bmat,0, option,stats,ptree)
	! call z_BPACK_PrintStructure(bmat,1, option,stats,ptree)


	! !**** solve phase
	! call EM_solve_CURV(bmat,option,msh,quant,ptree,stats)



	!**** print statistics
	call z_PrintStat(stats,ptree)
	retval(1)=stats%Time_Fill+stats%Time_Factor+stats%Time_Sol

	!**** deletion of quantities
	call delete_quant_EMCURV(quant)
	call z_delete_proctree(ptree)
	call z_delete_Hstat(stats)
	call z_delete_mesh(msh)
	call z_delete_kernelquant(ker)
	call z_BPACK_delete(bmat)




    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	if(parent/= MPI_COMM_NULL)then
		call MPI_REDUCE(retval, MPI_BOTTOM, 1, MPI_double_precision, MPI_MAX, Main_ID, parent,ierr)
		call MPI_Comm_disconnect(parent,ierr)
	endif

	call z_blacs_exit_wrp(1)
	call MPI_Finalize(ierr)


end PROGRAM ButterflyPACK_IE_2D





