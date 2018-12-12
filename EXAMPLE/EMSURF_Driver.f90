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


! This exmple works with double-complex precision data
#define DAT 0

#include "ButterflyPACK_config.fi"

PROGRAM ButterflyPACK_IE_3D
    use BPACK_DEFS
	use EMSURF_MODULE

	use BPACK_structure
	use BPACK_factor
	use BPACK_constr
	use BPACK_Solve_Mul
	use omp_lib
	use misc
    implicit none

	! include "mkl_vml.fi"

    real(kind=8) para
    real(kind=8) tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seed_myid(50)
	integer times(8)
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	complex(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:)

	character(len=:),allocatable  :: string
	character(len=1024)  :: strings
	character(len=6)  :: info_env
	integer :: length,edge
	integer :: ierr
	integer*8 oldmode,newmode
	type(Hoption)::option
	type(Hstat)::stats
	type(mesh)::msh
	type(Bmatrix)::bmat
	type(kernelquant)::ker
	type(quant_EMSURF),target::quant
	type(proctree)::ptree
	integer,allocatable:: groupmembers(:)
	integer nmpi
	CHARACTER (LEN=1000) DATA_DIR
	real(kind=8),allocatable::xyz(:,:)
	integer,allocatable::Permutation(:)
	integer Nunk_loc	
	
	! nmpi and groupmembers should be provided by the user
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo

	call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	deallocate(groupmembers)


	if(ptree%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "ButterflyPACK_IE_3D"
    write(*,*) "   "
	endif
	
	!**** initialize stats and option	
	call InitStat(stats)
	call SetDefaultOptions(option)


	!**** intialize the user-defined derived type quant
	! compute the quadrature rules
    quant%integral_points=6
    allocate (quant%ng1(quant%integral_points), quant%ng2(quant%integral_points), quant%ng3(quant%integral_points), quant%gauss_w(quant%integral_points))
    call gauss_points(quant)

    !*************************input******************************
	DATA_DIR='../EXAMPLE/EM3D_DATA/sphere_2300'

	quant%mesh_normal=1
	quant%scaling=1d0
	quant%wavelength=2.0
	quant%RCS_static=2
    quant%RCS_Nsample=1000
	quant%CFIE_alpha=1

	option%format= HMAT!  HODLR !
	option%near_para=2.01d0
	option%verbosity=2
	option%ILU=0
	option%forwardN15flag=0

	if(iargc()>=1)then
		CALL getarg(1, DATA_DIR)
	endif
	if(iargc()>=2)then
		call getarg(2,strings)
		read(strings,*)quant%wavelength
	endif
	if(iargc()>=3)then
		call getarg(3,strings)
		read(strings,*)option%precon
	endif
	if(iargc()>=4)then
		call getarg(4,strings)
		read(strings,*)option%xyzsort
	endif
	if(iargc()>=5)then
		call getarg(5,strings)
		read(strings,*)option%Nmin_leaf
	endif
	
    quant%omiga=2*pi/quant%wavelength/sqrt(mu0*eps0)
    quant%wavenum=2*pi/quant%wavelength
	! option%touch_para = 3* quant%minedgelength
	
   !***********************************************************************
	if(ptree%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'wavelength:',quant%wavelength
   write (*,*) ''
	endif
   !***********************************************************************

   
	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "geometry modeling for "//trim(DATA_DIR)//"......"
	call geo_modeling_SURF(quant,ptree%Comm,DATA_DIR)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "modeling finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()


	
	!**** initialization of the construction phase
	t1 = OMP_get_wtime()
	allocate(xyz(3,quant%Nunk))
	do ii=1, quant%Nunk
		xyz(:,ii) = quant%xyz(:,quant%maxnode+ii)
	enddo
    allocate(Permutation(quant%Nunk))
	call BPACK_construction_Element_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree,Zelem_EMSURF,quant,Coordinates=xyz)
	deallocate(Permutation) ! caller can use this permutation vector if needed 	
	deallocate(xyz)
	t2 = OMP_get_wtime()	
	
	
	
	!**** computation of the construction phase
	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Matrix construction......"
    call BPACK_construction_Element_Compute(bmat,option,stats,msh,ker,element_Zmn_user,ptree)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Matrix construction finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
 	t2 = OMP_get_wtime()


	!**** factorization phase
	if(option%precon/=NOPRECON)then
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Factor......"
	call BPACK_Factorization(bmat,option,stats,ptree,msh)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Factor finished"
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	if(option%ErrSol==1)then
		call BPACK_Test_Solve_error(bmat,msh%idxe-msh%idxs+1,option,ptree,stats)
	endif	
	end if

	!**** solve phase
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve......"
	call EM_solve_SURF(bmat,option,msh,quant,ptree,stats)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve finished"
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "

	
	!**** print statistics
	call PrintStat(stats,ptree)

	
	!**** deletion of quantities
	call delete_quant_EMSURF(quant)
	call delete_proctree(ptree)
	call delete_Hstat(stats)
	call delete_mesh(msh)
	call delete_kernelquant(ker)
	call BPACK_delete(bmat)


    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	call blacs_exit(1)
	call MPI_Finalize(ierr)

    ! ! ! ! pause

end PROGRAM ButterflyPACK_IE_3D






