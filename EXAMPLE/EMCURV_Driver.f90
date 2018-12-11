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


PROGRAM ButterflyPACK_IE_2D
    use BPACK_DEFS
    use EMCURV_MODULE
	
	use BPACK_structure
	use BPACK_Solve_Mul
	use BPACK_factor
	use BPACK_constr
	use BPACK_Utilities
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
	integer edge
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	complex(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:)
	
	character(len=:),allocatable  :: string
	character(len=1024)  :: strings	
	character(len=6)  :: info_env	
	integer :: length
	integer :: ierr
	integer*8 oldmode,newmode
	type(Hoption)::option	
	type(Hstat)::stats
	type(mesh)::msh
	type(kernelquant)::ker
	type(Bmatrix)::bmat
	integer,allocatable:: groupmembers(:)
	integer nmpi
	type(proctree)::ptree
	type(quant_EMCURV),target::quant
	CHARACTER (LEN=1000) DATA_DIR	
	integer:: randsize=50
	
	! nmpi and groupmembers should be provided by the user 
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo	
	
	! generate the process tree
	call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	deallocate(groupmembers)
	
	
	! set up and print MPI and thread number 
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'NUMBER_MPI=',nmpi
 	threads_num=1
    CALL getenv("OMP_NUM_THREADS", strings)
	strings = TRIM(strings)	
	if(LEN_TRIM(strings)>0)then
		read(strings , *) threads_num
	endif
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'OMP_NUM_THREADS=',threads_num
	call OMP_set_num_threads(threads_num)		
	

	! oldmode = vmlsetmode(VML_FTZDAZ_ON)
	! call vmlsetmode(VML_FTZDAZ_ON)	
	if(ptree%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "ButterflyPACK_IE_2D"
    write(*,*) "   "
	endif
	
	call InitStat(stats)
	call SetDefaultOptions(option)
	
 	! register the user-defined function and type in ker 
	ker%FuncZmn=>Zelem_EMCURV
	ker%QuantApp=>quant

	quant%RCS_static=1
    quant%RCS_Nsample=2000
	quant%model2d=10
	quant%wavelength=0.08
	
	msh%Nunk=5000	
	
	option%ErrSol=1
	option%format= HODLR! HMAT!  HODLR ! 
	option%near_para=0.01d0
	option%verbosity=2
	option%ILU=0  
	option%forwardN15flag=0 
    ! option%LRlevel=100

	if(iargc()>=1)then
		call getarg(1,strings)
		read(strings,*)option%LR_BLK_NUM	
	endif
	if(iargc()>=2)then	
		call getarg(2,strings)
		read(strings,*)quant%model2d
	endif
	if(iargc()>=3)then
		call getarg(3,strings)
		read(strings,*)msh%Nunk	
	endif
	if(iargc()>=4)then
		call getarg(4,strings)
		read(strings,*)quant%wavelength	
	endif
	if(iargc()>=5)then
		call getarg(5,strings)
		read(strings,*)option%tol_comp
	endif
	if(iargc()>=6)then
		call getarg(6,strings)
		read(strings,*)option%ErrFillFull
	endif
	if(iargc()>=7)then
		call getarg(7,strings)
		read(strings,*)option%RecLR_leaf		
	endif
	if(iargc()>=8)then
		call getarg(8,strings)
		read(strings,*)option%BACA_Batch	
	endif
	if(iargc()>=9)then
		call getarg(9,strings)
		read(strings,*)option%LRlevel	
	endif
	if(iargc()>=10)then
		call getarg(10,strings)
		read(strings,*)option%precon
	endif
	if(iargc()>=11)then
		call getarg(11,strings)
		read(strings,*)option%xyzsort
	endif	
	if(iargc()>=12)then
		call getarg(12,strings)
		read(strings,*)option%Nmin_leaf
	endif		
	
	quant%Nunk = msh%Nunk
    !*********************************************************

    quant%omiga=2*pi/quant%wavelength/sqrt(mu0*eps0)
    quant%wavenum=2*pi/quant%wavelength

   !***********************************************************************
   if(ptree%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'wavelength:',quant%wavelength
   write (*,*) ''
   endif
   !***********************************************************************
	
	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "geometry modeling......"
    call geo_modeling_CURV(quant,ptree%Comm)
	
	! generate the list of coordinates for clustering
	allocate(msh%xyz(2,quant%Nunk))
	do edge=1, quant%Nunk
		msh%xyz(:,edge) = quant%xyz(:,edge*2-1)
	enddo	
    option%touch_para = 3* quant%minedgelength
	
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "modeling finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()
	
	t1 = OMP_get_wtime()	
	call Cluster_partition(bmat,option,msh,ker,element_Zmn_user,ptree)
	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "constructing Hierarchical format......"
	call BPACK_structuring(bmat,option,msh,ptree,stats)
  
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Hierarchical format finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()
	
	
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Matrix construction......"
    call BPACK_construction_Element_Compute(bmat,option,stats,msh,ker,element_Zmn_user,ptree)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Matrix construction finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	
	
	if(option%precon/=NOPRECON)then
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Factor......"
    call BPACK_Factorization(bmat,option,stats,ptree,msh)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Factor finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "	
	end if

	if(option%ErrSol==1)then
		call BPACK_Test_Solve_error(bmat,msh%idxe-msh%idxs+1,option,ptree,stats)
	endif	
	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve......"
    call EM_solve_CURV(bmat,option,msh,quant,ptree,stats)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "	
	
	call PrintStat(stats,ptree)
	
	call delete_quant_EMCURV(quant)
	
	call delete_proctree(ptree)
	call delete_Hstat(stats)
	call delete_mesh(msh)
	call delete_kernelquant(ker)
	call BPACK_delete(bmat)
	
	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"
	
	call blacs_exit(1)
	call MPI_Finalize(ierr)
	
end PROGRAM ButterflyPACK_IE_2D




