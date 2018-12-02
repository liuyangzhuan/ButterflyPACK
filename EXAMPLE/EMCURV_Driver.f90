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

! Developers: Yang Liu, Xiaoye S. Li.
!             (Lawrence Berkeley National Lab, Computational Research Division).

PROGRAM HODLR_BUTTERFLY_SOLVER_2D
    use z_BPACK_DEFS
    use EMCURV_MODULE
	
	use z_BPACK_structure
	use z_BPACK_factor
	use z_BPACK_constr
	use z_BPACK_Utilities
	use omp_lib
	use z_misc
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
	type(z_Hoption)::option	
	type(z_Hstat)::stats
	type(z_mesh)::msh
	type(z_kernelquant)::ker
	type(z_hobf)::ho_bf
	type(z_Hmat)::h_mat
	integer,allocatable:: groupmembers(:)
	integer nmpi
	type(z_proctree)::ptree
	type(quant_EMCURV),target::quant
	CHARACTER (LEN=1000) DATA_DIR	
	
	
	! nmpi and groupmembers should be provided by the user 
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo	
	
	call z_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	deallocate(groupmembers)
	
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'NUMBER_MPI=',nmpi
	
 	threads_num=1
    CALL getenv("OMP_NUM_THREADS", strings)
	strings = TRIM(strings)	
	if(LEN_TRIM(strings)>0)then
		read(strings , *) threads_num
	endif
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'OMP_NUM_THREADS=',threads_num
	
	call OMP_set_num_threads(threads_num)		
		
		
	call DATE_AND_TIME(values=times)     ! Get the current time 
	seed_myid(1) = times(4) * (360000*times(5) + 6000*times(6) + 100*times(7) + times(8))
	
	! seed_myid(1) = myid*1000
	call RANDOM_SEED(PUT=seed_myid)
	
	! oldmode = vmlsetmode(VML_FTZDAZ_ON)
	! call vmlsetmode(VML_FTZDAZ_ON)
	if(ptree%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "HODLR_BUTTERFLY_SOLVER_2D"
    write(*,*) "   "
	endif
	
	call z_InitStat(stats)
	call z_SetDefaultOptions(option)
	
	time_tmp = 0

	
 	! register the user-defined function and type in ker 
	ker%FuncZmn=>Zelem_EMCURV
	ker%QuantApp=>quant

     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	! CALL getarg(1, DATA_DIR)
	
	
	! quant%model2d=10 ! Model type (1=strip; 2=corner reflector; 3=two opposite strips; 4=CR with RRS; 5=cylinder; 6=Rectangle Cavity); 7=half cylinder; 8=corrugated half cylinder; 9=corrugated corner reflector; 10=open polygon; 11=taller open polygon 
	! msh%Nunk=1280000
	! msh%Nunk=320000
	! msh%Nunk=80000
	! msh%Nunk=20000
	 ! msh%Nunk=5000 
    ! msh%Nunk=40000	
	! Refined_level=0
	
	
	
	
	! ker%Kernel = FULL	
	! allocate(ker%matZ_glo(msh%Nunk,msh%Nunk))
	! call RandomMat(msh%Nunk,msh%Nunk,msh%Nunk,ker%matZ_glo,0)
	! call MPI_Bcast(ker%matZ_glo,msh%Nunk*msh%Nunk,MPI_DOUBLE_COMPLEX,0,MPI_Comm_World,ierr)
	
	
	
	
	option%nogeo=0 
	! quant%wavelength=0.0006
	!quant%wavelength=0.0003
	! quant%wavelength=0.06

! quant%wavelength=0.08
! Discret=0.05
	quant%RCS_static=1
    quant%RCS_Nsample=100

	
	
	! quant%rank_approximate_para1=6.0
    ! quant%rank_approximate_para2=6.0
    ! quant%rank_approximate_para3=6.0
	
	
	option%Nmin_leaf=100
	! option%tol_comp=1d-4
	option%tol_Rdetect=3d-5	
	option%tol_LS=1d-12
	option%tol_itersol=1d-5
	option%n_iter=1000
	option%tol_rand=1d-3
	option%level_check=100
	option%precon= DIRECT ! HODLRPRECON ! NOPRECON !
	option%xyzsort=NATURAL !TM 
	option%lnoBP=40000
	option%TwoLayerOnly=1
    option%schulzorder=3
    option%schulzlevel=3000
	! option%LRlevel=100
	! option%ErrFillFull=0 
	! option%RecLR_leaf=ACA
	option%ErrSol=1
	! option%LR_BLK_NUM=2
	option%format= HMAT!  HODLR ! 
	option%near_para=0.01d0
	! option%verbosity=2
	option%ILU=0 
	option%forwardN15flag=0 
	
	call getarg(1,strings)
	read(strings,*)option%LR_BLK_NUM		
	call getarg(2,strings)
	read(strings,*)quant%model2d
	call getarg(3,strings)
	read(strings,*)msh%Nunk	
	call getarg(4,strings)
	read(strings,*)quant%wavelength	
	call getarg(5,strings)
	read(strings,*)option%tol_comp
	call getarg(6,strings)
	read(strings,*)option%ErrFillFull
	call getarg(7,strings)
	read(strings,*)option%RecLR_leaf		
	call getarg(8,strings)
	read(strings,*)option%BACA_Batch	
	call getarg(9,strings)
	read(strings,*)option%LRlevel

	
	
	! call MKL_set_num_threads(NUM_Threads)    ! this overwrites omp_set_num_threads for MKL functions 
	
    ! call OMP_set_dynamic(.true.)
    !call OMP_set_nested(.true.)
	
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
	
	! generate the list of points for clustering
	allocate(msh%xyz(2,quant%Nunk))
	do edge=1, quant%Nunk
		msh%xyz(:,edge) = quant%xyz(:,edge*2-1)
	enddo	
    option%touch_para = 3* quant%minedgelength
	
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "modeling finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1

	
if(option%format==HODLR)then		
	
	t1 = OMP_get_wtime()	
	call z_Cluster_partition(ho_bf,option,msh,ker,z_element_Zmn_user,ptree)
	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "constructing HODLR formatting......"
	call z_BPACK_structuring(ho_bf,option,msh,ptree,stats)
  
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "HODLR formatting finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1
	! stop 
	
    !pause
    
    !call compression_test()
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "HODLR construction......"
    ! call BPACK_construction(ho_bf,option,stats,msh,ker,element_Zmn_FULL,ptree)
    call z_BPACK_construction(ho_bf,option,stats,msh,ker,z_element_Zmn_user,ptree)
	! if(option%precon/=DIRECT)then
		! call copy_HOBF(ho_bf,ho_bf_copy)	
	! end if
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "HODLR construction finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	! write(*,*)t2-t1
	
	if(option%precon/=NOPRECON)then
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Cascading factorizing......"
    call z_BPACK_Factorization(ho_bf,option,stats,ptree,msh)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Cascading factorizing finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "	
	end if

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve......"
    call EM_solve_CURV(ho_bf,option,msh,quant,ptree,stats)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "	
	
	call z_PrintStat(stats,ptree)
	
	call delete_quant_EMCURV(quant)
	
	call z_delete_proctree(ptree)
	call z_delete_Hstat(stats)
	call z_delete_mesh(msh)
	call z_delete_kernelquant(ker)
	call z_HODLR_delete(ho_bf)
	
elseif(option%format==HMAT)then

	t1 = OMP_get_wtime()		
	call z_Cluster_partition(h_mat,option,msh,ker,z_element_Zmn_user,ptree)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "constructing H matrix formatting......"
	call z_BPACK_structuring(h_mat,option,msh,ptree,stats)	
  
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "H matrix formatting finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()	
	
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "H matrix construction......"
    call z_BPACK_construction(h_mat,option,stats,msh,ker,z_element_Zmn_user,ptree)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "H matrix construction finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	
	if(option%precon/=NOPRECON)then
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "H-LU factorizing......"
    call z_BPACK_Factorization(h_mat,option,stats,ptree,msh)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "H-LU factorizing finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "	
	end if	
	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve......"
    call EM_solve_CURV(h_mat,option,msh,quant,ptree,stats)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "	
	
	call z_PrintStat(stats,ptree)
	
	call delete_quant_EMCURV(quant)
	
	call z_delete_proctree(ptree)
	call z_delete_Hstat(stats)
	call z_delete_mesh(msh)
	call z_delete_kernelquant(ker)
	call z_Hmat_delete(h_mat)	
	
endif	
	
	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"
	
	call blacs_exit(1)
	call MPI_Finalize(ierr)
	
end PROGRAM HODLR_BUTTERFLY_SOLVER_2D



