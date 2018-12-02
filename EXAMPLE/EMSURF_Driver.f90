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

PROGRAM HODLR_BUTTERFLY_SOLVER_3D
    use z_BPACK_DEFS
	use EMSURF_MODULE
	
	use z_BPACK_structure
	use z_BPACK_factor
	use z_BPACK_constr
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
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	complex(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:)
	
	character(len=:),allocatable  :: string
	character(len=1024)  :: strings	
	character(len=6)  :: info_env	
	integer :: length,edge
	integer :: ierr
	integer*8 oldmode,newmode
	type(z_Hoption)::option	
	type(z_Hstat)::stats	
	type(z_mesh)::msh	
	type(z_hobf)::ho_bf,ho_bf_copy
	type(z_Hmat)::h_mat
	type(z_kernelquant)::ker
	type(quant_EMSURF),target::quant
	type(z_proctree)::ptree
	integer,allocatable:: groupmembers(:)	
	integer nmpi
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
    write(*,*) "HODLR_BUTTERFLY_SOLVER_3D"
    write(*,*) "   "
	endif
	call z_InitStat(stats)
	call z_SetDefaultOptions(option)
	
	time_tmp = 0
	
 	! register the user-defined function and type in ker 
	ker%FuncZmn=>Zelem_EMSURF
	ker%QuantApp=>quant	
	
    quant%integral_points=6
    allocate (quant%ng1(quant%integral_points), quant%ng2(quant%integral_points), quant%ng3(quant%integral_points), quant%gauss_w(quant%integral_points))
    call gauss_points(quant)

	
     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	CALL getarg(1, DATA_DIR)

	option%Nmin_leaf=40
	quant%mesh_normal=1
	! Refined_level=0
	
	! levelpara_control=0
	! ACA_tolerance_forward=1d-4
	option%tol_comp=1d-4
	! SVD_tolerance_factor=1d-4
	option%tol_Rdetect=3d-5
    ! Preset_level_butterfly=0
	option%nogeo=0
	quant%scaling=1d0
	quant%wavelength=2.0
	! Discret=0.05
	quant%RCS_static=2
    quant%RCS_Nsample=2000
    ! Optimizing_forward=0
    ! Fast_inverse=0
    ! Add_method_of_base_level=2
    ! quant%rank_approximate_para1=6.0
    ! quant%rank_approximate_para2=6.0
    ! quant%rank_approximate_para3=6.0
	option%tol_LS=1d-12
	! tfqmr_tolerance=1d-6
	option%tol_itersol=3d-3
	option%n_iter=1000
	option%tol_rand=5d-3
	! up_tolerance=1d-4
	! relax_lamda=1d0
	! SolvingMethod=1
	option%level_check=100
	! rank_tmp=7
	! schurinv=1
	! reducelevel_flag=0
	! directschur=1
	option%precon=DIRECT !HODLRPRECON !  HODLRPRECON ! NOPRECON !
	! verboselevel=2
	option%xyzsort=TM
	option%lnoBP=4000
	option%TwoLayerOnly=1
	quant%CFIE_alpha=1
    option%schulzorder=3
    option%schulzlevel=3000
	option%LRlevel=100
	option%ErrFillFull=0
	option%RecLR_leaf=BACA
	option%BACA_Batch=32
	
	option%ErrSol=1
	! option%LR_BLK_NUM=2
	option%format= HMAT!  HODLR ! 
	option%near_para=2.01d0
	option%verbosity=2
	option%ILU=0 	
	option%LR_BLK_NUM=1	
	option%forwardN15flag=0 
	
	! call MKL_set_num_threads(NUM_Threads)    ! this overwrites omp_set_num_threads for MKL functions 
	
    ! call OMP_set_dynamic(.true.)
    !call OMP_set_nested(.true.)
	

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
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "geometry modeling for "//trim(DATA_DIR)//"......"
    
	call geo_modeling_SURF(quant,ptree%Comm,DATA_DIR)
	
	! generate the list of points for clustering
	msh%Nunk=quant%Nunk
	allocate(msh%xyz(3,quant%Nunk))
    do edge=1, quant%Nunk
		msh%xyz(:,edge) = quant%xyz(:,quant%maxnode+edge)
    enddo	
	option%touch_para = 3* quant%minedgelength
	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "modeling finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1

if(option%format==HODLR)then		
	
		t1 = OMP_get_wtime()	
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "constructing HODLR formatting......"
		call z_Cluster_partition(ho_bf,option,msh,ker,z_element_Zmn_user,ptree)
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
		call z_BPACK_factorization(ho_bf,option,stats,ptree,msh)
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Cascading factorizing finished"
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "	
		end if

		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve......"
		call EM_solve_SURF(ho_bf,option,msh,quant,ptree,stats)
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve finished"
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "	
		

		call delete_quant_EMSURF(quant)
		
		call z_delete_proctree(ptree)
		call z_delete_Hstat(stats)
		call z_delete_mesh(msh)
		call z_delete_kernelquant(ker)
		call z_HODLR_delete(ho_bf)
	elseif(option%format==HMAT)then
		t1 = OMP_get_wtime()	
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "constructing H matrix formatting......"
		call z_Cluster_partition(h_mat,option,msh,ker,z_element_Zmn_user,ptree)
		call z_BPACK_structuring(h_mat,option,msh,ptree,stats)
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "H matrix formatting finished"
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
		t2 = OMP_get_wtime()
		! write(*,*)t2-t1
		! stop 
		
		!pause
		
		!call compression_test()
		t1 = OMP_get_wtime()	
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "H matrix construction......"
		call z_BPACK_construction(h_mat,option,stats,msh,ker,z_element_Zmn_user,ptree)
		! if(option%precon/=DIRECT)then
			! call copy_HOBF(h_mat,ho_bf_copy)	
		! end if    
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "H matrix construction finished"
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
		t2 = OMP_get_wtime()   
		! write(*,*)t2-t1
		
		if(option%precon/=NOPRECON)then								
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "H-LU factorizing......"
		call z_BPACK_factorization(h_mat,option,stats,ptree,msh)
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "H-LU factorizing finished"
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "	
		end if

		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve......"
		call EM_solve_SURF(h_mat,option,msh,quant,ptree,stats)
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve finished"
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "	
		

		call delete_quant_EMSURF(quant)
		
		call z_delete_proctree(ptree)
		call z_delete_Hstat(stats)
		call z_delete_mesh(msh)
		call z_delete_kernelquant(ker)
		call z_Hmat_delete(h_mat)	
	endif	
	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"
	
	call blacs_exit(1)
	call MPI_Finalize(ierr)
	
    ! ! ! ! pause

end PROGRAM HODLR_BUTTERFLY_SOLVER_3D






