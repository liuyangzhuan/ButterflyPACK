PROGRAM HODLR_BUTTERFLY_SOLVER_3D
    use MODULE_FILE
	use geometry_model
	use H_structure
	use cascading_factorization
	use EM_calculation
	use matrices_fill
	use omp_lib
	use MISC
    implicit none

	! include "mkl_vml.fi"	 
	
    real*8 para
    real*8 tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seed_myid(12)
	integer times(8)	
	real*8 t1,t2,x,y,z,r,theta,phi
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
	type(hobf)::ho_bf,ho_bf_copy
	type(kernelquant)::ker
	
 	threads_num=1
    CALL getenv("OMP_NUM_THREADS", strings)
	strings = TRIM(strings)	
	if(LEN_TRIM(strings)>0)then
		read(strings , *) threads_num
	endif
	write(*,*)'OMP_NUM_THREADS=',threads_num
	call OMP_set_num_threads(threads_num)
	
		
	call DATE_AND_TIME(values=times)     ! Get the current time 
	seed_myid(1) = times(4) * (360000*times(5) + 6000*times(6) + 100*times(7) + times(8))
	
	! seed_myid(1) = myid*1000
	call RANDOM_SEED(PUT=seed_myid)
	
	! oldmode = vmlsetmode(VML_FTZDAZ_ON)
	! call vmlsetmode(VML_FTZDAZ_ON)
	
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "HODLR_BUTTERFLY_SOLVER_3D"
    write(*,*) "FOR X64 COMPILER"
    write(*,*) "   "

	call InitStat(stats)
	
	! time_indexarray = 0
	! time_leastsquare = 0
	! time_buttermul = 0
	! time_buttermulinv = 0
	! time_kernelupdate = 0
	! time_memcopy = 0
	! time_gemm = 0
	! time_gemm1 = 0		   
    ! time_getvec = 0
	! time_resolve = 0
	
	time_tmp = 0
	
	msh%Origins=(/0d0,0d0,0d0/)
    msh%integral_points=6
    allocate (msh%ng1(msh%integral_points), msh%ng2(msh%integral_points), msh%ng3(msh%integral_points), msh%gauss_w(msh%integral_points))
    call gauss_points(msh)

	
     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	CALL getarg(1, DATA_DIR)
	ker%Kernel = EMSURF	
	
	option%Nmin_leaf=40
	msh%mesh_normal=1
	! Refined_level=0
	para=0.001
	! levelpara_control=0
	! ACA_tolerance_forward=1d-4
	option%tol_SVD=1d-4
	! SVD_tolerance_factor=1d-4
	option%tol_Rdetect=3d-5
    ! Preset_level_butterfly=0
	msh%scaling=1d0
	ker%wavelength=0.5
	! Discret=0.05
	ker%RCS_static=1
    ker%RCS_Nsample=2000
    ! Optimizing_forward=0
    ! Fast_inverse=0
    ! Add_method_of_base_level=2
    ker%rank_approximate_para1=6.0
    ker%rank_approximate_para2=6.0
    ker%rank_approximate_para3=6.0
	option%tol_LS=1d-12
	! tfqmr_tolerance=1d-6
	option%tol_itersol=3d-3
	option%N_iter=1000
	option%tol_rand=5d-3
	! up_tolerance=1d-4
	! relax_lamda=1d0
	! SolvingMethod=1
	option%level_check=100
	! rank_tmp=7
	! schurinv=1
	! reducelevel_flag=0
	! directschur=1
	option%precon=DIRECT
	! verboselevel=2
	option%xyzsort=1
	option%LnoBP=4000
	option%TwoLayerOnly=0
	ker%CFIE_alpha=1
	option%touch_para=3
    option%schulzorder=3
    option%schulzlevel=3000
	option%LRlevel=100
	
	! call MKL_set_num_threads(NUM_Threads)    ! this overwrites omp_set_num_threads for MKL functions 
	
    ! call OMP_set_dynamic(.true.)
    !call OMP_set_nested(.true.)
	

    !*********************************************************

    ker%omiga=2*pi/ker%wavelength/sqrt(mu0*eps0)
    ker%wavenum=2*pi/ker%wavelength

   !***********************************************************************
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'ker%wavelength:',ker%wavelength
   write (*,*) ''
   !***********************************************************************
	
	t1 = OMP_get_wtime()
    write(*,*) "geometry modeling......"
    call geo_modeling_SURF(msh,ker)
    write(*,*) "modeling finished"
    write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1

	t1 = OMP_get_wtime()	
    write(*,*) "constructing H_matrices formatting......"
    call H_matrix_structuring(ho_bf,para,option,msh)
	call BPlus_structuring(ho_bf,option,msh)
    write(*,*) "H_matrices formatting finished"
    write(*,*) "    "
	t2 = OMP_get_wtime()
	write(*,*)t2-t1
	! stop 
	
    !pause
    
    !call compression_test()
	t1 = OMP_get_wtime()	
    write(*,*) "H_matrices filling......"
    call matrices_filling(ho_bf,option,stats,msh,ker,element_Zmn_EMSURF)
	if(option%precon/=DIRECT)then
		call copy_HOBF(ho_bf,ho_bf_copy)	
	end if    
	write(*,*) "H_matrices filling finished"
    write(*,*) "    "
 	t2 = OMP_get_wtime()   
	write(*,*)t2-t1
	
    write(*,*) "Cascading factorizing......"
    call cascading_factorizing(ho_bf,option,stats)
    write(*,*) "Cascading factorizing finished"
    write(*,*) "    "	


    write(*,*) "EM_solve......"
    call EM_solve_SURF(ho_bf_copy,ho_bf,option,msh,ker)
    write(*,*) "EM_solve finished"
    write(*,*) "    "	
	

    write(*,*) "-------------------------------program end-------------------------------------"

    ! ! ! ! pause

end PROGRAM HODLR_BUTTERFLY_SOLVER_3D
