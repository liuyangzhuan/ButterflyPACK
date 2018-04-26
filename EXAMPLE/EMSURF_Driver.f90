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
	real*8,parameter :: cd = 299792458d0
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

	time_indexarray = 0
	time_leastsquare = 0
	time_buttermul = 0
	time_buttermulinv = 0
	time_kernelupdate = 0
	time_memcopy = 0
	time_gemm = 0
	time_gemm1 = 0		   
    time_getvec = 0
	time_resolve = 0
	time_halfbuttermul = 0
	time_tmp = 0
	
	Origins=(/0d0,0d0,0d0/)
    Bigvalue=10000000.0d0
    junit=(0d0,1d0)
    pi=4d0*atan(1d0)
    eps0=1d7/(4d0*pi*cd**2)
    mu0=pi*4d-7
    gamma=1.781072418d0
    impedence=sqrt(mu0/eps0)
    integral_points=6
    allocate (ng1(integral_points), ng2(integral_points), ng3(integral_points), gauss_w(integral_points))
    call gauss_points()

	
	! x = 1d0
	! y = 1d0
	! z = 0d0
	! call Cart2Sph(x,y,z,Origins,r,theta,phi)
	! write(*,*)r,theta,phi,sin(theta)
	! stop
	
     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	CALL getarg(1, DATA_DIR)
	Kernel = EMSURF	
	
	Nmin_leaf=40
	mesh_normal=1
	Refined_level=0
	para=0.001
	levelpara_control=0
	ACA_tolerance_forward=1d-4
	SVD_tolerance_forward=1d-4
	SVD_tolerance_factor=1d-4
	Rank_detection_factor=3d-5
    Preset_level_butterfly=0
	Scale=1d0
	wavelength=0.5
	Discret=0.05
	Static=1
    RCS_sample=2000
    Optimizing_forward=0
    Fast_inverse=0
    Add_method_of_base_level=2
    rank_approximate_para1=6.0
    rank_approximate_para2=6.0
    rank_approximate_para3=6.0
	LS_tolerance=1d-12
	tfqmr_tolerance=1d-6
	tfqmr_tolerance_solving=3d-3
	iter_tolerance=5d-3
	up_tolerance=1d-4
	relax_lamda=1d0
	SolvingMethod=1
	level_tmp=100
	rank_tmp=7
	schurinv=1
	reducelevel_flag=0
	directschur=1
	preconditioner=0
	verboselevel=2
	xyzsort=1
	LnoBP=4
	TwoLayerOnly=0
	CFIE_alpha=1
	explicitflag=1
	fullmatflag=0
	touch_para=3
	
	! open (90,file='input.txt')
	! read (90,*)

	! read (90,*) Nmin_leaf
	! read (90,*) mesh_normal
	! read (90,*) Refined_level
	! read (90,*) para
	! read (90,*) levelpara_control
	! read (90,*) ACA_tolerance_forward
	! read (90,*) SVD_tolerance_forward
	! read (90,*) SVD_tolerance_factor
	! read (90,*) Rank_detection_factor
    ! read (90,*) Preset_level_butterfly
	! read (90,*) Scale
	! read (90,*) wavelength
	! read (90,*) Discret
	! read (90,*) Static
    ! read (90,*) RCS_sample
    ! read (90,*) Optimizing_forward
    ! read (90,*) Fast_inverse
    ! read (90,*) Add_method_of_base_level
    ! read (90,*) rank_approximate_para1
    ! read (90,*) rank_approximate_para2
    ! read (90,*) rank_approximate_para3
	! read (90,*) threads_num
	! read (90,*) LS_tolerance
	! read (90,*) tfqmr_tolerance
	! read (90,*) tfqmr_tolerance_solving
	! read (90,*) iter_tolerance
	! read (90,*) up_tolerance
	! read (90,*) relax_lamda
	! read (90,*) SolvingMethod
	! read (90,*) level_tmp
	! read (90,*) rank_tmp
	! read (90,*) schurinv
	! read (90,*) reducelevel_flag
	! read (90,*) directschur
	! read (90,*) preconditioner
	! read (90,*) verboselevel
	! read (90,*) xyzsort
	! read (90,*) LnoBP
	! read (90,*) TwoLayerOnly
	! read (90,*) CFIE_alpha
	
	! close (90)
	!Nmin_leaf=250
	!para=0.01
	!tolerance=0.001
	!Scale=1.
	!alpha=0.5
    !wavelength=2.
    tolerance=ACA_tolerance_forward
    
	! call MKL_set_num_threads(NUM_Threads)    ! this overwrites omp_set_num_threads for MKL functions 
	
    ! call OMP_set_dynamic(.true.)
    !call OMP_set_nested(.true.)
	

    !*********************************************************

    omiga=2*pi/wavelength/sqrt(mu0*eps0)
    wavenum=2*pi/wavelength

   !***********************************************************************
   open (256,file='Info.txt')
   write (256,*) 'EFIE computing'
   write (256,*) 'wavelength:',wavelength
   close (256)
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'wavelength:',wavelength
   write (*,*) ''
   !***********************************************************************
	
	t1 = OMP_get_wtime()
    write(*,*) "geometry modeling......"
    call geo_modeling()
    write(*,*) "modeling finished"
    write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1

	t1 = OMP_get_wtime()	
    write(*,*) "constructing H_matrices formatting......"
    call H_matrix_structuring(para)
	call BPlus_structuring()
    write(*,*) "H_matrices formatting finished"
    write(*,*) "    "
	t2 = OMP_get_wtime()
	write(*,*)t2-t1
	! stop 
	
    !pause
    
    !call compression_test()
	t1 = OMP_get_wtime()	
    write(*,*) "H_matrices filling......"
    call matrices_filling(tolerance)
    write(*,*) "H_matrices filling finished"
    write(*,*) "    "
 	t2 = OMP_get_wtime()   
	write(*,*)t2-t1
	
    write(*,*) "Cascading factorizing......"
    call cascading_factorizing(tolerance)
    write(*,*) "Cascading factorizing finished"
    write(*,*) "    "	


    write(*,*) "EM_calculating......"
    call EM_calculating()
    write(*,*) "EM_calculating finished"
    write(*,*) "    "	
	
	! ! pause 
	
    !call butterfly_extraction_test()
    !pause

    !pause

    !call lowrank_test()

    !pause

!   write(*,*) "full matrix filling......"
!   call full_matrix()
!   write(*,*) "full matrix filling finished"
!   write(*,*) "    "
!!
!   write(*,*) "compare test......"
!   call compare_test0()
!   write(*,*) "compare test finished"
!   write(*,*) "    "

!    write(*,*) " generating random vectors......"
!    nn=20; mm=Maxedge
!    allocate (vectors_1(Maxedge,nn),vectors_2(Maxedge,nn),vectors_3(Maxedge,nn))
!    call random_vectors(Maxedge,nn)
!    write(*,*) "    "

!   write(*,*) " H multiplying vectors......"
!   vectors_2=(0.,0.)
!       call H_multiply_vector0(nn,1)
!   write(*,*) "    "


    ! ! ! write(*,*) "MLMDA LL Decompositing......"
    ! ! ! write(*,*) ''
    ! ! ! if (Fast_inverse==1) then
        ! ! ! write (*,*) 'Fast Inverse Option: Yes'
    ! ! ! else
        ! ! ! write (*,*) 'Fast Inverse Option: No'
    ! ! ! endif
    ! ! ! if (Add_method_of_base_level==1) then
        ! ! ! write (*,*) 'Addition method of base level: SVD'
    ! ! ! elseif (Add_method_of_base_level==2) then
        ! ! ! write (*,*) 'Addition method of base level: ACA'
    ! ! ! elseif (Add_method_of_base_level==3) then
        ! ! ! write (*,*) 'Addition method of base level: Regulation'
    ! ! ! endif

    ! ! !Primary_block=0
    ! ! call MLMDA_LUDecomposition()
    ! ! write(*,*) "MLMDA LL Decomposition finished"
    ! ! write(*,*) "    "
    ! ! ! pause

!   write(*,*) "full matrix LUD......"
!   call full_matrix_LUD()
!   !call full_matrix_inverse()
!   !pause
!   write(*,*) "full matrix LUD finished"
!   write(*,*) "    "

    ! ! ! ! write(*,*) "EM_calculating......"
    ! ! ! ! call EM_calculating()
    ! ! ! ! write(*,*) "EM_calculating finished"
    ! ! ! ! write(*,*) "    "

!    write(*,*) " compare test......"
!    call compare_test1()
!    write(*,*) "compare test finished"
!    write(*,*) "    "

!    write(*,*) "H multiplying vectors......"
!   vectors_3=(0.,0.)
!       call H_multiply_vector0(nn,2)
!   write(*,*) "    "
!
!    write(*,*) "comparing results......"
!    call compare_vectors(mm,nn)
!    write(*,*) "    "

    write(*,*) "-------------------------------program end-------------------------------------"

    ! ! ! ! pause

end PROGRAM HODLR_BUTTERFLY_SOLVER_3D
