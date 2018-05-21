PROGRAM HODLR_BUTTERFLY_SOLVER_2D
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
	type(kernelquant)::ker
	type(hobf)::ho_bf,ho_bf_copy
	integer,allocatable:: groupmembers(:)
	integer nmpi
	type(proctree)::ptree
	
	
	
	
	! nmpi and groupmembers should be provided by the user 
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo	
	
	call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	
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
    write(*,*) "HODLR_BUTTERFLY_SOLVER_2D"
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
	! time_halfbuttermul = 0
	time_tmp = 0
	
	msh%Origins=(/0d0,0d0,0d0/)
 

     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	CALL getarg(1, DATA_DIR)
	
	para=0.001
	
	ker%Kernel = EMCURV	
	
	msh%model2d=10 ! Model type (1=strip; 2=corner reflector; 3=two opposite strips; 4=CR with RRS; 5=cylinder; 6=Rectangle Cavity); 7=half cylinder; 8=corrugated half cylinder; 9=corrugated corner reflector; 10=open polygon; 11=taller open polygon 
	msh%Nunk=5000
	
	! Refined_level=0
	



	
	
	msh%scaling=1d0
	ker%wavelength=0.08
	! Discret=0.05
	ker%RCS_static=1
    ker%RCS_Nsample=2000

	
	
	ker%rank_approximate_para1=6.0
    ker%rank_approximate_para2=6.0
    ker%rank_approximate_para3=6.0
	
	
	option%Nmin_leaf=50
	option%tol_SVD=1d-4
	option%tol_Rdetect=3d-5	
	option%tol_LS=1d-12
	option%tol_itersol=1d-6
	option%N_iter=1000
	option%tol_rand=1d-3
	option%level_check=100
	option%precon=DIRECT
	option%xyzsort=3
	option%LnoBP=40000
	option%TwoLayerOnly=1
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
    call geo_modeling_CURV(msh,ker)
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
    call matrices_filling(ho_bf,option,stats,msh,ker,element_Zmn_EMCURV)
	if(option%precon/=DIRECT)then
		call copy_HOBF(ho_bf,ho_bf_copy)	
	end if
    write(*,*) "H_matrices filling finished"
    write(*,*) "    "
 	t2 = OMP_get_wtime()   
	write(*,*)t2-t1
	
	if(option%precon/=NOPRECON)then
    write(*,*) "Cascading factorizing......"
    call cascading_factorizing(ho_bf,option,stats)
    write(*,*) "Cascading factorizing finished"
    write(*,*) "    "	
	end if

    write(*,*) "EM_solve......"
    call EM_solve_CURV(ho_bf_copy,ho_bf,option,msh,ker)
    write(*,*) "EM_solve finished"
    write(*,*) "    "	
	
	
    write(*,*) "-------------------------------program end-------------------------------------"

    ! ! ! ! pause

end PROGRAM HODLR_BUTTERFLY_SOLVER_2D
