PROGRAM MLMDA_DIRECT_SOLVER_3D_CFIE

    use MODULE_FILE
	use geometry_model
	use H_structure
	use cascading_factorization
	use EM_calculation
	use matrices_fill
	use omp_lib
	use HODLR_randomMVP
    implicit none

    real*8 para,error
    real*8 tolerance
    integer Primary_block, nn, mm
    integer i,j,k, threads_num,ii,jj
	real*8,parameter :: cd = 299792458d0
	integer seed_myid(12)
	integer times(8)	
	real*8 t1,t2,t3,t4,x,y,z,r,theta,phi,tmp(3),Memory
	complex(kind=8),allocatable:: InputVec(:)
	complex(kind=8):: ctemp
	integer Ntunnel,kk,black_step,rankmax
	complex(kind=8),allocatable::Vout1(:,:),Vout2(:,:),Vin(:,:)
	character(len=1024)  :: strings
	
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
	
    !real scale


    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "MLMDA_DIRECT_SOLVER_3D_CFIE_NEW"
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
	
	
	Kernel = FULL

	Nmin_leaf=100
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
	wavelength=0.25
	Discret=0.05
	Static=2
    RCS_sample=1000
    Optimizing_forward=0
    Fast_inverse=0
    Add_method_of_base_level=2
    rank_approximate_para1=6.0
    rank_approximate_para2=6.0
    rank_approximate_para3=6.0
	LS_tolerance=1d-10
	tfqmr_tolerance=1d-6
	tfqmr_tolerance_solving=3d-3
	iter_tolerance=1d0
	up_tolerance=1
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
	LnoBP=600
	TwoLayerOnly=1
	CFIE_alpha=1
	explicitflag=1
	fullmatflag=1

	
	
	
	
	! open (90,file=trim(DATA_DIR)//'/input.txt')
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
	! read (90,*) explicitflag
	! read (90,*) fullmatflag
	
	! close (90)
	
	
	
	
	
	!Nmin_leaf=250
	!para=0.01
	!tolerance=0.001
	!Scale=1.
	!alpha=0.5
    !wavelength=2.
    tolerance=ACA_tolerance_forward
    ! call OMP_set_dynamic(.true.)
    !call OMP_set_nested(.true.)

    !*********************************************************

    omiga=2*pi/wavelength/sqrt(mu0*eps0)
    wavenum=2*pi/wavelength

	Ntunnel = 15600
	Maxedge = 3720
	Nin1 = 320*2 
	Nout1 = 610*2
	Nin2 = 320*2 
	Nout2 = 610*2	
	call assert(Nin1+Nout1+Nin2+Nout2==Maxedge,'The two surfaces have mismatched number of unknowns')
	
	! predefine the first three levels of tree due to the physical meanings
	allocate(basis_group_pre(7,2))
	basis_group_pre(1,1) = 1
	basis_group_pre(1,2) = Maxedge
	basis_group_pre(2,1) = 1
	basis_group_pre(2,2) = Nin1+Nout1	
	basis_group_pre(3,1) = Nin1+Nout1+1
	basis_group_pre(3,2) = Maxedge	
	basis_group_pre(4,1) = 1
	basis_group_pre(4,2) = Nin1	
	basis_group_pre(5,1) = 1+Nin1
	basis_group_pre(5,2) = Nout1+Nin1
	basis_group_pre(6,1) = 1+Nin1+Nout1
	basis_group_pre(6,2) = Nin2+Nin1+Nout1
	basis_group_pre(7,1) = 1+Nin2+Nin1+Nout1
	basis_group_pre(7,2) = Nout2+Nin2+Nin1+Nout1
	
	
	
	write(*,*)'Blackbox HODLR for scattering matrix compression'
	
	write(*,'(A10,I9,A11,I9)')' Ntunnel: ',Ntunnel,' Nsurface: ',Maxedge
	
	allocate(xyz(3,Maxedge))
	allocate(node_patch_of_edge(0:0,Maxedge))
	do kk=1,Maxedge
		node_patch_of_edge(0,kk)=kk
	enddo
	!!!!!! for simplicity, I'm duplicating locations of each point

	open(unit=521,file=trim(DATA_DIR)//'/edge_cen.out',status='old')
	do kk=1,Ntunnel/2
		read(521,*) tmp(1),tmp(2),tmp(3)
	end do          
	do kk=1,Maxedge/2
		read(521,*) xyz(1,2*kk-1),xyz(2,2*kk-1),xyz(3,2*kk-1)
		xyz(:,2*kk) = xyz(:,2*kk-1)
	end do
	close(521)
	
	
	
   ! !***********************************************************************
   ! open (256,file='Info.txt')
   ! write (256,*) 'CFIE computing'
   ! write (256,*) 'wavelength:',wavelength
   ! close (256)
   ! write (*,*) ''
   ! write (*,*) 'CFIE computing'
   ! write (*,*) 'wavelength:',wavelength
   ! write (*,*) ''
   ! !***********************************************************************
	
	! t1 = OMP_get_wtime()
    ! write(*,*) "geometry modeling......"
    ! call geo_modeling()
    ! write(*,*) "modeling finished"
    ! write(*,*) "    "
	! t2 = OMP_get_wtime()
	! ! write(*,*)t2-t1

	
	t1 = OMP_get_wtime()	
    write(*,*) "constructing H_matrices formatting......"
    call H_matrix_structuring(para)
	call BPlus_structuring()
    write(*,*) "H_matrices formatting finished"
    write(*,*) "    "
	t2 = OMP_get_wtime()
	write(*,*)t2-t1,'secnds'	
	! stop 
	
    !pause
    
	
	
	if(explicitflag ==1)then
		if(fullmatflag==1)then ! read the fullmat if it's already stored 
			t1 = OMP_get_wtime()
			write(*,*) "Generating fullmat ......"
			allocate(matZ_glo(Maxedge,Maxedge))
			matZ_glo = 0
			
			open(unit=888,file=trim(DATA_DIR)//'/fullmat.out',status='old')
			do ii=1,Maxedge
			do kk=1,Maxedge
				read(888,*)tmp(1),tmp(2)
				matZ_glo(kk,ii) = cmplx(tmp(1),tmp(2),dp) 
			end do
			end do
			close(unit=888)
			
			write(*,*) "Generating fullmat finished"
			t2 = OMP_get_wtime()   
			write(*,*)t2-t1, 'secnds'			
		else 			      ! generate the fullmat on-the-fly and store it in a file 
			write(*,*)'removed from previous versions'
		end if
	
		t1 = OMP_get_wtime()	
		write(*,*) "H_matrices filling......"
		call matrices_filling(tolerance)
		write(*,*) "H_matrices filling finished"
		write(*,*) "    "
		t2 = OMP_get_wtime()   
		write(*,*)t2-t1, 'secnds'			
	


		allocate(Vin(Maxedge,1))
		allocate(Vout1(Maxedge,1))
		allocate(Vout2(Maxedge,1))
		do ii=1,Maxedge
			Vin(ii,1) = random_complex_number()
		end do
		
		call MVM_Z_forward(Maxedge,1,Vin(:,1),Vout1(:,1),cascading_factors)
		
		do ii=1,Maxedge
			ctemp = 0d0
			do jj=1,Maxedge
				ctemp = ctemp + matZ_glo(new2old(ii),new2old(jj))*Vin(jj,1)
			end do
			Vout2(ii,1) = ctemp
		end do
		error = fnorm(Vout2-Vout1,Maxedge,1)/fnorm(Vout2,Maxedge,1)
		deallocate(Vin,Vout1,Vout2)
		
		write(*,*)error,'accuracy of construction'


		deallocate(matZ_glo)
		
	else if(explicitflag ==0)then
	
	
		if(fullmatflag==1)then ! read the fullmat and use it only for MVP
			t1 = OMP_get_wtime()
			write(*,*) "Reading fullmat ......"
			allocate(matZ_glo(Maxedge,Maxedge))
			matZ_glo = 0
			
			open(unit=888,file=trim(DATA_DIR)//'/fullmat.out',status='old')
			do ii=1,Maxedge
			do kk=1,Maxedge
				read(888,*)tmp(1),tmp(2)
				matZ_glo(kk,ii) = cmplx(tmp(1),tmp(2),dp) 
			end do
			end do
			close(unit=888)
			
			t2 = OMP_get_wtime()   
			write(*,*) "Reading fullmat finished ",t2-t1, 'secnds'

		! else if(fullmatflag==-1)then ! read the fullmat, explicitly contruct HODLR and use it only for MVP

			! t1 = OMP_get_wtime()
			! write(*,*) "Reading fullmat ......"
			! allocate(matZ_glo(Maxedge,Maxedge))
			! matZ_glo = 0
			
			! open(unit=888,file='fullmat.out',status='old')
			! do ii=1,Maxedge
			! do kk=1,Maxedge
				! read(888,*)tmp(1),tmp(2)
				! matZ_glo(kk,ii) = cmplx(tmp(1),tmp(2),dp) 
			! end do
			! end do
			! close(unit=888)
			
			! t2 = OMP_get_wtime()   
			! write(*,*) "Reading fullmat finished ",t2-t1, 'secnds'
			
			! t1 = OMP_get_wtime()	
			! write(*,*) "H_matrices filling......"
			! call matrices_filling(tolerance)
			! t2 = OMP_get_wtime()  
			! write(*,*) "H_matrices filling finished",t2-t1, 'secnds'

		else if(fullmatflag==0)then  ! true black-box MVP based HODLR construction
			write(*,*)'true black-box construction'
			! stop
		end if
		
		t1 = OMP_get_wtime()	
		write(*,*) "MVP-based HODLR construction......"		
		rankmax = 300
		call HODLR_MVP(rankmax,Memory,error)
		t2 = OMP_get_wtime()  
		write(*,*) "MVP-based HODLR construction finished",t2-t1, 'secnds. Error: ', error	
		
		if(fullmatflag==1)then	
			deallocate(matZ_glo)		
		end if	
	

	end if	
	
	
	
    ! write(*,*) "Cascading factorizing......"
    ! call cascading_factorizing(tolerance)
    ! write(*,*) "Cascading factorizing finished"
    ! write(*,*) "    "	


    ! write(*,*) "EM_calculating......"
    ! call EM_calculating()
    ! write(*,*) "EM_calculating finished"
    ! write(*,*) "    "	
	
	

    write(*,*) "-------------------------------program end-------------------------------------"

    ! ! ! ! pause

end PROGRAM MLMDA_DIRECT_SOLVER_3D_CFIE
