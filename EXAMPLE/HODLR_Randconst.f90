PROGRAM MLMDA_DIRECT_SOLVER_3D_CFIE

    use MODULE_FILE
	! use geometry_model
	use H_structure
	use cascading_factorization
	use matrices_fill
	use omp_lib
	use Bplus_compress_forward
	use HODLR_randomMVP
    implicit none

    real*8 para,error
    real*8 tolerance
    integer Primary_block, nn, mm
    integer i,j,k, threads_num,ii,jj
	integer seed_myid(12)
	integer times(8)	
	real*8 t1,t2,t3,t4,x,y,z,r,theta,phi,tmp(3),Memory
	complex(kind=8),allocatable:: InputVec(:)
	complex(kind=8):: ctemp
	integer Ntunnel,kk,black_step,rankmax
	complex(kind=8),allocatable::Vout1(:,:),Vout2(:,:),Vin(:,:)
	character(len=1024)  :: strings
	type(Hoption):: option
	type(Hstat)::stats	
	type(mesh)::msh	
	type(kernelquant)::ker
	integer:: explicitflag
	type(hobf)::ho_bf,ho_bf_copy
	integer Nin1,Nout1,Nin2,Nout2	
	type(proctree)::ptree
	integer,allocatable:: groupmembers(:)	
	
	
	
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
	
    !real scale


    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "MLMDA_DIRECT_SOLVER_3D_CFIE_NEW"
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

	msh%Origins=(/0d0,0d0,0d0/)
    
    msh%integral_points=6
    allocate (msh%ng1(msh%integral_points), msh%ng2(msh%integral_points), msh%ng3(msh%integral_points), msh%gauss_w(msh%integral_points))
    call gauss_points(msh)


     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	CALL getarg(1, DATA_DIR)
	
	
	ker%Kernel = FULL

	option%Nmin_leaf=100
	! Refined_level=0
	para=0.001
	! levelpara_control=0
	! ACA_tolerance_forward=1d-4
	option%tol_SVD=1d-4
	! SVD_tolerance_factor=1d-4
	option%tol_Rdetect=1d-4 !3d-5
    ! Preset_level_butterfly=0
	msh%scaling=1d0
	ker%wavelength=0.25
	! Discret=0.05
	ker%RCS_static=2
    ker%RCS_Nsample=1000
    ! Optimizing_forward=0
    ! Fast_inverse=0
    ! Add_method_of_base_level=2
    ker%rank_approximate_para1=6.0
    ker%rank_approximate_para2=6.0
    ker%rank_approximate_para3=6.0
	option%tol_LS=1d-10
	! tfqmr_tolerance=1d-6
	option%tol_itersol=3d-3
	option%N_iter=1000
	option%tol_rand=1d0
	! up_tolerance=1
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
	option%LnoBP=600
	option%TwoLayerOnly=1
	option%LRlevel=100
	ker%CFIE_alpha=1
	explicitflag=0
	! fullmatflag=1


	
	
	!Nmin_leaf=250
	!para=0.01
	!tolerance=0.001
	!msh%scaling=1.
	!alpha=0.5
    !ker%wavelength=2.
    ! tolerance=ACA_tolerance_forward
    ! call OMP_set_dynamic(.true.)
    !call OMP_set_nested(.true.)

    !*********************************************************

    ker%omiga=2*pi/ker%wavelength/sqrt(mu0*eps0)
    ker%wavenum=2*pi/ker%wavelength

	Ntunnel = 15600
	msh%Nunk = 3720
	Nin1 = 320*2 
	Nout1 = 610*2
	Nin2 = 320*2 
	Nout2 = 610*2	
	call assert(Nin1+Nout1+Nin2+Nout2==msh%Nunk,'The two surfaces have mismatched number of unknowns')
	
	! predefine the first three levels of tree due to the physical meanings
	allocate(basis_group_pre(7,2))
	basis_group_pre(1,1) = 1
	basis_group_pre(1,2) = msh%Nunk
	basis_group_pre(2,1) = 1
	basis_group_pre(2,2) = Nin1+Nout1	
	basis_group_pre(3,1) = Nin1+Nout1+1
	basis_group_pre(3,2) = msh%Nunk	
	basis_group_pre(4,1) = 1
	basis_group_pre(4,2) = Nin1	
	basis_group_pre(5,1) = 1+Nin1
	basis_group_pre(5,2) = Nout1+Nin1
	basis_group_pre(6,1) = 1+Nin1+Nout1
	basis_group_pre(6,2) = Nin2+Nin1+Nout1
	basis_group_pre(7,1) = 1+Nin2+Nin1+Nout1
	basis_group_pre(7,2) = Nout2+Nin2+Nin1+Nout1
	
	
	
	write(*,*)'Blackbox HODLR for scattering matrix compression'
	
	write(*,'(A10,I9,A11,I9)')' Ntunnel: ',Ntunnel,' Nsurface: ',msh%Nunk
	
	allocate(msh%xyz(3,msh%Nunk))
	allocate(msh%info_unk(0:0,msh%Nunk))
	do kk=1,msh%Nunk
		msh%info_unk(0,kk)=kk
	enddo
	!!!!!! for simplicity, I'm duplicating locations of each point

	open(unit=521,file=trim(DATA_DIR)//'/edge_cen.out',status='old')
	do kk=1,Ntunnel/2
		read(521,*) tmp(1),tmp(2),tmp(3)
	end do          
	do kk=1,msh%Nunk/2
		read(521,*) msh%xyz(1,2*kk-1),msh%xyz(2,2*kk-1),msh%xyz(3,2*kk-1)
		msh%xyz(:,2*kk) = msh%xyz(:,2*kk-1)
	end do
	close(521)
	
	
	
   ! !***********************************************************************
   ! write (*,*) ''
   ! write (*,*) 'CFIE computing'
   ! write (*,*) 'ker%wavelength:',ker%wavelength
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
    call H_matrix_structuring(ho_bf,para,option,msh)
	call BPlus_structuring(ho_bf,option,msh,ptree)
    write(*,*) "H_matrices formatting finished"
    write(*,*) "    "
	t2 = OMP_get_wtime()
	write(*,*)t2-t1,'secnds'	
	! stop 
	
    !pause
    


	t1 = OMP_get_wtime()
	write(*,*) "Generating fullmat ......"
	allocate(ker%matZ_glo(msh%Nunk,msh%Nunk))
	ker%matZ_glo = 0
	
	open(unit=888,file=trim(DATA_DIR)//'/fullmat.out',status='old')
	do ii=1,msh%Nunk
	do kk=1,msh%Nunk
		read(888,*)tmp(1),tmp(2)
		ker%matZ_glo(kk,ii) = cmplx(tmp(1),tmp(2),dp) 
	end do
	end do
	close(unit=888)
	
	write(*,*) "Generating fullmat finished"
	t2 = OMP_get_wtime()   
	write(*,*)t2-t1, 'secnds'	

	
	if(explicitflag ==1)then

		t1 = OMP_get_wtime()	
		write(*,*) "H_matrices filling......"
		call matrices_filling(ho_bf,option,stats,msh,ker,element_Zmn_FULL,ptree)
		if(option%precon/=DIRECT)then
			call copy_HOBF(ho_bf,ho_bf_copy)	
		end if		
		write(*,*) "H_matrices filling finished"
		write(*,*) "    "
		t2 = OMP_get_wtime()   
		write(*,*)t2-t1, 'secnds'			
	


		allocate(Vin(msh%Nunk,1))
		allocate(Vout1(msh%Nunk,1))
		allocate(Vout2(msh%Nunk,1))
		do ii=1,msh%Nunk
			Vin(ii,1) = random_complex_number()
		end do
		
		call MVM_Z_forward('N',msh%Nunk,1,1,ho_bf%Maxlevel+1,Vin,Vout1,ho_bf,ptree,stats)
		
		do ii=1,msh%Nunk
			ctemp = 0d0
			do jj=1,msh%Nunk
				ctemp = ctemp + ker%matZ_glo(new2old(ii),new2old(jj))*Vin(jj,1)
			end do
			Vout2(ii,1) = ctemp
		end do
		error = fnorm(Vout2-Vout1,msh%Nunk,1)/fnorm(Vout2,msh%Nunk,1)
		deallocate(Vin,Vout1,Vout2)
		
		write(*,*)error,'accuracy of construction'

		
	else if(explicitflag ==0)then

		t1 = OMP_get_wtime()	
		write(*,*) "MVP-based HODLR construction......"		
		rankmax = 300
		call HODLR_randomized(ho_bf,HODLR_MVP_randomized_Fullmat,msh%Nunk,rankmax,Memory,error,option,stats,ker,ptree)
		t2 = OMP_get_wtime()  
		write(*,*) "MVP-based HODLR construction finished",t2-t1, 'secnds. Error: ', error	

	end if	
	
	
	
    ! write(*,*) "Cascading factorizing......"
    ! call cascading_factorizing(ho_bf,option,stats,ptree)
    ! write(*,*) "Cascading factorizing finished"
    ! write(*,*) "    "	


    ! write(*,*) "EM_solve......"
    ! call EM_solve()
    ! write(*,*) "EM_solve finished"
    ! write(*,*) "    "	
	
	

    write(*,*) "-------------------------------program end-------------------------------------"

    ! ! ! ! pause

end PROGRAM MLMDA_DIRECT_SOLVER_3D_CFIE
