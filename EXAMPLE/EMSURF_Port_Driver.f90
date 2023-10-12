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
!> @brief This is an example that solves a 3D IE system for a cavity with ports in electromagnetics.
!> @details Note that instead of the use of precision dependent subroutine/module/type names "z_", one can also use the following \n
!> #define DAT 0 \n
!> #include "zButterflyPACK_config.fi" \n
!> which will macro replace precision-independent subroutine/module/type names "X" with "z_X" defined in SRC_DOUBLECOMLEX with double-complex precision

! This exmple works with double-complex precision data
PROGRAM ButterflyPACK_IE_3D
    use z_BPACK_DEFS
	use EMSURF_PORT_MODULE

	use z_BPACK_structure
	use z_BPACK_factor
	use z_BPACK_constr
	use z_BPACK_Solve_Mul
#ifdef HAVE_OPENMP
	use omp_lib
#endif
	use z_MISC_Utilities
	use z_BPACK_utilities
    implicit none

	! include "mkl_vml.fi"

    real(kind=8) para
    real(kind=8) tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj,pp,nn1
    integer i,j,k, threads_num
	integer seed_myid(50)
	integer times(8)
	real(kind=8) t1,t2,x,y,z,r,kc,betanm,norm1,normi,maxnorm
	complex(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:)

	character(len=:),allocatable  :: string
	character(len=1024)  :: strings,strings1
	character(len=6)  :: info_env
	integer :: length,edge
	integer :: ierr
	integer*8 oldmode,newmode
	type(z_Hoption)::option_post,option_sh,option_A,option_B
	type(z_Hstat)::stats_post,stats_sh,stats_A,stats_B
	type(z_mesh)::msh_post,msh_sh,msh_A,msh_B
	type(z_Bmatrix)::bmat_post, bmat_sh,bmat_A,bmat_B
	type(z_kernelquant)::ker_post, ker_sh,ker_A,ker_B
	type(quant_EMSURF),target::quant
	type(z_proctree)::ptree_post, ptree_sh,ptree_A,ptree_B
	integer,allocatable:: groupmembers(:)
	integer nmpi
	real(kind=8),allocatable::xyz(:,:)
	integer,allocatable::Permutation(:),order(:)
	integer Nunk_loc

	integer maxn, maxnev, maxncv, ldv, ith
	integer iparam(11), ipntr(14)
    logical,allocatable:: select(:)
    complex(kind=8),allocatable:: ax(:), mx(:), d(:), v(:,:), workd(:), workev(:), resid(:), workl(:), eigval(:), xloc(:,:), bloc(:,:), E_normal(:,:), eigvec_glo(:), eigvec_glo_ref(:)
    real(kind=8),allocatable:: rwork(:), rd(:,:),norms(:),norms_tmp(:),eigvec_glo_real_ref(:)

    character bmattype
	character(len=10) which
    integer ido, n, nev, ncv, lworkl, info, nconv, maxitr, ishfts, mode
    complex(kind=8) sigma
    real(kind=8) tol, retval(2), offset, rtemp1, rtemp2
    real(kind=8) dtheta,theta,phi,rcs
    real(kind=8) freq0,val0,valmin
    real(kind=8):: norm_thresh=1000d0
    real(kind=8):: dotproduct_thresh=0.4d0
	complex(kind=8) ctemp_loc,ctemp_1,ctemp
	logical rvec
	real(kind=8),external :: pdznorm2, dlapy2
	character(len=1024)  :: substring,substring1,substring2,substring3,substring4
	integer v_major,v_minor,v_bugfix
	integer Nfreq,Nmode,nth_mode,Nx,Ny

	integer nargs,flag
	integer parent,io

	logical :: exist, f_exist
	integer,parameter::Nx0=101,Ny0=101
	real(kind=8):: x0(Nx0),y0(Ny0),u(Nx0*Ny0),xx1(1000),yy1(1000), dx, v1(1000),vref(1000)
	character(len=1024)  :: string1,string2,string3,filename

# 125 "EMSURF_Port_Driver.f90"






	write(substring2 , *) 'pillbox'

	! nmpi and groupmembers should be provided by the user
	call MPI_Init(ierr)
	call MPI_COMM_GET_PARENT(parent, ierr) ! YL: this is needed if this function is spawned by a master process
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo

	call z_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree_A)
	deallocate(groupmembers)


	if(ptree_A%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "ButterflyPACK_IE_3D"
	call z_BPACK_GetVersionNumber(v_major,v_minor,v_bugfix)
	write(*,'(A23,I1,A1,I1,A1,I1,A1)') " ButterflyPACK Version:",v_major,".",v_minor,".",v_bugfix
    write(*,*) "   "
	endif

	!**** initialize stats and option
	call z_InitStat(stats_sh)
	call z_SetDefaultOptions(option_A)


	!**** intialize the user-defined derived type quant
	! compute the quadrature rules
    quant%integral_points=6
    allocate (quant%ng1(quant%integral_points), quant%ng2(quant%integral_points), quant%ng3(quant%integral_points), quant%gauss_w(quant%integral_points))
    call gauss_points(quant)

    !*************************input******************************
	quant%DATA_DIR='../EXAMPLE/EM3D_DATA/sphere_2300'

	quant%postprocess=1
	quant%mesh_normal=1
	quant%scaling=1d0
	quant%wavelength=2.0
	quant%freq=1/quant%wavelength/sqrt(BPACK_mu0*BPACK_eps0)
	quant%RCS_static=2
    quant%RCS_Nsample=1000
	quant%CFIE_alpha=1

	!**** default parameters for the eigen solvers
	quant%CMmode=0
	quant%SI=0
	quant%shift=(0d0, 0d0)
	quant%nev=1
	quant%tol_eig=1d-13
	quant%which='LM'


	option_A%ErrSol=1
	option_A%format=  HMAT!  HODLR !
	option_A%near_para=2.01d0
	option_A%verbosity=1
	option_A%ILU=0
	option_A%forwardN15flag=0
	option_A%LRlevel=100
	option_A%tol_itersol=1d-5
	option_A%sample_para=4d0

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
						if(trim(strings)=='--data_dir')then
							quant%data_dir=trim(strings1)
						else if	(trim(strings)=='--wavelength')then
							read(strings1,*)quant%wavelength
							quant%freq=1/quant%wavelength/sqrt(BPACK_mu0*BPACK_eps0)
						else if (trim(strings)=='--scaling')then
							read(strings1,*)quant%scaling
						else if (trim(strings)=='--norm_thresh')then
							read(strings1,*)norm_thresh
						else if (trim(strings)=='--freq')then
							read(strings1,*)quant%freq
							quant%wavelength=1/quant%freq/sqrt(BPACK_mu0*BPACK_eps0)
						else if	(trim(strings)=='--cmmode')then
							read(strings1,*)quant%CMmode
						else if	(trim(strings)=='--si')then
							read(strings1,*)quant%SI
						else if	(trim(strings)=='--postprocess')then
							read(strings1,*)quant%postprocess
						else if	(trim(strings)=='--nev')then
							read(strings1,*)quant%nev
						else if	(trim(strings)=='--tol_eig')then
							read(strings1,*)quant%tol_eig
						else if	(trim(strings)=='--model')then
							substring2=trim(strings1)
						else if	(trim(strings)=='--which')then
							quant%which=trim(strings1)
						else if	(trim(strings)=='--noport')then
							read(strings1,*)quant%noport
						else if	(trim(strings)=='--mesh_normal')then
							read(strings1,*)quant%mesh_normal
						else
							if(ptree_A%MyID==Main_ID)write(*,*)'ignoring unknown quant: ', trim(strings)
						endif
					else
						flag=0
					endif
				else
					flag=0
				endif
			enddo
		else if(trim(strings)=='-option')then ! options of ButterflyPACK
			call z_ReadOption(option_A,ptree_A,ii)
		else
			if(ptree_A%MyID==Main_ID)write(*,*)'ignoring unknown argument: ',trim(strings)
			ii=ii+1
		endif
	enddo

	call z_PrintOptions(option_A,ptree_A)


    quant%wavenum=2*BPACK_pi/quant%wavelength
	! option_A%touch_para = 3* quant%minedgelength


	!!!!!!! No port, the results should be identical as ie3deigen
	quant%Nport=0


	if(index(substring2,'rect_waveguide')>0)then
		!!!!!!! rectangular waveguide ports
		if(quant%noport==0)then
			quant%Nport=2
			allocate(quant%ports(quant%Nport))
			quant%ports(1)%origin=(/0.2d0,-0.1d0,0d0/)
			quant%ports(1)%x=(/-1d0,0d0,0d0/)
			quant%ports(1)%y=(/0d0,1d0,0d0/)
			quant%ports(1)%z=(/0d0,0d0,-1d0/)
			quant%ports(1)%type=1
			quant%ports(1)%a=0.4d0
			quant%ports(1)%b=0.2d0
			quant%ports(2)%origin=(/-0.2d0,-0.1d0,0.8d0/)
			quant%ports(2)%x=(/1d0,0d0,0d0/)
			quant%ports(2)%y=(/0d0,1d0,0d0/)
			quant%ports(2)%z=(/0d0,0d0,1d0/)
			quant%ports(2)%type=1
			quant%ports(2)%a=0.4d0
			quant%ports(2)%b=0.2d0
		endif
	endif


	if(index(substring2,'pillbox')>0)then
	!!!!!!! pillbox ports
	if(quant%noport==0)then
		quant%Nport=2
		allocate(quant%ports(quant%Nport))
		quant%ports(1)%origin=(/0d0,0d0,0d0/)
		quant%ports(1)%z=(/0d0,0d0,-1d0/)
		quant%ports(1)%x=(/-1d0,0d0,0d0/)
		quant%ports(1)%R=0.1
		quant%ports(1)%type=0
		quant%ports(2)%origin=(/0d0,0d0,0.1d0/)
		quant%ports(2)%z=(/0d0,0d0,1d0/)
		quant%ports(2)%x=(/1d0,0d0,0d0/)
		quant%ports(2)%R=0.1
		quant%ports(2)%type=0
	endif

	quant%Nobs=1000
	allocate(quant%obs_points(3,quant%Nobs))
	allocate(quant%obs_Efields(3,quant%Nobs))
	offset=0.001
	do ii=1,quant%Nobs
		quant%obs_points(1,ii)=0
		quant%obs_points(2,ii)=0
		quant%obs_points(3,ii)=(0.1-2*offset)/(quant%Nobs-1)*(ii-1)+offset
	enddo
	endif



	if(index(substring2,'cavity_wakefield')>0)then
	!!!!!!!!! cavity wakefield
	if(quant%noport==0)then
		quant%Nport=2
		allocate(quant%ports(quant%Nport))
		quant%ports(1)%origin=(/0d0,0d0,-0.1995d0/)
		quant%ports(1)%z=(/0d0,0d0,-1d0/)
		quant%ports(1)%x=(/1d0,0d0,0d0/)
		quant%ports(1)%R=0.037
		quant%ports(1)%type=0
		quant%ports(2)%origin=(/0d0,0d0,0.1995d0/)
		quant%ports(2)%z=(/0d0,0d0,1d0/)
		quant%ports(2)%x=(/1d0,0d0,0d0/)
		quant%ports(2)%R=0.037
		quant%ports(2)%type=0
	endif
	endif


	if(index(substring2,'cavity_rec')>0)then
	!!!!!!! cavity with 2 rectangular dumping ports and 2 circular beam ports
	if(quant%noport==0)then
		quant%Nport=2
		allocate(quant%ports(quant%Nport))
		quant%ports(1)%origin=(/0.035d0,0.2d0,0.01d0/)
		quant%ports(1)%x=(/-1d0,0d0,0d0/)
		quant%ports(1)%y=(/0d0,0d0,1d0/)
		quant%ports(1)%z=(/0d0,1d0,0d0/)
		quant%ports(1)%type=1
		quant%ports(1)%a=0.07d0
		quant%ports(1)%b=0.02d0

		quant%ports(2)%origin=(/-0.2d0,0.035d0,-0.03d0/)
		quant%ports(2)%x=(/0d0,-1d0,0d0/)
		quant%ports(2)%y=(/0d0,0d0,1d0/)
		quant%ports(2)%z=(/-1d0,0d0,0d0/)
		quant%ports(2)%type=1
		quant%ports(2)%a=0.07d0
		quant%ports(2)%b=0.02d0




		! quant%ports(1)%origin=(/0.035d0,0.2d0,0.01d0/)
		! quant%ports(1)%x=(/-1d0,0d0,0d0/)
		! quant%ports(1)%y=(/0d0,0d0,1d0/)
		! quant%ports(1)%z=(/0d0,1d0,0d0/)
		! quant%ports(1)%type=2
		! quant%ports(1)%a=0.07d0
		! quant%ports(1)%b=0.02d0
		! quant%ports(1)%Nx_arbi = 1400
		! quant%ports(1)%Ny_arbi = 400
		! quant%ports(1)%string_arbi="cavity_rec_port"
		! quant%ports(1)%nmode_arbi=4

		! quant%ports(2)%origin=(/-0.2d0,0.035d0,-0.03d0/)
		! quant%ports(2)%x=(/0d0,-1d0,0d0/)
		! quant%ports(2)%y=(/0d0,0d0,1d0/)
		! quant%ports(2)%z=(/-1d0,0d0,0d0/)
		! quant%ports(2)%type=2
		! quant%ports(2)%a=0.07d0
		! quant%ports(2)%b=0.02d0
		! quant%ports(2)%Nx_arbi = 1400
		! quant%ports(2)%Ny_arbi = 400
		! quant%ports(2)%string_arbi="cavity_rec_port"
		! quant%ports(2)%nmode_arbi=4




		! quant%ports(3)%origin=(/0d0,0d0,-0.1445476695d0/)
		! quant%ports(3)%z=(/0d0,0d0,-1d0/)
		! quant%ports(3)%x=(/1d0,0d0,0d0/)
		! quant%ports(3)%R=0.025
		! quant%ports(3)%type=0
		! quant%ports(4)%origin=(/0d0,0d0,0.1445476695d0/)
		! quant%ports(4)%z=(/0d0,0d0,1d0/)
		! quant%ports(4)%x=(/1d0,0d0,0d0/)
		! quant%ports(4)%R=0.025
		! quant%ports(4)%type=0
	endif

	quant%Nobs=1000
	allocate(quant%obs_points(3,quant%Nobs))
	allocate(quant%obs_Efields(3,quant%Nobs))
	offset=0.001
	do ii=1,quant%Nobs
		quant%obs_points(1,ii)=0
		quant%obs_points(2,ii)=0
		quant%obs_points(3,ii)=(0.289-2*offset)/(quant%Nobs-1)*(ii-1)+offset-0.144
	enddo
	endif


	if(index(substring2,'cavity_no_wg')>0)then
		!!!!!!! cavity with 2 circular beam ports
		quant%Nobs=1000
		allocate(quant%obs_points(3,quant%Nobs))
		allocate(quant%obs_Efields(3,quant%Nobs))
		offset=0.001
		do ii=1,quant%Nobs
			quant%obs_points(1,ii)=0
			quant%obs_points(2,ii)=0
			quant%obs_points(3,ii)=(0.289-2*offset)/(quant%Nobs-1)*(ii-1)+offset-0.144
		enddo
	endif


	if(index(substring2,'rfq_mirror')>0)then
	!!!!!!!! RFQ no port
	quant%Nport=0
	quant%Nobs=1000
	allocate(quant%obs_points(3,quant%Nobs))
	allocate(quant%obs_Efields(3,quant%Nobs))
	offset=0.001
	do ii=1,quant%Nobs
		quant%obs_points(1,ii)=0
		quant%obs_points(2,ii)=0
		quant%obs_points(3,ii)=(2.115227780-2*offset)/(quant%Nobs-1)*(ii-1)+offset - 1.05761389
	enddo
	endif




	if(ptree_A%MyID==Main_ID .and. quant%Nport>0)write(*,*)'Port No.  TypeTETM  n  m  kc  k  eta_nm'
	quant%Nunk_waveguidemode=0
	do pp=1,quant%Nport
		if(quant%ports(pp)%type==0)then
			call z_rrcurl(quant%ports(pp)%z,quant%ports(pp)%x,quant%ports(pp)%y)
			quant%ports(pp)%mmax=2
			quant%ports(pp)%nmax=2

			do nn=0,quant%ports(pp)%nmax
			do mm=1,quant%ports(pp)%mmax
				kc = r_TE_nm(nn+1,mm)/quant%ports(pp)%R
				if(quant%wavenum > kc)then
					quant%Nunk_waveguidemode = quant%Nunk_waveguidemode+2 ! 2 accounts for two polarizations for circular waveguides
				endif
				kc = r_TM_nm(nn+1,mm)/quant%ports(pp)%R
				if(quant%wavenum > kc)then
					quant%Nunk_waveguidemode = quant%Nunk_waveguidemode+2 ! 2 accounts for two polarizations for circular waveguides
				endif
			enddo
			enddo
		else if(quant%ports(pp)%type==1)then
			quant%ports(pp)%mmax=2
			quant%ports(pp)%nmax=2

			do nn=0,quant%ports(pp)%nmax
				do mm=0,quant%ports(pp)%mmax
					if(nn>0 .or. mm>0)then ! the lowest TE mode is 01 or 10
						kc = sqrt((nn*BPACK_pi/quant%ports(pp)%a)**2d0+(mm*BPACK_pi/quant%ports(pp)%b)**2d0)
						if(quant%wavenum > kc)then
							quant%Nunk_waveguidemode = quant%Nunk_waveguidemode+1
						endif
					endif
				enddo
			enddo
			do nn=0,quant%ports(pp)%nmax
				do mm=0,quant%ports(pp)%mmax
					if(nn>0 .and. mm>0)then ! the lowest TM mode is 11
						kc = sqrt((nn*BPACK_pi/quant%ports(pp)%a)**2d0+(mm*BPACK_pi/quant%ports(pp)%b)**2d0)
						if(quant%wavenum > kc)then
							quant%Nunk_waveguidemode = quant%Nunk_waveguidemode+1
						endif
					endif
				enddo
			enddo
		else if(quant%ports(pp)%type==2)then
			do nn=1,quant%ports(pp)%nmode_arbi
				write(string1 , *) quant%ports(pp)%Nx_arbi
				write(string2 , *) quant%ports(pp)%Ny_arbi
				write(string3 , *) nn
				open(22, file=trim(adjustl(quant%ports(pp)%string_arbi))//"_Nx_"//trim(adjustl(string1))//"_Ny_"//trim(adjustl(string2))//"_mode_"//trim(adjustl(string3))//".txt", status="old", action="read")
				read(22,*)quant%ports(pp)%kc_arbi(nn),dx,quant%ports(pp)%TETM_arbi(nn),quant%ports(pp)%A_n_arbi(nn)
				if(quant%wavenum > quant%ports(pp)%kc_arbi(nn))then
					quant%Nunk_waveguidemode = quant%Nunk_waveguidemode+1
				endif
				close(22)
			enddo
		endif
	enddo

	if(quant%Nunk_waveguidemode>0)then
		allocate(quant%waveguidemodes(quant%Nunk_waveguidemode,5))
		quant%waveguidemodes=0
	endif

	quant%Nunk_waveguidemode=0
	do pp=1,quant%Nport
		if(quant%ports(pp)%type==0)then
			call z_rrcurl(quant%ports(pp)%z,quant%ports(pp)%x,quant%ports(pp)%y)
			do nn=0,quant%ports(pp)%nmax
			do mm=1,quant%ports(pp)%mmax
				kc = r_TE_nm(nn+1,mm)/quant%ports(pp)%R
				quant%ports(pp)%A_TE_nm(nn+1,mm)=A_TE_nm_cir(nn+1,mm)
				quant%ports(pp)%impedance_TE_nm(nn+1,mm)=BPACK_Bigvalue
				if(quant%wavenum > kc)then
					quant%Nunk_waveguidemode = quant%Nunk_waveguidemode+1
					quant%waveguidemodes(quant%Nunk_waveguidemode,1)=pp      ! port ID
					quant%waveguidemodes(quant%Nunk_waveguidemode,2)=nn    ! nn
					quant%waveguidemodes(quant%Nunk_waveguidemode,3)=mm    ! mm
					quant%waveguidemodes(quant%Nunk_waveguidemode,4)=1       ! TETM
					quant%waveguidemodes(quant%Nunk_waveguidemode,5)=1       ! polarization
					quant%Nunk_waveguidemode = quant%Nunk_waveguidemode+1
					quant%waveguidemodes(quant%Nunk_waveguidemode,1)=pp      ! port ID
					quant%waveguidemodes(quant%Nunk_waveguidemode,2)=nn    ! nn
					quant%waveguidemodes(quant%Nunk_waveguidemode,3)=mm    ! mm
					quant%waveguidemodes(quant%Nunk_waveguidemode,4)=1       ! TETM
					quant%waveguidemodes(quant%Nunk_waveguidemode,5)=2       ! polarization
					betanm=sqrt(quant%wavenum**2d0-kc**2d0)
					quant%ports(pp)%impedance_TE_nm(nn+1,mm)=quant%wavenum*BPACK_impedence0/betanm
					if(ptree_A%MyID==Main_ID)write(*,*)pp,'CIR','TE',nn,mm,kc,quant%wavenum,quant%ports(pp)%impedance_TE_nm(nn+1,mm)
				endif
				kc = r_TM_nm(nn+1,mm)/quant%ports(pp)%R
				quant%ports(pp)%A_TM_nm(nn+1,mm)=A_TM_nm_cir(nn+1,mm)
				quant%ports(pp)%impedance_TM_nm(nn+1,mm)=BPACK_Bigvalue
				if(quant%wavenum > kc)then
					quant%Nunk_waveguidemode = quant%Nunk_waveguidemode+1
					quant%waveguidemodes(quant%Nunk_waveguidemode,1)=pp      ! port ID
					quant%waveguidemodes(quant%Nunk_waveguidemode,2)=nn    ! nn
					quant%waveguidemodes(quant%Nunk_waveguidemode,3)=mm    ! mm
					quant%waveguidemodes(quant%Nunk_waveguidemode,4)=2       ! TETM
					quant%waveguidemodes(quant%Nunk_waveguidemode,5)=1       ! polarization
					quant%Nunk_waveguidemode = quant%Nunk_waveguidemode+1
					quant%waveguidemodes(quant%Nunk_waveguidemode,1)=pp      ! port ID
					quant%waveguidemodes(quant%Nunk_waveguidemode,2)=nn    ! nn
					quant%waveguidemodes(quant%Nunk_waveguidemode,3)=mm    ! mm
					quant%waveguidemodes(quant%Nunk_waveguidemode,4)=2       ! TETM
					quant%waveguidemodes(quant%Nunk_waveguidemode,5)=2       ! polarization
					betanm=sqrt(quant%wavenum**2d0-kc**2d0)
					quant%ports(pp)%impedance_TM_nm(nn+1,mm)=BPACK_impedence0*betanm/quant%wavenum
					if(ptree_A%MyID==Main_ID)write(*,*)pp,'CIR','TM',nn,mm,kc,quant%wavenum,quant%ports(pp)%impedance_TM_nm(nn+1,mm)
				endif
			enddo
			enddo
		else if(quant%ports(pp)%type==1)then
			do nn=0,quant%ports(pp)%nmax
				do mm=0,quant%ports(pp)%mmax
					quant%ports(pp)%A_TE_nm(nn+1,mm+1)=0
					quant%ports(pp)%impedance_TE_nm(nn+1,mm+1)=BPACK_Bigvalue
					if(nn>0 .or. mm>0)then ! the lowest TE mode is 01 or 10
						kc = sqrt((nn*BPACK_pi/quant%ports(pp)%a)**2d0+(mm*BPACK_pi/quant%ports(pp)%b)**2d0)
						quant%ports(pp)%A_TE_nm(nn+1,mm+1)=1d0/sqrt(quant%ports(pp)%a/quant%ports(pp)%b*mm**2d0*A_nm_rec(nn+1,mm+1) + quant%ports(pp)%b/quant%ports(pp)%a*nn**2d0*B_nm_rec(nn+1,mm+1))
						if(quant%wavenum > kc)then
							quant%Nunk_waveguidemode = quant%Nunk_waveguidemode+1
							quant%waveguidemodes(quant%Nunk_waveguidemode,1)=pp
							quant%waveguidemodes(quant%Nunk_waveguidemode,2)=nn
							quant%waveguidemodes(quant%Nunk_waveguidemode,3)=mm
							quant%waveguidemodes(quant%Nunk_waveguidemode,4)=1
							quant%waveguidemodes(quant%Nunk_waveguidemode,5)=1
							betanm=sqrt(quant%wavenum**2d0-kc**2d0)
							quant%ports(pp)%impedance_TE_nm(nn+1,mm+1)=quant%wavenum*BPACK_impedence0/betanm
							if(ptree_A%MyID==Main_ID)write(*,*)pp,'RECT','TE',nn,mm,kc,quant%wavenum,quant%ports(pp)%impedance_TE_nm(nn+1,mm+1)
						endif
					endif
				enddo
			enddo

			do nn=0,quant%ports(pp)%nmax
				do mm=0,quant%ports(pp)%mmax
					quant%ports(pp)%A_TM_nm(nn+1,mm+1)=0
					quant%ports(pp)%impedance_TM_nm(nn+1,mm+1)=BPACK_Bigvalue
					if(nn>0 .and. mm>0)then ! the lowest TM mode is 11
						kc = sqrt((nn*BPACK_pi/quant%ports(pp)%a)**2d0+(mm*BPACK_pi/quant%ports(pp)%b)**2d0)
						quant%ports(pp)%A_TM_nm(nn+1,mm+1)=1d0/sqrt(quant%ports(pp)%a/quant%ports(pp)%b*mm**2d0*B_nm_rec(nn+1,mm+1) + quant%ports(pp)%b/quant%ports(pp)%a*nn**2d0*A_nm_rec(nn+1,mm+1))
						if(quant%wavenum > kc)then
							quant%Nunk_waveguidemode = quant%Nunk_waveguidemode+1
							quant%waveguidemodes(quant%Nunk_waveguidemode,1)=pp
							quant%waveguidemodes(quant%Nunk_waveguidemode,2)=nn
							quant%waveguidemodes(quant%Nunk_waveguidemode,3)=mm
							quant%waveguidemodes(quant%Nunk_waveguidemode,4)=2
							quant%waveguidemodes(quant%Nunk_waveguidemode,5)=1
							betanm=sqrt(quant%wavenum**2d0-kc**2d0)
							quant%ports(pp)%impedance_TM_nm(nn+1,mm+1)=BPACK_impedence0*betanm/quant%wavenum
							if(ptree_A%MyID==Main_ID)write(*,*)pp,'RECT','TM',nn,mm,kc,quant%wavenum,quant%ports(pp)%impedance_TM_nm(nn+1,mm+1)
						endif
					endif
				enddo
			enddo
		endif

	enddo




   !***********************************************************************
	if(ptree_A%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'frequency:',quant%freq
   write (*,*) 'wavelength:',quant%wavelength
   write (*,*) ''
	endif
   !***********************************************************************


	!***********************************************************************
	!**** construct compressed A
	!**** geometry generalization and discretization
	call geo_modeling_SURF(quant,ptree_A%Comm,quant%DATA_DIR)
	!**** register the user-defined function and type in ker
	ker_A%QuantApp => quant
	ker_A%FuncZmn => Zelem_EMSURF
	!**** initialization of the construction phase
	t1 = MPI_Wtime()
	allocate(xyz(3,quant%Nunk))
	do ii=1, quant%Nunk_int+quant%Nunk_port
		xyz(:,ii) = quant%xyz(:,quant%maxnode+ii)
	enddo
	! use port center as the geometry for the mode unknowns
	do ii=quant%Nunk_int + quant%Nunk_port+1,quant%Nunk
		pp = quant%waveguidemodes(ii-quant%Nunk_int - quant%Nunk_port,1)
		xyz(:,ii)=quant%ports(pp)%origin
	enddo

    allocate(Permutation(quant%Nunk))
	call z_BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat_A,option_A,stats_A,msh_A,ker_A,ptree_A,Coordinates=xyz)
	deallocate(Permutation) ! caller can use this permutation vector if needed
	deallocate(xyz)
	t2 = MPI_Wtime()
	call nxK_waveguidePrecompute(option_A,msh_A,quant,ptree_A,stats_A)
	!**** computation of the construction phase
    call z_BPACK_construction_Element(bmat_A,option_A,stats_A,msh_A,ker_A,ptree_A)
	!**** print statistics

	if(.not. (quant%SI==1 .and. abs(quant%shift)<1e-14))then
		call z_PrintStat(stats_A,ptree_A)
	endif


	!***********************************************************************
	!**** construct compressed A - sigma I or  A - sigma real(A)
	if(quant%SI==1)then
		call z_CopyOptions(option_A,option_sh)
		!**** create the process tree, can use larger number of mpis if needed
		allocate(groupmembers(nmpi))
		do ii=1,nmpi
			groupmembers(ii)=(ii-1)
		enddo
		call z_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree_sh)
		deallocate(groupmembers)

		if(abs(quant%shift)<1e-14)then  ! if zero shift, no need to compress the shifted matrix
			call z_CopyMesh(msh_A,msh_sh)
			call z_CopyStat(stats_A,stats_sh)
			call z_BPACK_copy(bmat_A,bmat_sh,ptree_sh)
		else
			!**** initialize stats_sh and option
			call z_InitStat(stats_sh)
			!**** register the user-defined function and type in ker
			ker_sh%QuantApp => quant
			ker_sh%FuncZmn => Zelem_EMSURF_Shifted
			!**** initialization of the construction phase
			t1 = MPI_Wtime()
			allocate(xyz(3,quant%Nunk))
			do ii=1, quant%Nunk_int+quant%Nunk_port
				xyz(:,ii) = quant%xyz(:,quant%maxnode+ii)
			enddo
			! use port center as the geometry for the mode unknowns
			do ii=quant%Nunk_int + quant%Nunk_port+1,quant%Nunk
				pp = quant%waveguidemodes(ii-quant%Nunk_int - quant%Nunk_port,1)
				xyz(:,ii)=quant%ports(pp)%origin
			enddo

			allocate(Permutation(quant%Nunk))
			call z_BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat_sh,option_sh,stats_sh,msh_sh,ker_sh,ptree_sh,Coordinates=xyz)
			deallocate(Permutation) ! caller can use this permutation vector if needed
			deallocate(xyz)
			t2 = MPI_Wtime()
			!**** computation of the construction phase
			call z_BPACK_construction_Element(bmat_sh,option_sh,stats_sh,msh_sh,ker_sh,ptree_sh)
		endif


		!**** factorization phase
		call z_BPACK_Factorization(bmat_sh,option_sh,stats_sh,ptree_sh,msh_sh)
		! !**** solve phase
		! call EM_solve_SURF(bmat_sh,option_sh,msh_sh,quant,ptree_sh,stats_sh)
		!**** print statistics
		call z_PrintStat(stats_sh,ptree_sh)
	endif

	!***********************************************************************
	allocate(xloc(Nunk_loc,quant%Nunk_waveguidemode))
	allocate(bloc(Nunk_loc,quant%Nunk_waveguidemode))
	call EM_solve_port_SURF(bmat_sh,option_sh,msh_sh,quant,ptree_sh,stats_sh,xloc,bloc)


	if(ptree_A%MyID==Main_ID)then
		write(*,*)'Postprocessing: '
	endif

	do ii=1,quant%Nunk_waveguidemode
		write(substring , *) ii
		write(substring1 , *) quant%freq
		call current_node_patch_mapping(0,trim(adjustl(substring2))//'_J_port_'//trim(adjustl(substring))//'_freq_'//trim(adjustl(substring1))//'.out',xloc(:,ii),msh_A,quant,ptree_A)
		call current_node_patch_mapping(1,trim(adjustl(substring2))//'_M_port_'//trim(adjustl(substring))//'_freq_'//trim(adjustl(substring1))//'.out',xloc(:,ii),msh_A,quant,ptree_A)
	enddo

	if(ptree_A%MyID==Main_ID .and. option_A%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	deallocate(xloc)

	if(quant%SI==1)then
		call z_delete_proctree(ptree_sh)
		call z_delete_Hstat(stats_sh)
		call z_delete_mesh(msh_sh)
		call z_delete_kernelquant(ker_sh)
		call z_BPACK_delete(bmat_sh)
	endif

	call z_delete_proctree(ptree_A)
	call z_delete_Hstat(stats_A)
	call z_delete_mesh(msh_A)
	call z_delete_kernelquant(ker_A)
	call z_BPACK_delete(bmat_A)

	!**** deletion of quantities
	call delete_quant_EMSURF(quant)


	call z_blacs_exit_wrp(1)
	call MPI_Finalize(ierr)

    ! ! ! ! pause

end PROGRAM ButterflyPACK_IE_3D







