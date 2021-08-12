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
	use EMSURF_PORT_MODULE

	use BPACK_structure
	use BPACK_factor
	use BPACK_constr
	use BPACK_Solve_Mul
	use omp_lib
	use MISC_Utilities
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
	type(Hoption)::option_post,option_sh,option_A,option_B
	type(Hstat)::stats_post,stats_sh,stats_A,stats_B
	type(mesh)::msh_post,msh_sh,msh_A,msh_B
	type(Bmatrix)::bmat_post, bmat_sh,bmat_A,bmat_B
	type(kernelquant)::ker_post, ker_sh,ker_A,ker_B
	type(quant_EMSURF),target::quant
	type(proctree)::ptree_post, ptree_sh,ptree_A,ptree_B
	integer,allocatable:: groupmembers(:)
	integer nmpi
	real(kind=8),allocatable::xyz(:,:)
	integer,allocatable::Permutation(:),order(:)
	integer Nunk_loc

	integer maxn, maxnev, maxncv, ldv, ith
	integer iparam(11), ipntr(14)
    logical,allocatable:: select(:)
    complex(kind=8),allocatable:: ax(:), mx(:), d(:), v(:,:), workd(:), workev(:), resid(:), workl(:), eigval(:), eigvec(:,:),E_normal(:,:), eigvec_glo(:), eigvec_glo_ref(:)
    real(kind=8),allocatable:: rwork(:), rd(:,:),norms(:),norms_tmp(:),eigvec_glo_real_ref(:)

    character bmattype
	character(len=10) which
    integer ido, n, nx, nev, ncv, lworkl, info, nconv, maxitr, ishfts, mode
    complex(kind=8) sigma
    real(kind=8) tol, retval(2), offset, rtemp1, rtemp2
    real(kind=8) dtheta,theta,phi,rcs
    real(kind=8) freq0,val0,valmin
    real(kind=8):: norm_thresh=1000d0
    real(kind=8):: dotproduct_thresh=0.4d0
	complex(kind=8) ctemp_loc,ctemp_1,ctemp
	logical rvec
	real(kind=8),external :: pdznorm2, dlapy2
	character(len=1024)  :: substring,substring1,substring2
	integer v_major,v_minor,v_bugfix
	integer Nfreq,Nmode,nth_mode

	integer nargs,flag
	integer parent,io

	logical :: exist, f_exist


	write(substring2 , *) 'pillbox'

	! nmpi and groupmembers should be provided by the user
	call MPI_Init(ierr)
	call MPI_COMM_GET_PARENT(parent, ierr) ! YL: this is needed if this function is spawned by a master process
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo

	call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree_A)
	deallocate(groupmembers)


	if(ptree_A%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "ButterflyPACK_IE_3D"
	call BPACK_GetVersionNumber(v_major,v_minor,v_bugfix)
	write(*,'(A23,I1,A1,I1,A1,I1,A1)') " ButterflyPACK Version:",v_major,".",v_minor,".",v_bugfix
    write(*,*) "   "
	endif

	!**** initialize stats and option
	call InitStat(stats_sh)
	call SetDefaultOptions(option_A)


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
	quant%freq=1/quant%wavelength/sqrt(mu0*eps0)
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
							quant%freq=1/quant%wavelength/sqrt(mu0*eps0)
						else if (trim(strings)=='--scaling')then
							read(strings1,*)quant%scaling
						else if (trim(strings)=='--norm_thresh')then
							read(strings1,*)norm_thresh
						else if (trim(strings)=='--freq')then
							read(strings1,*)quant%freq
							quant%wavelength=1/quant%freq/sqrt(mu0*eps0)
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
			call ReadOption(option_A,ptree_A,ii)
		else
			if(ptree_A%MyID==Main_ID)write(*,*)'ignoring unknown argument: ',trim(strings)
			ii=ii+1
		endif
	enddo

	call PrintOptions(option_A,ptree_A)


    quant%wavenum=2*pi/quant%wavelength
	! option_A%touch_para = 3* quant%minedgelength


	!!!!!!! No port, the results should be identical as ie3deigen
	quant%Nport=0


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
		quant%Nport=4
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

		quant%ports(3)%origin=(/0d0,0d0,-0.1445476695d0/)
		quant%ports(3)%z=(/0d0,0d0,-1d0/)
		quant%ports(3)%x=(/1d0,0d0,0d0/)
		quant%ports(3)%R=0.025
		quant%ports(3)%type=0
		quant%ports(4)%origin=(/0d0,0d0,0.1445476695d0/)
		quant%ports(4)%z=(/0d0,0d0,1d0/)
		quant%ports(4)%x=(/1d0,0d0,0d0/)
		quant%ports(4)%R=0.025
		quant%ports(4)%type=0
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
	do pp=1,quant%Nport
		quant%ports(pp)%mmax=1
		quant%ports(pp)%nmax=1

		if(quant%ports(pp)%type==0)then
			call curl(quant%ports(pp)%z,quant%ports(pp)%x,quant%ports(pp)%y)
			quant%ports(pp)%mmax=1
			quant%ports(pp)%nmax=1

			do nn=0,quant%ports(pp)%nmax
			do mm=1,quant%ports(pp)%mmax
				kc = r_TE_nm(nn+1,mm)/quant%ports(pp)%R
				quant%ports(pp)%A_TE_nm(nn+1,mm)=A_TE_nm_cir(nn+1,mm)
				quant%ports(pp)%impedance_TE_nm(nn+1,mm)=Bigvalue
				if(quant%wavenum > kc)then
					betanm=sqrt(quant%wavenum**2d0-kc**2d0)
					quant%ports(pp)%impedance_TE_nm(nn+1,mm)=quant%wavenum*impedence0/betanm
					if(ptree_A%MyID==Main_ID)write(*,*)pp,'CIR','TE',nn,mm,kc,quant%wavenum,quant%ports(pp)%impedance_TE_nm(nn+1,mm)
				endif
				kc = r_TM_nm(nn+1,mm)/quant%ports(pp)%R
				quant%ports(pp)%A_TM_nm(nn+1,mm)=A_TM_nm_cir(nn+1,mm)
				quant%ports(pp)%impedance_TM_nm(nn+1,mm)=Bigvalue
				if(quant%wavenum > kc)then
					betanm=sqrt(quant%wavenum**2d0-kc**2d0)
					quant%ports(pp)%impedance_TM_nm(nn+1,mm)=impedence0*betanm/quant%wavenum
					if(ptree_A%MyID==Main_ID)write(*,*)pp,'CIR','TM',nn,mm,kc,quant%wavenum,quant%ports(pp)%impedance_TM_nm(nn+1,mm)
				endif
			enddo
			enddo
		else if(quant%ports(pp)%type==1)then
			do nn=0,quant%ports(pp)%nmax
				do mm=0,quant%ports(pp)%mmax
					quant%ports(pp)%A_TE_nm(nn+1,mm+1)=0
					quant%ports(pp)%impedance_TE_nm(nn+1,mm+1)=Bigvalue
					if(nn>0 .or. mm>0)then ! the lowest TE mode is 01 or 10
						kc = sqrt((nn*pi/quant%ports(pp)%a)**2d0+(mm*pi/quant%ports(pp)%b)**2d0)
						quant%ports(pp)%A_TE_nm(nn+1,mm+1)=1d0/sqrt(quant%ports(pp)%a/quant%ports(pp)%b*mm**2d0*A_nm_rec(nn+1,mm+1) + quant%ports(pp)%b/quant%ports(pp)%a*nn**2d0*B_nm_rec(nn+1,mm+1))
						if(quant%wavenum > kc)then
							betanm=sqrt(quant%wavenum**2d0-kc**2d0)
							quant%ports(pp)%impedance_TE_nm(nn+1,mm+1)=quant%wavenum*impedence0/betanm
							if(ptree_A%MyID==Main_ID)write(*,*)pp,'RECT','TE',nn,mm,kc,quant%wavenum,quant%ports(pp)%impedance_TE_nm(nn+1,mm+1)
						endif
					endif
				enddo
			enddo

			do nn=0,quant%ports(pp)%nmax
				do mm=0,quant%ports(pp)%mmax
					quant%ports(pp)%A_TM_nm(nn+1,mm+1)=0
					quant%ports(pp)%impedance_TM_nm(nn+1,mm+1)=Bigvalue
					if(nn>0 .and. mm>0)then ! the lowest TM mode is 11
						kc = sqrt((nn*pi/quant%ports(pp)%a)**2d0+(mm*pi/quant%ports(pp)%b)**2d0)
						quant%ports(pp)%A_TM_nm(nn+1,mm+1)=1d0/sqrt(quant%ports(pp)%a/quant%ports(pp)%b*mm**2d0*B_nm_rec(nn+1,mm+1) + quant%ports(pp)%b/quant%ports(pp)%a*nn**2d0*A_nm_rec(nn+1,mm+1))
						if(quant%wavenum > kc)then
							betanm=sqrt(quant%wavenum**2d0-kc**2d0)
							quant%ports(pp)%impedance_TM_nm(nn+1,mm+1)=impedence0*betanm/quant%wavenum
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
	t1 = OMP_get_wtime()
	allocate(xyz(3,quant%Nunk))
	do ii=1, quant%Nunk_int+quant%Nunk_port
		xyz(:,ii) = quant%xyz(:,quant%maxnode+ii)
	enddo
	do ii=quant%Nunk_int+quant%Nunk_port+1,quant%Nunk
		xyz(:,ii) = quant%xyz(:,quant%maxnode+ii-quant%Nunk_port)
	enddo

    allocate(Permutation(quant%Nunk))
	call BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat_A,option_A,stats_A,msh_A,ker_A,ptree_A,Coordinates=xyz)
	deallocate(Permutation) ! caller can use this permutation vector if needed
	deallocate(xyz)
	t2 = OMP_get_wtime()
	!**** computation of the construction phase
    call BPACK_construction_Element(bmat_A,option_A,stats_A,msh_A,ker_A,ptree_A)
	!**** print statistics

	if(.not. (quant%SI==1 .and. abs(quant%shift)<1e-14))then
		call PrintStat(stats_A,ptree_A)
	endif


	!***********************************************************************
	!**** construct compressed A - sigma I or  A - sigma real(A)
	if(quant%SI==1)then
		call CopyOptions(option_A,option_sh)
		!**** create the process tree, can use larger number of mpis if needed
		allocate(groupmembers(nmpi))
		do ii=1,nmpi
			groupmembers(ii)=(ii-1)
		enddo
		call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree_sh)
		deallocate(groupmembers)

		if(abs(quant%shift)<1e-14)then  ! if zero shift, no need to compress the shifted matrix
			call CopyMesh(msh_A,msh_sh)
			call CopyStat(stats_A,stats_sh)
			call BPACK_copy(bmat_A,bmat_sh)
		else
			!**** initialize stats_sh and option
			call InitStat(stats_sh)
			!**** register the user-defined function and type in ker
			ker_sh%QuantApp => quant
			ker_sh%FuncZmn => Zelem_EMSURF_Shifted
			!**** initialization of the construction phase
			t1 = OMP_get_wtime()
			allocate(xyz(3,quant%Nunk))
			do ii=1, quant%Nunk_int+quant%Nunk_port
				xyz(:,ii) = quant%xyz(:,quant%maxnode+ii)
			enddo
			do ii=quant%Nunk_int+quant%Nunk_port+1,quant%Nunk
				xyz(:,ii) = quant%xyz(:,quant%maxnode+ii-quant%Nunk_port)
			enddo
			allocate(Permutation(quant%Nunk))
			call BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat_sh,option_sh,stats_sh,msh_sh,ker_sh,ptree_sh,Coordinates=xyz)
			deallocate(Permutation) ! caller can use this permutation vector if needed
			deallocate(xyz)
			t2 = OMP_get_wtime()
			!**** computation of the construction phase
			call BPACK_construction_Element(bmat_sh,option_sh,stats_sh,msh_sh,ker_sh,ptree_sh)
		endif


		!**** factorization phase
		call BPACK_Factorization(bmat_sh,option_sh,stats_sh,ptree_sh,msh_sh)
		! !**** solve phase
		! call EM_solve_SURF(bmat_sh,option_sh,msh_sh,quant,ptree_sh,stats_sh)
		!**** print statistics
		call PrintStat(stats_sh,ptree_sh)
	endif


	!***********************************************************************
	!**** construct compressed real(A)
	if(quant%CMmode==1)then ! solve the characteristic mode
		call CopyOptions(option_A,option_B)
		!**** register the user-defined function and type in ker
		ker_B%QuantApp => quant
		ker_B%FuncZmn => Zelem_EMSURF_Real
		!**** initialize stats_sh and option
		call InitStat(stats_B)
		!**** create the process tree, can use larger number of mpis if needed
		allocate(groupmembers(nmpi))
		do ii=1,nmpi
			groupmembers(ii)=(ii-1)
		enddo
		call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree_B)
		deallocate(groupmembers)
		!**** initialization of the construction phase
		t1 = OMP_get_wtime()
		allocate(xyz(3,quant%Nunk))
		do ii=1, quant%Nunk_int+quant%Nunk_port
			xyz(:,ii) = quant%xyz(:,quant%maxnode+ii)
		enddo
		do ii=quant%Nunk_int+quant%Nunk_port+1,quant%Nunk
			xyz(:,ii) = quant%xyz(:,quant%maxnode+ii-quant%Nunk_port)
		enddo
		allocate(Permutation(quant%Nunk))
		call BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat_B,option_B,stats_B,msh_B,ker_B,ptree_B,Coordinates=xyz)
		deallocate(Permutation) ! caller can use this permutation vector if needed
		deallocate(xyz)
		t2 = OMP_get_wtime()
		!**** computation of the construction phase
		call BPACK_construction_Element(bmat_B,option_B,stats_B,msh_B,ker_B,ptree_B)
		!**** factorization phase
		if(quant%SI==0)then
			call BPACK_Factorization(bmat_B,option_B,stats_B,ptree_B,msh_B)
		endif
		!**** print statistics
		call PrintStat(stats_B,ptree_B)
	endif

	!***********************************************************************
	!**** eigen solver
	allocate(eigval(quant%nev))
	allocate(eigvec(Nunk_loc,quant%nev))
	call BPACK_Eigen(bmat_A,option_A,ptree_A, stats_A,bmat_B,option_B,ptree_B, stats_B, bmat_sh,option_sh,ptree_sh,stats_sh, quant%Nunk,Nunk_loc, quant%nev,quant%tol_eig,quant%CMmode,quant%SI,quant%shift,quant%which,nconv,eigval,eigvec)


	! if(ptree_A%MyID==Main_ID)then
	! 	write(*,*)'Checking 1-norm of the eigenvectors: '
	! endif
	maxnorm=0
	allocate(norms(quant%nev))
	do nn=1,quant%nev
		norm1 = fnorm(eigvec(:,nn:nn), Nunk_loc, 1, '1')
		normi = fnorm(eigvec(:,nn:nn), Nunk_loc, 1, 'I')
		call MPI_ALLREDUCE(MPI_IN_PLACE, norm1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree_A%Comm, ierr)
		call MPI_ALLREDUCE(MPI_IN_PLACE, normi, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree_A%Comm, ierr)
		! if(ptree_A%MyID==Main_ID)then
		! 	write(*,*)norm1,normi,norm1/normi
		! endif
		norm1 =norm1/normi
		if(norm1>maxnorm)then
			nn1=nn
			maxnorm = norm1
		endif
		norms(nn)=norm1
	enddo
	allocate(order(quant%nev))
	allocate(norms_tmp(quant%nev))
	norms_tmp = norms
	call quick_sort(norms_tmp,order,quant%nev)
	deallocate(norms_tmp)
	do nn=quant%nev,1,-1
		if(ptree_A%MyID==Main_ID .and. norms(order(nn))>norm_thresh)then
			write(*,*)'eigval',dble(eigval(order(nn))),aimag(eigval(order(nn))),norms(order(nn)),order(nn)
		endif
	enddo

	!**** dump the frequency to file, which records the frequencies that have been tried
	if(ptree_A%MyID==Main_ID)then
		write(substring1 , *) quant%freq
		inquire(file=trim(adjustl(substring2))//"_freq_history.txt", exist=exist)
		if (exist) then
			open(200, status='SCRATCH', action="READWRITE")
			open(12, file=trim(adjustl(substring2))//"_freq_history.txt", status="old", action="readwrite")
			read(12,*)Nfreq
			f_exist=.false.
			do nn=1,Nfreq
				read(12,*)strings
				if(trim(adjustl(substring1))==trim(adjustl(strings)))f_exist=.true.
				write(200,*)trim(adjustl(strings))
			enddo
			if(.not. f_exist)then
				write(200,*)trim(adjustl(substring1))
				rewind(12)
				rewind(200)
				Nfreq = Nfreq+1
				write(12,*)Nfreq
				do nn=1,Nfreq
					read(200,*)strings
					write(12,*)trim(adjustl(strings))
				enddo
			endif
			close(12)
			close(200)
		else
			open(12, file=trim(adjustl(substring2))//"_freq_history.txt", status="new", action="write")
			write(12,*) 1
			write(12,*) trim(adjustl(substring1))
			close(12)
		end if
	endif

	allocate(eigvec_glo(quant%Nunk))
	allocate(eigvec_glo_ref(quant%Nunk))

	order=0  ! order maps the ith eigenvector to the nth mode found so far
	do ii=1,quant%nev
		if(norms(ii)>norm_thresh)then
			eigvec_glo=0
			eigvec_glo(msh_A%idxs:msh_A%idxe)=eigvec(:,ii)
			call MPI_ALLREDUCE(MPI_IN_PLACE,eigvec_glo,quant%Nunk,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree_A%Comm,ierr)
			eigvec_glo = eigvec_glo/dot_product(eigvec_glo,eigvec_glo)

			!**** compare with eigenvectors in file
			if(ptree_A%MyID==Main_ID)then
				write(substring1 , *) quant%freq
				inquire(file=trim(adjustl(substring2))//"_Nmodes.txt", exist=exist)
				if (exist) then
					open(12, file=trim(adjustl(substring2))//"_Nmodes.txt", status="old", action="readwrite")
					read(12,*)Nmode
					rewind(12)
					f_exist=.false.
					do nn=1,Nmode
						write(strings, *) nn
						open(13, file=trim(adjustl(substring2))//"_mode_vec"//trim(adjustl(strings))//".out", status="old", action="read")
						do i=1,quant%Nunk
							read(13,*)rtemp1,rtemp2
							eigvec_glo_ref(i) = rtemp1 + junit*rtemp2
						enddo
						close(13)
						write(*,*)abs(dot_product(eigvec_glo_ref,eigvec_glo)),'comparing eigenvector',ii, ' to existing mode',nn
						if(abs(dot_product(eigvec_glo_ref,eigvec_glo))>dotproduct_thresh)then
							f_exist=.true.
							order(ii)=nn
							exit
						endif
					enddo
					if(.not. f_exist)then
						Nmode = Nmode+1
						order(ii)=Nmode
						write(*,*)'New mode detected!'
						write(12,*)Nmode

						write(strings, *) Nmode
						open(13, file=trim(adjustl(substring2))//"_mode_vec"//trim(adjustl(strings))//".out", status="new", action="write")
						do nn=1,quant%Nunk
							write(13,*) dble(eigvec_glo(nn)),aimag(eigvec_glo(nn))
						enddo
						close(13)
					else
						write(strings, *) nth_mode 
						open(15, file=trim(adjustl(substring2))//'_EigVals_'//trim(adjustl(strings))//'.out', status="old", action="read")
						valmin=1e30
						do
							read(unit=15,fmt=*,iostat=io) freq0,val0
							if (io/=0) exit
							valmin = min(valmin,val0)
						enddo 
						close(15)
						if(valmin>abs(eigval(ii)))then
							open(13, file=trim(adjustl(substring2))//"_mode_vec"//trim(adjustl(strings))//".out", status="old", action="write")
							do nn=1,quant%Nunk
								write(13,*) dble(eigvec_glo(nn)),aimag(eigvec_glo(nn))
							enddo
							close(13)						
						endif					
					endif
					close(12)
				else
					Nmode=1
					order(ii)=Nmode
					open(12, file=trim(adjustl(substring2))//"_Nmodes.txt", status="new", action="write")
					write(12,*)Nmode
					close(12)

					write(strings, *) Nmode
					open(13, file=trim(adjustl(substring2))//"_mode_vec"//trim(adjustl(strings))//".out", status="new", action="write")
					do nn=1,quant%Nunk
						write(13,*) dble(eigvec_glo(nn)),aimag(eigvec_glo(nn))
					enddo
					close(13)
				end if

				!**** dump the eigenvalues of the found mode at the current frequency to file, this serves as the historical data for GPTune (note that records in the file can have duplicated frequencies)
				if(order(ii)>0)then
					write(substring , *) order(ii)
					write(substring1 , *) quant%freq
					inquire(file=trim(adjustl(substring2))//'_EigVals_'//trim(adjustl(substring))//'.out', exist=exist)
					if(exist)then
						open(15, file=trim(adjustl(substring2))//'_EigVals_'//trim(adjustl(substring))//'.out', status="old", position="append", action="write")
					else
						open(15, file=trim(adjustl(substring2))//'_EigVals_'//trim(adjustl(substring))//'.out', status="new", action="write")
					endif
					write(15,*)trim(adjustl(substring1)),abs(eigval(ii))
					close(15)
				endif

			endif
		endif
	enddo
	deallocate(eigvec_glo)
	deallocate(eigvec_glo_ref)
	call MPI_Bcast(order, quant%nev, MPI_integer, Main_ID, ptree_A%Comm, ierr)


	if(quant%postprocess==1)then

		if(ptree_A%MyID==Main_ID)then
			write(*,*)'Postprocessing: '
		endif

		do ii=1,quant%nev
			if(order(ii)>0)then
				write(substring , *) order(ii)
				write(substring1 , *) quant%freq
				call current_node_patch_mapping(0,trim(adjustl(substring2))//'_EigVecJ_'//trim(adjustl(substring))//'_freq_'//trim(adjustl(substring1))//'.out',eigvec(:,ii),msh_A,quant,ptree_A)
				call current_node_patch_mapping(1,trim(adjustl(substring2))//'_EigVecM_'//trim(adjustl(substring))//'_freq_'//trim(adjustl(substring1))//'.out',eigvec(:,ii),msh_A,quant,ptree_A)

				if(quant%CMmode==1)then
					if(ptree_A%MyID==Main_ID)open (100, file=trim(adjustl(substring2))//'VV_bistatic_'//'EigVec_'//trim(adjustl(substring))//'_freq_'//trim(adjustl(substring1))//'.txt')
					dtheta=180./quant%RCS_Nsample
					phi=0
					do i=0, quant%RCS_Nsample
						theta=i*dtheta
						ctemp_loc=(0.,0.)
						rcs=0
						!$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
						do edge=msh_A%idxs,msh_A%idxe
							call VV_polar_SURF(theta,phi,msh_A%new2old(edge),ctemp_1,eigvec(edge-msh_A%idxs+1,ii),quant)
							ctemp_loc=ctemp_loc+ctemp_1
						enddo
						!$omp end parallel do

						call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree_A%Comm,ierr)

						rcs=(abs(quant%wavenum*ctemp))**2/4/pi
						!rcs=rcs/quant%wavelength
						rcs=10*log10(rcs)
						if(ptree_A%MyID==Main_ID)write(100,*)theta,rcs
					enddo
					if(ptree_A%MyID==Main_ID)close(100)
				endif
			endif
		enddo

		call CopyOptions(option_A,option_post)
		!**** register the user-defined function and type in ker
		ker_post%QuantApp => quant
		ker_post%FuncZmn => Zelem_EMSURF_Post
		!**** initialize stats_post and option
		call InitStat(stats_post)
		!**** create the process tree, can use larger number of mpis if needed
		allocate(groupmembers(nmpi))
		do ii=1,nmpi
			groupmembers(ii)=(ii-1)
		enddo
		call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree_post)
		deallocate(groupmembers)
		!**** initialization of the construction phase
		t1 = OMP_get_wtime()
		allocate(xyz(3,quant%Nunk))
		do ii=1, quant%Nunk_int+quant%Nunk_port
			xyz(:,ii) = quant%xyz(:,quant%maxnode+ii)
		enddo
		do ii=quant%Nunk_int+quant%Nunk_port+1,quant%Nunk
			xyz(:,ii) = quant%xyz(:,quant%maxnode+ii-quant%Nunk_port)
		enddo
		allocate(Permutation(quant%Nunk))
		call BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat_post,option_post,stats_post,msh_post,ker_post,ptree_post,Coordinates=xyz)
		deallocate(Permutation) ! caller can use this permutation vector if needed
		deallocate(xyz)
		t2 = OMP_get_wtime()
		!**** computation of the construction phase
		call BPACK_construction_Element(bmat_post,option_post,stats_post,msh_post,ker_post,ptree_post)
		allocate(E_normal(Nunk_loc,quant%nev))
		E_normal=0
		call BPACK_Mult('N', Nunk_loc, quant%nev, eigvec, E_normal, bmat_post, ptree_post, option_post, stats_post)
		do nn=1,quant%nev
			if(order(nn)>0)then
				call EM_cavity_postprocess(option_A,msh_A,quant,ptree_A,stats_A,eigvec(:,nn),order(nn),norms(nn),eigval(nn),E_normal(:,nn),nn,substring2)
			endif
		enddo
		deallocate(E_normal)

		call delete_proctree(ptree_post)
		call delete_Hstat(stats_post)
		call delete_mesh(msh_post)
		call delete_kernelquant(ker_post)
		call BPACK_delete(bmat_post)
	endif


	retval(1) = abs(eigval(nn1))
	retval(2) = maxnorm

	if(ptree_A%MyID==Main_ID .and. option_A%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	deallocate(eigval)
	deallocate(eigvec)
	deallocate(order)
	deallocate(norms)

	if(quant%SI==1)then
		call delete_proctree(ptree_sh)
		call delete_Hstat(stats_sh)
		call delete_mesh(msh_sh)
		call delete_kernelquant(ker_sh)
		call BPACK_delete(bmat_sh)
	endif

	call delete_proctree(ptree_A)
	call delete_Hstat(stats_A)
	call delete_mesh(msh_A)
	call delete_kernelquant(ker_A)
	call BPACK_delete(bmat_A)

	if(quant%CMmode==1)then ! solve the characteristic mode
		call delete_proctree(ptree_B)
		call delete_Hstat(stats_B)
		call delete_mesh(msh_B)
		call delete_kernelquant(ker_B)
		call BPACK_delete(bmat_B)
	endif

	!**** deletion of quantities
	call delete_quant_EMSURF(quant)



	if(parent/= MPI_COMM_NULL)then
		call MPI_REDUCE(retval, MPI_BOTTOM, 2, MPI_double_precision, MPI_MAX, Main_ID, parent,ierr)
		call MPI_Comm_disconnect(parent,ierr)
	endif

	call blacs_exit(1)
	call MPI_Finalize(ierr)

    ! ! ! ! pause

end PROGRAM ButterflyPACK_IE_3D






