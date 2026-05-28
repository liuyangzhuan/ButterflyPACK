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
!> @brief This is an example that solves a 3D EFIE/MFIE/CFIE system for electromagnetics scattering.
!> @details This mixed-precision example constructs a double-complex forward matrix A,
!> constructs/factorizes a single-complex preconditioner M, and solves the left-preconditioned
!> system with double-complex TFQMR.


! This example works with double-complex outer iterations and a single-complex preconditioner.
PROGRAM ButterflyPACK_IE_3D
    use c_BPACK_DEFS
	use EMSURF_MODULE_MP, only: quant_EMSURF, delete_quant_EMSURF, gauss_points, &
		geo_modeling_SURF, Zelem_EMSURF, Zelem_EMSURF_block, EM_solve_SURF
	use EMSURF_MODULE, only: quant_EMSURF_DP => quant_EMSURF, &
		delete_quant_EMSURF_DP => delete_quant_EMSURF, gauss_points_DP => gauss_points, &
		geo_modeling_SURF_DP => geo_modeling_SURF, Zelem_EMSURF_DP => Zelem_EMSURF, &
		Zelem_EMSURF_block_DP => Zelem_EMSURF_block

	use c_BPACK_structure
	use c_BPACK_factor
	use c_BPACK_constr
#ifdef HAVE_OPENMP
	use omp_lib
#endif
	use c_MISC_Utilities
	use c_BPACK_utilities

	use z_BPACK_DEFS, only: z_Hoption, z_Hstat, z_mesh, z_Bmatrix, z_kernelquant, z_proctree
	use z_BPACK_structure, only: z_BPACK_delete, z_delete_mesh
	use z_BPACK_constr, only: z_BPACK_construction_Init, z_BPACK_construction_Element
	use z_MISC_Utilities, only: z_CreatePtree
	use z_BPACK_utilities, only: z_InitStat, z_SetDefaultOptions, z_ReadOption, z_PrintOptions, &
		z_PrintStat, z_delete_proctree, z_delete_Hstat, z_delete_kernelquant
    implicit none

    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num, provided
	integer seed_myid(50)
	integer times(8)
	real(kind=8) t1,t2

	character(len=:),allocatable  :: string
	character(len=1024)  :: strings,strings1
	character(len=6)  :: info_env
	integer :: length,edge
	integer :: ierr
	integer*8 oldmode,newmode
	type(c_Hoption)::option_m
	type(c_Hstat)::stats_m
	type(c_mesh)::msh_m
	type(c_Bmatrix)::bmat_m
	type(c_kernelquant)::ker_m
	type(quant_EMSURF),target::quant_m
	type(c_proctree)::ptree_m
	type(z_Hoption)::option_a
	type(z_Hstat)::stats_a
	type(z_mesh)::msh_a
	type(z_Bmatrix)::bmat_a
	type(z_kernelquant)::ker_a
	type(quant_EMSURF_DP),target::quant_a
	type(z_proctree)::ptree_a
	integer,allocatable:: groupmembers(:)
	integer nmpi
	real(kind=8),allocatable::xyz_a(:,:),xyz_m(:,:)
	integer,allocatable::Permutation_a(:),Permutation_m(:)
	integer Nunk_loc_a,Nunk_loc_m
	integer nargs,flag,opt_start,ii_z
	integer v_major,v_minor,v_bugfix

	! nmpi and groupmembers should be provided by the user

	! call MPI_Init(ierr)
	call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)

	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo

	call c_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree_m)
	call z_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree_a)
	deallocate(groupmembers)


	if(ptree_m%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "ButterflyPACK_IE_3D_MP"
	call c_BPACK_GetVersionNumber(v_major,v_minor,v_bugfix)
	write(*,'(A23,I1,A1,I1,A1,I1,A1)') " ButterflyPACK Version:",v_major,".",v_minor,".",v_bugfix
    write(*,*) "   "
	endif

	!**** initialize stats and options
	call c_InitStat(stats_m)
	call z_InitStat(stats_a)
	call c_SetDefaultOptions(option_m)
	call z_SetDefaultOptions(option_a)


	!**** initialize the user-defined derived types
	! compute the quadrature rules
    quant_m%integral_points=6
    allocate (quant_m%ng1(quant_m%integral_points), quant_m%ng2(quant_m%integral_points), &
		quant_m%ng3(quant_m%integral_points), quant_m%gauss_w(quant_m%integral_points))
    call gauss_points(quant_m)

    quant_a%integral_points=6
    allocate (quant_a%ng1(quant_a%integral_points), quant_a%ng2(quant_a%integral_points), &
		quant_a%ng3(quant_a%integral_points), quant_a%gauss_w(quant_a%integral_points))
    call gauss_points_DP(quant_a)

    !*************************input******************************
	quant_m%DATA_DIR='../EXAMPLE/EM3D_DATA/sphere_2300'
	quant_a%DATA_DIR='../EXAMPLE/EM3D_DATA/sphere_2300'

	quant_m%mesh_normal=1
	quant_a%mesh_normal=1
	quant_m%scaling=1d0
	quant_a%scaling=1d0
	quant_m%wavelength=2.0
	quant_a%wavelength=2.0d0
	quant_m%freq=1/quant_m%wavelength/sqrt(BPACK_mu0*BPACK_eps0)
	quant_a%freq=1/quant_a%wavelength/sqrt(BPACK_mu0*BPACK_eps0)
	quant_m%RCS_static=2
	quant_a%RCS_static=2
    quant_m%RCS_Nsample=1000
    quant_a%RCS_Nsample=1000
	quant_m%CFIE_alpha=1.0
	quant_a%CFIE_alpha=1.0d0

	option_m%ErrSol=1
	option_a%ErrSol=1
	option_m%format=  HODLR !HMAT!
	option_a%format=  HODLR !HMAT!
	option_m%near_para=2.01d0
	option_a%near_para=2.01d0
	option_m%verbosity=1
	option_a%verbosity=1
	option_m%ILU=0
	option_a%ILU=0
	option_m%forwardN15flag=0
	option_a%forwardN15flag=0
	option_m%LRlevel=100
	option_a%LRlevel=100
	option_m%tol_itersol=1d-5
	option_a%tol_itersol=1d-5
	option_m%sample_para=4d0
	option_a%sample_para=4d0
	option_m%knn=50
	option_a%knn=50

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
							quant_m%data_dir=trim(strings1)
							quant_a%data_dir=trim(strings1)
						else if	(trim(strings)=='--wavelength')then
							read(strings1,*)quant_a%wavelength
							quant_m%wavelength=quant_a%wavelength
							quant_a%freq=1/quant_a%wavelength/sqrt(BPACK_mu0*BPACK_eps0)
							quant_m%freq=quant_a%freq
						else if (trim(strings)=='--freq')then
							read(strings1,*)quant_a%freq
							quant_m%freq=quant_a%freq
							quant_a%wavelength=1/quant_a%freq/sqrt(BPACK_mu0*BPACK_eps0)
							quant_m%wavelength=quant_a%wavelength
						else
							if(ptree_m%MyID==Main_ID)write(*,*)'ignoring unknown quant: ', trim(strings)
						endif
					else
						flag=0
					endif
				else
					flag=0
				endif
			enddo
		else if(trim(strings)=='-option')then ! options of ButterflyPACK
			opt_start=ii
			call c_ReadOption(option_m,ptree_m,ii)
			ii_z=opt_start
			call z_ReadOption(option_a,ptree_a,ii_z)
		else
			if(ptree_m%MyID==Main_ID)write(*,*)'ignoring unknown argument: ',trim(strings)
			ii=ii+1
		endif
	enddo

	if(option_m%precon==NOPRECON)then
		if(ptree_m%MyID==Main_ID)write(*,*)'mixed-precision solve requires c_BPACK_Inv_Mult; switching M precon ', &
			'from NOPRECON to DIRECT'
		option_m%precon=DIRECT
	endif

    quant_m%wavenum=2*BPACK_pi/quant_m%wavelength
    quant_a%wavenum=2*BPACK_pi/quant_a%wavelength


   !***********************************************************************
	if(ptree_m%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'frequency:',quant_a%freq
   write (*,*) 'wavelength:',quant_a%wavelength
   write (*,*) ''
	endif
   !***********************************************************************


	!**** geometry generation and discretization
	call geo_modeling_SURF_DP(quant_a,ptree_a%Comm,quant_a%DATA_DIR)
	call geo_modeling_SURF(quant_m,ptree_m%Comm,quant_m%DATA_DIR)

	if(quant_a%Nunk/=quant_m%Nunk)then
		if(ptree_m%MyID==Main_ID)write(*,*)'double A and single M geometries have different Nunk:', quant_a%Nunk, quant_m%Nunk
		stop
	endif

	option_a%touch_para = 3* quant_a%minedgelength
	option_m%touch_para = 3* quant_m%minedgelength

	!**** register the user-defined functions and types
	ker_a%QuantApp => quant_a
	ker_a%FuncZmn => Zelem_EMSURF_DP
	ker_a%FuncZmnBlock => Zelem_EMSURF_block_DP

	ker_m%QuantApp => quant_m
	ker_m%FuncZmn => Zelem_EMSURF
	ker_m%FuncZmnBlock => Zelem_EMSURF_block

	!**** initialization of the construction phases
	allocate(xyz_a(3,quant_a%Nunk),xyz_m(3,quant_m%Nunk))
	do ii=1, quant_a%Nunk
		xyz_a(:,ii) = quant_a%xyz(:,quant_a%maxnode+ii)
		xyz_m(:,ii) = xyz_a(:,ii)
	enddo

	t1 = MPI_Wtime()
    allocate(Permutation_a(quant_a%Nunk))
	if(ptree_m%MyID==Main_ID .and. option_m%verbosity>=0)write(*,*) 'Constructing double-complex A'
	call z_PrintOptions(option_a,ptree_a)
	call z_BPACK_construction_Init(quant_a%Nunk,Permutation_a,Nunk_loc_a,bmat_a,option_a,stats_a,msh_a,ker_a,ptree_a,Coordinates=xyz_a)
	deallocate(Permutation_a)

    allocate(Permutation_m(quant_m%Nunk))
	if(ptree_m%MyID==Main_ID .and. option_m%verbosity>=0)write(*,*) 'Constructing/factorizing single-complex M'
	call c_PrintOptions(option_m,ptree_m)
	call c_BPACK_construction_Init(quant_m%Nunk,Permutation_m,Nunk_loc_m,bmat_m,option_m,stats_m,msh_m,ker_m,ptree_m,Coordinates=xyz_m)
	deallocate(Permutation_m)
	deallocate(xyz_a,xyz_m)
	t2 = MPI_Wtime()

	if(Nunk_loc_a/=Nunk_loc_m .or. msh_a%idxs/=msh_m%idxs .or. msh_a%idxe/=msh_m%idxe)then
		if(ptree_m%MyID==Main_ID)write(*,*)'double A and single M have different local ownership'
		stop
	endif
	do ii=msh_m%idxs,msh_m%idxe
		if(msh_a%new2old(ii)/=msh_m%new2old(ii))then
			if(ptree_m%MyID==Main_ID)write(*,*)'double A and single M have different permutations at local index',ii
			stop
		endif
	enddo


	!**** computation of the construction phases
    call z_BPACK_construction_Element(bmat_a,option_a,stats_a,msh_a,ker_a,ptree_a)
    call c_BPACK_construction_Element(bmat_m,option_m,stats_m,msh_m,ker_m,ptree_m)


	!**** factorization phase for the single-complex preconditioner
	call c_BPACK_Factorization(bmat_m,option_m,stats_m,ptree_m,msh_m)


	!**** solve phase
	call EM_solve_SURF(bmat_a,option_a,msh_a,stats_a,ptree_a,bmat_m,option_m,msh_m,quant_a,ptree_m,stats_m)


	!**** print statistics
	if(ptree_m%MyID==Main_ID .and. option_m%verbosity>=0)write(*,*) 'Double-complex A statistics'
	call z_PrintStat(stats_a,ptree_a)
	if(ptree_m%MyID==Main_ID .and. option_m%verbosity>=0)write(*,*) 'Single-complex M statistics'
	call c_PrintStat(stats_m,ptree_m)


	!**** deletion of quantities
	if(ptree_m%MyID==Main_ID .and. option_m%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	call delete_quant_EMSURF_DP(quant_a)
	call delete_quant_EMSURF(quant_m)
	call z_delete_proctree(ptree_a)
	call c_delete_proctree(ptree_m)
	call z_delete_Hstat(stats_a)
	call c_delete_Hstat(stats_m)
	call z_delete_mesh(msh_a)
	call c_delete_mesh(msh_m)
	call z_delete_kernelquant(ker_a)
	call c_delete_kernelquant(ker_m)
	call z_BPACK_delete(bmat_a)
	call c_BPACK_delete(bmat_m)

	call c_blacs_exit_wrp(1)
	call MPI_Finalize(ierr)

end PROGRAM ButterflyPACK_IE_3D
