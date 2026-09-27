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
!> @brief This is an example that solves 3D sound-soft acoustic scattering with the Chebyshev-based rectangular-polar Nystrom method (Bruno & Garza, JCP 2020) and ButterflyPACK's direct solver.
!> @details Analytic geometries: sphere (cubed-sphere patches), cube, disk and square plate. Closed surfaces use the combined-field equation, open surfaces the single-layer equation. \n
!> For the sphere, the driver checks the forward map against the exact eigenvalues of spherical harmonics and the far field against the Mie series. \n
!> Note that instead of the use of precision dependent subroutine/module/type names "z_", one can also use the following \n
!> #define DAT 0 \n
!> #include "zButterflyPACK_config.fi" \n
!> which will macro replace precision-independent subroutine/module/type names "X" with "z_X" defined in SRC_DOUBLECOMLEX with double-complex precision

! This exmple works with double-complex precision data
PROGRAM ButterflyPACK_Acoustic_SURF_CHEB
	use z_BPACK_DEFS
	use Acoustic_SURF_MODULE_CHEB

	use z_BPACK_utilities
	use z_BPACK_structure
	use z_BPACK_factor
	use z_BPACK_constr
	use z_BPACK_Solve_Mul
#ifdef HAVE_OPENMP
	use omp_lib
#endif
	use z_MISC_Utilities
	implicit none

	integer ii
	real(kind=8) t1, t2
	character(len=1024)  :: strings, strings1
	integer :: ierr
	type(z_Hoption)::option
	type(z_Hstat)::stats
	type(z_mesh)::msh
	type(z_Bmatrix)::bmat
	type(z_kernelquant)::ker
	type(quant_ACOUSTIC_CHEB), target::quant
	type(z_proctree)::ptree
	integer, allocatable:: groupmembers(:)
	integer nmpi, provided
	integer, allocatable::Permutation(:)
	integer Nunk_loc
	integer nargs, flag
	integer v_major, v_minor, v_bugfix

	call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)
	call MPI_Comm_size(MPI_Comm_World, nmpi, ierr)
	allocate (groupmembers(nmpi))
	do ii = 1, nmpi
		groupmembers(ii) = (ii - 1)
	enddo
	call z_CreatePtree(nmpi, groupmembers, MPI_Comm_World, ptree)
	deallocate (groupmembers)

	if (ptree%MyID == Main_ID) then
		write (*, *) "-------------------------------Program Start----------------------------------"
		write (*, *) "ButterflyPACK_Acoustic_SURF_CHEB"
		call z_BPACK_GetVersionNumber(v_major, v_minor, v_bugfix)
		write (*, '(A23,I1,A1,I1,A1,I1,A1)') " ButterflyPACK Version:", v_major, ".", v_minor, ".", v_bugfix
		write (*, *) "   "
	endif

	!**** initialize stats and option
	call z_InitStat(stats)
	call z_SetDefaultOptions(option)

	option%format = HODLR !HMAT!
	option%LRlevel = 0 ! low-rank blocks; -option --lrlevel 100 selects butterfly blocks
	option%tol_comp = 1d-6
	option%tol_rand = option%tol_comp
	option%tol_Rdetect = option%tol_comp*1d-1
	option%near_para = 2.01d0
	option%sample_para = 4d0
	option%verbosity = 1

	!**** read user-defined quantity parameters and ButterflyPACK options
	nargs = iargc()
	ii = 1
	do while (ii <= nargs)
		call getarg(ii, strings)
		if (trim(strings) == '-quant') then ! user-defined quantity parameters
			flag = 1
			do while (flag == 1)
				ii = ii + 1
				if (ii <= nargs) then
					call getarg(ii, strings)
					if (strings(1:2) == '--') then
						ii = ii + 1
						call getarg(ii, strings1)
						if (trim(strings) == '--geo') then
							select case (trim(strings1))
							case ('sphere')
								quant%geo = CHEB_SPHERE
							case ('cube')
								quant%geo = CHEB_CUBE
							case ('disk')
								quant%geo = CHEB_DISK
							case ('plate')
								quant%geo = CHEB_PLATE
							case default
								if (ptree%MyID == Main_ID) write (*, *) 'unknown geometry: ', trim(strings1)
								stop
							end select
						else if (trim(strings) == '--radius') then
							read (strings1, *) quant%radius
						else if (trim(strings) == '--split') then
							read (strings1, *) quant%split
						else if (trim(strings) == '--N') then
							read (strings1, *) quant%N
						else if (trim(strings) == '--nbeta') then
							read (strings1, *) quant%Nbeta
						else if (trim(strings) == '--delta') then
							read (strings1, *) quant%delta
						else if (trim(strings) == '--pxi') then
							read (strings1, *) quant%pxi
						else if (trim(strings) == '--pedge') then
							read (strings1, *) quant%pedge
						else if (trim(strings) == '--wavelength') then
							read (strings1, *) quant%wavelength
							quant%wavenum = 2d0*BPACK_pi/quant%wavelength
						else if (trim(strings) == '--wavenum') then
							read (strings1, *) quant%wavenum
							quant%wavelength = 2d0*BPACK_pi/quant%wavenum
						else if (trim(strings) == '--unknown') then
							if (trim(strings1) == 'phi') quant%unknown = UNK_PHI
							if (trim(strings1) == 'psi') quant%unknown = UNK_PSI
							if (trim(strings1) == 'auto') quant%unknown = -1
						else if (trim(strings) == '--inc_theta') then
							read (strings1, *) quant%inc_theta
						else if (trim(strings) == '--inc_phi') then
							read (strings1, *) quant%inc_phi
						else if (trim(strings) == '--eigtest') then
							read (strings1, *) quant%eigtest
						else if (trim(strings) == '--nfar') then
							read (strings1, *) quant%nfar
						else if (trim(strings) == '--naive_chord') then
							read (strings1, *) quant%naive_chord
						else if (trim(strings) == '--gmres_restart') then
							read (strings1, *) quant%gmres_restart
						else
							if (ptree%MyID == Main_ID) write (*, *) 'ignoring unknown quant: ', trim(strings)
						endif
					else
						flag = 0
					endif
				else
					flag = 0
				endif
			enddo
		else if (trim(strings) == '-option') then ! options of ButterflyPACK
			call z_ReadOption(option, ptree, ii)
		else
			if (ptree%MyID == Main_ID) write (*, *) 'ignoring unknown argument: ', trim(strings)
			ii = ii + 1
		endif
	enddo

	quant%dinc = (/sin(quant%inc_theta*BPACK_pi/180d0)*cos(quant%inc_phi*BPACK_pi/180d0), &
		&	sin(quant%inc_theta*BPACK_pi/180d0)*sin(quant%inc_phi*BPACK_pi/180d0), cos(quant%inc_theta*BPACK_pi/180d0)/)

	if (ptree%MyID == Main_ID) then
		write (*, *) ''
		write (*, *) 'Acoustic sound-soft scattering, Chebyshev rectangular-polar Nystrom'
		write (*, *) 'wavenumber:', quant%wavenum, ' wavelength:', quant%wavelength
		write (*, *) ''
	endif

	!**** geometry, Fejer nodes and near-field weights
	t1 = MPI_Wtime()
	call geo_modeling_CHEB(quant, ptree)
	call build_nodes_CHEB(quant, ptree)
	call precompute_nearfield_CHEB(quant, ptree)
	t2 = MPI_Wtime()
	if (ptree%MyID == Main_ID) write (*, *) 'geometry and near-field precomputation:', t2 - t1, 'Seconds'

	!**** register the user-defined function and type in ker
	ker%QuantApp => quant
	ker%FuncZmn => Zelem_Acoustic_CHEB

	!**** initialization of the construction phase
	allocate (Permutation(quant%Nunk))
	call z_PrintOptions(option, ptree)
	call z_BPACK_construction_Init(quant%Nunk, Permutation, Nunk_loc, bmat, option, stats, msh, ker, ptree, Coordinates=quant%xyz)
	deallocate (Permutation) ! caller can use this permutation vector if needed

	!**** computation of the construction phase
	call z_BPACK_construction_Element(bmat, option, stats, msh, ker, ptree)

	!**** forward-map test (sphere only), before the factorization overwrites the forward operator
	call forward_map_test_CHEB(bmat, option, msh, quant, ptree, stats)

	!**** factorization phase
	call z_BPACK_Factorization(bmat, option, stats, ptree, msh)

	!**** solve phase
	call solve_scattering_CHEB(bmat, option, msh, quant, ptree, stats)

	!**** print statistics
	call z_PrintStat(stats, ptree)

	!**** deletion of quantities
	call delete_quant_Acoustic_CHEB(quant)
	call z_delete_proctree(ptree)
	call z_delete_Hstat(stats)
	call z_delete_mesh(msh)
	call z_delete_kernelquant(ker)
	call z_BPACK_delete(bmat)

	if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) "-------------------------------program end-------------------------------------"

	call z_blacs_exit_wrp(1)
	call MPI_Finalize(ierr)

end PROGRAM ButterflyPACK_Acoustic_SURF_CHEB
