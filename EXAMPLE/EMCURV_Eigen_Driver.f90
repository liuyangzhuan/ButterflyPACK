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
!> @brief This is an example that uses PARPACK and ButterflyPACk to find the smallest eigen value or do characteristic mode analysis for a 2D EFIE system in electromagnetics.
!> @details Note that the use of the following \n
!> #define DAT 0 \n
!> #include "zButterflyPACK_config.fi" \n
!> which will macro replace subroutine, function, type names with those defined in SRC_DOUBLECOMLEX with double-complex precision

! This exmple works with double-complex precision data
PROGRAM ButterflyPACK_IE_2D
    use z_BPACK_DEFS
    use EMCURV_MODULE

	use z_BPACK_structure
	use z_BPACK_Solve_Mul
	use z_BPACK_factor
	use z_BPACK_constr
	use z_BPACK_Utilities
#ifdef HAVE_OPENMP
	use omp_lib
#endif
	use z_MISC_Utilities
    implicit none

	! include "mkl_vml.fi"

    real(kind=8) para
    real(kind=8) tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seed_myid(50)
	integer times(8)
	integer edge
	real(kind=8) t1,t2,t3, x,y,z,r,theta,phi
	complex(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:),eigs(:),eivect(:,:)

	character(len=:),allocatable  :: string
	character(len=1024)  :: strings,strings1
	character(len=6)  :: info_env
	integer :: length
	integer :: ierr
	integer*8 oldmode,newmode
	type(z_Hoption)::option_sh,option_A,option_B
	type(z_Hstat)::stats_sh,stats_A,stats_B
	type(z_mesh)::msh_sh,msh_A,msh_B
	type(z_kernelquant)::ker_sh,ker_A,ker_B
	type(z_Bmatrix)::bmat_sh,bmat_A,bmat_B
	integer,allocatable:: groupmembers(:)
	integer nmpi
	type(z_proctree)::ptree_sh,ptree_A,ptree_B
	type(quant_EMCURV),target::quant
	CHARACTER (LEN=1000) DATA_DIR
	integer:: randsize=50
	real(kind=8),allocatable::xyz(:,:)
	integer,allocatable::Permutation(:),tree(:)
	integer Nunk_loc,Maxlevel

	integer maxn, maxnev, maxncv, ldv
	integer iparam(11), ipntr(14)
    logical,allocatable:: select(:)
    complex(kind=8),allocatable:: ax(:), mx(:), d(:), v(:,:), workd(:), workev(:), resid(:), workl(:), eigval(:), eigvec(:,:)
    real(kind=8),allocatable:: rwork(:), rd(:,:)

    character bmattype
	character(len=10) which
    integer ido, n, nx, nev, ncv, lworkl, info, nconv, maxitr, ishfts, mode
    complex(kind=8) sigma
    real(kind=8):: tol, shift_r=0d0, shift_i=0d0
    logical rvec
	real(kind=8),external :: pdznorm2, dlapy2
	integer nargs,flag
	integer v_major,v_minor,v_bugfix

	integer nprow, npcol, myrow, mycol, myArows, myAcols
	integer desca(9)


	! nmpi and groupmembers should be provided by the user
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo

	! generate the process tree
	call z_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree_A)
	deallocate(groupmembers)


	if(ptree_A%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "ButterflyPACK_IE_2D"
	call z_BPACK_GetVersionNumber(v_major,v_minor,v_bugfix)
	write(*,'(A23,I1,A1,I1,A1,I1,A1)') " ButterflyPACK Version:",v_major,".",v_minor,".",v_bugfix
    write(*,*) "   "
	endif

	!**** initialize stats and option
	call z_InitStat(stats_A)
	call z_SetDefaultOptions(option_A)

	!**** intialize the user-defined derived type quant
	quant%RCS_static=1
    quant%RCS_Nsample=2000
	quant%model2d=10
	quant%wavelength=0.08
	quant%freq=1/quant%wavelength/sqrt(BPACK_mu0*BPACK_eps0)
	quant%Nunk=5000

	!**** default parameters for the eigen solvers
	quant%CMmode=0
	quant%SI=0
	quant%shift=(0d0, 0d0)
	quant%nev=1
	quant%tol_eig=1d-13
	quant%which='LM'



	option_A%ErrSol=1
	option_A%format=  HODLR !HMAT!  HODLR !
	option_A%near_para=0.01d0
	option_A%verbosity=2
	option_A%ILU=0
	option_A%forwardN15flag=0
        ! option_A%schulzlevel=0
        ! option_A%LRlevel=100
       ! option_A%level_check=1
    option_A%tol_itersol=1d-5
    ! option_A%sample_para=4d0

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
						if(trim(strings)=='--model2d')then
							read(strings1,*)quant%model2d
						else if	(trim(strings)=='--nunk')then
							read(strings1,*)quant%Nunk
						else if	(trim(strings)=='--wavelength')then
							read(strings1,*)quant%wavelength
							quant%freq=1/quant%wavelength/sqrt(BPACK_mu0*BPACK_eps0)
						else if (trim(strings)=='--freq')then
							read(strings1,*)quant%freq
							quant%wavelength=1/quant%freq/sqrt(BPACK_mu0*BPACK_eps0)
						else if	(trim(strings)=='--cmmode')then
							read(strings1,*)quant%CMmode
						else if	(trim(strings)=='--si')then
							read(strings1,*)quant%SI
						else if	(trim(strings)=='--nev')then
							read(strings1,*)quant%nev
						else if	(trim(strings)=='--tol_eig')then
							read(strings1,*)quant%tol_eig
						else if	(trim(strings)=='--shift_r')then
							read(strings1,*)shift_r
						else if	(trim(strings)=='--shift_i')then
							read(strings1,*)shift_i
						else if	(trim(strings)=='--which')then
							quant%which=trim(strings1)
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

	quant%shift = shift_r + shift_i*BPACK_junit

	call z_PrintOptions(option_A,ptree_A)

    quant%wavenum=2*BPACK_pi/quant%wavelength


   !***********************************************************************
   if(ptree_A%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'frequency:',quant%freq
   write (*,*) 'wavelength:',quant%wavelength
   write (*,*) ''
   endif
   !***********************************************************************


! generate a random matrix and test the dense eigen solvers in scalapack
#if 0
	call z_blacs_gridinfo_wrp(ptree_A%pgrp(1)%ctxt, nprow, npcol, myrow, mycol)
	if (myrow /= -1 .and. mycol /= -1) then
	myArows = z_numroc_wp(quant%Nunk, nbslpk, myrow, 0, nprow)
	myAcols = z_numroc_wp(quant%Nunk, nbslpk, mycol, 0, npcol)
	allocate(matZ1(max(1,myArows), max(1,myAcols)))
	allocate(eigvec(max(1,myArows), max(1,myAcols)))
	allocate(eigval(quant%Nunk))
	call z_descinit_wp(desca, quant%Nunk, quant%Nunk, nbslpk, nbslpk, 0, 0, ptree_A%pgrp(1)%ctxt, max(myArows, 1), info)
	! matZ1=1d0
	call z_RandomMat(max(1,myArows), max(1,myAcols),max(min(myArows,myAcols),1),matZ1,0)

	call z_pgeeigf90(matZ1, quant%Nunk, desca, eigval, eigvec, ptree_A%pgrp(1)%ctxt, ptree_A%pgrp(1)%ctxt_head)

	deallocate(matZ1)
	deallocate(eigvec)
	deallocate(eigval)
	endif
#endif



    !***********************************************************************
	!**** construct compressed A
	!**** geometry generalization and discretization
    call geo_modeling_CURV(quant,ptree_A%Comm)
	!**** register the user-defined function and type in ker_A
	ker_A%QuantApp => quant
	ker_A%FuncZmn => Zelem_EMCURV
	allocate(xyz(2,quant%Nunk))
	do edge=1, quant%Nunk
		xyz(:,edge) = quant%xyz(:,edge*2-1)
	enddo
	allocate(Permutation(quant%Nunk))
	call z_BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat_A,option_A,stats_A,msh_A,ker_A,ptree_A,Coordinates=xyz)
	deallocate(Permutation)
	deallocate(xyz)
	call z_BPACK_construction_Element(bmat_A,option_A,stats_A,msh_A,ker_A,ptree_A)
	!**** print statistics
	call z_PrintStat(stats_A,ptree_A)


! #define DENSEEIGEN

#ifdef DENSEEIGEN
	! call z_FULLMAT_Element(option_A,stats_A,msh_A,ker_A,ptree_sh)
	!***********************************************************************

	call z_BPACK_Convert2Dense(bmat_A,option_A,stats_A,msh_A,ker_A,ptree_A)
	allocate(eigval(quant%Nunk))
	allocate(eigvec(Nunk_loc,quant%Nunk))
	call z_BPACK_Eigen_Dense(bmat_A,option_A,stats_A,msh_A,ker_A,ptree_A,eigval,eigvec)
	deallocate(eigval)
	deallocate(eigvec)

	call z_delete_proctree(ptree_A)
	call z_delete_Hstat(stats_A)
	call z_delete_mesh(msh_A)
	call z_delete_kernelquant(ker_A)
	call z_BPACK_delete(bmat_A)

#else

    !***********************************************************************
	!**** construct compressed A - sigma I or  A - sigma real(A)
	if(quant%SI==1)then
		call z_CopyOptions(option_A,option_sh)
		!**** register the user-defined function and type in ker_sh
		ker_sh%QuantApp => quant
		ker_sh%FuncZmn => Zelem_EMCURV_Shifted
		!**** initialize stats and option
		call z_InitStat(stats_sh)
		!**** create the process tree, can use larger number of mpis if needed
		allocate(groupmembers(nmpi))
		do ii=1,nmpi
			groupmembers(ii)=(ii-1)
		enddo
		call z_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree_sh)
		deallocate(groupmembers)
		allocate(xyz(2,quant%Nunk))
		do edge=1, quant%Nunk
			xyz(:,edge) = quant%xyz(:,edge*2-1)
		enddo
		allocate(Permutation(quant%Nunk))
		call z_BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat_sh,option_sh,stats_sh,msh_sh,ker_sh,ptree_sh,Coordinates=xyz)
		deallocate(Permutation) ! caller can use this permutation vector if needed
		deallocate(xyz)
		!**** computation of the construction phase
		call z_BPACK_construction_Element(bmat_sh,option_sh,stats_sh,msh_sh,ker_sh,ptree_sh)
		!**** factorization phase
		call z_BPACK_Factorization(bmat_sh,option_sh,stats_sh,ptree_sh,msh_sh)
		!**** print statistics
		call z_PrintStat(stats_sh,ptree_sh)
	endif



    !***********************************************************************
	!**** construct compressed real(A)
	if(quant%CMmode==1)then ! solve the characteristic mode
		call z_CopyOptions(option_A,option_B)
		ker_B%QuantApp => quant
		ker_B%FuncZmn => Zelem_EMCURV_Real
		!**** initialize stats and option
		call z_InitStat(stats_B)
		!**** create the process tree, can use larger number of mpis if needed
		allocate(groupmembers(nmpi))
		do ii=1,nmpi
			groupmembers(ii)=(ii-1)
		enddo
		call z_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree_B)
		deallocate(groupmembers)
		allocate(xyz(2,quant%Nunk))
		do edge=1, quant%Nunk
			xyz(:,edge) = quant%xyz(:,edge*2-1)
		enddo
		allocate(Permutation(quant%Nunk))
		call z_BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat_B,option_B,stats_B,msh_B,ker_B,ptree_B,Coordinates=xyz)
		deallocate(Permutation)
		deallocate(xyz)
		call z_BPACK_construction_Element(bmat_B,option_B,stats_B,msh_B,ker_B,ptree_B)
		!**** factorization phase
		if(quant%SI==0)then
			call z_BPACK_Factorization(bmat_B,option_B,stats_B,ptree_B,msh_B)
		endif
		!**** print statistics
		call z_PrintStat(stats_B,ptree_B)
	endif
	!***********************************************************************

    !***********************************************************************
	!**** eigen solver
	allocate(eigval(quant%nev))
	allocate(eigvec(Nunk_loc,quant%nev))
	call z_BPACK_Eigen(bmat_A,option_A,ptree_A, stats_A,bmat_B,option_B,ptree_B, stats_B, bmat_sh,option_sh,ptree_sh,stats_sh, quant%Nunk,Nunk_loc, quant%nev,quant%tol_eig,quant%CMmode,quant%SI,quant%shift,quant%which,nconv,eigval,eigvec)
	deallocate(eigval)
	deallocate(eigvec)


    if(ptree_A%MyID==Main_ID .and. option_A%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

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

	if(quant%CMmode==1)then ! solve the characteristic mode
		call z_delete_proctree(ptree_B)
		call z_delete_Hstat(stats_B)
		call z_delete_mesh(msh_B)
		call z_delete_kernelquant(ker_B)
		call z_BPACK_delete(bmat_B)
	endif

#endif


	!**** deletion of quantities
	call delete_quant_EMCURV(quant)

	call z_blacs_exit_wrp(1)
	call MPI_Finalize(ierr)

end PROGRAM ButterflyPACK_IE_2D





