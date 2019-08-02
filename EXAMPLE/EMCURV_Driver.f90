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


PROGRAM ButterflyPACK_IE_2D
    use BPACK_DEFS
    use EMCURV_MODULE

	use BPACK_structure
	use BPACK_Solve_Mul
	use BPACK_factor
	use BPACK_constr
	use BPACK_Utilities
	use omp_lib
	use MISC_Utilities
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
	complex(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:)

	character(len=:),allocatable  :: string
	character(len=1024)  :: strings,strings1
	character(len=6)  :: info_env
	integer :: length
	integer :: ierr
	integer*8 oldmode,newmode
	type(Hoption),target::option,option1
	type(Hstat),target::stats,stats1
	type(mesh),target::msh,msh1
	type(kernelquant),target::ker,ker1
	type(Bmatrix),target::bmat,bmat1
	integer,allocatable:: groupmembers(:)
	integer nmpi
	type(proctree),target::ptree,ptree1
	type(quant_EMCURV),target::quant
	type(quant_bmat),target::quant1
	CHARACTER (LEN=1000) DATA_DIR
	integer:: randsize=50
	real(kind=8),allocatable::xyz(:,:)
	integer,allocatable::Permutation(:),tree(:)
	integer Nunk_loc,Maxlevel
	integer nargs,flag
	integer v_major,v_minor,v_bugfix

	! nmpi and groupmembers should be provided by the user
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo

	! generate the process tree
	call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	deallocate(groupmembers)


	if(ptree%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "ButterflyPACK_IE_2D"
	call BPACK_GetVersionNumber(v_major,v_minor,v_bugfix)
	write(*,'(A23,I1,A1,I1,A1,I1,A1)') " ButterflyPACK Version:",v_major,".",v_minor,".",v_bugfix
    write(*,*) "   "
	endif

	!**** initialize stats and option
	call InitStat(stats)
	call SetDefaultOptions(option)

	!**** intialize the user-defined derived type quant
	quant%RCS_static=1
    quant%RCS_Nsample=1
	quant%model2d=10
	quant%wavelength=0.08
	quant%freq=1/quant%wavelength/sqrt(mu0*eps0)
	quant%Nunk=5000

	option%ErrSol=1
	option%format=  HODLR !HMAT!  HODLR !
	option%near_para=0.01d0
	option%verbosity=2
	option%ILU=0
	option%forwardN15flag=0
        ! option%schulzlevel=0
        ! option%LRlevel=100
       ! option%level_check=1
    option%tol_itersol=1d-5
    ! option%sample_para=4d0

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
							quant%freq=1/quant%wavelength/sqrt(mu0*eps0)
						else if (trim(strings)=='--freq')then
							read(strings1,*)quant%freq
							quant%wavelength=1/quant%freq/sqrt(mu0*eps0)
						else
							if(ptree%MyID==Main_ID)write(*,*)'ignoring unknown quant: ', trim(strings)
						endif
					else
						flag=0
					endif
				else
					flag=0
				endif
			enddo
		else if(trim(strings)=='-option')then ! options of ButterflyPACK
			call ReadOption(option,ptree,ii)
		else
			if(ptree%MyID==Main_ID)write(*,*)'ignoring unknown argument: ',trim(strings)
			ii=ii+1
		endif
	enddo




    quant%wavenum=2*pi/quant%wavelength


   !***********************************************************************
   if(ptree%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'frequency:',quant%freq
   write (*,*) 'wavelength:',quant%wavelength
   write (*,*) ''
   endif
   !***********************************************************************

	!**** geometry generalization and discretization
    call geo_modeling_CURV(quant,ptree%Comm)

	option%touch_para = 3* quant%minedgelength

	!**** register the user-defined function and type in ker
	ker%QuantApp => quant
	ker%FuncZmn => Zelem_EMCURV

	!**** initialization of the construction phase
	allocate(xyz(2,quant%Nunk))
	do edge=1, quant%Nunk
		xyz(:,edge) = quant%xyz(:,edge*2-1)
	enddo
    allocate(Permutation(quant%Nunk))
	call PrintOptions(option,ptree)
	call BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree,Coordinates=xyz)
	deallocate(Permutation) ! caller can use this permutation vector if needed
	deallocate(xyz)

	!**** computation of the construction phase
    call BPACK_construction_Element(bmat,option,stats,msh,ker,ptree)


	select case(option%format)
	case(HODLR)
		Maxlevel=bmat%ho_bf%Maxlevel
	case(HMAT)
		Maxlevel=bmat%h_mat%Maxlevel
	end select
	if(option%lnoBP>Maxlevel)then	 ! haven't implement the following for Bplus.

		!*********** Construct the second HODLR by using the first HODLR as parallel element extraction
		call CopyOptions(option,option1)
		option1%nogeo=1   ! this indicates the second HOLDR construction requires no geometry information
		option1%xyzsort=NATURAL ! this indicates the second HOLDR construction requires no reordering
		option1%elem_extract=1 ! this indicates the second HOLDR construction uses user-defined parallel element extraction

		!**** register the user-defined function and type in ker
		ker1%FuncZmn=>Zelem_EMCURV
		ker1%FuncZmnBlock=>Zelem_block_Extraction
		ker1%QuantApp=>quant1

		quant1%bmat=>bmat
		quant1%msh=>msh
		quant1%ptree=>ptree
		quant1%stats=>stats
		quant1%option=>option

		msh1%Nunk = msh%Nunk

		allocate(xyz(2,quant%Nunk))
		do edge=1, quant%Nunk
			xyz(:,edge) = quant%xyz(:,edge*2-1)
		enddo
		allocate(msh1%xyz(2,quant%Nunk)) ! xyz is used to compute group centers to determine admissibility
		do edge=1,quant%Nunk
			msh1%xyz(:,edge) = xyz(:,msh%new2old(edge))
		enddo
		deallocate(xyz)

		call InitStat(stats1)

		! !**** generate the process tree for the second HODLR, can use larger number of MPIs if you want to
		! allocate(groupmembers(nmpi))
		! do ii=1,nmpi
			! groupmembers(ii)=(ii-1)
		! enddo
		! call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree1)
		! deallocate(groupmembers)


		!**** use the clustering tree from the first HODLR
		allocate (tree(2**Maxlevel))
		do ii=1,2**Maxlevel
			tree(ii)=msh%basis_group(2**Maxlevel+ii-1)%tail-msh%basis_group(2**Maxlevel+ii-1)%head+1
		enddo

		!**** initialization of the construction phase
		allocate(Permutation(msh1%Nunk))
		call BPACK_construction_Init(msh1%Nunk,Permutation,Nunk_loc,bmat1,option1,stats1,msh1,ker1,ptree,tree=tree)
		deallocate(Permutation) ! caller can use this permutation vector if needed
		deallocate(tree)

		!**** computation of the construction phase
		call BPACK_construction_Element(bmat1,option1,stats1,msh1,ker1,ptree)


		!**** deletion of quantities
		! call delete_proctree(ptree1)
		call delete_Hstat(stats1)
		call delete_mesh(msh1)
		call delete_kernelquant(ker1)
		call BPACK_delete(bmat1)
	endif



	!**** factorization phase
    call BPACK_Factorization(bmat,option,stats,ptree,msh)


	!**** solve phase
   ! call EM_solve_CURV(bmat,option,msh,quant,ptree,stats)



	!**** print statistics
	call PrintStat(stats,ptree)


	!**** deletion of quantities
	call delete_quant_EMCURV(quant)
	call delete_proctree(ptree)
	call delete_Hstat(stats)
	call delete_mesh(msh)
	call delete_kernelquant(ker)
	call BPACK_delete(bmat)




    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	call blacs_exit(1)
	call MPI_Finalize(ierr)

end PROGRAM ButterflyPACK_IE_2D




