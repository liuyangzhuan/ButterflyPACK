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
!> @brief This example generates a random LR product, or reads a full matrix from disk, and compress it using entry-valuation-based APIs
!> @details Note that instead of the use of precision dependent subroutine/module/type names "d_", one can also use the following \n
!> #define DAT 1 \n
!> #include "dButterflyPACK_config.fi" \n
!> which will macro replace precision-independent subroutine/module/type names "X" with "d_X" defined in SRC_DOUBLE with double precision

! This exmple works with double precision data
module APPLICATION_MODULE
use d_BPACK_DEFS
implicit none

	!**** define your application-related variables here
	type quant_app
		real(kind=8), allocatable :: matU_glo(:,:),matV_glo(:,:) ! Full Matrix: the random LR matrix to sample its entries
		real(kind=8), allocatable :: matZ_glo(:,:) ! Full Matrix: Full matrix read from files
		real(kind=8), allocatable :: locations(:,:) ! geometrical points
		integer:: rank
		integer:: Nunk
		integer:: Ndim
		real(kind=8):: lambda
		type(d_Bmatrix),pointer::bmat ! Use this metadata in matvec
		type(d_mesh),pointer::msh   ! Use this metadata in matvec
		type(d_proctree),pointer::ptree ! Use this metadata in matvec
		type(d_Hstat),pointer::stats ! Use this metadata in matvec
		type(d_Hoption),pointer::option ! Use this metadata in matvec
		CHARACTER (LEN=1000) DATA_DIR
		CHARACTER (LEN=1000) PERM_DIR
		CHARACTER (LEN=1000) LEAFS_DIR
		CHARACTER (LEN=1000) GEO_DIR
		integer:: tst=1
	end type quant_app

contains

	!**** user-defined subroutine to sample Z_mn as two LR products
	subroutine Zelem_LR(m,n,value_e,quant)
		use d_BPACK_DEFS
		implicit none

		class(*),pointer :: quant
		integer, INTENT(IN):: m,n
		real(kind=8)::value_e
		integer ii

		integer dimn

		select TYPE(quant)
		type is (quant_app)
			value_e = 0
			do ii=1,quant%rank
				value_e = value_e + quant%matU_glo(m,ii)*quant%matV_glo(ii,n)
			enddo
			if(m==n)then
				value_e = value_e + quant%lambda
			endif

			! value_e = quant%matZ_glo(m,n)
		class default
			write(*,*)"unexpected type"
			stop
		end select
	end subroutine Zelem_LR


	!**** user-defined subroutine to sample Z_mn as full matrix
	subroutine Zelem_FULL(m,n,value_e,quant)
		use d_BPACK_DEFS
		implicit none

		class(*),pointer :: quant
		integer, INTENT(IN):: m,n
		real(kind=8)::value_e
		integer ii

		integer dimn

		select TYPE(quant)
		type is (quant_app)
			value_e = quant%matZ_glo(m,n)
		class default
			write(*,*)"unexpected type"
			stop
		end select
	end subroutine Zelem_FULL

	subroutine HODLR_MVP_OneHODLR(trans,Mloc,Nloc,num_vect,Vin,Vout,quant)
		use d_BPACK_DEFS
		use d_MISC_DenseLA
		use d_MISC_Utilities
		use d_BPACK_Solve_Mul
		implicit none
		character trans
		real(kind=8) Vin(:,:),Vout(:,:)
		real(kind=8),allocatable:: Vin_tmp(:,:),Vout_tmp(:,:),Vin_tmp_2D(:,:),Vout_tmp_2D(:,:)
		integer ii,jj,nn,fl_transpose,kk,black_step
		integer, INTENT(in)::Mloc,Nloc,num_vect
		real(kind=8) n1,n2
		! type(d_mesh)::msh
		! type(d_proctree)::ptree
		integer idxs_o,idxe_o,N
		integer nproc,ctxt,info,nb1Dc, nb1Dr, level_p,pgno,num_blocks,ii_new,gg,proc,myi,myj,myAcols,myArows,nprow,npcol,myrow,mycol,Nrow,Ncol
		integer::descsVin(9),descsVout(9),descsMat2D(9),descsVin2D(9),descsVout2D(9)
		class(*),pointer :: quant
		type(d_Bmatrix),pointer::bmat
		! type(d_Hstat)::stats

		select TYPE(quant)
		type is (quant_app)
			pgno=1
			nproc = quant%ptree%pgrp(pgno)%nproc
			bmat=>quant%bmat
			call d_HODLR_Mult(trans,Nloc,num_vect,1,bmat%ho_bf%Maxlevel+1,Vin,Vout,bmat%ho_bf,quant%ptree,quant%option,quant%stats)
		end select

	end subroutine HODLR_MVP_OneHODLR

end module APPLICATION_MODULE





PROGRAM ButterflyPACK_FULL
    use d_BPACK_DEFS
    use APPLICATION_MODULE
	use d_BPACK_Solve_Mul

	use d_BPACK_structure
	use d_BPACK_factor
	use d_BPACK_constr
	use omp_lib
	use d_MISC_Utilities
	use d_BPACK_constr
	use d_BPACK_randomMVP
	use d_BPACK_utilities
    implicit none

    real(kind=8) tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seemyid(50)
	integer times(8)
	real(kind=8) t1,t2,x,y,z,r,theta,phi,error,memory,rtemp1,rtemp2
	real(kind=8),allocatable:: datain(:)

	character(len=:),allocatable  :: string
	character(len=1024)  :: strings,strings1
	character(len=6)  :: info_env
	integer :: length
	integer :: ierr
	integer*8 oldmode,newmode
	type(d_Hoption),target::option,option1
	type(d_Hstat),target::stats,stats1
	type(d_mesh),target::msh,msh1
	type(d_kernelquant),target::ker,ker1
	type(quant_app),target::quant,quant1
	type(d_Bmatrix),target::bmat,bmat1
	integer,allocatable:: groupmembers(:)
	integer nmpi
	integer level,Maxlevel
	type(d_proctree),target::ptree,ptree1
	integer,allocatable::Permutation(:)
	integer Nunk_loc
	integer,allocatable::tree1(:),tree(:)
	integer nargs,flag
	integer v_major,v_minor,v_bugfix,nleaf
	integer*8 blength
	real(kind=8),allocatable::matZ_glo_tmp(:,:)

	!**** nmpi and groupmembers should be provided by the user
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo

	!**** create the process tree
	call d_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	deallocate(groupmembers)


	if(ptree%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "ButterflyPACK_FULL"
	call d_BPACK_GetVersionNumber(v_major,v_minor,v_bugfix)
	write(*,'(A23,I1,A1,I1,A1,I1,A1)') " ButterflyPACK Version:",v_major,".",v_minor,".",v_bugfix
    write(*,*) "   "
	endif

	!**** initialize stats and option
	call d_InitStat(stats)
	call d_SetDefaultOptions(option)

	!**** set solver parameters
	option%nogeo=1  ! no geometry points available
	option%xyzsort=NATURAL ! no reordering will be perfomed (if matrix PSD, can use TM_GRAM as alternative reordering algorithm)
    ! option%verbosity=2
    ! option%LRlevel=100


	quant%DATA_DIR = '../EXAMPLE/FULLMAT_DATA/K05N4096.csv'
	quant%Nunk = 4096
	quant%rank = 2
	quant%Ndim = 0

	nargs = iargc()
	ii=1
	do while(ii<=nargs)
		call getarg(ii,strings)
		if(trim(strings)=='-quant')then ! user-defined quantity parameters !**** read the test number. 1: the kernel is product of two random matrices 2: kernel is a dense matrix stored in file
			flag=1
			do while(flag==1)
				ii=ii+1
				if(ii<=nargs)then
					call getarg(ii,strings)
					if(strings(1:2)=='--')then
						ii=ii+1
						call getarg(ii,strings1)
						if	(trim(strings)=='--tst')then
							read(strings1,*)quant%tst
						else if	(trim(strings)=='--data_dir')then
							quant%data_dir=trim(strings1)
						else if	(trim(strings)=='--perm_dir')then
							quant%perm_dir=trim(strings1)
						else if	(trim(strings)=='--leafs_dir')then
							quant%leafs_dir=trim(strings1)
						else if	(trim(strings)=='--geo_dir')then
							quant%geo_dir=trim(strings1)
							option%nogeo=0  ! geometry points available
						else if	(trim(strings)=='--nunk')then
							read(strings1,*)quant%Nunk
						else if	(trim(strings)=='--ndim')then
							read(strings1,*)quant%Ndim
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
			call d_ReadOption(option,ptree,ii)
		else
			if(ptree%MyID==Main_ID)write(*,*)'ignoring unknown argument: ',trim(strings)
			ii=ii+1
		endif
	enddo

	call d_PrintOptions(option,ptree)

!******************************************************************************!
! generate a LR matrix as two matrix product
	if(quant%tst==1)then
		allocate(tree(1))
		tree=quant%Nunk
		!**** Get matrix size and rank and create the matrix
		quant%lambda = 1d5
		allocate(quant%matU_glo(quant%Nunk,quant%rank))
		call d_RandomMat(quant%Nunk,quant%rank,quant%rank,quant%matU_glo,0)
		call MPI_Bcast(quant%matU_glo,quant%Nunk*quant%rank,MPI_DOUBLE_PRECISION,Main_ID,ptree%Comm,ierr)

		allocate(quant%matV_glo(quant%rank,quant%Nunk))
		call d_RandomMat(quant%rank,quant%Nunk,quant%rank,quant%matV_glo,0)
		call MPI_Bcast(quant%matV_glo,quant%Nunk*quant%rank,MPI_DOUBLE_PRECISION,Main_ID,ptree%Comm,ierr)
	   !***********************************************************************
	   if(ptree%MyID==Main_ID)then
	   write (*,*) ''
	   write (*,*) 'Random LR Kernel computing'
	   write (*,*) 'Matrix size:', quant%Nunk
       write (*,*) ''
	   endif
	   !***********************************************************************

		!**** register the user-defined function and type in ker
		ker%QuantApp => quant
		ker%FuncZmn => Zelem_LR

		!**** initialization of the construction phase
		allocate(Permutation(quant%Nunk))
		call d_BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree)
		deallocate(Permutation) ! caller can use this permutation vector if needed
	endif



!******************************************************************************!
! generate a LR matrix stored in a files
	if(quant%tst==2)then
		allocate(tree(1))
		tree=quant%Nunk

		!**** Get matrix size and rank and create the matrix
		allocate(quant%matZ_glo(quant%Nunk,quant%Nunk))


		!***** assuming reading one row every time
		allocate(datain(quant%Nunk))
		if(ptree%MyID==Main_ID)then


# 358 "FULLMAT_Driver.f90"
		open(10, file=quant%DATA_DIR)
		do ii=1,quant%Nunk
			read(10,*) datain(:)
			quant%matZ_glo(ii,:)=datain(:)
		enddo

		close(10)
		endif
		deallocate(datain)

# 384 "FULLMAT_Driver.f90"

		call d_assert(int(quant%Nunk,kind=8)*int(quant%Nunk,kind=8)< int(2d0**31d0,kind=8),'message size overflow in MPI!')

		call MPI_Bcast(quant%matZ_glo,quant%Nunk*quant%Nunk,MPI_DOUBLE_PRECISION,Main_ID,ptree%Comm,ierr)
		! do ii=1,10
		! call MPI_Bcast(quant%matZ_glo(1,(ii-1)*quant%Nunk/10+1),quant%Nunk*quant%Nunk/10,MPI_DOUBLE_PRECISION,Main_ID,ptree%Comm,ierr)
		! enddo
		! write(*,*)quant%matZ_glo(quant%Nunk,quant%Nunk)

		if(option%nogeo==0)then
			!***** assuming reading one row every time
			allocate(quant%locations(quant%Ndim,quant%Nunk))
			quant%locations=0
			allocate(datain(quant%Nunk))
			open(10, file=quant%GEO_DIR)
			do ii=1,quant%Ndim
				read(10,*) datain(:)
				quant%locations(ii,:)=datain(:)
			enddo
			close(10)
			deallocate(datain)
		endif

	   !***********************************************************************
	   if(ptree%MyID==Main_ID)then
	   write (*,*) ''
	   write (*,*) 'FullMat computing'
	   write (*,*) 'Matrix size:', quant%Nunk
	   write (*,*) ''
	   endif
	   !***********************************************************************

		!**** register the user-defined function and type in ker
		ker%QuantApp => quant
		ker%FuncZmn => Zelem_FULL

		!**** initialization of the construction phase
		allocate(Permutation(quant%Nunk))
		call d_BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree,Coordinates=quant%locations,tree=tree)
		deallocate(Permutation) ! caller can use this permutation vector if needed
		deallocate(tree)
	endif

!******************************************************************************!


	!**** computation of the construction phase
    call d_BPACK_construction_Element(bmat,option,stats,msh,ker,ptree)

	if(quant%tst==2)deallocate(quant%matZ_glo)


	!**** factorization phase
    call d_BPACK_Factorization(bmat,option,stats,ptree,msh)


	!**** solve phase
	if(option%ErrSol==1)then
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Test Solve ......"
		call d_BPACK_Test_Solve_error(bmat,msh%idxe-msh%idxs+1,option,ptree,stats)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Test Solve finished"
	endif
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "

	!**** print statistics
	call d_PrintStat(stats,ptree)


	!**** construct the second HODLR using the first HODLR as matvec
	call d_CopyOptions(option,option1)
	option1%nogeo=1
	option1%xyzsort=NATURAL
	ker1%FuncHMatVec=>HODLR_MVP_OneHODLR
	ker1%QuantApp=>quant1
	quant1%bmat=>bmat
	quant1%msh=>msh
	quant1%ptree=>ptree
	quant1%stats=>stats
	quant1%option=>option
	quant1%Nunk=quant%Nunk
	msh1%Nunk = msh%Nunk



	!**** initialize stats and option
	call d_InitStat(stats1)

	!**** create the process tree, can use larger number of mpis if needed
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo
	call d_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree1)
	deallocate(groupmembers)



	!**** use the clustering tree from the first HODLR
	select case(option%format)
	case(HODLR)
		Maxlevel=bmat%ho_bf%Maxlevel
	case(HMAT)
		Maxlevel=bmat%h_mat%Maxlevel
	end select
	allocate (tree1(2**Maxlevel))
	do ii=1,2**Maxlevel
		tree1(ii)=msh%basis_group(2**Maxlevel+ii-1)%tail-msh%basis_group(2**Maxlevel+ii-1)%head+1
	enddo


	!**** initialization of the construction phase
	allocate(Permutation(quant1%Nunk))
	call d_BPACK_construction_Init(quant1%Nunk,Permutation,Nunk_loc,bmat1,option1,stats1,msh1,ker1,ptree1,tree=tree1)
	deallocate(Permutation) ! caller can use this permutation vector if needed
	deallocate(tree1)


	!**** computation of the construction phase
	option%less_adapt=0
	call d_BPACK_construction_Matvec(bmat1,d_matvec_user,Memory,error,option1,stats1,ker1,ptree1,msh1)

	!**** print statistics
	call d_PrintStat(stats1,ptree1)


	!**** deletion of quantities
	call d_delete_proctree(ptree1)
	call d_delete_Hstat(stats1)
	call d_delete_mesh(msh1)
	call d_delete_kernelquant(ker1)
	call d_BPACK_delete(bmat1)

	if(allocated(quant%matU_glo))deallocate(quant%matU_glo)
	if(allocated(quant%matV_glo))deallocate(quant%matV_glo)
	if(allocated(quant%matZ_glo))deallocate(quant%matZ_glo)

	call d_delete_proctree(ptree)
	call d_delete_Hstat(stats)
	call d_delete_mesh(msh)
	call d_delete_kernelquant(ker)
	call d_BPACK_delete(bmat)

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	call blacs_exit(1)
	call MPI_Finalize(ierr)

end PROGRAM ButterflyPACK_FULL



