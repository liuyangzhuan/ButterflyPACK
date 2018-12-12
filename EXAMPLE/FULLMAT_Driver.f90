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

! This exmple works with double precision data
#define DAT 1

#include "ButterflyPACK_config.fi"

module APPLICATION_MODULE
use BPACK_DEFS
implicit none

	!**** define your application-related variables here
	type quant_app
		real(kind=8), allocatable :: matU_glo(:,:),matV_glo(:,:) ! Full Matrix: the random LR matrix to sample its entries
		real(kind=8), allocatable :: matZ_glo(:,:) ! Full Matrix: Full matrix read from files
		integer:: rank
		integer:: Nunk
		real(kind=8):: lambda
		type(Bmatrix),pointer::bmat ! Use this metadata in matvec
		type(mesh),pointer::msh   ! Use this metadata in matvec
		type(proctree),pointer::ptree ! Use this metadata in matvec
		type(Hstat),pointer::stats ! Use this metadata in matvec
		type(Hoption),pointer::option ! Use this metadata in matvec
	end type quant_app

contains

	!**** user-defined subroutine to sample Z_mn as two LR products
	subroutine Zelem_LR(m,n,value_e,quant)
		use BPACK_DEFS
		implicit none

		class(*),pointer :: quant
		integer, INTENT(IN):: m,n
		real(kind=8)::value_e
		integer ii

		real(kind=8) r_mn
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
		use BPACK_DEFS
		implicit none

		class(*),pointer :: quant
		integer, INTENT(IN):: m,n
		real(kind=8)::value_e
		integer ii

		real(kind=8) r_mn
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
		use BPACK_DEFS
		use DenseLA
		use misc
		use BPACK_Solve_Mul
		implicit none
		character trans
		real(kind=8) Vin(:,:),Vout(:,:)
		real(kind=8),allocatable:: Vin_tmp(:,:),Vout_tmp(:,:),Vin_tmp_2D(:,:),Vout_tmp_2D(:,:)
		real(kind=8) ctemp,a,b
		integer ii,jj,nn,fl_transpose,kk,black_step
		integer, INTENT(in)::Mloc,Nloc,num_vect
		real(kind=8) n1,n2,tmp(2)
		! type(mesh)::msh
		! type(proctree)::ptree
		integer idxs_o,idxe_o,N
		integer nproc,ctxt,info,nb1Dc, nb1Dr, level_p,pgno,num_blocks,ii_new,gg,proc,myi,myj,myAcols,myArows,nprow,npcol,myrow,mycol,Nrow,Ncol
		integer::descsVin(9),descsVout(9),descsMat2D(9),descsVin2D(9),descsVout2D(9)
		class(*),pointer :: quant
		type(Bmatrix),pointer::bmat
		! type(Hstat)::stats

		select TYPE(quant)
		type is (quant_app)
			pgno=1
			nproc = quant%ptree%pgrp(pgno)%nproc
			bmat=>quant%bmat
			call HODLR_Mult(trans,Nloc,num_vect,1,bmat%ho_bf%Maxlevel+1,Vin,Vout,bmat%ho_bf,quant%ptree,quant%option,quant%stats)
		end select

	end subroutine HODLR_MVP_OneHODLR

end module APPLICATION_MODULE





PROGRAM ButterflyPACK_FULL
    use BPACK_DEFS
    use APPLICATION_MODULE
	use BPACK_Solve_Mul

	use BPACK_structure
	use BPACK_factor
	use BPACK_constr
	use omp_lib
	use misc
	use BPACK_constr
	use BPACK_randomMVP
    implicit none

    real(kind=8) para
    real(kind=8) tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seemyid(50)
	integer times(8)
	real(kind=8) t1,t2,x,y,z,r,theta,phi,error,memory
	real(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:)
	real(kind=8),allocatable:: datain(:)

	character(len=:),allocatable  :: string
	character(len=1024)  :: strings
	character(len=6)  :: info_env
	integer :: length
	integer :: ierr
	integer*8 oldmode,newmode
	type(Hoption),target::option,option1
	type(Hstat),target::stats,stats1
	type(mesh),target::msh,msh1
	type(kernelquant),target::ker,ker1
	type(quant_app),target::quant,quant1
	type(Bmatrix),target::bmat,bmat1
	integer,allocatable:: groupmembers(:)
	integer nmpi
	integer level,Maxlevel,N_unk_loc
	type(proctree),target::ptree,ptree1
	CHARACTER (LEN=1000) DATA_DIR
	integer:: tst=1
	integer,allocatable::Permutation(:)
	integer Nunk_loc
	
	!**** nmpi and groupmembers should be provided by the user
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo

	!**** create the process tree
	call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	deallocate(groupmembers)


	if(ptree%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "ButterflyPACK_FULL"
    write(*,*) "   "
	endif

	!**** initialize stats and option
	call InitStat(stats)
	call SetDefaultOptions(option)
    
	!**** set solver parameters
	option%nogeo=1  ! no geometry points available
	option%xyzsort=NATURAL ! no reordering will be perfomed (if matrix PSD, can use TM_GRAM as alternative reordering algorithm)
    ! option%verbosity=2
    ! option%LRlevel=100
	

	!**** read the test number. 1: the kernel is product of two random matrices 2: kernel is a dense matrix stored in file
	if(iargc()>=1)then
		call getarg(1,strings)
		read(strings,*)tst
	endif

!******************************************************************************!
! generate a LR matrix as two matrix product
	if(tst==1)then
		!**** Get matrix size and rank and create the matrix
		quant%Nunk = 10000
		quant%rank = 2
		quant%lambda = 1d5
		allocate(quant%matU_glo(quant%Nunk,quant%rank))
		call RandomMat(quant%Nunk,quant%rank,quant%rank,quant%matU_glo,0)
		call MPI_Bcast(quant%matU_glo,quant%Nunk*quant%rank,MPI_DOUBLE_PRECISION,Main_ID,ptree%Comm,ierr)

		allocate(quant%matV_glo(quant%rank,quant%Nunk))
		call RandomMat(quant%rank,quant%Nunk,quant%rank,quant%matV_glo,0)
		call MPI_Bcast(quant%matV_glo,quant%Nunk*quant%rank,MPI_DOUBLE_PRECISION,Main_ID,ptree%Comm,ierr)
	   !***********************************************************************
	   if(ptree%MyID==Main_ID)then
	   write (*,*) ''
	   write (*,*) 'FullMat computing'
	   write (*,*) ''
	   endif
	   !***********************************************************************
	   
		!**** initialization of the construction phase
		allocate(Permutation(quant%Nunk))
		call BPACK_construction_Element_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree,Zelem_LR,quant)
		deallocate(Permutation) ! caller can use this permutation vector if needed 		   
	   
	endif



!******************************************************************************!
! generate a LR matrix stored in a files
	if(tst==2)then
		strings = '../EXAMPLE/FULLMAT_DATA/K05N4096.csv'
		if(iargc()>=2)then
			call getarg(2,strings)
		endif

		!**** Get matrix size and rank and create the matrix
		quant%Nunk = 4096
		allocate(quant%matZ_glo(quant%Nunk,quant%Nunk))
		allocate(datain(quant%Nunk))
		open(10, file=strings)
		do ii=1,quant%Nunk
			read(10,*) datain(:)
			quant%matZ_glo(:,ii)=datain
		enddo
		close(10)
		call MPI_Bcast(quant%matZ_glo,quant%Nunk*quant%Nunk,MPI_DOUBLE_PRECISION,Main_ID,ptree%Comm,ierr)
	   !***********************************************************************
	   if(ptree%MyID==Main_ID)then
	   write (*,*) ''
	   write (*,*) 'Random LR Kernel computing'
	   write (*,*) ''
	   endif
	   !***********************************************************************
	   
		!**** initialization of the construction phase
		allocate(Permutation(quant%Nunk))
		call BPACK_construction_Element_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree,Zelem_FULL,quant)
		deallocate(Permutation) ! caller can use this permutation vector if needed 		
	   
	endif

!******************************************************************************!


	!**** computation of the construction phase
	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Matrix construction......"
    call BPACK_construction_Element_Compute(bmat,option,stats,msh,ker,element_Zmn_user,ptree)

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Matrix construction finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
 	t2 = OMP_get_wtime()


	!**** factorization phase
	if(option%precon/=NOPRECON)then
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Factor ......"
    call BPACK_Factorization(bmat,option,stats,ptree,msh)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Factor finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	end if

	!**** solve phase
	if(option%ErrSol==1)then
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Test Solve ......"
		call BPACK_Test_Solve_error(bmat,msh%idxe-msh%idxs+1,option,ptree,stats)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Test Solve finished"
	endif
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	
	!**** print statistics
	call PrintStat(stats,ptree)


	!**** construct the second HODLR using the first HODLR as matvec
	call CopyOptions(option,option1)
	option1%nogeo=1
	option1%xyzsort=NATURAL
	ker1%FuncHMatVec=>HODLR_MVP_OneHODLR
	ker1%QuantApp=>quant1
	quant1%bmat=>bmat
	quant1%msh=>msh
	quant1%ptree=>ptree
	quant1%stats=>stats
	quant1%option=>option
	msh1%Nunk = msh%Nunk


	call InitStat(stats1)

	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo
	call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree1)
	deallocate(groupmembers)


	select case(option%format)
	case(HODLR)
		allocate (msh1%pretree(2**bmat%ho_bf%Maxlevel))
		do ii=1,2**bmat%ho_bf%Maxlevel
			msh1%pretree(ii)=msh%basis_group(2**bmat%ho_bf%Maxlevel+ii-1)%tail-msh%basis_group(2**bmat%ho_bf%Maxlevel+ii-1)%head+1
		enddo
	case(HMAT)
		allocate (msh1%pretree(2**bmat%h_mat%Maxlevel))
		do ii=1,2**bmat%h_mat%Maxlevel
			msh1%pretree(ii)=msh%basis_group(2**bmat%h_mat%Maxlevel+ii-1)%tail-msh%basis_group(2**bmat%h_mat%Maxlevel+ii-1)%head+1
		enddo
	end select

    call Cluster_partition(bmat1,option1,msh1,ker1,element_Zmn_user,ptree1)
	call BPACK_structuring(bmat1,option1,msh1,ptree1,stats)

	N_unk_loc = msh1%idxe-msh1%idxs+1
	t1 = OMP_get_wtime()
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "FastMATVEC-based Matrix construction......"
	call HODLR_randomized(bmat1%ho_bf,matvec_user,N_unk_loc,Memory,error,option1,stats1,ker1,ptree1,msh1)
	t2 = OMP_get_wtime()
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "FastMATVEC-based Matrix construction finished",t2-t1, 'secnds. Error: ', error

	!**** print statistics
	call PrintStat(stats1,ptree1)
	
	
	!**** deletion of quantities
	call delete_proctree(ptree1)
	call delete_Hstat(stats1)
	call delete_mesh(msh1)
	call delete_kernelquant(ker1)
	call BPACK_delete(bmat1)

	if(allocated(quant%matU_glo))deallocate(quant%matU_glo)
	if(allocated(quant%matV_glo))deallocate(quant%matV_glo)
	if(allocated(quant%matZ_glo))deallocate(quant%matZ_glo)

	call delete_proctree(ptree)
	call delete_Hstat(stats)
	call delete_mesh(msh)
	call delete_kernelquant(ker)
	call BPACK_delete(bmat)

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	call blacs_exit(1)
	call MPI_Finalize(ierr)

end PROGRAM ButterflyPACK_FULL


