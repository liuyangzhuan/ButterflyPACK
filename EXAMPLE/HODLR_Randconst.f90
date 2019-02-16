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


module APPLICATION_MODULE
use BPACK_DEFS
implicit none

	!**** define your application-related variables here
	type quant_app
		integer::Nunk ! matrix size
		complex(kind=8), allocatable :: matZ_glo(:,:) ! Full Matrix: Full matrix read from files
		complex(kind=8), allocatable :: matZ_loc(:,:) ! Local Matrix: Loccal matrix in a npx1 blasc grid
		integer,pointer :: N_p(:,:)=>null() ! column sizes of all processes sharing this hodlr
		type(Bmatrix),pointer::bmat ! Use this metadata in matvec
		type(mesh),pointer::msh   ! Use this metadata in matvec
		type(proctree),pointer::ptree ! Use this metadata in matvec
		type(Hstat),pointer::stats ! Use this metadata in matvec
		type(Hoption),pointer::option ! Use this metadata in matvec
	end type quant_app

contains

	!**** user-defined subroutine to sample Z_mn as full matrix
	subroutine Zelem_FULL(m,n,value_e,quant)
		use BPACK_DEFS
		implicit none

		class(*),pointer :: quant
		integer, INTENT(IN):: m,n
		complex(kind=8)::value_e
		integer ii

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
		complex(kind=8) Vin(:,:),Vout(:,:)
		complex(kind=8),allocatable:: Vin_tmp(:,:),Vout_tmp(:,:),Vin_tmp_2D(:,:),Vout_tmp_2D(:,:)
		complex(kind=8) ctemp,a,b
		integer ii,jj,nn,fl_transpose,kk,black_step
		integer, INTENT(in)::Mloc,Nloc,num_vect
		real(kind=8) n1,n2,tmp(2)
		integer idxs_o,idxe_o,N
		integer nproc,ctxt,info,nb1Dc, nb1Dr, level_p,pgno,num_blocks,ii_new,gg,proc,myi,myj,myAcols,myArows,nprow,npcol,myrow,mycol,Nrow,Ncol
		integer::descsVin(9),descsVout(9),descsMat2D(9),descsVin2D(9),descsVout2D(9)
		class(*),pointer :: quant
		type(Bmatrix),pointer::bmat

		select TYPE(quant)
		type is (quant_app)
			pgno=1
			nproc = quant%ptree%pgrp(pgno)%nproc

			bmat=>quant%bmat
			call BPACK_Mult(trans,Nloc,num_vect,Vin,Vout,bmat,quant%ptree,quant%option,quant%stats)
		end select

	end subroutine HODLR_MVP_OneHODLR

	subroutine HODLR_MVP_Fullmat(trans,Mloc,Nloc,num_vect,Vin,Vout,quant)
		use BPACK_DEFS
		use DenseLA
		use misc
		implicit none
		character trans
		complex(kind=8) Vin(:,:),Vout(:,:)
		complex(kind=8),allocatable:: Vin_tmp(:,:),Vout_tmp(:,:),Vin_tmp_2D(:,:),Vout_tmp_2D(:,:)
		complex(kind=8) ctemp,a,b
		integer ii,jj,nn,fl_transpose,kk,black_step
		integer, INTENT(in)::Mloc,Nloc,num_vect
		real(kind=8) n1,n2,tmp(2)
		integer idxs_o,idxe_o,N
		integer nproc,ctxt,info,nb1Dc, nb1Dr, level_p,pgno,num_blocks,ii_new,gg,proc,myi,myj,myAcols,myArows,nprow,npcol,myrow,mycol,Nrow,Ncol
		integer::descsVin(9),descsVout(9),descsMat2D(9),descsVin2D(9),descsVout2D(9)
		class(*),pointer :: quant

		select TYPE(quant)
		type is (quant_app)

		pgno=1
		nproc = quant%ptree%pgrp(pgno)%nproc
		N = quant%N_p(nproc,2)


		!!!!**** generate 2D grid blacs quantities
		ctxt = quant%ptree%pgrp(pgno)%ctxt
		call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
		if(myrow/=-1 .and. mycol/=-1)then
			myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
			myAcols = numroc_wp(num_vect, nbslpk, mycol, 0, npcol)
			call descinit( descsVin2D, N, num_vect, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
			call descinit( descsVout2D, N, num_vect, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
			myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
			myAcols = numroc_wp(N, nbslpk, mycol, 0, npcol)
			call descinit( descsMat2D, N, N, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
			allocate(Vin_tmp_2D(myArows,myAcols))
			allocate(Vout_tmp_2D(myArows,myAcols))
			Vout_tmp_2D=0
		else
			descsVin2D(2)=-1
			descsVout2D(2)=-1
			descsMat2D(2)=-1
		endif


		!!!!**** redistribution of input vectors
		call Redistribute1Dto2D(Vin,quant%N_p,0,pgno,Vin_tmp_2D,N,0,pgno,num_vect,quant%ptree)


		!!!!**** perform gemm on 2d grid
		if(myrow/=-1 .and. mycol/=-1)then
			call pgemmf90(trans,'N',N,num_vect,N,cone, quant%matZ_loc,1,1,descsMat2D,Vin_tmp_2D,1,1,descsVin2D,czero,Vout_tmp_2D,1,1,descsVout2D)
		endif


		!!!!**** redistribution of output vectors
		call Redistribute2Dto1D(Vout_tmp_2D,N,0,pgno,Vout,quant%N_p,0,pgno,num_vect,quant%ptree)


		!!!!**** deallocation buffers
		if(myrow/=-1 .and. mycol/=-1)then
			deallocate(Vin_tmp_2D)
			deallocate(Vout_tmp_2D)
		endif


		end select

	end subroutine HODLR_MVP_Fullmat


	subroutine CreateDistDenseMat(N,msh,ptree,quant)
		use BPACK_DEFS
		use DenseLA
		use misc
		implicit none
		complex(kind=8),allocatable:: Vin_tmp(:,:),Vout_tmp(:,:)
		complex(kind=8) ctemp,a,b
		integer ii,jj,nn,fl_transpose,kk,black_step
		integer, INTENT(in)::N
		real(kind=8) n1,n2,tmp(2)
		type(mesh)::msh
		type(proctree)::ptree
		type(quant_app) :: quant
		integer nproc,ctxt,nb1Dc, nb1Dr, level_p,pgno,num_blocks,ii_new,gg,proc,myi,myj,myAcols,myArows,nprow,npcol,myrow,mycol


		pgno=1
		nproc = ptree%pgrp(pgno)%nproc

		!!!!****** allocate index array for 1D HODLR layout
		level_p = ptree%nlevel-GetTreelevel(pgno)
		num_blocks = 2**level_p
		allocate(quant%N_p(nproc,2))
		quant%N_p(:,1) = N+1
		quant%N_p(:,2) = -N-1
		do ii=1,num_blocks
			ii_new=ii
			gg = 2**level_p+ii_new-1
			proc = ptree%pgrp(pgno*2**level_p+ii-1)%head - ptree%pgrp(pgno)%head
			quant%N_p(proc+1,1) = min(quant%N_p(proc+1,1),msh%basis_group(gg)%head)
			quant%N_p(proc+1,2) = max(quant%N_p(proc+1,2),msh%basis_group(gg)%tail)
		enddo


		!!!!****** assemble the full matrix in 1D blacs layout
		ctxt = ptree%pgrp(pgno)%ctxt
		call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
		if(myrow/=-1 .and. mycol/=-1)then
		myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(N, nbslpk, mycol, 0, npcol)
		allocate(quant%matZ_loc(myArows,myAcols))
		quant%matZ_loc=0
		do myi=1,myArows
			call l2g(myi,myrow,N,nprow,nbslpk,ii)
			do myj=1,myAcols
				call l2g(myj,mycol,N,npcol,nbslpk,jj)
				quant%matZ_loc(myi,myj) = quant%matZ_glo(msh%new2old(ii),msh%new2old(jj))
			enddo
		enddo
		endif

	end subroutine CreateDistDenseMat




end module APPLICATION_MODULE


PROGRAM ButterflyPACK_ScatteringMatrix_Matvec

    use BPACK_DEFS

	use BPACK_structure
	use BPACK_factor
	use BPACK_constr
	use omp_lib
	use Bplus_compress
	use BPACK_randomMVP
	use APPLICATION_MODULE

    implicit none

    real(kind=8) para,error,tmp1,tmp2,norm1,norm2
    real(kind=8) tolerance,rankrate
    integer Primary_block, nn, mm
    integer i,j,k, threads_num,ii,jj
	integer seed_myid(50)
	integer times(8)
	real(kind=8) t1,t2,t3,t4,x,y,z,r,theta,phi,tmp(3),Memory
	complex(kind=8),allocatable:: InputVec(:)
	complex(kind=8):: ctemp
	integer kk,black_step,rank0
	complex(kind=8),allocatable::Vout1(:,:),Vout2(:,:),Vin(:,:)
	character(len=1024)  :: strings
	type(Hoption),target:: option,option1
	type(Hstat),target::stats,stats1
	type(mesh),target::msh,msh1
	type(kernelquant),target::ker,ker1
	integer:: explicitflag
	type(Bmatrix),target::bmat,bmat1
	integer Nin1,Nout1,Nin2,Nout2
	type(proctree),target::ptree,ptree1
	integer,allocatable:: groupmembers(:)
	integer :: ierr
	integer :: nmpi
	CHARACTER (LEN=1000) DATA_DIR
	type(quant_app),target::quant,quant1
	integer N_unk_loc,Maxlevel
	integer,allocatable::tree(:),Permutation(:)
	real(kind=8),allocatable::xyz(:,:)
	integer Nunk_loc

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

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------Program Start----------------------------------"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "ButterflyPACK_ScatteringMatrix_Matvec"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "   "

	!**** initialize stats and option
	call InitStat(stats)
	call SetDefaultOptions(option)

	!**** intialize the user-defined derived type quant
	option%nogeo=0
	option%xyzsort=CKD


	explicitflag=0
	!**** initialize the user-defined derived type quant
	!*********** Construct the first HODLR by read-in the full matrix and (if explicitflag=0) use it as a matvec or (if explicitflag=1) use it as a fast entry evaluation
	DATA_DIR='../EXAMPLE/FULLMAT_DATA'
	if(iargc()>=1)then
		CALL getarg(1, DATA_DIR)
	endif


	!**** predefine the first three levels of tree due to the physical meanings
	quant%ptree=>ptree
	quant%Nunk = 3720
	Nin1 = 320*2
	Nout1 = 610*2
	Nin2 = 320*2
	Nout2 = 610*2
	call assert(Nin1+Nout1+Nin2+Nout2==quant%Nunk,'The two surfaces have mismatched number of unknowns')
	allocate(tree(4))
	tree(1) = Nin1
	tree(2) = Nout1
	tree(3) = Nin2
	tree(4) = Nout2


	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'Blackbox HODLR for scattering matrix compression'
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,'(A11,I9)')' Nsurface: ',quant%Nunk


	!**** generate the list of confidantes for clustering. For simplicity, duplicate locations of each point
	allocate(xyz(3,quant%Nunk))
	open(unit=521,file=trim(DATA_DIR)//'/Smatrix.geo',status='old')
	do kk=1,quant%Nunk/2
		read(521,*) xyz(1,2*kk-1),xyz(2,2*kk-1),xyz(3,2*kk-1)
		xyz(:,2*kk) = xyz(:,2*kk-1)
	end do
	close(521)


	!**** generate the full matrix used for entry evaluation function Zelem_FULL
	t1 = OMP_get_wtime()
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Generating fullmat ......"
	allocate(quant%matZ_glo(quant%Nunk,quant%Nunk))
	quant%matZ_glo = 0
	open(unit=888,file=trim(DATA_DIR)//'/Smatrix.mat',status='old')
	do ii=1,quant%Nunk
	do kk=1,quant%Nunk
		read(888,*)tmp(1),tmp(2)
		quant%matZ_glo(kk,ii) = cmplx(tmp(1),tmp(2),dp)
	end do
	end do
	close(unit=888)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Generating fullmat finished"
	t2 = OMP_get_wtime()
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)t2-t1, 'secnds'




	if(explicitflag ==1)then

		!**** register the user-defined function and type in ker
		ker%QuantApp => quant
		ker%FuncZmn => Zelem_FULL

		!**** initialization of the construction phase
	    allocate(Permutation(quant%Nunk))
		call BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree,Coordinates=xyz,tree=tree)
		deallocate(Permutation) ! caller can use this permutation vector if needed
		deallocate(xyz)
		deallocate(tree)

		!**** define other quantities in quant using information returned by BPACK_construction_Init
		call CreateDistDenseMat(quant%Nunk,msh,ptree,quant)

		!**** computation of the construction phase
		call BPACK_construction_Element(bmat,option,stats,msh,ker,element_Zmn_user,ptree)



		!**** check error of the entire construction
		N_unk_loc = msh%idxe-msh%idxs+1
		allocate(Vin(N_unk_loc,1))
		allocate(Vout1(N_unk_loc,1))
		allocate(Vout2(N_unk_loc,1))
		Vout2=0
		do ii=1,N_unk_loc
			call random_dp_number(Vin(ii,1))
		end do
		call BPACK_Mult('N',N_unk_loc,1,Vin,Vout1,bmat,ptree,option,stats)
		call matvec_user('N',N_unk_loc,N_unk_loc,1,Vin,Vout2,ker)
		tmp1 = fnorm(Vout2-Vout1,N_unk_loc,1)**2d0
		call MPI_ALLREDUCE(tmp1, norm1, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
		tmp2 = fnorm(Vout2,N_unk_loc,1)**2d0
		call MPI_ALLREDUCE(tmp2, norm2, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
		error = sqrt(norm1)/sqrt(norm2)
		deallocate(Vin,Vout1,Vout2)
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)error,'accuracy of construction'

	else if(explicitflag ==0)then

		!**** register the user-defined function and type in ker
		ker%QuantApp => quant
		ker%FuncHMatVec=>HODLR_MVP_Fullmat

		!**** initialization of the construction phase
	    allocate(Permutation(quant%Nunk))
		call BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree,Coordinates=xyz,tree=tree)
		deallocate(Permutation) ! caller can use this permutation vector if needed
		deallocate(xyz)
		deallocate(tree)


		!**** define other quantities in quant using information returned by BPACK_construction_Init
		call CreateDistDenseMat(quant%Nunk,msh,ptree,quant)


		!**** computation of the construction phase
		call BPACK_construction_Matvec(bmat,matvec_user,Memory,error,option,stats,ker,ptree,msh)


	end if



	!*********** Construct the second HODLR by using the first HODLR as a matvec

	call CopyOptions(option,option1)
	option1%nogeo=1   ! this indicates the second HOLDR construction requires no geometry information
	option1%xyzsort=NATURAL ! this indicates the second HOLDR construction requires no reordering

	!**** register the user-defined function and type in ker
	ker1%FuncZmn=>Zelem_FULL
	ker1%FuncHMatVec=>HODLR_MVP_OneHODLR
	ker1%QuantApp=>quant1

	quant1%bmat=>bmat
	quant1%msh=>msh
	quant1%ptree=>ptree
	quant1%stats=>stats
	quant1%option=>option
	quant1%Nunk=quant%Nunk

	msh1%Nunk = msh%Nunk
	call InitStat(stats1)


	!**** generate the process tree for the second HODLR, can use larger number of MPIs if you want to
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo
	call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree1)
	deallocate(groupmembers)


	!**** use the clustering tree from the first HODLR
	select case(option%format)
	case(HODLR)
		Maxlevel=bmat%ho_bf%Maxlevel
	case(HMAT)
		Maxlevel=bmat%h_mat%Maxlevel
	end select
	allocate (tree(2**Maxlevel))
	do ii=1,2**Maxlevel
		tree(ii)=msh%basis_group(2**Maxlevel+ii-1)%tail-msh%basis_group(2**Maxlevel+ii-1)%head+1
	enddo


	!**** initialization of the construction phase
	allocate(Permutation(quant1%Nunk))
	call BPACK_construction_Init(quant1%Nunk,Permutation,Nunk_loc,bmat1,option1,stats1,msh1,ker1,ptree1,tree=tree)
	deallocate(Permutation) ! caller can use this permutation vector if needed
	deallocate(tree)


	!**** computation of the construction phase
	call BPACK_construction_Matvec(bmat1,matvec_user,Memory,error,option1,stats1,ker1,ptree1,msh1)


	!**** deletion of quantities
	if(allocated(quant%matZ_glo))deallocate(quant%matZ_glo)
	if(allocated(quant%matZ_loc))deallocate(quant%matZ_loc)
	if(associated(quant%N_p))deallocate(quant%N_p)

	call delete_proctree(ptree)
	call delete_Hstat(stats)
	call delete_mesh(msh)
	call delete_kernelquant(ker)
	call BPACK_delete(bmat)

	call delete_proctree(ptree1)
	call delete_Hstat(stats1)
	call delete_mesh(msh1)
	call delete_kernelquant(ker1)
	call BPACK_delete(bmat1)


    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	call blacs_exit(1)
	call MPI_Finalize(ierr)


end PROGRAM ButterflyPACK_ScatteringMatrix_Matvec



