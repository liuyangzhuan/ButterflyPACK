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
!> @brief This example reads a frontal matrix from distributed files from disk, and compress it using entry-valuation or matvec-based APIs
!> @details Note that instead of the use of precision dependent subroutine/module/type names "z_", one can also use the following \n
!> #define DAT 0 \n
!> #include "zButterflyPACK_config.fi" \n
!> which will macro replace precision-independent subroutine/module/type names "X" with "z_X" defined in SRC_DOUBLECOMLEX with double-complex precision


! This exmple works with double-complex precision data
module APPLICATION_MODULE
use z_BPACK_DEFS
implicit none

	!**** define your application-related variables here
	type quant_app
		integer::Nunk ! matrix size
		complex(kind=8), allocatable :: matZ_loc(:,:) ! Local Matrix: Loccal matrix in a npx1 blasc z_grid
		integer,pointer :: N_p(:,:)=>null() ! column sizes of all processes sharing this hodlr
		type(z_Bmatrix),pointer::bmat ! Use this metadata in matvec
		type(z_mesh),pointer::msh   ! Use this metadata in matvec
		type(z_proctree),pointer::ptree ! Use this metadata in matvec
		type(z_Hstat),pointer::stats ! Use this metadata in matvec
		type(z_Hoption),pointer::option ! Use this metadata in matvec
		CHARACTER (LEN=1000) DATA_DIR
		integer:: explicitflag
	end type quant_app

contains




	subroutine HODLR_MVP_OneHODLR(trans,Mloc,Nloc,num_vect,Vin,Vout,quant)
		use z_BPACK_DEFS
		use z_MISC_DenseLA
		use z_MISC_Utilities
		use z_BPACK_Solve_Mul
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
		type(z_Bmatrix),pointer::bmat

		select TYPE(quant)
		type is (z_quant_bmat)
			pgno=1
			nproc = quant%ptree%pgrp(pgno)%nproc

			bmat=>quant%bmat
			call z_BPACK_Mult(trans,Nloc,num_vect,Vin,Vout,bmat,quant%ptree,quant%option,quant%stats)
		class default
			write(*,*)"unexpected quant type in HODLR_MVP_OneHODLR"
			stop
		end select

	end subroutine HODLR_MVP_OneHODLR

	subroutine HODLR_MVP_Fullmat(trans,Mloc,Nloc,num_vect,Vin,Vout,quant)
		use z_BPACK_DEFS
		use z_MISC_DenseLA
		use z_MISC_Utilities
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


		!!!!**** generate 2D z_grid blacs quantities
		ctxt = quant%ptree%pgrp(pgno)%ctxt
		call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
		if(myrow/=-1 .and. mycol/=-1)then
			myArows = z_numroc_wp(N, nbslpk, myrow, 0, nprow)
			myAcols = z_numroc_wp(num_vect, nbslpk, mycol, 0, npcol)
			call descinit( descsVin2D, N, num_vect, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
			call descinit( descsVout2D, N, num_vect, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
			myArows = z_numroc_wp(N, nbslpk, myrow, 0, nprow)
			myAcols = z_numroc_wp(N, nbslpk, mycol, 0, npcol)
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
		call z_Redistribute1Dto2D(Vin,quant%N_p,0,pgno,Vin_tmp_2D,N,0,pgno,num_vect,quant%ptree)


		!!!!**** perform gemm on 2d z_grid
		if(myrow/=-1 .and. mycol/=-1)then
			call z_pgemmf90(trans,'N',N,num_vect,N,BPACK_cone, quant%matZ_loc,1,1,descsMat2D,Vin_tmp_2D,1,1,descsVin2D,BPACK_czero,Vout_tmp_2D,1,1,descsVout2D)
		endif


		!!!!**** redistribution of output vectors
		call z_Redistribute2Dto1D(Vout_tmp_2D,N,0,pgno,Vout,quant%N_p,0,pgno,num_vect,quant%ptree)


		!!!!**** deallocation buffers
		if(myrow/=-1 .and. mycol/=-1)then
			deallocate(Vin_tmp_2D)
			deallocate(Vout_tmp_2D)
		endif


		end select

	end subroutine HODLR_MVP_Fullmat


	subroutine CreateDistDenseMat(N,msh,ptree,quant,option,DATA_DIR)
		use z_BPACK_DEFS
		use z_MISC_DenseLA
		use z_MISC_Utilities
		implicit none
		complex(kind=8),allocatable:: Vin_tmp(:,:),Vout_tmp(:,:)
		complex(kind=8) ctemp,a,b
		integer ii,jj,nn,fl_transpose,kk,black_step
		integer, INTENT(in)::N
		real(kind=8) n1,n2,t1,t2
		real(kind=8),allocatable::tmp(:)
		type(z_mesh)::msh
		type(z_proctree)::ptree
		type(z_Hoption)::option
		type(quant_app) :: quant
		integer nproc,ctxt,nb1Dc, nb1Dr, level_p,pgno,num_blocks,ii_new,gg,proc,myi,myj,myAcols,myArows,nprow,npcol,myrow,mycol
		CHARACTER (*) DATA_DIR
		character(len=1024)  :: smyrow,smycol


		pgno=1
		nproc = ptree%pgrp(pgno)%nproc

		!!!!****** allocate index array for 1D HODLR layout
		level_p = ptree%nlevel-z_GetTreelevel(pgno)
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


		!!!!****** assemble the full matrix in 2D blacs layout
		ctxt = ptree%pgrp(pgno)%ctxt
		call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)



		if(myrow/=-1 .and. mycol/=-1)then
		myArows = z_numroc_wp(N, nbslpk, myrow, 0, nprow)
		myAcols = z_numroc_wp(N, nbslpk, mycol, 0, npcol)
		allocate(quant%matZ_loc(myArows,myAcols))



		t1 = MPI_Wtime()
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Generating fullmat ......"
		write(smyrow , *) myrow
		write(smycol , *) mycol
		open(unit=888,file=trim(DATA_DIR)//'_pr'//trim(adjustl(smyrow))//'_pc'//trim(adjustl(smycol))//'.dat',status='old')
		read(888,*) ! skip the first line
		read(888,*) ! skip the second line !! M1,N1,Mloc1,Nloc1,Mbslpk1,Nbslpk1,nprow1,npcol1,myrow1,mycol1

		allocate(tmp(1:myAcols*2))
		do ii=1,myArows
		read(888,*)tmp
		do jj=1,myAcols
			quant%matZ_loc(ii,jj)=tmp(jj*2-1)+BPACK_junit*tmp(jj*2)
		enddo
		end do
		deallocate(tmp)
		close(888)


		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Generating fullmat finished"
		t2 = MPI_Wtime()
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)t2-t1, 'secnds'


		endif

	end subroutine CreateDistDenseMat




end module APPLICATION_MODULE


PROGRAM ButterflyPACK_FrontalMatrix_Matvec

    use z_BPACK_DEFS

	use z_BPACK_structure
	use z_BPACK_factor
	use z_BPACK_constr
#ifdef HAVE_OPENMP
	use omp_lib
#endif
	use z_Bplus_compress
	use z_BPACK_randomMVP
	use z_BPACK_utilities
	use APPLICATION_MODULE

    implicit none

    real(kind=8) para,error,tmp1,tmp2,norm1,norm2
    real(kind=8) tolerance,rankrate
    integer Primary_block, nn, mm
    integer i,j,k, threads_num,ii,jj
	integer seed_myid(50)
	integer times(8)
	real(kind=8) t1,t2,t3,t4,x,y,z,r,theta,phi,Memory
	real(kind=8),allocatable::tmp(:)
	complex(kind=8),allocatable:: InputVec(:)
	complex(kind=8):: ctemp
	integer kk,black_step,rank0
	complex(kind=8),allocatable::Vout1(:,:),Vout2(:,:),Vin(:,:)
	character(len=1024)  :: strings,strings1
	type(z_Hoption),target:: option,option1
	type(z_Hstat),target::stats,stats1
	type(z_mesh),target::msh,msh1
	type(z_kernelquant),target::ker,ker1
	type(z_Bmatrix),target::bmat,bmat1
	integer Nin1,Nout1,Nin2,Nout2
	type(z_proctree),target::ptree,ptree1
	integer,allocatable:: groupmembers(:)
	integer :: ierr
	integer :: nmpi
	type(quant_app),target::quant
	type(z_quant_bmat),target::quant1
	integer N_unk_loc,Maxlevel
	integer,allocatable::tree(:),Permutation(:)
	real(kind=8),allocatable::xyz(:,:)
	integer Nunk_loc
	integer nprow, npcol, myrow, mycol, nprow1, npcol1, myrow1, mycol1, M1,N1,Mloc1,Nloc1,Mbslpk1,Nbslpk1
	character(len=1024)  :: smyrow,smycol
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
	call z_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	deallocate(groupmembers)

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------Program Start----------------------------------"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "ButterflyPACK_FrontalMatrix_Matvec"
	call z_BPACK_GetVersionNumber(v_major,v_minor,v_bugfix)
	if(ptree%MyID==Main_ID)write(*,'(A23,I1,A1,I1,A1,I1,A1)') " ButterflyPACK Version:",v_major,".",v_minor,".",v_bugfix
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "   "

	!**** initialize stats and option
	call z_InitStat(stats)
	call z_SetDefaultOptions(option)

	!**** intialize the user-defined derived type quant
	quant%ptree=>ptree


	option%nogeo=1   ! this indicates the first HOLDR construction requires no geometry information
	option%xyzsort=NATURAL ! this indicates the first HOLDR construction requires no reordering
	option%Nmin_leaf=256
	option%rank0=16
	option%rankrate=2d0
	option%verbosity=2


	quant%explicitflag=0
	!**** initialize the user-defined derived type quant
	!*********** Construct the first HODLR by read-in the full matrix and (if explicitflag=0) use it as a matvec or (if explicitflag=1) use it as a fast entry evaluation
	quant%Nunk = 2500
	quant%DATA_DIR='../EXAMPLE/FULLMAT_DATA/Frontal_elastic/Frontal_elastic_2500'


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
						if	(trim(strings)=='--nunk')then
							read(strings1,*)quant%Nunk
						else if	(trim(strings)=='--data_dir')then
							quant%data_dir=trim(strings1)
						else if	(trim(strings)=='--explicitflag')then
							read(strings1,*)quant%explicitflag
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
			call z_ReadOption(option,ptree,ii)
		else
			if(ptree%MyID==Main_ID)write(*,*)'ignoring unknown argument: ',trim(strings)
			ii=ii+1
		endif
	enddo



	call blacs_gridinfo(ptree%pgrp(1)%ctxt, nprow, npcol, myrow, mycol)
	write(smyrow , *) myrow
	write(smycol , *) mycol
	open(unit=888,file=trim(quant%DATA_DIR)//'_pr'//trim(adjustl(smyrow))//'_pc'//trim(adjustl(smycol))//'.dat',status='old')
	read(888,*) ! skip the first line
	read(888,*)M1,N1,Mloc1,Nloc1,Mbslpk1,Nbslpk1,nprow1,npcol1,myrow1,mycol1
	call z_assert(M1==N1 .and. quant%Nunk == M1, 'input matrix dims not matching')
	call z_assert(Mbslpk1==Nbslpk1 .and. Mbslpk1==nbslpk, 'block sizes not matching')
	call z_assert(nprow1==nprow .and. npcol1==npcol, 'process grid dims not matching')
	close(888)


	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'Blackbox HODLR for frontal matrix compression'


	! allocate(tree(16))
	! tree = (/3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 4/)

	allocate(tree(1))
	tree(1) = quant%Nunk






	!**** register the user-defined function and type in ker
	ker%QuantApp => quant
	ker%FuncHMatVec=>HODLR_MVP_Fullmat

	!**** initialization of the construction phase
	allocate(Permutation(quant%Nunk))
	call z_BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree,tree=tree)
	deallocate(Permutation) ! caller can use this permutation vector if needed

	!**** define other quantities in quant using information returned by z_BPACK_construction_Init
	call CreateDistDenseMat(quant%Nunk,msh,ptree,quant,option,quant%DATA_DIR)


	!**** computation of the construction phase
	option%less_adapt=0
	call z_BPACK_construction_Matvec(bmat,z_matvec_user,Memory,error,option,stats,ker,ptree,msh)




	deallocate(tree)


	! call z_BPACK_Factorization(bmat,option,stats,ptree,msh)

	call z_PrintStat(stats,ptree)

	!*********** Construct the second HODLR by using the first HODLR as a matvec or entry extraction

	call z_CopyOptions(option,option1)
	option1%nogeo=1   ! this indicates the second HOLDR construction requires no geometry information
	option1%xyzsort=NATURAL ! this indicates the second HOLDR construction requires no reordering
	option1%elem_extract=1 ! this indicates the second HOLDR construction uses user-defined parallel element extraction


	!**** register the user-defined function and type in ker
	! ker1%FuncZmn=>Zelem_FULL
	ker1%FuncHMatVec=>HODLR_MVP_OneHODLR
	ker1%QuantApp=>quant1
	ker1%FuncZmnBlock=>z_Zelem_block_Extraction

	quant1%bmat=>bmat
	quant1%msh=>msh
	quant1%ptree=>ptree
	quant1%stats=>stats
	quant1%option=>option
	! quant1%Nunk=quant%Nunk



	msh1%Nunk = msh%Nunk
	call z_InitStat(stats1)


	!**** generate the process tree for the second HODLR, can use larger number of MPIs if you want to
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo
	call z_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree1)
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
	allocate(Permutation(msh1%Nunk))
	call z_PrintOptions(option1,ptree1)
	call z_BPACK_construction_Init(msh1%Nunk,Permutation,Nunk_loc,bmat1,option1,stats1,msh1,ker1,ptree1,tree=tree)
	deallocate(Permutation) ! caller can use this permutation vector if needed
	deallocate(tree)



	if(quant%explicitflag ==1)then
		!**** computation of the construction phase
		call z_BPACK_construction_Element(bmat1,option1,stats1,msh1,ker1,ptree1)
	else if(quant%explicitflag ==0)then
		!**** computation of the construction phase
		option1%less_adapt=0
		call z_BPACK_construction_Matvec(bmat1,z_matvec_user,Memory,error,option1,stats1,ker1,ptree1,msh1)
	endif


	call z_PrintStat(stats1,ptree)

	!**** deletion of quantities
	if(allocated(quant%matZ_loc))deallocate(quant%matZ_loc)
	if(associated(quant%N_p))deallocate(quant%N_p)

	call z_delete_proctree(ptree)
	call z_delete_Hstat(stats)
	call z_delete_mesh(msh)
	call z_delete_kernelquant(ker)
	call z_BPACK_delete(bmat)

	call z_delete_proctree(ptree1)
	call z_delete_Hstat(stats1)
	call z_delete_mesh(msh1)
	call z_delete_kernelquant(ker1)
	call z_BPACK_delete(bmat1)


    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	call blacs_exit(1)
	call MPI_Finalize(ierr)


end PROGRAM ButterflyPACK_FrontalMatrix_Matvec




