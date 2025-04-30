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
!> @brief This example reads a full matrix of kernel ridge system from disk, and compress it using entry-valuation-based APIs, and evaluate the prediction error.
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
		integer ntrain,ntest,diagflag ! size of training points and test points
		character(LEN=500)::trainfile_p,trainfile_tree,trainfile_l,trainfile_d,testfile_p,testfile_l !Kernel Regression: file pointers to train data, preordered tree, train labels, test data, and test labels
		real(kind=8)::lambda
		integer Nunk ! size of the matrix
		! real(kind=8),allocatable:: xyz(:,:)   ! coordinates of the points
		real(kind=8),allocatable:: matZ_glo(:,:)
		integer,allocatable:: perms(:)
		real(kind=8),allocatable:: diags(:)
	end type quant_app

contains



	!**** user-defined subroutine to sample Z_mn
	subroutine Zelem_FULL(m,n,value_e,quant)
		use d_BPACK_DEFS
		implicit none

		class(*),pointer :: quant
		integer, INTENT(IN):: m,n
		real(kind=8)::value_e

		real(kind=8) r_mn
		integer dimn

		select TYPE(quant)
		type is (quant_app)
			if(m<n)then
			value_e = quant%matZ_glo(m,n)
			else
			value_e = quant%matZ_glo(n,m)
			endif
			if(m==n)value_e = value_e +quant%lambda
		class default
			write(*,*)"unexpected type"
			stop
		end select
	end subroutine Zelem_FULL

end module APPLICATION_MODULE





PROGRAM ButterflyPACK_FullKRR
    use d_BPACK_DEFS
    use APPLICATION_MODULE

	use d_BPACK_structure
	use d_BPACK_factor
	use d_BPACK_constr
#ifdef HAVE_OPENMP
	use omp_lib
#endif
	use d_MISC_Utilities
	use d_BPACK_utilities

    implicit none

	! include "mkl_vml.fi"

    real(kind=8) para
    real(kind=8) tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer times(8)
	real(kind=8) t1,t2,x,y,z,r,theta,phi

	character(len=:),allocatable  :: string
	character(len=1024)  :: strings,strings1
	character(len=6)  :: info_env
	integer :: length
	integer :: edge_m,edge_n,flag,nargs
	integer :: ierr
	integer*8 oldmode,newmode
	type(d_Hoption)::option
	type(d_Hstat)::stats
	type(d_mesh)::msh
	type(d_kernelquant)::ker
	type(quant_app),target::quant
	type(d_Bmatrix)::bmat
	integer,allocatable:: groupmembers(:)
	integer nmpi
	type(d_proctree)::ptree
	CHARACTER (LEN=1000) DATA_DIR
	integer,allocatable::Permutation(:)
	integer Nunk_loc
	integer v_major,v_minor,v_bugfix
	real(kind=8),allocatable:: tmpline(:)

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
    write(*,*) "ButterflyPACK_FullKRR"
    write(*,*) "   "
	endif

   !***********************************************************************
   if(ptree%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'FULLKER computing'
   call d_BPACK_GetVersionNumber(v_major,v_minor,v_bugfix)
   write(*,'(A23,I1,A1,I1,A1,I1,A1)') " ButterflyPACK Version:",v_major,".",v_minor,".",v_bugfix
   write (*,*) ''
   endif
   !***********************************************************************



	!**** initialize stats and option
	call d_InitStat(stats)
	call d_SetDefaultOptions(option)

	!**** set solver parameters
	option%xyzsort=NATURAL
	option%nogeo=1
	option%Nmin_leaf=500
	option%ErrFillFull=0
	option%ErrSol=1


	!**** intialize the user-defined derived type quant
	quant%trainfile_p='../EXAMPLE/FULLMAT_DATA/FullMatKrr/kernel_train5k_test10k.txt'
	quant%trainfile_l='../EXAMPLE/FULLMAT_DATA/FullMatKrr/label_train5k_test10k.txt'
    quant%ntrain=5000
    quant%ntest=10000
    quant%lambda=0

	quant%diagflag=0
	nargs = iargc()
	ii=1
	do while(ii<=nargs)
		call getarg(ii,strings)
		if(trim(strings)=='-quant')then ! user-defined quantity parameters   !**** read data file directory, data set dimension, size of training and testing sets, RBF parameters and solver parameters
			flag=1
			do while(flag==1)
				ii=ii+1
				if(ii<=nargs)then
					call getarg(ii,strings)
					if(strings(1:2)=='--')then
						ii=ii+1
						call getarg(ii,strings1)
						if(trim(strings)=='--trainfile_p')then
							quant%trainfile_p=trim(strings1)
						else if	(trim(strings)=='--trainfile_l')then
							quant%trainfile_l=trim(strings1)
						else if	(trim(strings)=='--trainfile_d')then
							quant%trainfile_d=trim(strings1)
							quant%diagflag=1
						else if	(trim(strings)=='--ntrain')then
							read(strings1,*)quant%ntrain
						else if	(trim(strings)=='--ntest')then
							read(strings1,*)quant%ntest
						else if	(trim(strings)=='--lambda')then
							read(strings1,*)quant%lambda
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

	quant%Nunk = quant%ntrain

	call d_PrintOptions(option,ptree)


	t1 = MPI_Wtime()
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "reading fullmatrix......"
	allocate(quant%perms(quant%ntrain+quant%ntest))
	do ii=1,quant%ntrain+quant%ntest
		quant%perms(ii)=ii
	enddo
	! call d_rperm(quant%ntrain+quant%ntest,quant%perms)

	allocate(tmpline(quant%ntrain))
	open (90,file=quant%trainfile_p)
	allocate (quant%matZ_glo(quant%ntrain+quant%ntest,quant%ntrain))
	do edge_m=1,quant%ntrain+quant%ntest
	read (90,*) tmpline
	do edge_n=1,quant%ntrain
		quant%matZ_glo(quant%perms(edge_m),quant%perms(edge_n)) = tmpline(edge_n)
	enddo
	enddo
	close(90)
	deallocate(tmpline)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "reading fullmatrix finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = MPI_Wtime()

	!**** register the user-defined function and type in ker
	ker%QuantApp => quant
	ker%FuncZmn => Zelem_FULL

	!**** initialization of the construction phase
	allocate(Permutation(quant%Nunk))
	call d_BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree)
	deallocate(Permutation) ! caller can use this permutation vector if needed


	!**** computation of the construction phase
    call d_BPACK_construction_Element(bmat,option,stats,msh,ker,ptree)


	!**** factorization phase
    call d_BPACK_Factorization(bmat,option,stats,ptree,msh)


	!**** solve phase
    call FULLKER_solve(bmat,option,msh,quant,ptree,stats)



	!**** deletion of quantities
	call d_delete_proctree(ptree)
	call d_delete_Hstat(stats)
	call d_delete_mesh(msh)
	call d_delete_kernelquant(ker)
	call d_BPACK_delete(bmat)

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	call d_blacs_exit_wrp(1)
	call MPI_Finalize(ierr)

end PROGRAM ButterflyPACK_FullKRR



subroutine FULLKER_solve(bmat,option,msh,quant,ptree,stats)

    use d_BPACK_DEFS
	use APPLICATION_MODULE
#ifdef HAVE_OPENMP
	use omp_lib
#endif
	use d_BPACK_Solve_Mul
	use d_MISC_Utilities

    implicit none

    integer i, j, ii, jj, iii, jjj, ierr, ntest,Dimn,edge_m,edge_n,ncorrect
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, n_iter_max, iter, N_unk, N_unk_loc
    real(kind=8) theta, phi, dphi, rcs_V, rcs_H,error
    real(kind=8) T0,T1
    real(kind=8) n1,n2,rtemp
    real(kind=8) value_Z
    real(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:),vout(:,:),vout_tmp(:,:),vout_test(:,:),matrixtemp1(:,:),matrixtemp2(:,:),matrixtemp(:,:)
	real(kind=8):: rel_error
	type(d_Hoption)::option
	type(d_mesh)::msh
	type(quant_app)::quant
	type(d_proctree)::ptree
	type(d_Bmatrix)::bmat
	type(d_Hstat)::stats
	real(kind=8),allocatable:: current(:),voltage(:)
	real(kind=8), allocatable:: labels(:)
	integer, allocatable :: ipiv(:)
	real(kind=8) r_mn
	real(kind=8) label


    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Solve and Prediction......"

	N_unk=msh%Nunk
	N_unk_loc = msh%idxe-msh%idxs+1
	ntest=quant%ntest

	if(option%ErrSol==1)then
		call d_BPACK_Test_Solve_error(bmat,N_unk_loc,option,ptree,stats)
	endif


	!**** read training label as local RHS and solve for the weights


	allocate(labels(quant%ntrain+quant%ntest))
	allocate (x(N_unk_loc,1))
	x=0
	allocate (b(N_unk_loc,1))
	open (91,file=quant%trainfile_l)
	! read(91,*)N_unk,Dimn
	do ii=1,quant%ntrain+quant%ntest
		read(91,*)labels(quant%perms(ii))
	enddo
	close(91)
	do ii=1,N_unk_loc
		b(ii,1) = labels(msh%new2old(ii-1+msh%idxs))
	enddo
	allocate (vout_test(ntest,1))
	vout_test(:,1) = labels(1+quant%ntrain:quant%ntrain+quant%ntest)
	deallocate(labels)


	n1 = MPI_Wtime()

	call d_BPACK_Solution(bmat,x,b,N_unk_loc,1,option,ptree,stats)

	n2 = MPI_Wtime()
	stats%Time_Sol = stats%Time_Sol + n2-n1
	call MPI_ALLREDUCE(stats%Time_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'Solving:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Flop_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,'(A13,Es14.2)') 'Solve flops:',rtemp

	!**** prediction on the test sets

	if(quant%diagflag==1)then
		allocate(quant%diags(quant%ntrain))
		open (92,file=quant%trainfile_d)
		! read(91,*)N_unk,Dimn
		do ii=1,quant%ntrain
			read(92,*)quant%diags(quant%perms(ii))
		enddo
		close(92)
		do edge=1, N_unk_loc
			x(edge,1) = x(edge,1)*quant%diags(msh%new2old(edge+msh%idxs-1))
		enddo
		deallocate(quant%diags)
	endif



	T0 = MPI_Wtime()

	allocate (vout(ntest,1))
	allocate (vout_tmp(ntest,1))
	vout_tmp = 0
	do edge=1, N_unk_loc
		do edge_m=1,ntest
			value_Z = quant%matZ_glo(edge_m+quant%ntrain,msh%new2old(edge+msh%idxs-1))
			vout_tmp(edge_m,1) = vout_tmp(edge_m,1) + value_Z*x(edge,1)
		enddo
	enddo

	call MPI_REDUCE(vout_tmp, vout, ntest,MPI_double_precision, MPI_SUM, Main_ID, ptree%Comm,ierr)

		! allocate(matrixtemp2(N_unk,N_unk))
		! allocate(matrixtemp1(N_unk,N_unk))
		! matrixtemp2=quant%matZ_glo(msh%new2old,msh%new2old)
		! call d_GeneralInverse(N_unk,N_unk,matrixtemp2,matrixtemp1,1d-10)

		! allocate(ipiv(N_unk))
		! allocate(matrixtemp1(N_unk,N_unk))
		! matrixtemp1=quant%matZ_glo(msh%new2old,msh%new2old)
		! ipiv=0
		! call d_getrff90(matrixtemp1,ipiv)
		! call d_getrif90(matrixtemp1,ipiv)
		! deallocate(ipiv)




		! allocate(matrixtemp(ntest,N_unk))
		! matrixtemp=quant%matZ_glo(quant%ntrain+1:quant%ntrain+quant%ntest,msh%new2old)
		! call d_gemmf90(matrixtemp1,N_unk,b,N_unk,x,N_unk,'N','N',N_unk,1,N_unk,BPACK_cone,BPACK_czero)
		! allocate (vout(ntest,1))
		! vout=0
		! call d_gemmf90(matrixtemp,ntest,x,N_unk,vout,ntest,'N','N',ntest,1,N_unk,BPACK_cone,BPACK_czero)

	T1 = MPI_Wtime()
	if (ptree%MyID==Main_ID) then
		vout_tmp = vout-vout_test
		error=0
		do ii=1,ntest
			error = error + vout_tmp(ii,1)**2d0
		enddo
		error = error/ntest
		error = sqrt(error)

		write (*,*) ''
		write (*,*) 'Prediction time:',T1-T0,'Seconds'
		write (*,*) 'Prediction error:',error
		write (*,*) ''

	endif

	deallocate (vout)
	deallocate (vout_tmp)
	deallocate (vout_test)

	deallocate(x)
	deallocate(b)

	if(allocated(quant%matZ_glo))deallocate(quant%matZ_glo)
	if(allocated(quant%perms))deallocate(quant%perms)

	call MPI_barrier(ptree%Comm,ierr)

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Solve and Prediction finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "

    return

end subroutine FULLKER_solve


