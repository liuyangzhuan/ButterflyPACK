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
!> @brief This example generates a RBF kernel using training and testing data from disk, compress it using entry-valuation-based APIs, and evaluate the prediction error.
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
		real(kind=8) sigma, lambda ! Kernel Regression: RBF parameters
		integer dimn  ! dimension of the data sets
		integer ntrain,ntest ! size of training points and test points
		character(LEN=500)::trainfile_p,trainfile_tree,trainfile_l,testfile_p,testfile_l !Kernel Regression: file pointers to train data, preordered tree, train labels, test data, and test labels
		CHARACTER (LEN=1000) DATA_DIR
		integer Nunk ! size of the matrix
		real(kind=8),allocatable:: xyz(:,:)   ! coordinates of the points

	end type quant_app

contains


	!**** cutoff distance for gaussian kernel
	real(kind=8) function arg_thresh_Zmn(quant)
		use d_BPACK_DEFS
		implicit none

		type(quant_app)::quant
		arg_thresh_Zmn=-log(BPACK_SafeUnderflow)*2.0*quant%sigma**2
	end function arg_thresh_Zmn



	!**** user-defined subroutine to sample Z_mn
	subroutine Zelem_RBF(m,n,value_e,quant)
		use d_BPACK_DEFS
		implicit none

		class(*),pointer :: quant
		integer, INTENT(IN):: m,n
		real(kind=8)::value_e

		real(kind=8) r_mn
		integer dimn

		select TYPE(quant)
		type is (quant_app)
			dimn = size(quant%xyz,1)
			value_e=0
			r_mn=sum((quant%xyz(1:dimn,m)-quant%xyz(1:dimn,n))**2)
			if(r_mn<arg_thresh_Zmn(quant))then
				value_e = exp(-r_mn/2.0/quant%sigma**2)
			endif
			if(r_mn==0)value_e = value_e + quant%lambda
		class default
			write(*,*)"unexpected type"
			stop
		end select
	end subroutine Zelem_RBF

end module APPLICATION_MODULE





PROGRAM ButterflyPACK_KRR
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
    integer i,j,k, threads_num,iii
	integer times(8)
	real(kind=8) t1,t2,x,y,z,r,theta,phi

	character(len=:),allocatable  :: string
	character(len=1024)  :: strings,strings1
	character(len=6)  :: info_env
	integer :: length
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
	integer MPI_thread
	integer,allocatable::Permutation(:)
	integer Nunk_loc
	integer nargs,flag
	integer v_major,v_minor,v_bugfix

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
    write(*,*) "ButterflyPACK_KRR"
	call d_BPACK_GetVersionNumber(v_major,v_minor,v_bugfix)
	write(*,'(A23,I1,A1,I1,A1,I1,A1)') " ButterflyPACK Version:",v_major,".",v_minor,".",v_bugfix
    write(*,*) "   "
	endif

	!**** initialize stats and option
	call d_InitStat(stats)
	call d_SetDefaultOptions(option)


	!**** intialize the user-defined derived type quant
	quant%DATA_DIR='../EXAMPLE/KRR_DATA/susy_10Kn'
    quant%dimn=8
    quant%ntrain=10000
    quant%ntest=1000
	quant%sigma=0.1
	quant%lambda=10.0


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
						if(trim(strings)=='--dimn')then
							read(strings1,*)quant%dimn
						else if	(trim(strings)=='--data_dir')then
							quant%data_dir=trim(strings1)
						else if	(trim(strings)=='--ntrain')then
							read(strings1,*)quant%ntrain
						else if	(trim(strings)=='--ntest')then
							read(strings1,*)quant%ntest
						else if	(trim(strings)=='--sigma')then
							read(strings1,*)quant%sigma
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

	call d_PrintOptions(option,ptree)


	quant%trainfile_p=trim(quant%DATA_DIR)//'_train.csv'
	quant%trainfile_l=trim(quant%DATA_DIR)//'_train_label.csv'
	quant%testfile_p=trim(quant%DATA_DIR)//'_test.csv'
	quant%testfile_l=trim(quant%DATA_DIR)//'_test_label.csv'
	quant%Nunk = quant%ntrain
	if(ptree%MyID==Main_ID)write(*,*)'training set: ',trim(quant%trainfile_p)

	t1 = MPI_Wtime()
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "geometry modeling......"
	open (90,file=quant%trainfile_p)
	allocate (quant%xyz(quant%dimn,1:quant%Nunk))
	do ii=1,quant%Nunk
		read (90,*) quant%xyz(1:quant%dimn,ii)
	enddo
    close(90)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "modeling finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = MPI_Wtime()


	!**** register the user-defined function and type in ker
	ker%QuantApp => quant
	ker%FuncZmn => Zelem_RBF

    !**** initialization of the construction phase
    allocate(Permutation(quant%Nunk))
	call d_BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree,Coordinates=quant%xyz)
	deallocate(Permutation) ! caller can use this permutation vector if needed


	!**** computation of the construction phase
    call d_BPACK_construction_Element(bmat,option,stats,msh,ker,ptree)

	!**** factorization phase
    call d_BPACK_Factorization(bmat,option,stats,ptree,msh)


	!**** solve phase
    call RBF_solve(bmat,option,msh,quant,ptree,stats)


	!**** deletion of quantities
	if(allocated(quant%xyz))deallocate(quant%xyz)
	call d_delete_proctree(ptree)
	call d_delete_Hstat(stats)
	call d_delete_mesh(msh)
	call d_delete_kernelquant(ker)
	call d_BPACK_delete(bmat)

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"


	call blacs_exit(1)
	call MPI_Finalize(ierr)

end PROGRAM ButterflyPACK_KRR


subroutine RBF_solve(bmat,option,msh,quant,ptree,stats)

    use d_BPACK_DEFS
	use APPLICATION_MODULE
#ifdef HAVE_OPENMP
	use omp_lib
#endif
	use d_BPACK_Solve_Mul

    implicit none

    integer i, j, ii, jj, iii, jjj, ierr, ntest,Dimn,edge_m,edge_n,ncorrect
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, n_iter_max, iter, N_unk, N_unk_loc
    real(kind=8) theta, phi, dphi, rcs_V, rcs_H,rate
    real(kind=8) T0,T1
    real(kind=8) n1,n2,rtemp
    real(kind=8) value_Z
    real(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:),vout(:,:),vout_tmp(:,:)
	real(kind=8):: rel_error
	type(d_Hoption)::option
	type(d_mesh)::msh
	type(quant_app)::quant
	type(d_proctree)::ptree
	type(d_Bmatrix)::bmat
	type(d_Hstat)::stats
	real(kind=8),allocatable:: current(:),voltage(:)
	real(kind=8), allocatable:: labels(:)
	real(kind=8),allocatable:: xyz_test(:,:)
	real(kind=8) r_mn
	integer label

	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Solve and Prediction......"

	N_unk=msh%Nunk
	Dimn=quant%dimn
	N_unk_loc = msh%idxe-msh%idxs+1


	if(option%ErrSol==1)then
		call d_BPACK_Test_Solve_error(bmat,N_unk_loc,option,ptree,stats)
	endif


	!**** read training label as local RHS and solve for the weights



	allocate(labels(N_unk))
	allocate (x(N_unk_loc,1))
	x=0
	allocate (b(N_unk_loc,1))
	open (91,file=quant%trainfile_l)
	! read(91,*)N_unk,Dimn
	do ii=1,N_unk
		read(91,*)label
		labels(ii)=label
	enddo
	do ii=1,N_unk_loc
		b(ii,1) = labels(msh%new2old(ii-1+msh%idxs))
	enddo
	deallocate(labels)
	close(91)

	n1 = MPI_Wtime()

	call d_BPACK_Solution(bmat,x,b,N_unk_loc,1,option,ptree,stats)

	n2 = MPI_Wtime()
	stats%Time_Sol = stats%Time_Sol + n2-n1
	call MPI_ALLREDUCE(stats%Time_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'Solving:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Flop_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,'(A13Es14.2)') 'Solve flops:',rtemp

	!**** prediction on the test sets

	ntest=quant%ntest
	T0 = MPI_Wtime()
	open (92,file=quant%testfile_p)
	allocate (xyz_test(Dimn,ntest))
	do edge=1,ntest
		read (92,*) xyz_test(1:Dimn,edge)
	enddo
	close(92)

	allocate (vout(ntest,1))
	allocate (vout_tmp(ntest,1))
	vout_tmp = 0
	do edge=1, N_unk_loc
		do edge_m=1,ntest
			r_mn=sum((xyz_test(1:dimn,edge_m)-quant%xyz(1:dimn,msh%new2old(edge+msh%idxs-1)))**2)
			value_Z = exp(-r_mn/2.0/quant%sigma**2)
			vout_tmp(edge_m,1) = vout_tmp(edge_m,1) + value_Z*x(edge,1)
		enddo
	enddo

	call MPI_REDUCE(vout_tmp, vout, ntest,MPI_double_precision, MPI_SUM, Main_ID, ptree%Comm,ierr)
	T1 = MPI_Wtime()
	if (ptree%MyID==Main_ID) then
		do ii=1,ntest
			if(dble(vout(ii,1))>0)then
				vout(ii,1)=1
			else
				vout(ii,1)=-1
			endif
		enddo


		open (93,file=quant%testfile_l)
		do edge=1,ntest
			read (93,*) label
			vout_tmp(edge,1)=label
		enddo
		close(93)

		ncorrect=0
		do edge=1,ntest
			if(dble(vout_tmp(edge,1))*dble(vout(edge,1))>0)then
				ncorrect = ncorrect + 1
			endif
		enddo

		rate = dble(ncorrect)/dble(ntest)

		write (*,*) ''
		write (*,*) 'Prediction time:',T1-T0,'Seconds'
		write (*,*) 'Success rate:',rate
		write (*,*) ''

	endif

	deallocate (vout)
	deallocate (vout_tmp)

	deallocate(x)
	deallocate(b)
	deallocate(xyz_test)

	call MPI_barrier(ptree%Comm,ierr)

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Solve and Prediction finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "

    return

end subroutine RBF_solve


