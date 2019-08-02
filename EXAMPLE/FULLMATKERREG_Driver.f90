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
		integer ntrain,ntest ! size of training points and test points
		character(LEN=500)::trainfile_p,trainfile_tree,trainfile_l,testfile_p,testfile_l !Kernel Regression: file pointers to train data, preordered tree, train labels, test data, and test labels

		integer Nunk ! size of the matrix
		! real(kind=8),allocatable:: xyz(:,:)   ! coordinates of the points
		real(kind=8),allocatable:: matZ_glo(:,:)
		integer,allocatable:: perms(:)
	end type quant_app

contains



	!**** user-defined subroutine to sample Z_mn
	subroutine Zelem_FULL(m,n,value_e,quant)
		use BPACK_DEFS
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
			! if(m==n)value_e = value_e + 1d-4
		class default
			write(*,*)"unexpected type"
			stop
		end select
	end subroutine Zelem_FULL

end module APPLICATION_MODULE





PROGRAM ButterflyPACK_FullKRR
    use BPACK_DEFS
    use APPLICATION_MODULE

	use BPACK_structure
	use BPACK_factor
	use BPACK_constr
	use omp_lib
	use MISC_Utilities

    implicit none

	! include "mkl_vml.fi"

    real(kind=8) para
    real(kind=8) tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer times(8)
	real(kind=8) t1,t2,x,y,z,r,theta,phi

	character(len=:),allocatable  :: string
	character(len=1024)  :: strings
	character(len=6)  :: info_env
	integer :: length
	integer :: edge_m,edge_n
	integer :: ierr
	integer*8 oldmode,newmode
	type(Hoption)::option
	type(Hstat)::stats
	type(mesh)::msh
	type(kernelquant)::ker
	type(quant_app),target::quant
	type(Bmatrix)::bmat
	integer,allocatable:: groupmembers(:)
	integer nmpi
	type(proctree)::ptree
	CHARACTER (LEN=1000) DATA_DIR
	integer,allocatable::Permutation(:)
	integer Nunk_loc
	integer v_major,v_minor,v_bugfix

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
    write(*,*) "ButterflyPACK_FullKRR"
    write(*,*) "   "
	endif

	!**** initialize stats and option
	call InitStat(stats)
	call SetDefaultOptions(option)


    !**** read data file directory, data set dimension, size of training and testing sets, RBF parameters
	CALL getarg(1, DATA_DIR)
	quant%trainfile_p=trim(DATA_DIR)//'/QM7-gramian-Zm-1.00-Kv-0.50-Ke-0.05-q-5.0e-02.txt'
	quant%trainfile_l=trim(DATA_DIR)//'/QM7-atomization-energy.txt'
	call getarg(2,strings)
	read(strings,*)quant%ntrain
	quant%Nunk = quant%ntrain
	call getarg(3,strings)
	read(strings,*)quant%ntest
	call getarg(4,strings)
	read(strings,*)option%RecLR_leaf
	call getarg(5,strings)
	read(strings,*)option%tol_comp
	option%tol_rand=option%tol_comp
	option%tol_Rdetect=option%tol_comp*1d-1

    !**** set solver parameters
	option%xyzsort=NATURAL
	option%nogeo=1
	option%Nmin_leaf=500
	option%ErrFillFull=1
	option%ErrSol=1


   !***********************************************************************
   if(ptree%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'FULLKER computing'
   call BPACK_GetVersionNumber(v_major,v_minor,v_bugfix)
   write(*,'(A23,I1,A1,I1,A1,I1,A1)') " ButterflyPACK Version:",v_major,".",v_minor,".",v_bugfix
   write (*,*) ''
   endif
   !***********************************************************************

	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "reading fullmatrix......"
	allocate(quant%perms(quant%ntrain+quant%ntest))
	do ii=1,quant%ntrain+quant%ntest
		quant%perms(ii)=ii
	enddo
	! call rperm(quant%ntrain+quant%ntest,quant%perms)
	open (90,file=quant%trainfile_p)
	allocate (quant%matZ_glo(quant%ntrain+quant%ntest,quant%ntrain+quant%ntest))
	do edge_m=1,quant%ntrain+quant%ntest
	do edge_n=1,quant%ntrain+quant%ntest
		read (90,*) quant%matZ_glo(quant%perms(edge_m),quant%perms(edge_n))
	enddo
	enddo
	close(90)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "reading fullmatrix finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()

	!**** register the user-defined function and type in ker
	ker%QuantApp => quant
	ker%FuncZmn => Zelem_FULL

	!**** initialization of the construction phase
	allocate(Permutation(quant%Nunk))
	call BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree)
	deallocate(Permutation) ! caller can use this permutation vector if needed


	!**** computation of the construction phase
    call BPACK_construction_Element(bmat,option,stats,msh,ker,ptree)


	!**** factorization phase
    call BPACK_Factorization(bmat,option,stats,ptree,msh)


	!**** solve phase
    call FULLKER_solve(bmat,option,msh,quant,ptree,stats)



	!**** deletion of quantities
	call delete_proctree(ptree)
	call delete_Hstat(stats)
	call delete_mesh(msh)
	call delete_kernelquant(ker)
	call BPACK_delete(bmat)

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	call blacs_exit(1)
	call MPI_Finalize(ierr)

end PROGRAM ButterflyPACK_FullKRR



subroutine FULLKER_solve(bmat,option,msh,quant,ptree,stats)

    use BPACK_DEFS
	use APPLICATION_MODULE
	use omp_lib
	use BPACK_Solve_Mul
	use MISC_Utilities

    implicit none

    integer i, j, ii, jj, iii, jjj, ierr, ntest,Dimn,edge_m,edge_n,ncorrect
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, n_iter_max, iter, N_unk, N_unk_loc
    real(kind=8) theta, phi, dphi, rcs_V, rcs_H,error
    real T0,T1
    real(kind=8) n1,n2,rtemp
    real(kind=8) value_Z
    real(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:),vout(:,:),vout_tmp(:,:),vout_test(:,:),matrixtemp1(:,:),matrixtemp2(:,:),matrixtemp(:,:)
	real(kind=8):: rel_error
	type(Hoption)::option
	type(mesh)::msh
	type(quant_app)::quant
	type(proctree)::ptree
	type(Bmatrix)::bmat
	type(Hstat)::stats
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
		call BPACK_Test_Solve_error(bmat,N_unk_loc,option,ptree,stats)
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


	n1 = OMP_get_wtime()

	call BPACK_Solution(bmat,x,b,N_unk_loc,1,option,ptree,stats)

	n2 = OMP_get_wtime()
	stats%Time_Sol = stats%Time_Sol + n2-n1
	call MPI_ALLREDUCE(stats%Time_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'Solving:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Flop_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,'(A13Es14.2)') 'Solve flops:',rtemp

	!**** prediction on the test sets


	T1=secnds(0.0)

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
		! call GeneralInverse(N_unk,N_unk,matrixtemp2,matrixtemp1,1d-10)

		! allocate(ipiv(N_unk))
		! allocate(matrixtemp1(N_unk,N_unk))
		! matrixtemp1=quant%matZ_glo(msh%new2old,msh%new2old)
		! ipiv=0
		! call getrff90(matrixtemp1,ipiv)
		! call getrif90(matrixtemp1,ipiv)
		! deallocate(ipiv)




		! allocate(matrixtemp(ntest,N_unk))
		! matrixtemp=quant%matZ_glo(quant%ntrain+1:quant%ntrain+quant%ntest,msh%new2old)
		! call gemmf90(matrixtemp1,N_unk,b,N_unk,x,N_unk,'N','N',N_unk,1,N_unk,cone,czero)
		! allocate (vout(ntest,1))
		! vout=0
		! call gemmf90(matrixtemp,ntest,x,N_unk,vout,ntest,'N','N',ntest,1,N_unk,cone,czero)


	if (ptree%MyID==Main_ID) then
		vout_tmp = vout-vout_test
		error=0
		do ii=1,ntest
			error = error + vout_tmp(ii,1)**2d0
		enddo
		error = error/ntest
		error = sqrt(error)

		write (*,*) ''
		write (*,*) 'Prediction time:',secnds(T1),'Seconds'
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

