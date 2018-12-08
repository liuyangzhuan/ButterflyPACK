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

module APPLICATION_MODULE
use d_BPACK_DEFS
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
			! if(m==n)value_e = value_e + 1d-4
		class default
			write(*,*)"unexpected type"
			stop
		end select	
	end subroutine Zelem_FULL

end module APPLICATION_MODULE	





PROGRAM HODLR_BUTTERFLY_SOLVER_FULLKER
    use d_BPACK_DEFS
    use APPLICATION_MODULE
	
	use d_BPACK_structure
	use d_BPACK_factor
	use d_BPACK_constr
	use omp_lib
	use d_misc

    implicit none

	! include "mkl_vml.fi"	 
	
    real(kind=8) para
    real(kind=8) tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seed_myid(50)
	integer times(8)	
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	
	character(len=:),allocatable  :: string
	character(len=1024)  :: strings	
	character(len=6)  :: info_env	
	integer :: length
	integer :: edge_m,edge_n
	integer :: ierr
	integer*8 oldmode,newmode
	type(d_Hoption)::option	
	type(d_Hstat)::stats
	type(d_mesh)::msh
	type(d_kernelquant)::ker
	type(quant_app),target::quant
	type(d_hobf)::ho_bf,ho_bf_copy
	integer,allocatable:: groupmembers(:)
	integer nmpi
	type(d_proctree)::ptree
	CHARACTER (LEN=1000) DATA_DIR	
	
	
	!**** nmpi and groupmembers should be provided by the user 
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo	
	
	!**** create the process tree
	call d_createptree(nmpi,groupmembers,MPI_Comm_World,ptree)
	deallocate(groupmembers)
	
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'NUMBER_MPI=',nmpi
	
	!**** set number of threads
 	threads_num=1
    CALL getenv("OMP_NUM_THREADS", strings)
	strings = TRIM(strings)	
	if(LEN_TRIM(strings)>0)then
		read(strings , *) threads_num
	endif
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'OMP_NUM_THREADS=',threads_num
	call OMP_set_num_threads(threads_num)		
		

	if(ptree%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "HODLR_BUTTERFLY_SOLVER_FULLKER"
    write(*,*) "   "
	endif
	
	!**** initialize statistics variables  
	call d_initstat(stats)
	call d_setdefaultoptions(option)
	
	!**** register the user-defined function and type in ker 
	ker%FuncZmn=>Zelem_FULL
	ker%QuantApp=>quant
 

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
	
    !**** set solver parameters	
	option%xyzsort=NATURAL
	option%nogeo=1
	option%Nmin_leaf=500 
	! option%tol_comp=1d-4
	option%tol_Rdetect=3d-5	
	option%tol_LS=1d-12
	option%tol_itersol=1d-6
	option%n_iter=1000
	option%tol_rand=1d-3
	option%level_check=10000
	option%precon=DIRECT
	! option%xyzsort=TM
	option%lnoBP=40000
	option%TwoLayerOnly=1
    option%schulzorder=3
    option%schulzlevel=3000
	option%LRlevel=0
	option%ErrFillFull=1
	option%ErrSol=1
	! option%RecLR_leaf=RRQR


   !***********************************************************************
   if(ptree%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'FULLKER computing'
   write (*,*) ''
   endif
   !***********************************************************************
	
	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "reading fullmatrix......"
 		
	msh%Nunk=quant%Nunk	
	
	allocate(quant%perms(quant%ntrain+quant%ntest))
	
	do ii=1,quant%ntrain+quant%ntest
		quant%perms(ii)=ii
	enddo
	
	! call d_rperm(quant%ntrain+quant%ntest,quant%perms)
	
	
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
	! write(*,*)t2-t1

	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "constructing HODLR formatting......"
    call d_Cluster_partition(ho_bf,option,msh,ker,d_element_Zmn_user,ptree)
	call d_HODLR_structuring(ho_bf,option,msh,ptree,stats)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "HODLR formatting finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1

    
    !call compression_test()
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "HODLR construction......"
    call d_BPACK_construction(ho_bf,option,stats,msh,ker,d_element_Zmn_user,ptree)
	! call copy_HOBF(ho_bf,ho_bf_copy)	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "HODLR construction finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	! write(*,*)t2-t1
	
	if(option%precon/=NOPRECON)then
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Cascading factorizing......"
    call d_BPACK_factorization(ho_bf,option,stats,ptree,msh)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Cascading factorizing finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "	
	end if

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Solve and Prediction......"
    call FULLKER_solve(ho_bf,option,msh,quant,ptree,stats)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Solve and Prediction finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "	
	
	
	call d_delete_proctree(ptree)
	call d_delete_Hstat(stats)
	call d_delete_mesh(msh)
	call d_delete_kernelquant(ker)
	call d_HODLR_delete(ho_bf)
	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"
	
	call blacs_exit(1)
	call MPI_Finalize(ierr)
	
end PROGRAM HODLR_BUTTERFLY_SOLVER_FULLKER



subroutine FULLKER_solve(ho_bf_for,option,msh,quant,ptree,stats)
    
    use d_BPACK_DEFS
	use APPLICATION_MODULE
	use omp_lib
	use d_BPACK_Solve_Mul
	use d_MISC
    
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
	type(d_Hoption)::option
	type(d_mesh)::msh
	type(quant_app)::quant
	type(d_proctree)::ptree
	type(d_hobf)::ho_bf_for
	type(d_Hstat)::stats	
	real(kind=8),allocatable:: current(:),voltage(:)
	real(kind=8), allocatable:: labels(:)
	integer, allocatable :: ipiv(:)
	real(kind=8) r_mn
	real(kind=8) label

	
	N_unk=msh%Nunk
	N_unk_loc = msh%idxe-msh%idxs+1
	ntest=quant%ntest
	
	if(option%ErrSol==1)then
		call d_BPACK_Test_Solve_error(ho_bf_for,N_unk_loc,option,ptree,stats)
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
	
	call d_BPACK_solution(ho_bf_for,x,b,N_unk_loc,1,option,ptree,stats)
	
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
		! call d_gemmf90(matrixtemp1,N_unk,b,N_unk,x,N_unk,'N','N',N_unk,1,N_unk,cone,czero)
		! allocate (vout(ntest,1))		
		! vout=0
		! call d_gemmf90(matrixtemp,ntest,x,N_unk,vout,ntest,'N','N',ntest,1,N_unk,cone,czero) 
	 
	
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
	
    return
    
end subroutine FULLKER_solve

