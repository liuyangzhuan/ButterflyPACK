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

! Developers: Yang Liu, Xiaoye S. Li.
!             (Lawrence Berkeley National Lab, Computational Research Division).

module APPLICATION_MODULE
use d_BPACK_DEFS
implicit none

	!**** define your application-related variables here   
	type quant_app
		real(kind=8) sigma, lambda ! Kernel Regression: RBF parameters 
		integer dimn  ! dimension of the data sets
		integer ntrain,ntest ! size of training points and test points
		character(LEN=500)::trainfile_p,trainfile_tree,trainfile_l,testfile_p,testfile_l !Kernel Regression: file pointers to train data, preordered tree, train labels, test data, and test labels
		
		integer Nunk ! size of the matrix 
		real(kind=8),allocatable:: xyz(:,:)   ! coordinates of the points		
		
	end type quant_app

contains


	!**** cutoff distance for gaussian kernel 
	real(kind=8) function arg_thresh_Zmn(quant)
		use d_BPACK_DEFS
		implicit none 
		
		type(quant_app)::quant
		arg_thresh_Zmn=-log(SafeUnderflow)*2.0*quant%sigma**2 
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





PROGRAM HODLR_BUTTERFLY_SOLVER_RBF
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
    integer i,j,k, threads_num,iii
	integer seed_myid(50)
	integer times(8)	
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	
	character(len=:),allocatable  :: string
	character(len=1024)  :: strings	
	character(len=6)  :: info_env	
	integer :: length
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
	integer MPI_thread
	
	call MPI_Init(ierr)
	
	
	do iii=1,1
	! iii=0
	! do while(iii==0)
		
	
	!**** nmpi and groupmembers should be provided by the user 
	
	! call MPI_Init_thread( MPI_THREAD_SINGLE, MPI_thread,ierr); 
	
	call MPI_Query_thread(MPI_thread,ierr);
	select case (MPI_thread)
	case (MPI_THREAD_SINGLE)
		write(*,*)"MPI_Query_thread with MPI_THREAD_SINGLE"
	case (MPI_THREAD_SERIALIZED)
		write(*,*)"MPI_Query_thread with MPI_THREAD_SERIALIZED"
	case (MPI_THREAD_MULTIPLE)
		write(*,*)"MPI_Query_thread with MPI_THREAD_MULTIPLE"
	case (MPI_THREAD_FUNNELED)
		write(*,*)"MPI_Query_thread with MPI_THREAD_FUNNELED"
	end select
	

	
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
		
		
	!**** create a random seed	
	call DATE_AND_TIME(values=times)     ! Get the current time 
	seed_myid(1) = times(4) * (360000*times(5) + 6000*times(6) + 100*times(7) + times(8))
	! seed_myid(1) = myid*1000
	call RANDOM_SEED(PUT=seed_myid)
	

	if(ptree%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "HODLR_BUTTERFLY_SOLVER_RBF"
    write(*,*) "   "
	endif
	
	!**** initialize statistics variables  
	call d_initstat(stats)
	call d_setdefaultoptions(option)
	
	!**** register the user-defined function and type in ker 
	ker%FuncZmn=>Zelem_RBF
	ker%QuantApp=>quant
 

    !**** read data file directory, data set dimension, size of training and testing sets, RBF parameters
	CALL getarg(1, DATA_DIR)
	quant%trainfile_p=trim(DATA_DIR)//'_train.csv'			
	quant%trainfile_l=trim(DATA_DIR)//'_train_label.csv'		
	quant%testfile_p=trim(DATA_DIR)//'_test.csv'			
	quant%testfile_l=trim(DATA_DIR)//'_test_label.csv'		
	call getarg(2,strings)
	read(strings,*)quant%dimn
	call getarg(3,strings)
	read(strings,*)quant%ntrain
	quant%Nunk = quant%ntrain
	call getarg(4,strings)
	read(strings,*)quant%ntest
	call getarg(5,strings)
	read(strings,*)quant%sigma
	call getarg(6,strings)
	read(strings,*)quant%lambda	
	call getarg(7,strings)
	read(strings,*)option%xyzsort
	call getarg(8,strings)
	read(strings,*)option%RecLR_leaf	
	

    !**** set solver parameters	
	
	option%nogeo=0
	option%Nmin_leaf=200
	option%tol_comp=1d-2
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
	option%ErrFillFull=0
	option%ErrSol=1
	! option%RecLR_leaf=RRQR


   !***********************************************************************
   if(ptree%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'RBF computing'
   write (*,*) ''
   endif
   !***********************************************************************
	
	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "geometry modeling......"
    call geo_modeling_RBF(quant,ptree%Comm)
	
	msh%Nunk=quant%Nunk
	allocate(msh%xyz(quant%dimn,quant%Nunk))
	msh%xyz=quant%xyz
	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "modeling finished"
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
    call RBF_solve(ho_bf,option,msh,quant,ptree,stats)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Solve and Prediction finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "	
	
	if(allocated(quant%xyz))deallocate(quant%xyz)
	call d_delete_proctree(ptree)
	call d_delete_Hstat(stats)
	call d_delete_mesh(msh)
	call d_delete_kernelquant(ker)
	call d_HODLR_delete(ho_bf)
	
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"
	
	end do
	call blacs_exit(1)
	call MPI_Finalize(ierr)
	
end PROGRAM HODLR_BUTTERFLY_SOLVER_RBF


!**** read training sets 
subroutine geo_modeling_RBF(quant,MPIcomm)

    use d_BPACK_DEFS
	use APPLICATION_MODULE
    implicit none
	type(quant_app)::quant
    integer i,j,ii,jj,iii,jjj
    integer intemp
    integer node, patch, edge, flag
    integer node1, node2, num_node,Dimn
    integer node_temp(2)
    real(kind=8) dx, xx, yy, rr, theta,L,M,Am,tt,L1,L2,L3
    
    real(kind=8),allocatable :: node_xy_original(:,:)
    integer,allocatable :: num_edge_of_node(:)
    
    real(kind=8) a(3),b(3),c(3),r0, phi_start
	! type(d_proctree)::ptree
	integer MPIcomm,MyID,ierr
	call MPI_Comm_rank(MPIcomm,MyID,ierr)
	
	Dimn = quant%dimn
	
	open (90,file=quant%trainfile_p)
	allocate (quant%xyz(Dimn,0:quant%Nunk))
	do edge=1,quant%Nunk
		read (90,*) quant%xyz(1:Dimn,edge)
	enddo  			
    close(90)
    return
    
end subroutine geo_modeling_RBF


subroutine RBF_solve(ho_bf_for,option,msh,quant,ptree,stats)
    
    use d_BPACK_DEFS
	use APPLICATION_MODULE
	use omp_lib
	use d_BPACK_Solve_Mul
    
    implicit none
    
    integer i, j, ii, jj, iii, jjj, ierr, ntest,Dimn,edge_m,edge_n,ncorrect
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, n_iter_max, iter, N_unk, N_unk_loc
    real(kind=8) theta, phi, dphi, rcs_V, rcs_H,rate
    real T0,T1
    real(kind=8) n1,n2,rtemp	
    real(kind=8) value_Z
    real(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:),vout(:,:),vout_tmp(:,:)
	real(kind=8):: rel_error
	type(d_Hoption)::option
	type(d_mesh)::msh
	type(quant_app)::quant
	type(d_proctree)::ptree
	type(d_hobf)::ho_bf_for
	type(d_Hstat)::stats	
	real(kind=8),allocatable:: current(:),voltage(:)
	real(kind=8), allocatable:: labels(:)
	real(kind=8),allocatable:: xyz_test(:,:)
	real(kind=8) r_mn
	integer label

	N_unk=msh%Nunk
	Dimn=quant%dimn
	N_unk_loc = msh%idxe-msh%idxs+1	
	

	if(option%ErrSol==1)then
		call d_BPACK_Test_Solve_error(ho_bf_for,N_unk_loc,option,ptree,stats)
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
	
	n1 = OMP_get_wtime()
	
	call d_BPACK_solution(ho_bf_for,x,b,N_unk_loc,1,option,ptree,stats)
	
	n2 = OMP_get_wtime()
	stats%Time_Sol = stats%Time_Sol + n2-n1
	call MPI_ALLREDUCE(stats%Time_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'Solving:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Flop_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,'(A13Es14.2)') 'Solve flops:',rtemp		
	
	!**** prediction on the test sets	
	
	ntest=quant%ntest
	T1=secnds(0.0)
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
			r_mn=sum((xyz_test(1:dimn,edge_m)-msh%xyz(1:dimn,msh%new2old(edge+msh%idxs-1)))**2)
			value_Z = exp(-r_mn/2.0/quant%sigma**2)
			vout_tmp(edge_m,1) = vout_tmp(edge_m,1) + value_Z*x(edge,1)
		enddo
	enddo	
	
	call MPI_REDUCE(vout_tmp, vout, ntest,MPI_double_precision, MPI_SUM, Main_ID, ptree%Comm,ierr)
			
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
		write (*,*) 'Prediction time:',secnds(T1),'Seconds'
		write (*,*) 'Success rate:',rate
		write (*,*) ''
		
	endif		

	deallocate (vout)
	deallocate (vout_tmp)

	deallocate(x)
	deallocate(b)	
	deallocate(xyz_test)
	
	call MPI_barrier(ptree%Comm,ierr)
	
    return
    
end subroutine RBF_solve

