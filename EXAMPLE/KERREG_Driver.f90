module APPLICATION_MODULE
use MODULE_FILE
implicit none

	!**** define your application-related variables here   
	type quant_app
		real*8 sigma, lambda ! Kernel Regression: RBF parameters 		
		character(LEN=500)::trainfile_p,trainfile_tree,trainfile_l,testfile_p,testfile_l !Kernel Regression: file pointers to train data, preordered tree, train labels, test data, and test labels
	end type quant_app

contains


	!**** cutoff distance for gaussian kernel 
	real*8 function arg_thresh_Zmn(quant)
		use MODULE_FILE
		implicit none 
		
		type(quant_app)::quant
		arg_thresh_Zmn=-log(SafeUnderflow)*2.0*quant%sigma**2 
	end function arg_thresh_Zmn


	
	!**** user-defined subroutine to sample Z_mn
	subroutine Z_elem_RBF(ker,m,n,value_e,msh,quant)
		use MODULE_FILE
		implicit none 
		
		class(kernelquant)::ker ! this is required if F_Z_elem is a procedure pointer defined in type kernelquant
		class(*),pointer :: quant
		type(mesh)::msh
		integer, INTENT(IN):: m,n
		complex(kind=8)::value_e 

		real(kind=8) r_mn
		integer dimn
		
		select TYPE(quant)
		type is (quant_app)
			dimn = size(msh%xyz,1)
			value_e=0
			r_mn=sum((msh%xyz(1:dimn,m)-msh%xyz(1:dimn,n))**2)
			if(r_mn<arg_thresh_Zmn(quant))then
				value_e = exp(-r_mn/2.0/quant%sigma**2)
			endif
			if(r_mn==0)value_e = value_e + quant%lambda		
		class default
			write(*,*)"unexpected type"
			stop
		end select	
	end subroutine Z_elem_RBF

end module APPLICATION_MODULE	





PROGRAM HODLR_BUTTERFLY_SOLVER_RBF
    use MODULE_FILE
    use APPLICATION_MODULE
	! use geometry_model
	use H_structure
	use cascading_factorization
	use matrices_fill
	use omp_lib
	use MISC
	use matrices_fill
	
    implicit none

	! include "mkl_vml.fi"	 
	
    real*8 para
    real*8 tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seed_myid(12)
	integer times(8)	
	real*8 t1,t2,x,y,z,r,theta,phi
	complex(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:)
	
	character(len=:),allocatable  :: string
	character(len=1024)  :: strings	
	character(len=6)  :: info_env	
	integer :: length
	integer :: ierr
	integer*8 oldmode,newmode
	type(Hoption)::option	
	type(Hstat)::stats
	type(mesh)::msh
	type(kernelquant)::ker
	type(quant_app),target::quant
	type(hobf)::ho_bf,ho_bf_copy
	integer,allocatable:: groupmembers(:)
	integer nmpi
	type(proctree)::ptree
	
	! nmpi and groupmembers should be provided by the user 
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo	
	
	call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	
	if(ptree%MyID==Main_ID)write(*,*)'NUMBER_MPI=',nmpi
	
 	threads_num=1
    CALL getenv("OMP_NUM_THREADS", strings)
	strings = TRIM(strings)	
	if(LEN_TRIM(strings)>0)then
		read(strings , *) threads_num
	endif
	if(ptree%MyID==Main_ID)write(*,*)'OMP_NUM_THREADS=',threads_num
	
	call OMP_set_num_threads(threads_num)		
		
		
	call DATE_AND_TIME(values=times)     ! Get the current time 
	seed_myid(1) = times(4) * (360000*times(5) + 6000*times(6) + 100*times(7) + times(8))
	
	! seed_myid(1) = myid*1000
	call RANDOM_SEED(PUT=seed_myid)
	
	! oldmode = vmlsetmode(VML_FTZDAZ_ON)
	! call vmlsetmode(VML_FTZDAZ_ON)
	if(ptree%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "HODLR_BUTTERFLY_SOLVER_RBF"
    write(*,*) "   "
	endif
	
	call InitStat(stats)
	
	! time_indexarray = 0
	! time_leastsquare = 0
	! time_buttermul = 0
	! time_buttermulinv = 0
	! time_kernelupdate = 0
	! time_memcopy = 0
	! time_gemm = 0
	! time_gemm1 = 0		   
    ! time_getvec = 0
	! time_resolve = 0
	! time_halfbuttermul = 0
	time_tmp = 0
	
	msh%Origins=(/0d0,0d0,0d0/)
 
     
	! register the user-defined function and type in ker 
	ker%FuncZmn=>Z_elem_RBF
	ker%QuantZmn=>quant
 

     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	
	
	
	CALL getarg(1, DATA_DIR)
	quant%trainfile_p=trim(DATA_DIR)//'_train.csv'			
	quant%trainfile_l=trim(DATA_DIR)//'_train_label.csv'		
	quant%testfile_p=trim(DATA_DIR)//'_test.csv'			
	quant%testfile_l=trim(DATA_DIR)//'_test_label.csv'		
	
	
	para=0.001
		
	! msh%Nunk=100000

	msh%scaling=1d0
	quant%sigma=0.1d0
	quant%lambda=10.0d0


	
	option%Nmin_leaf=200
	option%tol_SVD=1d-4
	option%tol_Rdetect=3d-5	
	option%tol_LS=1d-12
	option%tol_itersol=1d-6
	option%N_iter=1000
	option%tol_rand=1d-3
	option%level_check=10000
	option%precon=DIRECT
	option%xyzsort=3
	option%LnoBP=40000
	option%TwoLayerOnly=1
	option%touch_para=3
    option%schulzorder=3
    option%schulzlevel=3000
	option%LRlevel=0
	option%ErrFillFull=0
	option%RecLR_leaf='Q'

	! call MKL_set_num_threads(NUM_Threads)    ! this overwrites omp_set_num_threads for MKL functions 
	
    ! call OMP_set_dynamic(.true.)
    !call OMP_set_nested(.true.)
	

    !*********************************************************

   !***********************************************************************
   if(ptree%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'RBF computing'
   write (*,*) ''
   endif
   !***********************************************************************
	
	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID)write(*,*) "geometry modeling......"
    call geo_modeling_RBF(msh,quant,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "modeling finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1

	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "constructing H_matrices formatting......"
    call H_matrix_structuring(ho_bf,para,option,msh,ptree)
	call BPlus_structuring(ho_bf,option,msh,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices formatting finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1
	! stop 
	
    !pause
    
    !call compression_test()
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices filling......"
    call matrices_filling(ho_bf,option,stats,msh,ker,element_Zmn_user,ptree)
	! if(option%precon/=DIRECT)then
		call copy_HOBF(ho_bf,ho_bf_copy)	
	! end if
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices filling finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	! write(*,*)t2-t1
	
	if(option%precon/=NOPRECON)then
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing......"
    call cascading_factorizing(ho_bf,option,stats,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	end if

    if(ptree%MyID==Main_ID)write(*,*) "Solve and Prediction......"
    call RBF_solve(ho_bf_copy,ho_bf,option,msh,quant,ptree,stats)
    if(ptree%MyID==Main_ID)write(*,*) "Solve and Prediction finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	
	
    if(ptree%MyID==Main_ID)write(*,*) "-------------------------------program end-------------------------------------"

    ! ! ! ! pause

end PROGRAM HODLR_BUTTERFLY_SOLVER_RBF



subroutine geo_modeling_RBF(msh,quant,ptree)

    use MODULE_FILE
	use APPLICATION_MODULE
    implicit none
    type(mesh)::msh
	type(quant_app)::quant
    integer i,j,ii,jj,iii,jjj
    integer intemp
    integer node, patch, edge, flag
    integer node1, node2, num_node,Maxedge,Dimn
    integer node_temp(2)
    real(kind=8) dx, xx, yy, rr, theta,L,M,Am,tt,L1,L2,L3
    
    real(kind=8),allocatable :: node_xy_original(:,:)
    integer,allocatable :: num_edge_of_node(:)
    
    real(kind=8) a(3),b(3),c(3),r0, phi_start
	type(proctree)::ptree
	
	! Maxedge	= msh%Nunk
	msh%Ncorner = 0
	
	open (90,file=quant%trainfile_p)
	read (90,*) Maxedge, Dimn
	
	msh%Nunk = Maxedge
	
	allocate (msh%xyz(Dimn,0:Maxedge), msh%info_unk(0:0,Maxedge))
	! write(*,*)Maxedge, Dimn,shape(msh%info_unk)
	do edge=1,Maxedge
		msh%info_unk(0,edge)=edge
		read (90,*) msh%xyz(1:Dimn,edge)
	enddo  			

    return
    
end subroutine geo_modeling_RBF


subroutine RBF_solve(ho_bf_for,ho_bf_inv,option,msh,quant,ptree,stats)
    
    use MODULE_FILE
	use APPLICATION_MODULE
	use omp_lib
	use HODLR_Solve
    ! use blas95
    implicit none
    
    integer i, j, ii, jj, iii, jjj, ierr, ntest,Dimn,edge_m,edge_n,ncorrect
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, N_iter_max, iter, N_unk, N_unk_loc
    real*8 theta, phi, dphi, rcs_V, rcs_H,rate
    real T0,T1
    real*8 n1,n2,rtemp	
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:),vout(:,:),vout_tmp(:,:)
	real*8:: rel_error
	type(Hoption)::option
	type(mesh)::msh
	type(quant_app)::quant
	type(proctree)::ptree
	type(hobf)::ho_bf_for,ho_bf_inv
	type(Hstat)::stats	
	complex(kind=8),allocatable:: current(:),voltage(:)
	complex(kind=8), allocatable:: labels(:)
	real*8,allocatable:: xyz_test(:,:)
	real(kind=8) r_mn
	integer label


	if(option%ErrSol==1)then
		call HODLR_Test_Solve_error(ho_bf_for,ho_bf_inv,option,ptree,stats)
	endif	
	
	! if(option%PRECON==DIRECT)then
		! msh%idxs = ho_bf_inv%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		! msh%idxe = ho_bf_inv%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)
	! else 
		! ! write(*,*)associated(ho_bf_for%levels(1)%BP_inverse),'dd' !%matrices_block(1)%N_p),'nima'
		! msh%idxs = ho_bf_for%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		! msh%idxe = ho_bf_for%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)	
	! endif
	
	N_unk=msh%Nunk
	N_unk_loc = msh%idxe-msh%idxs+1
	
	allocate(labels(N_unk))
	allocate (x(N_unk_loc,1))
	x=0	
	allocate (b(N_unk_loc,1))
	open (91,file=quant%trainfile_l)
	read(91,*)N_unk,Dimn
	do ii=1,N_unk
		read(91,*)label
		labels(ii)=label
	enddo
	do ii=1,N_unk_loc
		b(ii,1) = labels(msh%info_unk(0,ii-1+msh%idxs))
	enddo
	deallocate(labels)	
	
	
	n1 = OMP_get_wtime()
	
	call HODLR_Solution(ho_bf_for,ho_bf_inv,x,b,N_unk_loc,1,option,ptree,stats)
	
	n2 = OMP_get_wtime()
	stats%Time_Sol = stats%Time_Sol + n2-n1
	call MPI_ALLREDUCE(stats%Time_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,*) 'Solving:',rtemp,'Seconds'
	call MPI_ALLREDUCE(stats%Flop_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)write (*,'(A13Es14.2)') 'Solve flops:',rtemp		
	
	
	! prediction on the test sets
	T1=secnds(0.0)
	open (92,file=quant%testfile_p)
	read (92,*) ntest, Dimn
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
			r_mn=sum((xyz_test(1:dimn,edge_m)-msh%xyz(1:dimn,msh%info_unk(0,edge+msh%idxs-1)))**2)
			value_Z = exp(-r_mn/2.0/quant%sigma**2)
			vout_tmp(edge_m,1) = vout_tmp(edge_m,1) + value_Z*x(edge,1)
		enddo
	enddo	
	
	call MPI_REDUCE(vout_tmp, vout, ntest,MPI_double_complex, MPI_SUM, Main_ID, ptree%Comm,ierr)
			
	if (ptree%MyID==Main_ID) then
		do ii=1,ntest
			if(dble(vout(ii,1))>0)then
				vout(ii,1)=1
			else
				vout(ii,1)=-1
			endif
		enddo

		
		open (93,file=quant%testfile_l)
		read (93,*) ntest, Dimn
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

