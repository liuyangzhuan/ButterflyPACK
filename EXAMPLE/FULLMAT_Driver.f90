module APPLICATION_MODULE
use d_HODLR_DEFS
implicit none

	!**** define your application-related variables here   
	type quant_app
		real(kind=8), allocatable :: matU_glo(:,:),matV_glo(:,:) ! Full Matrix: the random LR matrix to sample its entries
		real(kind=8), allocatable :: matZ_glo(:,:) ! Full Matrix: Full matrix read from files		
		integer:: rank
		real(kind=8):: lambda
	end type quant_app

contains
	
	!**** user-defined subroutine to sample Z_mn as two LR products
	subroutine Z_elem_LR(m,n,value_e,quant)
		use d_HODLR_DEFS
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
	end subroutine Z_elem_LR


	!**** user-defined subroutine to sample Z_mn as full matrix
	subroutine Z_elem_FULL(m,n,value_e,quant)
		use d_HODLR_DEFS
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
	end subroutine Z_elem_FULL
	
	
	
end module APPLICATION_MODULE	





PROGRAM HODLR_BUTTERFLY_SOLVER
    use d_HODLR_DEFS
    use APPLICATION_MODULE
	use d_HODLR_Solve_Mul
	
	use d_HODLR_structure
	use d_HODLR_factor
	use d_HODLR_constr
	use omp_lib
	use d_misc
	use d_HODLR_constr
	
    implicit none

	! include "mkl_vml.fi"	 
	
    real(kind=8) para
    real(kind=8) tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seed_myid(50)
	integer times(8)	
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	real(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:)
	real(kind=8),allocatable:: datain(:)
	
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
	integer level,Maxlevel
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
	
	if(ptree%MyID==Main_ID)write(*,*)'NUMBER_MPI=',nmpi
	
	!**** set number of threads
 	threads_num=1
    CALL getenv("OMP_NUM_THREADS", strings)
	strings = TRIM(strings)	
	if(LEN_TRIM(strings)>0)then
		read(strings , *) threads_num
	endif
	if(ptree%MyID==Main_ID)write(*,*)'OMP_NUM_THREADS=',threads_num
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
	
	
	
!******************************************************************************!
! generate a LR matrix as two matrix product	
	
	
	
	! !**** register the user-defined function and type in ker 
	! ker%FuncZmn=>Z_elem_LR
	! ker%QuantZmn=>quant
 
    ! !**** Get matrix size and rank and create the matrix
	! msh%Nunk = 10000
	! quant%rank = 10
	! quant%lambda = 1d5
	! allocate(quant%matU_glo(msh%Nunk,quant%rank))
	! call RandomMat(msh%Nunk,quant%rank,quant%rank,quant%matU_glo,0)
	! call MPI_Bcast(quant%matU_glo,msh%Nunk*quant%rank,MPI_DOUBLE_COMPLEX,Main_ID,ptree%Comm,ierr)
	
	! allocate(quant%matV_glo(quant%rank,msh%Nunk))
	! call RandomMat(quant%rank,msh%Nunk,quant%rank,quant%matV_glo,0)	
	! call MPI_Bcast(quant%matV_glo,msh%Nunk*quant%rank,MPI_DOUBLE_COMPLEX,Main_ID,ptree%Comm,ierr)	
	

	
	
!******************************************************************************!
! generate a LR matrix stored in a files
	
	CALL getarg(1, strings)
	strings = TRIM(strings)	
	if(LEN_TRIM(strings)==0)then
		strings = './EXAMPLE/K05N4096.csv'	
	endif
	
	
	!**** register the user-defined function and type in ker 
	ker%FuncZmn=>Z_elem_FULL
	ker%QuantZmn=>quant

    !**** Get matrix size and rank and create the matrix
	msh%Nunk = 4096
	allocate(quant%matZ_glo(msh%Nunk,msh%Nunk))
	allocate(datain(msh%Nunk))
	open(10, file=strings)
	do ii=1,msh%Nunk
		read(10,*) datain(:)
		quant%matZ_glo(:,ii)=datain
	enddo
	close(10)
	call MPI_Bcast(quant%matZ_glo,msh%Nunk*msh%Nunk,MPI_DOUBLE_PRECISION,Main_ID,ptree%Comm,ierr)
	
		
!******************************************************************************!	
	
	
	
    !**** set solver parameters	
	
	option%nogeo=1
	option%Nmin_leaf=200
	option%tol_comp=1d-4
	option%tol_Rdetect=3d-5	
	option%tol_LS=1d-12
	option%tol_itersol=1d-6
	option%n_iter=1000
	option%tol_rand=1d-3
	option%level_check=10000
	option%precon=DIRECT
	option%xyzsort=NATURAL
	option%lnoBP=40000
	option%TwoLayerOnly=1
    option%schulzorder=3
    option%schulzlevel=3000
	option%LRlevel=0
	option%ErrFillFull=1
	option%ErrSol=1
	option%RecLR_leaf=RRQR
	
	
	
	CALL getarg(2, strings)
	strings = TRIM(strings)	
	if(LEN_TRIM(strings)>0)then
		read(strings,*)option%RecLR_leaf
	endif	
	
   !***********************************************************************
   if(ptree%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'Random LR Kernel computing'
   write (*,*) ''
   endif
   !***********************************************************************
	
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "constructing HODLR formatting......"
    call d_HODLR_structuring(ho_bf,option,msh,ker,d_element_Zmn_user,ptree)
	call d_BPlus_structuring(ho_bf,option,msh,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "HODLR formatting finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1

    
    !call compression_test()
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "HODLR construction......"
    call d_HODLR_construction(ho_bf,option,stats,msh,ker,d_element_Zmn_user,ptree)
	! call copy_HOBF(ho_bf,ho_bf_copy)	
    if(ptree%MyID==Main_ID)write(*,*) "HODLR construction finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	! write(*,*)t2-t1
	
	if(option%precon/=NOPRECON)then
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing......"
    call d_HODLR_Factorization(ho_bf,option,stats,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	end if
	
	if(option%ErrSol==1)then
    if(ptree%MyID==Main_ID)write(*,*) "Test Solve ......"
		call d_HODLR_Test_Solve_error(ho_bf,option,ptree,stats)
    if(ptree%MyID==Main_ID)write(*,*) "Test Solve finished"
	endif
	
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	
	deallocate(quant%matU_glo)
	deallocate(quant%matV_glo)
	
	
    if(ptree%MyID==Main_ID)write(*,*) "-------------------------------program end-------------------------------------"

	
end PROGRAM HODLR_BUTTERFLY_SOLVER

