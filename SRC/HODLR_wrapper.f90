#include "HODLR_config.fi"
module HODLR_wrapper
    use HODLR_DEFS
	
	use HODLR_structure
	use HODLR_factor
	use HODLR_constr
	use HODLR_randomMVP
	use omp_lib
	use misc
	use HODLR_Solve_Mul
    use iso_c_binding    
	

contains 



subroutine element_Zmn_user_C(edge_m,edge_n,value_e,msh,ker)
    use HODLR_DEFS
    implicit none
    integer edge_m, edge_n
    DT value_e
	type(mesh)::msh
	type(kernelquant)::ker	
	procedure(C_Zelem), POINTER :: proc
	
	value_e=0
	call c_f_procpointer(ker%C_FuncZmn, proc)
	call proc(msh%new2old(edge_m)-1,msh%new2old(edge_n)-1,value_e,ker%C_QuantApp)
	return
    
end subroutine element_Zmn_user_C


subroutine matvec_user_C(trans,M,N,num_vect,Vin,Vout,ker)

	integer, INTENT(IN):: M,N,num_vect
	DT::Vin(:,:),Vout(:,:) 
	type(kernelquant)::ker
	character trans
	
	procedure(C_HMatVec), POINTER :: proc
	
	call c_f_procpointer(ker%C_FuncHMatVec, proc)
	if(trans=='N')then  ! note that C_HMatVec takes Nin Nout, instead of takes M,N, so if statement is needed here
		call proc(trans//c_null_char,N,M,num_vect,Vin,Vout,ker%C_QuantApp)
	else
		call proc(trans//c_null_char,M,N,num_vect,Vin,Vout,ker%C_QuantApp)	
	endif
	return
	
end subroutine matvec_user_C	



subroutine Bmatvec_user_C(ker,block_o,trans,M,N,num_vect,Vin,Vout,a,b,ptree,stats,operand1)
	implicit none
	class(*)::ker	
	class(*),optional::operand1	
	type(matrixblock)::block_o
	character trans
	integer M, N, num_vect
	type(proctree)::ptree
	type(Hstat)::stats
	DT :: Vin(:,:), Vout(:,:),a,b
	procedure(C_BMatVec), POINTER :: proc
	
	select TYPE(ker)
	type is (kernelquant)
		call c_f_procpointer(ker%C_FuncBMatVec, proc)
		if(trans=='N')then  ! note that C_HMatVec takes Nin Nout, instead of takes M,N, so if statement is needed here
			call proc(trans//c_null_char,N,M,num_vect,Vin,Vout,ker%C_QuantApp,a,b)
		else
			call proc(trans//c_null_char,M,N,num_vect,Vin,Vout,ker%C_QuantApp,a,b)	
		endif			
	end select
end subroutine Bmatvec_user_C




!**** C interface of process tree construction
	!nmpi: number of MPIs for one hodlr
	!MPIcomm: MPI communicator from C caller
	!groupmembers: MPI ranks in MPIcomm for one hodlr
	!ptree_Cptr: the structure containing process tree
subroutine C_HODLR_Createptree(nmpi,groupmembers,MPIcomm,ptree_Cptr) bind(c, name="c_hodlr_createptree")	
	implicit none 
	integer nmpi
	integer MPIcomm
	integer:: groupmembers(*)
	type(c_ptr):: ptree_Cptr	
	type(proctree),pointer::ptree	

	allocate(ptree)
	call CreatePtree(nmpi,groupmembers,MPIcomm,ptree)
	ptree_Cptr=c_loc(ptree)
end subroutine C_HODLR_Createptree


!**** C interface of initializing statistics
	!nmpi: number of MPIs for one hodlr
	!MPIcomm: MPI communicator from C caller
	!groupmembers: MPI ranks in MPIcomm for one hodlr
	!stats_Cptr: the structure containing statistics
subroutine C_HODLR_Createstats(stats_Cptr) bind(c, name="c_hodlr_createstats")	
	implicit none 
	type(c_ptr), intent(out) :: stats_Cptr
	type(Hstat),pointer::stats

	allocate(stats)
	!**** initialize statistics variables  	
	call InitStat(stats)	
	stats_Cptr=c_loc(stats)
	
end subroutine C_HODLR_Createstats



!**** C interface of printing statistics
	!stats_Cptr: the structure containing statistics
	!ptree_Cptr: the structure containing process tree
subroutine C_HODLR_Printstats(stats_Cptr,ptree_Cptr) bind(c, name="c_hodlr_printstats")	
	implicit none 
	type(c_ptr) :: stats_Cptr
	type(c_ptr) :: ptree_Cptr
	type(Hstat),pointer::stats
	type(proctree),pointer::ptree	

	call c_f_pointer(stats_Cptr, stats)	
	call c_f_pointer(ptree_Cptr, ptree)
	!**** print statistics variables  	
	call PrintStat(stats,ptree)	
	
end subroutine C_HODLR_Printstats


!**** C interface of initializing option
	!option_Cptr: the structure containing option       
subroutine C_HODLR_Createoption(option_Cptr) bind(c, name="c_hodlr_createoption")	
	implicit none 
	type(c_ptr) :: option_Cptr
	type(Hoption),pointer::option	
	
	allocate(option)
	!**** set default hodlr options 	
	call SetDefaultOptions(option)
	
	option_Cptr=c_loc(option)
	
end subroutine C_HODLR_Createoption


!**** C interface of copy option
	!option_Cptr: the structure containing option       
	!option_Cptr1: the structure containing option       
subroutine C_HODLR_Copyoption(option_Cptr,option_Cptr1) bind(c, name="c_hodlr_copyoption")	
	implicit none 
	type(c_ptr) :: option_Cptr,option_Cptr1
	type(Hoption),pointer::option,option1	
	
	
	call c_f_pointer(option_Cptr, option)
	
	!****copy hodlr options 	
	allocate(option1)
	call CopyOptions(option,option1)
	
	option_Cptr1=c_loc(option1)
	
end subroutine C_HODLR_Copyoption



!**** C interface of set one entry in option
	!option_Cptr: the structure containing option       
subroutine C_HODLR_Setoption(option_Cptr,nam,val_Cptr) bind(c, name="c_hodlr_setoption")	
	implicit none 
	type(c_ptr) :: option_Cptr
	character(kind=c_char,len=1) :: nam(*)
	type(Hoption),pointer::option	
	! character::nam(:)	
	type(c_ptr),value :: val_Cptr
	integer,pointer::val_i
	real(kind=8),pointer::val_d
	integer strlen
	character(len=:),allocatable :: str
	integer valid_opt

	valid_opt = 0
	strlen=1
	do while(nam(strlen) /= c_null_char)
		strlen = strlen + 1
	enddo
	strlen = strlen -1
	allocate(character(len=strlen) :: str)
	str = transfer(nam(1:strlen), str)
	
	call c_f_pointer(option_Cptr, option)
           
	
!**** integer parameters 	
	if(trim(str)=='verbosity')then
	call c_f_pointer(val_Cptr, val_i)
	option%verbosity=val_i
	valid_opt=1
	endif
	if(trim(str)=='n_iter')then
	call c_f_pointer(val_Cptr, val_i)
	option%n_iter=val_i
	valid_opt=1
	endif
	if(trim(str)=='precon')then
	call c_f_pointer(val_Cptr, val_i)
	option%precon=val_i
	valid_opt=1
	endif
	if(trim(str)=='xyzsort')then
	call c_f_pointer(val_Cptr, val_i)
	option%xyzsort=val_i
	valid_opt=1
	endif
	if(trim(str)=='lnoBP')then
	call c_f_pointer(val_Cptr, val_i)
	option%lnoBP=val_i
	valid_opt=1
	endif
	if(trim(str)=='TwoLayerOnly')then
	call c_f_pointer(val_Cptr, val_i)
	option%TwoLayerOnly=val_i
	valid_opt=1
	endif
	if(trim(str)=='schulzorder')then
	call c_f_pointer(val_Cptr, val_i)
	option%schulzorder=val_i
	valid_opt=1
	endif
	if(trim(str)=='schulzlevel')then
	call c_f_pointer(val_Cptr, val_i)
	option%schulzlevel=val_i
	valid_opt=1
	endif
	if(trim(str)=='LRlevel')then
	call c_f_pointer(val_Cptr, val_i)
	option%LRlevel=val_i
	valid_opt=1	
	endif
	if(trim(str)=='ErrFillFull')then
	call c_f_pointer(val_Cptr, val_i)
	option%ErrFillFull=val_i
	valid_opt=1
	endif
	if(trim(str)=='BACA_Batch')then
	call c_f_pointer(val_Cptr, val_i)
	option%BACA_Batch=val_i
	valid_opt=1
	endif	
	if(trim(str)=='ErrSol')then
	call c_f_pointer(val_Cptr, val_i)
	option%ErrSol=val_i
	valid_opt=1
	endif
	if(trim(str)=='nogeo')then
	call c_f_pointer(val_Cptr, val_i)
	option%nogeo=val_i
	valid_opt=1
	endif
	if(trim(str)=='RecLR_leaf')then
	call c_f_pointer(val_Cptr, val_i)
	option%RecLR_leaf=val_i
	valid_opt=1
	endif
	if(trim(str)=='Nmin_leaf')then
	call c_f_pointer(val_Cptr, val_i)
	option%Nmin_leaf=val_i
	valid_opt=1
	endif	
	if(trim(str)=='LR_BLK_NUM')then
	call c_f_pointer(val_Cptr, val_i)
	option%LR_BLK_NUM=val_i
	valid_opt=1
	endif	
	if(trim(str)=='rank0')then
	call c_f_pointer(val_Cptr, val_i)
	option%rank0=val_i
	valid_opt=1
	endif	
	if(trim(str)=='itermax')then
	call c_f_pointer(val_Cptr, val_i)
	option%itermax=val_i
	valid_opt=1
	endif	
	if(trim(str)=='powiter')then
	call c_f_pointer(val_Cptr, val_i)
	option%powiter=val_i
	valid_opt=1
	endif
	
!**** double parameters 
	if(trim(str)=='tol_comp')then
	call c_f_pointer(val_Cptr, val_d)
	option%tol_comp=val_d
	valid_opt=1
	endif
	if(trim(str)=='tol_Rdetect')then
	call c_f_pointer(val_Cptr, val_d)
	option%tol_Rdetect=val_d
	valid_opt=1
	endif
	if(trim(str)=='tol_LS')then
	call c_f_pointer(val_Cptr, val_d)
	option%tol_LS=val_d
	valid_opt=1
	endif
	if(trim(str)=='tol_itersol')then
	call c_f_pointer(val_Cptr, val_d)
	option%tol_itersol=val_d
	valid_opt=1
	endif
	if(trim(str)=='tol_rand')then
	call c_f_pointer(val_Cptr, val_d)
	option%tol_rand=val_d
	valid_opt=1
	endif
	if(trim(str)=='touch_para')then
	call c_f_pointer(val_Cptr, val_d)
	option%touch_para=val_d
	valid_opt=1
	endif	
	if(trim(str)=='rankrate')then
	call c_f_pointer(val_Cptr, val_d)
	option%rankrate=val_d
	valid_opt=1
	endif	
    
	if(valid_opt==0)write(*,*)'invalid HODLR option: '//trim(str)
	
	deallocate(str)
	option_Cptr=c_loc(option)
	
end subroutine C_HODLR_Setoption



!**** C interface of HODLR construction
	!Npo: matrix size
	!Ndim: data set dimensionality (not used if nogeo=1)
	!Locations: coordinates used for clustering (not used if nogeo=1) 
	!Nmin: leafsize in HODLR tree 
	!tol: compression tolerance
	!nlevel: the number of top levels that have been ordered
	!tree: the order tree provided by the caller
	!Permutation: return the permutation vector new2old (indexed from 1)
	!Npo_loc: number of local row/column indices    
	!ho_bf_Cptr: the structure containing HODLR         
	!option_Cptr: the structure containing option         
	!stats_Cptr: the structure containing statistics         
	!msh_Cptr: the structure containing points and ordering information         
	!ker_Cptr: the structure containing kernel quantities
	!ptree_Cptr: the structure containing process tree
	!C_FuncZmn: the C_pointer to user-provided function to sample mn^th entry of the matrix
	!C_QuantApp: the C_pointer to user-defined quantities required to for entry evaluation and sampling
	!MPIcomm: user-provided MPI communicator
subroutine C_HODLR_Construct(Npo,Ndim,Locations,nlevel,tree,Permutation,Npo_loc,ho_bf_Cptr,option_Cptr,stats_Cptr,msh_Cptr,ker_Cptr,ptree_Cptr,C_FuncZmn,C_QuantApp,MPIcomm) bind(c, name="c_hodlr_construct")	
	implicit none 
	integer Npo,Ndim
	real(kind=8) Locations(*)
	
    real(kind=8) para
    real(kind=8) tolerance,h,lam
    integer Primary_block, nn, mm, MyID_old,Maxlevel,give,need
    integer i,j,k,ii,edge, threads_num,nth,Dimn,nmpi,ninc, acam
	real(kind=8),parameter :: cd = 299792458d0
	integer,allocatable:: groupmembers(:)
	integer nlevel,level
	integer Permutation(*),tree(*)
	integer Npo_loc
	! type(matricesblock), pointer :: blocks_i
	integer groupm
	integer MPIcomm
	type(c_ptr) :: ho_bf_Cptr
	type(c_ptr) :: option_Cptr
	type(c_ptr) :: stats_Cptr
	type(c_ptr) :: msh_Cptr
	type(c_ptr) :: ker_Cptr
	type(c_ptr) :: ptree_Cptr
	type(c_ptr), intent(in),target :: C_QuantApp
	type(c_funptr), intent(in),value,target :: C_FuncZmn
	
	type(Hoption),pointer::option	
	type(Hstat),pointer::stats
	type(mesh),pointer::msh
	type(kernelquant),pointer::ker
	type(hobf),pointer::ho_bf
	type(proctree),pointer::ptree	
	integer seed_myid(50)
	integer times(8)	
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	character(len=1024)  :: strings
	
	!**** allocate HODLR solver structures 
	allocate(ho_bf)
	! allocate(option)
	! allocate(stats)
	allocate(msh)
	allocate(ker)

	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	


	
	if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*)'NUMBER_MPI=',ptree%nproc
 	
 	threads_num=1
    CALL getenv("OMP_NUM_THREADS", strings)
	strings = TRIM(strings)	
	if(LEN_TRIM(strings)>0)then
		read(strings , *) threads_num
	endif
	if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*)'OMP_NUM_THREADS=',threads_num
	call OMP_set_num_threads(threads_num)
	
	!**** create a random seed		
	call DATE_AND_TIME(values=times)     ! Get the current time 
	seed_myid(1) = times(4) * (360000*times(5) + 6000*times(6) + 100*times(7) + times(8))
	! seed_myid(1) = myid*1000
	! call RANDOM_SEED(PUT=seed_myid)
	call init_random_seed()
	
	if(ptree%MyID==Main_ID .and. option%verbosity>0)then
    write(*,*) "HODLR_BUTTERFLY_SOLVER"
    write(*,*) "   "
	endif
	
	!**** register the user-defined function and type in ker 
	ker%C_QuantApp => C_QuantApp
	ker%C_FuncZmn => C_FuncZmn	
	

	msh%Nunk = Npo
		
		
	t1 = OMP_get_wtime()

	
	if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "User-supplied kernel:"
	Maxlevel=nlevel
	allocate(msh%pretree(2**Maxlevel))
	
	msh%pretree(1:2**Maxlevel) = tree(1:2**Maxlevel)	

	
	!**** make 0-element node a 1-element node 
 	
	! write(*,*)'before adjustment:',msh%pretree
	need = 0
	do ii=1,2**Maxlevel
		if(msh%pretree(ii)==0)need=need+1
	enddo
	do while(need>0) 
		give = ceiling_safe(need/dble(2**Maxlevel-need))
		do ii=1,2**Maxlevel
			nn = msh%pretree(ii)
			if(nn>1)then
				msh%pretree(ii) = msh%pretree(ii) - min(min(nn-1,give),need)
				need = need - min(min(nn-1,give),need)
			endif
		enddo
	enddo
	do ii=1,2**Maxlevel
		if(msh%pretree(ii)==0)msh%pretree(ii)=1
	enddo
	! write(*,*)'after adjustment:',msh%pretree
	tree(1:2**Maxlevel) = msh%pretree(1:2**Maxlevel) 
	
	
	!**** the geometry points are provided by user 
	if(option%nogeo==0)then 
		if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "User-supplied kernel requiring reorder:"
		Dimn = Ndim 
		allocate (msh%xyz(Dimn,0:msh%Nunk))
		ii=0
		do edge=1,msh%Nunk
			msh%xyz(1:Dimn,edge)= Locations(ii+1:ii+Dimn)
			ii = ii + Dimn
		enddo 
	endif
	
	
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "    "
	t2 = OMP_get_wtime()
	

	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "HODLR formatting......"
    call HODLR_structuring(ho_bf,option,msh,ker,element_Zmn_user_C,ptree)
	call BPlus_structuring(ho_bf,option,msh,ptree)
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "HODLR formatting finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "    "
	t2 = OMP_get_wtime()
	
	
    !call compression_test()
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "HODLR construction......"
    call HODLR_construction(ho_bf,option,stats,msh,ker,element_Zmn_user_C,ptree)
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "HODLR construction finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	! write(*,*)t2-t1

	
	!**** return the permutation vector 
	msh%idxs = ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
	msh%idxe = ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)		
	Npo_loc = msh%idxe-msh%idxs+1		
	if (ptree%MyID==Main_ID) then	
		do edge=1,Npo
			Permutation(edge) = msh%new2old(edge)
		enddo
	endif	
	
	!**** return the C address of hodlr structures to C caller
	ho_bf_Cptr=c_loc(ho_bf)
	option_Cptr=c_loc(option)
	stats_Cptr=c_loc(stats)
	msh_Cptr=c_loc(msh)
	ker_Cptr=c_loc(ker)
	ptree_Cptr=c_loc(ptree)


	
end subroutine C_HODLR_Construct






!**** C interface of HODLR construction via blackbox matvec
	!N: matrix size (in)
	!nlevel: the number of top levels that have been ordered (in)
	!tree: the order tree provided by the caller, if incomplete, the init routine will make it complete (inout)
	!Permutation: return the permutation vector new2old (indexed from 1) (out)
	!N_loc: number of local row/column indices (out)    
	!ho_bf_Cptr: the structure containing HODLR (out)        
	!option_Cptr: the structure containing option (in)        
	!stats_Cptr: the structure containing statistics (inout)         
	!msh_Cptr: the structure containing points and ordering information (out)        
	!ker_Cptr: the structure containing kernel quantities (out)
	!ptree_Cptr: the structure containing process tree (in)
subroutine C_HODLR_Construct_Matvec_Init(N,nlevel,tree,Permutation,N_loc,ho_bf_Cptr,option_Cptr,stats_Cptr,msh_Cptr,ker_Cptr,ptree_Cptr) bind(c, name="c_hodlr_construct_matvec_init")	
	implicit none 
	integer N
	
    real(kind=8) para
    real(kind=8) tolerance,h,lam
    integer Primary_block, nn, mm, MyID_old,Maxlevel,give,need
    integer i,j,k,ii,edge, threads_num,nth,Dimn,nmpi,ninc, acam
	real(kind=8),parameter :: cd = 299792458d0
	integer,allocatable:: groupmembers(:)
	integer nlevel,level
	integer Permutation(*),tree(*)
	integer N_loc
	! type(matricesblock), pointer :: blocks_i
	integer groupm
	type(c_ptr) :: ho_bf_Cptr
	type(c_ptr) :: option_Cptr
	type(c_ptr) :: stats_Cptr
	type(c_ptr) :: msh_Cptr
	type(c_ptr) :: ker_Cptr
	type(c_ptr) :: ptree_Cptr
	! type(c_ptr), intent(in),target :: C_QuantApp
	! type(c_funptr), intent(in),value,target :: C_FuncHMatVec
	
	type(Hoption),pointer::option	
	type(Hstat),pointer::stats
	type(mesh),pointer::msh
	type(kernelquant),pointer::ker
	type(hobf),pointer::ho_bf
	type(proctree),pointer::ptree	
	integer seed_myid(50)
	integer times(8)	
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	real(kind=8):: Memory=0d0,error
	character(len=1024)  :: strings
	integer N_unk_loc
	
	!**** allocate HODLR solver structures 
	allocate(ho_bf)
	! allocate(option)
	! allocate(stats)
	allocate(msh)
	allocate(ker)

	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	


	
	if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*)'NUMBER_MPI=',ptree%nproc
 	
 	threads_num=1
    CALL getenv("OMP_NUM_THREADS", strings)
	strings = TRIM(strings)	
	if(LEN_TRIM(strings)>0)then
		read(strings , *) threads_num
	endif
	if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*)'OMP_NUM_THREADS=',threads_num
	call OMP_set_num_threads(threads_num)
	
	!**** create a random seed		
	call DATE_AND_TIME(values=times)     ! Get the current time 
	seed_myid(1) = times(4) * (360000*times(5) + 6000*times(6) + 100*times(7) + times(8))
	! seed_myid(1) = myid*1000
	! call RANDOM_SEED(PUT=seed_myid)
	call init_random_seed()
	
	if(ptree%MyID==Main_ID .and. option%verbosity>0)then
    write(*,*) "HODLR_BUTTERFLY_SOLVER"
    write(*,*) "   "
	endif
	
	! !**** register the user-defined function and type in ker 
	! ker%C_QuantApp => C_QuantApp
	! ker%C_FuncHMatVec => C_FuncHMatVec	
	

	msh%Nunk = N
		
		
	t1 = OMP_get_wtime()

	
	if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "User-supplied kernel:"
	Maxlevel=nlevel
	allocate(msh%pretree(2**Maxlevel))
	
	msh%pretree(1:2**Maxlevel) = tree(1:2**Maxlevel)	

	
	!**** make 0-element node a 1-element node 
 	
	! write(*,*)'before adjustment:',msh%pretree
	need = 0
	do ii=1,2**Maxlevel
		if(msh%pretree(ii)==0)need=need+1
	enddo
	do while(need>0) 
		give = ceiling_safe(need/dble(2**Maxlevel-need))
		do ii=1,2**Maxlevel
			nn = msh%pretree(ii)
			if(nn>1)then
				msh%pretree(ii) = msh%pretree(ii) - min(min(nn-1,give),need)
				need = need - min(min(nn-1,give),need)
			endif
		enddo
	enddo
	do ii=1,2**Maxlevel
		if(msh%pretree(ii)==0)msh%pretree(ii)=1
	enddo
	! write(*,*)'after adjustment:',msh%pretree
	tree(1:2**Maxlevel) = msh%pretree(1:2**Maxlevel) 
	
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "    "
	t2 = OMP_get_wtime()
	

	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "HODLR formatting......"
    call HODLR_structuring(ho_bf,option,msh,ker,element_Zmn_user_C,ptree)
	call BPlus_structuring(ho_bf,option,msh,ptree)
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "HODLR formatting finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "    "
	t2 = OMP_get_wtime()
	
	!**** return the permutation vector 
	msh%idxs = ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
	msh%idxe = ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)		
	N_loc = msh%idxe-msh%idxs+1		
	if (ptree%MyID==Main_ID) then	
		do edge=1,N
			Permutation(edge) = msh%new2old(edge)
		enddo
	endif	
	
	!**** return the C address of hodlr structures to C caller
	ho_bf_Cptr=c_loc(ho_bf)
	option_Cptr=c_loc(option)
	stats_Cptr=c_loc(stats)
	msh_Cptr=c_loc(msh)
	ker_Cptr=c_loc(ker)
	ptree_Cptr=c_loc(ptree)
	
end subroutine C_HODLR_Construct_Matvec_Init



!**** C interface of HODLR construction via blackbox matvec
	!ho_bf_Cptr: the structure containing HODLR (inout)        
	!option_Cptr: the structure containing option (in)        
	!stats_Cptr: the structure containing statistics (inout)        
	!msh_Cptr: the structure containing points and ordering information (in)        
	!ker_Cptr: the structure containing kernel quantities (inout)
	!ptree_Cptr: the structure containing process tree (in)
	!C_FuncHMatVec: the C_pointer to user-provided function to multiply A and A* with vectors (in)
	!C_QuantApp: the C_pointer to user-defined quantities required to for entry evaluation and sampling (in)
subroutine C_HODLR_Construct_Matvec_Compute(ho_bf_Cptr,option_Cptr,stats_Cptr,msh_Cptr,ker_Cptr,ptree_Cptr,C_FuncHMatVec,C_QuantApp) bind(c, name="c_hodlr_construct_matvec_compute")	
	implicit none 

    real(kind=8) para
    real(kind=8) tolerance,h,lam
    integer Primary_block, nn, mm, MyID_old,Maxlevel,give,need
    integer i,j,k,ii,edge, threads_num,nth,Dimn,nmpi,ninc, acam
	real(kind=8),parameter :: cd = 299792458d0
	integer,allocatable:: groupmembers(:)
	integer level
	! type(matricesblock), pointer :: blocks_i
	integer groupm
	type(c_ptr) :: ho_bf_Cptr
	type(c_ptr) :: option_Cptr
	type(c_ptr) :: stats_Cptr
	type(c_ptr) :: msh_Cptr
	type(c_ptr) :: ker_Cptr
	type(c_ptr) :: ptree_Cptr
	type(c_ptr), intent(in),target :: C_QuantApp
	type(c_funptr), intent(in),value,target :: C_FuncHMatVec
	
	type(Hoption),pointer::option	
	type(Hstat),pointer::stats
	type(mesh),pointer::msh
	type(kernelquant),pointer::ker
	type(hobf),pointer::ho_bf
	type(proctree),pointer::ptree	
	integer seed_myid(50)
	integer times(8)	
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	real(kind=8):: Memory=0d0,error
	character(len=1024)  :: strings
	integer N_unk_loc
	
	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	call c_f_pointer(ho_bf_Cptr, ho_bf)
	call c_f_pointer(msh_Cptr, msh)
	call c_f_pointer(ker_Cptr, ker)	
	
	

	!**** register the user-defined function and type in ker 
	ker%C_QuantApp => C_QuantApp
	ker%C_FuncHMatVec => C_FuncHMatVec	
	

	t1 = OMP_get_wtime()	
     if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "FastMATVEC-based HODLR construction......"
	N_unk_loc = msh%idxe-msh%idxs+1
	call HODLR_randomized(ho_bf,matvec_user_C,N_unk_loc,Memory,error,option,stats,ker,ptree,msh)
	
 if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "FastMATVEC-based HODLR construction finished"
 if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	! write(*,*)t2-t1
	
	!**** return the C address of hodlr structures to C caller
	ho_bf_Cptr=c_loc(ho_bf)
	option_Cptr=c_loc(option)
	stats_Cptr=c_loc(stats)
	msh_Cptr=c_loc(msh)
	ker_Cptr=c_loc(ker)
	ptree_Cptr=c_loc(ptree)
	
end subroutine C_HODLR_Construct_Matvec_Compute


!**** C interface of BF construction via blackbox matvec
	!M,N: matrix size (in)
	!M_loc,N_loc: number of local row/column indices (out)    
	!bf_Cptr: the structure containing the block (out)        
	!option_Cptr: the structure containing option (in)        
	!stats_Cptr: the structure containing statistics (inout)        
	!msh_Cptr: the structure containing points and ordering information combined from mshr_Cptr and mshc_Cptr (out)      
	!mshr_Cptr: the structure containing points and ordering information for the row dimension (in)      
	!mshc_Cptr: the structure containing points and ordering information for the column dimension (in)      
	!ker_Cptr: the structure containing kernel quantities (out)
	!ptree_Cptr: the structure containing process tree (in)
subroutine C_BF_Construct_Matvec_Init(M,N,M_loc,N_loc,mshr_Cptr,mshc_Cptr,bf_Cptr,option_Cptr,stats_Cptr,msh_Cptr,ker_Cptr,ptree_Cptr) bind(c, name="c_bf_construct_matvec_init")	
	implicit none 
	integer M,N

    integer Maxlevel
    integer i,j,k,ii,edge, threads_num,nth,Dimn,nmpi,ninc, acam
	integer,allocatable:: groupmembers(:)
	integer level
	integer M_loc,N_loc
	type(c_ptr) :: bf_Cptr
	type(c_ptr) :: option_Cptr
	type(c_ptr) :: stats_Cptr
	type(c_ptr) :: msh_Cptr,mshr_Cptr,mshc_Cptr
	type(c_ptr) :: ker_Cptr
	type(c_ptr) :: ptree_Cptr

	type(Hoption),pointer::option	
	type(Hstat),pointer::stats
	type(mesh),pointer::msh,mshr,mshc
	type(kernelquant),pointer::ker
	type(matrixblock),pointer::blocks
	type(proctree),pointer::ptree	
	integer seed_myid(50)
	integer times(8)	
	real(kind=8) t1,t2
	character(len=1024)  :: strings
	
	!**** allocate HODLR solver structures 
	allocate(blocks)
	allocate(msh)
	allocate(ker)

	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	call c_f_pointer(mshr_Cptr, mshr)
	call c_f_pointer(mshc_Cptr, mshc)
	
	if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*)'NUMBER_MPI=',ptree%nproc
 	
 	threads_num=1
    CALL getenv("OMP_NUM_THREADS", strings)
	strings = TRIM(strings)	
	if(LEN_TRIM(strings)>0)then
		read(strings , *) threads_num
	endif
	if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*)'OMP_NUM_THREADS=',threads_num
	call OMP_set_num_threads(threads_num)
	
	!**** create a random seed		
	call DATE_AND_TIME(values=times)     ! Get the current time 
	seed_myid(1) = times(4) * (360000*times(5) + 6000*times(6) + 100*times(7) + times(8))
	! seed_myid(1) = myid*1000
	! call RANDOM_SEED(PUT=seed_myid)
	call init_random_seed()
	
	! if(ptree%MyID==Main_ID)then
    ! write(*,*) "HODLR_BUTTERFLY_SOLVER"
    ! write(*,*) "   "
	! endif
	
	! !**** register the user-defined function and type in ker 
	! ker%C_QuantApp => C_QuantApp
	! ker%C_FuncHMatVec => C_FuncHMatVec	
	
	call assert(mshr%Nunk==M,'mshr%Nunk\=M')
	call assert(mshc%Nunk==N,'mshc%Nunk\=N')
	call assert(mshc%Maxgroup==mshr%Maxgroup,'mshc%Maxgroup\=mshr%Maxgroup')
	
	msh%Nunk = N+M
	msh%Maxgroup=mshc%Maxgroup*2+1
	allocate (msh%basis_group(msh%Maxgroup))
	msh%basis_group(1)%head=1
	msh%basis_group(1)%tail=N+M
	msh%basis_group(1)%pgno=1
	call copy_basis_group(mshr%basis_group,1,mshr%Maxgroup,msh%basis_group,2,msh%Maxgroup,0)
	call copy_basis_group(mshc%basis_group,1,mshc%Maxgroup,msh%basis_group,3,msh%Maxgroup,mshr%Nunk)
	
	blocks%level=1
	blocks%col_group=3
	blocks%row_group=2
	blocks%pgno=1			
	blocks%headm=msh%basis_group(blocks%row_group)%head
	blocks%headn=msh%basis_group(blocks%col_group)%head
	blocks%M=M
	blocks%N=N
	blocks%style = 2
	
	
	Maxlevel = ceiling_safe(log(dble(msh%Maxgroup)) / log(2d0))
	if(blocks%level>option%LRlevel)then
		blocks%level_butterfly = 0 ! low rank below LRlevel
	else 				
		blocks%level_butterfly = Maxlevel - blocks%level   ! butterfly 
	endif	
	
	call ComputeParallelIndices(Maxlevel,blocks,blocks%pgno,ptree,msh,0)
		
		
		
	M_loc = blocks%M_loc	
	N_loc = blocks%N_loc	
		

	!**** return the C address of hodlr structures to C caller
	bf_Cptr=c_loc(blocks)
	option_Cptr=c_loc(option)
	stats_Cptr=c_loc(stats)
	msh_Cptr=c_loc(msh)
	ker_Cptr=c_loc(ker)
	ptree_Cptr=c_loc(ptree)
	
end subroutine C_BF_Construct_Matvec_Init





!**** C interface of BF construction via blackbox matvec
	!bf_Cptr: the structure containing the block (inout)        
	!option_Cptr: the structure containing option (in)        
	!stats_Cptr: the structure containing statistics (inout)        
	!msh_Cptr: the structure containing points and ordering information (in)    
	!ker_Cptr: the structure containing kernel quantities (inout)
	!ptree_Cptr: the structure containing process tree (in)
	!C_FuncBMatVec: the C_pointer to user-provided function to multiply A and A* with vectors (in)
	!C_QuantApp: the C_pointer to user-defined quantities required to for entry evaluation and sampling (in)	
subroutine C_BF_Construct_Matvec_Compute(bf_Cptr,option_Cptr,stats_Cptr,msh_Cptr,ker_Cptr,ptree_Cptr,C_FuncBMatVec,C_QuantApp) bind(c, name="c_bf_construct_matvec_compute")	
	implicit none 

    integer Maxlevel
    integer i,j,k,ii,edge, threads_num,nth,Dimn,nmpi,ninc, acam
	integer,allocatable:: groupmembers(:)
	integer level
	type(c_ptr) :: bf_Cptr
	type(c_ptr) :: option_Cptr
	type(c_ptr) :: stats_Cptr
	type(c_ptr) :: msh_Cptr
	type(c_ptr) :: ker_Cptr
	type(c_ptr) :: ptree_Cptr
	type(c_ptr), intent(in),target :: C_QuantApp
	type(c_funptr), intent(in),value,target :: C_FuncBMatVec

	type(Hoption),pointer::option	
	type(Hstat),pointer::stats
	type(mesh),pointer::msh
	type(kernelquant),pointer::ker
	type(matrixblock),pointer::blocks
	type(proctree),pointer::ptree	
	integer seed_myid(50)
	integer times(8)	
	real(kind=8) t1,t2,error
	
	!**** allocate HODLR solver structures 


	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	call c_f_pointer(msh_Cptr, msh)
	call c_f_pointer(bf_Cptr, blocks)
	call c_f_pointer(ker_Cptr, ker)
	
	
	! !**** register the user-defined function and type in ker 
	ker%C_QuantApp => C_QuantApp
	ker%C_FuncBMatVec => C_FuncBMatVec	
	
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "FastMATVEC-based BF construction......"	
	
	call BF_randomized(blocks%level_butterfly,option%rank0,option%rankrate,blocks,ker,Bmatvec_user_C,error,'CMatVec',option,stats,ptree,msh) 

	if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "FastMATVEC-based BF construction finished"
	
	!**** return the C address of hodlr structures to C caller
	bf_Cptr=c_loc(blocks)
	option_Cptr=c_loc(option)
	stats_Cptr=c_loc(stats)
	msh_Cptr=c_loc(msh)
	ker_Cptr=c_loc(ker)
	ptree_Cptr=c_loc(ptree)
	
end subroutine C_BF_Construct_Matvec_Compute



!**** C interface of HODLR factorization
	!ho_bf_for_Cptr: the structure containing HODLR                  
	!option_Cptr: the structure containing option         
	!stats_Cptr: the structure containing statistics         
	!ptree_Cptr: the structure containing process tree
subroutine C_HODLR_Factor(ho_bf_for_Cptr,option_Cptr,stats_Cptr,ptree_Cptr,msh_Cptr) bind(c, name="c_hodlr_factor")	
	implicit none 

	type(c_ptr), intent(inout) :: ho_bf_for_Cptr
	type(c_ptr), intent(in) :: ptree_Cptr	
	type(c_ptr), intent(in) :: msh_Cptr	
	type(c_ptr) :: option_Cptr
	type(c_ptr), intent(inout) :: stats_Cptr
	! type(c_ptr), intent(out) :: ho_bf_inv_Cptr
	

	type(Hoption),pointer::option	
	type(Hstat),pointer::stats
	type(hobf),pointer::ho_bf1
	! type(hobf),pointer::ho_bf_inv
	type(proctree),pointer::ptree	
	type(mesh),pointer::msh	

	! real(kind=8):: tol_fact
	
	call c_f_pointer(ho_bf_for_Cptr, ho_bf1)
	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	call c_f_pointer(msh_Cptr, msh)
	
	! allocate(ho_bf_tmp)
	! call copy_HOBF(ho_bf1,ho_bf_tmp)	! currently this subroutine only copies forward components 
	! ho_bf_inv=>ho_bf1
	! ho_bf1=>ho_bf_tmp
	
	
	if(option%precon/=NOPRECON)then
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "Cascading factorizing......"
    call HODLR_Factorization(ho_bf1,option,stats,ptree,msh)
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "Cascading factorizing finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "    "	
	end if
	
	! return the C address of hodlr structures to C caller
	! ho_bf_inv_Cptr=c_loc(ho_bf_inv)	
	ho_bf_for_Cptr=c_loc(ho_bf1)	

end subroutine C_HODLR_Factor


!**** C interface of HODLR solve
	!x: local solution vector         
	!b: local RHS        
	!Nloc: size of local RHS     
	!Nrhs: number of RHSs     
	!ho_bf_for_Cptr: the structure containing HODLR         
	!option_Cptr: the structure containing option         
	!stats_Cptr: the structure containing statistics         
	!ptree_Cptr: the structure containing process tree
subroutine C_HODLR_Solve(x,b,Nloc,Nrhs,ho_bf_for_Cptr,option_Cptr,stats_Cptr,ptree_Cptr) bind(c, name="c_hodlr_solve")	
	implicit none 

	integer Nloc,Nrhs
	DT::x(Nloc,Nrhs),b(Nloc,Nrhs)
	
	type(c_ptr), intent(in) :: ho_bf_for_Cptr
	type(c_ptr), intent(in) :: ptree_Cptr	
	type(c_ptr), intent(in) :: option_Cptr
	type(c_ptr), intent(in) :: stats_Cptr
	! type(c_ptr), intent(in) :: ho_bf_inv_Cptr
	

	type(Hoption),pointer::option	
	type(Hstat),pointer::stats
	type(hobf),pointer::ho_bf1
	! type(hobf),pointer::ho_bf_inv
	type(proctree),pointer::ptree	


	call c_f_pointer(ho_bf_for_Cptr, ho_bf1)
	! call c_f_pointer(ho_bf_inv_Cptr, ho_bf_inv)
	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	
	if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "Solve ......"
	
	if(option%ErrSol==1)then
		call HODLR_Test_Solve_error(ho_bf1,option,ptree,stats)
	endif		
	
	call HODLR_Solution(ho_bf1,x,b,Nloc,Nrhs,option,ptree,stats)

	if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "Solve finished"
	if(ptree%MyID==Main_ID .and. option%verbosity>0)write(*,*) "    "	
	
	! return the C address of hodlr structures to C caller
	
	! ho_bf_inv_Cptr=c_loc(ho_bf_inv)	
	! ho_bf_for_Cptr=c_loc(ho_bf_for)	

end subroutine C_HODLR_Solve


	

!**** C interface of butterfly-vector multiplication
	!xin: input vector  
    !Ninloc:size of local input vectors	
	!xout: output vector        
	!Noutloc:size of local output vectors	
	!Ncol: number of vectors     
	!bf_for_Cptr: the structure containing butterfly                
	!stats_Cptr: the structure containing statistics         
	!ptree_Cptr: the structure containing process tree
subroutine C_BF_Mult(trans,xin,xout,Ninloc,Noutloc,Ncol,bf_for_Cptr,stats_Cptr,ptree_Cptr,a,b) bind(c, name="c_bf_mult")	
	implicit none 
	real(kind=8) t1,t2
	integer Ninloc,Noutloc,Ncol
	DT::xin(Ninloc,Ncol),xout(Noutloc,Ncol),a,b
	
	character(kind=c_char,len=1) :: trans(*)
	type(c_ptr), intent(in) :: bf_for_Cptr
	type(c_ptr), intent(in) :: ptree_Cptr	
	type(c_ptr), intent(in) :: stats_Cptr
	! type(c_ptr), intent(in) :: ho_bf_inv_Cptr
	
	integer strlen
	character(len=:),allocatable :: str
	type(Hstat),pointer::stats
	type(matrixblock),pointer::blocks
	! type(hobf),pointer::ho_bf_inv
	type(proctree),pointer::ptree	

	t1 = OMP_get_wtime()

	strlen=1
	do while(trans(strlen) /= c_null_char)
		strlen = strlen + 1
	enddo
	strlen = strlen -1
	allocate(character(len=strlen) :: str)
	str = transfer(trans(1:strlen), str)	
	
	call c_f_pointer(bf_for_Cptr, blocks)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	
	if(trim(str)=='N')then
		call BF_block_MVP_dat(blocks,trim(str),Noutloc,Ninloc,Ncol,xin,xout,a,b,ptree,stats)	
	else
		call BF_block_MVP_dat(blocks,trim(str),Ninloc,Noutloc,Ncol,xin,xout,a,b,ptree,stats)
	endif
	
	t2 = OMP_get_wtime() 
	
	stats%Time_C_Mult = stats%Time_C_Mult + t2-t1
	stats%Flop_C_Mult = stats%Flop_C_Mult + stats%Flop_Tmp
	
	! write(*,*)t2-t1	
	deallocate(str)
end subroutine C_BF_Mult	
	
	
	
	

!**** C interface of HODLR-vector multiplication
	!xin: input vector  
    !Ninloc:size of local input vectors	
	!xout: output vector        
	!Noutloc:size of local output vectors	
	!Ncol: number of vectors     
	!ho_bf_for_Cptr: the structure containing HODLR         
	!option_Cptr: the structure containing option         
	!stats_Cptr: the structure containing statistics         
	!ptree_Cptr: the structure containing process tree
subroutine C_HODLR_Mult(trans,xin,xout,Ninloc,Noutloc,Ncol,ho_bf_for_Cptr,option_Cptr,stats_Cptr,ptree_Cptr) bind(c, name="c_hodlr_mult")	
	implicit none 
	real(kind=8) t1,t2
	integer Ninloc,Noutloc,Ncol
	DT::xin(Ninloc,Ncol),xout(Noutloc,Ncol)
	
	character(kind=c_char,len=1) :: trans(*)
	type(c_ptr), intent(in) :: ho_bf_for_Cptr
	type(c_ptr), intent(in) :: ptree_Cptr	
	type(c_ptr), intent(in) :: option_Cptr
	type(c_ptr), intent(in) :: stats_Cptr
	! type(c_ptr), intent(in) :: ho_bf_inv_Cptr
	
	integer strlen
	character(len=:),allocatable :: str
	type(Hoption),pointer::option	
	type(Hstat),pointer::stats
	type(hobf),pointer::ho_bf1
	! type(hobf),pointer::ho_bf_inv
	type(proctree),pointer::ptree	

	t1 = OMP_get_wtime()

	strlen=1
	do while(trans(strlen) /= c_null_char)
		strlen = strlen + 1
	enddo
	strlen = strlen -1
	allocate(character(len=strlen) :: str)
	str = transfer(trans(1:strlen), str)	
	
	call c_f_pointer(ho_bf_for_Cptr, ho_bf1)
	! call c_f_pointer(ho_bf_inv_Cptr, ho_bf_inv)
	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	
    ! if(ptree%MyID==Main_ID)write(*,*) "Multiply ......"
	
	call assert(Noutloc==Ninloc,"not square Z")
	call MVM_Z_forward(trim(str),Noutloc,Ncol,1,ho_bf1%Maxlevel+1,xin,xout,ho_bf1,ptree,stats)
	! need to use another Flop counter for this operation in future

    ! if(ptree%MyID==Main_ID)write(*,*) "Multiply finished"
    ! if(ptree%MyID==Main_ID)write(*,*) "    "	
	
	t2 = OMP_get_wtime() 
	
	stats%Time_C_Mult = stats%Time_C_Mult + t2-t1
	stats%Flop_C_Mult = stats%Flop_C_Mult + stats%Flop_Tmp
	
	! write(*,*)t2-t1	
	deallocate(str)
end subroutine C_HODLR_Mult




!**** C interface of HODLR(inverse)-vector multiplication
	!xin: input vector  
    !Ninloc:size of local input vectors	
	!xout: output vector        
	!Noutloc:size of local output vectors	
	!Ncol: number of vectors     
	!ho_bf_for_Cptr: the structure containing HODLR         
	!option_Cptr: the structure containing option         
	!stats_Cptr: the structure containing statistics         
	!ptree_Cptr: the structure containing process tree
subroutine C_HODLR_Inv_Mult(trans,xin,xout,Ninloc,Noutloc,Ncol,ho_bf_for_Cptr,option_Cptr,stats_Cptr,ptree_Cptr) bind(c, name="c_hodlr_inv_mult")	
	implicit none 
	real(kind=8) t1,t2
	integer Ninloc,Noutloc,Ncol
	DT::xin(Ninloc,Ncol),xout(Noutloc,Ncol)
	
	character(kind=c_char,len=1) :: trans(*)
	type(c_ptr), intent(in) :: ho_bf_for_Cptr
	type(c_ptr), intent(in) :: ptree_Cptr	
	type(c_ptr), intent(in) :: option_Cptr
	type(c_ptr), intent(in) :: stats_Cptr
	! type(c_ptr), intent(in) :: ho_bf_inv_Cptr
	
	integer strlen
	character(len=:),allocatable :: str
	type(Hoption),pointer::option	
	type(Hstat),pointer::stats
	type(hobf),pointer::ho_bf1
	! type(hobf),pointer::ho_bf_inv
	type(proctree),pointer::ptree	

	t1 = OMP_get_wtime()

	strlen=1
	do while(trans(strlen) /= c_null_char)
		strlen = strlen + 1
	enddo
	strlen = strlen -1
	allocate(character(len=strlen) :: str)
	str = transfer(trans(1:strlen), str)	
	
	call c_f_pointer(ho_bf_for_Cptr, ho_bf1)
	! call c_f_pointer(ho_bf_inv_Cptr, ho_bf_inv)
	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	
    ! if(ptree%MyID==Main_ID)write(*,*) "Multiply ......"
	
	call assert(Noutloc==Ninloc,"not square Z")
	call MVM_Z_factorized(trim(str),Noutloc,Ncol,xin,xout,ho_bf1,ptree,stats)
	! need to use another Flop counter for this operation in future

    ! if(ptree%MyID==Main_ID)write(*,*) "Multiply finished"
    ! if(ptree%MyID==Main_ID)write(*,*) "    "	
	
	t2 = OMP_get_wtime() 
	
	stats%Time_C_Mult = stats%Time_C_Mult + t2-t1
	stats%Flop_C_Mult = stats%Flop_C_Mult + stats%Flop_Tmp
	
	! write(*,*)t2-t1	
	deallocate(str)
end subroutine C_HODLR_Inv_Mult



!**** C interface of deleting statistics
	!stats_Cptr: the structure containing statistics
subroutine C_HODLR_Deletestats(stats_Cptr) bind(c, name="c_hodlr_deletestats")	
	implicit none 
	type(c_ptr), intent(inout) :: stats_Cptr
	type(Hstat),pointer::stats

	call c_f_pointer(stats_Cptr, stats)
	call delete_Hstat(stats)
	deallocate(stats)
	stats_Cptr=c_null_ptr
	
end subroutine C_HODLR_Deletestats


!**** C interface of deleting process tree
	!ptree_Cptr: the structure containing process tree
subroutine C_HODLR_Deleteproctree(ptree_Cptr) bind(c, name="c_hodlr_deleteproctree")	
	implicit none 
	type(c_ptr), intent(inout) :: ptree_Cptr
	type(proctree),pointer::ptree

	call c_f_pointer(ptree_Cptr, ptree)
	call delete_proctree(ptree)
	deallocate(ptree)
	ptree_Cptr=c_null_ptr
	
end subroutine C_HODLR_Deleteproctree


!**** C interface of deleting mesh 
	!msh_Cptr: the structure containing mesh
subroutine C_HODLR_Deletemesh(msh_Cptr) bind(c, name="c_hodlr_deletemesh")	
	implicit none 
	type(c_ptr), intent(inout) :: msh_Cptr
	type(mesh),pointer::msh

	call c_f_pointer(msh_Cptr, msh)
	call delete_mesh(msh)
	deallocate(msh)
	msh_Cptr=c_null_ptr
	
end subroutine C_HODLR_Deletemesh


!**** C interface of deleting kernelquant 
	!ker_Cptr: the structure containing kernelquant
subroutine C_HODLR_Deletekernelquant(ker_Cptr) bind(c, name="c_hodlr_deletekernelquant")	
	implicit none 
	type(c_ptr), intent(inout) :: ker_Cptr
	type(kernelquant),pointer::ker

	call c_f_pointer(ker_Cptr, ker)
	call delete_kernelquant(ker)
	deallocate(ker)
	ker_Cptr=c_null_ptr
	
end subroutine C_HODLR_Deletekernelquant


!**** C interface of deleting HOBF 
	!ho_bf_Cptr: the structure containing HOBF
subroutine C_HODLR_DeleteHOBF(ho_bf_Cptr) bind(c, name="c_hodlr_deletehobf")	
	implicit none 
	type(c_ptr), intent(inout) :: ho_bf_Cptr
	type(hobf),pointer::ho_bf

	call c_f_pointer(ho_bf_Cptr, ho_bf)
	call delete_HOBF(ho_bf)
	deallocate(ho_bf)
	ho_bf_Cptr=c_null_ptr
	
end subroutine C_HODLR_DeleteHOBF


!**** C interface of deleting a BF 
	!bf_Cptr: the structure containing BF
subroutine C_BF_DeleteBF(bf_Cptr) bind(c, name="c_bf_deletebf")	
	implicit none 
	type(c_ptr), intent(inout) :: bf_Cptr
	type(matrixblock),pointer::blocks

	call c_f_pointer(bf_Cptr, blocks)
	call delete_blocks(blocks,1)
	deallocate(blocks)
	bf_Cptr=c_null_ptr
end subroutine C_BF_DeleteBF



!**** C interface of deleting Hoption 
	!option_Cptr: the structure containing Hoption
subroutine C_HODLR_Deleteoption(option_Cptr) bind(c, name="c_hodlr_deleteoption")	
	implicit none 
	type(c_ptr), intent(inout) :: option_Cptr
	type(Hoption),pointer::option

	call c_f_pointer(option_Cptr, option)
	deallocate(option)
	option_Cptr=c_null_ptr
	
end subroutine C_HODLR_Deleteoption


end module HODLR_wrapper
