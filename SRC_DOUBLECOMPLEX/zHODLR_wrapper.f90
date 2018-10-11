#include "HODLR_config.fi"
module z_HODLR_wrapper
    use z_HODLR_DEFS
	
	use z_HODLR_structure
	use z_HODLR_factor
	use z_HODLR_constr
	use omp_lib
	use z_misc
	use z_HODLR_Solve_Mul
    use iso_c_binding    
	

contains 



subroutine element_Zmn_user_C(edge_m,edge_n,value_e,msh,ker)
    use z_HODLR_DEFS
    implicit none
    integer edge_m, edge_n
    DT value_e
	type(mesh)::msh
	type(kernelquant)::ker	
	procedure(C_Z_elem), POINTER :: proc
	
	value_e=0
	call c_f_procpointer(ker%C_FuncZmn, proc)
	call proc(msh%new2old(edge_m)-1,msh%new2old(edge_n)-1,value_e,ker%C_QuantZmn)
	return
    
end subroutine element_Zmn_user_C


!**** C interface of process tree construction
	!nmpi: number of MPIs for one hodlr
	!MPIcomm: MPI communicator from C caller
	!groupmembers: MPI ranks in MPIcomm for one hodlr
	!ptree_Cptr: the structure containing process tree
subroutine C_HODLR_Createptree(nmpi,groupmembers,MPIcomm,ptree_Cptr) bind(c, name="z_c_hodlr_createptree")	
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
subroutine C_HODLR_Createstats(stats_Cptr) bind(c, name="z_c_hodlr_createstats")	
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
subroutine C_HODLR_Printstats(stats_Cptr,ptree_Cptr) bind(c, name="z_c_hodlr_printstats")	
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
subroutine C_HODLR_Createoption(option_Cptr) bind(c, name="z_c_hodlr_createoption")	
	implicit none 
	type(c_ptr) :: option_Cptr
	type(Hoption),pointer::option	
	
	allocate(option)
	!**** set default hodlr options 	
	call SetDefaultOptions(option)
	
	option_Cptr=c_loc(option)
	
end subroutine C_HODLR_Createoption



!**** C interface of set one entry in option
	!option_Cptr: the structure containing option       
subroutine C_HODLR_Setoption(option_Cptr,nam,val_Cptr) bind(c, name="z_c_hodlr_setoption")	
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
	!Permutation: return the permutation vector new2old
	!Npo_loc: number of local row/column indices    
	!ho_bf_Cptr: the structure containing HODLR         
	!option_Cptr: the structure containing option         
	!stats_Cptr: the structure containing statistics         
	!msh_Cptr: the structure containing points and ordering information         
	!ker_Cptr: the structure containing kernel quantities
	!ptree_Cptr: the structure containing process tree
	!C_FuncZmn: the C_pointer to user-provided function to sample mn^th entry of the matrix
	!C_QuantZmn: the C_pointer to user-defined quantities required to sample mn^th entry of the matrix
	!MPIcomm: user-provided MPI communicator
subroutine C_HODLR_Construct(Npo,Ndim,Locations,nlevel,tree,Permutation,Npo_loc,ho_bf_Cptr,option_Cptr,stats_Cptr,msh_Cptr,ker_Cptr,ptree_Cptr,C_FuncZmn,C_QuantZmn,MPIcomm) bind(c, name="z_c_hodlr_construct")	
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
	type(c_ptr), intent(in),target :: C_QuantZmn
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
	


	
	if(ptree%MyID==Main_ID)write(*,*)'NUMBER_MPI=',ptree%nproc
 	
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
	! call RANDOM_SEED(PUT=seed_myid)
	call init_random_seed()
	
	if(ptree%MyID==Main_ID)then
    write(*,*) "HODLR_BUTTERFLY_SOLVER"
    write(*,*) "   "
	endif
	
	!**** register the user-defined function and type in ker 
	ker%C_QuantZmn => C_QuantZmn
	ker%C_FuncZmn => C_FuncZmn	
	

	msh%Nunk = Npo
		
		
	t1 = OMP_get_wtime()

	
	if(ptree%MyID==Main_ID)write(*,*) "User-supplied kernel:"
	Maxlevel=nlevel
	allocate(msh%pretree(2**Maxlevel))
	
	msh%pretree(1:2**Maxlevel) = tree(1:2**Maxlevel)	

	write(*,*)'before adjustment:',msh%pretree
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
	write(*,*)'after adjustment:',msh%pretree
	
	
	
	
	!**** the geometry points are provided by user 
	if(option%nogeo==0)then 
		if(ptree%MyID==Main_ID)write(*,*) "User-supplied kernel requiring reorder:"
		Dimn = Ndim 
		allocate (msh%xyz(Dimn,0:msh%Nunk))
		ii=0
		do edge=1,msh%Nunk
			msh%xyz(1:Dimn,edge)= Locations(ii+1:ii+Dimn)
			ii = ii + Dimn
		enddo 
	endif
	
	
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	

	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "HODLR formatting......"
    call HODLR_structuring(ho_bf,option,msh,ker,element_Zmn_user_C,ptree)
	call BPlus_structuring(ho_bf,option,msh,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "HODLR formatting finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	
	
    !call compression_test()
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "HODLR construction......"
    call z_HODLR_construction(ho_bf,option,stats,msh,ker,element_Zmn_user_C,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "HODLR construction finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
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



!**** C interface of HODLR factorization
	!ho_bf_for_Cptr: the structure containing HODLR                  
	!option_Cptr: the structure containing option         
	!stats_Cptr: the structure containing statistics         
	!ptree_Cptr: the structure containing process tree
subroutine C_HODLR_Factor(ho_bf_for_Cptr,option_Cptr,stats_Cptr,ptree_Cptr) bind(c, name="z_c_hodlr_factor")	
	implicit none 

	type(c_ptr), intent(inout) :: ho_bf_for_Cptr
	type(c_ptr), intent(in) :: ptree_Cptr	
	type(c_ptr) :: option_Cptr
	type(c_ptr), intent(inout) :: stats_Cptr
	! type(c_ptr), intent(out) :: ho_bf_inv_Cptr
	

	type(Hoption),pointer::option	
	type(Hstat),pointer::stats
	type(hobf),pointer::ho_bf1
	! type(hobf),pointer::ho_bf_inv
	type(proctree),pointer::ptree	

	! real(kind=8):: tol_fact
	
	call c_f_pointer(ho_bf_for_Cptr, ho_bf1)
	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	
	! allocate(ho_bf_tmp)
	! call copy_HOBF(ho_bf1,ho_bf_tmp)	! currently this subroutine only copies forward components 
	! ho_bf_inv=>ho_bf1
	! ho_bf1=>ho_bf_tmp
	
	
	if(option%precon/=NOPRECON)then
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing......"
    call HODLR_Factorization(ho_bf1,option,stats,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
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
subroutine C_HODLR_Solve(x,b,Nloc,Nrhs,ho_bf_for_Cptr,option_Cptr,stats_Cptr,ptree_Cptr) bind(c, name="z_c_hodlr_solve")	
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
	
    if(ptree%MyID==Main_ID)write(*,*) "Solve ......"
	
	if(option%ErrSol==1)then
		call HODLR_Test_Solve_error(ho_bf1,option,ptree,stats)
	endif		
	
	call HODLR_Solution(ho_bf1,x,b,Nloc,Nrhs,option,ptree,stats)

    if(ptree%MyID==Main_ID)write(*,*) "Solve finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	
	! return the C address of hodlr structures to C caller
	
	! ho_bf_inv_Cptr=c_loc(ho_bf_inv)	
	! ho_bf_for_Cptr=c_loc(ho_bf_for)	

end subroutine C_HODLR_Solve



!**** C interface of HODLR-vector multiplication
	!xin: input vector         
	!xout: output vector        
	!Nloc: size of local vectors     
	!Ncol: number of vectors     
	!ho_bf_for_Cptr: the structure containing HODLR         
	!option_Cptr: the structure containing option         
	!stats_Cptr: the structure containing statistics         
	!ptree_Cptr: the structure containing process tree
subroutine C_HODLR_Mult(trans,xin,xout,Nloc,Ncol,ho_bf_for_Cptr,option_Cptr,stats_Cptr,ptree_Cptr) bind(c, name="z_c_hodlr_mult")	
	implicit none 
	real(kind=8) t1,t2
	integer Nloc,Ncol
	DT::xin(Nloc,Ncol),xout(Nloc,Ncol)
	
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
	
    if(ptree%MyID==Main_ID)write(*,*) "Multiply ......"
	
	
	call MVM_Z_forward(trim(str),Nloc,Ncol,1,ho_bf1%Maxlevel+1,xin,xout,ho_bf1,ptree,stats)
	! need to use another Flop counter for this operation in future

    if(ptree%MyID==Main_ID)write(*,*) "Multiply finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	
	t2 = OMP_get_wtime() 
	
	stats%Time_C_Mult = stats%Time_C_Mult + t2-t1
	stats%Flop_C_Mult = stats%Flop_C_Mult + stats%Flop_Tmp
	
	! write(*,*)t2-t1	
	deallocate(str)
end subroutine C_HODLR_Mult



end module z_HODLR_wrapper
