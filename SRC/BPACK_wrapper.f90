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

#include "ButterflyPACK_config.fi"
module BPACK_wrapper
use BPACK_DEFS
use BPACK_structure
use BPACK_factor
use BPACK_constr
use BPACK_randomMVP
use omp_lib
use misc
use BPACK_Solve_Mul
use iso_c_binding

contains



subroutine element_Zmn_user_C(edge_m,edge_n,value_e,msh,option,ker)
    use BPACK_DEFS
    implicit none
    integer edge_m, edge_n
    DT value_e
	type(mesh)::msh
	type(Hoption)::option
	type(kernelquant)::ker
	procedure(C_Zelem), POINTER :: proc

	value_e=0
	call c_f_procpointer(ker%C_FuncZmn, proc)
	call proc(msh%new2old(edge_m)-1,msh%new2old(edge_n)-1,value_e,ker%C_QuantApp)
	value_e =value_e*option%scale_factor
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
subroutine C_BPACK_Createptree(nmpi,groupmembers,MPIcomm,ptree_Cptr) bind(c, name="c_bpack_createptree")
	implicit none
	integer nmpi
	integer MPIcomm
	integer:: groupmembers(*)
	type(c_ptr):: ptree_Cptr
	type(proctree),pointer::ptree

	allocate(ptree)
	call CreatePtree(nmpi,groupmembers,MPIcomm,ptree)
	ptree_Cptr=c_loc(ptree)
end subroutine C_BPACK_Createptree


!**** C interface of initializing statistics
	!nmpi: number of MPIs for one hodlr
	!MPIcomm: MPI communicator from C caller
	!groupmembers: MPI ranks in MPIcomm for one hodlr
	!stats_Cptr: the structure containing statistics
subroutine C_BPACK_Createstats(stats_Cptr) bind(c, name="c_bpack_createstats")
	implicit none
	type(c_ptr), intent(out) :: stats_Cptr
	type(Hstat),pointer::stats

	allocate(stats)
	!**** initialize statistics variables
	call InitStat(stats)
	stats_Cptr=c_loc(stats)

end subroutine C_BPACK_Createstats





!**** C interface of getting one entry in stats
	!stats_Cptr: the structure containing stats
subroutine C_BPACK_Getstats(stats_Cptr,nam,val_d) bind(c, name="c_bpack_getstats")
	implicit none
	real(kind=8)::val_d
	character(kind=c_char,len=1) :: nam(*)
	type(c_ptr) :: stats_Cptr
	type(Hstat),pointer::stats
	! character::nam(:)
	! type(c_ptr),value :: val_Cptr
	! integer,pointer::val_i
	! real(kind=8),pointer::val_d
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

	call c_f_pointer(stats_Cptr, stats)


	if(trim(str)=='Time_Fill')then
		val_d = stats%Time_Fill
	endif
	if(trim(str)=='Time_Factor')then
		val_d = stats%Time_Factor
	endif
	if(trim(str)=='Time_Solve')then
		val_d = stats%Time_Sol
	endif
	if(trim(str)=='Time_Sblock')then
		val_d = stats%Time_Sblock
	endif
	if(trim(str)=='Time_Inv')then
		val_d = stats%Time_Inv
	endif
	if(trim(str)=='Time_SMW')then
		val_d = stats%Time_SMW
	endif
	if(trim(str)=='Time_RedistB')then
		val_d = stats%Time_RedistB
	endif
	if(trim(str)=='Time_RedistV')then
		val_d = stats%Time_RedistV
	endif
	if(trim(str)=='Time_C_Mult')then
		val_d = stats%Time_C_Mult
	endif
	if(trim(str)=='Time_Direct_LU')then
		val_d = stats%Time_Direct_LU
	endif
	if(trim(str)=='Time_Add_Multiply')then
		val_d = stats%Time_Add_Multiply
	endif
	if(trim(str)=='Time_Multiply')then
		val_d = stats%Time_Multiply
	endif
	if(trim(str)=='Time_XLUM')then
		val_d = stats%Time_XLUM
	endif
	if(trim(str)=='Time_Split')then
		val_d = stats%Time_Split
	endif
	if(trim(str)=='Time_Comm')then
		val_d = stats%Time_Comm
	endif
	if(trim(str)=='Time_Idle')then
		val_d = stats%Time_Idle
	endif


	if(trim(str)=='Flop_Fill')then
		val_d = stats%Flop_Fill
	endif
	if(trim(str)=='Flop_Factor')then
		val_d = stats%Flop_Factor
	endif
	if(trim(str)=='Flop_Solve')then
		val_d = stats%Flop_Sol
	endif
	if(trim(str)=='Flop_C_Mult')then
		val_d = stats%Flop_C_Mult
	endif


	if(trim(str)=='Mem_Factor')then
		val_d = stats%Mem_Factor
	endif
	if(trim(str)=='Mem_Fill')then
		val_d = stats%Mem_Fill
	endif
	if(trim(str)=='Mem_Sblock')then
		val_d = stats%Mem_Sblock
	endif
	if(trim(str)=='Mem_SMW')then
		val_d = stats%Mem_SMW
	endif
	if(trim(str)=='Mem_Direct_inv')then
		val_d = stats%Mem_Direct_inv
	endif
	if(trim(str)=='Mem_Direct_for')then
		val_d = stats%Mem_Direct_for
	endif
	if(trim(str)=='Mem_int_vec')then
		val_d = stats%Mem_int_vec
	endif
	if(trim(str)=='Mem_Comp_for')then
		val_d = stats%Mem_Comp_for
	endif


	if(trim(str)=='Rank_max')then
		val_d = maxval(stats%rankmax_of_level_global)
	endif


	if(valid_opt==0)write(*,*)'invalid HODLR stats: '//trim(str)

	deallocate(str)


end subroutine C_BPACK_Getstats


!**** C interface of printing statistics
	!stats_Cptr: the structure containing statistics
	!ptree_Cptr: the structure containing process tree
subroutine C_BPACK_Printstats(stats_Cptr,ptree_Cptr) bind(c, name="c_bpack_printstats")
	implicit none
	type(c_ptr) :: stats_Cptr
	type(c_ptr) :: ptree_Cptr
	type(Hstat),pointer::stats
	type(proctree),pointer::ptree

	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	!**** print statistics variables
	call PrintStat(stats,ptree)

end subroutine C_BPACK_Printstats


!**** C interface of initializing option
	!option_Cptr: the structure containing option
subroutine C_BPACK_Createoption(option_Cptr) bind(c, name="c_bpack_createoption")
	implicit none
	type(c_ptr) :: option_Cptr
	type(Hoption),pointer::option

	allocate(option)
	!**** set default hodlr options
	call SetDefaultOptions(option)

	option_Cptr=c_loc(option)

end subroutine C_BPACK_Createoption


!**** C interface of copy option
	!option_Cptr: the structure containing option
	!option_Cptr1: the structure containing option
subroutine C_BPACK_Copyoption(option_Cptr,option_Cptr1) bind(c, name="c_bpack_copyoption")
	implicit none
	type(c_ptr) :: option_Cptr,option_Cptr1
	type(Hoption),pointer::option,option1


	call c_f_pointer(option_Cptr, option)

	!****copy hodlr options
	allocate(option1)
	call CopyOptions(option,option1)

	option_Cptr1=c_loc(option1)

end subroutine C_BPACK_Copyoption



!**** C interface of set one entry in option
	!option_Cptr: the structure containing option
subroutine C_BPACK_Setoption(option_Cptr,nam,val_Cptr) bind(c, name="c_bpack_setoption")
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
	if(trim(str)=='ILU')then
	call c_f_pointer(val_Cptr, val_i)
	option%ILU=val_i
	valid_opt=1
	endif
	if(trim(str)=='format')then
	call c_f_pointer(val_Cptr, val_i)
	option%format=val_i
	valid_opt=1
	endif
	if(trim(str)=='verbosity')then
	call c_f_pointer(val_Cptr, val_i)
	option%verbosity=val_i
	valid_opt=1
	endif
	if(trim(str)=='rmax')then
	call c_f_pointer(val_Cptr, val_i)
	option%rmax=val_i
	valid_opt=1
	endif
	if(trim(str)=='forwardN15flag')then
	call c_f_pointer(val_Cptr, val_i)
	option%forwardN15flag=val_i
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
	if(trim(str)=='near_para')then
	call c_f_pointer(val_Cptr, val_d)
	option%near_para=val_d
	valid_opt=1
	endif

	if(valid_opt==0)write(*,*)'invalid HODLR option: '//trim(str)

	deallocate(str)
	option_Cptr=c_loc(option)

end subroutine C_BPACK_Setoption



!**** C interface of matrix construction
	!Npo: matrix size
	!Ndim: data set dimensionality (not used if nogeo=1)
	!Locations: coordinates used for clustering (not used if nogeo=1)
	!Nmin: leafsize in HODLR tree
	!tol: compression tolerance
	!nlevel: the number of top levels that have been ordered
	!tree: the order tree provided by the caller
	!Permutation: return the permutation vector new2old (indexed from 1)
	!Npo_loc: number of local row/column indices
	!bmat_Cptr: the structure containing HODLR
	!option_Cptr: the structure containing option
	!stats_Cptr: the structure containing statistics
	!msh_Cptr: the structure containing points and ordering information
	!ker_Cptr: the structure containing kernel quantities
	!ptree_Cptr: the structure containing process tree
	!C_FuncZmn: the C_pointer to user-provided function to sample mn^th entry of the matrix
	!C_QuantApp: the C_pointer to user-defined quantities required to for entry evaluation and sampling
	!MPIcomm: user-provided MPI communicator
subroutine C_BPACK_Construct(Npo,Ndim,Locations,nlevel,tree,Permutation,Npo_loc,bmat_Cptr,option_Cptr,stats_Cptr,msh_Cptr,ker_Cptr,ptree_Cptr,C_FuncZmn,C_QuantApp,MPIcomm) bind(c, name="c_bpack_construct")
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
	type(c_ptr) :: bmat_Cptr
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
	type(Bmatrix),pointer::bmat
	type(proctree),pointer::ptree
	integer seed_myid(50)
	integer times(8)
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	character(len=1024)  :: strings

	!**** allocate HODLR solver structures
	allocate(bmat)

	! allocate(option)
	! allocate(stats)
	allocate(msh)
	allocate(ker)

	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)

	stats%Flop_Fill=0
	stats%Time_Fill=0


	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'NUMBER_MPI=',ptree%nproc

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
	! call RANDOM_SEED(PUT=seed_myid)
	call init_random_seed()

	if(ptree%MyID==Main_ID)then
    write(*,*) "HODLR_BUTTERFLY_SOLVER"
    write(*,*) "   "
	endif

	!**** register the user-defined function and type in ker
	ker%C_QuantApp => C_QuantApp
	ker%C_FuncZmn => C_FuncZmn


	msh%Nunk = Npo


	t1 = OMP_get_wtime()


	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "User-supplied kernel:"
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
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "User-supplied kernel requiring reorder:"
		Dimn = Ndim
		allocate (msh%xyz(Dimn,0:msh%Nunk))
		ii=0
		do edge=1,msh%Nunk
			msh%xyz(1:Dimn,edge)= Locations(ii+1:ii+Dimn)
			ii = ii + Dimn
		enddo
	endif


    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()


	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Hierarchical format......"
    call Cluster_partition(bmat,option,msh,ker,element_Zmn_user_C,ptree)
	call BPACK_structuring(bmat,option,msh,ptree,stats)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Hierarchical format finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()

	
	!**** computation of the construction phase
    call BPACK_construction_Element(bmat,option,stats,msh,ker,element_Zmn_user_C,ptree)


	!**** return the permutation vector
	select case(option%format)
	case(HODLR)
		msh%idxs = bmat%ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		msh%idxe = bmat%ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)
	case(HMAT)
		msh%idxs = bmat%h_mat%Local_blocks(1,1)%headm
		msh%idxe = bmat%h_mat%Local_blocks(1,1)%headm+bmat%h_mat%Local_blocks(1,1)%M-1
	end select

	Npo_loc = msh%idxe-msh%idxs+1
	if(ptree%MyID==Main_ID)then
		do edge=1,Npo
			Permutation(edge) = msh%new2old(edge)
		enddo
	endif

	!**** return the C address of hodlr structures to C caller
	bmat_Cptr=c_loc(bmat)
	option_Cptr=c_loc(option)
	stats_Cptr=c_loc(stats)
	msh_Cptr=c_loc(msh)
	ker_Cptr=c_loc(ker)
	ptree_Cptr=c_loc(ptree)



end subroutine C_BPACK_Construct






!**** C interface of matrix construction via blackbox matvec
	!N: matrix size (in)
	!nlevel: the number of top levels that have been ordered (in)
	!tree: the order tree provided by the caller, if incomplete, the init routine will make it complete (inout)
	!Permutation: return the permutation vector new2old (indexed from 1) (out)
	!N_loc: number of local row/column indices (out)
	!bmat_Cptr: the structure containing HODLR (out)
	!option_Cptr: the structure containing option (in)
	!stats_Cptr: the structure containing statistics (inout)
	!msh_Cptr: the structure containing points and ordering information (out)
	!ker_Cptr: the structure containing kernel quantities (out)
	!ptree_Cptr: the structure containing process tree (in)
subroutine c_bpack_construct_Matvec_Init(N,nlevel,tree,Permutation,N_loc,bmat_Cptr,option_Cptr,stats_Cptr,msh_Cptr,ker_Cptr,ptree_Cptr) bind(c, name="c_bpack_construct_matvec_init")
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
	type(c_ptr) :: bmat_Cptr
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
	type(Bmatrix),pointer::bmat
	type(proctree),pointer::ptree
	integer seed_myid(50)
	integer times(8)
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	real(kind=8):: Memory=0d0,error
	character(len=1024)  :: strings


	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)


	!**** allocate HODLR solver structures
	allocate(bmat)
	! allocate(option)
	! allocate(stats)
	allocate(msh)
	allocate(ker)


	stats%Flop_Fill=0
	stats%Time_Fill=0


	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'NUMBER_MPI=',ptree%nproc

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
	! call RANDOM_SEED(PUT=seed_myid)
	call init_random_seed()

	if(ptree%MyID==Main_ID)then
    write(*,*) "HODLR_BUTTERFLY_SOLVER"
    write(*,*) "   "
	endif

	! !**** register the user-defined function and type in ker
	! ker%C_QuantApp => C_QuantApp
	! ker%C_FuncHMatVec => C_FuncHMatVec


	msh%Nunk = N


	t1 = OMP_get_wtime()


	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "User-supplied kernel:"
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

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()


	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Hierarchical format......"
    call Cluster_partition(bmat,option,msh,ker,element_Zmn_user_C,ptree)
	call BPACK_structuring(bmat,option,msh,ptree,stats)
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Hierarchical format finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "
	t2 = OMP_get_wtime()

	!**** return the permutation vector
	select case(option%format)
	case(HODLR)
		msh%idxs = bmat%ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		msh%idxe = bmat%ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)
	case(HMAT)
		msh%idxs = bmat%h_mat%Local_blocks(1,1)%headm
		msh%idxe = bmat%h_mat%Local_blocks(1,1)%headm+bmat%h_mat%Local_blocks(1,1)%M-1
	end select

	N_loc = msh%idxe-msh%idxs+1
	if(ptree%MyID==Main_ID)then
		do edge=1,N
			Permutation(edge) = msh%new2old(edge)
		enddo
	endif

	!**** return the C address of hodlr structures to C caller
	bmat_Cptr=c_loc(bmat)
	option_Cptr=c_loc(option)
	stats_Cptr=c_loc(stats)
	msh_Cptr=c_loc(msh)
	ker_Cptr=c_loc(ker)
	ptree_Cptr=c_loc(ptree)

end subroutine c_bpack_construct_Matvec_Init



!**** C interface of matrix construction via blackbox matvec
	!bmat_Cptr: the structure containing HODLR (inout)
	!option_Cptr: the structure containing option (in)
	!stats_Cptr: the structure containing statistics (inout)
	!msh_Cptr: the structure containing points and ordering information (in)
	!ker_Cptr: the structure containing kernel quantities (inout)
	!ptree_Cptr: the structure containing process tree (in)
	!C_FuncHMatVec: the C_pointer to user-provided function to multiply A and A* with vectors (in)
	!C_QuantApp: the C_pointer to user-defined quantities required to for entry evaluation and sampling (in)
subroutine c_bpack_construct_Matvec_Compute(bmat_Cptr,option_Cptr,stats_Cptr,msh_Cptr,ker_Cptr,ptree_Cptr,C_FuncHMatVec,C_QuantApp) bind(c, name="c_bpack_construct_matvec_compute")
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
	type(c_ptr) :: bmat_Cptr
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
	type(Bmatrix),pointer::bmat
	type(proctree),pointer::ptree
	integer seed_myid(50)
	integer times(8)
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	real(kind=8):: Memory=0d0,error
	character(len=1024)  :: strings


	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	call c_f_pointer(bmat_Cptr, bmat)
	call c_f_pointer(msh_Cptr, msh)
	call c_f_pointer(ker_Cptr, ker)



	!**** register the user-defined function and type in ker
	ker%C_QuantApp => C_QuantApp
	ker%C_FuncHMatVec => C_FuncHMatVec


	!**** computation of the construction phase
	call BPACK_construction_Matvec(bmat,matvec_user_C,Memory,error,option,stats,ker,ptree,msh)


	!**** return the C address of hodlr structures to C caller
	bmat_Cptr=c_loc(bmat)
	option_Cptr=c_loc(option)
	stats_Cptr=c_loc(stats)
	msh_Cptr=c_loc(msh)
	ker_Cptr=c_loc(ker)
	ptree_Cptr=c_loc(ptree)

end subroutine c_bpack_construct_Matvec_Compute


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
	integer Maxgroup_rc

	!**** allocate HODLR solver structures
	allocate(blocks)
	allocate(msh)
	allocate(ker)

	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	call c_f_pointer(mshr_Cptr, mshr)
	call c_f_pointer(mshc_Cptr, mshc)

	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)'NUMBER_MPI=',ptree%nproc

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
	Maxgroup_rc = min(mshc%Maxgroup,mshr%Maxgroup)
	! call assert(mshc%Maxgroup==mshr%Maxgroup,'mshc%Maxgroup\=mshr%Maxgroup')

	msh%Nunk = N+M
	msh%Maxgroup=Maxgroup_rc*2+1
	allocate (msh%basis_group(msh%Maxgroup))
	msh%basis_group(1)%head=1
	msh%basis_group(1)%tail=N+M
	msh%basis_group(1)%pgno=1
	call copy_basis_group(mshr%basis_group,1,Maxgroup_rc,msh%basis_group,2,msh%Maxgroup,0)
	call copy_basis_group(mshc%basis_group,1,Maxgroup_rc,msh%basis_group,3,msh%Maxgroup,mshr%Nunk)

	blocks%level=1
	blocks%col_group=3
	blocks%row_group=2
	blocks%pgno=1
	blocks%headm=msh%basis_group(blocks%row_group)%head
	blocks%headn=msh%basis_group(blocks%col_group)%head
	blocks%M=M
	blocks%N=N
	blocks%style = 2

	Maxlevel = GetTreelevel(msh%Maxgroup)-1

	if(blocks%level>option%LRlevel)then
		blocks%level_butterfly = 0 ! low rank below LRlevel
	else
		blocks%level_butterfly = Maxlevel - blocks%level   ! butterfly
	endif

	call ComputeParallelIndices(blocks,blocks%pgno,ptree,msh,0)



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
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "FastMATVEC-based BF construction......"

	call BF_randomized(blocks%level_butterfly,option%rank0,option%rankrate,blocks,ker,Bmatvec_user_C,error,'CMatVec',option,stats,ptree,msh)

	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "FastMATVEC-based BF construction finished"

	!**** return the C address of hodlr structures to C caller
	bf_Cptr=c_loc(blocks)
	option_Cptr=c_loc(option)
	stats_Cptr=c_loc(stats)
	msh_Cptr=c_loc(msh)
	ker_Cptr=c_loc(ker)
	ptree_Cptr=c_loc(ptree)

end subroutine C_BF_Construct_Matvec_Compute



!**** C interface of HODLR factorization
	!bmat_Cptr: the structure containing HODLR
	!option_Cptr: the structure containing option
	!stats_Cptr: the structure containing statistics
	!ptree_Cptr: the structure containing process tree
subroutine C_BPACK_Factor(bmat_Cptr,option_Cptr,stats_Cptr,ptree_Cptr,msh_Cptr) bind(c, name="c_bpack_factor")
	implicit none

	type(c_ptr), intent(inout) :: bmat_Cptr
	type(c_ptr), intent(in) :: ptree_Cptr
	type(c_ptr), intent(in) :: msh_Cptr
	type(c_ptr) :: option_Cptr
	type(c_ptr), intent(inout) :: stats_Cptr


	type(Hoption),pointer::option
	type(Hstat),pointer::stats
	type(Bmatrix),pointer::bmat
	type(proctree),pointer::ptree
	type(mesh),pointer::msh

	! real(kind=8):: tol_fact

	call c_f_pointer(bmat_Cptr, bmat)
	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	call c_f_pointer(msh_Cptr, msh)

	stats%Flop_Factor=0
	stats%Time_Factor=0


    call BPACK_factorization(bmat,option,stats,ptree,msh)

	! return the C address of hodlr structures to C caller
	bmat_Cptr=c_loc(bmat)

end subroutine C_BPACK_Factor


!**** C interface of HODLR solve
	!x: local solution vector
	!b: local RHS
	!Nloc: size of local RHS
	!Nrhs: number of RHSs
	!bmat_Cptr: the structure containing HODLR
	!option_Cptr: the structure containing option
	!stats_Cptr: the structure containing statistics
	!ptree_Cptr: the structure containing process tree
subroutine C_BPACK_Solve(x,b,Nloc,Nrhs,bmat_Cptr,option_Cptr,stats_Cptr,ptree_Cptr) bind(c, name="c_bpack_solve")
	implicit none

	integer Nloc,Nrhs
	DT::x(Nloc,Nrhs),b(Nloc,Nrhs)

	type(c_ptr), intent(in) :: bmat_Cptr
	type(c_ptr), intent(in) :: ptree_Cptr
	type(c_ptr), intent(in) :: option_Cptr
	type(c_ptr), intent(in) :: stats_Cptr



	type(Hoption),pointer::option
	type(Hstat),pointer::stats
	type(Bmatrix),pointer::bmat

	type(proctree),pointer::ptree


	call c_f_pointer(bmat_Cptr, bmat)

	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)

	stats%Flop_Sol=0
	stats%Time_Sol=0

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Solve ......"

	if(option%ErrSol==1)then
		call BPACK_Test_Solve_error(bmat,Nloc,option,ptree,stats)
	endif

	call BPACK_Solution(bmat,x,b,Nloc,Nrhs,option,ptree,stats)

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Solve finished"
    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "


end subroutine C_BPACK_Solve




!**** C interface of butterfly-vector multiplication
	!xin: input vector
    !Ninloc:size of local input vectors
	!xout: output vector
	!Noutloc:size of local output vectors
	!Ncol: number of vectors
	!bf_for_Cptr: the structure containing butterfly
	!option_Cptr: the structure containing options
	!stats_Cptr: the structure containing statistics
	!ptree_Cptr: the structure containing process tree
subroutine C_BF_Mult(trans,xin,xout,Ninloc,Noutloc,Ncol,bf_for_Cptr,option_Cptr,stats_Cptr,ptree_Cptr) bind(c, name="c_bf_mult")
	implicit none
	real(kind=8) t1,t2
	integer Ninloc,Noutloc,Ncol
	DT::xin(Ninloc,Ncol),xout(Noutloc,Ncol)

	character(kind=c_char,len=1) :: trans(*)
	type(c_ptr), intent(in) :: bf_for_Cptr
	type(c_ptr), intent(in) :: ptree_Cptr
	type(c_ptr), intent(in) :: stats_Cptr
	type(c_ptr), intent(in) :: option_Cptr


	integer strlen
	character(len=:),allocatable :: str
	type(Hstat),pointer::stats
	type(Hoption),pointer::option
	type(matrixblock),pointer::blocks

	type(proctree),pointer::ptree

	t1 = OMP_get_wtime()

	strlen=1
	! do while(trans(strlen) /= c_null_char)
		! strlen = strlen + 1
	! enddo
	! strlen = strlen -1
	allocate(character(len=strlen) :: str)
	str = transfer(trans(1:strlen), str)

	call c_f_pointer(bf_for_Cptr, blocks)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(ptree_Cptr, ptree)

	stats%Flop_C_Mult=0
	stats%Time_C_Mult=0
	xout=0
	if(trim(str)=='N')then
		call BF_block_MVP_dat(blocks,trim(str),Noutloc,Ninloc,Ncol,xin,xout,cone,czero,ptree,stats)
	else
		call BF_block_MVP_dat(blocks,trim(str),Ninloc,Noutloc,Ncol,xin,xout,cone,czero,ptree,stats)
	endif

	t2 = OMP_get_wtime()

	xout = xout/option%scale_factor

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
	!bmat_Cptr: the structure containing HODLR
	!option_Cptr: the structure containing option
	!stats_Cptr: the structure containing statistics
	!ptree_Cptr: the structure containing process tree
subroutine C_BPACK_Mult(trans,xin,xout,Ninloc,Noutloc,Ncol,bmat_Cptr,option_Cptr,stats_Cptr,ptree_Cptr) bind(c, name="c_bpack_mult")
	implicit none
	real(kind=8) t1,t2
	integer Ninloc,Noutloc,Ncol
	DT::xin(Ninloc,Ncol),xout(Noutloc,Ncol)

	character(kind=c_char,len=1) :: trans(*)
	type(c_ptr), intent(in) :: bmat_Cptr
	type(c_ptr), intent(in) :: ptree_Cptr
	type(c_ptr), intent(in) :: option_Cptr
	type(c_ptr), intent(in) :: stats_Cptr


	integer strlen
	character(len=:),allocatable :: str
	type(Hoption),pointer::option
	type(Hstat),pointer::stats
	type(Bmatrix),pointer::bmat

	type(proctree),pointer::ptree

	t1 = OMP_get_wtime()

	strlen=1
	! do while(trans(strlen) /= c_null_char)
		! strlen = strlen + 1
	! enddo
	! strlen = strlen -1
	allocate(character(len=strlen) :: str)
	str = transfer(trans(1:strlen), str)

	call c_f_pointer(bmat_Cptr, bmat)

	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)

	stats%Flop_C_Mult=0
	stats%Time_C_Mult=0

    ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Multiply ......"

	call assert(Noutloc==Ninloc,"not square Z")
	call BPACK_Mult(trim(str),Noutloc,Ncol,xin,xout,bmat,ptree,option,stats)
	! need to use another Flop counter for this operation in future

    ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Multiply finished"
    ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "

	t2 = OMP_get_wtime()

	stats%Time_C_Mult = stats%Time_C_Mult + t2-t1
	stats%Flop_C_Mult = stats%Flop_C_Mult + stats%Flop_Tmp

	! write(*,*)t2-t1
	deallocate(str)
end subroutine C_BPACK_Mult




!**** C interface of HODLR(inverse)-vector multiplication
	!xin: input vector
    !Ninloc:size of local input vectors
	!xout: output vector
	!Noutloc:size of local output vectors
	!Ncol: number of vectors
	!bmat_Cptr: the structure containing HODLR
	!option_Cptr: the structure containing option
	!stats_Cptr: the structure containing statistics
	!ptree_Cptr: the structure containing process tree
subroutine C_BPACK_Inv_Mult(trans,xin,xout,Ninloc,Noutloc,Ncol,bmat_Cptr,option_Cptr,stats_Cptr,ptree_Cptr) bind(c, name="c_bpack_inv_mult")
	implicit none
	real(kind=8) t1,t2
	integer Ninloc,Noutloc,Ncol
	DT::xin(Ninloc,Ncol),xout(Noutloc,Ncol)

	character(kind=c_char,len=1) :: trans(*)
	type(c_ptr), intent(in) :: bmat_Cptr
	type(c_ptr), intent(in) :: ptree_Cptr
	type(c_ptr), intent(in) :: option_Cptr
	type(c_ptr), intent(in) :: stats_Cptr


	integer strlen
	character(len=:),allocatable :: str
	type(Hoption),pointer::option
	type(Hstat),pointer::stats
	type(Bmatrix),pointer::bmat

	type(proctree),pointer::ptree

	t1 = OMP_get_wtime()

	strlen=1
	! do while(trans(strlen) /= c_null_char)
		! strlen = strlen + 1
	! enddo
	! strlen = strlen -1
	allocate(character(len=strlen) :: str)
	str = transfer(trans(1:strlen), str)

	call c_f_pointer(bmat_Cptr, bmat)

	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)

	stats%Flop_C_Mult=0
	stats%Time_C_Mult=0

    ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Multiply ......"

	call assert(Noutloc==Ninloc,"not square Z")
	call BPACK_Inv_Mult(trim(str),Noutloc,Ncol,xin,xout,bmat,ptree,option,stats)
	! need to use another Flop counter for this operation in future

    ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "Multiply finished"
    ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "

	t2 = OMP_get_wtime()

	stats%Time_C_Mult = stats%Time_C_Mult + t2-t1
	stats%Flop_C_Mult = stats%Flop_C_Mult + stats%Flop_Tmp

	! write(*,*)t2-t1
	deallocate(str)
end subroutine C_BPACK_Inv_Mult



!**** C interface of deleting statistics
	!stats_Cptr: the structure containing statistics
subroutine C_BPACK_Deletestats(stats_Cptr) bind(c, name="c_bpack_deletestats")
	implicit none
	type(c_ptr), intent(inout) :: stats_Cptr
	type(Hstat),pointer::stats

	call c_f_pointer(stats_Cptr, stats)
	call delete_Hstat(stats)
	deallocate(stats)
	stats_Cptr=c_null_ptr

end subroutine C_BPACK_Deletestats


!**** C interface of deleting process tree
	!ptree_Cptr: the structure containing process tree
subroutine C_BPACK_Deleteproctree(ptree_Cptr) bind(c, name="c_bpack_deleteproctree")
	implicit none
	type(c_ptr), intent(inout) :: ptree_Cptr
	type(proctree),pointer::ptree

	call c_f_pointer(ptree_Cptr, ptree)
	call delete_proctree(ptree)
	deallocate(ptree)
	ptree_Cptr=c_null_ptr

end subroutine C_BPACK_Deleteproctree


!**** C interface of deleting mesh
	!msh_Cptr: the structure containing mesh
subroutine C_BPACK_Deletemesh(msh_Cptr) bind(c, name="c_bpack_deletemesh")
	implicit none
	type(c_ptr), intent(inout) :: msh_Cptr
	type(mesh),pointer::msh

	call c_f_pointer(msh_Cptr, msh)
	call delete_mesh(msh)
	deallocate(msh)
	msh_Cptr=c_null_ptr

end subroutine C_BPACK_Deletemesh


!**** C interface of deleting kernelquant
	!ker_Cptr: the structure containing kernelquant
subroutine C_BPACK_Deletekernelquant(ker_Cptr) bind(c, name="c_bpack_deletekernelquant")
	implicit none
	type(c_ptr), intent(inout) :: ker_Cptr
	type(kernelquant),pointer::ker

	call c_f_pointer(ker_Cptr, ker)
	call delete_kernelquant(ker)
	deallocate(ker)
	ker_Cptr=c_null_ptr

end subroutine C_BPACK_Deletekernelquant


!**** C interface of deleting HOBF
	!bmat_Cptr: the structure containing HOBF
subroutine C_BPACK_Delete(bmat_Cptr) bind(c, name="c_bpack_delete")
	implicit none
	type(c_ptr), intent(inout) :: bmat_Cptr
	type(Bmatrix),pointer::bmat

	call c_f_pointer(bmat_Cptr, bmat)
	call BPACK_delete(bmat)
	deallocate(bmat)
	bmat_Cptr=c_null_ptr

end subroutine C_BPACK_Delete


!**** C interface of deleting a BF
	!bf_Cptr: the structure containing BF
subroutine C_BF_DeleteBF(bf_Cptr) bind(c, name="c_bf_deletebf")
	implicit none
	type(c_ptr), intent(inout) :: bf_Cptr
	type(matrixblock),pointer::blocks

	call c_f_pointer(bf_Cptr, blocks)
	call BF_delete(blocks,1)
	deallocate(blocks)
	bf_Cptr=c_null_ptr
end subroutine C_BF_DeleteBF



!**** C interface of deleting Hoption
	!option_Cptr: the structure containing Hoption
subroutine C_BPACK_Deleteoption(option_Cptr) bind(c, name="c_bpack_deleteoption")
	implicit none
	type(c_ptr), intent(inout) :: option_Cptr
	type(Hoption),pointer::option

	call c_f_pointer(option_Cptr, option)
	deallocate(option)
	option_Cptr=c_null_ptr

end subroutine C_BPACK_Deleteoption


end module BPACK_wrapper

