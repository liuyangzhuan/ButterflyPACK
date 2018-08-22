module HODLR_C_Interface
    use MODULE_FILE
	! use geometry_model
	use H_structure
	use cascading_factorization
	use HODLR_construction
	use omp_lib
	use MISC
	use HODLR_Solve
    use iso_c_binding    
	
#include "HODLR_config.fi"
	
contains 



subroutine element_Zmn_user_C(edge_m,edge_n,value_e,msh,ker)
    use MODULE_FILE
    implicit none
    integer edge_m, edge_n
    DT value_e
	type(mesh)::msh
	type(kernelquant)::ker	
	procedure(C_Z_elem), POINTER :: proc
	
	value_e=0
	call c_f_procpointer(ker%C_FuncZmn, proc)
	call proc(msh%info_unk(0,edge_m)-1,msh%info_unk(0,edge_n)-1,value_e,ker%C_QuantZmn)
	return
    
end subroutine element_Zmn_user_C


!**** C interface of process tree construction
	!nmpi: number of MPIs for one hodlr
	!MPIcomm: MPI communicator from C caller
	!groupmembers: MPI ranks in MPIcomm for one hodlr
	!ptree_Cptr: the structure containing process tree
subroutine C_CreatePtree(nmpi,groupmembers,MPIcomm,ptree_Cptr) bind(c, name="c_createptree_")	
	implicit none 
	integer nmpi
	integer MPIcomm
	integer:: groupmembers(*)
	type(c_ptr):: ptree_Cptr	
	type(proctree),pointer::ptree	

	allocate(ptree)
	call CreatePtree(nmpi,groupmembers,MPIcomm,ptree)
	ptree_Cptr=c_loc(ptree)
end subroutine C_CreatePtree


!**** C interface of initializing statistics
	!nmpi: number of MPIs for one hodlr
	!MPIcomm: MPI communicator from C caller
	!groupmembers: MPI ranks in MPIcomm for one hodlr
	!stats_Cptr: the structure containing statistics
subroutine C_CreateStats(stats_Cptr) bind(c, name="c_createstats_")	
	implicit none 
	type(c_ptr), intent(out) :: stats_Cptr
	type(Hstat),pointer::stats

	allocate(stats)
	!**** initialize statistics variables  	
	call InitStat(stats)	
	stats_Cptr=c_loc(stats)
	
end subroutine C_CreateStats


!**** C interface of initializing option
	!option_Cptr: the structure containing option       
subroutine C_CreateOption(option_Cptr) bind(c, name="c_createoption_")	
	implicit none 
	type(c_ptr) :: option_Cptr
	type(Hoption),pointer::option	
	
	allocate(option)
	!**** set default hodlr options 	
	call SetDefaultOptions(option)
	
	option_Cptr=c_loc(option)
	
end subroutine C_CreateOption



!**** C interface of set one entry in option
	!option_Cptr: the structure containing option       
subroutine C_SetOption(option_Cptr,nam,val_Cptr) bind(c, name="c_setoption_")	
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
	if(trim(str)=='touch_para')then
	call c_f_pointer(val_Cptr, val_i)
	option%touch_para=val_i
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
	if(trim(str)=='ErrSol')then
	call c_f_pointer(val_Cptr, val_i)
	option%ErrSol=val_i
	valid_opt=1
	endif
	if(trim(str)=='preorder')then
	call c_f_pointer(val_Cptr, val_i)
	option%preorder=val_i
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
    
	if(valid_opt==0)write(*,*)'invalid HODLR option: '//trim(str)
	
	deallocate(str)
	option_Cptr=c_loc(option)
	
end subroutine C_SetOption



!**** C interface of HODLR construction
	!Npo: matrix size
	!Ndim: data set dimensionality (not used if preorder=1)
	!Locations: coordinates used for clustering (not used if preorder=1) 
	!Nmin: leafsize in HODLR tree (not used if preorder=1 and nlevel>0)
	!tol: compression tolerance
	!nth: number of threads
	!nmpi: number of processes calling this subroutine 
	!ninc: pick 1 out of ninc processes to create a new communicator for HODLR
	!preorder: 0: matrix not pre-ordered,requiring a clustering step 1: matrix pre-ordered 
	!nlevel: if preorder=1: nlevel=0: require a natural order tree; nlevel=1: an order tree is also provided
	!tree: the order tree provided by the caller
	!Permutation: return new2old if preorder=0; return a natural order if preorder=1
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
subroutine C_HODLR_Construct(Npo,Ndim,Locations,nlevel,tree,Permutation,Npo_loc,ho_bf_Cptr,option_Cptr,stats_Cptr,msh_Cptr,ker_Cptr,ptree_Cptr,C_FuncZmn,C_QuantZmn,MPIcomm) bind(c, name="c_hodlr_construct_")	
	implicit none 
	integer Npo,Ndim
	real(kind=8) Locations(*)
	
    real(kind=8) para
    real(kind=8) tolerance,h,lam
    integer Primary_block, nn, mm, MyID_old,Maxlevel
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
	integer seed_myid(12)
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
	call RANDOM_SEED(PUT=seed_myid)
	
	if(ptree%MyID==Main_ID)then
    write(*,*) "HODLR_BUTTERFLY_SOLVER"
    write(*,*) "   "
	endif
	
	!**** register the user-defined function and type in ker 
	ker%C_QuantZmn => C_QuantZmn
	ker%C_FuncZmn => C_FuncZmn	
	

	msh%Origins=(/0d0,0d0,0d0/)
	msh%scaling=1d0
	msh%Ncorner = 0
	msh%Nunk = Npo
		
		
	t1 = OMP_get_wtime()

	!**** the geometry points are provided by user and require reordering 
	if(option%preorder==0)then 
		if(ptree%MyID==Main_ID)write(*,*) "User-supplied kernel with geometrical points:"
		Dimn = Ndim 
		allocate (msh%xyz(Dimn,0:msh%Nunk))
		ii=0
		do edge=1,msh%Nunk
			msh%xyz(1:Dimn,edge)= Locations(ii+1:ii+Dimn)
			ii = ii + Dimn
		enddo  	
	else 
		if(nlevel==0)then  !**** the geometry points are not provided, the preorder tree is not provided, create a natural order tree 
			if(ptree%MyID==Main_ID)write(*,*) "Geometry-free kernel without tree:"
			level=0; i=1
			do while (int(msh%Nunk/i)>option%Nmin_leaf)
				level=level+1
				i=2**level
			enddo
			Maxlevel=level
			allocate(msh%pretree(2**Maxlevel))
			call CreateLeaf_Natural(Maxlevel,0,1,1,msh%Nunk,msh%pretree) 
			
		else   !**** the geometry points are not provided, but the preorder tree is provided 			
			
			if(ptree%MyID==Main_ID)write(*,*) "Geometry-free kernel with tree:"
			Maxlevel=nlevel
			allocate(msh%pretree(2**Maxlevel))
			msh%pretree(1:2**Maxlevel) = tree(1:2**Maxlevel)
		endif
	end if	
	
	!**** initialize the permutation vector (msh%info_unk(0,:) and old2new(:) coincide for user-supplied kernels)  
	allocate (msh%info_unk(0:0,msh%Nunk))
	do ii=1,msh%Nunk
		msh%info_unk(0,ii)=ii
	enddo 	
	
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	

	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "constructing H_matrices formatting......"
    call H_matrix_structuring(ho_bf,option,msh,ptree)
	call BPlus_structuring(ho_bf,option,msh,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices formatting finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	
	
    !call compression_test()
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices filling......"
    call matrices_construction(ho_bf,option,stats,msh,ker,element_Zmn_user_C,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices filling finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	! write(*,*)t2-t1

	
	!**** return the permutation vector if preorder==0
	msh%idxs = ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
	msh%idxe = ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)		
	Npo_loc = msh%idxe-msh%idxs+1		
	if(option%preorder==0)then 
	if (ptree%MyID==Main_ID) then	
		do edge=1,Npo
			Permutation(edge) = msh%info_unk(0,edge)
		enddo
	endif	
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
	!ho_bf_for_Cptr: the structure containing forward HODLR         
	!ho_bf_inv_Cptr: the structure containing factored HODLR         
	!option_Cptr: the structure containing option         
	!stats_Cptr: the structure containing statistics         
	!ptree_Cptr: the structure containing process tree
subroutine C_HODLR_Factor(ho_bf_for_Cptr,ho_bf_inv_Cptr,option_Cptr,stats_Cptr,ptree_Cptr) bind(c, name="c_hodlr_factor_")	
	implicit none 

	type(c_ptr), intent(inout) :: ho_bf_for_Cptr
	type(c_ptr), intent(in) :: ptree_Cptr	
	type(c_ptr) :: option_Cptr
	type(c_ptr), intent(inout) :: stats_Cptr
	type(c_ptr), intent(out) :: ho_bf_inv_Cptr
	

	type(Hoption),pointer::option	
	type(Hstat),pointer::stats
	type(hobf),pointer::ho_bf_for
	type(hobf),pointer::ho_bf_inv,ho_bf_tmp
	type(proctree),pointer::ptree	

	! real(kind=8):: tol_fact
	
	call c_f_pointer(ho_bf_for_Cptr, ho_bf_for)
	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	
	allocate(ho_bf_tmp)
	call copy_HOBF(ho_bf_for,ho_bf_tmp)	! currently this subroutine only copies forward components 
	ho_bf_inv=>ho_bf_for
	ho_bf_for=>ho_bf_tmp
	
	
	if(option%precon/=NOPRECON)then
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing......"
    call cascading_factorizing(ho_bf_inv,option,stats,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	end if
	
	! return the C address of hodlr structures to C caller
	ho_bf_inv_Cptr=c_loc(ho_bf_inv)	
	ho_bf_for_Cptr=c_loc(ho_bf_for)	

end subroutine C_HODLR_Factor


!**** C interface of HODLR solve
	!x: local solution vector         
	!b: local RHS        
	!Nloc: size of local RHS     
	!Nrhs: number of RHSs     
	!ho_bf_for_Cptr: the structure containing forward HODLR         
	!ho_bf_inv_Cptr: the structure containing factored HODLR         
	!option_Cptr: the structure containing option         
	!stats_Cptr: the structure containing statistics         
	!ptree_Cptr: the structure containing process tree
subroutine C_HODLR_Solve(x,b,Nloc,Nrhs,ho_bf_for_Cptr,ho_bf_inv_Cptr,option_Cptr,stats_Cptr,ptree_Cptr) bind(c, name="c_hodlr_solve_")	
	implicit none 

	integer Nloc,Nrhs
	DT::x(Nloc,Nrhs),b(Nloc,Nrhs)
	
	type(c_ptr), intent(in) :: ho_bf_for_Cptr
	type(c_ptr), intent(in) :: ptree_Cptr	
	type(c_ptr), intent(in) :: option_Cptr
	type(c_ptr), intent(in) :: stats_Cptr
	type(c_ptr), intent(in) :: ho_bf_inv_Cptr
	

	type(Hoption),pointer::option	
	type(Hstat),pointer::stats
	type(hobf),pointer::ho_bf_for
	type(hobf),pointer::ho_bf_inv
	type(proctree),pointer::ptree	


	call c_f_pointer(ho_bf_for_Cptr, ho_bf_for)
	call c_f_pointer(ho_bf_inv_Cptr, ho_bf_inv)
	call c_f_pointer(option_Cptr, option)
	call c_f_pointer(stats_Cptr, stats)
	call c_f_pointer(ptree_Cptr, ptree)
	
    if(ptree%MyID==Main_ID)write(*,*) "Solve ......"
	
	if(option%ErrSol==1)then
		call HODLR_Test_Solve_error(ho_bf_for,ho_bf_inv,option,ptree,stats)
	endif		
	
	call HODLR_Solution(ho_bf_for,ho_bf_inv,x,b,Nloc,Nrhs,option,ptree,stats)

    if(ptree%MyID==Main_ID)write(*,*) "Solve finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	
	! return the C address of hodlr structures to C caller
	
	! ho_bf_inv_Cptr=c_loc(ho_bf_inv)	
	! ho_bf_for_Cptr=c_loc(ho_bf_for)	

end subroutine C_HODLR_Solve


! subroutine H_Matrix_Apply(Npo,Ncol,Xin,Xout) bind(c, name="h_matrix_apply_")	
	! implicit none 
	! integer Npo,Ncol,Nmin, Ntot
	! real(kind=8) Xin(*),Xout(*)
	! ! real(kind=8) Xin(Npo*Ncol),Xout(Npo*Ncol)
	! real(kind=8) n1,n2
	
    ! integer i, j, ii, jj, iii, jjj, num_blocks, mm, nn
    ! integer level, blocks, edge, patch, node, group, groupm,Maxgroup_loc,g_start_glo,g_start_loc
    ! integer rank, index_near, m, n, length, flag, num_sample,vectors_x,vectors_y,vectors_start, Dimn
    ! real(kind=8) theta, phi, dphi, rcs_V, rcs_H, T0, T1, vecnorm, rtemp
    ! double complex value_Z, ctemp
    
    ! double complex, allocatable :: ctemp_vector(:), ctemp_vector1(:), ctemp_vector2(:), output(:,:)
    ! integer, allocatable :: intemp_vector1(:), intemp_vector2(:)     

	! type(matricesblock), pointer :: blocks_i,blocks_j
	! complex(kind=8), allocatable:: labels(:)		
	! integer segsize
	
	! if(MPI_Comm_H/=MPI_COMM_NULL)then
	
	
	! call MPI_verbose_barrier('1')
	! n1=OMP_get_wtime()		
	
	! call assert(Maxedge==Npo,'Npo incorrect')
	
	! blocks_i=>Local_matrices_block(1,1)
	! groupm = blocks_i%col_group
	! Maxgroup_loc = 2**(Maxlevel_for_blocks-blocks_i%level+1)-1
	! allocate(vectors_block_o(0:Maxgroup_loc))
	
	! do level=0,Maxlevel_for_blocks-blocks_i%level
		! g_start_loc = 2**level
		! g_start_glo = groupm*2**level
		! do group = 1,2**level
			! vectors_block_o(group+g_start_loc-1)%head=basis_group(group+g_start_glo-1)%head - basis_group(groupm)%head + 1
            ! vectors_block_o(group+g_start_loc-1)%tail=basis_group(group+g_start_glo-1)%tail - basis_group(groupm)%head + 1
            ! vectors_block_o(group+g_start_loc-1)%style=4
            ! vectors_block_o(group+g_start_loc-1)%group=group+g_start_glo-1
		! end do
	! end do	
	! allocate(vectors_block_o(1)%vector(vectors_block_o(1)%tail-vectors_block_o(1)%head+1,Ncol))	
	! vectors_block_o(1)%style=1
	! vectors_block_o(1)%vector=0.0	
	
	
	! num_blocks=2**Parallel_level_inverse_MPI
	! vectors_start=num_blocks-1		
	
	! ! i=MyID+1
	! ! vectors_x=vectors_start+i
	
	! vectors_x=1
	
	! call MPI_barrier(MPI_Comm_H,ierr)
	! n2=OMP_get_wtime()
	
	! if(MyID==0)write(*,*)"before kernel apply time: ",n2-n1 
		
	
	! call MPI_verbose_barrier('1')
	! n1=OMP_get_wtime()
	
	
	! allocate(vectors_block(0:Maxgroup_loc))
	
	! do j=1,num_blocks
	
		! blocks_i=>Local_matrices_block(j,1)
		
		! groupm = blocks_i%row_group
		! Maxgroup_loc = 2**(Maxlevel_for_blocks-blocks_i%level+1)-1
		
		
		! do level=0,Maxlevel_for_blocks-blocks_i%level
			! g_start_loc = 2**level
			! g_start_glo = groupm*2**level
			! do group = 1,2**level
				! vectors_block(group+g_start_loc-1)%head=basis_group(group+g_start_glo-1)%head - basis_group(groupm)%head + 1
				! vectors_block(group+g_start_loc-1)%tail=basis_group(group+g_start_glo-1)%tail - basis_group(groupm)%head + 1
				! vectors_block(group+g_start_loc-1)%style=4
				! vectors_block(group+g_start_loc-1)%group=group+g_start_glo-1
			! end do
		! end do	
		! allocate(vectors_block(1)%vector(vectors_block(1)%tail-vectors_block(1)%head+1,Ncol))	
		! vectors_block(1)%style=1
		! vectors_block(1)%vector=0.0			
		
		! if(MyID==j-1)then
			! ! do jj=1,Ncol
			! ! do ii=basis_group(groupm)%head,basis_group(groupm)%tail
				! ! vectors_block(1)%vector(ii-basis_group(groupm)%head+1,jj) = Xin((jj-1)*Maxedge+node_of_edge(0,ii))
			! ! enddo
			! ! enddo
			! do jj=1,Ncol
			! do ii=1,basis_group(groupm)%tail-basis_group(groupm)%head+1
				! vectors_block(1)%vector(ii,jj) = Xin((jj-1)*(basis_group(groupm)%tail-basis_group(groupm)%head+1)+ii)
			! enddo
			! enddo			
		! endif
		! call MPI_Bcast(vectors_block(1)%vector,(vectors_block(1)%tail-vectors_block(1)%head+1)*Ncol,MPI_double_complex,j-1,MPI_Comm_H,ierr)
		
		! vectors_y=1
		
		! ! write(*,*)vectors_x,vectors_y,vectors_block_o(vectors_x)%style,vectors_block_o(vectors_y)%style,'nima'
		! call Vector_add_multiply_o(vectors_x,'+',blocks_i,'N',vectors_y)	
		! ! call aggregate_vectors_o(vectors_x)				
	
		! deallocate(vectors_block(1)%vector)
	! enddo

	! call MPI_barrier(MPI_Comm_H,ierr)
	! n2=OMP_get_wtime()
	
	! if(MyID==0)write(*,*)"kernel apply time: ",n2-n1 
		


	! call MPI_verbose_barrier('1')
	! n1=OMP_get_wtime()		
	
	
	! blocks_i=>Local_matrices_block(1,1)
	! groupm = blocks_i%col_group
	! segsize = basis_group(groupm)%tail-basis_group(groupm)%head+1
	
	! do jj=1,Ncol
	! do ii=1,segsize
		! Xout((jj-1)*segsize+ii) = dble(vectors_block_o(1)%vector(ii,jj)) 
		! ! write(*,*)Xout((jj-1)*Maxedge+node_of_edge(0,ii))
	! enddo
	! enddo	

	! ! if(MyID/=0)then
		! ! Ntot = (vectors_block_o(1)%tail-vectors_block_o(1)%head+1)*Ncol	
		! ! call MPI_Send(vectors_block_o(1)%vector,Ntot,MPI_double_complex,Main_ID,MyID,MPI_Comm_H,ierr)
	! ! else 
		! ! do i=1,num_blocks
			! ! blocks_i=>Local_matrices_block(i,1)
			! ! groupm = blocks_i%row_group
			! ! Ntot = (basis_group(groupm)%tail - basis_group(groupm)%head + 1)*Ncol		
			! ! allocate(output(basis_group(groupm)%tail - basis_group(groupm)%head + 1,Ncol))	
			! ! if(i-1==0)then
				! ! output = vectors_block_o(1)%vector
			! ! else
				! ! call MPI_Recv(output,Ntot,MPI_double_complex,i-1,i-1,MPI_Comm_H,status,ierr)
			! ! endif
			
			! ! do jj=1,Ncol
			! ! do ii=basis_group(groupm)%head,basis_group(groupm)%tail
				! ! Xout((jj-1)*Maxedge+node_of_edge(0,ii)) = dble(output(ii-basis_group(groupm)%head+1,jj)) 
				! ! ! write(*,*)Xout((jj-1)*Maxedge+node_of_edge(0,ii))
			! ! enddo
			! ! enddo				
			
			! ! deallocate(output)
		! ! enddo
	! ! endif	
	
	! deallocate (vectors_block_o(1)%vector)
	! deallocate (vectors_block_o)
	! deallocate (vectors_block)
	
	! call MPI_barrier(MPI_Comm_H,ierr)
	! n2=OMP_get_wtime()
	
	! if(MyID==Main_ID)write(*,*)"after kernel apply time: ",n2-n1 
	! ! if(MyID==0)write(*,*)"output norm: ",sqrt(sum(Xout(1:segsize*Ncol)**2d0))
	
	! rtemp = sum(Xout(1:segsize*Ncol)**2d0)
	! call MPI_REDUCE(rtemp, vecnorm, 1,MPI_double, MPI_SUM, Main_ID, MPI_Comm_H,ierr)	
	! if(MyID==Main_ID)write(*,*)"output norm: ",sqrt(vecnorm)
	
	! ! if(MyID==0)then
	! ! do ii=1,Maxedge*Ncol
		! ! write(777,*)Xout(ii)
	! ! enddo
	
	! ! endif		
	
	! endif

! end subroutine H_Matrix_Apply



end module HODLR_C_Interface
