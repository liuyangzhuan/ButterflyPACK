module HODLR_C_Interface
    use MODULE_FILE
	! use geometry_model
	use H_structure
	use cascading_factorization
	use matrices_fill
	use omp_lib
	use MISC
	use HODLR_Solve
    use iso_c_binding    
	
contains 



subroutine element_Zmn_user_C(edge_m,edge_n,value_e,msh,ker)
    
    use MODULE_FILE
    implicit none
    
    integer edge_m, edge_n
    complex(kind=8) value_e
	type(mesh)::msh
	type(kernelquant)::ker	
	
	procedure(C_Z_elem), POINTER :: proc
	
	value_e=0
	call c_f_procpointer(ker%C_FuncZmn, proc)
	! write(*,*)'ddd'
	call proc(msh%info_unk(0,edge_m)-1,msh%info_unk(0,edge_n)-1,value_e,ker%C_QuantZmn)
	! call proc(msh%info_unk(0,edge_m)-1,msh%info_unk(0,edge_n)-1,value_e)
	return
    
end subroutine element_Zmn_user_C





subroutine C_HODLR_Fill(Npo,Ndim,Locations,Nmin,tol,nth,nmpi,ninc,preorder,tree,Permutation,Npo_loc,ho_bf_Cptr,option_Cptr,stats_Cptr,msh_Cptr,ker_Cptr,ptree_Cptr,C_FuncZmn,C_QuantZmn,MPIcomm) bind(c, name="c_hodlr_fill_")	
	implicit none 
	integer Npo,Ndim,Nmin
	real*8 Locations(Npo*Ndim)
	
    real*8 para
    real*8 tolerance,tol,h,lam
    integer Primary_block, nn, mm, MyID_old
    integer i,j,k,ii,edge, threads_num,nth,Dimn,nmpi,ninc, acam
	real(kind=8),parameter :: cd = 299792458d0
	integer,allocatable:: groupmembers(:)
	integer Permutation(*),tree(*)
	integer Npo_loc
	! type(matricesblock), pointer :: blocks_i
	integer groupm,preorder
	integer MPIcomm
	type(c_ptr), intent(out) :: ho_bf_Cptr
	type(c_ptr), intent(out) :: option_Cptr
	type(c_ptr), intent(out) :: stats_Cptr
	type(c_ptr), intent(out) :: msh_Cptr
	type(c_ptr), intent(out) :: ker_Cptr
	type(c_ptr), intent(out) :: ptree_Cptr
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
	real*8 t1,t2,x,y,z,r,theta,phi
	complex(kind=8):: Ctmp

	allocate(ho_bf)
	allocate(option)
	allocate(stats)
	allocate(msh)
	allocate(ker)
	allocate(ptree)
	

	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)*ninc
	enddo

	call CreatePtree(nmpi,groupmembers,MPIcomm,ptree)
	deallocate(groupmembers)
	
	if(ptree%MyID==Main_ID)write(*,*)'NUMBER_MPI=',nmpi
 	threads_num=nth
	if(ptree%MyID==Main_ID)write(*,*)'OMP_NUM_THREADS=',threads_num
	call OMP_set_num_threads(threads_num)		
		
		
	call DATE_AND_TIME(values=times)     ! Get the current time 
	seed_myid(1) = times(4) * (360000*times(5) + 6000*times(6) + 100*times(7) + times(8))
	! seed_myid(1) = myid*1000
	call RANDOM_SEED(PUT=seed_myid)
	
	if(ptree%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "HODLR_BUTTERFLY_SOLVER"
    write(*,*) "   "
	endif
	
	call InitStat(stats)

	time_tmp = 0
	
	msh%Origins=(/0d0,0d0,0d0/)
 

     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	
	ker%C_QuantZmn => C_QuantZmn
	ker%C_FuncZmn => C_FuncZmn
		
	para=0.001
	
	! ker%rank_approximate_para1=6.0
    ! ker%rank_approximate_para2=6.0
    ! ker%rank_approximate_para3=6.0
	
	
	option%Nmin_leaf=Nmin
	option%tol_SVD=tol
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

   !***********************************************************************
   if(ptree%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'User-supplied kernel with geometry:'
   ! write (*,*) 'wavelength:',ker%wavelength
   write (*,*) ''
   endif
   !***********************************************************************
	
	
	
	
	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID)write(*,*) "geometry modeling......"
	
	Dimn = Ndim 
	msh%Ncorner = 0
	msh%Nunk = Npo

	allocate (msh%xyz(Dimn,0:msh%Nunk), msh%info_unk(0:0,msh%Nunk))
	! write(*,*)msh%Nunk, Dimn,shape(msh%info_unk)
	ii=0
	do edge=1,msh%Nunk
		msh%info_unk(0,edge)=edge
		msh%xyz(1:Dimn,edge)= Locations(ii+1:ii+Dimn)
		ii = ii + Dimn
	enddo  	
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
	
	
	
	

	! call c_f_procpointer(C_FuncZmn, proc)
	! ! write(*,*)'ddd'
	! call proc(0,0,Ctmp,C_QuantZmn)	
	! write(*,*)'1 1 element is: ', Ctmp
	
	! call element_Zmn_user_C(1,1,Ctmp,msh,ker)
	! write(*,*)'1 1 element is: ', Ctmp
	! call element_Zmn_user_C(3,5,Ctmp,msh,ker)
	! write(*,*)'3 5 element is: ', Ctmp
	! stop
	
    !pause
    
    !call compression_test()
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices filling......"
    call matrices_filling(ho_bf,option,stats,msh,ker,element_Zmn_user_C,ptree)
	! if(option%precon/=DIRECT)then
		! call copy_HOBF(ho_bf,ho_bf_copy)	
	! end if
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices filling finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	! write(*,*)t2-t1

	msh%idxs = ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
	msh%idxe = ho_bf%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)	
		
	Npo_loc = msh%idxe-msh%idxs+1		
	

	if (ptree%MyID==Main_ID) then	
		do edge=1,Npo
			Permutation(edge) = msh%info_unk(0,edge)
		enddo
	endif	
	
	! return the C address of hodlr structures to C caller
	ho_bf_Cptr=c_loc(ho_bf)
	option_Cptr=c_loc(option)
	stats_Cptr=c_loc(stats)
	msh_Cptr=c_loc(msh)
	ker_Cptr=c_loc(ker)
	ptree_Cptr=c_loc(ptree)


	
end subroutine C_HODLR_Fill




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

	! real*8:: tol_fact
	
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



subroutine C_HODLR_Solve(x,b,Nloc,Nrhs,ho_bf_for_Cptr,ho_bf_inv_Cptr,option_Cptr,stats_Cptr,ptree_Cptr) bind(c, name="c_hodlr_solve_")	
	implicit none 

	integer Nloc,Nrhs
	complex(kind=8)::x(Nloc,Nrhs),b(Nloc,Nrhs)
	
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
	! real*8 Xin(*),Xout(*)
	! ! real*8 Xin(Npo*Ncol),Xout(Npo*Ncol)
	! real*8 n1,n2
	
    ! integer i, j, ii, jj, iii, jjj, num_blocks, mm, nn
    ! integer level, blocks, edge, patch, node, group, groupm,Maxgroup_loc,g_start_glo,g_start_loc
    ! integer rank, index_near, m, n, length, flag, num_sample,vectors_x,vectors_y,vectors_start, Dimn
    ! real*8 theta, phi, dphi, rcs_V, rcs_H, T0, T1, vecnorm, rtemp
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
