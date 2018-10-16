PROGRAM HODLR_BUTTERFLY_SOLVER_3D
    use z_HODLR_DEFS
	use EMSURF_MODULE
	
	use z_HODLR_structure
	use z_HODLR_factor
	use z_HODLR_constr
	use omp_lib
	use z_misc
    implicit none

	! include "mkl_vml.fi"	 
	
    real(kind=8) para
    real(kind=8) tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seed_myid(50)
	integer times(8)	
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	complex(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:)
	
	character(len=:),allocatable  :: string
	character(len=1024)  :: strings	
	character(len=6)  :: info_env	
	integer :: length,edge
	integer :: ierr
	integer*8 oldmode,newmode
	type(z_Hoption)::option	
	type(z_Hstat)::stats	
	type(z_mesh)::msh	
	type(z_hobf)::ho_bf,ho_bf_copy
	type(z_kernelquant)::ker
	type(quant_EMSURF),target::quant
	type(z_proctree)::ptree
	integer,allocatable:: groupmembers(:)	
	integer nmpi
	CHARACTER (LEN=1000) DATA_DIR	
	
	! nmpi and groupmembers should be provided by the user 
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo	
	
	call z_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	deallocate(groupmembers)
	
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
    write(*,*) "HODLR_BUTTERFLY_SOLVER_3D"
    write(*,*) "   "
	endif
	call z_InitStat(stats)
	call z_SetDefaultOptions(option)
	
	time_tmp = 0
	
 	! register the user-defined function and type in ker 
	ker%FuncZmn=>Z_elem_EMSURF
	ker%QuantZmn=>quant	
	
    quant%integral_points=6
    allocate (quant%ng1(quant%integral_points), quant%ng2(quant%integral_points), quant%ng3(quant%integral_points), quant%gauss_w(quant%integral_points))
    call gauss_points(quant)

	
     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	CALL getarg(1, DATA_DIR)

	option%Nmin_leaf=40
	quant%mesh_normal=1
	! Refined_level=0
	
	! levelpara_control=0
	! ACA_tolerance_forward=1d-4
	option%tol_comp=1d-4
	! SVD_tolerance_factor=1d-4
	option%tol_Rdetect=3d-5
    ! Preset_level_butterfly=0
	option%nogeo=0
	quant%scaling=1d0
	quant%wavelength=1.0
	! Discret=0.05
	quant%RCS_static=2
    quant%RCS_Nsample=2000
    ! Optimizing_forward=0
    ! Fast_inverse=0
    ! Add_method_of_base_level=2
    ! quant%rank_approximate_para1=6.0
    ! quant%rank_approximate_para2=6.0
    ! quant%rank_approximate_para3=6.0
	option%tol_LS=1d-12
	! tfqmr_tolerance=1d-6
	option%tol_itersol=3d-3
	option%n_iter=1000
	option%tol_rand=5d-3
	! up_tolerance=1d-4
	! relax_lamda=1d0
	! SolvingMethod=1
	option%level_check=100
	! rank_tmp=7
	! schurinv=1
	! reducelevel_flag=0
	! directschur=1
	option%precon=DIRECT
	! verboselevel=2
	option%xyzsort=TM
	option%lnoBP=4000
	option%TwoLayerOnly=0
	quant%CFIE_alpha=1
    option%schulzorder=3
    option%schulzlevel=3000
	option%LRlevel=100
	option%ErrFillFull=0
	option%RecLR_leaf=ACA
	! call MKL_set_num_threads(NUM_Threads)    ! this overwrites omp_set_num_threads for MKL functions 
	
    ! call OMP_set_dynamic(.true.)
    !call OMP_set_nested(.true.)
	

    !*********************************************************

    quant%omiga=2*pi/quant%wavelength/sqrt(mu0*eps0)
    quant%wavenum=2*pi/quant%wavelength

   !***********************************************************************
	if(ptree%MyID==Main_ID)then							  
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'wavelength:',quant%wavelength
   write (*,*) ''
	endif		
   !***********************************************************************
	
	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID)write(*,*) "geometry modeling......"
    
	call geo_modeling_SURF(quant,ptree%Comm,DATA_DIR)
	
	! generate the list of points for clustering
	msh%Nunk=quant%Nunk
	allocate(msh%xyz(3,quant%Nunk))
    do edge=1, quant%Nunk
		msh%xyz(:,edge) = quant%xyz(:,quant%maxnode+edge)
    enddo	
	option%touch_para = 3* quant%minedgelength
	
    if(ptree%MyID==Main_ID)write(*,*) "modeling finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1

	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "constructing HODLR formatting......"
    call z_HODLR_structuring(ho_bf,option,msh,ker,z_element_Zmn_user,ptree)
	call z_BPlus_structuring(ho_bf,option,msh,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "HODLR formatting finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1
	! stop 
	
    !pause
    
    !call compression_test()
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "HODLR construction......"
    call z_HODLR_construction(ho_bf,option,stats,msh,ker,z_element_Zmn_user,ptree)
	! if(option%precon/=DIRECT)then
		! call copy_HOBF(ho_bf,ho_bf_copy)	
	! end if    
	if(ptree%MyID==Main_ID)write(*,*) "HODLR construction finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	! write(*,*)t2-t1
	
	if(option%precon/=NOPRECON)then								
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing......"
    call z_HODLR_Factorization(ho_bf,option,stats,ptree,msh)
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	end if

    if(ptree%MyID==Main_ID)write(*,*) "EM_solve......"
    call EM_solve_SURF(ho_bf,option,msh,quant,ptree,stats)
    if(ptree%MyID==Main_ID)write(*,*) "EM_solve finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	

	call delete_quant_EMSURF(quant)
	
	call z_delete_proctree(ptree)
	call z_delete_Hstat(stats)
	call z_delete_mesh(msh)
	call z_delete_kernelquant(ker)
	call z_delete_HOBF(ho_bf)
	
    if(ptree%MyID==Main_ID)write(*,*) "-------------------------------program end-------------------------------------"
	
	call blacs_exit(1)
	call MPI_Finalize(ierr)
	
    ! ! ! ! pause

end PROGRAM HODLR_BUTTERFLY_SOLVER_3D


subroutine EM_solve_SURF(ho_bf_inv,option,msh,quant,ptree,stats)
    use EMSURF_MODULE
    use z_HODLR_DEFS
	! use RCS_Bi
	! use RCS_Mono
	! use element_vinc
	use z_HODLR_Solve_Mul
	
    
    implicit none
    
    integer i, j, ii, jj, iii, jjj,ierr
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, n_iter_max, iter ,N_unk, N_unk_loc
    real(kind=8) theta, phi, dphi, rcs_V, rcs_H
    real T0
    real(kind=8) n1,n2,rtemp
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:)
	real(kind=8):: rel_error
	type(z_Hoption)::option
	type(z_hobf)::ho_bf_inv
	type(z_mesh)::msh
	type(quant_EMSURF)::quant
	type(z_proctree)::ptree
	type(z_Hstat)::stats
	complex(kind=8),allocatable:: current(:,:),voltage(:,:)
	
	
	if(option%ErrSol==1)then
		call z_hodlr_test_solve_error(ho_bf_inv,option,ptree,stats)
	endif
	
	
	! if(option%PRECON==DIRECT)then
		! msh%idxs = ho_bf_inv%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		! msh%idxe = ho_bf_inv%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)
	! else 
		! msh%idxs = ho_bf_for%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		! msh%idxe = ho_bf_for%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)	
	! endif
	
	! N_unk=msh%Nunk
	N_unk_loc = msh%idxe-msh%idxs+1
	
    if (quant%RCS_static==2) then
    
        theta=90
        phi=0
        
        allocate (current(N_unk_loc,2))
		Current=0
        allocate (voltage(N_unk_loc,2))

		
        !$omp parallel do default(shared) private(edge,value_Z)
        do edge=msh%idxs, msh%idxe
            call element_Vinc_VV_SURF(theta,phi,msh%new2old(edge),value_Z,quant)
			voltage(edge-msh%idxs+1,1)=value_Z
			call element_Vinc_HH_SURF(theta,phi,msh%new2old(edge),value_Z,quant)
            voltage(edge-msh%idxs+1,2)=value_Z
        enddo    
        !$omp end parallel do
        
        T0=secnds(0.0)
        		
		call z_hodlr_solution(ho_bf_inv,Current,Voltage,N_unk_loc,2,option,ptree,stats)
					
		
        if(ptree%MyID==Main_ID)write (*,*) ''
        if(ptree%MyID==Main_ID)write (*,*) 'Solving:',secnds(T0),'Seconds'
        if(ptree%MyID==Main_ID)write (*,*) ''

		T0=secnds(0.0)
        call RCS_bistatic_SURF(Current,msh,quant,ptree)
		
		! call current_node_patch_mapping('V',curr(:,1),msh)    		
		! call current_node_patch_mapping('H',curr(:,2),msh)          

        if(ptree%MyID==Main_ID)write (*,*) ''
        if(ptree%MyID==Main_ID)write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
        if(ptree%MyID==Main_ID)write (*,*) ''
		deallocate(Current)
		deallocate(Voltage)
	
    elseif (quant%RCS_static==1) then
    
        allocate (current(N_unk_loc,1))

        
        num_sample=quant%RCS_Nsample
		theta=90.
        dphi=180./num_sample
		allocate (b(N_unk_loc,num_sample+1))
		allocate (x(N_unk_loc,num_sample+1))        
		x=0
		
		
        if(ptree%MyID==Main_ID)open (100, file='bistaticH.out')

        n1=OMP_get_wtime()
		
        do j=0, num_sample 
            phi=j*dphi
			!$omp parallel do default(shared) private(edge,value_Z)
			do edge=msh%idxs, msh%idxe
				call element_Vinc_HH_SURF(theta,phi,msh%new2old(edge),value_Z,quant)
				b(edge-msh%idxs+1,j+1)=value_Z
			enddo    
			!$omp end parallel do
        enddo
		
		call z_hodlr_solution(ho_bf_inv,x,b,N_unk_loc,num_sample+1,option,ptree,stats)
			
		do j=0, num_sample 			
			phi=j*dphi
			
			Current(:,1)=x(:,j+1)
			
            call RCS_monostatic_HH_SURF(theta,phi,rcs_H,Current(:,1),msh,quant,ptree)
!             !$omp parallel do default(shared) private(i)
!             do i=1, N_unk
!                 current(i)=vectors_block(0)%vector(i,2)
!             enddo
!             !$omp end parallel do
!             call RCS_monostatic_HH_SURF(theta,phi,rcs_H)
            
            if(ptree%MyID==Main_ID)write (100,*) phi,rcs_H !,rcs_H
            
            ! deallocate (vectors_block)
            
        enddo
        
		n2 = OMP_get_wtime()
		stats%Time_Sol = stats%Time_Sol + n2-n1
		call MPI_ALLREDUCE(stats%Time_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
		
        if(ptree%MyID==Main_ID)then
			close(100)
			write (*,*) ''
			write (*,*) 'Solving:',rtemp,'Seconds'
			write (*,*) ''   
		endif
		
		call MPI_ALLREDUCE(stats%Flop_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
		if(ptree%MyID==Main_ID)write (*,'(A13Es14.2)') 'Solve flops:',rtemp	
	
		
		deallocate(Current)
		
    endif
        
    return
    
end subroutine EM_solve_SURF



