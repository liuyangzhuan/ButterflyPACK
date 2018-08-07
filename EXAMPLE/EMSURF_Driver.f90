PROGRAM HODLR_BUTTERFLY_SOLVER_3D
    use MODULE_FILE
	! use geometry_model
	use H_structure
	use cascading_factorization
	use matrices_fill
	use omp_lib
	use MISC
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
	type(hobf)::ho_bf,ho_bf_copy
	type(kernelquant)::ker
	type(proctree)::ptree
	integer,allocatable:: groupmembers(:)	
	integer nmpi
	
	
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
    write(*,*) "HODLR_BUTTERFLY_SOLVER_3D"
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
	
	time_tmp = 0
	
	msh%Origins=(/0d0,0d0,0d0/)
    msh%integral_points=6
    allocate (msh%ng1(msh%integral_points), msh%ng2(msh%integral_points), msh%ng3(msh%integral_points), msh%gauss_w(msh%integral_points))
    call gauss_points(msh)

	
     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	CALL getarg(1, DATA_DIR)
	ker%Kernel = EMSURF	
	
	option%Nmin_leaf=40
	msh%mesh_normal=1
	! Refined_level=0
	para=0.001
	! levelpara_control=0
	! ACA_tolerance_forward=1d-4
	option%tol_SVD=1d-4
	! SVD_tolerance_factor=1d-4
	option%tol_Rdetect=3d-5
    ! Preset_level_butterfly=0
	msh%scaling=1d0
	ker%wavelength=1.0
	! Discret=0.05
	ker%RCS_static=2
    ker%RCS_Nsample=2000
    ! Optimizing_forward=0
    ! Fast_inverse=0
    ! Add_method_of_base_level=2
    ker%rank_approximate_para1=6.0
    ker%rank_approximate_para2=6.0
    ker%rank_approximate_para3=6.0
	option%tol_LS=1d-12
	! tfqmr_tolerance=1d-6
	option%tol_itersol=3d-3
	option%N_iter=1000
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
	option%xyzsort=3
	option%LnoBP=4000
	option%TwoLayerOnly=0
	ker%CFIE_alpha=1
	option%touch_para=3
    option%schulzorder=3
    option%schulzlevel=3000
	option%LRlevel=100
	option%ErrFillFull=0
	option%RecLR_leaf='A'	
	! call MKL_set_num_threads(NUM_Threads)    ! this overwrites omp_set_num_threads for MKL functions 
	
    ! call OMP_set_dynamic(.true.)
    !call OMP_set_nested(.true.)
	

    !*********************************************************

    ker%omiga=2*pi/ker%wavelength/sqrt(mu0*eps0)
    ker%wavenum=2*pi/ker%wavelength

   !***********************************************************************
	if(ptree%MyID==Main_ID)then							  
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'wavelength:',ker%wavelength
   write (*,*) ''
	endif		
   !***********************************************************************
	
	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID)write(*,*) "geometry modeling......"
    call geo_modeling_SURF(msh,ker,ptree)
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
    call matrices_filling(ho_bf,option,stats,msh,ker,element_Zmn_EMSURF,ptree)
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

    if(ptree%MyID==Main_ID)write(*,*) "EM_solve......"
    call EM_solve_SURF(ho_bf_copy,ho_bf,option,msh,ker,ptree,stats)
    if(ptree%MyID==Main_ID)write(*,*) "EM_solve finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	

    if(ptree%MyID==Main_ID)write(*,*) "-------------------------------program end-------------------------------------"

    ! ! ! ! pause

end PROGRAM HODLR_BUTTERFLY_SOLVER_3D





subroutine geo_modeling_SURF(msh,ker,ptree)

    use MODULE_FILE
	use misc
    implicit none
    
    integer i,j,ii,jj,iii,jjj
    integer intemp
    integer node, patch, edge, flag
    integer node1, node2
    integer node_temp(2)
	integer Dimn
    real*8 a(3),b(3),c(3),r0
	type(mesh)::msh
	type(kernelquant)::ker
	type(proctree)::ptree
	
	integer Maxedge
	
    Dimn=3
    
    open(11,file=trim(DATA_DIR)//'/node.geo')
    open(111,file=trim(DATA_DIR)//'/elem.geo')
    
    read(11,*)msh%maxnode
    read(111,*)msh%maxpatch
    Maxedge=msh%maxpatch*3/2
	
    
	
    allocate(msh%xyz(3,msh%maxnode+Maxedge))
    allocate(msh%node_of_patch(0:3,msh%maxpatch),msh%info_unk(0:6,maxedge+1000))
    allocate(msh%normal_of_patch(3,msh%maxpatch))
    
    
    !************msh%xyz****************
    i=1
    do while(i<=msh%maxnode)
        read(11,*)intemp,msh%xyz(1:3,i)
        msh%xyz(1:3,i)=msh%xyz(1:3,i)/msh%scaling
        i=i+1
    enddo
    close(11)
    
    i=1
    if (msh%mesh_normal==1) then
        do while(i<=msh%maxpatch)
            read(111,*)intemp,msh%node_of_patch(1:3,i)
            i=i+1 
        enddo
    elseif (msh%mesh_normal==-1) then
        do while(i<=msh%maxpatch)
            read(111,*)intemp,msh%node_of_patch(3,i),msh%node_of_patch(2,i),msh%node_of_patch(1,i)
            i=i+1 
        enddo
    endif
    close(111)
    
    !************msh%normal_of_patch****************
    
    !$omp parallel do default(shared) private(patch,i,a,b,c,r0)
    do patch=1,msh%maxpatch
        do i=1,3
            a(i)=(msh%xyz(i,msh%node_of_patch(2,patch))-msh%xyz(i,msh%node_of_patch(1,patch)))
            b(i)=(msh%xyz(i,msh%node_of_patch(3,patch))-msh%xyz(i,msh%node_of_patch(1,patch)))
        enddo
        call curl(a,b,c)
        r0=sqrt(c(1)**2+c(2)**2+c(3)**2)
        c(1)=c(1)/r0
        c(2)=c(2)/r0
        c(3)=c(3)/r0
        msh%normal_of_patch(1:3,patch)=c(1:3)	    
    enddo
    !$omp end parallel do
    
    !************msh%info_unk****************

    edge=0
    do i=1,msh%maxpatch-1
        do j=i+1,msh%maxpatch
            flag=0;node1=0;node2=0;iii=1
            do ii=1,3
                do jj=1,3
	     	         if(msh%node_of_patch(ii,i)==msh%node_of_patch(jj,j))then
                        flag=flag+1
                        node_temp(iii)=msh%node_of_patch(ii,i)
                        iii=iii+1
                    endif
                enddo
            enddo
            if(flag==2)then
                edge=edge+1
                if(node_temp(1)<node_temp(2)) then
                    msh%info_unk(1,edge)=node_temp(1)
                    msh%info_unk(2,edge)=node_temp(2)
                else
                    msh%info_unk(1,edge)=node_temp(2)
                    msh%info_unk(2,edge)=node_temp(1)
                endif
                msh%info_unk(3,edge)=i
                msh%info_unk(4,edge)=j       ! notice that : i<j  
                msh%info_unk(0,edge)=0
            endif
        enddo
    enddo
    
    Maxedge=edge
    
    !$omp parallel do default(shared) private(edge,node_temp,jj,iii,jjj)
    do edge=1,maxedge
	    node_temp(1)=0
	    node_temp(2)=0	    
	    do jj=3,4
             do iii=1,3
                 do jjj=1,2
        	            if(msh%node_of_patch(iii,msh%info_unk(jj,edge))==msh%info_unk(jjj,edge)) then
                           node_temp(jj-2)=node_temp(jj-2)+iii               
        	            endif
                 enddo
             enddo
         enddo
         node_temp(1)=6-node_temp(1)
         node_temp(2)=6-node_temp(2)
         msh%info_unk(5,edge)=msh%node_of_patch(node_temp(1),msh%info_unk(3,edge))
         msh%info_unk(6,edge)=msh%node_of_patch(node_temp(2),msh%info_unk(4,edge))
    enddo
    !$omp end parallel do
    
    node=msh%maxnode
    do edge=1, Maxedge
        node=node+1
        msh%info_unk(0,edge)=node
        do i=1,3
            msh%xyz(i,node)=1./2.*(msh%xyz(i,msh%info_unk(1,edge))+msh%xyz(i,msh%info_unk(2,edge)))
        enddo
    enddo
	
	msh%maxedgelength = 0
	do edge=1,Maxedge
		msh%maxedgelength = max(msh%maxedgelength,sqrt(sum(abs(msh%xyz(:,msh%info_unk(1,edge))-msh%xyz(:,msh%info_unk(2,edge)))**2)))
	end do	

	msh%minedgelength = BigValue
	do edge=1,Maxedge
		msh%minedgelength = min(msh%minedgelength,sqrt(sum(abs(msh%xyz(:,msh%info_unk(1,edge))-msh%xyz(:,msh%info_unk(2,edge)))**2)))
	end do	
	
	! write(*,*)	msh%xyz(1,1:100),sum(msh%xyz(1,:))
	! stop
	msh%Nunk = Maxedge

    if(ptree%MyID==Main_ID)write (*,*) ''
    if(ptree%MyID==Main_ID)write (*,*) 'Maxedge:',Maxedge
	if(ptree%MyID==Main_ID)write (*,*) 'minedgelength:',msh%minedgelength
	if(ptree%MyID==Main_ID)write (*,*) 'wavelength/minedgelength:',ker%wavelength/msh%minedgelength
	if(ptree%MyID==Main_ID)write (*,*) 'maxedgelength:',msh%maxedgelength
	if(ptree%MyID==Main_ID)write (*,*) 'wavelength/maxedgelength:',ker%wavelength/msh%maxedgelength

    if(ptree%MyID==Main_ID)write (*,*) '' 
    
    return
    
end subroutine geo_modeling_SURF




subroutine EM_solve_SURF(ho_bf_for,ho_bf_inv,option,msh,ker,ptree,stats)
    
    use MODULE_FILE
	use RCS_Bi
	use RCS_Mono
	use element_vinc
	use HODLR_Solve
	
    ! use blas95
    implicit none
    
    integer i, j, ii, jj, iii, jjj,ierr
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, N_iter_max, iter ,N_unk, N_unk_loc
    real*8 theta, phi, dphi, rcs_V, rcs_H
    real T0
    real*8 n1,n2,rtemp
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:)
	real*8:: rel_error
	type(Hoption)::option
	type(hobf)::ho_bf_for,ho_bf_inv
	type(mesh)::msh
	type(kernelquant)::ker
	type(proctree)::ptree
	type(Hstat)::stats
	complex(kind=8),allocatable:: current(:,:),voltage(:,:)
	
	
	if(option%ErrSol==1)then
		call HODLR_Test_Solve_error(ho_bf_for,ho_bf_inv,option,msh,ker,ptree,stats)
	endif
	
	
	if(option%PRECON==DIRECT)then
		msh%idxs = ho_bf_inv%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		msh%idxe = ho_bf_inv%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)
	else 
		msh%idxs = ho_bf_for%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		msh%idxe = ho_bf_for%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)	
	endif
	
	N_unk=msh%Nunk
	N_unk_loc = msh%idxe-msh%idxs+1
	
    if (ker%RCS_static==2) then
    
        theta=90
        phi=0
        
        allocate (current(N_unk_loc,2))
		Current=0
        allocate (voltage(N_unk_loc,2))

		
        !$omp parallel do default(shared) private(edge,value_Z)
        do edge=msh%idxs, msh%idxe
            call element_Vinc_VV_SURF(theta,phi,edge,value_Z,msh,ker)
			voltage(edge-msh%idxs+1,1)=value_Z
			call element_Vinc_HH_SURF(theta,phi,edge,value_Z,msh,ker)
            voltage(edge-msh%idxs+1,2)=value_Z
        enddo    
        !$omp end parallel do
        
        T0=secnds(0.0)
        		
		call HODLR_Solution(ho_bf_for,ho_bf_inv,Current,Voltage,N_unk_loc,2,option,ptree,stats)
					
		
        if(ptree%MyID==Main_ID)write (*,*) ''
        if(ptree%MyID==Main_ID)write (*,*) 'Solving:',secnds(T0),'Seconds'
        if(ptree%MyID==Main_ID)write (*,*) ''

		T0=secnds(0.0)
        call RCS_bistatic_SURF(Current,msh,ker,ptree)
		
		! call current_node_patch_mapping('V',curr(:,1),msh)    		
		! call current_node_patch_mapping('H',curr(:,2),msh)          

        if(ptree%MyID==Main_ID)write (*,*) ''
        if(ptree%MyID==Main_ID)write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
        if(ptree%MyID==Main_ID)write (*,*) ''
		deallocate(Current)
		deallocate(Voltage)
	
    elseif (ker%RCS_static==1) then
    
        allocate (current(N_unk_loc,1))

        
        num_sample=ker%RCS_Nsample
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
				call element_Vinc_HH_SURF(theta,phi,edge,value_Z,msh,ker)
				b(edge-msh%idxs+1,j+1)=value_Z
			enddo    
			!$omp end parallel do
        enddo
		
		call HODLR_Solution(ho_bf_for,ho_bf_inv,x,b,N_unk_loc,num_sample+1,option,ptree,stats)
			
		do j=0, num_sample 			
			phi=j*dphi
			
			Current(:,1)=x(:,j+1)
			
            call RCS_monostatic_HH_SURF(theta,phi,rcs_H,Current(:,1),msh,ker,ptree)
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

