module EM_calculation
use RCS_Bi
use RCS_Mono
use element_vinc
use HODLR_Solve
contains 


subroutine EM_solve_SURF(ho_bf_for,ho_bf_inv,option,msh,ker,ptree,stats)
    
    use MODULE_FILE
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
					
		
        write (*,*) ''
        write (*,*) 'Solving:',secnds(T0),'Seconds'
        write (*,*) ''

		T0=secnds(0.0)
        call RCS_bistatic_SURF(Current,msh,ker,ptree)
		
		! call current_node_patch_mapping('V',curr(:,1),msh)    		
		! call current_node_patch_mapping('H',curr(:,2),msh)          

        write (*,*) ''
        write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
        write (*,*) ''
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
		
		deallocate(Current)
		
    endif
        
    return
    
end subroutine EM_solve_SURF



subroutine EM_solve_CURV(ho_bf_for,ho_bf_inv,option,msh,ker,ptree,stats)
    
    use MODULE_FILE
    ! use blas95
    implicit none
    
    integer i, j, ii, jj, iii, jjj, ierr
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, N_iter_max, iter, N_unk, N_unk_loc
    real*8 theta, phi, dphi, rcs_V, rcs_H
    real T0
    real*8 n1,n2,rtemp	
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:)
	real*8:: rel_error
	type(Hoption)::option
	type(mesh)::msh
	type(kernelquant)::ker
	type(proctree)::ptree
	type(hobf)::ho_bf_for,ho_bf_inv
	type(Hstat)::stats	
	complex(kind=8),allocatable:: current(:),voltage(:)
	

	if(option%PRECON==DIRECT)then
		msh%idxs = ho_bf_inv%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		msh%idxe = ho_bf_inv%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)
	else 
		! write(*,*)associated(ho_bf_for%levels(1)%BP_inverse),'dd' !%matrices_block(1)%N_p),'nima'
		msh%idxs = ho_bf_for%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		msh%idxe = ho_bf_for%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)	
	endif
	
	N_unk=msh%Nunk
	N_unk_loc = msh%idxe-msh%idxs+1
	
    if (ker%RCS_static==2) then
    
        phi=180d0
        
        allocate (current(N_unk_loc))
		Current=0
        allocate (voltage(N_unk_loc))
								  
        !$omp parallel do default(shared) private(edge,value_Z)
        do edge=msh%idxs, msh%idxe
            call element_Vinc_VV_CURV(phi,edge,value_Z,msh,ker)
            voltage(edge-msh%idxs+1)=value_Z
        enddo    
        !$omp end parallel do
        
        n1 = OMP_get_wtime()
        
		call HODLR_Solution(ho_bf_for,ho_bf_inv,Current,Voltage,N_unk_loc,1,option,ptree,stats)
		
		n2 = OMP_get_wtime()
		stats%Time_Sol = stats%Time_Sol + n2-n1
		call MPI_ALLREDUCE(stats%Time_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
		if(ptree%MyID==Main_ID)then
			write (*,*) ''
			write (*,*) 'Solving:',rtemp,'Seconds'
			write (*,*) ''
		endif
		
		T0=secnds(0.0)
        call RCS_bistatic_CURV(Current,msh,ker,ptree)

		if(ptree%MyID==Main_ID)then
			write (*,*) ''
			write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
			write (*,*) ''
		endif
    
		deallocate(current)
		deallocate(voltage)
	
    elseif (ker%RCS_static==1) then
    
        allocate (current(N_unk_loc))
        num_sample=ker%RCS_Nsample
        dphi=180./num_sample
        
		allocate (b(N_unk_loc,num_sample+1))
		allocate (x(N_unk_loc,num_sample+1))
		x=0
		
        if(ptree%MyID==Main_ID)open (100, file='RCS_monostatic.txt')

        n1=OMP_get_wtime()       
        do j=0, num_sample 
            phi=j*dphi
			!$omp parallel do default(shared) private(edge,value_Z)
			do edge=msh%idxs, msh%idxe
				call element_Vinc_VV_CURV(phi,edge,value_Z,msh,ker)
				b(edge-msh%idxs+1,j+1)=value_Z
			enddo    
			!$omp end parallel do
		enddo
		
		call HODLR_Solution(ho_bf_for,ho_bf_inv,x,b,N_unk_loc,num_sample+1,option,ptree,stats)
			
			
		do j=0, num_sample 	
			phi=j*dphi
			Current=x(:,j+1)
            call RCS_monostatic_VV_CURV(phi,rcs_V,Current,msh,ker,ptree)
!             !$omp parallel do default(shared) private(i)
!             do i=1, N_unk
!                 current(i)=vectors_block(0)%vector(i,2)
!             enddo
!             !$omp end parallel do
!             call RCS_monostatic_HH(theta,phi,rcs_H)
            
            if(ptree%MyID==Main_ID)write (100,*) j,phi,rcs_V !,rcs_H
            
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

		
	! call MPI_ALLREDUCE(stats%Time_RedistV,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	! if(ptree%MyID==Main_ID)write (*,*) '     Time_RedistV:', rtemp			
		
		deallocate(b)
		deallocate(x)
		deallocate(Current)
		
    endif
        
    return
    
end subroutine EM_solve_CURV

end module EM_calculation