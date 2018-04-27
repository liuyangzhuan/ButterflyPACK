module EM_calculation
use RCS_Bi
use RCS_Mono
use element_vinc
use HODLR_Solve
contains 

subroutine EM_calculating()
    
    use MODULE_FILE
    implicit none
	if(Kernel==EMCURV)then
		call EM_calculating_CURV()
	elseif(Kernel==EMSURF)then
		call EM_calculating_SURF()
	else
		write(*,*)'unknown Kernel for EM_calculating'
		stop
	endif	
	
end subroutine EM_calculating  


subroutine EM_calculating_SURF()
    
    use MODULE_FILE
    ! use blas95
    implicit none
    
    integer i, j, ii, jj, iii, jjj
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, N_iter_max, iter 
    real*8 theta, phi, dphi, rcs_V, rcs_H
    real T0
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:)
	real*8:: rel_error
	
    if (Static==2) then
    
        theta=90
        phi=0
        
        allocate (current(Maxedge))
        allocate (voltage(Maxedge))

		allocate(current2com(Maxedge,2))
		
        !$omp parallel do default(shared) private(edge,value_Z)
        do edge=1, Maxedge
            call element_Vinc_VV_SURF(theta,phi,edge,value_Z)
            voltage(edge)=value_Z
        enddo    
        !$omp end parallel do
        
        T0=secnds(0.0)
        		
		call HODLR_Solution(cascading_factors_copy,cascading_factors,Current,Voltage,Maxedge,1)

		current2com(:,1) = Current
		
		
        !$omp parallel do default(shared) private(edge,value_Z)
        do edge=1, Maxedge
            call element_Vinc_HH_SURF(theta,phi,edge,value_Z)
            voltage(edge)=value_Z
        enddo    
        !$omp end parallel do
        
        T0=secnds(0.0)
        

		call HODLR_Solution(cascading_factors_copy,cascading_factors,Current,Voltage,Maxedge,1)		
		
		current2com(:,2) = Current		
		
        write (*,*) ''
        write (*,*) 'Solving:',secnds(T0),'Seconds'
        write (*,*) ''

		T0=secnds(0.0)
        call RCS_bistatic_SURF()
        

        write (*,*) ''
        write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
        write (*,*) ''
		deallocate(Current)
		deallocate(Voltage)
	
    elseif (Static==1) then
    
        allocate (current(Maxedge))

        
        num_sample=RCS_sample
		theta=90.
        dphi=180./num_sample
		allocate (b(Maxedge,num_sample+1))
		allocate (x(Maxedge,num_sample+1))        
		
		
		
        open (100, file='bistaticH.out')

        T0=secnds(0.0)        
        do j=0, num_sample 
            phi=j*dphi
			!$omp parallel do default(shared) private(edge,value_Z)
			do edge=1, Maxedge
				call element_Vinc_HH_SURF(theta,phi,edge,value_Z)
				b(edge,j+1)=value_Z
			enddo    
			!$omp end parallel do
        enddo
		
		call HODLR_Solution(cascading_factors_copy,cascading_factors,x,b,Maxedge,num_sample+1)
			
		do j=0, num_sample 			
			phi=j*dphi
			
			Current=x(:,j+1)
			
            call RCS_monostatic_HH_SURF(theta,phi,rcs_H)
!             !$omp parallel do default(shared) private(i)
!             do i=1, Maxedge
!                 current(i)=vectors_block(0)%vector(i,2)
!             enddo
!             !$omp end parallel do
!             call RCS_monostatic_HH_SURF(theta,phi,rcs_H)
            
            write (100,*) phi,rcs_H !,rcs_H
            
            ! deallocate (vectors_block)
            
        enddo
        
        close(100)
        write (*,*) ''
        write (*,*) 'Solving:',secnds(T0),'Seconds'
        write (*,*) ''   
		deallocate(Current)
		
    endif
        
    return
    
end subroutine EM_calculating_SURF



subroutine EM_calculating_CURV()
    
    use MODULE_FILE
    ! use blas95
    implicit none
    
    integer i, j, ii, jj, iii, jjj
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, N_iter_max, iter 
    real*8 theta, phi, dphi, rcs_V, rcs_H
    real T0
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:)
	real*8:: rel_error
	
    if (Static==2) then
    
        phi=180d0
        
        allocate (current(Maxedge))
        allocate (voltage(Maxedge))

		allocate(current2com(Maxedge,1))								  
        !$omp parallel do default(shared) private(edge,value_Z)
        do edge=1, Maxedge
            call element_Vinc_VV_CURV(phi,edge,value_Z)
            voltage(edge)=value_Z
        enddo    
        !$omp end parallel do
        
        T0=secnds(0.0)
        
		call HODLR_Solution(cascading_factors_copy,cascading_factors,Current,Voltage,Maxedge,1)
		
		current2com(:,1) = Current		
		
        write (*,*) ''
        write (*,*) 'Solving:',secnds(T0),'Seconds'
        write (*,*) ''

		T0=secnds(0.0)
        call RCS_bistatic_CURV()

        write (*,*) ''
        write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
        write (*,*) ''
    
		deallocate(current)
		deallocate(voltage)
	
    elseif (Static==1) then
    
        allocate (current(Maxedge))
        num_sample=RCS_sample
        dphi=180./num_sample
        
		allocate (b(Maxedge,num_sample+1))
		allocate (x(Maxedge,num_sample+1))
		
        open (100, file='RCS_monostatic.txt')

        T0=secnds(0.0)        
        do j=0, num_sample 
            phi=j*dphi
			!$omp parallel do default(shared) private(edge,value_Z)
			do edge=1, Maxedge
				call element_Vinc_VV_CURV(phi,edge,value_Z)
				b(edge,j+1)=value_Z
			enddo    
			!$omp end parallel do
		enddo
		
		call HODLR_Solution(cascading_factors_copy,cascading_factors,x,b,Maxedge,num_sample+1)
			
			
		do j=0, num_sample 	
			phi=j*dphi
			Current=x(:,j+1)
            call RCS_monostatic_VV_CURV(phi,rcs_V)
!             !$omp parallel do default(shared) private(i)
!             do i=1, Maxedge
!                 current(i)=vectors_block(0)%vector(i,2)
!             enddo
!             !$omp end parallel do
!             call RCS_monostatic_HH(theta,phi,rcs_H)
            
            write (100,*) j,phi,rcs_V !,rcs_H
            
            ! deallocate (vectors_block)
            
        enddo
        
        close(100)
        write (*,*) ''
        write (*,*) 'Solving:',secnds(T0),'Seconds'
        write (*,*) ''     

		deallocate(b)
		deallocate(x)
		deallocate(Current)
		
    endif
        
    return
    
end subroutine EM_calculating_CURV

end module EM_calculation