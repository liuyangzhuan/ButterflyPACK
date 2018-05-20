module EM_calculation
use RCS_Bi
use RCS_Mono
use element_vinc
use HODLR_Solve
contains 


subroutine EM_solve_SURF(ho_bf_for,ho_bf_inv,option,msh,ker)
    
    use MODULE_FILE
    ! use blas95
    implicit none
    
    integer i, j, ii, jj, iii, jjj
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, N_iter_max, iter ,N_unk
    real*8 theta, phi, dphi, rcs_V, rcs_H
    real T0
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:)
	real*8:: rel_error
	type(Hoption)::option
	type(hobf)::ho_bf_for,ho_bf_inv
	type(mesh)::msh
	type(kernelquant)::ker
	
	complex(kind=8),allocatable:: current(:,:),voltage(:,:)
	N_unk=msh%Nunk
	
    if (ker%RCS_static==2) then
    
        theta=90
        phi=0
        
        allocate (current(N_unk,2))
        allocate (voltage(N_unk,2))

		
        !$omp parallel do default(shared) private(edge,value_Z)
        do edge=1, N_unk
            call element_Vinc_VV_SURF(theta,phi,edge,value_Z,msh,ker)
			voltage(edge,1)=value_Z
			call element_Vinc_HH_SURF(theta,phi,edge,value_Z,msh,ker)
            voltage(edge,2)=value_Z
        enddo    
        !$omp end parallel do
        
        T0=secnds(0.0)
        		
		call HODLR_Solution(ho_bf_for,ho_bf_inv,Current,Voltage,N_unk,2,option)
					
		
        write (*,*) ''
        write (*,*) 'Solving:',secnds(T0),'Seconds'
        write (*,*) ''

		T0=secnds(0.0)
        call RCS_bistatic_SURF(Current,msh,ker)
        

        write (*,*) ''
        write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
        write (*,*) ''
		deallocate(Current)
		deallocate(Voltage)
	
    elseif (ker%RCS_static==1) then
    
        allocate (current(N_unk,1))

        
        num_sample=ker%RCS_Nsample
		theta=90.
        dphi=180./num_sample
		allocate (b(N_unk,num_sample+1))
		allocate (x(N_unk,num_sample+1))        
		
		
		
        open (100, file='bistaticH.out')

        T0=secnds(0.0)        
        do j=0, num_sample 
            phi=j*dphi
			!$omp parallel do default(shared) private(edge,value_Z)
			do edge=1, N_unk
				call element_Vinc_HH_SURF(theta,phi,edge,value_Z,msh,ker)
				b(edge,j+1)=value_Z
			enddo    
			!$omp end parallel do
        enddo
		
		call HODLR_Solution(ho_bf_for,ho_bf_inv,x,b,N_unk,num_sample+1,option)
			
		do j=0, num_sample 			
			phi=j*dphi
			
			Current(:,1)=x(:,j+1)
			
            call RCS_monostatic_HH_SURF(theta,phi,rcs_H,Current(:,1),msh,ker)
!             !$omp parallel do default(shared) private(i)
!             do i=1, N_unk
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
    
end subroutine EM_solve_SURF



subroutine EM_solve_CURV(ho_bf_for,ho_bf_inv,option,msh,ker)
    
    use MODULE_FILE
    ! use blas95
    implicit none
    
    integer i, j, ii, jj, iii, jjj
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, N_iter_max, iter, N_unk
    real*8 theta, phi, dphi, rcs_V, rcs_H
    real T0
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:)
	real*8:: rel_error
	type(Hoption)::option
	type(mesh)::msh
	type(kernelquant)::ker
	
	type(hobf)::ho_bf_for,ho_bf_inv
	complex(kind=8),allocatable:: current(:),voltage(:)
	N_unk=msh%Nunk

    if (ker%RCS_static==2) then
    
        phi=180d0
        
        allocate (current(N_unk))
        allocate (voltage(N_unk))
								  
        !$omp parallel do default(shared) private(edge,value_Z)
        do edge=1, N_unk
            call element_Vinc_VV_CURV(phi,edge,value_Z,msh,ker)
            voltage(edge)=value_Z
        enddo    
        !$omp end parallel do
        
        T0=secnds(0.0)
        
		call HODLR_Solution(ho_bf_for,ho_bf_inv,Current,Voltage,N_unk,1,option)
		
        write (*,*) ''
        write (*,*) 'Solving:',secnds(T0),'Seconds'
        write (*,*) ''

		T0=secnds(0.0)
        call RCS_bistatic_CURV(Current,msh,ker)

        write (*,*) ''
        write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
        write (*,*) ''
    
		deallocate(current)
		deallocate(voltage)
	
    elseif (ker%RCS_static==1) then
    
        allocate (current(N_unk))
        num_sample=ker%RCS_Nsample
        dphi=180./num_sample
        
		allocate (b(N_unk,num_sample+1))
		allocate (x(N_unk,num_sample+1))
		
        open (100, file='RCS_monostatic.txt')

        T0=secnds(0.0)        
        do j=0, num_sample 
            phi=j*dphi
			!$omp parallel do default(shared) private(edge,value_Z)
			do edge=1, N_unk
				call element_Vinc_VV_CURV(phi,edge,value_Z,msh,ker)
				b(edge,j+1)=value_Z
			enddo    
			!$omp end parallel do
		enddo
		
		call HODLR_Solution(ho_bf_for,ho_bf_inv,x,b,N_unk,num_sample+1,option)
			
			
		do j=0, num_sample 	
			phi=j*dphi
			Current=x(:,j+1)
            call RCS_monostatic_VV_CURV(phi,rcs_V,Current,msh,ker)
!             !$omp parallel do default(shared) private(i)
!             do i=1, N_unk
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
    
end subroutine EM_solve_CURV

end module EM_calculation