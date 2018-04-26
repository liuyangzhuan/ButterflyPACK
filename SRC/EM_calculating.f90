module EM_calculation
use Butterfly_rightmultiply
use Butterfly_inversion
use Butterfly_compress_forward
use RCS_Bi
use RCS_Mono
use element_vinc
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
    complex(kind=8),allocatable:: Voltage_pre(:)
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
        

		if(preconditioner==1)then
			call initialize_r0(Maxedge)
			allocate(Voltage_pre(Maxedge))
			call MVM_Z_factorized(Maxedge,Voltage,Voltage_pre)
			N_iter_max = 100
			iter = 0
			rel_error = tfqmr_tolerance_solving
			call ztfqmr_MLMDA(N_iter_max,Maxedge,Voltage_pre,Current,rel_error,iter)
			deallocate(Voltage_pre)
			deallocate(r0_initial)
		else 			
			call MVM_Z_factorized(Maxedge,Voltage,Current)
		end if		
		current2com(:,1) = Current
		
		
        !$omp parallel do default(shared) private(edge,value_Z)
        do edge=1, Maxedge
            call element_Vinc_HH_SURF(theta,phi,edge,value_Z)
            voltage(edge)=value_Z
        enddo    
        !$omp end parallel do
        
        T0=secnds(0.0)
        

		if(preconditioner==1)then
			call initialize_r0(Maxedge)
			allocate(Voltage_pre(Maxedge))
			call MVM_Z_factorized(Maxedge,Voltage,Voltage_pre)
			N_iter_max = 100
			iter = 0
			rel_error = tfqmr_tolerance_solving
			call ztfqmr_MLMDA(N_iter_max,Maxedge,Voltage_pre,Current,rel_error,iter)
			deallocate(Voltage_pre)
			deallocate(r0_initial)
		else 			
			call MVM_Z_factorized(Maxedge,Voltage,Current)
		end if		
		current2com(:,2) = Current		
		
        write (*,*) ''
        write (*,*) 'Solving:',secnds(T0),'Seconds'
        write (*,*) ''

		T0=secnds(0.0)
        call RCS_bistatic_SURF()
        

        write (*,*) ''
        write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
        write (*,*) ''
    
    elseif (Static==1) then
    
        allocate (current(Maxedge))
        allocate (voltage(Maxedge))
        
        
        num_sample=RCS_sample
		theta=90.
        dphi=180./num_sample
        
        open (100, file='bistaticH.out')

        T0=secnds(0.0)        
        do j=0, num_sample 
        

            phi=j*dphi
			!$omp parallel do default(shared) private(edge,value_Z)
			do edge=1, Maxedge
				call element_Vinc_HH_SURF(theta,phi,edge,value_Z)
				voltage(edge)=value_Z
			enddo    
			!$omp end parallel do
            
			if(preconditioner==1)then
				call initialize_r0(Maxedge)				
				allocate(Voltage_pre(Maxedge))
				call MVM_Z_factorized(Maxedge,Voltage,Voltage_pre)
				N_iter_max = 100
				iter = 0
				rel_error = tfqmr_tolerance_solving
				call ztfqmr_MLMDA(N_iter_max,Maxedge,Voltage_pre,Current,rel_error,iter)
				deallocate(Voltage_pre)
				deallocate(r0_initial)
			else 			
				call MVM_Z_factorized(Maxedge,Voltage,Current)
			end if
			
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
    complex(kind=8),allocatable:: Voltage_pre(:)
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
        

		if(preconditioner==1)then
			call initialize_r0(Maxedge)
			allocate(Voltage_pre(Maxedge))
			call MVM_Z_factorized(Maxedge,Voltage,Voltage_pre)
			N_iter_max = 100
			iter = 0
			rel_error = tfqmr_tolerance_solving
			call ztfqmr_MLMDA(N_iter_max,Maxedge,Voltage_pre,Current,rel_error,iter)
			deallocate(Voltage_pre)
			deallocate(r0_initial)
		else 			
			call MVM_Z_factorized(Maxedge,Voltage,Current)
		end if		
		current2com(:,1) = Current		
		
        write (*,*) ''
        write (*,*) 'Solving:',secnds(T0),'Seconds'
        write (*,*) ''

		T0=secnds(0.0)
        call RCS_bistatic_CURV()

        write (*,*) ''
        write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
        write (*,*) ''
    
    elseif (Static==1) then
    
        allocate (current(Maxedge))
        allocate (voltage(Maxedge))
        
        
        num_sample=RCS_sample
        dphi=180./num_sample
        
        open (100, file='RCS_monostatic.txt')

        T0=secnds(0.0)        
        do j=0, num_sample 
        

            phi=j*dphi
			!$omp parallel do default(shared) private(edge,value_Z)
			do edge=1, Maxedge
				call element_Vinc_VV_CURV(phi,edge,value_Z)
				voltage(edge)=value_Z
			enddo    
			!$omp end parallel do
            
			if(preconditioner==1)then
				call initialize_r0(Maxedge)				
				allocate(Voltage_pre(Maxedge))
				call MVM_Z_factorized(Maxedge,Voltage,Voltage_pre)
				N_iter_max = 100
				iter = 0
				rel_error = tfqmr_tolerance_solving
				call ztfqmr_MLMDA(N_iter_max,Maxedge,Voltage_pre,Current,rel_error,iter)
				deallocate(Voltage_pre)
				deallocate(r0_initial)
			else 			
				call MVM_Z_factorized(Maxedge,Voltage,Current)
			end if
			
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
    endif
        
    return
    
end subroutine EM_calculating_CURV



  subroutine ztfqmr_MLMDA(ntotal,nn,b,x,err,iter)
    implicit none
	integer level_c,rowblock
    integer,intent(in)::ntotal
    integer::iter,itmax,it,nn	
    complex(kind=dp),dimension(1:nn)::x,bb,b,ytmp
    real(kind=dp)::err,rerr
    complex(kind=dp),dimension(1:nn)::w,yo,ayo,ye,aye,r,d,v
    real(kind=dp)::ta,we,cm
    complex(kind=dp)::we_local,we_sum,rerr_local,rerr_sum,err_local,err_sum
    complex(kind=dp)::ta_local,ta_sum,bmag_local,bmag_sum1,dumb_ali(6)
    complex(kind=dp)::etha,rho,rho_local,amgis,amgis_local,ahpla,dum,dum_local,beta
    real(kind=dp)::bmag
    real(kind=dp)::tim1,tim2
    integer::kk,ll
    ! Variables for storing current
    integer::count_unk,srcbox_id
    complex(kind=dp),dimension(:),allocatable::curr_coefs_dum_local,curr_coefs_dum_global
    real(kind=dp)::mem_est
	character:: trans

    itmax=iter

    ! ! ! if (myid == main_id) then
       ! ! ! call cpu_time(tim1)
       ! ! ! open(unit=32,file='iterations.out',status='unknown')
    ! ! ! end if
    
    if (iter.eq.0) itmax=ntotal
    bb=b 
    
    !  set initial values
    !
    d=cmplx(0.0_dp,0.0_dp,dp)
    ! write(*,*)'1'
	! call SmartMultifly(trans,nn,level_c,rowblock,1,x,r)    
    call MVM_Z_forward(nn,x,ytmp,cascading_factors_copy)
	call MVM_Z_factorized(nn,ytmp,r)
	
    r=bb-r !residual from the initial guess
    w=r
    yo=r
        ! ! write(*,*)'2'
    	! ! if(isnan(sum(abs(yo)**2)))then
			! ! write(*,*)'shitddd'
			! ! stop
		! ! end if		
	! call SmartMultifly(trans,nn,level_c,rowblock,1,yo,ayo)    
    call MVM_Z_forward(nn,yo,ytmp,cascading_factors_copy)
	call MVM_Z_factorized(nn,ytmp,ayo)
	
    v=ayo
    we=0.0_dp
    etha=cmplx(0.0_dp,0.0_dp,dp)
    
    ta_local=dot_product(r,r)
    rho_local=dot_product(r0_initial,r)
    bmag_local=dot_product(bb,bb)
    
    dumb_ali(1:3)=(/ta_local,rho_local,bmag_local/)
    ! call MPI_ALLREDUCE(dumb_ali(1:3),dumb_ali(4:6),3,MPI_DOUBLE_COMPLEX,&
         ! MPI_SUM,MPI_Comm_world,ierr)
	dumb_ali(4:6) = dumb_ali(1:3)	 
    ta_sum=dumb_ali(4);rho=dumb_ali(5);bmag_sum1=dumb_ali(6)
    ta=sqrt(abs(ta_sum))
    bmag=sqrt(abs(bmag_sum1))
    rerr=ta/bmag
    
    iters: do it=1,itmax
       amgis_local=dot_product(r0_initial,v)
       ! ! call MPI_ALLREDUCE(amgis_local,amgis,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
            ! ! MPI_Comm_world,ierr)
	   amgis = 	amgis_local
       ahpla=rho/amgis
       ye=yo-ahpla*v
           ! write(*,*)'3'
       ! call SmartMultifly(trans,nn,level_c,rowblock,1,ye,aye)
	   call MVM_Z_forward(nn,ye,ytmp,cascading_factors_copy)
	   call MVM_Z_factorized(nn,ytmp,aye)
	
       !  start odd (2n-1) m loop
       d=yo+(we*we*etha/ahpla)*d
       w=w-ahpla*ayo
       we_local=dot_product(w,w)
       ! call MPI_ALLREDUCE(we_local,we_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
            ! MPI_Comm_world,ierr)
		we_sum = we_local	
       we=sqrt(abs(we_sum))/ta
       cm=1.0d0/sqrt(1.0d0+we*we)
       ta=ta*we*cm
       etha=ahpla*cm*cm
       x=x+etha*d
       !  check if the result has converged.
       !a        if (err*bmag .gt. ta*sqrt(2.*it)) then
       !
       !  start even (2n)  m loop
       d=ye+(we*we*etha/ahpla)*d
       w=w-ahpla*aye
       we_local=dot_product(w,w)
       ! call MPI_ALLREDUCE(we_local,we_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
            ! MPI_Comm_world,ierr)
		we_sum = we_local	
       we=sqrt(abs(we_sum))/ta
       cm=1.0d0/sqrt(1.0d0+we*we)
       ta=ta*we*cm
       etha=ahpla*cm*cm
       x=x+etha*d
      
       
       !  check if the result has converged.
       if (mod(it,1)==0 .or. rerr<1.0_dp*err) then
    ! write(*,*)'4'
		  ! call SmartMultifly(trans,nn,level_c,rowblock,1,x,r)
		  call MVM_Z_forward(nn,x,ytmp,cascading_factors_copy)
		  call MVM_Z_factorized(nn,ytmp,r)
	
          r=bb-r
          rerr_local=dot_product(r,r)
          ! call MPI_ALLREDUCE(rerr_local,rerr_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
               ! MPI_Comm_world,ierr)
			rerr_sum = rerr_local   
          rerr=sqrt(abs(rerr_sum))/bmag
          
          ! ! if (myid==main_id) then
             print*,'# ofiter,error:',it,rerr
             ! write(32,*)'# ofiter,error:',it,rerr ! iterations file
          ! ! end if
          
          if (err > rerr) then
             err=rerr
             iter=it

             ! ! ! if (myid == main_id) then
                ! ! ! print*,'Total number of iterations and achieved residual :::',it,err
             ! ! ! end if
             
             return
          endif
       end if
       !  make preparations for next iteration
       dum_local=dot_product( r0_initial,w)
       ! call MPI_ALLREDUCE(dum_local,dum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,&
            ! MPI_Comm_world,ierr)
		dum = dum_local	
       beta=dum/rho
       rho=dum
       yo=w+beta*ye
           ! write(*,*)'5'
       ! call SmartMultifly(trans,nn,level_c,rowblock,1,yo,ayo)
		call MVM_Z_forward(nn,yo,ytmp,cascading_factors_copy)
		call MVM_Z_factorized(nn,ytmp,ayo)
	
       !MAGIC
       v=ayo+beta*( aye+beta*v )
    enddo iters
    ! write(*,*)'6'
    ! call SmartMultifly(trans,nn,level_c,rowblock,1,x,r)
    call MVM_Z_forward(nn,x,ytmp,cascading_factors_copy)
	call MVM_Z_factorized(nn,ytmp,r)	   
    !MAGIC
    r=bb-r
    err_local=dot_product(r,r)
    ! call MPI_ALLREDUCE(err_local,err_sum,1,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_Comm_world,&
         ! ierr)
	err_sum = err_local	 
    err=sqrt(abs(err_sum))/bmag
    iter=itmax



   print*,'Iterative solver is terminated without convergence!!!',it,err
   stop
	
    return
  end subroutine ztfqmr_MLMDA




end module EM_calculation