module HODLR_Solve
use Butterfly_rightmultiply
use Butterfly_inversion
use Butterfly_compress_forward
contains 


subroutine HODLR_Solution(hobf_forward,hobf_inverse,x,b,Ns,num_vectors,option)
    
    use MODULE_FILE
    implicit none
    
    integer i, j, ii, jj, iii, jjj
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, N_iter_max, iter, Ns, num_vectors 
    real*8 theta, phi, dphi, rcs_V, rcs_H
    real T0
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:)
	real*8:: rel_error
	type(Hoption)::option
	
	! type(cascadingfactors)::cascading_factors_forward(:),cascading_factors_inverse(:)
	type(hobf)::hobf_forward,hobf_inverse
	complex(kind=8)::x(Ns,num_vectors),b(Ns,num_vectors)
	complex(kind=8),allocatable::r0_initial(:)

	if(option%precon/=DIRECT)then
        allocate(r0_initial(1:Ns))	
		do ii=1,Ns
		   r0_initial(ii)= random_complex_number()
		end do			
		
		do ii=1,num_vectors
			iter = 0
			call HODLR_Ztfqmr(option%precon,option%N_iter,Ns,b(:,ii),x(:,ii),option%tol_itersol,iter,r0_initial,hobf_forward,hobf_inverse)
		end do
		
		deallocate(r0_initial)
	else 			
		call MVM_Z_factorized(Ns,num_vectors,b,x,hobf_inverse)
	end if	
 
    return
    
end subroutine HODLR_Solution


  subroutine HODLR_Ztfqmr(precond,ntotal,nn,b,x,err,iter,r0_initial,hobf_forward,hobf_inverse)
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
	! type(cascadingfactors)::cascading_factors_forward(:),cascading_factors_inverse(:)
	type(hobf)::hobf_forward,hobf_inverse
	complex(kind=8)::r0_initial(:)
	integer precond
	
    itmax=iter

	! ! call MVM_Z_factorized(nn,1,b,bb,hobf_inverse)	
	! bb=b 
	call HODLR_ApplyPrecon(precond,nn,b,bb,hobf_inverse)	
	
	
    ! ! ! if (myid == main_id) then
       ! ! ! call cpu_time(tim1)
       ! ! ! open(unit=32,file='iterations.out',status='unknown')
    ! ! ! end if
    
    if (iter.eq.0) itmax=ntotal
    
    
    !  set initial values
    !
    d=cmplx(0.0_dp,0.0_dp,dp)
    ! write(*,*)'1'
	! call SmartMultifly(trans,nn,level_c,rowblock,1,x,r)    
    call MVM_Z_forward('N',nn,1,1,hobf_forward%Maxlevel+1,x,ytmp,hobf_forward)
	! ! call MVM_Z_factorized(nn,1,ytmp,r,hobf_inverse)
	! r=ytmp
	call HODLR_ApplyPrecon(precond,nn,ytmp,r,hobf_inverse)	
	
    r=bb-r !residual from the initial guess
    w=r
    yo=r
        ! ! write(*,*)'2'
    	! ! if(isnan(sum(abs(yo)**2)))then
			! ! write(*,*)'shitddd'
			! ! stop
		! ! end if		
	! call SmartMultifly(trans,nn,level_c,rowblock,1,yo,ayo)    
    call MVM_Z_forward('N',nn,1,1,hobf_forward%Maxlevel+1,yo,ytmp,hobf_forward)
	! ! call MVM_Z_factorized(nn,1,ytmp,ayo,hobf_inverse)
	! ayo=ytmp
	call HODLR_ApplyPrecon(precond,nn,ytmp,ayo,hobf_inverse)	
	
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
	   call MVM_Z_forward('N',nn,1,1,hobf_forward%Maxlevel+1,ye,ytmp,hobf_forward)
	   ! ! call MVM_Z_factorized(nn,1,ytmp,aye,hobf_inverse)
	   ! aye=ytmp
		call HODLR_ApplyPrecon(precond,nn,ytmp,aye,hobf_inverse)	
	
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
		  call MVM_Z_forward('N',nn,1,1,hobf_forward%Maxlevel+1,x,ytmp,hobf_forward)
		  ! ! call MVM_Z_factorized(nn,1,ytmp,r,hobf_inverse)
		  ! r=ytmp
		  call HODLR_ApplyPrecon(precond,nn,ytmp,r,hobf_inverse)	
	
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
		call MVM_Z_forward('N',nn,1,1,hobf_forward%Maxlevel+1,yo,ytmp,hobf_forward)
		! ! call MVM_Z_factorized(nn,1,ytmp,ayo,hobf_inverse)
		! ayo=ytmp
		call HODLR_ApplyPrecon(precond,nn,ytmp,ayo,hobf_inverse)	
	
       !MAGIC
       v=ayo+beta*( aye+beta*v )
    enddo iters
    ! write(*,*)'6'
    ! call SmartMultifly(trans,nn,level_c,rowblock,1,x,r)
    call MVM_Z_forward('N',nn,1,1,hobf_forward%Maxlevel+1,x,ytmp,hobf_forward)
	
	! ! call MVM_Z_factorized(nn,1,ytmp,r,hobf_inverse)	 
	! r=ytmp
	call HODLR_ApplyPrecon(precond,nn,ytmp,r,hobf_inverse)	
	
	
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
  end subroutine HODLR_Ztfqmr



  subroutine HODLR_ApplyPrecon(precond,nn,x,y,hobf_inverse)
    implicit none
	integer nn
	complex(kind=dp),dimension(1:nn)::x,y
	integer precond
	type(hobf),optional::hobf_inverse
	if(precond==NOPRECON)then
		y=x
	else if (precond==SAIPRECON)then 
		call MVM_Z_factorized(nn,1,x,y,hobf_inverse)	 
	endif
	end subroutine HODLR_ApplyPrecon

end module HODLR_Solve