#include "HODLR_config.fi"
module z_HODLR_Solve_Mul

use z_Bplus_compress


contains 


subroutine HODLR_Solution(hobf_inverse,x,b,Ns_loc,num_vectors,option,ptree,stats)
    
    use z_HODLR_DEFS
    implicit none
    
    integer i, j, ii, jj, iii, jjj
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, n_iter_max, iter, Ns_loc, num_vectors 
    real(kind=8) theta, phi, dphi, rcs_V, rcs_H
    real T0
	real(kind=8):: rel_error
	type(Hoption)::option
	type(proctree)::ptree
	type(Hstat)::stats
	! type(cascadingfactors)::cascading_factors_forward(:),cascading_factors_inverse(:)
	type(hobf)::hobf_inverse
	DT::x(Ns_loc,num_vectors),b(Ns_loc,num_vectors)
	DT,allocatable::r0_initial(:)
	real(kind=8) n1,n2,rtemp

	n1 = OMP_get_wtime()
	
	if(option%precon/=DIRECT)then
        allocate(r0_initial(1:Ns_loc))	
		do ii=1,Ns_loc
		   call random_dp_number(r0_initial(ii))
		end do			
		
		do ii=1,num_vectors
			iter = 0
			rel_error = option%tol_itersol
			call HODLR_Ztfqmr(option%precon,option%n_iter,Ns_loc,b(:,ii),x(:,ii),rel_error,iter,r0_initial,hobf_inverse,ptree,stats)
		end do
		
		deallocate(r0_initial)
	else 			
		call MVM_Z_factorized(Ns_loc,num_vectors,b,x,hobf_inverse,ptree,stats)
	end if	
 
 
	n2 = OMP_get_wtime()
	stats%Time_Sol = stats%Time_Sol + n2-n1
 
    return
    
end subroutine HODLR_Solution


  subroutine HODLR_Ztfqmr(precond,ntotal,nn_loc,b,x,err,iter,r0_initial,hobf_inverse,ptree,stats)
    implicit none
	integer level_c,rowblock,ierr
    integer,intent(in)::ntotal
    integer::iter,itmax,it,nn_loc	
    DT,dimension(1:nn_loc)::x,bb,b,ytmp
    real(kind=dp)::err,rerr
    DT,dimension(1:nn_loc)::w,yo,ayo,ye,aye,r,d,v
    real(kind=dp)::ta,we,cm
    DT::we_local,we_sum,rerr_local,rerr_sum,err_local,err_sum
    DT::ta_local,ta_sum,bmag_local,bmag_sum1,dumb_ali(6)
    DT::etha,rho,rho_local,amgis,amgis_local,ahpla,dum,dum_local,beta
    real(kind=dp)::bmag
    real(kind=dp)::tim1,tim2
    integer::kk,ll
    ! Variables for storing current
    integer::count_unk,srcbox_id
    DT,dimension(:),allocatable::curr_coefs_dum_local,curr_coefs_dum_global
    real(kind=dp)::mem_est
	character:: trans
	type(Hstat)::stats
		
	! type(cascadingfactors)::cascading_factors_forward(:),cascading_factors_inverse(:)
	type(hobf)::hobf_inverse
	DT::r0_initial(:)
	integer precond
	type(proctree)::ptree
	
    itmax=iter

	call HODLR_ApplyPrecon(precond,nn_loc,b,bb,ptree,hobf_inverse,stats)	
	
	
    ! ! ! if (myid == main_id) then
       ! ! ! call cpu_time(tim1)
       ! ! ! open(unit=32,file='iterations.out',status='unknown')
    ! ! ! end if
    
    if (iter.eq.0) itmax=ntotal
    
    
    !  set initial values
    !
    d=0d0
    ! write(*,*)'1'
	! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,x,r)    
    call MVM_Z_forward('N',nn_loc,1,1,hobf_inverse%Maxlevel+1,x,ytmp,hobf_inverse,ptree,stats)
	stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
	call HODLR_ApplyPrecon(precond,nn_loc,ytmp,r,ptree,hobf_inverse,stats)	
	
    r=bb-r !residual from the initial guess
    w=r
    yo=r
        ! ! write(*,*)'2'
    	! ! if(isnan(sum(abs(yo)**2)))then
			! ! write(*,*)'shitddd'
			! ! stop
		! ! end if		
	! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,yo,ayo)    
    call MVM_Z_forward('N',nn_loc,1,1,hobf_inverse%Maxlevel+1,yo,ytmp,hobf_inverse,ptree,stats)
	stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
	call HODLR_ApplyPrecon(precond,nn_loc,ytmp,ayo,ptree,hobf_inverse,stats)	
	
    v=ayo
    we=0.0_dp
    etha=0d0
    
    ta_local=dot_product(r,r)
    rho_local=dot_product(r0_initial,r)
    bmag_local=dot_product(bb,bb)
    
    dumb_ali(1:3)=(/ta_local,rho_local,bmag_local/)
    call MPI_ALLREDUCE(dumb_ali(1:3),dumb_ali(4:6),3,MPI_DT,&
         MPI_SUM,ptree%Comm,ierr)
    ta_sum=dumb_ali(4);rho=dumb_ali(5);bmag_sum1=dumb_ali(6)
    ta=sqrt(abs(ta_sum))
    bmag=sqrt(abs(bmag_sum1))
    rerr=ta/bmag
    
    iters: do it=1,itmax
       amgis_local=dot_product(r0_initial,v)
       call MPI_ALLREDUCE(amgis_local,amgis,1,MPI_DT,MPI_SUM,&
            ptree%Comm,ierr)
       ahpla=rho/amgis
       ye=yo-ahpla*v
           ! write(*,*)'3'
       ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,ye,aye)
	   call MVM_Z_forward('N',nn_loc,1,1,hobf_inverse%Maxlevel+1,ye,ytmp,hobf_inverse,ptree,stats)
	   stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
	   call HODLR_ApplyPrecon(precond,nn_loc,ytmp,aye,ptree,hobf_inverse,stats)	
	
       !  start odd (2n-1) m loop
       d=yo+(we*we*etha/ahpla)*d
       w=w-ahpla*ayo
       we_local=dot_product(w,w)
       call MPI_ALLREDUCE(we_local,we_sum,1,MPI_DT,MPI_SUM,&
            ptree%Comm,ierr)	
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
       call MPI_ALLREDUCE(we_local,we_sum,1,MPI_DT,MPI_SUM,&
            ptree%Comm,ierr)	
       we=sqrt(abs(we_sum))/ta
       cm=1.0d0/sqrt(1.0d0+we*we)
       ta=ta*we*cm
       etha=ahpla*cm*cm
       x=x+etha*d
      
       
       !  check if the result has converged.
       if (mod(it,1)==0 .or. rerr<1.0_dp*err) then
    ! write(*,*)'4'
		  ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,x,r)
		  call MVM_Z_forward('N',nn_loc,1,1,hobf_inverse%Maxlevel+1,x,ytmp,hobf_inverse,ptree,stats)
		  stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
		  call HODLR_ApplyPrecon(precond,nn_loc,ytmp,r,ptree,hobf_inverse,stats)	
	
          r=bb-r
          rerr_local=dot_product(r,r)
          call MPI_ALLREDUCE(rerr_local,rerr_sum,1,MPI_DT,MPI_SUM,&
               ptree%Comm,ierr)
          rerr=sqrt(abs(rerr_sum))/bmag
          
          if (ptree%MyID==Main_ID) then
             print*,'# ofiter,error:',it,rerr
             ! write(32,*)'# ofiter,error:',it,rerr ! iterations file
          end if
          
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
       call MPI_ALLREDUCE(dum_local,dum,1,MPI_DT,MPI_SUM,&
            ptree%Comm,ierr)
       beta=dum/rho
       rho=dum
       yo=w+beta*ye
           ! write(*,*)'5'
       ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,yo,ayo)
		call MVM_Z_forward('N',nn_loc,1,1,hobf_inverse%Maxlevel+1,yo,ytmp,hobf_inverse,ptree,stats)
		stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
		call HODLR_ApplyPrecon(precond,nn_loc,ytmp,ayo,ptree,hobf_inverse,stats)	
	
       !MAGIC
       v=ayo+beta*( aye+beta*v )
    enddo iters
    ! write(*,*)'6'
    ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,x,r)
    call MVM_Z_forward('N',nn_loc,1,1,hobf_inverse%Maxlevel+1,x,ytmp,hobf_inverse,ptree,stats)
	stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
	call HODLR_ApplyPrecon(precond,nn_loc,ytmp,r,ptree,hobf_inverse,stats)	
	
	
    !MAGIC
    r=bb-r
    err_local=dot_product(r,r)
    call MPI_ALLREDUCE(err_local,err_sum,1,MPI_DT,MPI_SUM,ptree%Comm,&
         ierr)
    err=sqrt(abs(err_sum))/bmag
    iter=itmax



   print*,'Iterative solver is terminated without convergence!!!',it,err
   stop
	
    return
  end subroutine HODLR_Ztfqmr



  subroutine HODLR_ApplyPrecon(precond,nn_loc,x,y,ptree,hobf_inverse,stats)
    implicit none
	integer nn_loc
	DT,dimension(1:nn_loc)::x,y
	integer precond
	type(hobf)::hobf_inverse
	type(proctree)::ptree
	type(Hstat)::stats
	
	if(precond==NOPRECON)then
		y=x
	else if (precond==HODLRPRECON)then 
		call MVM_Z_factorized(nn_loc,1,x,y,hobf_inverse,ptree,stats)	 
	endif
	end subroutine HODLR_ApplyPrecon

	
	
subroutine HODLR_Test_Solve_error(ho_bf_inv,option,ptree,stats)
    
    use z_HODLR_DEFS
    
    implicit none
    
    integer i, j, ii, jj, iii, jjj, ierr
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, n_iter_max, iter, N_unk, N_unk_loc
    real(kind=8) theta, phi, dphi, rcs_V, rcs_H
    real T0
    real(kind=8) n1,n2,rtemp	
    DT value_Z
    DT,allocatable:: Voltage_pre(:),x(:,:),xtrue(:,:),b(:,:),btrue(:,:)
	real(kind=8):: rel_error,rtemp1,rtemp2,rtemp3,rtemp4,norm1,norm2,norm3,norm4
	type(Hoption)::option
	! type(mesh)::msh
	! type(kernelquant)::ker
	type(proctree)::ptree
	type(hobf)::ho_bf_inv
	type(Hstat)::stats	
	DT,allocatable:: current(:),voltage(:)
	integer idxs,idxe

	! if(option%PRECON==DIRECT)then
		idxs = ho_bf_inv%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		idxe = ho_bf_inv%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)
	! else 
		! ! write(*,*)associated(ho_bf_for%levels(1)%BP_inverse),'dd' !%matrices_block(1)%N_p),'nima'
		! msh%idxs = ho_bf_for%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		! msh%idxe = ho_bf_for%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)	
	! endif
	
	! N_unk=msh%Nunk
	N_unk_loc = idxe-idxs+1	
	
	allocate (x(N_unk_loc,1))
	x=0		
	allocate (xtrue(N_unk_loc,1))
	call RandomMat(N_unk_loc,1,1,xtrue,0)	
	allocate (btrue(N_unk_loc,1))	
	btrue=0	
	allocate (b(N_unk_loc,1))	
	b=0
	call MVM_Z_forward('N',N_unk_loc,1,1,ho_bf_inv%Maxlevel+1,xtrue,btrue,ho_bf_inv,ptree,stats)
	stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
	call HODLR_Solution(ho_bf_inv,x,btrue,N_unk_loc,1,option,ptree,stats)
	call MVM_Z_forward('N',N_unk_loc,1,1,ho_bf_inv%Maxlevel+1,x,b,ho_bf_inv,ptree,stats)
	stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
	
	rtemp1 = fnorm(xtrue-x,N_unk_loc,1)**2d0;
	rtemp2 = fnorm(xtrue,N_unk_loc,1)**2d0;
	
	call MPI_ALLREDUCE(rtemp1, norm1, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
	call MPI_ALLREDUCE(rtemp2, norm2, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)then
		write(*,*)'||X_t-H\(H*X_t)||_F/||X_t||_F: ',sqrt(norm1/norm2)
	endif
	
	rtemp3 = fnorm(btrue-b,N_unk_loc,1)**2d0;
	rtemp4 = fnorm(btrue,N_unk_loc,1)**2d0;
	
	call MPI_ALLREDUCE(rtemp3, norm3, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
	call MPI_ALLREDUCE(rtemp4, norm4, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
	if(ptree%MyID==Main_ID)then
		write(*,*)'||B-H*(H\B)||_F/||B||_F: ',sqrt(norm3/norm4)
	endif	
	
	
	deallocate(x)
	deallocate(xtrue)
	deallocate(btrue)
	deallocate(b)
	
end subroutine HODLR_Test_Solve_error
	

subroutine MVM_Z_factorized(Ns,num_vectors,Vin,Vout,ho_bf1,ptree,stats)

    use z_HODLR_DEFS
    
    implicit none
    
	integer Ns
	integer level_c,rowblock,head,tail
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors,pp
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real(kind=8) a,b,c,d
    DT ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o
	
    ! type(vectorsblock), pointer :: random1, random2
    
    real(kind=8),allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	
	DT,allocatable::vec_old(:,:),vec_new(:,:)
	! complex(kind=8)::Vin(:),Vout(:)
	DT::Vin(Ns,num_vectors),Vout(Ns,num_vectors)
	type(hobf)::ho_bf1
    type(proctree)::ptree
	type(Hstat)::stats
		
	idx_start_glo = ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)	 
	   		
	! get the right multiplied vectors
	ctemp1=1.0d0 ; ctemp2=0.0d0
	! allocate(vec_old(Ns,num_vectors))
	allocate(vec_new(Ns,num_vectors))
	Vout = Vin
	! write(*,*)'ddddd',Ns,num_vectors
	! write(*,*)'begin'
	
	do level = ho_bf1%Maxlevel+1,1,-1
		vec_new = 0
		do ii = ho_bf1%levels(level)%Bidxs,ho_bf1%levels(level)%Bidxe
			pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level)%BP_inverse(ii)%pgno)%head + 1
			head = ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%N_p(pp,1) + ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%headn -1
			tail = head + ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%N_loc -1
			idx_start_loc = head-idx_start_glo+1
			idx_end_loc = tail-idx_start_glo+1					
		
			if(level==ho_bf1%Maxlevel+1)then	
				call fullmat_block_MVP_dat(ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vectors,&
				&Vout(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2)
				stats%Flop_Sol = stats%Flop_Sol + flops_zgemm(idx_end_loc-idx_start_loc+1,num_vectors,idx_end_loc-idx_start_loc+1)
			else 
				stats%Flop_Tmp=0
				call Bplus_block_MVP_inverse_dat(ho_bf1,level,ii,'N',idx_end_loc-idx_start_loc+1,num_vectors,Vout(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ptree,stats)
				stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
				
			endif
		end do
		Vout = vec_new
	end do
	! Vout = vec_new(1:Ns,1)
	! deallocate(vec_old)
	deallocate(vec_new)	
	 
    return                

end subroutine MVM_Z_factorized

subroutine MVM_Z_forward(trans,Ns,num_vectors,level_start,level_end,Vin,Vout,ho_bf1,ptree,stats)

    use z_HODLR_DEFS
    
    implicit none
    
	character trans,trans_tmp
	integer Ns, level_start, level_end
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real(kind=8) a,b,c,d
    DT ctemp, ctemp1, ctemp2
	! type(matrixblock),pointer::block_o
	type(blockplus),pointer::bplus_o
	type(proctree)::ptree
    ! type(vectorsblock), pointer :: random1, random2
    type(Hstat)::stats
	
    real(kind=8),allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_m,idx_end_m,idx_start_n,idx_end_n,pp,head,tail,idx_start_loc,idx_end_loc
	
	DT,allocatable::vec_old(:,:),vec_new(:,:)
	! complex(kind=8)::Vin(:,:),Vout(:,:)
	DT::Vin(Ns,num_vectors),Vout(Ns,num_vectors)
	type(hobf)::ho_bf1
 
	idx_start_glo = ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1) 
	
	trans_tmp = trans
	if(trans=='C')then
		trans_tmp = 'T'
		Vin=conjg(cmplx(Vin,kind=8))
	endif
	
	! get the right multiplied vectors
	ctemp1=1.0d0 ; ctemp2=1.0d0
	! allocate(vec_old(Ns,num_vectors))
	allocate(vec_new(Ns,num_vectors))
	! vec_old(1:Ns,1:num_vectors) = Vin
	vec_new = 0
	stats%Flop_Tmp =0
	
	do level = level_start,level_end !ho_bf1%Maxlevel+1
		do ii =ho_bf1%levels(level)%Bidxs,ho_bf1%levels(level)%Bidxe
			
			pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level)%BP_inverse(ii)%pgno)%head + 1
			head = ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%N_p(pp,1) + ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%headn -1
			tail = head + ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%N_loc -1
			idx_start_loc = head-idx_start_glo+1
			idx_end_loc = tail-idx_start_glo+1					
			
			if(level==ho_bf1%Maxlevel+1)then	
				call fullmat_block_MVP_dat(ho_bf1%levels(level)%BP(ii)%LL(1)%matrices_block(1),trans_tmp,idx_end_loc-idx_start_loc+1,num_vectors,&
				&Vin(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2)
				stats%Flop_Tmp = stats%Flop_Tmp + flops_zgemm(idx_end_loc-idx_start_loc+1,num_vectors,idx_end_loc-idx_start_loc+1)				
			else
				call Bplus_block_MVP_twoforward_dat(ho_bf1,level,ii,trans_tmp,idx_end_loc-idx_start_loc+1,num_vectors,Vin(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2,ptree,stats)
			endif
		
		end do				
	end do
	
	Vout = vec_new(1:Ns,1:num_vectors)
	! deallocate(vec_old)
	deallocate(vec_new)	

	if(trans=='C')then
		Vout=conjg(cmplx(Vout,kind=8))
		Vin=conjg(cmplx(Vin,kind=8))
	endif
	
    return                

end subroutine MVM_Z_forward

end module z_HODLR_Solve_Mul