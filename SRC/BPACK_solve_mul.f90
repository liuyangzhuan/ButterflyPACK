! “ButterflyPACK” Copyright (c) 2018, The Regents of the University of California, through
! Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.

! If you have questions about your rights to use or distribute this software, please contact
! Berkeley Lab's Intellectual Property Office at  IPO@lbl.gov.

! NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the
! U.S. Government consequently retains certain rights. As such, the U.S. Government has been
! granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable
! worldwide license in the Software to reproduce, distribute copies to the public, prepare
! derivative works, and perform publicly and display publicly, and to permit other to do so.

! Developers: Yang Liu
!             (Lawrence Berkeley National Lab, Computational Research Division).

#include "ButterflyPACK_config.fi"
module BPACK_Solve_Mul

use Bplus_compress


contains


subroutine BPACK_Solution(bmat,x,b,Ns_loc,num_vectors,option,ptree,stats)

    use BPACK_DEFS
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
	type(Bmatrix)::bmat
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
			call BPACK_Ztfqmr(option%precon,option%n_iter,Ns_loc,b(:,ii),x(:,ii),rel_error,iter,r0_initial,bmat,ptree,option,stats)
		end do

		deallocate(r0_initial)
	else
		call BPACK_Inv_Mult('N',Ns_loc,num_vectors,b,x,bmat,ptree,option,stats)
	end if


	n2 = OMP_get_wtime()
	stats%Time_Sol = stats%Time_Sol + n2-n1

    return

end subroutine BPACK_Solution


  subroutine BPACK_Ztfqmr(precond,ntotal,nn_loc,b,x,err,iter,r0_initial,bmat,ptree,option,stats)
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
	type(Bmatrix)::bmat
	DT::r0_initial(:)
	integer precond
	type(proctree)::ptree
	type(Hoption)::option

    itmax=iter

	call BPACK_ApplyPrecon(precond,nn_loc,b,bb,ptree,bmat,option,stats)


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
    call BPACK_Mult('N',nn_loc,1,x,ytmp,bmat,ptree,option,stats)
	stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
	call BPACK_ApplyPrecon(precond,nn_loc,ytmp,r,ptree,bmat,option,stats)

    r=bb-r !residual from the initial guess
    w=r
    yo=r
        ! ! write(*,*)'2'
    	! ! if(isnan(sum(abs(yo)**2)))then
			! ! write(*,*)'shitddd'
			! ! stop
		! ! end if
	! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,yo,ayo)
    call BPACK_Mult('N',nn_loc,1,yo,ytmp,bmat,ptree,option,stats)
	stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
	call BPACK_ApplyPrecon(precond,nn_loc,ytmp,ayo,ptree,bmat,option,stats)

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
	   call BPACK_Mult('N',nn_loc,1,ye,ytmp,bmat,ptree,option,stats)
	   stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
	   call BPACK_ApplyPrecon(precond,nn_loc,ytmp,aye,ptree,bmat,option,stats)

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
		  call BPACK_Mult('N',nn_loc,1,x,ytmp,bmat,ptree,option,stats)
		  stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
		  call BPACK_ApplyPrecon(precond,nn_loc,ytmp,r,ptree,bmat,option,stats)

          r=bb-r
          rerr_local=dot_product(r,r)
          call MPI_ALLREDUCE(rerr_local,rerr_sum,1,MPI_DT,MPI_SUM,&
               ptree%Comm,ierr)
          rerr=sqrt(abs(rerr_sum))/bmag

          if (ptree%MyID==Main_ID .and. option%verbosity>=0) then
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
		call BPACK_Mult('N',nn_loc,1,yo,ytmp,bmat,ptree,option,stats)
		stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
		call BPACK_ApplyPrecon(precond,nn_loc,ytmp,ayo,ptree,bmat,option,stats)

       !MAGIC
       v=ayo+beta*( aye+beta*v )
    enddo iters
    ! write(*,*)'6'
    ! call SmartMultifly(trans,nn_loc,level_c,rowblock,1,x,r)
    call BPACK_Mult('N',nn_loc,1,x,ytmp,bmat,ptree,option,stats)
	stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
	call BPACK_ApplyPrecon(precond,nn_loc,ytmp,r,ptree,bmat,option,stats)


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
  end subroutine BPACK_Ztfqmr



  subroutine BPACK_ApplyPrecon(precond,nn_loc,x,y,ptree,bmat,option,stats)
    implicit none
	integer nn_loc
	DT,dimension(1:nn_loc)::x,y
	integer precond
	type(Bmatrix)::bmat
	type(proctree)::ptree
	type(Hstat)::stats
	type(Hoption)::option

	if(precond==NOPRECON)then
		y=x
	else if (precond==HODLRPRECON)then
		call BPACK_Inv_Mult('N',nn_loc,1,x,y,bmat,ptree,option,stats)
	endif
	end subroutine BPACK_ApplyPrecon



subroutine BPACK_Test_Solve_error(bmat,N_unk_loc,option,ptree,stats)

    use BPACK_DEFS

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
	type(Bmatrix)::bmat
	type(Hstat)::stats
	DT,allocatable:: current(:),voltage(:)
	integer idxs,idxe


	! idxs = bmat%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
	! idxe = bmat%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)


	! ! N_unk=msh%Nunk
	! N_unk_loc = idxe-idxs+1

	allocate (x(N_unk_loc,1))
	x=0
	allocate (xtrue(N_unk_loc,1))
	xtrue=0
	call RandomMat(N_unk_loc,1,1,xtrue,0)
	allocate (btrue(N_unk_loc,1))
	btrue=0
	allocate (b(N_unk_loc,1))
	b=0
	call BPACK_Mult('N',N_unk_loc,1,xtrue,btrue,bmat,ptree,option,stats)
	stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp
	call BPACK_Solution(bmat,x,btrue,N_unk_loc,1,option,ptree,stats)
	call BPACK_Mult('N',N_unk_loc,1,x,b,bmat,ptree,option,stats)
	stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp

	rtemp1 = fnorm(xtrue-x,N_unk_loc,1)**2d0;
	rtemp2 = fnorm(xtrue,N_unk_loc,1)**2d0;

	call MPI_ALLREDUCE(rtemp1, norm1, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
	call MPI_ALLREDUCE(rtemp2, norm2, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)then
		write(*,*)'||X_t-H\(H*X_t)||_F/||X_t||_F: ',sqrt(norm1/norm2)
	endif

	rtemp3 = fnorm(btrue-b,N_unk_loc,1)**2d0;
	rtemp4 = fnorm(btrue,N_unk_loc,1)**2d0;

	call MPI_ALLREDUCE(rtemp3, norm3, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
	call MPI_ALLREDUCE(rtemp4, norm4, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)then
		write(*,*)'||B-H*(H\B)||_F/||B||_F: ',sqrt(norm3/norm4)
	endif


	deallocate(x)
	deallocate(xtrue)
	deallocate(btrue)
	deallocate(b)

end subroutine BPACK_Test_Solve_error


subroutine BPACK_Inv_Mult(trans,Ns,num_vectors,Vin,Vout,bmat,ptree,option,stats)
    implicit none

	integer Ns
	character trans
	integer num_vectors
	DT::Vin(Ns,num_vectors),Vout(Ns,num_vectors)
	type(Bmatrix)::bmat
    type(proctree)::ptree
	type(Hstat)::stats
	type(Hoption)::option

	select case(option%format)
    case(HODLR)
		call HODLR_Inv_Mult(trans,Ns,num_vectors,Vin,Vout,bmat%ho_bf,ptree,option,stats)
    case(HMAT)
		call Hmat_Inv_Mult(trans,Ns,num_vectors,Vin,Vout,bmat%h_mat,ptree,option,stats)
	end select

end subroutine BPACK_Inv_Mult


subroutine BPACK_Mult(trans,Ns,num_vectors,Vin,Vout,bmat,ptree,option,stats)
    use BPACK_DEFS
    implicit none
	character trans
	integer Ns
	integer num_vectors
	DT::Vin(Ns,num_vectors),Vout(Ns,num_vectors)
	type(Bmatrix)::bmat
	type(proctree)::ptree
    type(Hstat)::stats
    type(Hoption)::option


	select case(option%format)
    case(HODLR)
		call HODLR_Mult(trans,Ns,num_vectors,1,bmat%ho_bf%Maxlevel+1,Vin,Vout,bmat%ho_bf,ptree,option,stats)
    case(HMAT)
		call Hmat_Mult(trans,Ns,num_vectors,Vin,Vout,bmat%h_mat,ptree,option,stats)
	end select


end subroutine BPACK_Mult




subroutine HODLR_Inv_Mult(trans,Ns,num_vectors,Vin,Vout,ho_bf1,ptree,option,stats)

    use BPACK_DEFS

    implicit none

	integer Ns
	integer level_c,rowblock,head,tail
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors,pp
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character trans,trans_tmp
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
	type(Hoption)::option

	idx_start_glo = ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)

	stats%Flop_Tmp=0

	trans_tmp = trans
	if(trans=='C')then
		trans_tmp = 'T'
		Vin=conjg(cmplx(Vin,kind=8))
	endif

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
				call Full_block_MVP_dat(ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1),trans_tmp,idx_end_loc-idx_start_loc+1,num_vectors,&
				&Vout(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2)
				stats%Flop_Sol = stats%Flop_Sol + flops_zgemm(idx_end_loc-idx_start_loc+1,num_vectors,idx_end_loc-idx_start_loc+1)
			else
				call Bplus_block_MVP_inverse_dat(ho_bf1,level,ii,trans_tmp,idx_end_loc-idx_start_loc+1,num_vectors,Vout(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ptree,stats)

			endif
		end do
		Vout = vec_new
	end do
	! Vout = vec_new(1:Ns,1)
	! deallocate(vec_old)
	deallocate(vec_new)

    stats%Flop_Sol = stats%Flop_Sol + stats%Flop_Tmp

		! do ii=1,Ns
		! write(131,*)abs(Vout(ii,1))
		! enddo

	if(trans=='C')then
		Vout=conjg(cmplx(Vout,kind=8))
		Vin=conjg(cmplx(Vin,kind=8))
	endif
	Vout = Vout*option%scale_factor

    return

end subroutine HODLR_Inv_Mult






subroutine HODLR_Mult(trans,Ns,num_vectors,level_start,level_end,Vin,Vout,ho_bf1,ptree,option,stats)
    use BPACK_DEFS
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
    type(Hoption)::option

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
				call Full_block_MVP_dat(ho_bf1%levels(level)%BP(ii)%LL(1)%matrices_block(1),trans_tmp,idx_end_loc-idx_start_loc+1,num_vectors,&
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

	Vout = Vout/option%scale_factor

	! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)"output norm: ",fnorm(Vout,Ns,num_vectors)**2d0


    return

end subroutine HODLR_Mult


subroutine Hmat_Inv_Mult(trans,Ns,num_vectors,Vin,Vout,h_mat,ptree,option,stats)
    implicit none

	integer Ns,ii
	character trans,trans_tmp
	integer num_vectors
	DT::Vin(Ns,num_vectors),Vout(Ns,num_vectors)
	type(Hmat)::h_mat
    type(proctree)::ptree
	type(Hstat)::stats
	type(Hoption)::option

	trans_tmp = trans
	if(trans=='C')then
		trans_tmp = 'T'
		Vin=conjg(cmplx(Vin,kind=8))
	endif

	Vout = Vin

	if(trans=='N')then
		! write(*,*)fnorm(Vout,Ns,num_vectors),'before L solve'
		call Hmat_Lsolve_Toplevel(h_mat,trans_tmp,Vout,Ns,num_vectors,ptree,stats)
		! write(*,*)fnorm(Vout,Ns,num_vectors),'before U solve',abs(Vout)
		call Hmat_Usolve_Toplevel(h_mat,trans_tmp,Vout,Ns,num_vectors,ptree,stats)
		! write(*,*)fnorm(Vout,Ns,num_vectors),'after LU solve'
		! do ii=1,Ns
		! write(130,*)abs(Vout(ii,1))
		! enddo
	else
		call Hmat_Usolve_Toplevel(h_mat,trans_tmp,Vout,Ns,num_vectors,ptree,stats)
		call Hmat_Lsolve_Toplevel(h_mat,trans_tmp,Vout,Ns,num_vectors,ptree,stats)
	endif

	if(trans=='C')then
		Vout=conjg(cmplx(Vout,kind=8))
		Vin=conjg(cmplx(Vin,kind=8))
	endif

	Vout = Vout*option%scale_factor

end subroutine Hmat_Inv_Mult


subroutine Hmat_Mult(trans,Ns,num_vectors,Vin,Vout,h_mat,ptree,option,stats)

    use BPACK_DEFS

    implicit none

	character trans,trans_tmp
	integer Ns
	integer level_c,rowblock
    integer j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real(kind=8) vecnorm,n1,n2
    DT ctemp, ctemp1, ctemp2
	! type(matrixblock),pointer::block_o
	type(blockplus),pointer::bplus_o
	type(proctree)::ptree
    ! type(vectorsblock), pointer :: random1, random2
    type(Hstat)::stats
    type(Hoption)::option

    real(kind=8),allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_m,idx_end_m,idx_start_n,idx_end_n,pp,head,tail,idx_start_loc,idx_end_loc,Nmax,Nmsg
	type(matrixblock), pointer :: blocks_i,blocks_j
	DT,allocatable::vec_old(:,:),vec_new(:,:),vin_tmp(:,:),vout_tmp(:,:),vec_buffer(:,:)

	! DT::Vin(:,:),Vout(:,:)
	DT::Vin(Ns,num_vectors),Vout(Ns,num_vectors)
	type(Hmat)::h_mat
	integer ierr
	integer,allocatable::status_all(:,:),srequest_all(:)
	integer :: status(MPI_Status_size)

	if(ptree%MyID/=MPI_COMM_NULL)then

		trans_tmp = trans
		if(trans=='C')then
			trans_tmp = 'T'
			Vin=conjg(cmplx(Vin,kind=8))
		endif

		call MPI_barrier(ptree%Comm,ierr)


		n1=OMP_get_wtime()


		blocks_i=>h_mat%Local_blocks_copy(1,1)
		num_blocks=2**h_mat%Dist_level

		Nmsg=num_blocks

		call MPI_ALLREDUCE(blocks_i%N, Nmax, 1,MPI_INTEGER, MPI_MAX, ptree%Comm,ierr)
		allocate(srequest_all(Nmsg))
		allocate(status_all(MPI_status_size,Nmsg))
		srequest_all=MPI_request_null
		do j=1,Nmsg
			call MPI_Isend(Vin,blocks_i%N*num_vectors,MPI_DT,j-1,ptree%MyID,ptree%Comm,srequest_all(j),ierr)
		enddo
		allocate(vec_buffer(Nmax,num_vectors))

		vout=0.0
		do while(Nmsg>0)
			call MPI_Recv(vec_buffer,Nmax*num_vectors,MPI_DT,MPI_ANY_SOURCE,MPI_ANY_TAG,ptree%Comm,status,ierr)
			j = status(MPI_TAG)+1
			blocks_i=>h_mat%Local_blocks_copy(j,1)
			allocate(vin_tmp(blocks_i%N,num_vectors))
			vin_tmp = vec_buffer(1:blocks_i%N,1:num_vectors)
			call Hmat_block_MVP_dat(blocks_i,trans_tmp,blocks_i%headm,blocks_i%headn,num_vectors,vin_tmp,Vout,cone,ptree,stats)
			deallocate(vin_tmp)

			Nmsg = Nmsg -1
		enddo

		deallocate(vec_buffer)
		call MPI_waitall(num_blocks,srequest_all,status_all,ierr)
		deallocate(srequest_all)
		deallocate(status_all)


		if(trans=='C')then
			Vout=conjg(cmplx(Vout,kind=8))
			Vin=conjg(cmplx(Vin,kind=8))
		endif

		Vout = Vout/option%scale_factor

		call MPI_barrier(ptree%Comm,ierr)
		n2=OMP_get_wtime()

		vecnorm = fnorm(Vout,Ns,num_vectors)**2d0
		call MPI_AllREDUCE(MPI_IN_PLACE, vecnorm, 1,MPI_double, MPI_SUM, ptree%Comm,ierr)
		! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*)"output norm: ",sqrt(vecnorm)
	endif

    return

end subroutine Hmat_Mult


subroutine Hmat_Lsolve_Toplevel(h_mat,trans,xloc,nloc,nvec,ptree,stats)
    implicit none

	type(Hmat)::h_mat
	character::trans
	type(proctree)::ptree
	type(Hstat)::stats

    integer i, j, k, ii, jj, kk, iii, jjj, num_blocks, mm, nn,groupm,groupn
    integer vectors_start, vectors_x, vectors_y, id_l,nvec,nloc

    type(matrixblock), pointer :: blocks_l
    integer Nreq, Nmod, Bufsize
	DT,allocatable::recv_buf(:),vin(:,:)
	DT::xloc(:,:)
    integer :: status(MPI_Status_size)
	integer,allocatable::status_all(:,:),request_all(:)
	integer idx_start
	integer ierr

	if(trans=='N')then
		num_blocks=2**h_mat%Dist_level
		vectors_start=num_blocks-1

		blocks_l=>h_mat%Local_blocks(ptree%MyID+1,1)
		groupm = blocks_l%row_group

		vectors_x=vectors_start+ptree%MyID+1

		Bufsize = nloc*nvec
		call MPI_ALLREDUCE(MPI_IN_PLACE, Bufsize, 1,MPI_integer, MPI_MAX, ptree%Comm,ierr)
		allocate(recv_buf(Bufsize))

		Nmod = ptree%MyID
		do while(Nmod>0)
			call MPI_Recv(recv_buf,Bufsize,MPI_DT,MPI_ANY_SOURCE,MPI_ANY_TAG,ptree%Comm,status,ierr)
			vectors_y = status(MPI_TAG)
			j=vectors_y-vectors_start
			blocks_l=>h_mat%Local_blocks(j,1)
			groupn = blocks_l%col_group
			nn=blocks_l%N
			allocate(vin(nn,nvec))
			j=0
			do jjj=1,nvec
			do iii=1,nn
				j=j+1
				vin(iii,jjj)=recv_buf(j)
			enddo
			enddo
			call Hmat_block_MVP_dat(blocks_l,'N',blocks_l%headm,blocks_l%headn,nvec,vin,xloc,-cone,ptree,stats)
			deallocate(vin)
			Nmod = Nmod-1
		enddo
		deallocate(recv_buf)


		blocks_l=>h_mat%Local_blocks(ptree%MyID+1,1)
		idx_start = blocks_l%headm
		call Hmat_Lsolve(blocks_l,'N',idx_start,nvec,xloc,ptree,stats)

		Nreq=num_blocks-1-ptree%MyID
		if(Nreq>0)then
			allocate(status_all(MPI_status_size,Nreq))
			allocate(request_all(Nreq))
		end if
		do i=1, Nreq
			call MPI_Isend(xloc,nloc*nvec,MPI_DT,i+ptree%MyID,vectors_x,ptree%Comm,request_all(i),ierr)
		enddo
		if(Nreq>0)then
			call MPI_waitall(Nreq,request_all,status_all,ierr)
			deallocate(status_all)
			deallocate(request_all)
		endif
	else
		write(*,*)'XxL^-1 with MPI is not yet implemented'
		stop
	endif

	call MPI_barrier(ptree%Comm,ierr)
    return

end subroutine Hmat_Lsolve_Toplevel





subroutine Hmat_Usolve_Toplevel(h_mat,trans,xloc,nloc,nvec,ptree,stats)

    implicit none

	type(Hmat)::h_mat
	character::trans
	type(proctree)::ptree
	type(Hstat)::stats

    integer i, j, k, ii, jj, kk, iii, jjj, num_blocks, mm, nn,idx_start,nvec,nloc
    integer vectors_start, vectors_x, vectors_y, id_u,groupm,groupn,ierr

    type(matrixblock), pointer :: blocks_u
	integer :: status(MPI_Status_size)
    integer Nreq, Nmod, Bufsize
	DT,allocatable::recv_buf(:),vin(:,:)
	DT::xloc(:,:)
	integer,allocatable::status_all(:,:),request_all(:)

	if(trans=='N')then

		num_blocks=2**h_mat%Dist_level
		vectors_start=num_blocks-1

		Bufsize = nloc*nvec
		call MPI_ALLREDUCE(MPI_IN_PLACE, Bufsize, 1,MPI_integer, MPI_MAX, ptree%Comm,ierr)
		allocate(recv_buf(Bufsize))


		vectors_x=vectors_start+ptree%MyID+1
		Nmod = num_blocks -1 - ptree%MyID

		blocks_u=>h_mat%Local_blocks(ptree%MyID+1,1)

		do while(Nmod>0)
			call MPI_Recv(recv_buf,Bufsize,MPI_DT,MPI_ANY_SOURCE,MPI_ANY_TAG,ptree%Comm,status,ierr)
			vectors_y = status(MPI_TAG)
			j=vectors_y-vectors_start
			blocks_u=>h_mat%Local_blocks(j,1)

			nn=blocks_u%N
			allocate(vin(nn,nvec))
			j=0
			do jjj=1,nvec
			do iii=1,nn
				j=j+1
				vin(iii,jjj)=recv_buf(j)
			enddo
			enddo

			call Hmat_block_MVP_dat(blocks_u,'N',blocks_u%headm,blocks_u%headn,nvec,vin,xloc,-cone,ptree,stats)

			deallocate(vin)

			Nmod = Nmod-1
		enddo
		deallocate(recv_buf)

		blocks_u=>h_mat%Local_blocks(ptree%MyID+1,1)
		idx_start = blocks_u%headm
		call Hmat_Usolve(blocks_u,'N',idx_start,nvec,xloc,ptree,stats)

		Nreq=ptree%MyID
		if(Nreq>0)then
			allocate(status_all(MPI_status_size,Nreq))
			allocate(request_all(Nreq))
		end if

		do i=Nreq,1,-1
			call MPI_Isend(xloc,nloc*nvec,MPI_DT,i-1,vectors_x,ptree%Comm,request_all(i),ierr)
		enddo

		if(Nreq>0)then
			call MPI_waitall(Nreq,request_all,status_all,ierr)
			deallocate(status_all)
			deallocate(request_all)
		endif

	else
		write(*,*)'XxU^-1 with MPI is not yet implemented'
		stop
	endif

	call MPI_barrier(ptree%Comm,ierr)
    return
end subroutine Hmat_Usolve_Toplevel









end module BPACK_Solve_Mul