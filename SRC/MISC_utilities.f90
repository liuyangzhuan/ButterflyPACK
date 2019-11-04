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
#ifdef Intel
#ifndef USEVSL
#define USEVSL
include "mkl_vsl.f90"
#endif
#endif


module MISC_Utilities
use BPACK_DEFS
#ifdef Intel
USE IFPORT
#endif
use omp_lib
use MISC_DenseLA
use BPACK_linkedlist


integer, parameter :: int64 = selected_int_kind(18)

contains

! #ifndef mymacro(x)
! #define mymacro(x) print *, "Now giving information about ", "x" ; \
                   ! call mysub( x, size(x,1), size(x,2) ) ; \
                   ! print *, "About to do function on ", "x"; \
                   ! call dofunction(x) ; \
                   ! print *, "x"," is a nice matrix! Huzzah!"
! #endif


subroutine linspaceI(startI,endI,N,array)
implicit none
integer startI,endI,N
integer array(1:N)
real(kind=8)::rtemp
integer i

rtemp=dble(endI-startI+1)/dble(N)

array(1)=int(dble(endI-startI+1)/dble(N)/2.) + startI - 1
if (array(1)==startI - 1) then
	array(1)=startI
endif

do i=2, N
	array(i)=array(1)+int(dble(i-1)*rtemp)
enddo

end subroutine linspaceI


  subroutine copysubmat(A,ais,ajs,B,bis,bjs,m,n,trans)
	implicit none
	integer ais,ajs,bis,bjs,m,n
	DT::A(:,:),B(:,:)
	real(kind=8)::n1,n2
	integer ii,jj,ijind
	character trans

	! n1 = OMP_get_wtime()
	if(trans=='N')then
		do ii =1,m
		do jj =1,n
			B(ii+bis-1,jj+bjs-1) = A(ii+ais-1,jj+ajs-1)
		end do
		end do
	elseif(trans=='T')then
		do ii =1,m
		do jj =1,n
			B(jj+bis-1,ii+bjs-1) = A(ii+ais-1,jj+ajs-1)
		end do
		end do
	endif
	! n2 = OMP_get_wtime()
	! time_memcopy = time_memcopy + n2-n1

 end subroutine copysubmat




  subroutine copymatT(A,B,m,n)
	implicit none
	integer m,n
	DT::A(:,:),B(:,:)
	real(kind=8)::n1,n2
	integer ii,jj,ijind

	! n1 = OMP_get_wtime()


	! !$omp parallel do default(shared) private(ii,jj)
		do ii =1,m
		do jj =1,n
			B(jj,ii) = A(ii,jj)
		end do
		end do
	! !$omp end parallel do

	! n2 = OMP_get_wtime()
	! time_memcopy = time_memcopy + n2-n1

 end subroutine copymatT


function isnanMat(A,m,n)
	implicit none
	logical:: isnanMat
	DT::A(:,:)
	integer m,n,ii,jj
	isnanMat = .false.
	do ii =1,m
	do jj =1,n
		isnanMat = isnanMat .or. isnan(abs(A(ii,jj)))
	end do
	end do
 end function isnanMat


! function fnorm(A,m,n)
	! implicit none
	! real(kind=8):: fnorm
	! DT::A(m,n)
	! integer m,n,ii,jj
	! fnorm = 0
	! do ii =1,m
	! do jj =1,n
		! fnorm = fnorm + abs(A(ii,jj))**2d0
	! end do
	! end do
	! fnorm = sqrt(fnorm)
 ! end function fnorm


subroutine LR_ReCompression(matU,matV,M,N,rmax,rank,SVD_tolerance,Flops)


    use BPACK_DEFS
    implicit none

	integer N,M,rmax
	real(kind=8),optional:: Flops
	real(kind=8):: flop
	real(kind=8) SVD_tolerance
	integer i, j, ii, jj, indx, rank_1, rank_2
	integer rank, ranknew,ldaU,ldaV

	DT::matU(:,:),matV(:,:)
	DT, allocatable :: QQ1(:,:), RR1(:,:),QQ2(:,:), RR2(:,:), UUsml(:,:), VVsml(:,:),tau_Q(:),mattemp(:,:),matU1(:,:),matV1(:,:)
	real(kind=8), allocatable :: Singularsml(:)
	integer,allocatable :: jpvt1(:),jpvt2(:)

	if(present(Flops))Flops=0d0

	ldaU = size(matU,1)
	ldaV = size(matV,1)

	call assert(rmax<=min(M,N),'rmax too big in LR_ReCompression')

	allocate(QQ1(M,rmax))
	! call copymatN(matU(1:M,1:rmax),QQ1,M,rmax)
	QQ1 = matU(1:M,1:rmax)
	allocate (tau_Q(rmax))
	! allocate (jpvt1(rmax))
	! jpvt1=0
	! call geqp3f90(QQ1,jpvt1,tau_Q)
	call geqrff90(QQ1,tau_Q,flop=flop)
	if(present(Flops))Flops = Flops + flop

	allocate (RR1(rmax,rmax))
	RR1=0d0
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rmax
		do i=1, j
			RR1(i,j)=QQ1(i,j)
		enddo
	enddo
	! !$omp end parallel do
	call un_or_gqrf90(QQ1,tau_Q,M,rmax,rmax,flop=flop)
	deallocate(tau_Q)
	! deallocate(jpvt1)
	if(present(Flops))Flops = Flops + flop

	allocate(QQ2(N,rmax))
	call copymatT(matV(1:rmax,1:N),QQ2,rmax,N)
	allocate (tau_Q(rmax))
	! call geqrff90(QQ2,tau_Q)
	! allocate (jpvt2(rmax))
	! jpvt2=0
	call geqrff90(QQ2,tau_Q,flop=flop)

	if(present(Flops))Flops = Flops + flop

	allocate (RR2(rmax,rmax))
	RR2=0d0
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rmax
		do i=1, j
			RR2(i,j)=QQ2(i,j)
		enddo
	enddo
	! !$omp end parallel do

	call un_or_gqrf90(QQ2,tau_Q,N,rmax,rmax,flop=flop)
	! deallocate(jpvt2)
	deallocate(tau_Q)
	if(present(Flops))Flops = Flops + flop

	allocate(mattemp(rmax,rmax))
	mattemp=0
	! call zgemm('N','T',rmax,rmax,rmax, cone, RR1, rmax,RR2,rmax,czero,mattemp,rmax)
	call gemmf90(RR1,rmax,RR2,rmax,mattemp,rmax,'N','T',rmax,rmax,rmax, cone,czero,flop=flop)
	if(present(Flops))Flops = Flops + flop
	allocate(UUsml(rmax,rmax),VVsml(rmax,rmax),Singularsml(rmax))
	call SVD_Truncate(mattemp,rmax,rmax,rmax,UUsml,VVsml,Singularsml,SVD_tolerance,ranknew,flop=flop)
	if(present(Flops))Flops = Flops + flop
	do i=1,ranknew
		UUsml(1:rmax,i) = UUsml(1:rmax,i) * Singularsml(i)
	enddo
	! call zgemm('N','N',M,ranknew,rmax, cone, QQ1, M,UUsml,rmax,czero,matU,ldaU)
	call gemmf90(QQ1,M,UUsml,rmax,matU,ldaU,'N','N',M,ranknew,rmax, cone,czero,flop=flop)
	if(present(Flops))Flops = Flops + flop
	! call zgemm('N','T',ranknew,N,rmax, cone, VVsml, rmax,QQ2,N,czero,matV,ldaV)
	call gemmf90(VVsml,rmax,QQ2,N,matV,ldaV,'N','T',ranknew,N,rmax, cone,czero,flop=flop)
	if(present(Flops))Flops = Flops + flop

	rank = ranknew


	deallocate(mattemp,RR1,QQ1,UUsml,VVsml,Singularsml)
	deallocate(QQ2,RR2)


end subroutine LR_ReCompression


subroutine LR_FnormUp(matU,matV,M,N,ruv,rup,ldV,normUV,normUVupdate,tolerance,Flops)


    use BPACK_DEFS
    implicit none

	integer N,M,ruv,rup,ldV
	real(kind=8) normUVupdate,normUV,inner_UV
	real(kind=8),optional:: Flops
	real(kind=8):: flop
	real(kind=8) tolerance
	integer i, j, k, ii, jj, indx, rank_1, rank_2
	integer rank, ranknew,mn, ranknew1, ranknew2
	DT::ctemp,inner_V,inner_U
	DT::matU(:,:),matV(:,:)
	DT, allocatable :: QQ1(:,:), RR1(:,:),QQ2(:,:), RR2(:,:), UUsml(:,:), VVsml(:,:),tau_Q(:),mattemp(:,:),matU1(:,:),matV1(:,:),UU(:,:),VV(:,:)
	real(kind=8), allocatable :: Singularsml(:),Singular(:)
	integer,allocatable :: jpvt(:),jpvt1(:),jpvt2(:)

	if(present(Flops))Flops=0

		! inner_UV=0
		! !$omp parallel do default(shared) private(j,i,k,ctemp,inner_V,inner_U) reduction(+:inner_UV)
		! do j=1,ruv
			! do k=1,rup
			! inner_U = dot_product(matU(:,j),matU(:,ruv+k))
			! inner_V = dot_product(matV(j,:),matV(ruv+k,:))
			! inner_UV=inner_UV+2*dble(inner_U*inner_V)
			! enddo
		! enddo
		! !$omp end parallel do


		inner_UV=0
		if(ruv>0)then
		allocate(UU(ruv,rup))
		UU=0
		call gemmf77('T','N',ruv,rup,M,cone,matU,M,matU(1,ruv+1),M,czero,UU,ruv)
		allocate(VV(ruv,rup))
		VV=0
		call gemmf77('N','T',ruv,rup,N,cone,matV,ldV,matV(ruv+1,1),ldV,czero,VV,ruv)
		do j=1,ruv
			do k=1,rup
				inner_UV = inner_UV + 2*dble(UU(j,k)*VV(j,k))
			enddo
		enddo
		deallocate(UU)
		deallocate(VV)
		endif


		normUV = sqrt(normUV**2d0+normUVupdate**2d0+inner_UV)

		if(present(Flops))Flops= Flops + (M+N)*ruv*rup

end subroutine LR_FnormUp

subroutine LR_Fnorm(matU,matV,M,N,rmax,norm,tolerance,Flops)


    use BPACK_DEFS
    implicit none

	integer N,M,rmax
	real(kind=8) norm
	real(kind=8),optional:: Flops
	real(kind=8):: flop
	real(kind=8) tolerance
	integer i, j, ii, jj, indx, rank_1, rank_2
	integer rank, ranknew,mn, ranknew1, ranknew2

	DT::matU(:,:),matV(:,:)
	DT, allocatable :: QQ1(:,:), RR1(:,:),QQ2(:,:), RR2(:,:), UUsml(:,:), VVsml(:,:),tau_Q(:),mattemp(:,:),matU1(:,:),matV1(:,:),UU(:,:),VV(:,:)
	real(kind=8), allocatable :: Singularsml(:),Singular(:)
	integer,allocatable :: jpvt(:),jpvt1(:),jpvt2(:)

	if(present(Flops))Flops=0

	! mn= min(M,rmax)


	! allocate(QQ1(M,rmax))
	! call copymatN(matU(1:M,1:rmax),QQ1,M,rmax)

	! allocate(UU(M,mn))
	! allocate(VV(mn,rmax))
	! allocate(Singular(mn))
	! call gesvd_robust(QQ1,Singular,UU,VV,M,rmax,mn)
	! do ii=1,mn
		! VV(ii,:) = VV(ii,:)*Singular(ii)
	! enddo
	! ! allocate(QQ2(rmax,N))
	! ! call copymatN(matV(1:rmax,1:N),QQ2,rmax,N)

	! allocate(mattemp(mn,N))
	! mattemp=0
	! call zgemm('N','N',mn,N,rmax, cone, VV, mn,matV,size(matV,1),czero,mattemp,mn)

	! norm = fnorm(mattemp,mn,N)

	! deallocate(QQ1)
	! deallocate(UU)
	! deallocate(VV)
	! deallocate(Singular)
	! deallocate(mattemp)




	allocate(QQ1(M,rmax))
	! call copymatN(matU(1:M,1:rmax),QQ1,M,rmax)
	QQ1 = matU(1:M,1:rmax)
	allocate (tau_Q(rmax))
	call geqrff90(QQ1,tau_Q,flop=flop)
	if(present(Flops))Flops= Flops + flop
	deallocate(tau_Q)

	allocate (RR1(rmax,rmax))
	RR1=0d0
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rmax
		do i=1, j
			RR1(i,j)=QQ1(i,j)
		enddo
	enddo
	! !$omp end parallel do

	allocate(QQ2(N,rmax))
	call copymatT(matV(1:rmax,1:N),QQ2,rmax,N)
	allocate (tau_Q(rmax))
	call geqrff90(QQ2,tau_Q,flop=flop)
	if(present(Flops))Flops= Flops + flop
	deallocate(tau_Q)

	allocate (RR2(rmax,rmax))
	RR2=0d0
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rmax
		do i=1, j
			RR2(i,j)=QQ2(i,j)
		enddo
	enddo
	! !$omp end parallel do

	allocate(mattemp(rmax,rmax))
	mattemp=0
	! call zgemm('N','T',rmax,rmax,rmax, cone, RR1, rmax,,rmax,czero,mattemp,rmax)
	call gemmf90(RR1,rmax,RR2,rmax,mattemp,rmax,'N','T',rmax,rmax,rmax, cone,czero,flop=flop)
	if(present(Flops))Flops= Flops + flop
	norm = fnorm(mattemp,rmax,rmax)

	deallocate(mattemp,RR1,QQ1)
	deallocate(QQ2,RR2)




	! allocate(QQ1(M,rmax))
	! QQ1 = matU(1:M,1:rmax)
	! allocate (tau_Q(rmax))
	! allocate (jpvt1(rmax))
	! jpvt1=0
	! call geqp3modf90(QQ1,jpvt1,tau_Q,tolerance,SafeUnderflow,ranknew1,flop=flop)
	! if(present(Flops))Flops= Flops + flop
	! deallocate(tau_Q)


	! allocate(QQ2(N,rmax))
	! call copymatT(matV(1:rmax,1:N),QQ2,rmax,N)
	! allocate (tau_Q(rmax))
	! allocate (jpvt2(rmax))
	! jpvt2=0
	! call geqp3modf90(QQ2,jpvt2,tau_Q,tolerance,SafeUnderflow,ranknew2,flop=flop)
	! if(present(Flops))Flops= Flops + flop
	! deallocate(tau_Q)

	! if(ranknew1>0 .and. ranknew2>0)then
		! allocate (RR1(ranknew1,rmax))
		! RR1=0d0
		! ! !$omp parallel do default(shared) private(i,j)
		! do j=1, ranknew1
			! do i=1, j
				! RR1(i,jpvt1(j))=QQ1(i,j)
			! enddo
		! enddo
		! ! !$omp end parallel do

		! allocate (RR2(ranknew2,rmax))
		! RR2=0d0
		! ! !$omp parallel do default(shared) private(i,j)
		! do j=1, ranknew2
			! do i=1, j
				! RR2(i,jpvt2(j))=QQ2(i,j)
			! enddo
		! enddo
		! ! !$omp end parallel do

		! allocate(mattemp(ranknew1,ranknew2))
		! mattemp=0
		! ! call zgemm('N','T',ranknew1,ranknew2,rmax, cone, RR1, ranknew1,RR2,ranknew2,czero,mattemp,ranknew1)
		! call gemmf90(RR1, ranknew1,RR2,ranknew2,mattemp,ranknew1,'N','T',ranknew1,ranknew2,rmax, cone,czero,flop=flop)


		! norm = fnorm(mattemp,rmax,rmax)

		! deallocate(mattemp,RR1,QQ1)
		! deallocate(QQ2,RR2)
		! deallocate(jpvt1,jpvt2)

	! else
		! norm=0
		! deallocate(QQ2,QQ1)
		! deallocate(jpvt1,jpvt2)
	! endif



	! allocate(QQ1(M,rmax))
	! QQ1 = matU(1:M,1:rmax)
	! allocate (tau_Q(rmax))
	! allocate (jpvt(rmax))
	! jpvt=0
	! call geqp3modf90(QQ1,jpvt,tau_Q,tolerance,SafeUnderflow,ranknew,flop=flop)
	! if(present(Flops))Flops= Flops + flop
	! ! call geqp3f90(QQ1,jpvt,tau_Q)
	! ranknew=rmax


	! if(ranknew>0)then
		! allocate (RR1(ranknew,rmax))
		! RR1=0d0
		! ! !$omp parallel do default(shared) private(i,j)
		! do j=1, ranknew
			! do i=1, j
				! RR1(i,j)=QQ1(i,j)
			! enddo
		! enddo
		! ! !$omp end parallel do

		! call un_or_gqrf90(QQ1,tau_Q,M,rmax,rmax,flop=flop)
		! if(present(Flops))Flops= Flops + flop

		! allocate(QQ2(rmax,N))
		! do ii=1,rmax
			! QQ2(ii,:) = matV(jpvt(ii),:)
		! enddo
		! allocate(mattemp(ranknew,N))
		! mattemp=0
		! ! call zgemm('N','N',ranknew,N,rmax, cone, RR1, ranknew,QQ2,rmax,czero,mattemp,ranknew)
		! call gemmf90(RR1, ranknew,QQ2,rmax,mattemp,ranknew,'N','N',ranknew,N,rmax, cone,czero,flop=flop)

		! norm = fnorm(mattemp,ranknew,rmax)

		! deallocate(mattemp,RR1)
		! deallocate(QQ2)
	! else
		! norm=0
	! endif

	! deallocate(tau_Q)
	! deallocate(jpvt)
	! deallocate(QQ1)


end subroutine LR_Fnorm




subroutine GetRank(M,N,mat,rank,eps,flop)
!
!
implicit none
integer M,N,rank,mn,i,mnl,ii,jj
DT::mat(:,:)
real(kind=8) eps
DT, allocatable :: UU(:,:), VV(:,:),Atmp(:,:),A_tmp(:,:),tau(:),mat1(:,:)
real(kind=8),allocatable :: Singular(:)
integer,allocatable :: jpvt(:)
real(kind=8),optional::flop

allocate(mat1(M,N))

mat1 = mat
mn=min(M,N)
mnl=max(M,N)
allocate(UU(M,mn))
allocate(VV(mn,N))
allocate(Singular(mn))


if(isnan(fnorm(mat,M,N)))then
	write(*,*)'input matrix NAN in GetRank'
	stop
end if


call gesvd_robust(mat1,Singular,UU,VV,M,N,mn,flop)
! write(*,*)Singular,'hh'
if(Singular(1)<SafeUnderflow)then
	rank = 1
	deallocate(UU,VV,Singular)
else
	if(isnan(sum(Singular)))then
		deallocate(UU,VV,Singular)
		write(*,*)'gesvd_robust wrong in GetRank, switching to QR'

		allocate(Atmp(mnl,mn))
		allocate(A_tmp(mn,mn))

		if(M>=N)then
			Atmp = mat
		else
			call copymatT(mat,Atmp,M,N)
		end if

		allocate(jpvt(mn))
		allocate(tau(mn))



		! RRQR
		jpvt = 0
		call geqp3f90(Atmp,jpvt,tau,flop)
		if(isnan(fnorm(Atmp,mnl,mn)))then
			write(*,*)'Q or R has NAN in GetRank'
			stop
		end if
		A_tmp = 0
		!$omp parallel do default(shared) private(ii,jj)
		do ii=1, mn
			do jj=ii, mn
				A_tmp(ii,jj)=Atmp(ii,jj)
			enddo
		enddo
		!$omp end parallel do


		rank = mn
		do i=1,mn
			if (abs(A_tmp(i,i))/abs(A_tmp(1,1))/mnl<=eps) then
				rank=i
				if(abs(A_tmp(i,i))<SafeUnderflow)rank = i -1
				exit
			end if
		end do

		deallocate(jpvt)
		deallocate(tau)
		deallocate(Atmp)
		deallocate(A_tmp)


	else

		! write(*,*)Singular,'hh'

		rank = mn
		do i=1,mn
			if (Singular(i)/Singular(1)<=eps .or. Singular(i)<1d-60) then
				rank=i
				if(Singular(i)<Singular(1)*eps/10 .or. Singular(i)<1d-60)rank = i -1
				exit
			end if
		end do


		deallocate(UU,VV,Singular)

	end if
endif
deallocate(mat1)


end subroutine GetRank



subroutine PComputeRange(M_p,N,mat,rank,eps,ptree,pgno,Flops)

implicit none
integer M,N,M_loc,rank,mn,i,mnl,ii,jj,rrflag
DT::mat(:,:)
real(kind=8) eps
DT, allocatable :: UU(:,:), VV(:,:),Atmp(:,:),A_tmp(:,:),tau(:),mat1D(:,:),mat2D(:,:)
real(kind=8),allocatable :: Singular(:)
real(kind=8) :: flop
integer,allocatable :: jpvt(:)
integer, allocatable :: ipiv(:),jpiv(:),JPERM(:)
integer:: M_p(:,:)
integer,allocatable:: M_p_1D(:,:)
type(proctree)::ptree
integer pgno,proc,nproc,nb1Dc,nb1Dr,ctxt1D,ctxt,idxs_o,idxe_o,ierr
integer myArows,myAcols,info,nprow,npcol,myrow,mycol,taun
integer::descsMat1D(9),descsMat2D(9)
real(kind=8),optional:: Flops

rank=0

if(present(Flops))Flops=0d0

if(IOwnPgrp(ptree,pgno))then


	proc = ptree%MyID - ptree%pgrp(pgno)%head
	nproc = ptree%pgrp(pgno)%nproc
	M_loc = M_p(proc+1,2)-M_p(proc+1,1)+1
	M = M_p(nproc,2)

	if(nproc==1)then
		call ComputeRange(M,N,mat,rank,1,eps,Flops=flop)
		if(present(Flops))Flops = Flops + flop
	else

		mn=min(M,N)


		!!!!**** generate 2D grid blacs quantities
		ctxt = ptree%pgrp(pgno)%ctxt
		call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
		if(myrow/=-1 .and. mycol/=-1)then
			myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
			myAcols = numroc_wp(N, nbslpk, mycol, 0, npcol)
			call descinit( descsMat2D, M, N, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
			allocate(mat2D(myArows,myAcols))
			mat2D=0
		else
			descsMat2D(2)=-1
			allocate(mat2D(1,1))
			mat2D=0
		endif


		!!!!**** redistribution of input matrix
		call Redistribute1Dto2D(mat,M_p,0,pgno,mat2D,M,0,pgno,N,ptree)


		! Compute RRQR
		if(myrow/=-1 .and. mycol/=-1)then
			allocate(ipiv(myAcols))
			ipiv=0
			taun=numroc_wp(mn, nbslpk, mycol, 0, npcol)
			allocate(tau(taun))
			tau=0
			allocate(jpiv(N))
			jpiv=0
			allocate(JPERM(N))
			JPERM=0
			call pgeqpfmodf90(M, N, mat2D, 1, 1, descsMat2D, ipiv, tau, JPERM, jpiv, rank,eps, SafeUnderflow,flop=flop)
			if(present(Flops))Flops = Flops + flop/dble(nprow*npcol)

			if(rank>0)then
				call pun_or_gqrf90(ctxt,mat2D,tau,M,rank,rank,descsMat2D,1,1,flop=flop)
				if(present(Flops))Flops = Flops + flop/dble(nprow*npcol)
			else
				rank=1
				if(myArows>0 .and. myAcols>0)mat2D=0d0
			endif

			deallocate(ipiv)
			deallocate(tau)
			deallocate(jpiv)
			deallocate(JPERM)

			if(isnan(fnorm(mat2D,myArows,myAcols)))then
				write(*,*)'Q or R has NAN in PComputeRange'
				stop
			end if
		endif

		call MPI_ALLREDUCE(MPI_IN_PLACE, rank, 1,MPI_integer, MPI_MAX, ptree%pgrp(pgno)%Comm,ierr)

		!!!!**** redistribution of output matrix
		call Redistribute2Dto1D(mat2D,M,0,pgno,mat,M_p,0,pgno,N,ptree)


		deallocate(mat2D)
	endif
endif

end subroutine PComputeRange




subroutine ComputeRange(M,N,mat,rank,rrflag,eps,Flops)

implicit none
integer M,N,rank,mn,i,mnl,ii,jj,rrflag
DT::mat(:,:)
real(kind=8) eps
DT, allocatable :: UU(:,:), VV(:,:),Atmp(:,:),A_tmp(:,:),tau(:)
real(kind=8),allocatable :: Singular(:)
integer,allocatable :: jpvt(:)
real(kind=8),optional:: Flops
real(kind=8):: flop
if(present(Flops))Flops=0d0

	mn=min(M,N)
	if(isnan(fnorm(mat,M,N)))then
		write(*,*)'input matrix NAN in ComputeRange'
		stop
	end if


	allocate(jpvt(N))
	jpvt = 0
	allocate(tau(mn))
	if(rrflag==1)then
		! RRQR
		call geqp3modf90(mat,jpvt,tau,eps,SafeUnderflow,rank,flop=flop)
		if(present(Flops))Flops = Flops + flop

		if(rank>0)then
			call un_or_gqrf90(mat,tau,M,rank,rank,flop=flop)
			if(present(Flops))Flops = Flops + flop
		else
			rank=1
			mat(1:M,1)=0d0
		endif
	else
		call geqp3f90(mat,jpvt,tau,flop=flop)
		if(present(Flops))Flops = Flops + flop
		call un_or_gqrf90(mat,tau,M,N,mn,flop=flop)
		if(present(Flops))Flops = Flops + flop
		rank=mn
	endif
	if(isnan(fnorm(mat,M,N)))then
		write(*,*)'Q or R has NAN in ComputeRange'
		stop
	end if


	deallocate(jpvt)
	deallocate(tau)

end subroutine ComputeRange



subroutine CheckRandomizedLR(M,N,mat,tolerance)


implicit none
integer M,N,rank,mn,i,j,k,rankmax_c,rankmax_r,rankmax_min,flag0,rank_new,rmax
DT::mat(M,N),ctemp
DT,allocatable::mat1(:,:),mat2(:,:)
real(kind=8) tolerance,Smax
DT, allocatable :: UU(:,:), VV(:,:),matrix_U(:,:), matrix_V(:,:),U_new(:,:),V_new(:,:),U_new1(:,:),V_new1(:,:),test_in(:,:),test_out1(:,:),test_out2(:,:),test_out3(:,:)
real(kind=8),allocatable :: Singular(:)
integer,allocatable:: select_row(:), select_column(:)
DT,allocatable::MatrixSubselection(:,:)

allocate(mat1(M,N))
allocate(mat2(M,N))
call GetRank(M,N,mat,rank,tolerance)
rankmax_r = rank*3
write(*,*)rankmax_r,min(m,n)
if(rankmax_r>min(m,n))rankmax_r = min(m,n)
rankmax_c = rankmax_r
rankmax_min = min(rankmax_r,rankmax_c)

! rankmax_r = M
! rankmax_c = N
! rankmax_min = min(rankmax_r,rankmax_c)






! method 1: SeudoSkeleton


allocate(select_row(rankmax_r),select_column(rankmax_c))
call linspaceI(1,M,rankmax_r,select_row)
call linspaceI(1,N,rankmax_c,select_column)

allocate (MatrixSubselection(rankmax_r,rankmax_c))
MatrixSubselection = mat(select_row,select_column)

allocate (matrix_U(M,rankmax_c),matrix_V(N,rankmax_r))
matrix_U = mat(1:M,select_column)
do i=1,N
do j=1,rankmax_r
matrix_V(i,j) = mat(select_row(j),i)
end do
end do


deallocate(select_column,select_row)


 ! write(*,*)tolerance
allocate (UU(rankmax_r,rankmax_min),VV(rankmax_min,rankmax_c),Singular(rankmax_min))
call gesvd_robust(MatrixSubselection,Singular,UU,VV,rankmax_r,rankmax_c,rankmax_min)
deallocate (MatrixSubselection)

flag0=0 ; i=0
do while (flag0==0 .and. i<rankmax_min)
	i=i+1
	if (Singular(i)/Singular(1)<=tolerance/10) then
		flag0=1
	endif
enddo

rank_new=i
! if(rank_new>rank)rank_new=rank

write(*,*)'rank=',rank,'rank_new',rank_new

! allocate(U_new1(M,rank_new))
! allocate(V_new1(rank_new,N))

! do i =1,M
! do j =1,rank_new
	! U_new1(i,j) = UU(i,j)*Singular(j)
! end do
! end do
! V_new1 = VV(1:rank_new,1:N)

! call gemm_omp(U_new1,V_new1,mat2,M,N,rank_new)
! write(*,*)Singular(1:rank_new+2)
Smax = Singular(1)

allocate(U_new(M,rank_new))
allocate(V_new(rank_new,N))

do j=1, rank_new
	do i=1, M
		ctemp=0
		do k=1, rankmax_c
			ctemp=ctemp+matrix_U(i,k)*conjg(cmplx(VV(j,k),kind=8))
		enddo
		U_new(i,j)=ctemp
	enddo
enddo


do j=1, rank_new
	do i=1, N
		ctemp=0
		do k=1, rankmax_r
			ctemp=ctemp+conjg(cmplx(UU(k,j),kind=8))*matrix_V(i,k)
		enddo
		V_new(j,i)=ctemp/Singular(j)
	enddo
enddo

deallocate (matrix_U,VV)
deallocate (matrix_V,UU,Singular)

! call gemm_omp(U_new,V_new,mat1,M,N,rank_new)
call gemmf90(U_new,M,V_new,rank_new,mat1,M,'N','N',M,N,rank_new,cone,czero)

write(*,*)M,N
do j=1,N
do i=1,M
	write(211,*)dble(mat1(i,j))
	write(212,*)aimag(cmplx(mat1(i,j),kind=8))
end do
end do


! write(*,*)'F-norm residual:', fnorm(mat-mat1,M,N)/fnorm(mat,M,N),fnorm(mat-mat2,M,N)/fnorm(mat,M,N)
write(*,*)'F-norm residual:', fnorm(mat-mat1,M,N)/fnorm(mat,M,N)
allocate(test_in(N,1))
allocate(test_out1(M,1))
allocate(test_out2(M,1))
! allocate(test_out3(M,1))
do i=1,N
	call random_dp_number(test_in(i,1))
end do

! call gemm_omp(mat,test_in,test_out1,M,1,N)
call gemmf90(mat,M,test_in,N,test_out1,M,'N','N',M,1,N,cone,czero)
! call gemm_omp(mat1,test_in,test_out2,M,1,N)
call gemmf90(mat1,M,test_in,N,test_out2,M,'N','N',M,1,N,cone,czero)

! call gemm_omp(mat2,test_in,test_out3,M,1,N)
! write(*,*)'testing vector error:', fnorm(test_out1-test_out2,M,1)/fnorm(test_out1,M,1),fnorm(test_out1-test_out3,M,1)/fnorm(test_out1,M,1)
write(*,*)'testing vector error:', fnorm(test_out1-test_out2,M,1)/fnorm(test_out1,M,1)






! method 2: ACA

rmax = min(M,N)
allocate(matrix_U(M,rmax))
allocate(matrix_V(rmax,N))

call ACA_CompressionFull(mat,matrix_U,matrix_V,M,N,rmax,rank_new,tolerance*0.3d0,tolerance)

! call gemm_omp(matrix_U(1:M,1:rank_new),matrix_V(1:rank_new,1:N),mat2,M,N,rank_new)
call gemmf90(matrix_U,M,matrix_V,rmax,mat2,M,'N','N',M,N,rank_new,cone,czero)


write(*,*)'F-norm residual:', fnorm(mat-mat2,M,N)/fnorm(mat,M,N),' rank:',rank_new


deallocate(mat1)
deallocate(mat2)

end subroutine CheckRandomizedLR



! generate random permutation of 1:N
subroutine rperm(N, p)

integer, intent(in):: N
integer:: p(N)

integer:: i
integer:: k, j, ipj, itemp, m
real(kind=8), dimension(100) :: u

call assert(N<1d9,'In rperm, N too large')
! write(*,*)p(1)
p = (/ (i, i=1,N) /)

! Generate up to 100 U(0,1) numbers at a time.
do i=1,N,100
  m = min(N-i+1, 100)
  call random_number(u)
  do j=1,m
	ipj = i+j-1
	k = int(u(j)*(N-ipj+1)) + ipj
	itemp = p(ipj)
	p(ipj) = p(k)
	p(k) = itemp
  end do
end do
return

end subroutine rperm



subroutine init_random_seed()
implicit none
integer, allocatable :: seed(:)
integer :: i, n, un, istat, dt(8), pid
integer*8 :: t

call random_seed(size = n)
allocate(seed(n))
! First try if the OS provides a random number generator

! Fallback to XOR:ing the current time and pid. The PID is
! useful in case one launches multiple instances of the same
! program in parallel.
call system_clock(t)
if (t == 0) then
  call date_and_time(values=dt)
  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
	   + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
	   + dt(3) * 24_int64 * 60 * 60 * 1000 &
	   + dt(5) * 60 * 60 * 1000 &
	   + dt(6) * 60 * 1000 + dt(7) * 1000 &
	   + dt(8)
end if

! pid = getpid()
pid = 0

t = ieor(t, int(pid, kind(t)))
do i = 1, n
  seed(i) = lcg(t)
end do
call random_seed(put=seed)
contains
! This simple PRNG might not be good enough for real work, but is
! sufficient for seeding a better PRNG.
function lcg(s)
  integer :: lcg
  integer(int64) :: s
  if (s == 0) then
	 s = 104729
  else
	 s = mod(s, 4294967296_int64)
  end if
  s = mod(s * 279470273_int64, 4294967291_int64)
  lcg = int(mod(s, int(huge(0), int64)), kind(0))
end function lcg
end subroutine init_random_seed


subroutine random_dp_number(val)
implicit none
real(kind=8):: a=0,b=0,c=0,d=0
integer seed
class(*) val

select type(val)
type is (complex(kind=8))

	! Uniform distribution
	call random_number(a)
	call random_number(b)
	call random_number(c)
	call random_number(d)
	if (c<0.5d0) then
		a=-a
	endif
	if (d<0.5d0) then
		b=-b
	endif
	val=a+junit*b


	! ! ! Normal distribution
	! ! call random_number(a)
	! ! seed = a*10000000
	! ! random_dp_number =  c8_normal_01 ( seed )
type is (real(kind=8))
	! Uniform distribution
	call random_number(a)
	val = a*2d0-1d0

end select
return
end subroutine random_dp_number



complex(kind=8) function c8_normal_01(seed)

!*****************************************************************************80
!
!! C8_NORMAL_01 returns a unit pseudonormal C8.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
!    Output, complex (kind = 8 ) C8_NORMAL_01, a sample of the PDF.
!
!    Output, integer ( kind = 4 ) SEED, a seed for the random number
!    generator.
!
  implicit none

  real(kind=8), parameter :: r8_pi = 3.141592653589793D+00
  ! real(kind=8) r8_uniform_01
  integer ( kind = 4 ) seed
  real(kind=8) v1
  real(kind=8) v2
  real(kind=8) x_c
  real(kind=8) x_r

  v1 = r8_uniform_01 ( seed )
  v2 = r8_uniform_01 ( seed )

  x_r = sqrt ( - 2.0D+00 * log ( v1 ) ) * cos ( 2.0D+00 * r8_pi * v2 )
  x_c = sqrt ( - 2.0D+00 * log ( v1 ) ) * sin ( 2.0D+00 * r8_pi * v2 )

  c8_normal_01 = cmplx ( x_r, x_c, kind = 8 )

  return
end function c8_normal_01

real(kind=8) function r8_normal_01(seed)
!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real(kind=8) R8_NORMAL_01, a normally distributed
!    random value.
!
  implicit none

  real(kind=8) r1
  real(kind=8) r2

  real(kind=8), parameter :: r8_pi = 3.141592653589793D+00
  ! real(kind=8) r8_uniform_01
  integer ( kind = 4 ) seed
  real(kind=8) x

  r1 = r8_uniform_01 ( seed )
  r2 = r8_uniform_01 ( seed )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  r8_normal_01 = x

  return
end function r8_normal_01


real(kind=8) function r8_uniform_01(seed)

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real(kind=8) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end function r8_uniform_01


! ****************************************************************************************** !
! ***************                          assert function                               !*********** !
! ***************                    									          !*********** !
! ****************************************************************************************** !
! output an error message is the first argument is false
subroutine assert(statement,msg)
	logical::statement
	character(*)::msg
	if(.not. statement)then
		write(*,*)msg
		stop
	end if
end subroutine assert

function floor_safe(input)
real(kind=8)::input
integer floor_safe
integer input_nint
input_nint = NINT(input)
if(abs(input_nint-input)<1d-13)then
	floor_safe = input_nint
else
	floor_safe = floor(input)
end if

end function floor_safe

function ceiling_safe(input)
real(kind=8)::input
integer ceiling_safe
integer input_nint
input_nint = NINT(input)
if(abs(input_nint-input)<1d-13)then
	ceiling_safe = input_nint
else
	ceiling_safe = ceiling(input)
end if

end function ceiling_safe

function INT_safe(input)
real(kind=8)::input
integer INT_safe
integer input_nint
input_nint = NINT(input)
if(abs(input_nint-input)<1d-13)then
	INT_safe = input_nint
else
	INT_safe = INT(input)
end if
end function INT_safe


subroutine cscalar(a,b,c)
	  implicit none
	real(kind=8) b(3)
	complex(kind=8) a(3),c
	complex(kind=8) ax,ay,az
	real(kind=8) bx,by,bz
	ax=a(1);ay=a(2);az=a(3)
	bx=b(1);by=b(2);bz=b(3)
	c=ax*bx+ay*by+az*bz
	  return
end subroutine cscalar

subroutine scalar(a,b,c)
	  implicit none
	real(kind=8) a(3),b(3),c
	real(kind=8) ax,ay,az
	real(kind=8) bx,by,bz
	ax=a(1);ay=a(2);az=a(3)
	bx=b(1);by=b(2);bz=b(3)
	c=ax*bx+ay*by+az*bz
	  return
end subroutine scalar


subroutine curl(a,b,c)
	  implicit none
	real(kind=8) a(3),b(3),c(3)
	real(kind=8) ax,ay,az
	real(kind=8) bx,by,bz
	real(kind=8) cx,cy,cz
	ax=a(1);ay=a(2);az=a(3)
	bx=b(1);by=b(2);bz=b(3)

	cx=ay*bz-by*az
	cy=-ax*bz+bx*az
	cz=ax*by-bx*ay
	c(1)=cx;c(2)=cy;c(3)=cz
	  return
end subroutine curl


subroutine ccurl(a,b,c)
	  implicit none
	real(kind=8) a(3)
	complex(kind=8) b(3),c(3)
	real(kind=8) ax,ay,az
	complex(kind=8) bx,by,bz
	complex(kind=8) cx,cy,cz
	ax=a(1);ay=a(2);az=a(3)
	bx=b(1);by=b(2);bz=b(3)

	cx=ay*bz-by*az
	cy=-ax*bz+bx*az
	cz=ax*by-bx*ay
	c(1)=cx;c(2)=cy;c(3)=cz
	  return
end subroutine ccurl

 real(kind=8) function norm_vector(vector,n)
	  implicit none
	  integer n
	  DT:: vector(n)
	  integer i
	  real(kind=8) sum

	  sum=0.0d0
      do i=1,n
	      sum=sum+dble(vector(i)*conjg(cmplx(vector(i),kind=8)))
      enddo

      norm_vector=sum

      return
end function norm_vector


subroutine LeastSquare(m,n,k,A,b,x,eps_r,Flops)
	!
	!
	implicit none

	integer m,n,k,mn_min
	DT:: A(m,n),x(n,k),b(m,k)
	real(kind=8),allocatable::Singular(:)
	DT,allocatable:: Atmp(:,:),xtmp(:,:),tau(:),A_tmp(:,:),UU(:,:),VV(:,:),UU_h(:,:),VV_h(:,:),VV_h2(:,:),z_arb(:,:),x_arb(:,:),matrixtemp(:,:),A_tmp_rank(:,:),xtmp_rank(:,:),xtmp_rank3(:,:),A_tmp_rank2(:,:),xtmp_rank2(:,:)
	real(kind=8):: eps_r
	integer ii,jj,i,j,flag0,rank
	integer,allocatable:: jpvt(:)
	real(kind=8),optional::Flops
	real(kind=8)::flop

	DT:: alpha,beta
	alpha=1d0
	beta=0d0

	if(present(Flops))Flops=0

	allocate(xtmp(n,k))
	allocate(tau(n))
	allocate(jpvt(n))
	allocate(A_tmp(n,n))
	allocate(Atmp(m,n))
	allocate(UU(m,n))
	allocate(VV(n,n))
	allocate(Singular(n))

	if(m<n)write(*,*)m,n
	call assert(m>=n,'m should not be less than n for least square')


	if(fnorm(b,m,k)<SafeUnderflow)then
		write(*,*)'warning: RHS zero in least square. |b|= ',fnorm(b,m,k),'size b: ',m,k,'size A',m,n
		x = 0
	else

	Atmp = A


	! SVD
    call gesvd_robust(Atmp,Singular,UU,VV,m,n,n,flop)
	if(present(Flops))Flops=Flops+flop

! !!!!!!!  If SVD fails, uncomment the following If statement, but the code might become slow
	! if(isnan(sum(Singular)))then

		! write(*,*)'gesvd wrong in LeastSquare, switching to QR'

		! ! call GetRank(m,n,Atmp,rank,Rank_detection_factor)
		! ! write(*,*)rank,'kao kao'

		! ! stop
		! Atmp = A


		! ! RRQR
		! jpvt = 0
		! call geqp3f90(Atmp,jpvt,tau)
		! if(isnan(fnorm(Atmp,m,n)))then
			! write(*,*)'Q or R has NAN in LeastSquare'
			! stop
		! end if
		! call un_or_mqrf90(Atmp,tau,b,'L','C',m,n,n)
		! A_tmp = 0
		! ! !$omp parallel do default(shared) private(ii,jj)
		! do ii=1, n
			! do jj=ii, n
				! A_tmp(ii,jj)=Atmp(ii,jj)
			! enddo
		! enddo
		! ! !$omp end parallel do
		! ! !$omp parallel do default(shared) private(ii,jj)
		! do jj=1, k
			! do ii=1, n
				! xtmp(ii,jj)=b(ii,jj)
			! enddo
		! enddo
		! ! !$omp end parallel do

		! flag0=0 ; i=0

		! rank = n
		! do i=1,n
			! if (abs(A_tmp(i,i))/abs(A_tmp(1,1))/m<=eps_r) then
				! rank=i
				! if(abs(A_tmp(i,i))<SafeUnderflow)rank = i -1
				! ! write(*,*)rank,n
				! exit
			! end if
		! end do

		! allocate(A_tmp_rank(rank,rank))
		! do jj=1, rank
			! do ii=1, rank
				! A_tmp_rank(ii,jj)=A_tmp(ii,jj)
			! enddo
		! enddo

		! allocate(xtmp_rank(rank,k))
		! do ii = 1,rank
			! do jj =1,k
				! xtmp_rank(ii,jj) = xtmp(ii,jj)
			! end do
		! end do


		! call trsmf90(A_tmp_rank,xtmp_rank,'L','U','N','N',rank,rank)

		! xtmp = 0

		! ! do ii = 1,n
		! ! do jj =1,k
			! ! xtmp(ii,jj) = random_dp_number()
		! ! end do
		! ! end do

		! do ii = 1,rank
			! do jj =1,k
				! xtmp(ii,jj) = xtmp_rank(ii,jj)
			! end do
		! end do

		! do ii=1,n
			! x(jpvt(ii),1:k) = xtmp(ii,1:k)
		! end do

		! if(isnan(fnorm(x,n,k)))then
			! write(*,*)'trisolve has NAN in LeastSquare'
			! stop
		! end if

		! deallocate(A_tmp_rank,xtmp_rank)
	! else
		if(Singular(1)<SafeUnderflow)then
			write(*,*)'warning: Matrix zero in least square'
			rank=1
			x = 0
		else
			rank = n
			do i=1,n
				if (Singular(i)/Singular(1)<=eps_r) then
					rank=i
					if(Singular(i)<Singular(1)*eps_r/10)rank = i -1
					exit
				end if
			end do


			allocate(UU_h(rank,m))
			allocate(VV_h(n,rank))
			allocate(matrixtemp(rank,k))

			do ii=1,rank
			do jj =1,m
				UU_h(ii,jj) = conjg(cmplx(UU(jj,ii),kind=8))
			end do
			end do

			do ii=1,n
			do jj=1,rank
				VV_h(ii,jj) = conjg(cmplx(VV(jj,ii),kind=8))/Singular(jj)
			end do
			end do

			call gemmf90(UU_h, rank, b,m, matrixtemp,rank,'N','N',rank,k,m,alpha,beta,flop)
			if(present(Flops))Flops=Flops+flop
			call gemmf90(VV_h, n, matrixtemp,rank, x,n,'N','N',n,k,rank,alpha,beta,flop)
			if(present(Flops))Flops=Flops+flop

			deallocate(UU_h,VV_h,matrixtemp)
		end if
	! end if


	do ii = 1,k
	if(isnan(abs(sum(x(:,ii)))))then
		! do jj =1,rank
		! write(*,*)jj,'hh',A_tmp_rank(jj,:)
		! end do
		write(*,*)'hh',rank,Singular,fnorm(A,m,n)
		stop
	end if
	end do

	! deallocate(A_tmp_rank,xtmp_rank)



	 ! A = Atmp

	 end if

	deallocate(Singular)
	deallocate(xtmp)
	deallocate(tau)
	deallocate(jpvt)
	deallocate(A_tmp)
	deallocate(Atmp)
	deallocate(UU)
	deallocate(VV)

end subroutine LeastSquare




subroutine GeneralInverse(m,n,A,A_inv,eps_r,Flops)
	!
	!
	implicit none

	integer m,n,mn_min
	DT:: A(:,:),A_inv(:,:)
	real(kind=8),allocatable::Singular(:)
	DT,allocatable:: Atmp(:,:),tau(:),UU(:,:),VV(:,:),UU_h(:,:),VV_h(:,:),matrixtemp(:,:),A_tmp_rank(:,:),xtmp_rank(:,:),xtmp_rank3(:,:),A_tmp_rank2(:,:),xtmp_rank2(:,:)
	real(kind=8):: eps_r
	integer ii,jj,i,j,flag0,rank
	integer,allocatable:: jpvt(:)
	real(kind=8),optional::Flops
	real(kind=8)::flop

	DT:: alpha,beta
	alpha=1d0
	beta=0d0

	if(present(Flops))Flops=0

	A_inv=0

	allocate(Atmp(m,n))
	Atmp=0
	! SVD
	Atmp = A
	mn_min = min(m,n)
	allocate(Singular(mn_min))
	allocate(UU(m,mn_min))
	allocate(VV(mn_min,n))

	Singular=0
	UU=0
	VV=0

    call gesvd_robust(Atmp,Singular,UU,VV,m,n,mn_min,flop=flop)
	if(present(Flops))Flops=Flops+flop

	if(Singular(1)<SafeUnderflow)then
		write(*,*)'Warning: A zero in GeneralInverse'
		A_inv = 0
	else
		rank = mn_min
		do i=1,mn_min
			if (Singular(i)/Singular(1)<=eps_r) then
				rank=i
				if(Singular(i)<Singular(1)*eps_r/10)rank = i -1
				exit
			end if
		end do


		allocate(UU_h(rank,m))
		allocate(VV_h(n,rank))
		UU_h=0
		VV_h=0

		do ii=1,rank
		do jj =1,m
			UU_h(ii,jj) = conjg(cmplx(UU(jj,ii),kind=8))
		end do
		end do

		do ii=1,n
		do jj=1,rank
			VV_h(ii,jj) = conjg(cmplx(VV(jj,ii),kind=8))/Singular(jj)
		end do
		end do

		call gemmf90(VV_h,n, UU_h,rank, A_inv,n,'N','N',n,m,rank,alpha,beta,flop=flop)
		if(present(Flops))Flops=Flops+flop

		deallocate(UU_h,VV_h)
	endif

	deallocate(Singular)
	! deallocate(jpvt)
	deallocate(UU)
	deallocate(VV)
	deallocate(Atmp)

end subroutine GeneralInverse


	subroutine RandomizedSVD(matRcol,matZRcol,matRrow,matZcRrow,matU,matV,Singular,rankmax_r,rankmax_c,rmax,rank,tolerance,SVD_tolerance,Flops)
	!
	!
    use BPACK_DEFS
    implicit none

    integer i, j, ii, jj, indx, rank_1, rank_2
    integer blocks, index_j, group_nn, rank_blocks
    integer group_m,group_n,size_of_groupm,size_of_groupn
    real(kind=8) norm_Z,norm_U,norm_V,tolerance,SVD_tolerance
    integer edgefine_m,edgefine_n, level_butterfly, level_blocks
    integer edge_m, edge_n, header_m, header_n
    integer rank,rank1,rank2,rank12, ranknew, row, column, rankmax,rankmax_c,rankmax_r,rankmax_min,rmax,idxs_r,idxs_c
    DT value_Z,value_UV,maxvalue
    DT inner_U,inner_V,ctemp
    real(kind=8) inner_UV
    integer,allocatable:: select_column(:), select_row(:)
	DT::matU(rankmax_r,rmax),matV(rmax,rankmax_c),matRcol(rankmax_c,rmax),matZRcol(rankmax_r,rmax),matRrow(rankmax_r,rmax),matZcRrow(rankmax_c,rmax)
	DT,allocatable::matZRcol1(:,:),matZcRrow1(:,:),tau(:)
	DT,allocatable::QQ1(:,:),RR1(:,:),QQ2(:,:),RR2(:,:),RrowcZRcol(:,:)
	real(kind=8)::Singular(rmax)
    DT,allocatable:: row_R(:),column_R(:),matM(:,:),matrixtemp(:,:)
    real(kind=8),allocatable:: norm_row_R(:),norm_column_R(:)

	DT, allocatable :: RrowcQ1(:,:),RrowcQ1inv(:,:),Q2cRcol(:,:),Q2cRcolinv(:,:), QQ2tmp(:,:), RR2tmp(:,:), UUsml(:,:), VVsml(:,:),tau_Q(:),mattemp(:,:),matU1(:,:),matV1(:,:)
	real(kind=8), allocatable :: Singularsml(:)
	integer,allocatable:: jpvt(:)
	real(kind=8),optional::Flops
	real(kind=8)::flop


	if(present(Flops))Flops=0

	allocate(matZRcol1(rankmax_r,rmax))
	allocate(matZcRrow1(rankmax_c,rmax))
	allocate(tau(rmax))
	allocate(jpvt(rmax))
	allocate(QQ1(rankmax_r,rmax))
	allocate(RR1(rmax,rmax))
	allocate(QQ2(rankmax_c,rmax))
	allocate(RR2(rmax,rmax))
	allocate(RrowcZRcol(rmax,rmax))


	rankmax_min = min(rankmax_r,rankmax_c)
    matU = 0
	matV = 0
	Singular = 0

	QQ1 = matZRcol
	jpvt = 0
	tau=0
	call geqp3f90(QQ1,jpvt,tau,flop=flop)
	if(present(Flops))Flops=Flops+flop
	!  call geqrff90(QQ1,tau)
	RR1=0d0
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rmax
		do i=1, j
			RR1(i,j)=QQ1(i,j)
		enddo
	enddo
	! !$omp end parallel do
	call un_or_gqrf90(QQ1,tau,rankmax_r,rmax,rmax,flop=flop)
	if(present(Flops))Flops=Flops+flop

	! write(*,*)fnorm(QQ1,rankmax_r,rmax),rankmax_r,rmax,fnorm(RR1,rmax,rmax),rmax,rmax,'really'

	if(abs(RR1(1,1))<SafeUnderflow)then
		rank1=1
	else
		rank1 = rmax
		do i=1,rmax
			if (abs(RR1(i,i))/abs(RR1(1,1))/rankmax_r<=tolerance) then   ! be careful with the tolerance here
				rank1=i
				if(abs(RR1(i,i))<SafeUnderflow)rank1 = i -1
				exit
			end if
		end do
	endif
! write(*,*)shape(QQ2),shape(matZcRrow)

	QQ2 = matZcRrow
	jpvt = 0
	tau = 0
	call geqp3f90(QQ2,jpvt,tau,flop=flop)
	if(present(Flops))Flops=Flops+flop
	!  call geqrff90(QQ2,tau)
	RR2=0d0
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rmax
		do i=1, j
			RR2(i,j)=QQ2(i,j)
		enddo
	enddo
	! !$omp end parallel do
	call un_or_gqrf90(QQ2,tau,rankmax_c,rmax,rmax,flop=flop)
	if(present(Flops))Flops=Flops+flop

	if(abs(RR2(1,1))<SafeUnderflow)then
		rank2=1
	else
		rank2 = rmax
		do i=1,rmax
			if (abs(RR2(i,i))/abs(RR2(1,1))/rankmax_c<=tolerance) then  ! be careful with the tolerance here
				rank2=i
				if(abs(RR2(i,i))<SafeUnderflow)rank2 = i -1
				exit
			end if
		end do
	endif
	! write(*,*)rank2,rank1,rmax,'ha',rankmax_r,rankmax_c
	! write(111,*)QQ1,QQ2
	! stop

	!  do ii=1,rankmax_r
	!  do jj=1,rmax
		!  write(133,*)dble(QQ1(ii,jj)),aimag(QQ1(ii,jj))
	!  enddo
	!  enddo

	!  do ii=1,rankmax_c
	!  do jj=1,rmax
		!  write(134,*)dble(QQ2(ii,jj)),aimag(QQ2(ii,jj))
	!  enddo
	!  enddo

	allocate(RrowcQ1(rmax,rank1))
	RrowcQ1=0
	allocate(RrowcQ1inv(rank1,rmax))
	RrowcQ1inv=0
	! call gemmHN_omp(matRrow,QQ1(1:rankmax_r,1:rank1),RrowcQ1,rmax,rank1,rankmax_r)
	call gemmf90(matRrow,rankmax_r,QQ1,rankmax_r,RrowcQ1,rmax,'C','N',rmax,rank1,rankmax_r, cone,czero,flop=flop)
	if(present(Flops))Flops=Flops+flop

	call GeneralInverse(rmax,rank1,RrowcQ1,RrowcQ1inv,tolerance,Flops=flop)
	if(present(Flops))Flops=Flops+flop
	deallocate(RrowcQ1)

	! allocate(ipiv(rmax))
	! call getrff90(RrowcQ1,ipiv)
	! call getrif90(RrowcQ1,ipiv)
	! RrowcQ1inv = RrowcQ1
	! deallocate(ipiv)
	! deallocate(RrowcQ1)





	allocate(Q2cRcol(rank2,rmax))
	Q2cRcol=0
	allocate(Q2cRcolinv(rmax,rank2))
	Q2cRcolinv=0
	! call gemmHN_omp(QQ2(1:rankmax_c,1:rank2),matRcol,Q2cRcol,rank2,rmax,rankmax_c)
	call gemmf90(QQ2,rankmax_c,matRcol,rankmax_c,Q2cRcol,rank2,'C','N',rank2,rmax,rankmax_c, cone,czero,flop=flop)
	if(present(Flops))Flops=Flops+flop
	call GeneralInverse(rank2,rmax,Q2cRcol,Q2cRcolinv,tolerance,Flops=flop)
	if(present(Flops))Flops=Flops+flop
	deallocate(Q2cRcol)

	! allocate(ipiv(rmax))
	! call getrff90(Q2cRcol,ipiv)
	! call getrif90(Q2cRcol,ipiv)
	! Q2cRcolinv = Q2cRcol
	! deallocate(ipiv)
	! deallocate(Q2cRcol)


	! call gemmHN_omp(matRrow,matZRcol,RrowcZRcol,rmax,rmax,rankmax_r)
	call gemmf90(matRrow,rankmax_r,matZRcol,rankmax_r,RrowcZRcol,rmax,'C','N',rmax,rmax,rankmax_r, cone,czero,flop=flop)
	if(present(Flops))Flops=Flops+flop

	allocate(matrixtemp(rmax,rank2))
	matrixtemp=0
	allocate(matM(rank1,rank2))
	matM=0
	! call gemm_omp(RrowcZRcol,Q2cRcolinv,matrixtemp,rmax,rank2,rmax)
	call gemmf90(RrowcZRcol,rmax,Q2cRcolinv,rmax,matrixtemp,rmax,'N','N',rmax,rank2,rmax,cone,czero,flop=flop)
	if(present(Flops))Flops=Flops+flop

	! call gemm_omp(RrowcQ1inv,matrixtemp,matM,rank1,rank2,rmax)
	call gemmf90(RrowcQ1inv,rank1,matrixtemp,rmax,matM,rank1,'N','N',rank1,rank2,rmax,cone,czero,flop=flop)
	if(present(Flops))Flops=Flops+flop

	deallocate(matrixtemp,RrowcQ1inv,Q2cRcolinv)

	! allocate(matrixtemp(rankmax_r,rank2))
	! allocate(matM(rank1,rank2))
	! call gemm_omp(matZRcol,Q2cRcolinv,matrixtemp,rankmax_r,rank2,rmax)
	! call gemmHN_omp(QQ1(1:rankmax_r,1:rank1),matrixtemp,matM,rank1,rankmax_r,rank2)
	! deallocate(matrixtemp,RrowcQ1inv,Q2cRcolinv)

	rank12 = min(rank1,rank2)
	allocate(UUsml(rank1,rank12),VVsml(rank12,rank2),Singularsml(rank12))
	UUsml=0
	VVsml=0
	Singularsml=0
	call SVD_Truncate(matM,rank1,rank2,rank12,UUsml,VVsml,Singularsml,SVD_tolerance,rank,flop=flop)
	if(present(Flops))Flops=Flops+flop

	! write(111,*)UUsml(1:rank1,1:rank)
	! stop

	! call gemm_omp(QQ1(1:rankmax_r,1:rank1),UUsml(1:rank1,1:rank),matU(1:rankmax_r,1:rank),rankmax_r,rank,rank1)
	call gemmf90(QQ1,rankmax_r,UUsml,rank1,matU,rankmax_r,'N','N',rankmax_r,rank,rank1,cone,czero,flop=flop)
	if(present(Flops))Flops=Flops+flop

	! call gemmNH_omp(VVsml(1:rank,1:rank2),QQ2(1:rankmax_c,1:rank2),matV(1:rank,1:rankmax_c),rank,rankmax_c,rank2)
	call gemmf90(VVsml,rank12,QQ2,rankmax_c,matV,rmax,'N','H',rank,rankmax_c,rank2, cone,czero,flop=flop)
	if(present(Flops))Flops=Flops+flop

	Singular(1:rank) = Singularsml(1:rank)
	deallocate(UUsml,VVsml,Singularsml,matM)


	deallocate(matZRcol1)
	deallocate(matZcRrow1)
	deallocate(tau)
	deallocate(jpvt)
	deallocate(QQ1)
	deallocate(RR1)
	deallocate(QQ2)
	deallocate(RR2)
	deallocate(RrowcZRcol)



	! allocate(matM(rankmax_r,rankmax_c))
	! if(rankmax_r==rmax)then
	! call GetRank(rmax,rmax,matRrow,rank,tolerance)
	! write(*,*)rmax,rank,'ga'
		! call GeneralInverse(rmax,rmax,matRrow,matInv,tolerance)
		! allocate(matrixtemp(rmax,rmax))
		! matrixtemp = matInv
		! do ii=1,rmax
		! do jj=1,rmax
			! matInv(ii,jj) = conjg(cmplx(matrixtemp(jj,ii),kind=8))
		! end do
		! end do
		! deallocate(matrixtemp)
		! call gemmNH_omp(matInv,matZcRrow,matM,rmax,rankmax_c,rmax)
	! else if(rankmax_c==rmax)then
	! call GetRank(rmax,rmax,matRcol,rank,tolerance)
	! write(*,*)rmax,rank,'ga'
		! call GeneralInverse(rmax,rmax,matRcol,matInv,tolerance)
		! call gemm_omp(matZRcol,matInv,matM,rankmax_r,rankmax_c,rmax)
	! end if

	! write(*,*)fnorm(matM,rankmax_r,rankmax_c),'woao'

	! rank12 = min(rankmax_r,rankmax_c)
	! allocate(UUsml(rankmax_r,rank12),VVsml(rank12,rankmax_c),Singularsml(rank12))
	! call SVD_Truncate(matM,rankmax_r,rankmax_c,rank12,UUsml,VVsml,Singularsml,SVD_tolerance,rank)
	! matU(1:rankmax_r,1:rank) = UUsml(1:rankmax_r,1:rank)
	! matV(1:rank,1:rankmax_c) = VVsml(1:rank,1:rankmax_c)
	! Singular(1:rank) = Singularsml(1:rank)
	! deallocate(UUsml,VVsml,Singularsml,matM)

    return

end subroutine RandomizedSVD




subroutine RandomSubMat(ms,me,ns,ne,k,A,Oflag)
	implicit none
	integer ms,me,ns,ne,k,m,n
	DT:: A(:,:)
	DT,allocatable::matrix_small(:,:)
	integer:: Oflag

	m=me-ms+1
	n=ne-ns+1
	call assert(m>0,'m<=0 in RandomSubMat')
	call assert(n>0,'n<=0 in RandomSubMat')

	allocate(matrix_small(m,n))
	call RandomMat(m,n,k,matrix_small,Oflag)
	A(ms:me,ns:ne)=matrix_small
	deallocate(matrix_small)

end subroutine RandomSubMat


subroutine RandomMat(m,n,k,A,Oflag)


	!
	!
#ifdef USEVSL
	use mkl_vsl_type
	use mkl_vsl
#endif

	implicit none

	integer m,n,k,mn_min,ktmp
	class(*):: A(:,:)
	real(kind=8):: c
	real(kind=8),allocatable::Singular(:),Ar(:,:),Ai(:,:)
	complex(kind=8):: ctemp
	complex(kind=8),allocatable:: UU(:,:),VV(:,:)
	integer ii,jj,kk,i,j,flag0,rank
	integer:: Oflag
#ifdef USEVSL
	type(VSL_STREAM_STATE) ::stream(2)
#endif
	integer brng, method, seedd,ierror




 	ktmp = k
	if(ktmp>min(m,n))then
		! write(*,*)'k is not properly set in RandomMat'
		ktmp=min(m,n)
	end if
	! call assert(m<=n,'m>n')

	call assert(k<=min(m,n),'k too large in RandomMat')


#ifdef USEVSL

	select type(A)
	type is (complex(kind=8))

		allocate(Ar(m,n))
		allocate(Ai(m,n))
		brng   = VSL_BRNG_MCG31
		method = VSL_RNG_METHOD_UNIFORM_STD
		call random_number(c)
		seedd = NINT(1000*c)
		ierror=vslnewstream( stream(1), brng,  seedd )
		call random_number(c)
		seedd = NINT(1000*c)
		ierror=vslnewstream( stream(2), brng,  seedd )
		ierror=vdrnguniform( method, stream(1), M*N, Ar, -1d0, 1d0)
		ierror=vdrnguniform( method, stream(2), M*N, Ai, -1d0, 1d0)
		! ierror=vdrnggaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream, M*N, mati, 1d0, 1d0)
		ierror=vsldeletestream(stream(1))
		ierror=vsldeletestream(stream(2))
		A = Ar + junit*Ai
		deallocate(Ar)
		deallocate(Ai)
	type is (real(kind=8))

		allocate(Ar(m,n))
		brng   = VSL_BRNG_MCG31
		method = VSL_RNG_METHOD_UNIFORM_STD
		call random_number(c)
		seedd = NINT(1000*c)
		ierror=vslnewstream( stream(1), brng,  seedd )
		call random_number(c)
		seedd = NINT(1000*c)
		ierror=vdrnguniform( method, stream(1), M*N, Ar, -1d0, 1d0)
		ierror=vsldeletestream(stream(1))
		A = Ar
		deallocate(Ar)

	end select
#else
	do ii=1,m
	do jj=1,n
		call random_dp_number(A(ii,jj))
	end do
	end do


! select type(A)
	! type is (complex(kind=8))

	! mn_min = min(m,n)
	! allocate(Singular(mn_min))
	! allocate(UU(m,mn_min))
	! allocate(VV(mn_min,n))

	! call gesvd_robust(A,Singular,UU,VV,m,n,mn_min)
	! Singular = Singular + 1d0  ! this is to make sure the singular values are not the same
	! ! Singular = 1d0  ! this is to make sure the singular values are not the same



	! ! if(Oflag==1)then
		! ! call assert(mn_min==n,'mn_min not correct')
		! ! A = UU(1:m,1:mn_min)
	! ! else
		! ! !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
		! do jj=1, n
			! do ii=1, m
				! ctemp=(0d0,0d0)
				! do kk=1, k
					! ctemp=ctemp+UU(ii,kk)*VV(kk,jj)*Singular(kk)
				! enddo
				! A(ii,jj)=ctemp
			! enddo
		! enddo
		! ! !$omp end parallel do
	! ! end if

	! deallocate(Singular)
	! deallocate(UU)
	! deallocate(VV)
! end select

#endif


end subroutine RandomMat



subroutine ID_Selection(Mat,select_column,select_row,m,n,rank,tolerance)

    use BPACK_DEFS
    implicit none

    integer i, j, ii, jj, indx, rank_1, rank_2
    integer blocks, index_j, group_nn, rank_blocks
    integer group_m,group_n,size_of_groupm,size_of_groupn
    real(kind=8)tolerance
    integer edgefine_m,edgefine_n, level_butterfly, level_blocks
    integer edge_m, edge_n, header_m, header_n
    integer rank, row, column, rankmax,n,m,rank_c,rank_r
    DT value_Z,value_UV,maxvalue,Mat(m,n)
	DT,allocatable::Mat1(:,:),Mat1T(:,:)
    DT inner_U,inner_V,ctemp
    real(kind=8) inner_UV,flop
    integer select_column(n), select_row(m)

    DT,allocatable:: row_R(:),column_R(:),matrixtemp_V(:,:),matrixtemp_U(:,:),tau(:)
    real(kind=8),allocatable:: norm_row_R(:),norm_column_R(:)
	integer,allocatable :: jpvt(:)

	select_column = 0
	select_row = 0

	allocate(Mat1(m,n))
	Mat1 = Mat
	allocate(Mat1T(n,m))
	call copymatT(Mat,Mat1T,m,n)

	allocate(jpvt(max(m,n)))
	allocate(tau(max(m,n)))

	jpvt=0
	call geqp3modf90(Mat1T,jpvt,tau,tolerance,SafeUnderflow,rank_r)
	select_row(1:rank_r) = jpvt(1:rank_r)
	if(rank_r==0)then
		rank_r=1
		select_row(1)=1
	endif


	jpvt=0
	call geqp3modf90(Mat1,jpvt,tau,tolerance,SafeUnderflow,rank_c)
	select_column(1:rank_c) = jpvt(1:rank_c)
	if(rank_c==0)then
		rank_c=1
		select_column(1)=1
	endif

	rank = min(rank_c,rank_r)


	deallocate(jpvt)
	deallocate(tau)
	deallocate(Mat1)
	deallocate(Mat1T)

    return

end subroutine ID_Selection





subroutine ACA_CompressionFull(mat,matU,matV,rankmax_r,rankmax_c,rmax,rank,tolerance,SVD_tolerance)


    use BPACK_DEFS
    implicit none

    integer i, j, ii, jj, indx, rank_1, rank_2
    integer blocks, index_j, group_nn, rank_blocks
    integer group_m,group_n,size_of_groupm,size_of_groupn
    real(kind=8) norm_Z,norm_U,norm_V,tolerance,SVD_tolerance
    integer edgefine_m,edgefine_n, level_butterfly, level_blocks
    integer edge_m, edge_n, header_m, header_n
    integer rank, ranknew, row, column, rankmax,rankmax_c,rankmax_r,rankmax_min,rmax
    DT value_Z,value_UV,maxvalue
    DT inner_U,inner_V,ctemp
    real(kind=8) inner_UV
    integer,allocatable:: select_column(:), select_row(:)
	DT::mat(rankmax_r,rankmax_c),matU(rankmax_r,rmax),matV(rmax,rankmax_c)

    DT,allocatable:: row_R(:),column_R(:)
    real(kind=8),allocatable:: norm_row_R(:),norm_column_R(:)

	DT, allocatable :: QQ1(:,:), RR1(:,:),QQ2(:,:), RR2(:,:),QQ2tmp(:,:), RR2tmp(:,:), UUsml(:,:), VVsml(:,:),tau_Q(:),mattemp(:,:),matU1(:,:),matV1(:,:)
	real(kind=8), allocatable :: Singularsml(:)


	rankmax_min = min(rankmax_r,rankmax_c)
    norm_Z=0
	select_column = 0
	select_row = 0

    allocate(row_R(rankmax_c),column_R(rankmax_r))
    allocate(norm_row_R(rankmax_c),norm_column_R(rankmax_r))

    select_row(1)=1

    ! !$omp parallel do default(shared) private(j,value_Z)
    do j=1, rankmax_c
        value_Z=mat(select_row(1),j)
        row_R(j)=value_Z
        norm_row_R(j)=dble(value_Z*conjg(cmplx(value_Z,kind=8)))
    enddo
    ! !$omp end parallel do

    select_column(1)=maxloc(norm_row_R,1)
    maxvalue=row_R(select_column(1))

    ! !$omp parallel do default(shared) private(j)
    do j=1, rankmax_c
        row_R(j)=row_R(j)/maxvalue
    enddo
    ! !$omp end parallel do
    ! !$omp parallel do default(shared) private(j)
    do j=1, rankmax_c
        matV(1,j)=row_R(j)
    enddo
    ! !$omp end parallel do

    ! !$omp parallel do default(shared) private(i,value_Z)
    do i=1,rankmax_r
        value_Z=mat(i,select_column(1))
        column_R(i)=value_Z
        norm_column_R(i)=dble(value_Z*conjg(cmplx(value_Z,kind=8)))
    enddo
    ! !$omp end parallel do

    norm_column_R(select_row(1))=0

    ! !$omp parallel do default(shared) private(i)
    do i=1,rankmax_r
        matU(i,1)=column_R(i)
    enddo
    ! !$omp end parallel do

    norm_U=norm_vector(column_R,rankmax_r)
    norm_V=norm_vector(row_R,rankmax_c)
    norm_Z=norm_Z+norm_U*norm_V

	! if(rankmax<2)write(*,*)'rankmax'
    select_row(2)=maxloc(norm_column_R,1)

    rank=1
	! write(*,*)column_R,row_R
	! write(*,*)norm_Z*ACA_tolerance_forward**2,norm_U*norm_V,'hehe',ACA_tolerance_forward
    do while (norm_Z*tolerance**2<norm_U*norm_V .and. rank<rankmax_min)

        ! !$omp parallel do default(shared) private(j,i,value_Z,value_UV)
        do j=1,rankmax_c
            value_Z=mat(select_row(rank+1),j)
            value_UV=0
            do i=1,rank
                value_UV=value_UV+matU(select_row(rank+1),i)*matV(i,j)
            enddo
            row_R(j)=value_Z-value_UV
            norm_row_R(j)=dble(row_R(j)*conjg(cmplx(row_R(j),kind=8)))
        enddo
        ! !$omp end parallel do

        do i=1,rank
        	norm_row_R(select_column(i))=0
        enddo

        select_column(rank+1)=maxloc(norm_row_R,1)
        maxvalue=row_R(select_column(rank+1))

        ! !$omp parallel do default(shared) private(j)
        do j=1,rankmax_c
    	row_R(j)=row_R(j)/maxvalue
        enddo
        ! !$omp end parallel do
        ! !$omp parallel do default(shared) private(j)
        do j=1,rankmax_c
            matV(rank+1,j)=row_R(j)
        enddo
        ! !$omp end parallel do

        ! !$omp parallel do default(shared) private(i,j,value_Z,value_UV)
        do i=1,rankmax_r
            value_Z=mat(i,select_column(rank+1))
            value_UV=0
            do j=1,rank
                value_UV=value_UV+matU(i,j)*matV(j,select_column(rank+1))
            enddo
            column_R(i)=value_Z-value_UV
            norm_column_R(i)=dble(column_R(i)*conjg(cmplx(column_R(i),kind=8)))
        enddo
        ! !$omp end parallel do

        do i=1,rank+1
            norm_column_R(select_row(i))=0
        enddo

        ! !$omp parallel do default(shared) private(i)
        do i=1,rankmax_r
            matU(i,rank+1)=column_R(i)
        enddo
        ! !$omp end parallel do

        norm_U=norm_vector(column_R,rankmax_r)
        norm_V=norm_vector(row_R,rankmax_c)
        inner_UV=0
        do j=1,rank
            inner_U=0
            inner_V=0
	        ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_U)
            do i=1,rankmax_r
                ctemp=matU(i,rank+1)*matU(i,j)
                inner_U=inner_U+ctemp*conjg(cmplx(ctemp,kind=8))
            enddo
            ! !$omp end parallel do
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_V)
            do i=1,rankmax_c
                ctemp=matV(rank+1,i)*matV(j,i)
                inner_V=inner_V+ctemp*conjg(cmplx(ctemp,kind=8))
            enddo
            ! !$omp end parallel do
            inner_UV=inner_UV+dble(2*sqrt(inner_U*inner_V))
        enddo
        norm_Z=norm_Z+inner_UV+norm_U*norm_V

        rank=rank+1
		if(rank>rmax)then
			write(*,*)'increase rmax',rank,rmax
			stop
		end if
        if (rank<rankmax_min) then
            select_row(rank+1)=maxloc(norm_column_R,1)
        endif

    enddo

	! write(*,*)select_row(1:rank),select_column(1:rank)

    deallocate(row_R,column_R)
    deallocate(norm_row_R,norm_column_R)


! ACA followed by SVD

	allocate(QQ1(rankmax_r,rank))
	! call copymatN(matU(1:rankmax_r,1:rank),QQ1,rankmax_r,rank)
	QQ1 = matU(1:rankmax_r,1:rank)
	allocate (tau_Q(rank))
	call geqrff90(QQ1,tau_Q)

	allocate (RR1(rank,rank))
	RR1=0d0
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rank
		do i=1, j
			RR1(i,j)=QQ1(i,j)
		enddo
	enddo
	! !$omp end parallel do
	call un_or_gqrf90(QQ1,tau_Q,rankmax_r,rank,rank)
	deallocate(tau_Q)


	allocate(QQ2tmp(rankmax_c,rank))
	call copymatT(matV(1:rank,1:rankmax_c),QQ2tmp,rank,rankmax_c)
	allocate (tau_Q(rank))
	call geqrff90(QQ2tmp,tau_Q)

	allocate (RR2tmp(rank,rank))
	RR2tmp=0d0
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rank
		do i=1, j
			RR2tmp(i,j)=QQ2tmp(i,j)
		enddo
	enddo
	! !$omp end parallel do
	call un_or_gqrf90(QQ2tmp,tau_Q,rankmax_c,rank,rank)
	deallocate(tau_Q)

	allocate(QQ2(rank,rankmax_c))
	call copymatT(QQ2tmp,QQ2,rankmax_c,rank)
	allocate(RR2(rank,rank))
	call copymatT(RR2tmp,RR2,rank,rank)



	! allocate(matU1(rankmax_r,rank))
	! allocate(matV1(rank,rankmax_c))
	! call gemm_omp(QQ1,RR1,matU1,rankmax_r,rank,rank)

	! call gemm_omp(RR2,QQ2,matV1,rank,rankmax_c,rank)

	! write(*,*)fnorm(matU1-matU(1:rankmax_r,1:rank),rankmax_r,rank),fnorm(matV1-matV(1:rank,1:rankmax_c),rank,rankmax_c)



	deallocate(QQ2tmp,RR2tmp)
	allocate(mattemp(rank,rank))
	! call gemm_omp(RR1,RR2,mattemp,rank,rank,rank)
	call gemmf90(RR1,rank,RR2,rank,mattemp,rank,'N','N',rank,rank,rank,cone,czero)



	allocate(UUsml(rank,rank),VVsml(rank,rank),Singularsml(rank))
	call SVD_Truncate(mattemp,rank,rank,rank,UUsml,VVsml,Singularsml,SVD_tolerance,ranknew)

	! call gemm_omp(QQ1,UUsml(1:rank,1:ranknew),matU(1:rankmax_r,1:ranknew),rankmax_r,ranknew,rank)
	call gemmf90(QQ1,rankmax_r,UUsml,rank,matU,rankmax_r,'N','N',rankmax_r,ranknew,rank,cone,czero)
	! call gemm_omp(VVsml(1:ranknew,1:rank),QQ2,matV(1:ranknew,1:rankmax_c),ranknew,rankmax_c,rank)
	call gemmf90(VVsml,rank,QQ2,rank,matV,rmax,'N','N',ranknew,rankmax_c,rank,cone,czero)
	! write(*,*)'aca rank:',rank,'after svd',ranknew

	rank = ranknew
	do i=1,ranknew
		matU(1:rankmax_r,i) = matU(1:rankmax_r,i)*Singularsml(i)
	end do
	deallocate(mattemp,RR1,RR2,QQ1,QQ2,UUsml,VVsml,Singularsml)

	deallocate(select_column)
	deallocate(select_row)

    return

end subroutine ACA_CompressionFull


recursive subroutine RRQR_LQ(mat,mm,nn,mn,UU,VV,tolerance,rank,lr,flops)
!
!
implicit none
integer mm,nn,mn,rank,ii,jj,rank_new
real(kind=8):: tolerance
DT::mat(mm,nn),UU(mm,mn),VV(mn,nn)
DT,allocatable::mat0(:,:),UUt(:,:),VVt(:,:),tau(:),RR1(:,:)
real(kind=8):: Singular(mn)
integer::i,j,flag
real(kind=8),optional::flops
real(kind=8)::flop
character::lr
integer, allocatable :: jpiv(:)

if(lr=='L')then ! LQ
	allocate(mat0(nn,mm))
	allocate(UUt(nn,mn))
	allocate(VVt(mn,mm))
	call copymatT(mat,mat0,mm,nn)
	call RRQR_LQ(mat0,nn,mm,mn,UUt,VVt,tolerance,rank,'R',flops)
	call copymatT(UUt,VV,nn,mn)
	call copymatT(VVt,UU,mn,mm)
	deallocate(mat0)
	deallocate(UUt)
	deallocate(VVt)
elseif(lr=='R')then ! QR
	if(present(flops))flops=0
	allocate(mat0(mm,nn))
	mat0 = mat
	allocate(jpiv(nn))
	jpiv=0
	allocate(tau(mn))
	tau=0
	call geqp3modf90(mat0,jpiv,tau,tolerance,SafeUnderflow,rank,flop=flop)
	if(present(flops))flops = flops + flop
	if(rank>0)then
		allocate(RR1(1:rank,nn))
		RR1=0
		do j=1, nn
			do i=1, min(j,rank)
				RR1(i,j)=mat0(i,j)
			enddo
		enddo

		call un_or_gqrf90(mat0,tau,mm,rank,rank,flop=flop)
		if(present(flops))flops = flops + flop
		UU(1:mm,1:rank)=mat0(1:mm,1:rank)
		do j=1,nn
			VV(1:rank,jpiv(j))=RR1(1:rank,j)
		enddo
		deallocate(RR1)
	else
		rank=1
		UU=0
		VV=0
	endif

	deallocate(mat0)
	deallocate(jpiv)
	deallocate(tau)
endif


end subroutine RRQR_LQ


subroutine SVD_Truncate(mat,mm,nn,mn,UU,VV,Singular,tolerance,rank,flop)
!
!
implicit none
integer mm,nn,mn,rank,ii,jj
real(kind=8):: tolerance
DT::mat(mm,nn),UU(mm,mn),VV(mn,nn)
DT,allocatable::mat0(:,:)
real(kind=8):: Singular(mn)
integer::i,flag
real(kind=8),optional::flop

UU=0
VV=0
Singular=0

allocate(mat0(mm,nn))
mat0 = mat
call gesvd_robust(mat0,Singular,UU,VV,mm,nn,mn,flop=flop)
rank=mn

if(Singular(1)<SafeUnderflow)then
	rank = 1
	UU=0
	VV=0
	Singular=0
else
	rank = mn
	do i=1,mn
		if (Singular(i)/Singular(1)<=tolerance) then
			rank=i
			if(Singular(i)<Singular(1)*tolerance/10)rank = i -1
			exit
		end if
	end do
endif
deallocate(mat0)



end subroutine SVD_Truncate




subroutine PSVD_Truncate(mm,nn,mat,descMat,UU,VV,descUU,descVV,Singular,tolerance,rank,ctxt,flop)
implicit none
integer mm,nn,mnmin,rank,ii,jj
real(kind=8):: tolerance
DT::mat(:,:),UU(:,:),VV(:,:)
DT,allocatable::mat0(:,:)
real(kind=8):: Singular(:)
integer::i,flag
real(kind=8),optional::flop
integer::descMat(9),descUU(9),descVV(9)
integer::ctxt,nprow, npcol, myrow, mycol,iproc,myi

call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
mnmin = min(mm,nn)

call pgesvdf90('V', 'V', mm, nn, mat, 1, 1, descMat, Singular, UU, 1, 1, descUU, VV, 1, 1, descVV,flop=flop)

rank = mnmin
if(Singular(1)>SafeUnderflow)then
	rank = mnmin
	do i=1,mnmin
		if (Singular(i)/Singular(1)<=tolerance) then
			rank=i
			if(Singular(i)<Singular(1)*tolerance/10)rank = i -1
			exit
		end if
	end do
else
	rank=1
	Singular(1)=0
	VV=0
	UU=0
endif

end subroutine PSVD_Truncate


! sort the first dimension of the array
subroutine PIKSRT_DBLE_Multi(N,M,ARR)
  implicit none
  integer j,i,N,M
  real(kind=8) ARR(N,M)
  real(kind=8),allocatable::a(:)
  allocate(a(M))
  do j=2, N
    a=ARR(j,:)
    do i=j-1,1,-1
      if (ARR(i,1)<=a(1)) goto 10
      ARR(i+1,:)=ARR(i,:)
    end do
	i=0
10  ARR(i+1,:)=a
  end do
  deallocate(a)
  return
end subroutine PIKSRT_DBLE_Multi


! sort the first dimension of the array
subroutine PIKSRT_INT_Multi(N,M,ARR)
  implicit none
  integer j,i,N,M
  integer ARR(N,M)
  integer,allocatable::a(:)
  allocate(a(M))
  do j=2, N
    a=ARR(j,:)
    do i=j-1,1,-1
      if (ARR(i,1)<=a(1)) goto 11
      ARR(i+1,:)=ARR(i,:)
    end do
	i=0
11  ARR(i+1,:)=a
  end do
  deallocate(a)
  return
end subroutine PIKSRT_INT_Multi


!*** remove the duplicates in an integer array
subroutine remove_dup_int(array,nin,nout)
	implicit none
	integer array(nin)
	real(kind=8) array1(nin)
	integer order(nin)
	integer nin,nout,last,ii

	array1(1:nin)=array(1:nin)
	call quick_sort(array1,order,nin)

	nout=0
	last=0
	do ii=1,nin
		if(NINT(array1(ii))/=last)then
			nout=nout+1
			last = NINT(array1(ii))
			array(nout)=last
		endif
	enddo
end subroutine remove_dup_int


subroutine quick_sort(list, order, n)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.

IMPLICIT NONE
integer n
real(kind=8) list(n)
INTEGER  order(n)

! Local variable
INTEGER :: i

DO i = 1, n
  order(i) = i
END DO

CALL quick_sort_1(1, n)

CONTAINS

RECURSIVE subroutine quick_sort_1(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
real(kind=8)                :: reference, temp
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end)

ELSE
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1

  DO
    ! Scan list from left end until element >= reference is found
    DO
      i = i + 1
      IF (list(i) >= reference) EXIT
    END DO
    ! Scan list from right end until element <= reference is found
    DO
      j = j - 1
      IF (list(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL quick_sort_1(left_end, j)
  IF (i < right_end) CALL quick_sort_1(i, right_end)
END IF

END subroutine quick_sort_1


subroutine interchange_sort(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
real(kind=8)                :: temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i); list(i) = list(j); list(j) = temp
      itemp = order(i); order(i) = order(j); order(j) = itemp
    END IF
  END DO
END DO

END subroutine interchange_sort

END subroutine quick_sort


subroutine Cart2Sph(xin,yin,zin,origin,r,theta,phi)
implicit none
real(kind=8),intent(in)::xin,yin,zin,origin(3)
real(kind=8),intent(out)::r,theta,phi
real(kind=8):: x,y,z

	x = xin-origin(1)
	y = yin-origin(2)
	z = zin-origin(3)

	r = sqrt(x**2 + y**2 + z**2)
	theta = acos(z/r)
	if(r==0)theta=0
	if(x==0 .and. y>=0)then
		phi = pi/2
	else if(x==0 .and. y<0)then
		phi = 3*pi/2

	else
		if(y==0 .and. x>=0)then
			phi = 0
		else if(y==0 .and. x<0)then
			phi = pi
		else

			phi = atan(y/x)
			if(phi>0 .and. x<0)phi = phi+pi
			if(phi<0 .and. x>0)phi = phi + 2*pi
			if(phi<0 .and. x<0)phi = phi + pi
		end if
	end if
end subroutine Cart2Sph





complex(kind=8) function Hankel02_Func(x)

use BPACK_DEFS
implicit none

    real(kind=8) x
    complex(kind=8) y

    Hankel02_Func=BesselJ0_func(x)-junit*BesselY0_func(x)

    return

end function Hankel02_Func

real(kind=8) function BesselJ0_func(x)

    implicit none

    real(kind=8) x, z, ax
    real(kind=8) y, rtemp1, rtemp2, xx

    ax=abs(x)
    if (ax<8d0) then
        y=x*x
        rtemp1=57568490574.0d0+y*(-13362690354.0d0+y*(651619640.7d0+y*(-11214424.18d0+y*(77392.33017d0+y*(-184.9052456d0)))))
        rtemp2=57568490411.0d0+y*(1029532985.0d0+y*(9494680.718d0+y*(59272.64853d0+y*(267.8532712d0+y*1.0d0))))
        BesselJ0_func=rtemp1/rtemp2
    else
        z=8.0d0/ax
        y=z*z

        xx=ax-0.785398164d0
        rtemp1=1.0d0+y*(-0.1098628627d-2+y*(0.2734510407d-4+y*(-0.2073370639d-5+y*0.2093887211d-6)))
        rtemp2=-0.1562499995d-1+y*(0.1430488765d-3+y*(-0.6911147651d-5+y*(0.7621095161d-6-y*0.934935152d-7)))
        BesselJ0_func=sqrt(0.636619772d0/ax)*(cos(xx)*rtemp1-z*sin(xx)*rtemp2)
    endif

    return

end function BesselJ0_func

real(kind=8) function BesselY0_func(x)

    implicit none

    real(kind=8) x, z, ax
    real(kind=8) y, rtemp1, rtemp2, xx

    if (x<8.0d0) then
        y=x*x
        rtemp1=-2957821389.0d0+y*(7062834065.0d0+y*(-512359803.6d0+y*(10879881.29d0+y*(-86327.92757d0+y*228.4622733d0))))
        rtemp2=40076544269.0d0+y*(745249964.8d0+y*(7189466.438d0+y*(47447.26470d0+y*(226.1030244d0+y*1.0d0))))
        BesselY0_func=(rtemp1/rtemp2)+0.636619772d0*BesselJ0_func(x)*LOG(x)
    else
        z=8.0d0/x
        y=z*z

        xx=x-0.785398164d0
        rtemp1=1.0d0+y*(-0.1098628627d-2+y*(0.2734510407d-4+y*(-0.2073370639d-5+y*0.2093887211d-6)))
        rtemp2=-0.1562499995d-1+y*(0.1430488765d-3+y*(-0.6911147651d-5+y*(0.7621095161d-6-y*0.934935152d-7)))
        BesselY0_func=sqrt(0.636619772d0/x)*(sin(xx)*rtemp1+z*cos(xx)*rtemp2)
    endif

    return

end function BesselY0_func

!**** create a array treeleaf holding size of each leaf box of a tree with nlevel levels (0<=level<=nlevel)
recursive subroutine CreateLeaf_Natural(nlevel,level,group,idxs,idxe,treeleaf)
implicit none
integer nlevel,level,group,idxs,idxe
integer treeleaf(2**nlevel)

	if(level==nlevel)then
		treeleaf(group-2**nlevel+1)=idxe-idxs+1
	else
		call CreateLeaf_Natural(nlevel,level+1,2*group,idxs,int((idxs+idxe)/2d0),treeleaf)
		call CreateLeaf_Natural(nlevel,level+1,2*group+1,int((idxs+idxe)/2d0)+1,idxe,treeleaf)
	endif
end subroutine CreateLeaf_Natural


subroutine NumberingPtree(ptree)
	implicit none
	type(proctree)::ptree
	integer :: level,group

	ptree%pgrp(1)%head=0 ; ptree%pgrp(1)%tail=ptree%nproc-1; ptree%pgrp(1)%nproc=ptree%nproc
	do level=0, ptree%nlevel-1
		do group=2**level, 2**(level+1)-1
			if (level<ptree%nlevel-1) then
				if (ptree%pgrp(group)%nproc==1) then
					ptree%pgrp(2*group)%head=ptree%pgrp(group)%head
					ptree%pgrp(2*group)%tail=ptree%pgrp(group)%tail
					ptree%pgrp(2*group)%nproc=ptree%pgrp(group)%nproc
					ptree%pgrp(2*group+1)%head=ptree%pgrp(group)%head
					ptree%pgrp(2*group+1)%tail=ptree%pgrp(group)%tail
					ptree%pgrp(2*group+1)%nproc=ptree%pgrp(group)%nproc
				else
					ptree%pgrp(2*group)%head=ptree%pgrp(group)%head
					ptree%pgrp(2*group)%tail=int((ptree%pgrp(group)%head+ptree%pgrp(group)%tail)/2)
					ptree%pgrp(2*group)%nproc=ptree%pgrp(2*group)%tail-ptree%pgrp(2*group)%head+1

					ptree%pgrp(2*group+1)%head=ptree%pgrp(2*group)%tail+1
					ptree%pgrp(2*group+1)%tail=ptree%pgrp(group)%tail
					ptree%pgrp(2*group+1)%nproc=ptree%pgrp(2*group+1)%tail-ptree%pgrp(2*group+1)%head+1
				endif
			end if
		end do
	end do
end subroutine NumberingPtree


!**** computation of the local butterfly block ranges owned by this MPI rank
	!ptree: process tree
	!pgno: the process group number that shares this butterfly
	!level_butterfly: number of butterfly levels
	!level: the specified level that requires the computation 0<=level<=level_butterfly+1
	!idx_r: starting index of local row blocks
	!inc_r: increments of local row blocks
	!nr: number of local row blocks
	!idx_c: starting index of local column blocks
	!inc_c: increments of local column blocks
	!nc: number of local column blocks
	!dir: 'R': row-wise ordering (2^level row blocks and 2^(level_butterfly-level) column blocks) or 'C': column-wise ordering (2^(level-1) row blocks and 2^(level_butterfly-level+1) column blocks)
subroutine GetLocalBlockRange(ptree,pgno,level,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,dir)
	implicit none
	type(proctree)::ptree
	integer :: level_p,level,ll,group,level_butterfly,idx_r,inc_r,nr,idx_c,inc_c,nc,pgno,pgno1,found,nleaf,ith
	character::dir
	if(IOwnPgrp(ptree,pgno))then
		found=0
		level_p=0
		pgno1=pgno
		ith=1
		do while(found==0)
			if(ptree%pgrp(pgno1)%nproc==1)then
				found=1
				exit
			endif
			if(IOwnPgrp(ptree,pgno1*2))then
				pgno1=pgno1*2
				ith=ith*2
			elseif(IOwnPgrp(ptree,pgno1*2+1))then
				pgno1=pgno1*2+1
				ith=ith*2+1
			endif
			level_p=level_p+1
		enddo
		ith = ith - 2**level_p + 1
		call assert(level_p<=level_butterfly,'too many processes sharing this block')
		nleaf = 2**(level_butterfly-level_p)
		if(dir=='R')then
			idx_c=(ith-1)*nleaf+1
			nc=nleaf
			inc_c=1
			idx_r=1
			nr=1
			inc_r=1
			do ll=1,level
				if(ll/=level_butterfly+1)then
					idx_r = 2*idx_r-mod(idx_c,2) ! odd/even column indices mapped to odd/even row indices
					idx_c = floor((idx_c-1)/2d0)+1
					if(nc>1)then ! if at previous level column blocks are contiguous, double #row and half #column
						nc=nc/2
						nr=nr*2
					else         ! otherwise double row increments
						inc_r=inc_r*2
					endif
				endif
			enddo
		elseif(dir=='C')then
			idx_r=(ith-1)*nleaf+1
			nr=nleaf
			inc_r=1
			idx_c=1
			nc=1
			inc_c=1
			do ll=level_butterfly,level,-1
				if(ll/=0)then
					idx_c = 2*idx_c-mod(idx_r,2) ! odd/even row indices mapped to odd/even column indices
					idx_r = floor((idx_r-1)/2d0)+1
					if(nr>1)then ! if at previous level row blocks are contiguous, double #column and half #row
						nr=nr/2
						nc=nc*2
					else
						inc_c=inc_c*2 ! otherwise double column increments
					endif
				endif
			enddo
		endif
	else
		idx_r=0
		inc_r=0
		nr=0
		idx_c=0
		inc_c=0
		nc=0
	endif
end subroutine GetLocalBlockRange



!**** computation of the local butterfly block ranges owned by this MPI rank
	!ptree: process tree
	!pgno: the process group number that shares this butterfly
	!level_butterfly: number of butterfly levels
	!level: the level where the block resides 0<=level<=level_butterfly+1
	!index_i: row index of the block
	!index_j: column index of the block
	!dir: 'R': row-wise ordering (2^level row blocks and 2^(level_butterfly-level) column blocks) or 'C': column-wise ordering (2^(level-1) row blocks and 2^(level_butterfly-level+1) column blocks)
	!pid: mpi rank that handles the target block
subroutine GetBlockPID(ptree,pgno,level,level_butterfly,index_i,index_j,dir,pid)
	implicit none
	type(proctree)::ptree
	integer :: pid,level_p,num_blocks,level,ll,group,level_butterfly,idx_r,idx_c,idx,index_i,index_j,pgno,pgno1,found,nleaf,ith
	character::dir

	level_p = ptree%nlevel-GetTreelevel(pgno)
	idx_r=index_i
	idx_c=index_j

	if(dir=='R')then
		do ll=level,1,-1
			if(ll/=level_butterfly+1)then
				idx_c = 2*idx_c-mod(idx_r,2) ! odd/even row indices mapped to odd/even column indices
				idx_r = floor_safe((idx_r-1)/2d0)+1
			endif
		enddo
		idx=idx_c
	elseif(dir=='C')then
		do ll=level,level_butterfly
			if(ll/=0)then
				idx_r = 2*idx_r-mod(idx_c,2) ! odd/even column indices mapped to odd/even row indices
				idx_c = floor_safe((idx_c-1)/2d0)+1
			endif
		enddo
		idx=idx_r
	endif
	idx = idx + 2**level_butterfly -1 ! index of this leaf-level group in the corresponding level_butterfly-level tree
	do ll=1,level_butterfly-level_p
		idx = floor_safe(idx/2d0)
	enddo
	idx = idx - 2**level_p +1 ! index of this group's ancestor group in the corresponding level_p-level tree
	pgno1 = pgno*2**level_p ! index of the first process group at pgno's level + level_p of the process tree
	pgno1 = pgno1 + idx -1 ! index of the process group that handles group idx in the process tree



	pid = ptree%pgrp(pgno1)%head
end subroutine GetBlockPID


subroutine CreatePtree(nmpi,groupmembers,MPI_Comm_base,ptree)
	implicit none
	integer nmpi,MPI_Comm_base,MPI_Group_base,MPI_Group_H,MPI_Group_parent,MPI_Comm_parent,groupmembers(nmpi)
	type(proctree)::ptree,ptreecol,ptreerow
	integer :: ierr,Maxgrp,Maxgrpcol,Maxgrprow,MyID_old,level,group,icontxt,ii,jj,kk,pid
	integer :: nprow,npcol,myrow,mycol,nproc,nproc_tmp,nproc_tmp1,nsprow,nsprow1,nspcol,nlevelrow,nlevelcol
	integer,allocatable::pmap(:,:),groupmembers_sml(:),MPI_Group_H_sml(:)

	integer,external :: sys2blacs_handle

	call MPI_Comm_rank(MPI_Comm_base,MyID_old,ierr)
	call MPI_Comm_group(MPI_Comm_base,MPI_Group_base,ierr)
	call MPI_Group_incl(MPI_Group_base, nmpi, groupmembers, MPI_Group_H, ierr)
	call MPI_Comm_Create(MPI_Comm_base,MPI_Group_H,ptree%Comm,ierr)

	if(ptree%Comm/=MPI_COMM_NULL)then
		call MPI_Comm_size(ptree%Comm,ptree%nproc,ierr)
		call MPI_Comm_rank(ptree%Comm,ptree%MyID,ierr)

		call assert(groupmembers(ptree%MyID+1)==MyID_old,'it is assumed the new ID of the Ith proc in groupmembers is I-1')

		ptree%nlevel = ceiling_safe(log(dble(ptree%nproc)) / log(2d0))+1
		Maxgrp=2**(ptree%nlevel)-1
		allocate (ptree%pgrp(Maxgrp))
		call NumberingPtree(ptree)
		allocate (MPI_Group_H_sml(Maxgrp))
		MPI_Group_H_sml = MPI_GROUP_NULL
		do group=1,Maxgrp
			nproc=ptree%pgrp(group)%nproc
			if(group==1)then
				MPI_Comm_parent=ptree%Comm
				MPI_Group_parent=MPI_Group_H
			else
				MPI_Comm_parent=ptree%pgrp(floor_safe(group/2d0))%Comm
				MPI_Group_parent=MPI_Group_H_sml(floor_safe(group/2d0))
			endif

			! ! create the 2D grids as square as possible
			nprow = floor_safe(sqrt(dble(nproc)))
			npcol = floor_safe(nproc/dble(nprow))

			! the following guarantees column dimension is at most one more level than row dimension, this makes parallel ACA implementation easier
			nlevelrow = ceiling_safe(log(dble(nprow)) / log(2d0))+1
			nlevelcol = ceiling_safe(log(dble(npcol)) / log(2d0))+1
			if(nlevelcol>nlevelrow+1)then
				npcol = 2**nprow
			endif

			! ! trail to power of 2 grids, the following two lines can be removed
			! nproc_tmp = 2**floor_safe(log10(dble(nproc))/log10(2d0))
			! nprow  = 2**floor_safe(log10(sqrt(dble(nproc_tmp)))/log10(2d0))
			! npcol = floor_safe(nproc/dble(nprow))


			ptree%pgrp(group)%nprow=nprow
			ptree%pgrp(group)%npcol=npcol
			ptree%pgrp(group)%ctxt=-1
			ptree%pgrp(group)%ctxt1D=-1
			ptree%pgrp(group)%ctxt_head=-1
			ptree%pgrp(group)%Comm=MPI_COMM_NULL


			if(MPI_Comm_parent/=MPI_COMM_NULL)then

				! create the local communicator for this tree node
				allocate(groupmembers_sml(ptree%pgrp(group)%nproc))
				do ii=1,ptree%pgrp(group)%nproc
					groupmembers_sml(ii)=ii-1
				enddo
				if(mod(group,2)==1 .and. group>1)then
				if(ptree%pgrp(floor_safe(group/2d0))%nproc>1)then ! this makes sure the two childs are not the same pgrp as the parent
				groupmembers_sml=groupmembers_sml+ptree%pgrp(group-1)%nproc
				endif
				endif
				call MPI_Group_incl(MPI_Group_parent, ptree%pgrp(group)%nproc, groupmembers_sml, MPI_Group_H_sml(group), ierr)
				call MPI_Comm_Create(MPI_Comm_parent,MPI_Group_H_sml(group),ptree%pgrp(group)%Comm,ierr)
				deallocate(groupmembers_sml)

				! write(*,*)ptree%MyID,group,ptree%pgrp(group)%Comm,ptree%pgrp(group)%Comm==MPI_COMM_NULL,ierr

				! call MPI_barrier(MPI_Comm_parent,ierr)

			endif


			if(ptree%pgrp(group)%Comm/=MPI_COMM_NULL)then

				allocate(pmap(nprow,npcol))
				do jj=1,npcol
				do ii=1,nprow   ! 'row major here'
					kk=npcol*(ii-1)+jj
					pmap(ii,jj)=kk-1
				enddo
				enddo

				! the context involving 2D grids
				ptree%pgrp(group)%ctxt = sys2blacs_handle(ptree%pgrp(group)%Comm)
				call blacs_gridmap( ptree%pgrp(group)%ctxt, pmap, nprow, nprow, npcol )
				deallocate(pmap)

				allocate(pmap(nproc,1))
				do kk=1,nproc
					pmap(kk,1)=kk-1
				enddo

				! the context involving 1D grids non-cyclic
				ptree%pgrp(group)%ctxt1D = sys2blacs_handle(ptree%pgrp(group)%Comm)
				call blacs_gridmap( ptree%pgrp(group)%ctxt1D, pmap, nproc, nproc, 1 )
				! the context involving head proc only
				ptree%pgrp(group)%ctxt_head = sys2blacs_handle(ptree%pgrp(group)%Comm)
				pid=0
				call blacs_gridmap( ptree%pgrp(group)%ctxt_head, pid, 1, 1, 1 )
				deallocate(pmap)
				! call MPI_barrier(ptree%pgrp(group)%Comm,ierr)

			endif

			! create the hierarchical process grids used for parallel ACA
			nsprow = ptree%pgrp(group)%nprow
			nspcol = ptree%pgrp(group)%npcol
			ptreecol%nproc = nspcol
			ptreecol%nlevel = ceiling_safe(log(dble(ptreecol%nproc)) / log(2d0))+1
			Maxgrpcol=2**(ptreecol%nlevel)-1
			allocate (ptreecol%pgrp(Maxgrpcol))
			call NumberingPtree(ptreecol)
			ptreerow%nproc = nsprow
			ptreerow%nlevel = ceiling_safe(log(dble(ptreerow%nproc)) / log(2d0))+1
			Maxgrprow=2**(ptreerow%nlevel)-1
			allocate (ptreerow%pgrp(Maxgrprow))
			call NumberingPtree(ptreerow)

			allocate(ptree%pgrp(group)%gd)
			ptree%pgrp(group)%gd%nsprow=nsprow
			ptree%pgrp(group)%gd%nspcol=nspcol
			ptree%pgrp(group)%gd%hprow=0
			ptree%pgrp(group)%gd%hpcol=0
			ptree%pgrp(group)%gd%gprow=1
			ptree%pgrp(group)%gd%gpcol=1
			call CreateNewGrid(ptree%pgrp(group)%gd,0,ptree,ptreerow,ptreecol,nspcol,MPI_Group_H_sml(group),ptree%pgrp(group)%Comm)

			deallocate(ptreecol%pgrp)
			deallocate(ptreerow%pgrp)

			! call MPI_barrier(ptree%Comm,ierr)
		enddo

		! call MPI_barrier(ptree%Comm,ierr)
		do group=1,Maxgrp
			if(MPI_Group_H_sml(group)/=MPI_GROUP_NULL)call MPI_Group_Free(MPI_Group_H_sml(group),ierr)
		enddo
		deallocate(MPI_Group_H_sml)

		! call MPI_barrier(ptree%Comm,ierr)

		! if(ptree%MyID==Main_ID)then
		! do group=1,Maxgrp
		! write(*,'(A5,I5,A9,I5,A6,I5,A6,I5,A6,I5,A13,I5,A13,I5)')'myid',ptree%MyID,'group no',group,'nproc',ptree%pgrp(group)%nproc,'nprow',ptree%pgrp(group)%nprow,'npcol',ptree%pgrp(group)%npcol,'grid%nsprow',ptree%pgrp(group)%gd%nsprow,'grid%nspcol',ptree%pgrp(group)%gd%nspcol
		! enddo
		! endif

	end if

	! call MPI_barrier(ptree%Comm,ierr)

	call MPI_Group_Free(MPI_Group_base,ierr)
	call MPI_Group_Free(MPI_Group_H,ierr)

	! stop
end subroutine CreatePtree



!******** Create a new square grid gd.
! Note: Ideally, MPI_Comm_parent should involve smaller process counts with increasing recursion level, but I haven't figured out the local numbering in groupmembers_sml and pmap
recursive subroutine CreateNewGrid(gd,cridx,ptree,ptreerow,ptreecol,nspcol_parent,MPI_Group_parent,MPI_Comm_parent)
	implicit none
	type(grid)::gd
	integer cridx,Iown,nproc
	type(proctree)::ptree,ptreerow,ptreecol
	integer nspcol_parent,ii,jj,kk,nsproc
	integer MPI_Group_parent,MPI_Group_H_sml,ierr,MPI_Comm_parent
	integer,allocatable::pmap(:,:),groupmembers_sml(:)
	integer,external :: sys2blacs_handle

	gd%ctxt=-1
	gd%Comm=MPI_COMM_NULL
	MPI_Group_H_sml=MPI_GROUP_NULL
	if(MPI_Comm_parent/=MPI_COMM_NULL)then

		! create the local communicator for this grid
		nsproc = gd%nsprow*gd%nspcol
		Iown=0
		allocate(groupmembers_sml(nsproc))
		do jj=1,gd%nspcol
		do ii=1,gd%nsprow   ! 'row major here'
			kk=nspcol_parent*(gd%hprow+ii-1)+jj+gd%hpcol
			groupmembers_sml(jj+(ii-1)*gd%nspcol)=kk-1
			! if(kk+ptree%pgrp(group)%head-1==31)Iown=1
		enddo
		enddo
		! if(ptree%MyID==0 )write(*,*)'size',nsproc,'cridx',cridx
		call MPI_Comm_size(MPI_Comm_parent,nproc,ierr)
		call MPI_Group_incl(MPI_Group_parent, nsproc, groupmembers_sml, MPI_Group_H_sml, ierr)
		call MPI_Comm_Create(MPI_Comm_parent,MPI_Group_H_sml,gd%Comm,ierr)
		deallocate(groupmembers_sml)
		! call MPI_Group_Free(MPI_Group_H_sml,ierr)
	endif

	if(gd%Comm/=MPI_COMM_NULL)then
		allocate(pmap(gd%nsprow,gd%nspcol))
		do jj=1,gd%nspcol
		do ii=1,gd%nsprow   ! 'row major here'
			kk=gd%nspcol*(ii-1)+jj
			pmap(ii,jj)=kk-1
		enddo
		enddo

		! the context involving 2D grids
		gd%ctxt = sys2blacs_handle(gd%Comm)
		call blacs_gridmap( gd%ctxt, pmap, gd%nsprow, gd%nsprow, gd%nspcol )
		deallocate(pmap)
	endif

	if(cridx<ptreerow%nlevel+ptreecol%nlevel-2)then
		allocate(gd%gdc(2))
		if(mod(cridx+1,2)==1)then
			do ii=1,2
				gd%gdc(ii)%gprow=gd%gprow
				gd%gdc(ii)%hprow=0
				gd%gdc(ii)%nsprow=gd%nsprow
				gd%gdc(ii)%gpcol=gd%gpcol*2+ii-1
				gd%gdc(ii)%hpcol=0
				if(ii==2 .and. ptreecol%pgrp(gd%gpcol)%nproc>1)then
					gd%gdc(ii)%hpcol= ptreecol%pgrp(gd%gpcol*2)%nproc
				endif
				gd%gdc(ii)%nspcol=ptreecol%pgrp(gd%gdc(ii)%gpcol)%nproc
			enddo
		else
			do ii=1,2
				gd%gdc(ii)%gpcol=gd%gpcol
				gd%gdc(ii)%hpcol=0
				gd%gdc(ii)%nspcol=gd%nspcol
				gd%gdc(ii)%gprow=gd%gprow*2+ii-1
				gd%gdc(ii)%hprow=0
				if(ii==2 .and. ptreerow%pgrp(gd%gprow)%nproc>1)then
					gd%gdc(ii)%hprow= ptreerow%pgrp(gd%gprow*2)%nproc
				endif
				gd%gdc(ii)%nsprow=ptreerow%pgrp(gd%gdc(ii)%gprow)%nproc
			enddo
		endif

		do ii=1,2
			call CreateNewGrid(gd%gdc(ii),cridx+1,ptree,ptreerow,ptreecol,gd%nspcol,MPI_Group_H_sml,gd%Comm)
		enddo
	endif
	if(MPI_Group_H_sml/=MPI_GROUP_NULL)call MPI_Group_Free(MPI_Group_H_sml,ierr)

end subroutine CreateNewGrid

! redistribute array 1D block array dat_i distributed among process group pgno_i to 1D block array dat_o distributed among process group pgno_o, M_p_i/M_p_o denote the starting index of each process, head_i/head_o denote the global index of the first element (among all processes) in the dat_i/dat_o
subroutine Redistribute1Dto1D(dat_i,M_p_i,head_i,pgno_i,dat_o,M_p_o,head_o,pgno_o,N,ptree)
implicit none
DT::dat_i(:,:),dat_o(:,:)
integer pgno_i,pgno_o,N
integer M_p_i(:,:),M_p_o(:,:)
integer nproc_i, nproc_o,idxs_i,idxs_o,idxe_i,idxe_o,ii,jj,iii,jjj
type(proctree)::ptree
type(commquant1D),allocatable::sendquant(:),recvquant(:)
integer,allocatable::S_req(:),R_req(:)
integer,allocatable:: statuss(:,:),statusr(:,:)
integer tag,Nreqs,Nreqr,recvid,sendid,ierr,head_i,head_o,sizes,sizer,offs,offr


if(pgno_i==pgno_o .and. ptree%pgrp(pgno_i)%nproc==1)then
	idxs_i=M_p_i(1,1)+head_i
	idxe_i=M_p_i(1,2)+head_i
	idxs_o=M_p_o(1,1) + head_o
	idxe_o=M_p_o(1,2) + head_o
	if(idxs_o<=idxe_i .and. idxe_o>=idxs_i)then
		offs=max(idxs_i,idxs_o) - idxs_i
		sizes = min(idxe_i,idxe_o) - max(idxs_i,idxs_o) + 1
		offr=max(idxs_i,idxs_o) - idxs_o
		sizer = min(idxe_i,idxe_o) - max(idxs_i,idxs_o) + 1
		dat_o(offr+1:offr+sizer,1:N) = dat_i(offs+1:offs+sizes,1:N)
	endif
else

	nproc_i = ptree%pgrp(pgno_i)%nproc
	nproc_o = ptree%pgrp(pgno_o)%nproc
	tag = pgno_o

	allocate(statuss(MPI_status_size,nproc_o))
	allocate(statusr(MPI_status_size,nproc_i))
	allocate(S_req(nproc_o))
	allocate(R_req(nproc_i))


	allocate(sendquant(nproc_o))
	do ii=1,nproc_o
		sendquant(ii)%size=0
	enddo

	allocate(recvquant(nproc_i))
	do ii=1,nproc_i
		recvquant(ii)%size=0
	enddo

	if(IOwnPgrp(ptree,pgno_i))then
		ii = ptree%myid-ptree%pgrp(pgno_i)%head+1
		idxs_i=M_p_i(ii,1)+head_i
		idxe_i=M_p_i(ii,2)+head_i

		do jj=1,nproc_o
			idxs_o=M_p_o(jj,1) + head_o
			idxe_o=M_p_o(jj,2) + head_o
			if(idxs_o<=idxe_i .and. idxe_o>=idxs_i)then
				sendquant(jj)%offset=max(idxs_i,idxs_o) - idxs_i
				sendquant(jj)%size = min(idxe_i,idxe_o) - max(idxs_i,idxs_o) + 1
				if(sendquant(jj)%size>0)then
				allocate(sendquant(jj)%dat(sendquant(jj)%size,N))
				sendquant(jj)%dat = dat_i(sendquant(jj)%offset+1:sendquant(jj)%offset+sendquant(jj)%size,1:N)
				endif
			endif
		enddo
	endif

	if(IOwnPgrp(ptree,pgno_o))then
		jj = ptree%myid-ptree%pgrp(pgno_o)%head+1
		idxs_o=M_p_o(jj,1) + head_o
		idxe_o=M_p_o(jj,2) + head_o

		do ii=1,nproc_i
			idxs_i=M_p_i(ii,1)+head_i
			idxe_i=M_p_i(ii,2)+head_i
			if(idxs_o<=idxe_i .and. idxe_o>=idxs_i)then
				recvquant(ii)%offset=max(idxs_i,idxs_o) - idxs_o
				recvquant(ii)%size = min(idxe_i,idxe_o) - max(idxs_i,idxs_o) + 1
				allocate(recvquant(ii)%dat(recvquant(ii)%size,N))
				recvquant(ii)%dat = 0
			endif
		enddo
	endif


	! post receive
	Nreqr=0
	do ii=1,nproc_i
		if(recvquant(ii)%size>0)then
			jj = ptree%myid-ptree%pgrp(pgno_o)%head+1
			sendid = ii+ptree%pgrp(pgno_i)%head-1
			if(ptree%MyID/=sendid)then
				Nreqr = Nreqr+1
				call MPI_Irecv(recvquant(ii)%dat,recvquant(ii)%size*N,MPI_DT,sendid,tag,ptree%Comm,R_req(Nreqr),ierr)
			endif
		endif
	enddo

	! post send
	Nreqs=0
	do jj=1,nproc_o
		if(sendquant(jj)%size>0)then
			ii = ptree%myid-ptree%pgrp(pgno_i)%head+1
			recvid = jj+ptree%pgrp(pgno_o)%head-1
			if(ptree%MyID==recvid)then
				recvquant(ii)%dat=sendquant(jj)%dat ! make the direct copy if I own both send and receive pieces
			else
				Nreqs = Nreqs+1
				call MPI_Isend(sendquant(jj)%dat,sendquant(jj)%size*N,MPI_DT,recvid,tag,ptree%Comm,S_req(Nreqs),ierr)
			endif
		endif
	enddo

	if(Nreqs>0)then
		call MPI_waitall(Nreqs,S_req,statuss,ierr)
	endif
	if(Nreqr>0)then
		call MPI_waitall(Nreqr,R_req,statusr,ierr)
	endif


	! copy data from receive buffer
	do ii=1,nproc_i
		if(recvquant(ii)%size>0)then
			dat_o(recvquant(ii)%offset+1:recvquant(ii)%offset+recvquant(ii)%size,1:N) = recvquant(ii)%dat
		endif
	enddo

	! deallocation
	deallocate(S_req)
	deallocate(R_req)
	deallocate(statuss)
	deallocate(statusr)
	do jj=1,nproc_o
		if(allocated(sendquant(jj)%dat))deallocate(sendquant(jj)%dat)
	enddo
	deallocate(sendquant)
	do ii=1,nproc_i
		if(allocated(recvquant(ii)%dat))deallocate(recvquant(ii)%dat)
	enddo
	deallocate(recvquant)
endif

end subroutine Redistribute1Dto1D



! redistribute 1D block array dat_i distributed among process group pgno_i to 2D block array dat_o distributed among process group pgno_o
subroutine Redistribute1Dto2D(dat_i,M_p_i,head_i,pgno_i,dat_o,M,head_o,pgno_o,N,ptree)
implicit none
DT::dat_i(:,:),dat_o(:,:)
DT,allocatable::dat_1D(:,:)
integer pgno_i,pgno_o,N,M
integer M_p_i(:,:)
integer nproc_i, nproc_o,idxs_i,idxs_o,idxe_i,idxe_o,ii,jj,iii,jjj
type(proctree)::ptree
integer,allocatable::S_req(:),R_req(:)
integer,allocatable:: statuss(:,:),statusr(:,:)
integer,allocatable:: M_p_1D(:,:)
integer tag,Nreqs,Nreqr,recvid,sendid,ierr,head_i,head_o,sizes,sizer,offs,offr
integer ctxt1D,nproc,nprow,npcol,myrow,mycol,nb1Dc,nb1Dr,myArows,myAcols
integer::desc1D(9),desc2D(9)
integer::ctxt,info

ctxt1D = ptree%pgrp(pgno_o)%ctxt1D
nproc = ptree%pgrp(pgno_o)%nproc
nb1Dc=N
nb1Dr=ceiling_safe(M/dble(nproc))
allocate(M_p_1D(nproc,2))
do ii=1,nproc
	M_p_1D(ii,1) = (ii-1)*nb1Dr+1
	M_p_1D(ii,2) = ii*nb1Dr
enddo
M_p_1D(nproc,2) = M
jj = ptree%myid-ptree%pgrp(pgno_o)%head+1
myArows = M_p_1D(jj,2)-M_p_1D(jj,1)+1
myAcols = N
call descinit( desc1D, M, N, nb1Dr, nb1Dc, 0, 0, ctxt1D, max(myArows,1), info )
! write(*,*)ptree%MyID,M,N,myArows,myAcols,'1D',nproc
if(myArows>0 .and. myAcols>0)then
allocate(dat_1D(myArows,myAcols))
dat_1D=0
endif
ctxt = ptree%pgrp(pgno_o)%ctxt
call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
if(myrow/=-1 .and. mycol/=-1)then
	myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
	myAcols = numroc_wp(N, nbslpk, mycol, 0, npcol)
	call descinit( desc2D, M, N, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
else
	desc2D(2)=-1
endif


call Redistribute1Dto1D(dat_i,M_p_i,head_i,pgno_i,dat_1D,M_p_1D,head_o,pgno_o,N,ptree)

! write(*,*)ptree%MyID,M,N,myArows,myAcols,size(dat_1D,1),size(dat_1D,2),size(dat_o,1),size(dat_o,2),'2D'
call pgemr2df90(M, N, dat_1D, 1, 1, desc1D, dat_o, 1, 1, desc2D, ctxt1D)

deallocate(M_p_1D)
if(allocated(dat_1D))deallocate(dat_1D)

end subroutine Redistribute1Dto2D



! redistribute 2D block array dat_i distributed among process group pgno_i to 1D block array dat_o distributed among process group pgno_o
subroutine Redistribute2Dto1D(dat_i,M,head_i,pgno_i,dat_o,M_p_o,head_o,pgno_o,N,ptree)
implicit none
DT::dat_i(:,:),dat_o(:,:)
DT,allocatable::dat_1D(:,:)
integer pgno_i,pgno_o,N,M
integer M_p_o(:,:)
integer nproc_i, nproc_o,idxs_i,idxs_o,idxe_i,idxe_o,ii,jj,iii,jjj
type(proctree)::ptree
integer,allocatable::S_req(:),R_req(:)
integer,allocatable:: statuss(:,:),statusr(:,:)
integer,allocatable:: M_p_1D(:,:)
integer tag,Nreqs,Nreqr,recvid,sendid,ierr,head_i,head_o,sizes,sizer,offs,offr
integer ctxt1D,nproc,nprow,npcol,myrow,mycol,nb1Dc,nb1Dr,myArows,myAcols
integer::desc1D(9),desc2D(9)
integer::ctxt,info

ctxt1D = ptree%pgrp(pgno_i)%ctxt1D
nproc = ptree%pgrp(pgno_i)%nproc
nb1Dc=N
nb1Dr=ceiling_safe(M/dble(nproc))
allocate(M_p_1D(nproc,2))
do ii=1,nproc
	M_p_1D(ii,1) = (ii-1)*nb1Dr+1
	M_p_1D(ii,2) = ii*nb1Dr
enddo
M_p_1D(nproc,2) = M
jj = ptree%myid-ptree%pgrp(pgno_i)%head+1
myArows = M_p_1D(jj,2)-M_p_1D(jj,1)+1
myAcols = N
call descinit( desc1D, M, N, nb1Dr, nb1Dc, 0, 0, ctxt1D, max(myArows,1), info )
! write(*,*)ptree%MyID,M,N,myArows,myAcols,'1D',nproc
if(myArows>0 .and. myAcols>0)then
allocate(dat_1D(myArows,myAcols))
dat_1D=0
endif

ctxt = ptree%pgrp(pgno_i)%ctxt
call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
if(myrow/=-1 .and. mycol/=-1)then
	myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
	myAcols = numroc_wp(N, nbslpk, mycol, 0, npcol)
	call descinit( desc2D, M, N, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
else
	myArows=0
	myAcols=0
	desc2D(2)=-1
endif

! write(*,*)ptree%MyID,M,N,myArows,myAcols,size(dat_i,1),size(dat_i,2),size(dat_1D,1),size(dat_1D,2),'2D2D',isnanMat(dat_i,size(dat_i,1),size(dat_i,2)),isnanMat(dat_1D,size(dat_1D,1),size(dat_1D,2)),myrow,mycol,pgno_i,ctxt
call pgemr2df90(M, N, dat_i, 1, 1, desc2D, dat_1D, 1, 1, desc1D, ctxt1D)
! write(*,*)ptree%MyID,M,N,myArows,myAcols,size(dat_1D,1),size(dat_1D,2),size(dat_o,1),size(dat_o,2),'1D1D',isnanMat(dat_1D,size(dat_1D,1),size(dat_1D,2)),isnanMat(dat_o,size(dat_o,1),size(dat_o,2)),myrow,mycol
call Redistribute1Dto1D(dat_1D,M_p_1D,head_i,pgno_i,dat_o,M_p_o,head_o,pgno_o,N,ptree)

deallocate(M_p_1D)
if(allocated(dat_1D))deallocate(dat_1D)

end subroutine Redistribute2Dto1D




! get the level (indexed from 1) of a node in a tree. gno is the node number starting from root (1). Note that the process tree levels are indexed from 1, the basis_group levels are indexed from 0
integer function GetTreelevel(gno)
	implicit none
	integer gno,ii,level
	ii=gno
	level=0
	do while(ii/=0)
		ii = INT(ii/2d0)
		level = level+1
	enddo
	GetTreelevel = level
end function GetTreelevel

! check if I share this process group
logical function IOwnPgrp(ptree,pgno)
	implicit none
	integer pgno
	type(proctree)::ptree
	IOwnPgrp = ptree%MyID>=ptree%pgrp(pgno)%head .and. ptree%MyID<=ptree%pgrp(pgno)%tail
end function IOwnPgrp


! convert global index to local index in block-cyclic distribution
subroutine g2l(i,n,np,nb,p,il)
   implicit none
   integer :: i    ! global array index, input
   integer :: n    ! global array dimension, input
   integer :: np   ! processor array dimension, input
   integer :: nb   ! block size, input
   integer :: p    ! processor array index, output
   integer :: il   ! local array index, output
   integer :: im1
   im1 = i-1
   p   = mod((im1/nb),np)
   il  = (im1/(np*nb))*nb + mod(im1,nb) + 1
   return
end subroutine g2l



! convert local index to global index in block-cyclic distribution
subroutine l2g(il,p,n,np,nb,i)
   implicit none
   integer :: il   ! local array index, input
   integer :: p    ! processor array index, input
   integer :: n    ! global array dimension, input
   integer :: np   ! processor array dimension, input
   integer :: nb   ! block size, input
   integer :: i    ! global array index, output
   integer :: ilm1
   ilm1 = il-1
   i    = (((ilm1/nb) * np) + p)*nb + mod(ilm1,nb) + 1
   return
end subroutine l2g


integer function lcm(a,b)
integer:: a,b
	lcm = a*b / gcd(a,b)
end function lcm

integer function gcd(a,b)
integer :: a,b,t,as,bs
	as=a
	bs=b
	do while (bs/=0)
		t = bs
		bs = mod(as,bs)
		as = t
	end do
	gcd = abs(as)
end function gcd


integer function numroc_wp(n, nb, iproc, isrcproc, nprocs)
integer :: n, nb, iproc, isrcproc, nprocs
integer :: numroc ! blacs routine
numroc_wp = numroc(n, nb, iproc, isrcproc, nprocs)
end function numroc_wp

! get my row and column index in a 2D grid that is column major consisting of contiguous proc ids
subroutine Gridinfo_2D(pmap,MyID,myrow,mycol)
implicit none
integer pmap(3)
integer myrow,mycol,nprow,npcol,idstart,idend,MyID
nprow=pmap(1)
npcol=pmap(2)
idstart=pmap(3)
idend=pmap(3)+nprow*npcol-1
myrow=-1
mycol=-1
if(MyID>=idstart .and. MyID<=idend)then
	mycol = (MyID-idstart)/nprow
	myrow = mod(MyID-idstart,nprow)
endif
end subroutine Gridinfo_2D




end module MISC_Utilities
