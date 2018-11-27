#include "HODLR_config.fi"
module DenseLA
use BPACK_DEFS
use omp_lib



! interface geqrff90
  ! module procedure dgeqrff90
  ! module procedure zgeqrff90
! end interface


! interface geqp3f90
  ! module procedure dgeqp3f90
  ! module procedure zgeqp3f90
! end interface


! interface geqp3modf90
  ! module procedure dgeqp3modf90
  ! module procedure zgeqp3modf90
! end interface


! interface un_or_mqrf90
  ! module procedure ormqrf90
  ! module procedure unmqrf90
! end interface


! interface un_or_gqrf90
  ! module procedure ungqrf90
  ! module procedure orgqrf90
! end interface


! interface getrff90
  ! module procedure dgetrff90
  ! module procedure zgetrff90
! end interface

! interface getrsf90
  ! module procedure dgetrsf90
  ! module procedure zgetrsf90
! end interface

! interface getrif90
  ! module procedure dgetrif90
  ! module procedure zgetrif90
! end interface

! interface trsmf90
  ! module procedure dtrsmf90
  ! module procedure ztrsmf90
! end interface

! interface pun_or_mqrf90
  ! module procedure pzunmqrf90
  ! module procedure pdormqrf90
! end interface


! interface pgeqpfmodf90
  ! module procedure pzgeqpfmodf90
  ! module procedure pdgeqpfmodf90
! end interface


! interface pgetrif90
  ! module procedure pzgetrif90
  ! module procedure pdgetrif90
! end interface


! interface pgesvdf90
  ! module procedure pzgesvdf90
  ! module procedure pdgesvdf90
! end interface

contains



real(kind=8) function fnorm(Matrix,M,N,norm)
class(*) Matrix(:,:)
integer M,N
character,optional:: norm
select type(Matrix)
type is (real(kind=8))
	fnorm = dlangef90(Matrix,M,N,norm)	
type is (complex(kind=8))
	fnorm = zlangef90(Matrix,M,N,norm)	
end select

end function fnorm



real(kind=8) function dlangef90(Matrix,M,N,norm)	

! 

implicit none
real(kind=8) Matrix(:,:)
integer M,N,lda
character,optional:: norm
character::opt
real(kind=8),allocatable:: WORK(:)
real(kind=8) dlange

opt = 'F'
if(present(norm))opt=norm
if(opt=='I' .or. opt=='i')then
allocate(WORK(M))
WORK=0
end if

lda = size(Matrix,1)
dlangef90 = dlange( opt, M, N, Matrix, lda, WORK )

if(opt=='I' .or. opt=='i')then
deallocate(WORK)
end if

end function dlangef90	


real(kind=8) function zlangef90(Matrix,M,N,norm)	

! 

implicit none
complex(kind=8) Matrix(:,:)
integer M,N,lda
character,optional:: norm
character::opt
complex(kind=8),allocatable:: WORK(:)
real(kind=8) zlange

opt = 'F'
if(present(norm))opt=norm
if(opt=='I' .or. opt=='i')then
allocate(WORK(M))
WORK=0
end if

lda = size(Matrix,1)
zlangef90 = zlange( opt, M, N, Matrix, lda, WORK )

if(opt=='I' .or. opt=='i')then
deallocate(WORK)
end if

end function zlangef90	


subroutine gesvd_robust(Matrix,Singular,UU,VV,mm,nn,mn_min,flop)
    use BPACK_DEFS
	implicit none 
	integer mm,nn,mn_min
	DT Matrix(:,:),UU(:,:),VV(:,:)
	real(kind=8) Singular(:)
	real(kind=8),optional::flop
	

	if(mm==1)then
		UU(1,1)=1d0
		Singular(1) = fnorm(Matrix,mm,nn)
		if(Singular(1)<SafeUnderflow)then
			VV = 0d0
		else 
			VV = Matrix/Singular(1)
		endif
	elseif(nn==1)then
		VV(1,1)=1d0
		Singular(1) = fnorm(Matrix,mm,nn)
		if(Singular(1)<SafeUnderflow)then
			UU = 0d0
		else 
			UU = Matrix/Singular(1)
		endif			
	else
		if(fnorm(Matrix,mm,nn)<SafeUnderflow)then
			Singular=0
			UU=0
			VV=0
		else 
			! write(*,*)'ga',fnorm(Matrix,mm,nn),shape(Matrix),shape(UU),shape(VV)
			Singular=0
			UU=0
			VV=0
			
			
			call gesddf90(Matrix,Singular,UU,VV,flop=flop)
			
			!!!!!! gesvd (QR iteration) can occasionally fail compared to gesdd (DC) 
			! call gesvdf90(Matrix,Singular,UU,VV,flop=flop)  
			
		endif
	endif
	
	
end subroutine gesvd_robust
		

		
		

subroutine gesvdf90(Matrix,Singular,UU,VV,flop)
implicit none
class(*)::Matrix(:,:),UU(:,:),VV(:,:)
real(kind=8) Singular(:)
real(kind=8),optional::flop


select type(Matrix)
type is (real(kind=8))
	select type(UU)
	type is (real(kind=8))
	select type(VV)
	type is (real(kind=8))
		call dgesvdf90(Matrix,Singular,UU,VV,flop)
	end select
	end select
	
type is (complex(kind=8))
	select type(UU)
	type is (complex(kind=8))
	select type(VV)
	type is (complex(kind=8))
		call zgesvdf90(Matrix,Singular,UU,VV,flop)
	end select
	end select
	
end select

end subroutine gesvdf90
		
		
		
subroutine zgesvdf90(Matrix,Singular,UU,VV,flop)	

! 

implicit none
complex(kind=8)Matrix(:,:),UU(:,:),VV(:,:)
real(kind=8) Singular(:)
integer m,n,mn_min
real(kind=8),optional::flop

integer LWORK,INFO
complex(kind=8):: TEMP(1)
real(kind=8),allocatable::RWORK(:)

complex(kind=8),allocatable:: WORK(:)

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)


allocate(RWORK(5*mn_min))
RWORK=0
LWORK=-1
call ZGESVD('S','S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, TEMP, LWORK, RWORK, INFO)

LWORK=NINT(dble(TEMP(1)*2.001))  ! increase this 2.001 factor when not converging
allocate(WORK(LWORK))     
WORK=0

! write(*,*)'sssss', TEMP(1),LWORK,fnorm(Matrix,m,n)

call ZGESVD('S','S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, WORK, LWORK, RWORK, INFO)

if(INFO/=0)then
	write(*,*)'SVD failed!!',INFO
	stop
endif

deallocate(WORK,RWORK)

if(present(flop))flop = flops_zgesvd(m,n)

end subroutine zgesvdf90	


subroutine dgesvdf90(Matrix,Singular,UU,VV,flop)	

! 

implicit none
real(kind=8) Matrix(:,:),UU(:,:),VV(:,:)
real(kind=8) Singular(:)
integer m,n,mn_min

integer LWORK,INFO
real(kind=8):: TEMP(1)
real(kind=8),optional::flop

real(kind=8),allocatable:: WORK(:)

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)

LWORK=-1
call DGESVD('S','S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, TEMP, LWORK, INFO)

LWORK=NINT(dble(TEMP(1)*2.001))  ! increase this 2.001 factor when not converging
allocate(WORK(LWORK))     
WORK=0

! write(*,*)'sssss', TEMP(1),LWORK,fnorm(Matrix,m,n)

call DGESVD('S','S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, WORK, LWORK, INFO)

if(INFO/=0)then
	write(*,*)'SVD failed!!',INFO
	stop
endif

deallocate(WORK)
if(present(flop))flop = flops_dgesvd(m,n)
end subroutine dgesvdf90




subroutine gesddf90(Matrix,Singular,UU,VV,flop)
implicit none
class(*)::Matrix(:,:),UU(:,:),VV(:,:)
real(kind=8) Singular(:)
real(kind=8),optional::flop


select type(Matrix)
type is (real(kind=8))
	select type(UU)
	type is (real(kind=8))
	select type(VV)
	type is (real(kind=8))
		call dgesddf90(Matrix,Singular,UU,VV,flop)
	end select
	end select
	
type is (complex(kind=8))
	select type(UU)
	type is (complex(kind=8))
	select type(VV)
	type is (complex(kind=8))
		call zgesddf90(Matrix,Singular,UU,VV,flop)
	end select
	end select
	
end select

end subroutine gesddf90


subroutine zgesddf90(Matrix,Singular,UU,VV,flop)	

! 

implicit none
complex(kind=8)Matrix(:,:),UU(:,:),VV(:,:)
real(kind=8) Singular(:)
integer m,n,mn_min

integer LWORK,INFO
complex(kind=8):: TEMP(1)
real(kind=8),allocatable::RWORK(:)
integer,allocatable::IWORK(:)

complex(kind=8),allocatable:: WORK(:)
real(kind=8),optional::flop

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)

allocate(RWORK(max(1,min(m,n)*max(5*min(m,n)+7,2*max(m,n)+2*min(m,n)+1))))
allocate(IWORK(8*mn_min))
RWORK=0
LWORK=-1
call ZGESDD('S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, TEMP, LWORK, RWORK, IWORK, INFO)


LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK))     
WORK=0

! write(*,*)'sssss', TEMP(1),LWORK,fnorm(Matrix,m,n)

call ZGESDD('S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, WORK, LWORK, RWORK, IWORK, INFO)

if(INFO/=0)then
	write(*,*)'SDD failed!!',INFO
	stop
endif

deallocate(WORK,RWORK,IWORK)
if(present(flop))flop = flops_zgesdd(m,n)
end subroutine zgesddf90	



subroutine dgesddf90(Matrix,Singular,UU,VV,flop)	

! 

implicit none
real(kind=8) Matrix(:,:),UU(:,:),VV(:,:)
real(kind=8) Singular(:)
integer m,n,mn_min

integer LWORK,INFO
real(kind=8):: TEMP(1)
real(kind=8),allocatable::RWORK(:)
integer,allocatable::IWORK(:)

complex(kind=8),allocatable:: WORK(:)
real(kind=8),optional::flop

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)

allocate(IWORK(8*mn_min))
LWORK=-1
call DGESDD('S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, TEMP, LWORK, IWORK, INFO)


LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK))     
WORK=0

! write(*,*)'sssss', TEMP(1),LWORK,fnorm(Matrix,m,n)

call DGESDD('S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, WORK, LWORK, IWORK, INFO)

if(INFO/=0)then
	write(*,*)'SDD failed!!',INFO
	stop
endif

deallocate(WORK,IWORK)
if(present(flop))flop = flops_dgesdd(m,n)
end subroutine dgesddf90	




subroutine geqrff90(Matrix,tau,flop)
implicit none
class(*)::Matrix(:,:),tau(:)
real(kind=8),optional::flop


select type(Matrix)
type is (real(kind=8))
	select type(tau)
	type is (real(kind=8))
		call dgeqrff90(Matrix,tau,flop)
	end select
	
type is (complex(kind=8))
	select type(tau)
	type is (complex(kind=8))
		call zgeqrff90(Matrix,tau,flop)
	end select
end select

end subroutine geqrff90




subroutine zgeqrff90(Matrix,tau,flop)	

! 

implicit none
complex(kind=8)Matrix(:,:),tau(:)
integer m,n,mn_min
real(kind=8),optional::flop

integer LWORK,INFO
complex(kind=8):: TEMP(1)

complex(kind=8),allocatable:: WORK(:)

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)


LWORK=-1
call ZGEQRF(m, n, Matrix, m, tau, TEMP, LWORK, INFO)

LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK))     
WORK=0
call ZGEQRF(m, n, Matrix, m, tau, WORK, LWORK, INFO)


if(INFO/=0)then
	write(*,*)'zgeqrff90 failed!!',INFO
	stop
endif

deallocate(WORK)


! call geqrf(Matrix,tau)
if(present(flop))flop = flops_zgetrf(m,n)
end subroutine zgeqrff90




subroutine dgeqrff90(Matrix,tau,flop)	

! 

implicit none
real(kind=8) Matrix(:,:),tau(:)
integer m,n,mn_min

integer LWORK,INFO
real(kind=8):: TEMP(1)

real(kind=8),allocatable:: WORK(:)
real(kind=8),optional::flop

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)


LWORK=-1
call DGEQRF(m, n, Matrix, m, tau, TEMP, LWORK, INFO)

LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK))     
WORK=0
call DGEQRF(m, n, Matrix, m, tau, WORK, LWORK, INFO)


if(INFO/=0)then
	write(*,*)'dgeqrff90 failed!!',INFO
	stop
endif

deallocate(WORK)


! call geqrf(Matrix,tau)
if(present(flop))flop = flops_dgetrf(m,n)
end subroutine dgeqrff90





subroutine geqp3f90(Matrix,jpvt,tau,flop)	

! 

implicit none
class(*):: Matrix(:,:),tau(:)
integer jpvt(:)
real(kind=8),optional::flop

select type(Matrix)
type is (real(kind=8))
	select type(tau)
	type is (real(kind=8))
		call dgeqp3f90(Matrix,jpvt,tau,flop)	
	end select
	
type is (complex(kind=8))
	select type(tau)
	type is (complex(kind=8))
		call zgeqp3f90(Matrix,jpvt,tau,flop)	
	end select
end select

end subroutine geqp3f90


subroutine zgeqp3f90(Matrix,jpvt,tau,flop)	

! 

implicit none
complex(kind=8)Matrix(:,:),tau(:)
integer jpvt(:)
integer m,n,mn_min

integer LWORK,INFO
complex(kind=8):: TEMP(1)
real(kind=8),allocatable::RWORK(:)

complex(kind=8),allocatable:: WORK(:)

real(kind=8),optional::flop

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)

allocate(RWORK(2*max(m,n)))
RWORK=0
LWORK=-1
call ZGEQP3(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, RWORK, INFO)
LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK)) 
WORK=0    
call ZGEQP3(m, n, Matrix, m, jpvt, tau, WORK, LWORK, RWORK, INFO)

if(INFO/=0)then
	write(*,*)'geqp3f90 failed!!',INFO
	stop
endif

deallocate(WORK)
deallocate(RWORK)

! call geqp3(Matrix,jpvt,tau)
if(present(flop))flop = flops_zgetrf(m,n)
end subroutine zgeqp3f90
	


subroutine dgeqp3f90(Matrix,jpvt,tau,flop)	

! 

implicit none
real(kind=8) Matrix(:,:),tau(:)
integer jpvt(:)
integer m,n,mn_min

integer LWORK,INFO
real(kind=8):: TEMP(1)

real(kind=8),allocatable:: WORK(:)
real(kind=8),optional::flop

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)

LWORK=-1
call DGEQP3(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, INFO)
LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK)) 
WORK=0    
call DGEQP3(m, n, Matrix, m, jpvt, tau, WORK, LWORK, INFO)

if(INFO/=0)then
	write(*,*)'geqp3f90 failed!!',INFO
	stop
endif

deallocate(WORK)

! call geqp3(Matrix,jpvt,tau)
if(present(flop))flop = flops_dgetrf(m,n)
end subroutine dgeqp3f90	
	
	
subroutine geqp3modf90(Matrix,jpvt,tau,rtol,atol,rank,flop)	

! 

implicit none
class(*):: Matrix(:,:),tau(:)
integer jpvt(:)
real(kind=8)::rtol,atol
integer rank
real(kind=8),optional::flop

select type(Matrix)
type is (real(kind=8))
	select type(tau)
	type is (real(kind=8))
		call dgeqp3modf90(Matrix,jpvt,tau,rtol,atol,rank,flop)	
	end select
	
type is (complex(kind=8))
	select type(tau)
	type is (complex(kind=8))
		call zgeqp3modf90(Matrix,jpvt,tau,rtol,atol,rank,flop)	
	end select
end select

end subroutine geqp3modf90	
	
	

subroutine dgeqp3modf90(Matrix,jpvt,tau,rtol,atol,rank,flop)	

! 

implicit none
real(kind=8) Matrix(:,:),tau(:)
integer jpvt(:)
integer m,n,mn_min

real(kind=8)::rtol,atol
integer LWORK,INFO,rank
real(kind=8):: TEMP(1)
real(kind=8),optional::flop
real(kind=8),allocatable:: WORK(:)

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)
rank=0


LWORK=-1
! call DGEQP3(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, INFO)
call DGEQP3mod( m, n, Matrix, m, jpvt, tau, TEMP, LWORK,INFO, RANK, RTOL, ATOL)


LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK)) 
WORK=0    
! call DGEQP3(m, n, Matrix, m, jpvt, tau, WORK, LWORK, INFO)
call DGEQP3mod( m, n, Matrix, m, jpvt, tau, WORK, LWORK,INFO, RANK, RTOL, ATOL)


if(INFO/=0)then
	write(*,*)'geqp3modf90 failed!!',INFO
	stop
endif

deallocate(WORK)

! call geqp3(Matrix,jpvt,tau)
if(present(flop))flop = flops_dgeqpfmod(m,n,rank)
end subroutine dgeqp3modf90		
	

subroutine zgeqp3modf90(Matrix,jpvt,tau,rtol,atol,rank,flop)	

! 

implicit none
complex(kind=8)Matrix(:,:),tau(:)
integer jpvt(:)
integer m,n,mn_min

real(kind=8)::rtol,atol
integer LWORK,INFO,rank
complex(kind=8):: TEMP(1)
real(kind=8),allocatable::RWORK(:)

complex(kind=8),allocatable:: WORK(:)
real(kind=8),optional::flop

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)
rank=0

allocate(RWORK(2*max(m,n)))
RWORK=0
LWORK=-1
! call ZGEQP3(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, RWORK, INFO)
call ZGEQP3mod( m, n, Matrix, m, jpvt, tau, TEMP, LWORK, RWORK,INFO, RANK, RTOL, ATOL)


LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK)) 
WORK=0    
! call ZGEQP3(m, n, Matrix, m, jpvt, tau, WORK, LWORK, RWORK, INFO)
call ZGEQP3mod( m, n, Matrix, m, jpvt, tau, WORK, LWORK, RWORK,INFO, RANK, RTOL, ATOL)


if(INFO/=0)then
	write(*,*)'geqp3modf90 failed!!',INFO
	stop
endif

deallocate(WORK)
deallocate(RWORK)

! call geqp3(Matrix,jpvt,tau)
if(present(flop))flop = flops_zgeqpfmod(m,n,rank)
end subroutine zgeqp3modf90	
	

subroutine un_or_mqrf90(a,tau,c,side,trans,m,n,k,flop)	

! 

implicit none
class(*) a(:,:),tau(:),c(:,:)
integer m,n,k
character:: side, trans
real(kind=8),optional::flop

select type(a)
type is (real(kind=8))
	select type(tau)
	type is (real(kind=8))
	select type(c)
	type is (real(kind=8))
		call ormqrf90(a,tau,c,side,trans,m,n,k,flop)	
	end select
	end select
type is (complex(kind=8))
	select type(tau)
	type is (complex(kind=8))
	select type(c)
	type is (complex(kind=8))
		call unmqrf90(a,tau,c,side,trans,m,n,k,flop)	
	end select
	end select
end select


end	subroutine un_or_mqrf90
	
	
	
subroutine ormqrf90(a,tau,c,side,trans,m,n,k,flop)	

! 

implicit none
real(kind=8) a(:,:),tau(:),c(:,:)
integer m,n,k,lda,ldc, mn_min

character:: side, trans,trans1

integer LWORK,INFO
real(kind=8):: TEMP(1)
real(kind=8),allocatable:: WORK(:)
real(kind=8),optional::flop

trans1=trans
if(trans1=='C')trans1='T'  

ldc=size(c,1)
lda=size(a,1)

LWORK=-1
call DORMQR(side, trans1, m, n, k, a, lda, tau, c, ldc, TEMP, LWORK, INFO)
LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK))
WORK=0     
call DORMQR(side, trans1, m, n, k, a, lda, tau, c, ldc, WORK, LWORK, INFO)

if(INFO/=0)then
	write(*,*)'ormqrf90 failed!!',INFO
	stop
endif

deallocate(WORK)

! call unmqr(a,tau,c,side,trans)
if(present(flop))flop = flops_dunmqr(side,m,n,k)
end subroutine ormqrf90

	
	
subroutine unmqrf90(a,tau,c,side,trans,m,n,k,flop)	

! 

implicit none
complex(kind=8)a(:,:),tau(:),c(:,:)
integer m,n,k,lda,ldc, mn_min
real(kind=8),optional::flop
character:: side, trans

integer LWORK,INFO
complex(kind=8):: TEMP(1)
complex(kind=8),allocatable:: WORK(:)


ldc=size(c,1)
lda=size(a,1)

LWORK=-1
call ZUNMQR(side, trans, m, n, k, a, lda, tau, c, ldc, TEMP, LWORK, INFO)
LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK))
WORK=0     
call ZUNMQR(side, trans, m, n, k, a, lda, tau, c, ldc, WORK, LWORK, INFO)

if(INFO/=0)then
	write(*,*)'unmqrf90 failed!!',INFO
	stop
endif

deallocate(WORK)

! call unmqr(a,tau,c,side,trans)
if(present(flop))flop = flops_zunmqr(side,m,n,k)
end subroutine unmqrf90

	
	

subroutine un_or_gqrf90(Matrix,tau,m,n,k,flop)
class(*) Matrix(:,:),tau(:)
real(kind=8),optional::flop
integer m,n,k

select type(Matrix)
type is (real(kind=8))
	select type(tau)
	type is (real(kind=8))
		call orgqrf90(Matrix,tau,m,n,k,flop)
	end select
	
type is (complex(kind=8))
	select type(tau)
	type is (complex(kind=8))
		call ungqrf90(Matrix,tau,m,n,k,flop)
	end select
end select


end subroutine un_or_gqrf90



subroutine pun_or_gqrf90(Matrix,tau,m,n,k,desca,ia,ja,flop)
class(*) Matrix(:,:),tau(:)
real(kind=8),optional::flop
integer m,n,k,ia,ja
integer desca(9)

select type(Matrix)
type is (real(kind=8))
	select type(tau)
	type is (real(kind=8))
		call porgqrf90(Matrix,tau,m,n,k,desca,ia,ja,flop)
	end select
	
type is (complex(kind=8))
	select type(tau)
	type is (complex(kind=8))
		call pungqrf90(Matrix,tau,m,n,k,desca,ia,ja,flop)
	end select
end select


end subroutine pun_or_gqrf90

	
	

subroutine orgqrf90(Matrix,tau,m,n,k,flop)	

! 

implicit none
real(kind=8) Matrix(:,:),tau(:)
integer m,n,k,lda

integer LWORK,INFO
real(kind=8):: TEMP(1)

real(kind=8),allocatable:: WORK(:)
real(kind=8),optional::flop

lda=size(Matrix,1)


LWORK=-1
call DORGQR(m, n, k, Matrix, lda, tau, TEMP, LWORK, INFO)

LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK))     
WORK=0
call DORGQR(m, n, k, Matrix, lda, tau, WORK, LWORK, INFO)


if(INFO/=0)then
	write(*,*)'orgqrf90 failed!!',INFO
	stop
endif

deallocate(WORK)


! call ungqr(Matrix,tau)
if(present(flop))flop = flops_dungqr(m,n,k)
end subroutine orgqrf90	
	


subroutine porgqrf90(Matrix,tau,m,n,k,desca,ia,ja,flop)

! 

implicit none
real(kind=8) Matrix(:,:),tau(:)
integer m,n,k,lda

integer LWORK,INFO
real(kind=8):: TEMP(1)
integer desca(9)
integer ia,ja
real(kind=8),allocatable:: WORK(:)
real(kind=8),optional::flop


LWORK=-1
call PDORGQR(m, n, k, Matrix, ia, ja, desca, tau, TEMP, LWORK, INFO)

LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK))     
WORK=0
call PDORGQR(m, n, k, Matrix, ia, ja, desca, tau, WORK, LWORK, INFO)


if(INFO/=0)then
	write(*,*)'porgqrf90 failed!!',INFO
	stop
endif

deallocate(WORK)


! call ungqr(Matrix,tau)
if(present(flop))flop = flops_dungqr(m,n,k)
end subroutine porgqrf90		
	
	
	
subroutine ungqrf90(Matrix,tau,m,n,k,flop)	

! 

implicit none
complex(kind=8)Matrix(:,:),tau(:)
integer m,n,k,lda
real(kind=8),optional::flop
integer LWORK,INFO
complex(kind=8):: TEMP(1)

complex(kind=8),allocatable:: WORK(:)

lda=size(Matrix,1)

LWORK=-1
call ZUNGQR(m, n, k, Matrix, lda, tau, TEMP, LWORK, INFO)

LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK))     
WORK=0
call ZUNGQR(m, n, k, Matrix, lda, tau, WORK, LWORK, INFO)


if(INFO/=0)then
	write(*,*)'ungqrf90 failed!!',INFO
	stop
endif

deallocate(WORK)

if(present(flop))flop = flops_zungqr(m,n,k)
! call ungqr(Matrix,tau)

end subroutine ungqrf90	
	

	
subroutine pungqrf90(Matrix,tau,m,n,k,desca,ia,ja,flop)	

! 

implicit none
complex(kind=8)Matrix(:,:),tau(:)
integer m,n,k,lda
real(kind=8),optional::flop
integer LWORK,INFO
complex(kind=8):: TEMP(1)
integer desca(9)
integer ia,ja
complex(kind=8),allocatable:: WORK(:)

LWORK=-1
call PZUNGQR(m, n, k, Matrix, ia, ja, desca, tau, TEMP, LWORK, INFO)

LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK))     
WORK=0
call PZUNGQR(m, n, k, Matrix, ia, ja, desca, tau, WORK, LWORK, INFO)


if(INFO/=0)then
	write(*,*)'pungqrf90 failed!!',INFO
	stop
endif

deallocate(WORK)

if(present(flop))flop = flops_zungqr(m,n,k)
! call ungqr(Matrix,tau)

end subroutine pungqrf90	
	
	

subroutine getrff90(Matrix,ipiv,flop)	
class(*) Matrix(:,:)
integer ipiv(:)
real(kind=8),optional::flop

select type(Matrix)
type is (real(kind=8))
	call dgetrff90(Matrix,ipiv,flop)
type is (complex(kind=8))
	call zgetrff90(Matrix,ipiv,flop)
end select

end subroutine getrff90
	
		

subroutine dgetrff90(Matrix,ipiv,flop)	

! 

implicit none
real(kind=8) Matrix(:,:)
integer ipiv(:)
integer m,n,mn_min

integer INFO
real(kind=8),optional::flop

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)
	
call DGETRF( m, n, Matrix, m, ipiv, INFO )
	

if(INFO/=0)then
	write(*,*)'getrff90 failed!!',INFO
	stop
endif
	
	
! call getrf(Matrix,ipiv)	
if(present(flop))flop = flops_dgeqrf(m,n)
end subroutine dgetrff90		
	
	
subroutine zgetrff90(Matrix,ipiv,flop)	

! 

implicit none
complex(kind=8) Matrix(:,:)
integer ipiv(:)
integer m,n,mn_min
real(kind=8),optional::flop
integer INFO

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)
	
call ZGETRF( m, n, Matrix, m, ipiv, INFO )
	

if(INFO/=0)then
	write(*,*)'getrff90 failed!!',INFO
	stop
endif
	
	
! call getrf(Matrix,ipiv)	
if(present(flop))flop = flops_zgeqrf(m,n)
end subroutine zgetrff90	
	

	
subroutine getrsf90(Matrix,ipiv,B,trans,flop)	
class(*) Matrix(:,:),B(:,:)
integer ipiv(:)
character:: trans
real(kind=8),optional::flop

select type(Matrix)
type is (real(kind=8))
	select type(B)
	type is (real(kind=8))
	call dgetrsf90(Matrix,ipiv,B,trans,flop)		
	end select
type is (complex(kind=8))
	select type(B)
	type is (complex(kind=8))
	call zgetrsf90(Matrix,ipiv,B,trans,flop)		
	end select
end select

end subroutine getrsf90
		
	


subroutine dgetrsf90(Matrix,ipiv,B,trans,flop)	

! 

implicit none
real(kind=8) Matrix(:,:),B(:,:)
integer ipiv(:)
integer m,n,mn_min,nrhs
character:: trans,trans1
real(kind=8),optional::flop
integer INFO

n=size(Matrix,1)
nrhs=size(B,2)

trans1=trans
if(trans1=='C')trans1='T' 

call DGETRS( trans1, n, nrhs, Matrix, n, ipiv, B, n, INFO )


if(INFO/=0)then
	write(*,*)'getrsf90 failed!!',INFO
	stop
endif
	
	
! call getrs(Matrix,ipiv,B,trans)	
if(present(flop))flop = flops_dgetrs(n,nrhs)	
end subroutine dgetrsf90		
	
subroutine zgetrsf90(Matrix,ipiv,B,trans,flop)	

! 

implicit none
complex(kind=8) Matrix(:,:),B(:,:)
integer ipiv(:)
integer m,n,mn_min,nrhs
character:: trans
real(kind=8),optional::flop
integer INFO

n=size(Matrix,1)
nrhs=size(B,2)


call ZGETRS( trans, n, nrhs, Matrix, n, ipiv, B, n, INFO )


if(INFO/=0)then
	write(*,*)'getrsf90 failed!!',INFO
	stop
endif
	
	
! call getrs(Matrix,ipiv,B,trans)	
if(present(flop))flop = flops_zgetrs(n,nrhs)	
end subroutine zgetrsf90		
	
	
	

subroutine getrif90(Matrix,ipiv,flop)	

implicit none
class(*) Matrix(:,:)
integer ipiv(:)
real(kind=8),optional::flop	
	
select type(Matrix)
type is (real(kind=8))
	call dgetrif90(Matrix,ipiv,flop)			
type is (complex(kind=8))
	call zgetrif90(Matrix,ipiv,flop)		
end select	
	
end subroutine getrif90

subroutine dgetrif90(Matrix,ipiv,flop)	

! 

implicit none
real(kind=8) Matrix(:,:)
integer ipiv(:)
integer m,n,mn_min

integer LWORK,INFO
real(kind=8):: TEMP(1)

real(kind=8),allocatable:: WORK(:)
real(kind=8),optional::flop
m=size(Matrix,1)
n=size(Matrix,2)

if(m>n)then
	write(*,*)'size of Matrix incorrect'
	stop
endif

mn_min = min(m,n)


LWORK=-1
call DGETRI(m, Matrix, m, ipiv, TEMP, LWORK, INFO)

LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK))     
WORK=0
call DGETRI(m, Matrix, m, ipiv, WORK, LWORK, INFO)


if(INFO/=0)then
	write(*,*)'getrif90 failed!!',INFO
	stop
endif

deallocate(WORK)


! call getri(Matrix,ipiv)	
if(present(flop))flop = flops_dgetri(mn_min)
end subroutine dgetrif90	
	
	
subroutine zgetrif90(Matrix,ipiv,flop)	

! 

implicit none
complex(kind=8)Matrix(:,:)
integer ipiv(:)
integer m,n,mn_min

integer LWORK,INFO
complex(kind=8):: TEMP(1)
real(kind=8),optional::flop
complex(kind=8),allocatable:: WORK(:)

m=size(Matrix,1)
n=size(Matrix,2)
if(m>n)then
	write(*,*)'size of Matrix incorrect'
	stop
endif
mn_min = min(m,n)


LWORK=-1
call ZGETRI(m, Matrix, m, ipiv, TEMP, LWORK, INFO)

LWORK=NINT(dble(TEMP(1)*2.001))
allocate(WORK(LWORK))     
WORK=0
call ZGETRI(m, Matrix, m, ipiv, WORK, LWORK, INFO)


if(INFO/=0)then
	write(*,*)'getrif90 failed!!',INFO
	stop
endif

deallocate(WORK)


! call getri(Matrix,ipiv)	
if(present(flop))flop = flops_zgetri(mn_min)
end subroutine zgetrif90		
	


	

subroutine trsmf90(Matrix,B,side,uplo,transa,diag,m,n,flop)	

implicit none
class(*) Matrix(:,:),B(:,:)
integer m,n
character side,uplo,transa,diag
real(kind=8),optional::flop

select type(Matrix)
type is (real(kind=8))
	select type(B)
	type is (real(kind=8))
	call dtrsmf90(Matrix,B,side,uplo,transa,diag,m,n,flop)	
	end select
type is (complex(kind=8))
	select type(B)
	type is (complex(kind=8))
	call ztrsmf90(Matrix,B,side,uplo,transa,diag,m,n,flop)	
	end select
end select


end subroutine trsmf90
	


subroutine dtrsmf90(Matrix,B,side,uplo,transa,diag,m,n,flop)	

! 

implicit none
real(kind=8) Matrix(:,:),B(:,:),alpha
integer m,n,lda,ldb
character side,uplo,transa,diag
real(kind=8),optional::flop
integer INFO

alpha=1d0

ldb=size(B,1)
lda=size(Matrix,1)

call DTRSM(side, uplo, transa, diag, m, n, alpha, Matrix, lda, B, ldb)


! call trsm(Matrix,B,side,uplo,transa,diag)	
if(present(flop))flop = flops_dtrsm(side,m,n)
end subroutine dtrsmf90		


subroutine ztrsmf90(Matrix,B,side,uplo,transa,diag,m,n,flop)	

! 

implicit none
complex(kind=8) Matrix(:,:),B(:,:),alpha
integer m,n,lda,ldb
character side,uplo,transa,diag
real(kind=8),optional::flop
integer INFO

alpha=1d0

ldb=size(B,1)
lda=size(Matrix,1)

call ZTRSM(side, uplo, transa, diag, m, n, alpha, Matrix, lda, B, ldb)


! call trsm(Matrix,B,side,uplo,transa,diag)	
if(present(flop))flop = flops_ztrsm(side,m,n)
end subroutine ztrsmf90		



subroutine gemmf90(MatA,lda,MatB,ldb,MatC,ldc,transa,transb,m,n,k,al,be,flop)
integer m,n,k,lda,ldb,ldc
class(*) MatA(:,:),MatB(:,:),MatC(:,:)
class(*),optional:: al,be
character transa,transb
real(kind=8),optional::flop
real(kind=8)alphar,betar
complex(kind=8)alphac,betac

! use default value if not presentd
alphar=1d0
betar=0d0
alphac=1d0
betac=0d0

select type(MatA)
type is (real(kind=8))
	select type(MatB)
	type is (real(kind=8))
	select type(MatC)
	type is (real(kind=8))
	
	if(present(al))then
	select type(al)
	type is (real(kind=8))
		alphar=al
	end select
	endif

	if(present(be))then
	select type(be)
	type is (real(kind=8))
		betar=be
	end select
	endif	
	call dgemmf90(MatA,lda,MatB,ldb,MatC,ldc,transa,transb,m,n,k,alphar,betar,flop)
	
	end select		
	end select
	
type is (complex(kind=8))
	select type(MatB)
	type is (complex(kind=8))
	select type(MatC)
	type is (complex(kind=8))
	
	if(present(al))then
	select type(al)
	type is (complex(kind=8))
		alphac=al
	end select
	endif

	if(present(be))then
	select type(be)
	type is (complex(kind=8))
		betac=be
	end select
	endif	
	call zgemmf90(MatA,lda,MatB,ldb,MatC,ldc,transa,transb,m,n,k,alphac,betac,flop)
	
	end select		
	end select
end select


end subroutine gemmf90



subroutine dgemmf90(MatA,lda,MatB,ldb,MatC,ldc,transa,transb,m,n,k,alpha,beta,flop)	

! 

implicit none
integer m,n,k,lda,ldb,ldc
real(kind=8) MatA(:,:),MatB(:,:),MatC(:,:)
real(kind=8):: alpha,beta
character transa,transb
real(kind=8),optional::flop

call dgemm(transa, transb, m, n, k, alpha, MatA, lda, MatB, ldb, beta, MatC, ldc)

! call gemmf90(MatA,MatB,MatC,transa,transb,alpha,beta)	
if(present(flop))flop = flops_dgemm(m,n,k)
end subroutine dgemmf90



subroutine zgemmf90(MatA,lda,MatB,ldb,MatC,ldc,transa,transb,m,n,k,alpha,beta,flop)	

! 

implicit none
integer m,n,k,lda,ldb,ldc
complex(kind=8) MatA(:,:),MatB(:,:),MatC(:,:)
complex(kind=8):: alpha,beta
character transa,transb
real(kind=8),optional::flop

call zgemm(transa, transb, m, n, k, alpha, MatA, lda, MatB, ldb, beta, MatC, ldc)

! call gemmf90(MatA,MatB,MatC,transa,transb,alpha,beta)	
if(present(flop))flop = flops_zgemm(m,n,k)				
end subroutine zgemmf90





subroutine pun_or_mqrf90(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc,flop)
implicit none
character side,trans
integer m,n,k,ia,ja,ic,jc
class(*) MatA(:,:),MatC(:,:),tau(:)
integer desca(9),descc(9)
real(kind=8),optional::flop


select type(MatA)
type is (real(kind=8))
	select type(MatC)
	type is (real(kind=8))
	select type(tau)
	type is (real(kind=8))	
	call pdormqrf90(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc,flop)	
	end select
	end select
type is (complex(kind=8))
	select type(MatC)
	type is (complex(kind=8))
	select type(tau)
	type is (complex(kind=8))	
	call pzunmqrf90(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc,flop)
	end select
	end select
end select

end subroutine pun_or_mqrf90




subroutine pdormqrf90(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc,flop)
implicit none
character side,trans
integer m,n,k,ia,ja,ic,jc
real(kind=8) MatA(:,:),MatC(:,:),tau(:)
integer desca(9),descc(9)
real(kind=8),allocatable:: WORK(:)
integer LWORK,INFO
real(kind=8):: TEMP(1)
real(kind=8),optional::flop
LWORK=-1
call pdormqr(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, TEMP, lwork, info)
lwork=NINT(dble(TEMP(1)*2.001))
allocate(WORK(lwork))     
WORK=0		
call pdormqr(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, WORK, lwork, info)
deallocate(WORK)

if(present(flop))flop = flops_dunmqr(side,m,n,k)

end subroutine pdormqrf90


subroutine pzunmqrf90(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc,flop)
implicit none
character side,trans
integer m,n,k,ia,ja,ic,jc
complex(kind=8) MatA(:,:),MatC(:,:),tau(:)
integer desca(9),descc(9)
complex(kind=8),allocatable:: WORK(:)
integer LWORK,INFO
complex(kind=8):: TEMP(1)
real(kind=8),optional::flop

LWORK=-1
call pzunmqr(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, TEMP, lwork, info)
lwork=NINT(dble(TEMP(1)*2.001))
allocate(WORK(lwork))     
WORK=0		
call pzunmqr(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, WORK, lwork, info)

deallocate(WORK)
if(present(flop))flop = flops_zunmqr(side,m,n,k)
end subroutine pzunmqrf90






subroutine pgeqpfmodf90(M, N, Matrix, ia, ja, desca, ipiv, tau, JPERM, jpiv, rank,rtol, atol,flop)
implicit none
integer M,N,ia,ja
class(*) Matrix(:,:),tau(:)
integer ipiv(:),jpiv(:),JPERM(:)
integer rank
integer desca(9)
real(kind=8)::rtol,atol
real(kind=8),optional::flop

select type(Matrix)
type is (real(kind=8))
	select type(tau)
	type is (real(kind=8))	
	call pdgeqpfmodf90(M, N, Matrix, ia, ja, desca, ipiv, tau, JPERM, jpiv, rank,rtol, atol,flop)	
	end select
type is (complex(kind=8))
	select type(tau)
	type is (complex(kind=8))	
	call pzgeqpfmodf90(M, N, Matrix, ia, ja, desca, ipiv, tau, JPERM, jpiv, rank,rtol, atol,flop)
	end select
end select

end subroutine pgeqpfmodf90



subroutine pzgeqpfmodf90(M, N, Matrix, ia, ja, desca, ipiv, tau, JPERM, jpiv, rank,rtol, atol,flop)
implicit none
integer M,N,ia,ja
complex(kind=8) Matrix(:,:),tau(:)
integer ipiv(:),jpiv(:),JPERM(:)
integer rank
integer desca(9)
real(kind=8)::rtol,atol
integer LWORK,LRWORK,INFO,ierr
complex(kind=8):: TEMP(1)
real(kind=8),allocatable::RWORK(:)
complex(kind=8),allocatable:: WORK(:)
real(kind=8):: RTEMP(1)
real(kind=8),optional::flop

LWORK=-1
LRWORK=-1
call PZGEQPFmod(M, N, Matrix, 1, 1, desca, ipiv, tau, TEMP, lwork, RTEMP, lrwork, info, JPERM, jpiv, rank,rtol, atol)
lwork=NINT(dble(TEMP(1)*2.001))
allocate(WORK(lwork))     
WORK=0
lrwork=NINT(dble(RTEMP(1)*2.001))
allocate(RWORK(lrwork))     
RWORK=0		
call PZGEQPFmod(M, N, Matrix, 1, 1, desca, ipiv, tau, WORK, lwork, RWORK, lrwork, info, JPERM, jpiv, rank,rtol, atol)

deallocate(WORK)
deallocate(RWORK)

if(present(flop))flop = flops_zgeqpfmod(m,n,rank)
end subroutine pzgeqpfmodf90


subroutine pdgeqpfmodf90(M, N, Matrix, ia, ja, desca, ipiv, tau, JPERM, jpiv, rank,rtol, atol,flop)
implicit none
integer M,N,ia,ja
real(kind=8) Matrix(:,:),tau(:)
integer ipiv(:),jpiv(:),JPERM(:)
integer rank
integer desca(9)
real(kind=8)::rtol,atol
integer LWORK,INFO,ierr
real(kind=8):: TEMP(1)
real(kind=8),allocatable:: WORK(:)

real(kind=8),optional::flop
LWORK=-1

call PDGEQPFmod(M, N, Matrix, 1, 1, desca, ipiv, tau, TEMP, lwork, info, JPERM, jpiv, rank,rtol, atol)
lwork=NINT(dble(TEMP(1)*2.001))
allocate(WORK(lwork))     
WORK=0	
call PDGEQPFmod(M, N, Matrix, 1, 1, desca, ipiv, tau, WORK, lwork, info, JPERM, jpiv, rank,rtol, atol)

deallocate(WORK)
if(present(flop))flop = flops_dgeqpfmod(m,n,rank)

end subroutine pdgeqpfmodf90




subroutine pgemr2df90(M, N, MatA, ia, ja, desca, MatB, ib, jb, descb, ictxt)
implicit none
integer M,N,ia,ja,ib,jb
class(*) MatA(:,:),MatB(:,:)
integer desca(9),descb(9)
integer ictxt

select type(MatA)
type is (real(kind=8))
	select type(MatB)
	type is (real(kind=8))
	call pdgemr2d(M, N, MatA, ia, ja, desca, MatB, ib, jb, descb, ictxt)	
	end select
type is (complex(kind=8))
	select type(MatB)
	type is (complex(kind=8))
	call pzgemr2d(M, N, MatA, ia, ja, desca, MatB, ib, jb, descb, ictxt)
	end select
end select
end subroutine pgemr2df90	




subroutine pgemmf90(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc,flop)
implicit none
character transa,transb
integer m,n,k,ia,ja,ib,jb,ic,jc
class(*) alpha,beta,a(:,:),b(:,:),c(:,:)
integer desca(9),descb(9),descc(9)
real(kind=8),optional::flop
select type(a)
type is (real(kind=8))
	select type(b)
	type is (real(kind=8))
	select type(c)
	type is (real(kind=8))
	select type(alpha)
	type is (real(kind=8))	
	select type(beta)
	type is (real(kind=8))		
	call pdgemm(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
	if(present(flop))flop = flops_dgemm(m,n,k)
	end select
	end select	
	end select	
	end select
type is (complex(kind=8))
	select type(b)
	type is (complex(kind=8))
	select type(c)
	type is (complex(kind=8))
	select type(alpha)
	type is (complex(kind=8))
	select type(beta)
	type is (complex(kind=8))	
	call pzgemm(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
	if(present(flop))flop = flops_zgemm(m,n,k)
	end select
	end select	
	end select	
	end select
end select
end subroutine pgemmf90


subroutine ptrsmf90(side, uplo, transa, diag, m, n, alpha, a, ia, ja, desca, b, ib, jb, descb,flop)
implicit none
character side, uplo,transa,diag
integer m,n,ia,ja,ib,jb
integer desca(9),descb(9)
class(*) a(:,:),b(:,:),alpha
real(kind=8),optional::flop

select type(a)
type is (real(kind=8))
	select type(b)
	type is (real(kind=8))
	select type(alpha)
	type is (real(kind=8))	
	call pdtrsm(side, uplo, transa, diag, m, n, alpha, a, ia, ja, desca, b, ib, jb, descb)
	if(present(flop))flop = flops_dtrsm(side, m,n)
	end select
	end select	
type is (complex(kind=8))
	select type(b)
	type is (complex(kind=8))
	select type(alpha)
	type is (complex(kind=8))
	call pztrsm(side, uplo, transa, diag, m, n, alpha, a, ia, ja, desca, b, ib, jb, descb)
	if(present(flop))flop = flops_ztrsm(side, m,n)
	end select
	end select		
end select
end subroutine ptrsmf90


subroutine pgetrff90(m, n, a, ia, ja, desca, ipiv, info,flop)
implicit none
integer m,n,ia,ja
class(*) a(:,:)
integer desca(9)
integer ipiv(:)
integer info
real(kind=8),optional::flop

select type(a)
type is (real(kind=8))
	call pdgetrf(m, n, a, ia, ja, desca, ipiv, info)
	if(present(flop))flop = flops_dgetrf(m,n)	
type is (complex(kind=8))
	call pzgetrf(m, n, a, ia, ja, desca, ipiv, info)	
	if(present(flop))flop = flops_zgetrf(m,n)		
end select
end subroutine pgetrff90	


subroutine pgetrif90(n, a, ia, ja, desca, ipiv,flop)
implicit none
integer n,ia,ja
class(*):: a(:,:)
integer desca(9)
integer ipiv(:)
real(kind=8),optional::flop

select type(a)
type is (real(kind=8))
	call pdgetrif90(n, a, ia, ja, desca, ipiv,flop)
type is (complex(kind=8))
	call pzgetrif90(n, a, ia, ja, desca, ipiv,flop)
end select

end subroutine pgetrif90


subroutine pdgetrif90(n, a, ia, ja, desca, ipiv,flop)
implicit none
integer n,ia,ja
real(kind=8):: a(:,:)
real(kind=8):: TEMP(1)
integer desca(9)
integer ipiv(:)
integer info
integer TEMPI(1)
integer lwork,liwork
integer, allocatable :: iwork(:)
real(kind=8),allocatable:: work(:)
real(kind=8),optional::flop

call pdgetri(n,a,1,1,desca,ipiv,TEMP,-1,TEMPI,-1,info)
LWORK=NINT(dble(TEMP(1)*2.001))
allocate(work(lwork))
work=0
liwork=TEMPI(1)
allocate(iwork(liwork))
iwork=0	
call pdgetri(n,a,1,1,desca,ipiv,work,lwork,iwork,liwork,info)	
deallocate(iwork)
deallocate(work)

if(present(flop))flop = flops_dgetri(n)

end subroutine pdgetrif90





subroutine pzgetrif90(n, a, ia, ja, desca, ipiv,flop)
implicit none
integer n,ia,ja
complex(kind=8):: a(:,:)
complex(kind=8):: TEMP(1)
integer desca(9)
integer ipiv(:)
integer info
integer TEMPI(1)
integer lwork,liwork
integer, allocatable :: iwork(:)
complex(kind=8),allocatable:: work(:)
real(kind=8),optional::flop

call pzgetri(n,a,1,1,desca,ipiv,TEMP,-1,TEMPI,-1,info)
LWORK=NINT(dble(TEMP(1)*2.001))
allocate(work(lwork))
work=0
liwork=TEMPI(1)
allocate(iwork(liwork))
iwork=0	
call pzgetri(n,a,1,1,desca,ipiv,work,lwork,iwork,liwork,info)
	
deallocate(iwork)
deallocate(work)

if(present(flop))flop = flops_zgetri(n)

end subroutine pzgetrif90



subroutine pgesvdf90(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt,flop)
implicit none

character jobu,jobvt
integer m,n,ia,ja,iu,ju,ivt,jvt
class(*):: a(:,:),u(:,:),vt(:,:)
integer desca(9),descu(9),descvt(9)
real(kind=8):: s(:)
real(kind=8),optional::flop


select type(a)
type is (real(kind=8))
	select type(u)
	type is (real(kind=8))
	select type(vt)
	type is (real(kind=8))	
	call pdgesvdf90(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt,flop)
	end select
	end select	
type is (complex(kind=8))
	select type(u)
	type is (complex(kind=8))
	select type(vt)
	type is (complex(kind=8))
	call pzgesvdf90(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt,flop)
	end select
	end select		
end select

end subroutine pgesvdf90



subroutine pdgesvdf90(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt,flop)
implicit none

character jobu,jobvt
integer m,n,ia,ja,iu,ju,ivt,jvt
real(kind=8):: a(:,:),u(:,:),vt(:,:)
integer desca(9),descu(9),descvt(9)
real(kind=8):: s(:)
real(kind=8):: TEMP(1)
integer LWORK,mnmax,mnmin
real(kind=8),allocatable:: WORK(:)	
integer info
real(kind=8),optional::flop
mnmax = max(m,n)
mnmin = min(m,n)

lwork=-1
call pdgesvd(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, TEMP, lwork, info)

lwork=NINT(dble(TEMP(1)*2.001))
allocate(WORK(lwork))     
WORK=0
call pdgesvd(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, WORK, lwork, info)

deallocate(WORK)	
if(present(flop))flop = flops_dgesvd(m,n)

end subroutine pdgesvdf90



subroutine pzgesvdf90(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt,flop)
implicit none

character jobu,jobvt
integer m,n,ia,ja,iu,ju,ivt,jvt
complex(kind=8):: a(:,:),u(:,:),vt(:,:)
integer desca(9),descu(9),descvt(9)
real(kind=8):: s(:)
complex(kind=8):: TEMP(1)
integer LWORK,mnmax,mnmin
complex(kind=8),allocatable:: WORK(:)	
real(kind=8),allocatable::RWORK(:)
integer info
real(kind=8),optional::flop

mnmax = max(m,n)
mnmin = min(m,n)

allocate(rwork(1+4*mnmax))
rwork=0
lwork=-1
call pzgesvd(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, TEMP, lwork, rwork, info)

lwork=NINT(dble(TEMP(1)*2.001))
allocate(WORK(lwork))     
WORK=0
call pzgesvd(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, WORK, lwork, rwork, info)


deallocate(WORK,rwork)	
if(present(flop))flop = flops_zgesvd(m,n)

end subroutine pzgesvdf90


real(kind=8) function flops_zgesdd(m, n)
	implicit none 
	integer m,n
	if(m>n)then	
		flops_zgesdd = 4.*(8.*m*n*n + 4./3.*n*n*n)
	else
		flops_zgesdd = 4.*(8.*n*m*m + 4./3.*m*m*m)
	endif
end function flops_zgesdd

real(kind=8) function flops_dgesdd(m, n)
	implicit none 
	integer m,n
	if(m>n)then	
		flops_dgesdd = 8.*m*n*n + 4./3.*n*n*n
	else
		flops_dgesdd = 8.*n*m*m + 4./3.*m*m*m
	endif
end function flops_dgesdd


real(kind=8) function flops_zgesvd(m, n)
	implicit none 
	integer m,n
	if(m>n)then	
		flops_zgesvd = 4.*(12.*m*n*n + 16./3.*n*n*n)
	else
		flops_zgesvd = 4.*(12.*n*m*m + 16./3.*m*m*m)
	endif
end function flops_zgesvd


real(kind=8) function flops_dgesvd(m, n)
	implicit none 
	integer m,n
	if(m>n)then	
		flops_dgesvd = 12.*m*n*n + 16./3.*n*n*n
	else
		flops_dgesvd = 12.*n*m*m + 16./3.*m*m*m
	endif
end function flops_dgesvd


real(kind=8) function flops_dgeqpfmod(m, n, k)
	implicit none 
	integer m,n,k
	if(m>n)then	
		flops_dgeqpfmod = 2.*m*n*n - 2./3.*n*n*n - (2.*(m-k)*(n-k)*(n-k) - 2./3.*(n-k)*(n-k)*(n-k))
	else
		flops_dgeqpfmod = 2.*n*m*m - 2./3.*m*m*m - (2.*(n-k)*(m-k)*(m-k) - 2./3.*(m-k)*(m-k)*(m-k))
	endif
end function flops_dgeqpfmod

real(kind=8) function flops_zgeqpfmod(m, n, k)
	implicit none 
	integer m,n,k
	if(m>n)then	
		flops_zgeqpfmod = 4.*(2.*m*n*n - 2./3.*n*n*n - (2.*(m-k)*(n-k)*(n-k) - 2./3.*(n-k)*(n-k)*(n-k)))
	else
		flops_zgeqpfmod = 4.*(2.*n*m*m - 2./3.*m*m*m - (2.*(n-k)*(m-k)*(m-k) - 2./3.*(m-k)*(m-k)*(m-k)))
	endif
end function flops_zgeqpfmod



real(kind=8) function fmuls_geqrf(m, n)
	implicit none 
	integer m,n
	if(m>n)then	
		fmuls_geqrf = m*n*n - 1./3.*n*n*n +   m*n + 0.5*n*n + 23./6.*n
	else
		fmuls_geqrf = n*m*m - 1./3.*m*m*m + 2*n*m - 0.5*m*m + 23./6.*m
	endif
end function fmuls_geqrf
real(kind=8) function fadds_geqrf(m, n)
	implicit none 
	integer m,n
	if(m>n)then	
		fadds_geqrf = m*n*n - 1./3.*n*n*n + 0.5*n*n       + 5./6.*n
	else
		fadds_geqrf = n*m*m - 1./3.*m*m*m + n*m - 0.5*m*m + 5./6.*m
	endif
end function fadds_geqrf
real(kind=8) function flops_zgeqrf(m, n)
	implicit none 
	integer m,n
	flops_zgeqrf = 6.*fmuls_geqrf(m, n) + 2.*fadds_geqrf(m, n)
end function flops_zgeqrf
real(kind=8) function flops_dgeqrf(m, n)
	implicit none 
	integer m,n
	flops_dgeqrf = fmuls_geqrf(m, n) + fadds_geqrf(m, n)
end function flops_dgeqrf



real(kind=8) function fmuls_ungqr(m, n, k)
	implicit none 
	integer m,n,k
    fmuls_ungqr = 2.*m*n*k - (m + n)*k*k + 2./3.*k*k*k + 2.*n*k - k*k - 5./3.*k
end function fmuls_ungqr
real(kind=8) function fadds_ungqr(m, n, k)
	implicit none 
	integer m,n,k
    fadds_ungqr = 2.*m*n*k - (m + n)*k*k + 2./3.*k*k*k + n*k - m*k + 1./3.*k
end function fadds_ungqr
real(kind=8) function flops_zungqr(m, n,k)
	implicit none 
	integer m,n,k
	flops_zungqr = 6.*fmuls_ungqr(m, n,k) + 2.*fadds_ungqr(m, n,k)
end function flops_zungqr
real(kind=8) function flops_dungqr(m, n,k)
	implicit none 
	integer m,n,k
	flops_dungqr = fmuls_ungqr(m, n,k) + fadds_ungqr(m, n,k)
end function flops_dungqr



real(kind=8) function fmuls_unmqr(side, m, n, k)
	integer m,n,k
	character side
	if(side=='L')then
		fmuls_unmqr = 2.*n*m*k - n*k*k + 2.*n*k
	else
		fmuls_unmqr = 2.*n*m*k - m*k*k + m*k + n*k - 0.5*k*k + 0.5*k
	endif
end function fmuls_unmqr

real(kind=8) function fadds_unmqr(side, m, n, k)
	integer m,n,k
	character side
	if(side=='L')then
		fadds_unmqr = 2.*n*m*k - n*k*k + n*k
	else
		fadds_unmqr = 2.*n*m*k - m*k*k + m*k
	endif
end function fadds_unmqr

real(kind=8) function flops_zunmqr(side, m, n, k)
	integer m,n,k
	character side
	flops_zunmqr = 6.*fmuls_unmqr(side, m, n, k) + 2.*fadds_unmqr(side, m, n, k)
end function flops_zunmqr


real(kind=8) function flops_dunmqr(side, m, n, k)
	integer m,n,k
	character side
	flops_dunmqr = fmuls_unmqr(side, m, n, k) + fadds_unmqr(side, m, n, k)
end function flops_dunmqr


real(kind=8) function fmuls_getrf(m, n)
	implicit none 
	integer m,n
	
	if(m>n)then	
		fmuls_getrf = 0.5*m*n*n - 1./6.*n*n*n + 0.5*m*n - 0.5*n*n + 2./3.*n
	else
		fmuls_getrf = 0.5*n*m*m - 1./6.*m*m*m + 0.5*n*m - 0.5*m*m + 2./3.*m
	endif	
end function fmuls_getrf
real(kind=8) function fadds_getrf(m, n)
	implicit none 
	integer m,n
	
	if(m>n)then	
		fadds_getrf = 0.5*m*n*n - 1./6.*n*n*n - 0.5*m*n + 1./6.*n
	else
		fadds_getrf = 0.5*n*m*m - 1./6.*m*m*m - 0.5*n*m + 1./6.*m
	endif	
end function fadds_getrf
real(kind=8) function flops_zgetrf(m, n)
	implicit none 
	integer m,n
	flops_zgetrf = 6.*fmuls_getrf(m, n) + 2.*fadds_getrf(m, n)
end function flops_zgetrf
real(kind=8) function flops_dgetrf(m, n)
	implicit none 
	integer m,n
	flops_dgetrf = fmuls_getrf(m, n) + fadds_getrf(m, n)
end function flops_dgetrf




real(kind=8) function fmuls_getrs(n, nrhs)
	implicit none
	integer n,nrhs
    fmuls_getrs =  nrhs*n*n
end function fmuls_getrs	
real(kind=8) function fadds_getrs(n, nrhs)
	implicit none
	integer n,nrhs
    fadds_getrs =  nrhs*n*(n - 1)
end function fadds_getrs	
real(kind=8) function flops_zgetrs(n, nrhs)
	implicit none 
	integer n,nrhs
	flops_zgetrs = 6.*fmuls_getrs(n,nrhs) + 2.*fadds_getrs(n,nrhs)
end function flops_zgetrs
real(kind=8) function flops_dgetrs(n, nrhs)
	implicit none 
	integer n,nrhs
	flops_dgetrs = fmuls_getrs(n,nrhs) + fadds_getrs(n,nrhs)
end function flops_dgetrs



real(kind=8) function fmuls_getri(n)
	implicit none 
	integer n
    fmuls_getri = 2./3.*n*n*n + 0.5*n*n + 5./6.*n
end function fmuls_getri	
real(kind=8) function fadds_getri(n)
	implicit none 
	integer n
    fadds_getri = 2./3.*n*n*n - 1.5*n*n + 5./6.*n
end function fadds_getri	
real(kind=8) function flops_zgetri(n)
	implicit none 
	integer n
	flops_zgetri = 6.*fmuls_getri(n) + 2.*fadds_getri(n)
end function flops_zgetri
real(kind=8) function flops_dgetri(n)
	implicit none 
	integer n
	flops_dgetri = fmuls_getri(n) + fadds_getri(n)
end function flops_dgetri


real(kind=8) function fmuls_trsm(side, m, n)
	integer m,n
	character side
	if(side=='L')then
		fmuls_trsm = 0.5*n*m*(m + 1)
	elseif(side=='R')then
		fmuls_trsm = 0.5*m*n*(n + 1)
	endif
end function fmuls_trsm	
real(kind=8) function fadds_trsm(side, m, n)
	integer m,n
	character side
	if(side=='L')then
		fadds_trsm = 0.5*n*m*(m - 1)
	elseif(side=='R')then
		fadds_trsm = 0.5*m*n*(n - 1)
	endif
end function fadds_trsm	
real(kind=8) function flops_ztrsm(side, m, n)
	integer m,n
	character side
	flops_ztrsm = 6.*fmuls_trsm(side, m, n) + 2.*fadds_trsm(side, m, n)
end function flops_ztrsm
real(kind=8) function flops_dtrsm(side, m, n)
	integer m,n
	character side
	flops_dtrsm = fmuls_trsm(side, m, n) + fadds_trsm(side, m, n)
end function flops_dtrsm


real(kind=8) function fmuls_gemm(m, n, k)
	implicit none 
	integer m,n,k
    fmuls_gemm = m*n*k
end function fmuls_gemm
real(kind=8) function fadds_gemm(m, n, k)
	implicit none 
	integer m,n,k
    fadds_gemm = m*n*k
end function fadds_gemm
real(kind=8) function flops_zgemm(m, n,k)
	implicit none 
	integer m,n,k
	flops_zgemm = 6.*fmuls_gemm(m, n,k) + 2.*fadds_gemm(m, n,k)
end function flops_zgemm
real(kind=8) function flops_dgemm(m, n,k)
	implicit none 
	integer m,n,k
	flops_dgemm = fmuls_gemm(m, n,k) + fadds_gemm(m, n,k)
end function flops_dgemm



end module DenseLA