! include "mkl_vsl.f90"

module misc
use MODULE_FILE
#ifdef Intel
USE IFPORT    !! uncomment this when using ifort
#endif  
use omp_lib

integer, parameter :: int64 = selected_int_kind(18) 

contains

subroutine linspaceI(startI,endI,N,array)
implicit none
integer startI,endI,N
integer array(1:N)
real*8::rtemp
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


subroutine check_NAN(A,M,N,nan)
! ! use lapack95
! ! use blas95
implicit none 
	
	logical::nan
	integer::M,N
	complex(kind=8)::A(M,N),ctemp
	integer ii,jj
	
	ctemp = 0
	do ii=1,M
	do jj=1,N
		ctemp = ctemp + A(ii,jj)
	end do
	end do
	nan = isnan(abs(ctemp))

 end subroutine check_NAN


 subroutine gemm_omp(A,B,C,m,n,k)
 	! ! use lapack95
	! ! use blas95
	implicit none 
	
	integer m,n,k
	complex(kind=8)::A(m,n),B(n,k),C(m,k)
	complex(kind=8),allocatable::A1(:,:),B1(:,:)
	complex(kind=8)::ctemp
	real*8::n1,n2
	integer ii,jj,kk
	
	allocate(A1(m,n))
	allocate(B1(n,k))
	
	! n1 = OMP_get_wtime()	

	if(isnan(fnorm(A,m,n)) .or. isnan(fnorm(B,n,k)))then
		write(*,*)fnorm(A,m,n),fnorm(B,n,k),'diaod'
		stop
	end if
	
	A1 = A
	B1 = B
	
	!$omp parallel do default(shared) private(ii,jj,kk,ctemp)	
	do ii =1,m	
	do jj =1,k	
	ctemp = 0
	do kk=1,n
		ctemp = ctemp + A1(ii,kk)*B1(kk,jj)	
	end do
	C(ii,jj) = ctemp
	end do
	end do
	!$omp end parallel do
	 
	! ! call gemmf90(A,B,C,'N','N')
	
							 
 
	! n2 = OMP_get_wtime()
	! time_gemm = time_gemm + n2-n1
	
	deallocate(A1)
	deallocate(B1)
 end subroutine gemm_omp

 
 subroutine gemmHN_omp(A,B,C,m,n,k)
 	! ! use lapack95
	! ! use blas95
	implicit none 
	
	integer m,n,k
	complex(kind=8)::A(n,m),B(n,k),C(m,k)
	complex(kind=8),allocatable::A1(:,:),B1(:,:)	
	complex(kind=8)::ctemp
	real*8::n1,n2
	integer ii,jj,kk
	allocate(A1(n,m))
	allocate(B1(n,k))	
	
	! n1 = OMP_get_wtime()	
	A1 = A
	B1 = B
	
	!$omp parallel do default(shared) private(ii,jj,kk,ctemp)	
	do jj =1,k	
	do ii =1,m
		   
	ctemp = 0
	do kk=1,n
		ctemp = ctemp + conjg(A1(kk,ii))*B1(kk,jj)	
	end do
	C(ii,jj) = ctemp
	end do
	end do
	!$omp end parallel do
	
	! ! call gemmf90(A,B,C,'N','N')
	
	! n2 = OMP_get_wtime()
	! time_gemm = time_gemm + n2-n1
	
	deallocate(A1)
	deallocate(B1)	
 end subroutine gemmHN_omp 
 

 subroutine gemmNH_omp(A,B,C,m,n,k)
 	! ! use lapack95
	! ! use blas95
	implicit none 
	
	integer m,n,k
	complex(kind=8)::A(m,n),B(k,n),C(m,k)
	complex(kind=8),allocatable::A1(:,:),B1(:,:)	
	complex(kind=8)::ctemp
	real*8::n1,n2
	integer ii,jj,kk
	
	allocate(A1(m,n))
	allocate(B1(k,n))
	
	! n1 = OMP_get_wtime()	
	A1 = A
	B1 = B
	!$omp parallel do default(shared) private(ii,jj,kk,ctemp)	
	do jj =1,k
	do ii =1,m
		   
	ctemp = 0
	do kk=1,n
		ctemp = ctemp + A1(ii,kk)*conjg(B1(jj,kk))	
	end do
	C(ii,jj) = ctemp
	end do
	end do
	!$omp end parallel do
	
	! ! call gemmf90(A,B,C,'N','N')
	
	! n2 = OMP_get_wtime()
	! time_gemm = time_gemm + n2-n1
	deallocate(A1)
	deallocate(B1)		
 end subroutine gemmNH_omp 
 
 subroutine gemmTN_omp(A,B,C,m,n,k)
 	! ! use lapack95
	! ! use blas95
	implicit none 
	
	integer m,n,k
	complex(kind=8)::A(n,m),B(n,k),C(m,k)
	complex(kind=8),allocatable::A1(:,:),B1(:,:)		
	complex(kind=8)::ctemp
	real*8::n1,n2
	integer ii,jj,kk

	allocate(A1(n,m))
	allocate(B1(n,k))	
	
	! n1 = OMP_get_wtime()	
	A1 = A
	B1 = B
	!$omp parallel do default(shared) private(ii,jj,kk,ctemp)	
	do jj =1,k	
	do ii =1,m
		   
	ctemp = 0
	do kk=1,n
		ctemp = ctemp + A1(kk,ii)*B1(kk,jj)	
	end do
	C(ii,jj) = ctemp
	end do
	end do
	!$omp end parallel do
	
	! ! call gemmf90(A,B,C,'N','N')
	
	! n2 = OMP_get_wtime()
	! time_gemm = time_gemm + n2-n1
	deallocate(A1)
	deallocate(B1)		
 end subroutine gemmTN_omp 
 
 

 subroutine gemmNT_omp(A,B,C,m,n,k)
 	! ! use lapack95
	! ! use blas95
	implicit none 
	
	integer m,n,k
	complex(kind=8)::A(m,n),B(k,n),C(m,k)
	complex(kind=8),allocatable::A1(:,:),B1(:,:)		
	complex(kind=8)::ctemp
	real*8::n1,n2
	integer ii,jj,kk
	allocate(A1(m,n))
	allocate(B1(k,n))	
	
	! n1 = OMP_get_wtime()	
	A1 = A
	B1 = B
	!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
	do jj =1,k	
	do ii =1,m
		   
	ctemp = 0
	do kk=1,n
		ctemp = ctemp + A1(ii,kk)*B1(jj,kk)	
	end do
	C(ii,jj) = ctemp
	end do
	end do
	!$omp end parallel do
	
	! ! call gemmf90(A,B,C,'N','N')
	
	! n2 = OMP_get_wtime()
	! time_gemm = time_gemm + n2-n1

	deallocate(A1)
	deallocate(B1)		
 end subroutine gemmNT_omp  
 
 
 
! #ifndef mymacro(x)

! #define mymacro(x) print *, "Now giving information about ", "x" ; \
                   ! call mysub( x, size(x,1), size(x,2) ) ; \
                   ! print *, "About to do function on ", "x"; \
                   ! call dofunction(x) ; \
                   ! print *, "x"," is a nice matrix! Huzzah!"
! #endif				   
 
 
  subroutine copymatN_omp(A,B,m,n)
	implicit none 
	integer m,n
	complex(kind=8)::A(:,:),B(:,:)
	character:: chara
	real*8::n1,n2
	integer ii,jj,ijind
	
	! n1 = OMP_get_wtime()
	

	! !$omp parallel do default(shared) private(ii,jj)	
		do ii =1,m
		do jj =1,n
			B(ii,jj) = A(ii,jj)
		end do
		end do
	! !$omp end parallel do	
	

	! n2 = OMP_get_wtime()
	! time_memcopy = time_memcopy + n2-n1
	
 end subroutine copymatN_omp
 
 
  subroutine copymatT_omp(A,B,m,n)
	implicit none 
	integer m,n
	complex(kind=8)::A(:,:),B(:,:) 
	real*8::n1,n2
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
	
 end subroutine copymatT_omp
 

function fnorm(A,m,n)
	implicit none 
	real*8:: fnorm 
	complex(kind=8)::A(m,n)
	integer m,n,ii,jj
	fnorm = 0
	do ii =1,m
	do jj =1,n
		fnorm = fnorm + abs(A(ii,jj))**2d0
	end do
	end do
	fnorm = sqrt(fnorm)
 end function fnorm
 


subroutine GetRank(M,N,mat,rank,eps)
! ! use lapack95
! ! use blas95
implicit none 
integer M,N,rank,mn,i,mnl,ii,jj
complex(kind=8)::mat(M,N)
real*8 eps
complex(kind=8), allocatable :: UU(:,:), VV(:,:),Atmp(:,:),A_tmp(:,:),tau(:),mat1(:,:)
real*8,allocatable :: Singular(:)
integer,allocatable :: jpvt(:)

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


call gesvd_robust(mat1,Singular,UU,VV,M,N,mn)
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
			call copymatT_omp(mat,Atmp,M,N)
		end if
				
		allocate(jpvt(mn))
		allocate(tau(mn))
		
		
		
		! RRQR
		jpvt = 0
		call geqp3f90(Atmp,jpvt,tau)
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
				if(abs(A_tmp(i,i))<1d-300)rank = i -1
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




subroutine ComputeRange(M,N,mat,rank,rrflag,eps)

! use lapack95
! use blas95
implicit none 
integer M,N,rank,mn,i,mnl,ii,jj,rrflag
complex(kind=8)::mat(:,:)
real*8 eps
complex(kind=8), allocatable :: UU(:,:), VV(:,:),Atmp(:,:),A_tmp(:,:),tau(:)
real*8,allocatable :: Singular(:)
integer,allocatable :: jpvt(:)

mn=min(M,N)
mnl=max(M,N)

if(isnan(fnorm(mat,M,N)))then
	write(*,*)'input matrix NAN in ComputeRange'
	stop
end if


	allocate(Atmp(mnl,mn))
	allocate(A_tmp(mn,mn))
	
	if(M>=N)then
		Atmp = mat
	else 
		call copymatT_omp(mat,Atmp,M,N)
	end if
			
	allocate(jpvt(mn))
	allocate(tau(mn))
	
	
	! RRQR
	jpvt = 0
	call geqp3f90(Atmp,jpvt,tau)
	if(isnan(fnorm(Atmp,mnl,mn)))then
		write(*,*)'Q or R has NAN in ComputeRange'
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
	
	call ungqrf90(Atmp,tau)
	
	if(M>=N)then
		mat = Atmp
	else 
		call copymatT_omp(Atmp,mat,N,M)
	end if	
	
	rank = mn
	
	if(rrflag==1)then
		do i=1,mn
			if (abs(A_tmp(i,i))/abs(A_tmp(1,1))/mnl<=eps) then
				rank=i
				if(abs(A_tmp(i,i))<1d-300)rank = i -1
				exit
			end if
		end do
	endif
	
	deallocate(jpvt)
	deallocate(tau)
	deallocate(Atmp)
	deallocate(A_tmp)

end subroutine ComputeRange



subroutine CheckRandomizedLR(M,N,mat,tolerance)
! use lapack95
! use blas95
implicit none 
integer M,N,rank,mn,i,j,k,rankmax_c,rankmax_r,rankmax_min,flag0,rank_new,rmax
complex(kind=8)::mat(M,N),ctemp
complex(kind=8),allocatable::mat1(:,:),mat2(:,:)								
real*8 tolerance,Smax
complex(kind=8), allocatable :: UU(:,:), VV(:,:),matrix_U(:,:), matrix_V(:,:),U_new(:,:),V_new(:,:),U_new1(:,:),V_new1(:,:),test_in(:,:),test_out1(:,:),test_out2(:,:),test_out3(:,:)
real*8,allocatable :: Singular(:)
integer,allocatable:: select_row(:), select_column(:)
complex(kind=8),allocatable::MatrixSubselection(:,:)

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

! call gemm_omp(U_new1,V_new1,mat2,M,rank_new,N)
! write(*,*)Singular(1:rank_new+2)
Smax = Singular(1)
        
allocate(U_new(M,rank_new))
allocate(V_new(rank_new,N))
		
do j=1, rank_new
	do i=1, M
		ctemp=0
		do k=1, rankmax_c
			ctemp=ctemp+matrix_U(i,k)*conjg(VV(j,k))
		enddo
		U_new(i,j)=ctemp
	enddo
enddo

		
do j=1, rank_new
	do i=1, N
		ctemp=0
		do k=1, rankmax_r
			ctemp=ctemp+conjg(UU(k,j))*matrix_V(i,k)
		enddo
		V_new(j,i)=ctemp/Singular(j)
	enddo
enddo

deallocate (matrix_U,VV)
deallocate (matrix_V,UU,Singular)

call gemm_omp(U_new,V_new,mat1,M,rank_new,N)


write(*,*)M,N
do j=1,N
do i=1,M
	write(211,*)dble(mat1(i,j))
	write(212,*)aimag(mat1(i,j))
end do
end do


! write(*,*)'F-norm residual:', fnorm(mat-mat1,M,N)/fnorm(mat,M,N),fnorm(mat-mat2,M,N)/fnorm(mat,M,N)
write(*,*)'F-norm residual:', fnorm(mat-mat1,M,N)/fnorm(mat,M,N)
allocate(test_in(N,1))
allocate(test_out1(M,1))
allocate(test_out2(M,1))
! allocate(test_out3(M,1))
do i=1,N
	test_in(i,1) = random_complex_number()
end do

call gemm_omp(mat,test_in,test_out1,M,N,1)
call gemm_omp(mat1,test_in,test_out2,M,N,1)
! call gemm_omp(mat2,test_in,test_out3,M,N,1)
! write(*,*)'testing vector error:', fnorm(test_out1-test_out2,M,1)/fnorm(test_out1,M,1),fnorm(test_out1-test_out3,M,1)/fnorm(test_out1,M,1)
write(*,*)'testing vector error:', fnorm(test_out1-test_out2,M,1)/fnorm(test_out1,M,1)






! method 2: ACA

rmax = min(M,N)
allocate(matrix_U(M,rmax))
allocate(matrix_V(rmax,N))

call ACA_CompressionFull(mat,matrix_U,matrix_V,M,N,rmax,rank_new,tolerance*0.3d0,tolerance)

call gemm_omp(matrix_U(1:M,1:rank_new),matrix_V(1:rank_new,1:N),mat2,M,rank_new,N)
write(*,*)'F-norm residual:', fnorm(mat-mat2,M,N)/fnorm(mat,M,N),' rank:',rank_new
		
		
deallocate(mat1)
deallocate(mat2)			
		
end subroutine CheckRandomizedLR




! Form the first mxn elements in the memory occupied by mut to a new mxn matrix new_mut
subroutine DirectCopy(mut,new_mut,m,n)
integer m,n
complex(kind=8):: mut(m,n),new_mut(m,n)
new_mut = mut(1:m,1:n)
end subroutine DirectCopy



  ! generate random permutation of 1:N  
  subroutine rperm(N, p)

    integer, intent(in):: N
    integer:: p(N)

    integer:: i
    integer:: k, j, ipj, itemp, m
    real*8, dimension(100) :: u

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
   pid = getpid()
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



  real*8 function random_real_number(mag)
	implicit none 
	real*8 a,b,c,d,mag
	integer seed
	
	! Uniform distribution
	call random_number(a)	
	a = (a-0.5d0)*2d0*mag
	
	random_real_number=a

	
	return
	end function random_real_number





  complex(kind=8) function random_complex_number()
	implicit none 
	real*8 a,b,c,d
	integer seed
	
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
	random_complex_number=a+junit*b
	
	
	! ! ! Normal distribution
	! ! call random_number(a)
	! ! seed = a*10000000
	! ! random_complex_number =  c8_normal_01 ( seed )
	
	return
	end function random_complex_number

	

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

  real*8, parameter :: r8_pi = 3.141592653589793D+00
  ! real*8 r8_uniform_01
  integer ( kind = 4 ) seed
  real*8 v1
  real*8 v2
  real*8 x_c
  real*8 x_r

  v1 = r8_uniform_01 ( seed )
  v2 = r8_uniform_01 ( seed )

  x_r = sqrt ( - 2.0D+00 * log ( v1 ) ) * cos ( 2.0D+00 * r8_pi * v2 )
  x_c = sqrt ( - 2.0D+00 * log ( v1 ) ) * sin ( 2.0D+00 * r8_pi * v2 )

  c8_normal_01 = cmplx ( x_r, x_c, kind = 8 )

  return
end function c8_normal_01

real*8 function r8_normal_01(seed)
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
!    Output, real*8 R8_NORMAL_01, a normally distributed
!    random value.
!
  implicit none

  real*8 r1
  real*8 r2

  real*8, parameter :: r8_pi = 3.141592653589793D+00
  ! real*8 r8_uniform_01
  integer ( kind = 4 ) seed
  real*8 x

  r1 = r8_uniform_01 ( seed )
  r2 = r8_uniform_01 ( seed )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  r8_normal_01 = x

  return
end function r8_normal_01


real*8 function r8_uniform_01(seed)

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
!    Output, real*8 R8_UNIFORM_01, a new pseudorandom variate,
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
real*8::input
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
real*8::input
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
real*8::input
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
	real*8 b(3)
	complex(kind=8) a(3),c
	complex(kind=8) ax,ay,az
	real*8 bx,by,bz
	ax=a(1);ay=a(2);az=a(3)
	bx=b(1);by=b(2);bz=b(3)
	c=ax*bx+ay*by+az*bz
	  return
end subroutine cscalar

subroutine scalar(a,b,c)
	  implicit none
	real*8 a(3),b(3),c
	real*8 ax,ay,az
	real*8 bx,by,bz
	ax=a(1);ay=a(2);az=a(3)
	bx=b(1);by=b(2);bz=b(3)
	c=ax*bx+ay*by+az*bz
	  return
end subroutine scalar


subroutine curl(a,b,c)
	  implicit none
	real*8 a(3),b(3),c(3)
	real*8 ax,ay,az
	real*8 bx,by,bz
	real*8 cx,cy,cz
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
	real*8 a(3)
	complex(kind=8) b(3),c(3)
	real*8 ax,ay,az
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

 real*8 function norm_vector(vector,n)
	  implicit none
	  integer n
	  complex(kind=8) vector(n)
	  integer i
	  real*8 sum
	  
	  sum=0.0d0
      do i=1,n
	      sum=sum+vector(i)*conjg(vector(i))
      enddo
      
      norm_vector=sum
      
      return 
end function

subroutine matrix_real_inverse(a,n)
           implicit none
		   integer i,j,k,n
		   real*8     a(n,n)
	       real*8     h,z
		   real*8,allocatable :: b(:),c(:)
		   integer,allocatable :: p(:),q(:)
           allocate(b(n),c(n),p(n),q(n))

	       do  6 i=1,n
	        p(i)=0
	        q(i)=0
	        b(i)=0
	        c(i)=0
6        continue

	       do  100 k=1,n

	        h=0
	       do  10 i=k,n
	       do  10  j=k,n
	        if(abs(a(i,j)).le.abs(h))  goto 10
	        h=a(i,j)
	        p(k)=i
	        q(k)=j
10       continue
	       if(abs(h).le.1d-35)  goto  150
	       if(p(k).eq.k) goto  40

	       do  30 j=1,n
	        z=a(p(k),j)
	        a(p(k),j)=a(k,j)
30        a(k,j)=z
40       if(q(k).eq.k) goto 60

	       do  50 i=1,n
	        z=a(i,q(k))
	        a(i,q(k))=a(i,k)
50        a(i,k)=z

60       do  80 j=1,n
	       if(j-k) 65,70,65
65       c(j)=a(j,k)
	       b(j)=-a(k,j)/h
	       goto 75
70       c(j)=(1.,0.)
	       b(j)=1./h
75       a(k,j)=(0.,0.)
	       a(j,k)=(0.,0.)
80       continue

	        do 90 i=1,n
	        do 90 j=1,n
90        a(i,j)=a(i,j)+c(i)*b(j)
100       continue

	       do 140 k=n,1,-1
	       if(p(k).eq.k)  goto  120

	       do 110 i=1,n
	        z=a(i,p(k))
	        a(i,p(k))=a(i,k)
110      a(i,k)=z
120       if(q(k).eq.k) goto 140
	        do  130 j=1,n 
	         z=a(q(k),j)
	         a(q(k),j)=a(k,j)
130        a(k,j)=z
140       continue

            deallocate(b,c,p,q)

	        return

150	       write(*,160)
160	       format(1x,'matrix is singular')
	         stop 00
end subroutine matrix_real_inverse



subroutine LeastSquare(m,n,k,A,b,x,eps_r)
	! ! use lapack95
	! ! use blas95
	implicit none
	
	integer m,n,k,mn_min
	complex(kind=8):: A(m,n),x(n,k),b(m,k)
	real*8,allocatable::Singular(:)
	complex(kind=8),allocatable:: Atmp(:,:),xtmp(:,:),tau(:),A_tmp(:,:),UU(:,:),VV(:,:),UU_h(:,:),VV_h(:,:),VV_h2(:,:),z_arb(:,:),x_arb(:,:),matrixtemp(:,:),A_tmp_rank(:,:),xtmp_rank(:,:),xtmp_rank3(:,:),A_tmp_rank2(:,:),xtmp_rank2(:,:)
	real*8:: eps_r
	integer ii,jj,i,j,flag0,rank
	integer,allocatable:: jpvt(:)

	complex(kind=8):: alpha,beta
	alpha=1d0
	beta=0d0

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
	
	
	if(fnorm(b,m,k)<1d-80)then
		write(*,*)'warning: RHS zero in least square'
		x = 0
	else 
	
	Atmp = A
	
	
	! SVD
    call gesvd_robust(Atmp,Singular,UU,VV,m,n,n)
	
	
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
		! call unmqrf90(Atmp,tau,b,'L','C')
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
		
		
		! call trsmf90(A_tmp_rank,xtmp_rank,'L','U','N','N')	
		
		! xtmp = 0
		
		! ! do ii = 1,n
		! ! do jj =1,k
			! ! xtmp(ii,jj) = random_complex_number()
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
				UU_h(ii,jj) = conjg(UU(jj,ii))	
			end do
			end do
			
			do ii=1,n
			do jj=1,rank
				VV_h(ii,jj) = conjg(VV(jj,ii))/Singular(jj)	
			end do
			end do 
			
			call gemmf90(UU_h, b, matrixtemp,'N','N',alpha,beta)  
			call gemmf90(VV_h, matrixtemp, x,'N','N',alpha,beta)  
			
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
	
	
	

subroutine GeneralInverse(m,n,A,A_inv,eps_r)
	! ! use lapack95
	! ! use blas95
	implicit none
	
	integer m,n,mn_min
	complex(kind=8):: A(m,n),A_inv(n,m)
	real*8,allocatable::Singular(:)
	complex(kind=8),allocatable:: Atmp(:,:),tau(:),UU(:,:),VV(:,:),UU_h(:,:),VV_h(:,:),matrixtemp(:,:),A_tmp_rank(:,:),xtmp_rank(:,:),xtmp_rank3(:,:),A_tmp_rank2(:,:),xtmp_rank2(:,:)
	real*8:: eps_r
	integer ii,jj,i,j,flag0,rank
	integer,allocatable:: jpvt(:) 

	complex(kind=8):: alpha,beta
	alpha=1d0
	beta=0d0
	
	allocate(Atmp(m,n))
	
	! SVD
	Atmp = A
	mn_min = min(m,n)
	allocate(Singular(mn_min))
	allocate(UU(m,mn_min))
	allocate(VV(mn_min,n))
	
	! write(*,*)fnorm(Atmp,m,n),'ggg'
	
    call gesvd_robust(Atmp,Singular,UU,VV,m,n,mn_min)
	
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

		do ii=1,rank
		do jj =1,m
			UU_h(ii,jj) = conjg(UU(jj,ii))	
		end do
		end do
		
		do ii=1,n
		do jj=1,rank
			VV_h(ii,jj) = conjg(VV(jj,ii))/Singular(jj)	
		end do
		end do 
		
		call gemmf90(VV_h, UU_h, A_inv,'N','N',alpha,beta)  
		
		deallocate(UU_h,VV_h)
	endif
	
	deallocate(Singular)		
	! deallocate(jpvt)	
	deallocate(UU)	
	deallocate(VV)	
	deallocate(Atmp)	
	
end subroutine GeneralInverse	
	

	subroutine RandomizedSVD(matRcol,matZRcol,matRrow,matZcRrow,matU,matV,Singular,rankmax_r,rankmax_c,rmax,rank,tolerance,SVD_tolerance)
	! ! use lapack95
	! ! use blas95
    use MODULE_FILE
    implicit none

    integer i, j, ii, jj, indx, rank_1, rank_2
    integer blocks, index_j, group_nn, rank_blocks
    integer group_m,group_n,size_of_groupm,size_of_groupn
    real*8 norm_Z,norm_U,norm_V,tolerance,SVD_tolerance
    integer edgefine_m,edgefine_n, level_butterfly, level_blocks
    integer edge_m, edge_n, header_m, header_n
    integer rank,rank1,rank2,rank12, ranknew, row, column, rankmax,rankmax_c,rankmax_r,rankmax_min,rmax,idxs_r,idxs_c
    complex(kind=8) value_Z,value_UV,maxvalue
    complex(kind=8) inner_U,inner_V,ctemp
    real*8 inner_UV
    integer,allocatable:: select_column(:), select_row(:)
	complex(kind=8)::matU(rankmax_r,rmax),matV(rmax,rankmax_c),matRcol(rankmax_c,rmax),matZRcol(rankmax_r,rmax),matRrow(rankmax_r,rmax),matZcRrow(rankmax_c,rmax)
	complex(kind=8),allocatable::matZRcol1(:,:),matZcRrow1(:,:),tau(:)
	complex(kind=8),allocatable::QQ1(:,:),RR1(:,:),QQ2(:,:),RR2(:,:),RrowcZRcol(:,:)
	real*8::Singular(rmax)
    complex(kind=8),allocatable:: row_R(:),column_R(:),matM(:,:),matrixtemp(:,:)
    real*8,allocatable:: norm_row_R(:),norm_column_R(:)
	
	complex(kind=8), allocatable :: RrowcQ1(:,:),RrowcQ1inv(:,:),Q2cRcol(:,:),Q2cRcolinv(:,:), QQ2tmp(:,:), RR2tmp(:,:), UUsml(:,:), VVsml(:,:),tau_Q(:),mattemp(:,:),matU1(:,:),matV1(:,:)	
	real*8, allocatable :: Singularsml(:)
	integer,allocatable:: jpvt(:)
	
	allocate(select_column(rankmax_c))
	allocate(select_row(rankmax_r))
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
	call geqp3f90(QQ1,jpvt,tau)
	!  call geqrff90(QQ1,tau)
	RR1=(0,0)
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rmax
		do i=1, j
			RR1(i,j)=QQ1(i,j)
		enddo
	enddo
	! !$omp end parallel do	
	call ungqrf90(QQ1,tau)
	
	! write(*,*)fnorm(QQ1,rankmax_r,rmax),rankmax_r,rmax,fnorm(RR1,rmax,rmax),rmax,rmax,'really'
	
	if(abs(RR1(1,1))<SafeUnderflow)then
		rank1=1
	else
		rank1 = rmax
		do i=1,rmax
			if (abs(RR1(i,i))/abs(RR1(1,1))/rankmax_r<=tolerance) then   ! be careful with the tolerance here
				rank1=i
				if(abs(RR1(i,i))<1d-300)rank1 = i -1
				exit
			end if
		end do
	endif
! write(*,*)shape(QQ2),shape(matZcRrow)
	
	QQ2 = matZcRrow
	jpvt = 0
	call geqp3f90(QQ2,jpvt,tau)
	!  call geqrff90(QQ2,tau)
	RR2=(0,0)
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rmax
		do i=1, j
			RR2(i,j)=QQ2(i,j)
		enddo
	enddo
	! !$omp end parallel do	
	call ungqrf90(QQ2,tau)
	
	if(abs(RR2(1,1))<SafeUnderflow)then
		rank2=1
	else	
		rank2 = rmax
		do i=1,rmax
			if (abs(RR2(i,i))/abs(RR2(1,1))/rankmax_c<=tolerance) then  ! be careful with the tolerance here
				rank2=i
				if(abs(RR2(i,i))<1d-300)rank2 = i -1
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
	allocate(RrowcQ1inv(rank1,rmax))
	call gemmHN_omp(matRrow,QQ1(1:rankmax_r,1:rank1),RrowcQ1,rmax,rankmax_r,rank1)
	
	call GeneralInverse(rmax,rank1,RrowcQ1,RrowcQ1inv,tolerance)
	deallocate(RrowcQ1)

	! allocate(ipiv(rmax))
	! call getrff90(RrowcQ1,ipiv)
	! call getrif90(RrowcQ1,ipiv)	
	! RrowcQ1inv = RrowcQ1
	! deallocate(ipiv)		
	! deallocate(RrowcQ1)
	

	
	
	
	allocate(Q2cRcol(rank2,rmax))
	allocate(Q2cRcolinv(rmax,rank2))
	call gemmHN_omp(QQ2(1:rankmax_c,1:rank2),matRcol,Q2cRcol,rank2,rankmax_c,rmax)
	
	call GeneralInverse(rank2,rmax,Q2cRcol,Q2cRcolinv,tolerance)
	deallocate(Q2cRcol)

	! allocate(ipiv(rmax))
	! call getrff90(Q2cRcol,ipiv)
	! call getrif90(Q2cRcol,ipiv)	
	! Q2cRcolinv = Q2cRcol
	! deallocate(ipiv)		
	! deallocate(Q2cRcol)	
	
	
	
	call gemmHN_omp(matRrow,matZRcol,RrowcZRcol,rmax,rankmax_r,rmax)	
	allocate(matrixtemp(rmax,rank2))
	allocate(matM(rank1,rank2))
	call gemm_omp(RrowcZRcol,Q2cRcolinv,matrixtemp,rmax,rmax,rank2)
	call gemm_omp(RrowcQ1inv,matrixtemp,matM,rank1,rmax,rank2)
	deallocate(matrixtemp,RrowcQ1inv,Q2cRcolinv)
	

	! allocate(matrixtemp(rankmax_r,rank2))
	! allocate(matM(rank1,rank2))
	! call gemm_omp(matZRcol,Q2cRcolinv,matrixtemp,rankmax_r,rmax,rank2)
	! call gemmHN_omp(QQ1(1:rankmax_r,1:rank1),matrixtemp,matM,rank1,rankmax_r,rank2)
	! deallocate(matrixtemp,RrowcQ1inv,Q2cRcolinv)	
	
	rank12 = min(rank1,rank2)
	allocate(UUsml(rank1,rank12),VVsml(rank12,rank2),Singularsml(rank12))
	call SVD_Truncate(matM,rank1,rank2,rank12,UUsml,VVsml,Singularsml,SVD_tolerance,rank)
	
	! write(111,*)UUsml(1:rank1,1:rank)
	! stop	
	
	call gemm_omp(QQ1(1:rankmax_r,1:rank1),UUsml(1:rank1,1:rank),matU(1:rankmax_r,1:rank),rankmax_r,rank1,rank)
	call gemmNH_omp(VVsml(1:rank,1:rank2),QQ2(1:rankmax_c,1:rank2),matV(1:rank,1:rankmax_c),rank,rank2,rankmax_c)
	Singular(1:rank) = Singularsml(1:rank)
	deallocate(UUsml,VVsml,Singularsml,matM)
	
	deallocate(select_column)
	deallocate(select_row)
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
			! matInv(ii,jj) = conjg(matrixtemp(jj,ii))
		! end do
		! end do
		! deallocate(matrixtemp)
		! call gemmNH_omp(matInv,matZcRrow,matM,rmax,rmax,rankmax_c)
	! else if(rankmax_c==rmax)then
	! call GetRank(rmax,rmax,matRcol,rank,tolerance)
	! write(*,*)rmax,rank,'ga'	
		! call GeneralInverse(rmax,rmax,matRcol,matInv,tolerance)
		! call gemm_omp(matZRcol,matInv,matM,rankmax_r,rmax,rankmax_c)
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
	
	

subroutine RandomMat(m,n,k,A,Oflag)
	! ! use lapack95
	! ! use blas95
	! use mkl_vsl_type
	! use mkl_vsl				 
	implicit none
	
	integer m,n,k,mn_min,ktmp
	complex(kind=8):: A(m,n)
	real(kind=8):: c			 
	real*8,allocatable::Singular(:),Ar(:,:),Ai(:,:)
	complex(kind=8):: ctemp
	complex(kind=8),allocatable:: UU(:,:),VV(:,:)
	integer ii,jj,kk,i,j,flag0,rank
	integer:: Oflag
	
	! type(VSL_STREAM_STATE) ::stream(2)
	integer brng, method, seedd,ierror								   
 	ktmp = k
	if(ktmp>min(m,n))then
		! write(*,*)'k is not properly set in RandomMat'
		ktmp=min(m,n)
	end if
	! call assert(m<=n,'m>n')
	
	call assert(k<=min(m,n),'k too large in RandomMat')
	


	! allocate(Ar(m,n))
	! allocate(Ai(m,n))
	! brng   = VSL_BRNG_MCG31
	! method = VSL_RNG_METHOD_UNIFORM_STD
	! call random_number(c)
	! seedd = NINT(1000*c)
	! ierror=vslnewstream( stream(1), brng,  seedd )
	! call random_number(c)
	! seedd = NINT(1000*c)
	! ierror=vslnewstream( stream(2), brng,  seedd )
	! ierror=vdrnguniform( method, stream(1), M*N, Ar, -1d0, 1d0) 
	! ierror=vdrnguniform( method, stream(2), M*N, Ai, -1d0, 1d0) 
	! ! ierror=vdrnggaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream, M*N, mati, 1d0, 1d0) 
	! ierror=vsldeletestream(stream(1))	
	! ierror=vsldeletestream(stream(2))	
	! A = Ar + junit*Ai
	! deallocate(Ar)
	! deallocate(Ai)



	do ii=1,m
	do jj=1,n
		A(ii,jj)=random_complex_number()
	end do
	end do


	
	! mn_min = min(m,n)
	! allocate(Singular(mn_min))
	! allocate(UU(m,mn_min))
	! allocate(VV(mn_min,n))	
	
	! call gesvd_robust(A,Singular,UU,VV,m,n,mn_min)
	! ! Singular = Singular + 1d0  ! this is to make sure the singular values are not the same
	! Singular = 1d0  ! this is to make sure the singular values are not the same
	
	! if(Oflag==1)then
		! ! call assert(mn_min==n,'mn_min not correct')
		! A = UU(1:m,1:mn_min)
	! else 
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
	! end if
	
	! deallocate(Singular)
	! deallocate(UU)
	! deallocate(VV)
	
	
	
end subroutine RandomMat	
	


subroutine KernelUpdate(matA,matB,matC,null,m,n,k,nulldim,up_tol,eps_r,error,rankmax)
! ! use lapack95
! ! use blas95
implicit none 
integer m,n,i,j,k,nulldim,nulldim_new,rank,ii,jj,kk
real*8::up_tol,eps_r
complex(kind=8)::matA(m,n),matC(n,k),matB(m,k),matB0(m,k),matBtmp(m,k),null(n,n),matC_old(n,k)
complex(kind=8),allocatable::matD(:,:),matDtmp(:,:),matE(:,:),UU(:,:),VV(:,:),nulltmp(:,:),nullnew(:,:),matrixtemp(:,:)
real*8,allocatable::Singular(:)
complex(kind=8)::ctemp
real*8::n1,n2
real*8 error
integer,optional:: rankmax

! n1 = OMP_get_wtime()
matB0 = matB

call gemm_omp(matA,matC,matBtmp,m,n,k)
error = fnorm(matB0-matBtmp,m,k)
! error = 0

if(nulldim>0)then				   
   call gemm_omp(matA,matC,matBtmp,m,n,k)  
   if(fnorm(matB0,m,k)*up_tol<fnorm(matB0-matBtmp,m,k))then
		call copymatN_omp(matC,matC_old,n,k)
		allocate(matD(m,nulldim))
		allocate(matDtmp(m,nulldim))
		allocate(matE(nulldim,k))
		call gemm_omp(matA,null(1:n,1:nulldim),matD,m,n,nulldim)		
		matB0 = matB0 - matBtmp
		matDtmp = matD
		! stop
		
		allocate(Singular(nulldim))
		allocate(UU(m,nulldim))
		allocate(VV(nulldim,nulldim))

		
		! write(*,*)fnorm(null,n,nulldim),fnorm(matD,m,nulldim)
		
		call gesvd_robust(matD,Singular,UU,VV,m,nulldim,nulldim)
		
		rank = nulldim
		do i=1,nulldim
			if (Singular(i)/Singular(1)<=eps_r) then
				rank=i
				if(Singular(i)<Singular(1)*eps_r/10)rank = i -1
				exit
			end if
		end do	
		
		if(present(rankmax))rank = min(rank,rankmax)
		
		allocate(matrixtemp(rank,k))

		
		! !$omp parallel do default(shared) private(ii,jj,kk,ctemp)	
		do ii =1,rank
		do jj =1,k
		ctemp = 0
		do kk=1,m
			ctemp = ctemp + conjg(UU(kk,ii))*matB0(kk,jj)	
		end do
		matrixtemp(ii,jj) = ctemp
		end do
		end do
		! !$omp end parallel do	

		! !$omp parallel do default(shared) private(ii,jj,kk,ctemp)	
		do ii =1,nulldim
		do jj =1,k
		ctemp = 0
		do kk=1,rank
			ctemp = ctemp + conjg(VV(kk,ii))*matrixtemp(kk,jj)/Singular(kk)		
		end do
		matE(ii,jj) = ctemp
		end do
		end do
		! !$omp end parallel do	
		
		
		
		
		call gemm_omp(null(1:n,1:nulldim),matE,matC,n,nulldim,k) 
		matC = matC + matC_old
		
		
		! write(*,*)fnorm(matC,n,k),fnorm(matE,nulldim,k),fnorm(matDtmp,m,nulldim),fnorm(matB,m,k),'ker'
		
		if(rank==nulldim)then
			null(1:n,1:nulldim) = 0
			nulldim = 0
		else 
			nulldim_new = nulldim-rank
			allocate(nulltmp(nulldim,nulldim_new))
			allocate(nullnew(n,nulldim_new))

			! !$omp parallel do default(shared) private(ii,kk)	
			do ii =1,nulldim
			do kk=1,nulldim_new
				nulltmp(ii,kk) =  conjg(VV(kk+rank,ii))		
			end do
			end do
			! !$omp end parallel do	
			call gemm_omp(null(1:n,1:nulldim),nulltmp,nullnew,n,nulldim,nulldim_new)
			null(1:n,1:nulldim) = 0
			nulldim = nulldim_new
			call copymatN_omp(nullnew,null(1:n,1:nulldim),n,nulldim_new)
			deallocate(nulltmp)
			deallocate(nullnew)
		end if
	
	
		
	
		! write(*,*)fnorm(null,n,nulldim),'hah',rank,nulldim
	
		deallocate(matrixtemp)
		deallocate(matD,matE,matDtmp)
		deallocate(Singular,UU,VV)
		
		call gemm_omp(matA,matC,matBtmp,m,n,k)  
		error = -fnorm(matB-matBtmp,m,k)
		! error = 0
		
		! if(error>1d-4)then
			! call GetRank(m,n,matA,rank,1d-4)
			! write(*,*)m,n,k,rank,error
			! call GetRank(m,k,matB,rank,1d-4)	
			! write(*,*)m,n,k,rank
			! stop
		! end if
   end if
   
end if

! n2 = OMP_get_wtime()
! time_kernelupdate = time_kernelupdate + n2-n1	

end subroutine KernelUpdate	


subroutine ACA_SubsetSelection(select_column,select_row,rankmax_r,rankmax_c,rank)

    use MODULE_FILE
    implicit none

    integer i, j, ii, jj, indx, rank_1, rank_2
    integer blocks, index_j, group_nn, rank_blocks
    integer group_m,group_n,size_of_groupm,size_of_groupn
    real*8 norm_Z,norm_U,norm_V,tolerance
    integer edgefine_m,edgefine_n, level_butterfly, level_blocks
    integer edge_m, edge_n, header_m, header_n
    integer rank, row, column, rankmax,rankmax_c,rankmax_r,rankmax_min
    complex(kind=8) value_Z,value_UV,maxvalue
    complex(kind=8) inner_U,inner_V,ctemp
    real*8 inner_UV
    integer select_column(rankmax_c), select_row(rankmax_r)

    complex(kind=8),allocatable:: row_R(:),column_R(:)
    real*8,allocatable:: norm_row_R(:),norm_column_R(:)
	rankmax_min = min(rankmax_r,rankmax_c)
    norm_Z=0
	select_column = 0
	select_row = 0
	
    allocate(row_R(rankmax_c),column_R(rankmax_r))
    allocate(norm_row_R(rankmax_c),norm_column_R(rankmax_r))

    select_row(1)=1

    ! !$omp parallel do default(shared) private(j,value_Z)
    do j=1, rankmax_c
        value_Z=MatrixSubselection(select_row(1),j)
        row_R(j)=value_Z
        norm_row_R(j)=value_Z*conjg(value_Z)
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
        matrixtemp_V(j,1)=row_R(j)
    enddo
    ! !$omp end parallel do

    ! !$omp parallel do default(shared) private(i,value_Z)
    do i=1,rankmax_r
        value_Z=MatrixSubselection(i,select_column(1))
        column_R(i)=value_Z
        norm_column_R(i)=value_Z*conjg(value_Z)
    enddo
    ! !$omp end parallel do

    norm_column_R(select_row(1))=0

    ! !$omp parallel do default(shared) private(i)
    do i=1,rankmax_r
        matrixtemp_U(i,1)=column_R(i)
    enddo
    ! !$omp end parallel do

    norm_U=norm_vector(column_R,rankmax_r)
    norm_V=norm_vector(row_R,rankmax_c)
    norm_Z=norm_Z+norm_U*norm_V

	! if(rankmax<2)write(*,*)'rankmax'
    select_row(2)=maxloc(norm_column_R,1)

    rank=1
	! write(*,*)sum(column_R),sum(row_R),norm_U,norm_V,'hehe'
    do while (norm_Z*ACA_tolerance_forward**2<norm_U*norm_V .and. rank<rankmax_min)

        ! !$omp parallel do default(shared) private(j,i,value_Z,value_UV)
        do j=1,rankmax_c
            value_Z=MatrixSubselection(select_row(rank+1),j)
            value_UV=0
            do i=1,rank
                value_UV=value_UV+matrixtemp_U(select_row(rank+1),i)*matrixtemp_V(j,i)
            enddo
            row_R(j)=value_Z-value_UV
            norm_row_R(j)=row_R(j)*conjg(row_R(j))
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
            matrixtemp_V(j,rank+1)=row_R(j)
        enddo
        ! !$omp end parallel do

        ! !$omp parallel do default(shared) private(i,j,value_Z,value_UV)
        do i=1,rankmax_r
            value_Z=MatrixSubselection(i,select_column(rank+1))
            value_UV=0
            do j=1,rank
                value_UV=value_UV+matrixtemp_U(i,j)*matrixtemp_V(select_column(rank+1),j)
            enddo
            column_R(i)=value_Z-value_UV
            norm_column_R(i)=column_R(i)*conjg(column_R(i))
        enddo
        ! !$omp end parallel do

        do i=1,rank+1
            norm_column_R(select_row(i))=0
        enddo

        ! !$omp parallel do default(shared) private(i)
        do i=1,rankmax_r
            matrixtemp_U(i,rank+1)=column_R(i)
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
                ctemp=matrixtemp_U(i,rank+1)*matrixtemp_U(i,j)
                inner_U=inner_U+ctemp*conjg(ctemp)
            enddo
            ! !$omp end parallel do
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_V)
            do i=1,rankmax_c
                ctemp=matrixtemp_V(i,rank+1)*matrixtemp_V(i,j)
                inner_V=inner_V+ctemp*conjg(ctemp)
            enddo
            ! !$omp end parallel do
            inner_UV=inner_UV+2*sqrt(inner_U*inner_V)
        enddo
        norm_Z=norm_Z+inner_UV+norm_U*norm_V

        rank=rank+1
        if (rank<rankmax_min) then
            select_row(rank+1)=maxloc(norm_column_R,1)
        endif

    enddo

	! write(*,*)select_row(1:rank),select_column(1:rank)
	
    deallocate(row_R,column_R)
    deallocate(norm_row_R,norm_column_R)

    return

end subroutine ACA_SubsetSelection





subroutine ACA_CompressionFull(mat,matU,matV,rankmax_r,rankmax_c,rmax,rank,tolerance,SVD_tolerance)
	! use lapack95
	! use blas95
    use MODULE_FILE
    implicit none

    integer i, j, ii, jj, indx, rank_1, rank_2
    integer blocks, index_j, group_nn, rank_blocks
    integer group_m,group_n,size_of_groupm,size_of_groupn
    real*8 norm_Z,norm_U,norm_V,tolerance,SVD_tolerance
    integer edgefine_m,edgefine_n, level_butterfly, level_blocks
    integer edge_m, edge_n, header_m, header_n
    integer rank, ranknew, row, column, rankmax,rankmax_c,rankmax_r,rankmax_min,rmax
    complex(kind=8) value_Z,value_UV,maxvalue
    complex(kind=8) inner_U,inner_V,ctemp
    real*8 inner_UV
    integer,allocatable:: select_column(:), select_row(:)
	complex(kind=8)::mat(rankmax_r,rankmax_c),matU(rankmax_r,rmax),matV(rmax,rankmax_c)
	
    complex(kind=8),allocatable:: row_R(:),column_R(:)
    real*8,allocatable:: norm_row_R(:),norm_column_R(:)
	
	complex(kind=8), allocatable :: QQ1(:,:), RR1(:,:),QQ2(:,:), RR2(:,:),QQ2tmp(:,:), RR2tmp(:,:), UUsml(:,:), VVsml(:,:),tau_Q(:),mattemp(:,:),matU1(:,:),matV1(:,:)	
	real*8, allocatable :: Singularsml(:)
	
		
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
        norm_row_R(j)=value_Z*conjg(value_Z)
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
        norm_column_R(i)=value_Z*conjg(value_Z)
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
            norm_row_R(j)=row_R(j)*conjg(row_R(j))
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
            norm_column_R(i)=column_R(i)*conjg(column_R(i))
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
                inner_U=inner_U+ctemp*conjg(ctemp)
            enddo
            ! !$omp end parallel do
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_V)
            do i=1,rankmax_c
                ctemp=matV(rank+1,i)*matV(j,i)
                inner_V=inner_V+ctemp*conjg(ctemp)
            enddo
            ! !$omp end parallel do
            inner_UV=inner_UV+2*sqrt(inner_U*inner_V)
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
	call copymatN_omp(matU(1:rankmax_r,1:rank),QQ1,rankmax_r,rank)
	allocate (tau_Q(rank))
	call geqrff90(QQ1,tau_Q)
	
	allocate (RR1(rank,rank))
	RR1=(0,0)
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rank
		do i=1, j
			RR1(i,j)=QQ1(i,j)
		enddo
	enddo
	! !$omp end parallel do	
	call ungqrf90(QQ1,tau_Q)
	deallocate(tau_Q)


	allocate(QQ2tmp(rankmax_c,rank))
	call copymatT_omp(matV(1:rank,1:rankmax_c),QQ2tmp,rank,rankmax_c)
	allocate (tau_Q(rank))
	call geqrff90(QQ2tmp,tau_Q)
	
	allocate (RR2tmp(rank,rank))
	RR2tmp=(0,0)
	! !$omp parallel do default(shared) private(i,j)
	do j=1, rank
		do i=1, j
			RR2tmp(i,j)=QQ2tmp(i,j)
		enddo
	enddo
	! !$omp end parallel do	
	call ungqrf90(QQ2tmp,tau_Q)
	deallocate(tau_Q)
	
	allocate(QQ2(rank,rankmax_c))
	call copymatT_omp(QQ2tmp,QQ2,rankmax_c,rank)	
	allocate(RR2(rank,rank))
	call copymatT_omp(RR2tmp,RR2,rank,rank)	
	
	
	
	! allocate(matU1(rankmax_r,rank))
	! allocate(matV1(rank,rankmax_c))
	! call gemm_omp(QQ1,RR1,matU1,rankmax_r,rank,rank)
	
	! call gemm_omp(RR2,QQ2,matV1,rank,rank,rankmax_c)
	
	! write(*,*)fnorm(matU1-matU(1:rankmax_r,1:rank),rankmax_r,rank),fnorm(matV1-matV(1:rank,1:rankmax_c),rank,rankmax_c)
	
	
	
	deallocate(QQ2tmp,RR2tmp)
	allocate(mattemp(rank,rank))
	call gemm_omp(RR1,RR2,mattemp,rank,rank,rank)
	
	
	
	
	allocate(UUsml(rank,rank),VVsml(rank,rank),Singularsml(rank))
	call SVD_Truncate(mattemp,rank,rank,rank,UUsml,VVsml,Singularsml,SVD_tolerance,ranknew)
	
	call gemm_omp(QQ1,UUsml(1:rank,1:ranknew),matU(1:rankmax_r,1:ranknew),rankmax_r,rank,ranknew)
	
	call gemm_omp(VVsml(1:ranknew,1:rank),QQ2,matV(1:ranknew,1:rankmax_c),ranknew,rank,rankmax_c)
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


subroutine SVD_Truncate(mat,mm,nn,mn,UU,VV,Singular,tolerance,rank)
! ! use lapack95
! ! use blas95
implicit none 
integer mm,nn,mn,rank,ii,jj
real*8:: tolerance
complex(kind=8)::mat(mm,nn),UU(mm,mn),VV(mn,nn)
complex(kind=8),allocatable::mat0(:,:)					 
real*8:: Singular(mn)
integer::i,flag

allocate(mat0(mm,nn))
mat0 = mat
call gesvd_robust(mat0,Singular,UU,VV,mm,nn,mn)
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



! sort the first dimension of the array
subroutine PIKSRT_DBLE_Multi(N,M,ARR)
  implicit none
  integer j,i,N,M
  real*8 ARR(N,M)
  real*8,allocatable::a(:)
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


subroutine quick_sort(list, order, n)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.

IMPLICIT NONE
integer n
real*8 list(n)
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
real*8                :: reference, temp
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
real*8                :: temp

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
real*8,intent(in)::xin,yin,zin,origin(3)
real*8,intent(out)::r,theta,phi
real*8:: x,y,z
	
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


subroutine gesvd_robust(Matrix,Singular,UU,VV,mm,nn,mn_min)
    use MODULE_FILE
	implicit none 
	integer mm,nn,mn_min
	complex(kind=8)Matrix(mm,nn),UU(mm,mn_min),VV(mn_min,nn)
	real*8 Singular(mn_min)

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
			write(*,*)'ga',fnorm(Matrix,mm,nn),shape(Matrix),shape(UU),shape(VV)
			call gesvdf90(Matrix,Singular,UU,VV)
			write(*,*)'gani'
		endif
	endif
	
	
end subroutine gesvd_robust
		
	
subroutine gesvdf90(Matrix,Singular,UU,VV)	

! ! use lapack95

implicit none
complex(kind=8)Matrix(:,:),UU(:,:),VV(:,:)
real*8 Singular(:)
integer m,n,mn_min

integer LWORK,INFO
complex(kind=8):: TEMP(1)
real*8,allocatable::RWORK(:)

complex(kind=8),allocatable:: WORK(:)

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)


allocate(RWORK(5*m))
LWORK=-1
call ZGESVD('S','S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, TEMP, LWORK, RWORK, INFO)

LWORK=TEMP(1)*1.001
allocate(WORK(LWORK))     

call ZGESVD('S','S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, WORK, LWORK, RWORK, INFO)

if(INFO/=0)then
	write(*,*)'SVD failed!!',INFO
	stop
endif

deallocate(WORK,RWORK)

end subroutine gesvdf90	


subroutine geqrff90(Matrix,tau)	

! ! use lapack95

implicit none
complex(kind=8)Matrix(:,:),tau(:)
integer m,n,mn_min

integer LWORK,INFO
complex(kind=8):: TEMP(1)

complex(kind=8),allocatable:: WORK(:)

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)


LWORK=-1
call ZGEQRF(m, n, Matrix, m, tau, TEMP, LWORK, INFO)

LWORK=TEMP(1)*1.001
allocate(WORK(LWORK))     

call ZGEQRF(m, n, Matrix, m, tau, WORK, LWORK, INFO)


if(INFO/=0)then
	write(*,*)'geqrff90 failed!!',INFO
	stop
endif

deallocate(WORK)


! call geqrf(Matrix,tau)

end subroutine geqrff90


subroutine geqp3f90(Matrix,jpvt,tau)	

! ! use lapack95

implicit none
complex(kind=8)Matrix(:,:),tau(:)
integer jpvt(:)
integer m,n,mn_min

integer LWORK,INFO
complex(kind=8):: TEMP(1)
real*8,allocatable::RWORK(:)

complex(kind=8),allocatable:: WORK(:)

m=size(Matrix,1)
n=size(Matrix,2)
mn_min = min(m,n)

allocate(RWORK(2*m))

LWORK=-1
call ZGEQP3(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, RWORK, INFO)
LWORK=TEMP(1)*1.001
allocate(WORK(LWORK))     
call ZGEQP3(m, n, Matrix, m, jpvt, tau, WORK, LWORK, RWORK, INFO)

if(INFO/=0)then
	write(*,*)'geqp3f90 failed!!',INFO
	stop
endif

deallocate(WORK)
deallocate(RWORK)

! call geqp3(Matrix,jpvt,tau)

end subroutine geqp3f90
	

subroutine unmqrf90(a,tau,c,side,trans)	

! ! use lapack95

implicit none
complex(kind=8)a(:,:),tau(:),c(:,:)
integer m,n,k,lda,ldc, mn_min

character side, trans

integer LWORK,INFO
complex(kind=8):: TEMP(1)
complex(kind=8),allocatable:: WORK(:)

m=size(c,1)
n=size(c,2)
mn_min = min(m,n)
ldc=m

if(side == 'L')then
	k=m
	lda=m
end if

if(side == 'R')then
	k=n
	lda=n
end if


LWORK=-1
call ZUNMQR(side, trans, m, n, k, a, lda, tau, c, ldc, TEMP, LWORK, INFO)
LWORK=TEMP(1)*1.001
allocate(WORK(LWORK))     
call ZUNMQR(side, trans, m, n, k, a, lda, tau, c, ldc, WORK, LWORK, INFO)

if(INFO/=0)then
	write(*,*)'unmqrf90 failed!!',INFO
	stop
endif

deallocate(WORK)

! call unmqr(a,tau,c,side,trans)

end subroutine unmqrf90

	
	
subroutine ungqrf90(Matrix,tau)	

! ! use lapack95

implicit none
complex(kind=8)Matrix(:,:),tau(:)
integer m,n,mn_min

integer LWORK,INFO
complex(kind=8):: TEMP(1)

complex(kind=8),allocatable:: WORK(:)

m=size(Matrix,1)
n=size(Matrix,2)
call assert(m>=n,'size of Matrix incorrect')

mn_min = min(m,n)


LWORK=-1
call ZUNGQR(m, n, n, Matrix, m, tau, TEMP, LWORK, INFO)

LWORK=TEMP(1)*1.001
allocate(WORK(LWORK))     

call ZUNGQR(m, n, n, Matrix, m, tau, WORK, LWORK, INFO)


if(INFO/=0)then
	write(*,*)'ungqrf90 failed!!',INFO
	stop
endif

deallocate(WORK)


! call ungqr(Matrix,tau)

end subroutine ungqrf90	
	
	

subroutine getrff90(Matrix,ipiv)	

! ! use lapack95

implicit none
complex(kind=8) Matrix(:,:)
integer ipiv(:)
integer m,n,mn_min

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
	
end subroutine getrff90	
	
	
	
subroutine getrsf90(Matrix,ipiv,B,trans)	

! ! use lapack95

implicit none
complex(kind=8) Matrix(:,:),B(:,:)
integer ipiv(:)
integer m,n,mn_min,nrhs
character trans

integer INFO

n=size(Matrix,1)
nrhs=size(B,2)


call ZGETRS( trans, n, nrhs, Matrix, n, ipiv, B, n, INFO )


if(INFO/=0)then
	write(*,*)'getrsf90 failed!!',INFO
	stop
endif
	
	
! call getrs(Matrix,ipiv,B,trans)	
	
end subroutine getrsf90		
	
	
	
	
subroutine getrif90(Matrix,ipiv)	

! ! use lapack95

implicit none
complex(kind=8)Matrix(:,:)
integer ipiv(:)
integer m,n,mn_min

integer LWORK,INFO
complex(kind=8):: TEMP(1)

complex(kind=8),allocatable:: WORK(:)

m=size(Matrix,1)
n=size(Matrix,2)
call assert(m<=n,'size of Matrix incorrect')

mn_min = min(m,n)


LWORK=-1
call ZGETRI(m, Matrix, m, ipiv, TEMP, LWORK, INFO)

LWORK=TEMP(1)*1.001
allocate(WORK(LWORK))     

call ZGETRI(m, Matrix, m, ipiv, WORK, LWORK, INFO)


if(INFO/=0)then
	write(*,*)'getrif90 failed!!',INFO
	stop
endif

deallocate(WORK)


! call getri(Matrix,ipiv)	

end subroutine getrif90		
	


subroutine trsmf90(Matrix,B,side,uplo,transa,diag)	

! ! use blas95

implicit none
complex(kind=8) Matrix(:,:),B(:,:),alpha
integer m,n,lda,ldb
character side,uplo,transa,diag

integer INFO

alpha=1d0

m=size(B,1)
n=size(B,2)

ldb=m

if(side=='L' .or. side=='l')lda=m
if(side=='R' .or. side=='r')lda=n

call ZTRSM(side, uplo, transa, diag, m, n, alpha, Matrix, lda, B, ldb)


! call trsm(Matrix,B,side,uplo,transa,diag)	
	
end subroutine trsmf90		



subroutine gemmf90(MatA,MatB,MatC,transa,transb,alpha,beta)	

! ! use blas95

implicit none
complex(kind=8) MatA(:,:),MatB(:,:),MatC(:,:),alpha,beta
integer m,n,k,lda,ldb,ldc
character transa,transb

m=size(MatC,1)
n=size(MatC,2)
ldc=m

if(transa=='N' .or. transa=='n')then
k=size(MatA,2)
lda=m
endif

if(transa=='T' .or. transa=='t')then
k=size(MatA,1)
lda=k
endif


if(transb=='N' .or. transb=='n')then
ldb=k
endif

if(transb=='T' .or. transb=='t')then
ldb=n
endif

call zgemm(transa, transb, m, n, k, alpha, MatA, lda, MatB, ldb, beta, MatC, ldc)

! call gemmf90(MatA,MatB,MatC,transa,transb,alpha,beta)	

end subroutine gemmf90

end module misc
