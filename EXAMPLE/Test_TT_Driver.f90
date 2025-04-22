!> @file
!> @brief This example tests the functionality of QTT-SVD compression and contraction when applied to a 2D IE operator
!> @details Note that instead of the use of precision dependent subroutine/module/type names "z_", one can also use the following \n
!> #define DAT 0 \n
!> #include "zButterflyPACK_config.fi" \n
!> which will macro replace precision-independent subroutine/module/type names "X" with "z_X" defined in SRC_DOUBLECOMLEX with double-complex precision

#define DAT 0
#include "zButterflyPACK_config.fi"



!====================================================================
program main
  use MISC_Utilities
  implicit none



  ! Local declarations
  integer :: fullvec
  real*8 :: tol
  real*8 :: wavelen, ppw, ds, waven
  real*8 :: zdist
  integer :: N, s, d, num_legs, n_i,brand_size,i
  integer :: io, jo, is, js
  real*8 :: r
  real*8, dimension(3) :: oo, ss
  integer :: idx

  complex*16, allocatable :: tensor_in(:,:,:,:)

  ! Flattened or reshaped arrays
  integer, allocatable :: shape_hd(:)
  integer, allocatable :: perm(:)
  integer, allocatable :: shapes(:)
  integer, allocatable :: ms(:),ns(:)
  complex*16, allocatable :: tensor_hd(:), tensor_hd1(:)
  integer :: tot_size

  type(TTtype) :: tt         ! The TT
  type(TTtype) :: tmat        ! Interpreted as TT-matrix
  complex*16, allocatable :: brand(:), matrix(:,:)
  type(TTtype) :: tb_tt, qtc1, qtc
  complex*16, allocatable :: tb(:), tc(:), tc1(:)

  !=================================================================
  ! Setup parameters
  !=================================================================
  fullvec = 1
  tol = 1d-5

  wavelen = 0.0078125d0
  ppw     = 4.0d0
  ds      = wavelen / ppw
  waven   = 2d0*3.141592653589793d0 / wavelen
  zdist   = 1d0

  N = 8
  d = 2
  n_i=N*N
  num_legs = d

  !=================================================================
  ! Build the 4D operator (N x N x N x N)
  !    complex kernel: exp(-i*waven*r)/r
  !=================================================================
  allocate(tensor_in(N,N,N,N))
  tensor_in = (0d0, 0d0)  ! init to (0 + 0i)

  do io=1,N
     do jo=1,N
        oo = (/ io*ds, jo*ds, 0d0 /)
        do is=1,N
           do js=1,N
              ss = (/ is*ds, js*ds, 0d0 /)
              r = sqrt( (oo(1)-ss(1))**2 + (oo(2)-ss(2))**2 + (oo(3)-ss(3))**2 )
              if (r /= 0d0) then
                 tensor_in(io,jo,is,js) = exp( -dcmplx(0d0,1d0)*waven*r ) / r
              else
                 tensor_in(io,jo,is,js) = dcmplx(0d0,0d0)
              end if
           end do
        end do
     end do
  end do


  ! Flatten "tensor_in" => 1D array "tensor_hd"
  allocate(tensor_hd(N*N*N*N))
  tensor_hd = reshape(tensor_in, shape(tensor_hd))


  allocate(shapes(num_legs))
  shapes = N*N


  ! Build  TT decomposition
  tt%d = num_legs
  tt%mpo = 1
  allocate(tt%n(num_legs))
  tt%n = shapes
  shapes = N
  allocate(tt%m_n(num_legs,2))
  tt%m_n(:,1)=shapes
  tt%m_n(:,2)=shapes

  deallocate(shapes)
!   write(*,*)real(tensor_hd(1:8)),'good so far'
  call TT_Compress_SVD(tensor_hd, tol, tt)
  call TT_Decompress(tt,tensor_hd1)
  print *,"Ranks in tt:",tt%r

     write(*,*)'TT relative error: ', sqrt(sum(abs(tensor_hd-tensor_hd1)**2d0))/sqrt(sum(abs(tensor_hd)**2d0))
     write(*,*)'TT compression ratio: ', N*N*N*N/dble(product(shape(tt%core)))


  ! Multiply by either a full vector or a TT-vector
  if (fullvec == 1) then
     allocate(tb(n_i))
     tb = 1.d0

     allocate(tc(n_i))  ! result
     call TT_Apply_Fullvec(tt, tb, tc)

      allocate(matrix(product(tt%m_n(:,1)),product(tt%m_n(:,2))))
      matrix = reshape(tensor_hd1,[product(tt%m_n(:,1)),product(tt%m_n(:,2))])
      allocate(tc1(product(tt%m_n(:,1))))  ! result
      tc1=0
      call gemmf90(matrix, product(tt%m_n(:,1)), tb, product(tt%m_n(:,2)), tc1, product(tt%m_n(:,1)), 'N', 'N', product(tt%m_n(:,1)), 1, product(tt%m_n(:,2)), BPACK_cone, BPACK_czero)


     write(*,*)'fnorm of output:', sqrt(sum(abs(tc)**2d0)),sqrt(sum(abs(tc1)**2d0)),sqrt(sum(abs(tc-tc1)**2d0))





     qtc%mpo = 0
     qtc%d = num_legs
     allocate(qtc%n(num_legs))
     do i=1,num_legs
       qtc%n(i) = N
     end do


     call TT_Compress_SVD(tc, 1d-30, qtc)
     print *,"Ranks in qtc:", qtc%r


     call TT_Delete(tt)
     call TT_Delete(qtc)
     deallocate(tb)
     deallocate(tc)
     deallocate(tensor_hd)
     deallocate(tensor_in)


  else

     brand_size = N*N
     allocate(brand(brand_size))
     brand = 1.d0

     tb_tt%mpo = 0
     tb_tt%d = tt%d
     allocate(tb_tt%n(tb_tt%d))
     do i=1,tb_tt%d
       tb_tt%n(i) = N
     end do
     call TT_Compress_SVD(brand, tol, tb_tt)
    !  print *,"Ranks in tb_tt:"
    !  do i=1,tb_tt%d+1
    !    print *, tb_tt%r(i)
    !  end do

     call TT_Apply_TTvec(tt, tb_tt, qtc1)

     call TT_Decompress(qtc1, tc)
     write(*,*)'fnorm of output: ', sqrt(sum(abs(tc)**2d0))

     print *,"Ranks in qtc1:",qtc1%r

     call TT_Delete(tt)
     call TT_Delete(tb_tt)
     call TT_Delete(qtc1)
     deallocate(brand)
     deallocate(tc)
     deallocate(tensor_hd)
     deallocate(tensor_in)

  end if


end program main

