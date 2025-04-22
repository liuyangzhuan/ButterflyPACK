!> @file
!> @brief This example tests the functionality of QTT-SVD compression and contraction when applied to a 2D/3D IE operator
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
   integer :: NDIM
   real*8 :: tol
   real*8 :: wavelen, ppw, ds, waven
   real*8 :: zdist
   integer :: s, d, num_twos, n_i,brand_size,i
   integer :: io, jo, ko, is, js, ks
   real*8 :: r
   real*8, dimension(3) :: oo, ss
   integer :: idx, Nmax,ZFPflag
   integer, allocatable:: N(:)

   complex*16, allocatable :: tensor_in4(:,:,:,:),tensor_in6(:,:,:,:,:,:)

   ! Flattened or reshaped arrays
   integer, allocatable :: shape_hd(:)
   integer, allocatable :: perm_o1(:),perm_o2(:),perm_o3(:),perm_s1(:),perm_s2(:),perm_s3(:)
   integer, allocatable :: shapes(:)
   integer, allocatable :: ms(:),ns(:)
   complex*16, allocatable :: tensor_hd(:), tensor_hd1(:), matrix(:,:)
   integer :: tot_size

   type(TTtype) :: qtt         ! The TT
   type(TTtype) :: tmat        ! Interpreted as TT-matrix
   complex*16, allocatable :: brand(:)
   type(TTtype) :: tb_tt, qtc1, qtc
   complex*16, allocatable :: tb(:), tc(:), tc1(:)
   real*8,allocatable:: real_data(:), imag_data(:)

   !=================================================================
   ! Setup parameters
   !=================================================================
   fullvec = 1
   tol = 1d-4
   NDIM = 3
   ZFPflag = 1

   wavelen = 0.0625
   ppw     = 4.0d0
   ds      = wavelen / ppw
   waven   = 2d0*3.141592653589793d0 / wavelen
   allocate(N(NDIM*2))

   if(NDIM==2)then
      zdist   = 1d0
      N = [5,5,6,6]
      Nmax = next_power_of_2(maxval(N))
      allocate(perm_o1(Nmax))
      call rperm(Nmax, perm_o1)
      allocate(perm_o2(Nmax))
      call rperm(Nmax, perm_o2)
      allocate(perm_s1(Nmax))
      call rperm(Nmax, perm_s1)
      allocate(perm_s2(Nmax))
      call rperm(Nmax, perm_s2)

      !=================================================================
      ! Build the 4D operator (N x N x N x N)
      !    complex kernel: exp(-i*waven*r)/r
      !=================================================================
      allocate(tensor_in4(N(1),N(2),N(3),N(4)))
      tensor_in4 = (0d0, 0d0)  ! init to (0 + 0i)

      do io=1,N(1)
         do jo=1,N(2)
            oo = (/ perm_o1(io)*ds, perm_o2(jo)*ds, 0d0 /)
            do is=1,N(3)
               do js=1,N(4)
                  ss = (/ perm_s1(is)*ds, perm_s2(js)*ds, zdist /)
                  r = sqrt( (oo(1)-ss(1))**2 + (oo(2)-ss(2))**2 + (oo(3)-ss(3))**2 )
                  if (r /= 0d0) then
                     tensor_in4(io,jo,is,js) = exp( -dcmplx(0d0,1d0)*0*r ) / r
                  else
                     tensor_in4(io,jo,is,js) = dcmplx(0d0,0d0)
                  end if
               end do
            end do
         end do
      end do
      ! Flatten "tensor_in" => 1D array "tensor_hd"
      allocate(tensor_hd(product(N)))
      tensor_hd = reshape(tensor_in4, shape(tensor_hd))
      qtt%d_org = NDIM
   else if (NDIM==3)then
      zdist   = 2d0
      ! N = [5,5,6,6,6,4]
      N = [3,5,5,3,5,5]
      ! Nmax = next_power_of_2(maxval(N))
      Nmax = 16
      ! N = Nmax
      allocate(perm_o1(Nmax))
      allocate(perm_o2(Nmax))
      allocate(perm_o3(Nmax))
      allocate(perm_s1(Nmax))
      allocate(perm_s2(Nmax))
      allocate(perm_s3(Nmax))

      ! !!!!!!!!!!!! randomly choose N(i) indices out of Nmax
      ! call rperm(Nmax, perm_o1)
      ! call rperm(Nmax, perm_o2)
      ! call rperm(Nmax, perm_o3)
      ! call rperm(Nmax, perm_s1)
      ! call rperm(Nmax, perm_s2)
      ! call rperm(Nmax, perm_s3)

      !!!!!!!!!!! N(i) indices by linspace(1,Nmax,N(i))
      call linspaceI(1, Nmax, N(1), perm_o1(1:N(1)))
      call linspaceI(1, Nmax, N(2), perm_o2(1:N(2)))
      call linspaceI(1, Nmax, N(3), perm_o3(1:N(3)))
      call linspaceI(1, Nmax, N(4), perm_s1(1:N(4)))
      call linspaceI(1, Nmax, N(5), perm_s2(1:N(5)))
      call linspaceI(1, Nmax, N(6), perm_s3(1:N(6)))

      ! !!!!!!!!!!! The first N(i) indices of 1:Nmax
      ! call linspaceI(1, Nmax, Nmax, perm_o1)
      ! call linspaceI(1, Nmax, Nmax, perm_o2)
      ! call linspaceI(1, Nmax, Nmax, perm_o3)
      ! call linspaceI(1, Nmax, Nmax, perm_s1)
      ! call linspaceI(1, Nmax, Nmax, perm_s2)
      ! call linspaceI(1, Nmax, Nmax, perm_s3)


      !=================================================================
      ! Build the 6D operator
      !    complex kernel: exp(-i*waven*r)/r
      !=================================================================
      allocate(tensor_in6(N(1),N(2),N(3),N(4),N(5),N(6)))
      tensor_in6 = (0d0, 0d0)  ! init to (0 + 0i)

      do io=1,N(1)
         do jo=1,N(2)
            do ko=1,N(3)
               oo = (/ perm_o1(io)*ds, perm_o2(jo)*ds, perm_o3(ko)*ds /)
               do is=1,N(4)
                  do js=1,N(5)
                     do ks=1,N(6)
                        ss = (/ perm_s1(is)*ds, perm_s2(js)*ds, perm_s3(ks)*ds + zdist /)
                        r = sqrt( (oo(1)-ss(1))**2 + (oo(2)-ss(2))**2 + (oo(3)-ss(3))**2 )
                        if (r /= 0d0) then
                           tensor_in6(io,jo,ko,is,js,ks) = exp( -dcmplx(0d0,1d0)*waven*r ) / r
                        else
                           tensor_in6(io,jo,ko,is,js,ks) = dcmplx(0d0,0d0)
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do

      ! Flatten "tensor_in" => 1D array "tensor_hd"
      allocate(tensor_hd(product(N)))
      tensor_hd = reshape(tensor_in6, shape(tensor_hd))

      ! allocate(real_data(product(N)))
      ! allocate(imag_data(product(N)))
      ! read(777,*)real_data,imag_data
      ! tensor_hd = real_data + BPACK_junit*imag_data


      deallocate(perm_o1)
      deallocate(perm_o2)
      deallocate(perm_o3)
      deallocate(perm_s1)
      deallocate(perm_s2)
      deallocate(perm_s3)


      qtt%d_org = NDIM
   endif


   qtt%mpo = 1
   allocate(qtt%m_n_org(qtt%d_org,2))
   qtt%m_n_org(:,1)=N(1:NDIM)
   qtt%m_n_org(:,2)=N(NDIM+1:NDIM*2)

   call QTT_Compress_SVD(tensor_hd, tol, qtt, ZFPflag)
   call QTT_Decompress(qtt,tensor_hd1)
   print *,"Original tensor dimensions:",qtt%m_n_org(:,1)," by ",qtt%m_n_org(:,2)
   print *,"Ranks in qtt:",qtt%r
   write(*,*)'QTT relative error: ', sqrt(sum(abs(tensor_hd-tensor_hd1)**2d0))/sqrt(sum(abs(tensor_hd)**2d0))
   if(allocated(qtt%core))then
      write(*,*)'QTT compression ratio: ', product(N)/dble(product(shape(qtt%core)))
   else
      write(*,*)'QTT compression ratio: ', product(N)*16/dble((SIZEOF(qtt%coreZFP%buffer_r)+SIZEOF(qtt%coreZFP%buffer_i)))
   endif

   allocate(tb(product(qtt%m_n_org(:,2))))
   tb = 1.d0
   allocate(tc(product(qtt%m_n_org(:,1))))  ! result
   call QTT_Apply_Fullvec(qtt, tb, tc)

   allocate(matrix(product(qtt%m_n_org(:,1)),product(qtt%m_n_org(:,2))))
   matrix = reshape(tensor_hd1,[product(qtt%m_n_org(:,1)),product(qtt%m_n_org(:,2))])
   allocate(tc1(product(qtt%m_n_org(:,1))))  ! result
   call gemmf90(matrix, product(qtt%m_n_org(:,1)), tb, product(qtt%m_n_org(:,2)), tc1, product(qtt%m_n_org(:,1)), 'N', 'N', product(qtt%m_n_org(:,1)), 1, product(qtt%m_n_org(:,2)), BPACK_cone, BPACK_czero)

   write(*,*)'fnorm of output:', sqrt(sum(abs(tc)**2d0)), sqrt(sum(abs(tc1)**2d0)), sqrt(sum(abs(tc-tc1)**2d0))
   call TT_Delete(qtt)

   deallocate(tb)
   deallocate(tc)
   deallocate(tensor_hd)
   deallocate(tensor_hd1)

   if(NDIM==2)then
      deallocate(tensor_in4)
   else if (NDIM==3)then
      deallocate(tensor_in6)
   endif

end program main

