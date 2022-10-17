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

!> @file MISC_DenseLA.f90
!> @brief Low-level wrappers for calling BLAS/LAPACK/ScaLAPACK functions


#include "ButterflyPACK_config.fi"
module MISC_DenseLA
   use BPACK_DEFS
#ifdef HAVE_OPENMP
   use omp_lib
#endif

contains

   DTR function fnorm(Matrix, M, N, norm)
      DT Matrix(:, :)
      integer M, N
      character, optional:: norm
#if DAT==0
      fnorm = zlangef90(Matrix, M, N, norm)
#elif DAT==1
      fnorm = dlangef90(Matrix, M, N, norm)
#elif DAT==2
      fnorm = clangef90(Matrix, M, N, norm)
#elif DAT==3
      fnorm = slangef90(Matrix, M, N, norm)
#endif

   end function fnorm

   integer function numroc_wp(n, nb, iproc, isrcproc, nprocs)
      integer :: n, nb, iproc, isrcproc, nprocs
      integer :: numroc ! blacs routine
      numroc_wp = numroc(n, nb, iproc, isrcproc, nprocs)
   end function numroc_wp

   pure subroutine g2l(i, n, np, nb, p, il)
      implicit none
      integer,intent(in) :: i    ! global array index, input
      integer,intent(in) :: n    ! global array dimension, input
      integer,intent(in) :: np   ! processor array dimension, input
      integer,intent(in) :: nb   ! block size, input
      integer,intent(out) :: p    ! processor array index, output
      integer,intent(out) :: il   ! local array index, output
      integer :: im1
      im1 = i - 1
      p = mod((im1/nb), np)
      il = (im1/(np*nb))*nb + mod(im1, nb) + 1
      return
   end subroutine g2l

! convert local index to global index in block-cyclic distribution
   pure subroutine l2g(il, p, n, np, nb, i)
      implicit none
      integer,intent(in) :: il   ! local array index, input
      integer,intent(in) :: p    ! processor array index, input
      integer,intent(in) :: n    ! global array dimension, input
      integer,intent(in) :: np   ! processor array dimension, input
      integer,intent(in) :: nb   ! block size, input
      integer,intent(out) :: i    ! global array index, output
      integer :: ilm1
      ilm1 = il - 1
      i = (((ilm1/nb)*np) + p)*nb + mod(ilm1, nb) + 1
      return
   end subroutine l2g

   real(kind=8) function dlangef90(Matrix, M, N, norm)

!

      implicit none
      real(kind=8) Matrix(:, :)
      integer M, N, lda
      character, optional:: norm
      character::opt
      real(kind=8), allocatable:: WORK(:)
      real(kind=8) dlange

      opt = 'F'
      if (present(norm)) opt = norm
      if (opt == 'I' .or. opt == 'i') then
         allocate (WORK(M))
         WORK = 0
      end if

      lda = size(Matrix, 1)
      dlangef90 = dlange(opt, M, N, Matrix, lda, WORK)

      if (opt == 'I' .or. opt == 'i') then
         deallocate (WORK)
      end if

   end function dlangef90


   real(kind=4) function slangef90(Matrix, M, N, norm)

!

      implicit none
      real(kind=4) Matrix(:, :)
      integer M, N, lda
      character, optional:: norm
      character::opt
      real(kind=4), allocatable:: WORK(:)
      real(kind=4) slange

      opt = 'F'
      if (present(norm)) opt = norm
      if (opt == 'I' .or. opt == 'i') then
         allocate (WORK(M))
         WORK = 0
      end if

      lda = size(Matrix, 1)
      slangef90 = slange(opt, M, N, Matrix, lda, WORK)

      if (opt == 'I' .or. opt == 'i') then
         deallocate (WORK)
      end if

   end function slangef90


   real(kind=8) function zlangef90(Matrix, M, N, norm)
!
      implicit none
      complex(kind=8) Matrix(:, :)
      integer M, N, lda
      character, optional:: norm
      character::opt
      complex(kind=8), allocatable:: WORK(:)
      real(kind=8) zlange

      opt = 'F'
      if (present(norm)) opt = norm
      if (opt == 'I' .or. opt == 'i') then
         allocate (WORK(M))
         WORK = 0
      else
         allocate(WORK(1))
      end if

      lda = size(Matrix, 1)
      zlangef90 = zlange(opt, M, N, Matrix, lda, WORK)

      ! if (opt == 'I' .or. opt == 'i') then
         deallocate (WORK)
      ! end if

   end function zlangef90


   real(kind=4) function clangef90(Matrix, M, N, norm)
!
      implicit none
      complex(kind=4) Matrix(:, :)
      integer M, N, lda
      character, optional:: norm
      character::opt
      complex(kind=4), allocatable:: WORK(:)
      real(kind=4) clange

      opt = 'F'
      if (present(norm)) opt = norm
      if (opt == 'I' .or. opt == 'i') then
         allocate (WORK(M))
         WORK = 0
      else
         allocate(WORK(1))
      end if

      lda = size(Matrix, 1)
      clangef90 = clange(opt, M, N, Matrix, lda, WORK)

      ! if (opt == 'I' .or. opt == 'i') then
         deallocate (WORK)
      ! end if

   end function clangef90


   real(kind=8) function pfnorm(M, N, Matrix, ia, ja, desca, norm)
      integer M, N, ia, ja
      DT Matrix(:, :)
      integer desca(9)
      character, optional:: norm

#if DAT==0
         pfnorm = pzlangef90(M, N, Matrix, ia, ja, desca, norm)
#elif DAT==1
         pfnorm = pdlangef90(M, N, Matrix, ia, ja, desca, norm)
#elif DAT==2
         pfnorm = pclangef90(M, N, Matrix, ia, ja, desca, norm)
#elif DAT==3
         pfnorm = pslangef90(M, N, Matrix, ia, ja, desca, norm)
#endif
   end function pfnorm


   real(kind=8) function pdlangef90(M, N, Matrix, ia, ja, desca, norm)
      implicit none
      integer M, N, ia, ja
      real(kind=8) Matrix(:, :)
      integer desca(9)
      character, optional:: norm
      real(kind=8), allocatable:: WORK(:)
      real(kind=8) pdlange
      character::opt
      integer IROFFA,ICOFFA

      IROFFA = MOD( ia-1, nbslpk )
      ICOFFA = MOD( ja-1, nbslpk )
      opt = 'F'
      if (present(norm)) opt = norm
      if (opt == 'I' .or. opt == 'i' .or. opt == '1') then
         allocate (WORK(max(M+IROFFA,N+ICOFFA)))
         WORK = 0
      end if
      pdlangef90 = pdlange(opt, M, N, Matrix, ia, ja, desca, WORK)
      if (opt == 'I' .or. opt == 'i' .or. opt == '1') then
         deallocate (WORK)
      end if
   end function pdlangef90


   real(kind=4) function pslangef90(M, N, Matrix, ia, ja, desca, norm)
      implicit none
      integer M, N, ia, ja
      real(kind=4) Matrix(:, :)
      integer desca(9)
      character, optional:: norm
      real(kind=4), allocatable:: WORK(:)
      real(kind=4) pslange
      character::opt
      integer IROFFA,ICOFFA

      IROFFA = MOD( ia-1, nbslpk )
      ICOFFA = MOD( ja-1, nbslpk )
      opt = 'F'
      if (present(norm)) opt = norm
      if (opt == 'I' .or. opt == 'i' .or. opt == '1') then
         allocate (WORK(max(M+IROFFA,N+ICOFFA)))
         WORK = 0
      end if
      pslangef90 = pslange(opt, M, N, Matrix, ia, ja, desca, WORK)
      if (opt == 'I' .or. opt == 'i' .or. opt == '1') then
         deallocate (WORK)
      end if
   end function pslangef90


   real(kind=8) function pzlangef90(M, N, Matrix, ia, ja, desca, norm)
      implicit none
      integer M, N, ia, ja
      complex(kind=8) Matrix(:, :)
      integer desca(9)
      character, optional:: norm
      complex(kind=8), allocatable:: WORK(:)
      real(kind=8) pzlange
      character::opt
      integer IROFFA,ICOFFA

      IROFFA = MOD( ia-1, nbslpk )
      ICOFFA = MOD( ja-1, nbslpk )
      opt = 'F'
      if (present(norm)) opt = norm
      if (opt == 'I' .or. opt == 'i' .or. opt == '1') then
         allocate (WORK(max(M+IROFFA,N+ICOFFA)))
         WORK = 0
      end if
      pzlangef90 = pzlange(opt, M, N, Matrix, ia, ja, desca, WORK)
      if (opt == 'I' .or. opt == 'i' .or. opt == '1') then
         deallocate (WORK)
      end if
   end function pzlangef90

   real(kind=4) function pclangef90(M, N, Matrix, ia, ja, desca, norm)
      implicit none
      integer M, N, ia, ja
      complex(kind=4) Matrix(:, :)
      integer desca(9)
      character, optional:: norm
      complex(kind=4), allocatable:: WORK(:)
      real(kind=4) pclange
      character::opt
      integer IROFFA,ICOFFA

      IROFFA = MOD( ia-1, nbslpk )
      ICOFFA = MOD( ja-1, nbslpk )
      opt = 'F'
      if (present(norm)) opt = norm
      if (opt == 'I' .or. opt == 'i' .or. opt == '1') then
         allocate (WORK(max(M+IROFFA,N+ICOFFA)))
         WORK = 0
      end if
      pclangef90 = pclange(opt, M, N, Matrix, ia, ja, desca, WORK)
      if (opt == 'I' .or. opt == 'i' .or. opt == '1') then
         deallocate (WORK)
      end if
   end function pclangef90


   subroutine gesvd_robust(Matrix, Singular, UU, VV, mm, nn, mn_min, flop)
      implicit none
      integer mm, nn, mn_min
      DT Matrix(:, :), UU(:, :), VV(:, :),Matrix0(mm,nn)
      DTR Singular(:)
      real(kind=8), optional::flop
      integer INFO

      if (mm == 1) then
         UU(1, 1) = 1d0
         Singular(1) = fnorm(Matrix, mm, nn)
         if (Singular(1) < BPACK_SafeUnderflow) then
            VV = 0d0
         else
            VV = Matrix/Singular(1)
         endif
         if (present(flop)) flop = mm*nn
      elseif (nn == 1) then
         VV(1, 1) = 1d0
         Singular(1) = fnorm(Matrix, mm, nn)
         if (Singular(1) < BPACK_SafeUnderflow) then
            UU = 0d0
         else
            UU = Matrix/Singular(1)
         endif
         if (present(flop)) flop = mm*nn
      else
         if (fnorm(Matrix, mm, nn) < BPACK_SafeUnderflow) then
            Singular = 0
            UU = 0
            VV = 0
            if (present(flop))flop=0
         else
            ! write(*,*)'ga',fnorm(Matrix,mm,nn),shape(Matrix),shape(UU),shape(VV)
            Singular = 0
            UU = 0
            VV = 0
            Matrix0 =Matrix(1:mm,1:nn)
            call gesddf90(Matrix0, Singular, UU, VV, INFO, flop=flop)

            if (INFO /= 0) then
               !!!!!! gesvd (QR iteration) can occasionally fail compared to gesdd (DC)
               call gesvdf90(Matrix,Singular,UU,VV,flop=flop)
            endif
         endif
      endif

   end subroutine gesvd_robust

   subroutine gesvdf90(Matrix, Singular, UU, VV, flop)
      implicit none
      DT::Matrix(:, :), UU(:, :), VV(:, :)
      DTR Singular(:)
      real(kind=8), optional::flop

#if DAT==0
      call zgesvdf90(Matrix, Singular, UU, VV, flop)
#elif DAT==1
      call dgesvdf90(Matrix, Singular, UU, VV, flop)
#elif DAT==2
      call cgesvdf90(Matrix, Singular, UU, VV, flop)
#elif DAT==3
      call sgesvdf90(Matrix, Singular, UU, VV, flop)
#endif


   end subroutine gesvdf90

   subroutine zgesvdf90(Matrix, Singular, UU, VV, flop)

!

      implicit none
      complex(kind=8) Matrix(:, :), UU(:, :), VV(:, :)
      real(kind=8) Singular(:)
      integer m, n, mn_min
      real(kind=8), optional::flop

      integer LWORK, INFO
      complex(kind=8):: TEMP(1)
      real(kind=8), allocatable::RWORK(:)

      complex(kind=8), allocatable:: WORK(:)

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      allocate (RWORK(5*mn_min))
      RWORK = 0
      LWORK = -1
      call ZGESVD('S', 'S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, TEMP, LWORK, RWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))  ! increase this 2.001 factor when not converging
      allocate (WORK(LWORK))
      WORK = 0

! write(*,*)'sssss', TEMP(1),LWORK,fnorm(Matrix,m,n)

      call ZGESVD('S', 'S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, WORK, LWORK, RWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'SVD failed!!', INFO
         stop
      endif

      deallocate (WORK, RWORK)

      if (present(flop)) flop = flops_zgesvd(m, n)

   end subroutine zgesvdf90


   subroutine cgesvdf90(Matrix, Singular, UU, VV, flop)
   !
         implicit none
         complex(kind=4) Matrix(:, :), UU(:, :), VV(:, :)
         real(kind=4) Singular(:)
         integer m, n, mn_min
         real(kind=8), optional::flop

         integer LWORK, INFO
         complex(kind=4):: TEMP(1)
         real(kind=4), allocatable::RWORK(:)

         complex(kind=4), allocatable:: WORK(:)

         m = size(Matrix, 1)
         n = size(Matrix, 2)
         mn_min = min(m, n)

         allocate (RWORK(5*mn_min))
         RWORK = 0
         LWORK = -1
         call CGESVD('S', 'S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, TEMP, LWORK, RWORK, INFO)

         LWORK = NINT(dble(TEMP(1)*2.001 + 1))  ! increase this 2.001 factor when not converging
         allocate (WORK(LWORK))
         WORK = 0

   ! write(*,*)'sssss', TEMP(1),LWORK,fnorm(Matrix,m,n)

         call CGESVD('S', 'S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, WORK, LWORK, RWORK, INFO)

         if (INFO /= 0) then
            write (*, *) 'SVD failed!!', INFO
            stop
         endif

         deallocate (WORK, RWORK)

         if (present(flop)) flop = flops_zgesvd(m, n)
   end subroutine cgesvdf90


   subroutine dgesvdf90(Matrix, Singular, UU, VV, flop)

!

      implicit none
      real(kind=8) Matrix(:, :), UU(:, :), VV(:, :)
      real(kind=8) Singular(:)
      integer m, n, mn_min

      integer LWORK, INFO
      real(kind=8):: TEMP(1)
      real(kind=8), optional::flop

      real(kind=8), allocatable:: WORK(:)

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      LWORK = -1
      call DGESVD('S', 'S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))  ! increase this 2.001 factor when not converging
      allocate (WORK(LWORK))
      WORK = 0

! write(*,*)'sssss', TEMP(1),LWORK,fnorm(Matrix,m,n)

      call DGESVD('S', 'S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'SVD failed!!', INFO
         stop
      endif

      deallocate (WORK)
      if (present(flop)) flop = flops_dgesvd(m, n)
   end subroutine dgesvdf90


   subroutine sgesvdf90(Matrix, Singular, UU, VV, flop)
      !
      implicit none
      real(kind=4) Matrix(:, :), UU(:, :), VV(:, :)
      real(kind=4) Singular(:)
      integer m, n, mn_min

      integer LWORK, INFO
      real(kind=4):: TEMP(1)
      real(kind=8), optional::flop

      real(kind=4), allocatable:: WORK(:)

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      LWORK = -1
      call SGESVD('S', 'S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))  ! increase this 2.001 factor when not converging
      allocate (WORK(LWORK))
      WORK = 0

! write(*,*)'sssss', TEMP(1),LWORK,fnorm(Matrix,m,n)

      call SGESVD('S', 'S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'SVD failed!!', INFO
         stop
      endif

      deallocate (WORK)
      if (present(flop)) flop = flops_dgesvd(m, n)
   end subroutine sgesvdf90


   subroutine gesddf90(Matrix, Singular, UU, VV, INFO, flop)
      implicit none
      DT::Matrix(:, :), UU(:, :), VV(:, :)
      DTR Singular(:)
      real(kind=8), optional::flop
      integer INFO

#if DAT==0
      call zgesddf90(Matrix, Singular, UU, VV, INFO, flop)
#elif DAT==1
      call dgesddf90(Matrix, Singular, UU, VV, INFO, flop)
#elif DAT==2
      call cgesddf90(Matrix, Singular, UU, VV, INFO, flop)
#elif DAT==3
      call sgesddf90(Matrix, Singular, UU, VV, INFO, flop)
#endif




   end subroutine gesddf90

   subroutine zgesddf90(Matrix, Singular, UU, VV, INFO, flop)
!
      implicit none
      complex(kind=8) Matrix(:, :), UU(:, :), VV(:, :)
      real(kind=8) Singular(:)
      integer m, n, mn_min

      integer LWORK, INFO
      complex(kind=8):: TEMP(1)
      real(kind=8), allocatable::RWORK(:)
      integer, allocatable::IWORK(:)

      complex(kind=8), allocatable:: WORK(:)
      real(kind=8), optional::flop

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      allocate (RWORK(max(1, min(m, n)*max(5*min(m, n) + 7, 2*max(m, n) + 2*min(m, n) + 1))))
      allocate (IWORK(8*mn_min))
      RWORK = 0
      LWORK = -1
      call ZGESDD('S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, TEMP, LWORK, RWORK, IWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0

! write(*,*)'sssss', TEMP(1),LWORK,fnorm(Matrix,m,n)

      call ZGESDD('S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, WORK, LWORK, RWORK, IWORK, INFO)

      ! if (INFO /= 0) then
      !    write (*, *) 'SDD failed!!', INFO
      !    stop
      ! endif

      deallocate (WORK, RWORK, IWORK)
      if (present(flop)) flop = flops_zgesdd(m, n)
   end subroutine zgesddf90

   subroutine cgesddf90(Matrix, Singular, UU, VV, INFO, flop)
      !
      implicit none
      complex(kind=4) Matrix(:, :), UU(:, :), VV(:, :)
      real(kind=4) Singular(:)
      integer m, n, mn_min

      integer LWORK, INFO
      complex(kind=4):: TEMP(1)
      real(kind=4), allocatable::RWORK(:)
      integer, allocatable::IWORK(:)

      complex(kind=4), allocatable:: WORK(:)
      real(kind=8), optional::flop

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      allocate (RWORK(max(1, min(m, n)*max(5*min(m, n) + 7, 2*max(m, n) + 2*min(m, n) + 1))))
      allocate (IWORK(8*mn_min))
      RWORK = 0
      LWORK = -1
      call CGESDD('S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, TEMP, LWORK, RWORK, IWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0

! write(*,*)'sssss', TEMP(1),LWORK,fnorm(Matrix,m,n)

      call CGESDD('S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, WORK, LWORK, RWORK, IWORK, INFO)

      ! if (INFO /= 0) then
      !    write (*, *) 'SDD failed!!', INFO
      !    stop
      ! endif

      deallocate (WORK, RWORK, IWORK)
      if (present(flop)) flop = flops_zgesdd(m, n)
   end subroutine cgesddf90


   subroutine dgesddf90(Matrix, Singular, UU, VV, INFO, flop)

!

      implicit none
      real(kind=8) Matrix(:, :), UU(:, :), VV(:, :)
      real(kind=8) Singular(:)
      integer m, n, mn_min

      integer LWORK, INFO
      real(kind=8):: TEMP(1)
      real(kind=8), allocatable::RWORK(:)
      integer, allocatable::IWORK(:)

      real(kind=8), allocatable:: WORK(:)
      real(kind=8), optional::flop

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      allocate (IWORK(8*mn_min))
      LWORK = -1
      call DGESDD('S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, TEMP, LWORK, IWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0

! write(*,*)'sssss', TEMP(1),LWORK,fnorm(Matrix,m,n)

      call DGESDD('S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, WORK, LWORK, IWORK, INFO)

      ! if (INFO /= 0) then
      !    write (*, *) 'SDD failed!!', INFO
      !    stop
      ! endif

      deallocate (WORK, IWORK)
      if (present(flop)) flop = flops_dgesdd(m, n)
   end subroutine dgesddf90

   subroutine sgesddf90(Matrix, Singular, UU, VV, INFO, flop)
      !
      implicit none
      real(kind=4) Matrix(:, :), UU(:, :), VV(:, :)
      real(kind=4) Singular(:)
      integer m, n, mn_min

      integer LWORK, INFO
      real(kind=4):: TEMP(1)
      real(kind=4), allocatable::RWORK(:)
      integer, allocatable::IWORK(:)

      real(kind=4), allocatable:: WORK(:)
      real(kind=8), optional::flop

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      allocate (IWORK(8*mn_min))
      LWORK = -1
      call SGESDD('S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, TEMP, LWORK, IWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0

! write(*,*)'sssss', TEMP(1),LWORK,fnorm(Matrix,m,n)

      call SGESDD('S', m, n, Matrix, m, Singular, UU, m, VV, mn_min, WORK, LWORK, IWORK, INFO)

      ! if (INFO /= 0) then
      !    write (*, *) 'SDD failed!!', INFO
      !    stop
      ! endif

      deallocate (WORK, IWORK)
      if (present(flop)) flop = flops_dgesdd(m, n)
   end subroutine sgesddf90

   subroutine geqrff90(Matrix, tau, flop)
      implicit none
      DT::Matrix(:, :), tau(:)
      real(kind=8), optional::flop

#if DAT==0
            call zgeqrff90(Matrix, tau, flop)
#elif DAT==1
            call dgeqrff90(Matrix, tau, flop)
#elif DAT==2
            call cgeqrff90(Matrix, tau, flop)
#elif DAT==3
            call sgeqrff90(Matrix, tau, flop)
#endif


   end subroutine geqrff90

   subroutine zgeqrff90(Matrix, tau, flop)

!

      implicit none
      complex(kind=8) Matrix(:, :), tau(:)
      integer m, n, mn_min
      real(kind=8), optional::flop

      integer LWORK, INFO
      complex(kind=8):: TEMP(1)

      complex(kind=8), allocatable:: WORK(:)

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      LWORK = -1
      call ZGEQRF(m, n, Matrix, m, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call ZGEQRF(m, n, Matrix, m, tau, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'zgeqrff90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call geqrf(Matrix,tau)
      if (present(flop)) flop = flops_zgetrf(m, n)
   end subroutine zgeqrff90

   subroutine cgeqrff90(Matrix, tau, flop)
      !
      implicit none
      complex(kind=4) Matrix(:, :), tau(:)
      integer m, n, mn_min
      real(kind=8), optional::flop

      integer LWORK, INFO
      complex(kind=4):: TEMP(1)

      complex(kind=4), allocatable:: WORK(:)

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      LWORK = -1
      call CGEQRF(m, n, Matrix, m, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call CGEQRF(m, n, Matrix, m, tau, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'cgeqrff90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call geqrf(Matrix,tau)
      if (present(flop)) flop = flops_zgetrf(m, n)
   end subroutine cgeqrff90

   subroutine dgeqrff90(Matrix, tau, flop)

!

      implicit none
      real(kind=8) Matrix(:, :), tau(:)
      integer m, n, mn_min

      integer LWORK, INFO
      real(kind=8):: TEMP(1)

      real(kind=8), allocatable:: WORK(:)
      real(kind=8), optional::flop

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      LWORK = -1
      call DGEQRF(m, n, Matrix, m, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call DGEQRF(m, n, Matrix, m, tau, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'dgeqrff90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call geqrf(Matrix,tau)
      if (present(flop)) flop = flops_dgetrf(m, n)
   end subroutine dgeqrff90

   subroutine sgeqrff90(Matrix, tau, flop)
      !
      implicit none
      real(kind=4) Matrix(:, :), tau(:)
      integer m, n, mn_min

      integer LWORK, INFO
      real(kind=4):: TEMP(1)

      real(kind=4), allocatable:: WORK(:)
      real(kind=8), optional::flop

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      LWORK = -1
      call SGEQRF(m, n, Matrix, m, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call SGEQRF(m, n, Matrix, m, tau, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'sgeqrff90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call geqrf(Matrix,tau)
      if (present(flop)) flop = flops_dgetrf(m, n)
   end subroutine sgeqrff90


   subroutine geqp3f90(Matrix, jpvt, tau, flop)

!

      implicit none
      DT:: Matrix(:, :), tau(:)
      integer jpvt(:)
      real(kind=8), optional::flop


#if DAT==0
            call zgeqp3f90(Matrix, jpvt, tau, flop)
#elif DAT==1
            call dgeqp3f90(Matrix, jpvt, tau, flop)
#elif DAT==2
            call cgeqp3f90(Matrix, jpvt, tau, flop)
#elif DAT==3
            call sgeqp3f90(Matrix, jpvt, tau, flop)
#endif


   end subroutine geqp3f90

   subroutine zgeqp3f90(Matrix, jpvt, tau, flop)

!

      implicit none
      complex(kind=8) Matrix(:, :), tau(:)
      integer jpvt(:)
      integer m, n, mn_min

      integer LWORK, INFO
      complex(kind=8):: TEMP(1)
      real(kind=8), allocatable::RWORK(:)

      complex(kind=8), allocatable:: WORK(:)

      real(kind=8), optional::flop

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      allocate (RWORK(2*max(m, n)))
      RWORK = 0
      LWORK = -1
      call ZGEQP3(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, RWORK, INFO)
      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call ZGEQP3(m, n, Matrix, m, jpvt, tau, WORK, LWORK, RWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'geqp3f90 failed!!', INFO
         stop
      endif

      deallocate (WORK)
      deallocate (RWORK)

! call geqp3(Matrix,jpvt,tau)
      if (present(flop)) flop = flops_zgetrf(m, n)
   end subroutine zgeqp3f90

   subroutine cgeqp3f90(Matrix, jpvt, tau, flop)
      !
      implicit none
      complex(kind=4) Matrix(:, :), tau(:)
      integer jpvt(:)
      integer m, n, mn_min

      integer LWORK, INFO
      complex(kind=4):: TEMP(1)
      real(kind=4), allocatable::RWORK(:)

      complex(kind=4), allocatable:: WORK(:)

      real(kind=8), optional::flop

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      allocate (RWORK(2*max(m, n)))
      RWORK = 0
      LWORK = -1
      call CGEQP3(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, RWORK, INFO)
      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call CGEQP3(m, n, Matrix, m, jpvt, tau, WORK, LWORK, RWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'geqp3f90 failed!!', INFO
         stop
      endif

      deallocate (WORK)
      deallocate (RWORK)

! call geqp3(Matrix,jpvt,tau)
      if (present(flop)) flop = flops_zgetrf(m, n)
   end subroutine cgeqp3f90

   subroutine dgeqp3f90(Matrix, jpvt, tau, flop)

!

      implicit none
      real(kind=8) Matrix(:, :), tau(:)
      integer jpvt(:)
      integer m, n, mn_min

      integer LWORK, INFO
      real(kind=8):: TEMP(1)

      real(kind=8), allocatable:: WORK(:)
      real(kind=8), optional::flop

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      LWORK = -1
      call DGEQP3(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, INFO)
      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call DGEQP3(m, n, Matrix, m, jpvt, tau, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'geqp3f90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call geqp3(Matrix,jpvt,tau)
      if (present(flop)) flop = flops_dgetrf(m, n)
   end subroutine dgeqp3f90

   subroutine sgeqp3f90(Matrix, jpvt, tau, flop)
      !
      implicit none
      real(kind=4) Matrix(:, :), tau(:)
      integer jpvt(:)
      integer m, n, mn_min

      integer LWORK, INFO
      real(kind=4):: TEMP(1)

      real(kind=4), allocatable:: WORK(:)
      real(kind=8), optional::flop

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      LWORK = -1
      call SGEQP3(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, INFO)
      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call SGEQP3(m, n, Matrix, m, jpvt, tau, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'geqp3f90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call geqp3(Matrix,jpvt,tau)
      if (present(flop)) flop = flops_dgetrf(m, n)
   end subroutine sgeqp3f90


   subroutine geqp3modf90(Matrix, jpvt, tau, rtol, atol, rank, flop)

!

      implicit none
      DT:: Matrix(:, :), tau(:)
      integer jpvt(:)
      real(kind=8)::rtol, atol
      integer rank
      real(kind=8), optional::flop
      integer m,n,mn_min

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)
      if(mn_min==0)then
         rank=0
      else
#if DAT==0
         call zgeqp3modf90(Matrix, jpvt, tau, rtol, atol, rank, flop)
#elif DAT==1
         call dgeqp3modf90(Matrix, jpvt, tau, rtol, atol, rank, flop)
#elif DAT==2
         call cgeqp3modf90(Matrix, jpvt, tau, real(rtol), real(atol), rank, flop)
#elif DAT==3
         call sgeqp3modf90(Matrix, jpvt, tau, real(rtol), real(atol), rank, flop)
#endif

      endif

   end subroutine geqp3modf90

   subroutine dgeqp3modf90(Matrix, jpvt, tau, rtol, atol, rank, flop)

!

      implicit none
      real(kind=8) Matrix(:, :), tau(:)
      integer jpvt(:)
      integer m, n, mn_min

      real(kind=8)::rtol, atol
      integer LWORK, INFO, rank
      real(kind=8):: TEMP(1)
      real(kind=8), optional::flop
      real(kind=8), allocatable:: WORK(:)

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)
      rank = 0

      LWORK = -1
! call DGEQP3(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, INFO)
      call DGEQP3mod(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, INFO, RANK, RTOL, ATOL)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
! call DGEQP3(m, n, Matrix, m, jpvt, tau, WORK, LWORK, INFO)
      call DGEQP3mod(m, n, Matrix, m, jpvt, tau, WORK, LWORK, INFO, RANK, RTOL, ATOL)

      if (INFO /= 0) then
         write (*, *) 'geqp3modf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call geqp3(Matrix,jpvt,tau)
      if (present(flop)) flop = flops_dgeqpfmod(m, n, rank)
   end subroutine dgeqp3modf90


   subroutine sgeqp3modf90(Matrix, jpvt, tau, rtol, atol, rank, flop)
      !
      implicit none
      real(kind=4) Matrix(:, :), tau(:)
      integer jpvt(:)
      integer m, n, mn_min

      real(kind=4)::rtol, atol
      integer LWORK, INFO, rank
      real(kind=4):: TEMP(1)
      real(kind=8), optional::flop
      real(kind=4), allocatable:: WORK(:)

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)
      rank = 0

      LWORK = -1
! call DGEQP3(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, INFO)
      call SGEQP3mod(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, INFO, RANK, RTOL, ATOL)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
! call DGEQP3(m, n, Matrix, m, jpvt, tau, WORK, LWORK, INFO)
      call SGEQP3mod(m, n, Matrix, m, jpvt, tau, WORK, LWORK, INFO, RANK, RTOL, ATOL)

      if (INFO /= 0) then
         write (*, *) 'geqp3modf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call geqp3(Matrix,jpvt,tau)
      if (present(flop)) flop = flops_dgeqpfmod(m, n, rank)
   end subroutine sgeqp3modf90


   subroutine zgeqp3modf90(Matrix, jpvt, tau, rtol, atol, rank, flop)
!
      implicit none
      complex(kind=8) Matrix(:, :), tau(:)
      integer jpvt(:)
      integer m, n, mn_min

      real(kind=8)::rtol, atol
      integer LWORK, INFO, rank
      complex(kind=8):: TEMP(1)
      real(kind=8), allocatable::RWORK(:)

      complex(kind=8), allocatable:: WORK(:)
      real(kind=8), optional::flop

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)
      rank = 0

      allocate (RWORK(2*max(m, n)))
      RWORK = 0
      LWORK = -1
! call ZGEQP3(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, RWORK, INFO)
      call ZGEQP3mod(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, RWORK, INFO, RANK, RTOL, ATOL)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
! call ZGEQP3(m, n, Matrix, m, jpvt, tau, WORK, LWORK, RWORK, INFO)
      call ZGEQP3mod(m, n, Matrix, m, jpvt, tau, WORK, LWORK, RWORK, INFO, RANK, RTOL, ATOL)

      if (INFO /= 0) then
         write (*, *) 'geqp3modf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)
      deallocate (RWORK)

! call geqp3(Matrix,jpvt,tau)
      if (present(flop)) flop = flops_zgeqpfmod(m, n, rank)
   end subroutine zgeqp3modf90

   subroutine cgeqp3modf90(Matrix, jpvt, tau, rtol, atol, rank, flop)
      !
      implicit none
      complex(kind=4) Matrix(:, :), tau(:)
      integer jpvt(:)
      integer m, n, mn_min

      real(kind=4)::rtol, atol
      integer LWORK, INFO, rank
      complex(kind=4):: TEMP(1)
      real(kind=4), allocatable::RWORK(:)

      complex(kind=4), allocatable:: WORK(:)
      real(kind=8), optional::flop

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)
      rank = 0

      allocate (RWORK(2*max(m, n)))
      RWORK = 0
      LWORK = -1
! call ZGEQP3(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, RWORK, INFO)
      call CGEQP3mod(m, n, Matrix, m, jpvt, tau, TEMP, LWORK, RWORK, INFO, RANK, RTOL, ATOL)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
! call ZGEQP3(m, n, Matrix, m, jpvt, tau, WORK, LWORK, RWORK, INFO)
      call CGEQP3mod(m, n, Matrix, m, jpvt, tau, WORK, LWORK, RWORK, INFO, RANK, RTOL, ATOL)

      if (INFO /= 0) then
         write (*, *) 'geqp3modf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)
      deallocate (RWORK)

! call geqp3(Matrix,jpvt,tau)
      if (present(flop)) flop = flops_zgeqpfmod(m, n, rank)
   end subroutine cgeqp3modf90

   subroutine un_or_mqrf90(a, tau, c, side, trans, m, n, k, flop)

!

      implicit none
      DT a(:, :), tau(:), c(:, :)
      integer m, n, k
      character:: side, trans
      real(kind=8), optional::flop

#if DAT==0
         call zunmqrf90(a, tau, c, side, trans, m, n, k, flop)
#elif DAT==1
         call dormqrf90(a, tau, c, side, trans, m, n, k, flop)
#elif DAT==2
         call cunmqrf90(a, tau, c, side, trans, m, n, k, flop)
#elif DAT==3
         call sormqrf90(a, tau, c, side, trans, m, n, k, flop)
#endif


   end subroutine un_or_mqrf90

   subroutine dormqrf90(a, tau, c, side, trans, m, n, k, flop)
!
      implicit none
      real(kind=8) a(:, :), tau(:), c(:, :)
      integer m, n, k, lda, ldc, mn_min

      character:: side, trans, trans1

      integer LWORK, INFO
      real(kind=8):: TEMP(1)
      real(kind=8), allocatable:: WORK(:)
      real(kind=8), optional::flop

      trans1 = trans
      if (trans1 == 'C') trans1 = 'T'

      ldc = size(c, 1)
      lda = size(a, 1)

      LWORK = -1
      call DORMQR(side, trans1, m, n, k, a, lda, tau, c, ldc, TEMP, LWORK, INFO)
      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call DORMQR(side, trans1, m, n, k, a, lda, tau, c, ldc, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'dormqrf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call unmqr(a,tau,c,side,trans)
      if (present(flop)) flop = flops_dunmqr(side, m, n, k)
   end subroutine dormqrf90

   subroutine sormqrf90(a, tau, c, side, trans, m, n, k, flop)
      !
      implicit none
      real(kind=4) a(:, :), tau(:), c(:, :)
      integer m, n, k, lda, ldc, mn_min

      character:: side, trans, trans1

      integer LWORK, INFO
      real(kind=4):: TEMP(1)
      real(kind=4), allocatable:: WORK(:)
      real(kind=8), optional::flop

      trans1 = trans
      if (trans1 == 'C') trans1 = 'T'

      ldc = size(c, 1)
      lda = size(a, 1)

      LWORK = -1
      call SORMQR(side, trans1, m, n, k, a, lda, tau, c, ldc, TEMP, LWORK, INFO)
      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call SORMQR(side, trans1, m, n, k, a, lda, tau, c, ldc, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'sormqrf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call unmqr(a,tau,c,side,trans)
      if (present(flop)) flop = flops_dunmqr(side, m, n, k)
   end subroutine sormqrf90


   subroutine zunmqrf90(a, tau, c, side, trans, m, n, k, flop)
!
      implicit none
      complex(kind=8) a(:, :), tau(:), c(:, :)
      integer m, n, k, lda, ldc, mn_min
      real(kind=8), optional::flop
      character:: side, trans

      integer LWORK, INFO
      complex(kind=8):: TEMP(1)
      complex(kind=8), allocatable:: WORK(:)

      ldc = size(c, 1)
      lda = size(a, 1)

      LWORK = -1
      call ZUNMQR(side, trans, m, n, k, a, lda, tau, c, ldc, TEMP, LWORK, INFO)
      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call ZUNMQR(side, trans, m, n, k, a, lda, tau, c, ldc, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'zunmqrf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call unmqr(a,tau,c,side,trans)
      if (present(flop)) flop = flops_zunmqr(side, m, n, k)
   end subroutine zunmqrf90

   subroutine cunmqrf90(a, tau, c, side, trans, m, n, k, flop)
      !
      implicit none
      complex(kind=4) a(:, :), tau(:), c(:, :)
      integer m, n, k, lda, ldc, mn_min
      real(kind=8), optional::flop
      character:: side, trans

      integer LWORK, INFO
      complex(kind=4):: TEMP(1)
      complex(kind=4), allocatable:: WORK(:)

      ldc = size(c, 1)
      lda = size(a, 1)

      LWORK = -1
      call CUNMQR(side, trans, m, n, k, a, lda, tau, c, ldc, TEMP, LWORK, INFO)
      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call CUNMQR(side, trans, m, n, k, a, lda, tau, c, ldc, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'cunmqrf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call unmqr(a,tau,c,side,trans)
      if (present(flop)) flop = flops_zunmqr(side, m, n, k)
   end subroutine cunmqrf90


   subroutine un_or_gqrf90(Matrix, tau, m, n, k, flop)
      DT Matrix(:, :), tau(:)
      real(kind=8), optional::flop
      integer m, n, k


#if DAT==0
            call zungqrf90(Matrix, tau, m, n, k, flop)
#elif DAT==1
            call dorgqrf90(Matrix, tau, m, n, k, flop)
#elif DAT==2
            call cungqrf90(Matrix, tau, m, n, k, flop)
#elif DAT==3
            call sorgqrf90(Matrix, tau, m, n, k, flop)
#endif



   end subroutine un_or_gqrf90

   subroutine pun_or_gqrf90(ctxt, Matrix, tau, m, n, k, desca, ia, ja, flop)
      DT Matrix(:, :), tau(:)
      real(kind=8), optional::flop
      integer m, n, k, ia, ja
      integer desca(9)
      character ROWBTOP, COLBTOP
      integer ctxt

      CALL PB_TOPGET(ctxt, 'Broadcast', 'Rowwise', ROWBTOP)
      CALL PB_TOPGET(ctxt, 'Broadcast', 'Columnwise', COLBTOP)

#if DAT==0
      call zpungqrf90(Matrix, tau, m, n, k, desca, ia, ja, flop)
#elif DAT==1
      call dporgqrf90(Matrix, tau, m, n, k, desca, ia, ja, flop)
#elif DAT==2
      call cpungqrf90(Matrix, tau, m, n, k, desca, ia, ja, flop)
#elif DAT==3
      call sporgqrf90(Matrix, tau, m, n, k, desca, ia, ja, flop)
#endif



      CALL PB_TOPSET(ctxt, 'Broadcast', 'Rowwise', ROWBTOP)
      CALL PB_TOPSET(ctxt, 'Broadcast', 'Columnwise', COLBTOP)

   end subroutine pun_or_gqrf90

   subroutine dorgqrf90(Matrix, tau, m, n, k, flop)
!
      implicit none
      real(kind=8) Matrix(:, :), tau(:)
      integer m, n, k, lda

      integer LWORK, INFO
      real(kind=8):: TEMP(1)

      real(kind=8), allocatable:: WORK(:)
      real(kind=8), optional::flop

      lda = size(Matrix, 1)

      LWORK = -1
      call DORGQR(m, n, k, Matrix, lda, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call DORGQR(m, n, k, Matrix, lda, tau, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'dorgqrf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call ungqr(Matrix,tau)
      if (present(flop)) flop = flops_dungqr(m, n, k)
   end subroutine dorgqrf90

   subroutine sorgqrf90(Matrix, tau, m, n, k, flop)
      !
      implicit none
      real(kind=4) Matrix(:, :), tau(:)
      integer m, n, k, lda

      integer LWORK, INFO
      real(kind=4):: TEMP(1)

      real(kind=4), allocatable:: WORK(:)
      real(kind=8), optional::flop

      lda = size(Matrix, 1)

      LWORK = -1
      call SORGQR(m, n, k, Matrix, lda, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call SORGQR(m, n, k, Matrix, lda, tau, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'sorgqrf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call ungqr(Matrix,tau)
      if (present(flop)) flop = flops_dungqr(m, n, k)
   end subroutine sorgqrf90

   subroutine dporgqrf90(Matrix, tau, m, n, k, desca, ia, ja, flop)
!
      implicit none
      real(kind=8) Matrix(:, :), tau(:)
      integer m, n, k, lda

      integer LWORK, INFO
      real(kind=8):: TEMP(1)
      integer desca(9)
      integer ia, ja
      real(kind=8), allocatable:: WORK(:)
      real(kind=8), optional::flop

      LWORK = -1
      call PDORGQR(m, n, k, Matrix, ia, ja, desca, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call PDORGQR(m, n, k, Matrix, ia, ja, desca, tau, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'dporgqrf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call ungqr(Matrix,tau)
      if (present(flop)) flop = flops_dungqr(m, n, k)
   end subroutine dporgqrf90

   subroutine sporgqrf90(Matrix, tau, m, n, k, desca, ia, ja, flop)
      !
      implicit none
      real(kind=4) Matrix(:, :), tau(:)
      integer m, n, k, lda

      integer LWORK, INFO
      real(kind=4):: TEMP(1)
      integer desca(9)
      integer ia, ja
      real(kind=4), allocatable:: WORK(:)
      real(kind=8), optional::flop

      LWORK = -1
      call PSORGQR(m, n, k, Matrix, ia, ja, desca, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call PSORGQR(m, n, k, Matrix, ia, ja, desca, tau, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'sporgqrf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call ungqr(Matrix,tau)
      if (present(flop)) flop = flops_dungqr(m, n, k)
   end subroutine sporgqrf90


   subroutine zungqrf90(Matrix, tau, m, n, k, flop)
!
      implicit none
      complex(kind=8) Matrix(:, :), tau(:)
      integer m, n, k, lda
      real(kind=8), optional::flop
      integer LWORK, INFO
      complex(kind=8):: TEMP(1)

      complex(kind=8), allocatable:: WORK(:)

      lda = size(Matrix, 1)

      LWORK = -1
      call ZUNGQR(m, n, k, Matrix, lda, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call ZUNGQR(m, n, k, Matrix, lda, tau, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'zungqrf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

      if (present(flop)) flop = flops_zungqr(m, n, k)
! call ungqr(Matrix,tau)

   end subroutine zungqrf90


   subroutine cungqrf90(Matrix, tau, m, n, k, flop)
      !
      implicit none
      complex(kind=4) Matrix(:, :), tau(:)
      integer m, n, k, lda
      real(kind=8), optional::flop
      integer LWORK, INFO
      complex(kind=4):: TEMP(1)

      complex(kind=4), allocatable:: WORK(:)

      lda = size(Matrix, 1)

      LWORK = -1
      call CUNGQR(m, n, k, Matrix, lda, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call CUNGQR(m, n, k, Matrix, lda, tau, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'cungqrf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

      if (present(flop)) flop = flops_zungqr(m, n, k)
! call ungqr(Matrix,tau)

   end subroutine cungqrf90

   subroutine zpungqrf90(Matrix, tau, m, n, k, desca, ia, ja, flop)
!
      implicit none
      complex(kind=8) Matrix(:, :), tau(:)
      integer m, n, k, lda
      real(kind=8), optional::flop
      integer LWORK, INFO
      complex(kind=8):: TEMP(1)
      integer desca(9)
      integer ia, ja
      complex(kind=8), allocatable:: WORK(:)

      LWORK = -1
      call PZUNGQR(m, n, k, Matrix, ia, ja, desca, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call PZUNGQR(m, n, k, Matrix, ia, ja, desca, tau, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'zpungqrf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

      if (present(flop)) flop = flops_zungqr(m, n, k)
! call ungqr(Matrix,tau)

   end subroutine zpungqrf90


   subroutine cpungqrf90(Matrix, tau, m, n, k, desca, ia, ja, flop)
      !
      implicit none
      complex(kind=4) Matrix(:, :), tau(:)
      integer m, n, k, lda
      real(kind=8), optional::flop
      integer LWORK, INFO
      complex(kind=4):: TEMP(1)
      integer desca(9)
      integer ia, ja
      complex(kind=4), allocatable:: WORK(:)

      LWORK = -1
      call PCUNGQR(m, n, k, Matrix, ia, ja, desca, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call PCUNGQR(m, n, k, Matrix, ia, ja, desca, tau, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'cpungqrf90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

      if (present(flop)) flop = flops_zungqr(m, n, k)
! call ungqr(Matrix,tau)

   end subroutine cpungqrf90

   subroutine getrff90(Matrix, ipiv, flop)
      DT Matrix(:, :)
      integer ipiv(:)
      real(kind=8), optional::flop


#if DAT==0
      call zgetrff90(Matrix, ipiv, flop)
#elif DAT==1
      call dgetrff90(Matrix, ipiv, flop)
#elif DAT==2
      call cgetrff90(Matrix, ipiv, flop)
#elif DAT==3
      call sgetrff90(Matrix, ipiv, flop)
#endif


   end subroutine getrff90

#if DAT==0
   subroutine zgetrff90(Matrix, ipiv, flop)
!
      implicit none
      complex(kind=8) Matrix(:, :)
      integer ipiv(:)
      integer m, n, mn_min, ii
      real(kind=8), optional::flop
      integer INFO
      real(kind=8)::norm1,dlamch

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      call ZGETRF(m, n, Matrix, m, ipiv, INFO)

      if (INFO /= 0) then
         write (*, *) 'getrff90 failed!!', INFO
         stop
      endif

! call getrf(Matrix,ipiv)
      if (present(flop)) flop = flops_zgeqrf(m, n)


      norm1 = fnorm(Matrix,m, n,'1')
      do ii=1,mn_min
         if(abs(Matrix(ii, ii))<dlamch('E'))then
            Matrix(ii, ii) = sqrt(dlamch('E'))*norm1
         endif
      enddo

   end subroutine zgetrff90
#elif DAT==1
   subroutine dgetrff90(Matrix, ipiv, flop)
   !
      implicit none
      real(kind=8) Matrix(:, :)
      integer ipiv(:)
      integer m, n, mn_min, ii

      integer INFO
      real(kind=8), optional::flop
      real(kind=8)::norm1,dlamch

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      call DGETRF(m, n, Matrix, m, ipiv, INFO)

      if (INFO /= 0) then
         write (*, *) 'getrff90 failed!!', INFO
         stop
      endif

! call getrf(Matrix,ipiv)
      if (present(flop)) flop = flops_dgeqrf(m, n)

      norm1 = fnorm(Matrix,m, n,'1')
      do ii=1,mn_min
         if(abs(Matrix(ii, ii))<dlamch('E'))then
            Matrix(ii, ii) = sign(1d0,dble(Matrix(ii, ii)))*sqrt(dlamch('E'))*norm1
         endif
      enddo
   end subroutine dgetrff90

#elif DAT==2
   subroutine cgetrff90(Matrix, ipiv, flop)
!
      implicit none
      complex(kind=4) Matrix(:, :)
      integer ipiv(:)
      integer m, n, mn_min, ii
      real(kind=8), optional::flop
      integer INFO
      real(kind=4)::norm1,slamch

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      call CGETRF(m, n, Matrix, m, ipiv, INFO)

      if (INFO /= 0) then
         write (*, *) 'getrff90 failed!!', INFO
         stop
      endif

! call getrf(Matrix,ipiv)
      if (present(flop)) flop = flops_zgeqrf(m, n)


      norm1 = fnorm(Matrix,m, n,'1')
      do ii=1,mn_min
         if(abs(Matrix(ii, ii))<slamch('E'))then
            Matrix(ii, ii) = sqrt(slamch('E'))*norm1
         endif
      enddo
   end subroutine cgetrff90
#elif DAT==3
   subroutine sgetrff90(Matrix, ipiv, flop)
   !
      implicit none
      real(kind=4) Matrix(:, :)
      integer ipiv(:)
      integer m, n, mn_min, ii

      integer INFO
      real(kind=8), optional::flop
      real(kind=4)::norm1,slamch

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      mn_min = min(m, n)

      call SGETRF(m, n, Matrix, m, ipiv, INFO)

      if (INFO /= 0) then
         write (*, *) 'getrff90 failed!!', INFO
         stop
      endif

! call getrf(Matrix,ipiv)
      if (present(flop)) flop = flops_dgeqrf(m, n)

      norm1 = fnorm(Matrix,m, n,'1')
      do ii=1,mn_min
         if(abs(Matrix(ii, ii))<slamch('E'))then
            Matrix(ii, ii) = sign(1d0,dble(Matrix(ii, ii)))*sqrt(slamch('E'))*norm1
         endif
      enddo
   end subroutine sgetrff90
#endif

   subroutine getrsf90(Matrix, ipiv, B, trans, flop)
      DT Matrix(:, :), B(:, :)
      integer ipiv(:)
      character:: trans
      real(kind=8), optional::flop


#if DAT==0
      call zgetrsf90(Matrix, ipiv, B, trans, flop)
#elif DAT==1
      call dgetrsf90(Matrix, ipiv, B, trans, flop)
#elif DAT==2
      call cgetrsf90(Matrix, ipiv, B, trans, flop)
#elif DAT==3
      call sgetrsf90(Matrix, ipiv, B, trans, flop)
#endif



   end subroutine getrsf90

   subroutine dgetrsf90(Matrix, ipiv, B, trans, flop)

!

      implicit none
      real(kind=8) Matrix(:, :), B(:, :)
      integer ipiv(:)
      integer m, n, mn_min, nrhs
      character:: trans, trans1
      real(kind=8), optional::flop
      integer INFO

      n = size(Matrix, 1)
      nrhs = size(B, 2)

      trans1 = trans
      if (trans1 == 'C') trans1 = 'T'

      call DGETRS(trans1, n, nrhs, Matrix, n, ipiv, B, n, INFO)

      if (INFO /= 0) then
         write (*, *) 'getrsf90 failed!!', INFO
         stop
      endif

! call getrs(Matrix,ipiv,B,trans)
      if (present(flop)) flop = flops_dgetrs(n, nrhs)
   end subroutine dgetrsf90

   subroutine sgetrsf90(Matrix, ipiv, B, trans, flop)
      !
      implicit none
      real(kind=4) Matrix(:, :), B(:, :)
      integer ipiv(:)
      integer m, n, mn_min, nrhs
      character:: trans, trans1
      real(kind=8), optional::flop
      integer INFO

      n = size(Matrix, 1)
      nrhs = size(B, 2)

      trans1 = trans
      if (trans1 == 'C') trans1 = 'T'

      call SGETRS(trans1, n, nrhs, Matrix, n, ipiv, B, n, INFO)

      if (INFO /= 0) then
         write (*, *) 'getrsf90 failed!!', INFO
         stop
      endif

! call getrs(Matrix,ipiv,B,trans)
      if (present(flop)) flop = flops_dgetrs(n, nrhs)
   end subroutine sgetrsf90


   subroutine zgetrsf90(Matrix, ipiv, B, trans, flop)
!
      implicit none
      complex(kind=8) Matrix(:, :), B(:, :)
      integer ipiv(:)
      integer m, n, mn_min, nrhs
      character:: trans
      real(kind=8), optional::flop
      integer INFO

      n = size(Matrix, 1)
      nrhs = size(B, 2)

      call ZGETRS(trans, n, nrhs, Matrix, n, ipiv, B, n, INFO)

      if (INFO /= 0) then
         write (*, *) 'getrsf90 failed!!', INFO
         stop
      endif

! call getrs(Matrix,ipiv,B,trans)
      if (present(flop)) flop = flops_zgetrs(n, nrhs)
   end subroutine zgetrsf90

   subroutine cgetrsf90(Matrix, ipiv, B, trans, flop)
      !
      implicit none
      complex(kind=4) Matrix(:, :), B(:, :)
      integer ipiv(:)
      integer m, n, mn_min, nrhs
      character:: trans
      real(kind=8), optional::flop
      integer INFO

      n = size(Matrix, 1)
      nrhs = size(B, 2)

      call CGETRS(trans, n, nrhs, Matrix, n, ipiv, B, n, INFO)

      if (INFO /= 0) then
         write (*, *) 'getrsf90 failed!!', INFO
         stop
      endif

! call getrs(Matrix,ipiv,B,trans)
      if (present(flop)) flop = flops_zgetrs(n, nrhs)
   end subroutine cgetrsf90

   subroutine getrif90(Matrix, ipiv, flop)

      implicit none
      DT Matrix(:, :)
      integer ipiv(:)
      real(kind=8), optional::flop

#if DAT==0
      call zgetrif90(Matrix, ipiv, flop)
#elif DAT==1
      call dgetrif90(Matrix, ipiv, flop)
#elif DAT==2
      call cgetrif90(Matrix, ipiv, flop)
#elif DAT==3
      call sgetrif90(Matrix, ipiv, flop)
#endif


   end subroutine getrif90

   subroutine dgetrif90(Matrix, ipiv, flop)
!
      implicit none
      real(kind=8) Matrix(:, :)
      integer ipiv(:)
      integer m, n, mn_min

      integer LWORK, INFO
      real(kind=8):: TEMP(1)

      real(kind=8), allocatable:: WORK(:)
      real(kind=8), optional::flop
      m = size(Matrix, 1)
      n = size(Matrix, 2)

      if (m > n) then
         write (*, *) 'size of Matrix incorrect'
         stop
      endif

      mn_min = min(m, n)

      LWORK = -1
      call DGETRI(m, Matrix, m, ipiv, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call DGETRI(m, Matrix, m, ipiv, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'getrif90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call getri(Matrix,ipiv)
      if (present(flop)) flop = flops_dgetri(mn_min)
   end subroutine dgetrif90

   subroutine sgetrif90(Matrix, ipiv, flop)
      !
      implicit none
      real(kind=4) Matrix(:, :)
      integer ipiv(:)
      integer m, n, mn_min

      integer LWORK, INFO
      real(kind=4):: TEMP(1)

      real(kind=4), allocatable:: WORK(:)
      real(kind=8), optional::flop
      m = size(Matrix, 1)
      n = size(Matrix, 2)

      if (m > n) then
         write (*, *) 'size of Matrix incorrect'
         stop
      endif

      mn_min = min(m, n)

      LWORK = -1
      call SGETRI(m, Matrix, m, ipiv, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call SGETRI(m, Matrix, m, ipiv, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'getrif90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call getri(Matrix,ipiv)
      if (present(flop)) flop = flops_dgetri(mn_min)
   end subroutine sgetrif90

   subroutine zgetrif90(Matrix, ipiv, flop)
!
      implicit none
      complex(kind=8) Matrix(:, :)
      integer ipiv(:)
      integer m, n, mn_min

      integer LWORK, INFO
      complex(kind=8):: TEMP(1)
      real(kind=8), optional::flop
      complex(kind=8), allocatable:: WORK(:)

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      if (m > n) then
         write (*, *) 'size of Matrix incorrect'
         stop
      endif
      mn_min = min(m, n)

      LWORK = -1
      call ZGETRI(m, Matrix, m, ipiv, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call ZGETRI(m, Matrix, m, ipiv, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'getrif90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call getri(Matrix,ipiv)
      if (present(flop)) flop = flops_zgetri(mn_min)
   end subroutine zgetrif90

   subroutine cgetrif90(Matrix, ipiv, flop)
      !
      implicit none
      complex(kind=4) Matrix(:, :)
      integer ipiv(:)
      integer m, n, mn_min

      integer LWORK, INFO
      complex(kind=4):: TEMP(1)
      real(kind=8), optional::flop
      complex(kind=4), allocatable:: WORK(:)

      m = size(Matrix, 1)
      n = size(Matrix, 2)
      if (m > n) then
         write (*, *) 'size of Matrix incorrect'
         stop
      endif
      mn_min = min(m, n)

      LWORK = -1
      call CGETRI(m, Matrix, m, ipiv, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call CGETRI(m, Matrix, m, ipiv, WORK, LWORK, INFO)

      if (INFO /= 0) then
         write (*, *) 'getrif90 failed!!', INFO
         stop
      endif

      deallocate (WORK)

! call getri(Matrix,ipiv)
      if (present(flop)) flop = flops_zgetri(mn_min)
   end subroutine cgetrif90

   subroutine trsmf90(Matrix, B, side, uplo, transa, diag, m, n, flop)

      implicit none
      DT Matrix(:, :), B(:, :)
      integer m, n
      character side, uplo, transa, diag
      real(kind=8), optional::flop

#if DAT==0
      call ztrsmf90(Matrix, B, side, uplo, transa, diag, m, n, flop)
#elif DAT==1
      call dtrsmf90(Matrix, B, side, uplo, transa, diag, m, n, flop)
#elif DAT==2
      call ctrsmf90(Matrix, B, side, uplo, transa, diag, m, n, flop)
#elif DAT==3
      call strsmf90(Matrix, B, side, uplo, transa, diag, m, n, flop)
#endif

   end subroutine trsmf90

   subroutine dtrsmf90(Matrix, B, side, uplo, transa, diag, m, n, flop)
!
      implicit none
      real(kind=8) Matrix(:, :), B(:, :), alpha
      integer m, n, lda, ldb
      character side, uplo, transa, diag
      real(kind=8), optional::flop
      integer INFO

      alpha = 1d0

      ldb = size(B, 1)
      lda = size(Matrix, 1)

      call DTRSM(side, uplo, transa, diag, m, n, alpha, Matrix, lda, B, ldb)

! call trsm(Matrix,B,side,uplo,transa,diag)
      if (present(flop)) flop = flops_dtrsm(side, m, n)
   end subroutine dtrsmf90

   subroutine strsmf90(Matrix, B, side, uplo, transa, diag, m, n, flop)
      !
      implicit none
      real(kind=4) Matrix(:, :), B(:, :), alpha
      integer m, n, lda, ldb
      character side, uplo, transa, diag
      real(kind=8), optional::flop
      integer INFO

      alpha = 1d0

      ldb = size(B, 1)
      lda = size(Matrix, 1)

      call STRSM(side, uplo, transa, diag, m, n, alpha, Matrix, lda, B, ldb)

! call trsm(Matrix,B,side,uplo,transa,diag)
      if (present(flop)) flop = flops_dtrsm(side, m, n)
   end subroutine strsmf90

   subroutine ztrsmf90(Matrix, B, side, uplo, transa, diag, m, n, flop)
!
      implicit none
      complex(kind=8) Matrix(:, :), B(:, :), alpha
      integer m, n, lda, ldb
      character side, uplo, transa, diag
      real(kind=8), optional::flop
      integer INFO

      alpha = 1d0

      ldb = size(B, 1)
      lda = size(Matrix, 1)

! write(*,*)shape(Matrix),lda,shape(B),ldb
      call ZTRSM(side, uplo, transa, diag, m, n, alpha, Matrix, lda, B, ldb)

! call trsm(Matrix,B,side,uplo,transa,diag)
      if (present(flop)) flop = flops_ztrsm(side, m, n)
   end subroutine ztrsmf90

   subroutine ctrsmf90(Matrix, B, side, uplo, transa, diag, m, n, flop)
      !
      implicit none
      complex(kind=4) Matrix(:, :), B(:, :), alpha
      integer m, n, lda, ldb
      character side, uplo, transa, diag
      real(kind=8), optional::flop
      integer INFO

      alpha = 1d0

      ldb = size(B, 1)
      lda = size(Matrix, 1)

! write(*,*)shape(Matrix),lda,shape(B),ldb
      call CTRSM(side, uplo, transa, diag, m, n, alpha, Matrix, lda, B, ldb)

! call trsm(Matrix,B,side,uplo,transa,diag)
      if (present(flop)) flop = flops_ztrsm(side, m, n)
   end subroutine ctrsmf90

#ifdef HAVE_MKL
   subroutine gemm_batch_mkl(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size, flop)
      character*1::transa_array(:),transb_array(:)
      DT::alpha_array(:),beta_array(:)
      integer::group_size(:),m_array(:),n_array(:),k_array(:),lda_array(:),ldb_array(:),ldc_array(:)
      integer(c_intptr_t)::a_array(:),b_array(:),c_array(:)
      real(kind=8), optional::flop
      integer m,n,k,group_count,ii

#if DAT==0
      call zgemm_batch(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size)
#elif DAT==1
      call dgemm_batch(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size)
#elif DAT==2
      call cgemm_batch(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size)
#elif DAT==3
      call sgemm_batch(transa_array, transb_array, m_array, n_array, k_array, alpha_array, a_array, lda_array, b_array, ldb_array, beta_array, c_array, ldc_array, group_count, group_size)
#endif

      if (present(flop))then
      flop=0
      do ii=1,group_count
         m=m_array(ii)
         n=n_array(ii)
         k=k_array(ii)
         flop = flop + flops_gemm(m, n, k)
      enddo
      endif
   end subroutine gemm_batch_mkl
#endif

   subroutine gemmf90(MatA, lda, MatB, ldb, MatC, ldc, transa, transb, m, n, k, alpha, beta, flop)
      implicit none
      integer m, n, k, lda, ldb, ldc
      DT MatA(lda, *), MatB(ldb, *), MatC(ldc, *)
      DT alpha, beta
      character transa, transb
      real(kind=8), optional::flop

      call gemmf77(transa, transb, m, n, k, alpha, MatA, lda, MatB, ldb, beta, MatC, ldc)

! call gemmf90(MatA,MatB,MatC,transa,transb,alpha,beta)
      if (present(flop)) flop = flops_gemm(m, n, k)
   end subroutine gemmf90

   subroutine pun_or_mqrf90(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, flop)
      implicit none
      character side, trans
      integer m, n, k, ia, ja, ic, jc
      DT MatA(:, :), MatC(:, :), tau(:)
      integer desca(9), descc(9)
      real(kind=8), optional::flop


#if DAT==0
      call pzunmqrf90(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, flop)
#elif DAT==1
      call pdormqrf90(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, flop)
#elif DAT==2
      call pcunmqrf90(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, flop)
#elif DAT==3
      call psormqrf90(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, flop)
#endif



   end subroutine pun_or_mqrf90

   subroutine pdormqrf90(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, flop)
      implicit none
      character side, trans
      integer m, n, k, ia, ja, ic, jc
      real(kind=8) MatA(:, :), MatC(:, :), tau(:)
      integer desca(9), descc(9)
      real(kind=8), allocatable:: WORK(:)
      integer LWORK, INFO
      real(kind=8):: TEMP(1)
      real(kind=8), optional::flop
      LWORK = -1
      call pdormqr(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, TEMP, lwork, info)
      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      call pdormqr(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, WORK, lwork, info)
      deallocate (WORK)

      if (present(flop)) flop = flops_dunmqr(side, m, n, k)

   end subroutine pdormqrf90

   subroutine psormqrf90(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, flop)
      implicit none
      character side, trans
      integer m, n, k, ia, ja, ic, jc
      real(kind=4) MatA(:, :), MatC(:, :), tau(:)
      integer desca(9), descc(9)
      real(kind=4), allocatable:: WORK(:)
      integer LWORK, INFO
      real(kind=4):: TEMP(1)
      real(kind=8), optional::flop
      LWORK = -1
      call psormqr(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, TEMP, lwork, info)
      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      call psormqr(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, WORK, lwork, info)
      deallocate (WORK)

      if (present(flop)) flop = flops_dunmqr(side, m, n, k)

   end subroutine psormqrf90

   subroutine pzunmqrf90(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, flop)
      implicit none
      character side, trans
      integer m, n, k, ia, ja, ic, jc
      complex(kind=8) MatA(:, :), MatC(:, :), tau(:)
      integer desca(9), descc(9)
      complex(kind=8), allocatable:: WORK(:)
      integer LWORK, INFO
      complex(kind=8):: TEMP(1)
      real(kind=8), optional::flop

      LWORK = -1
      call pzunmqr(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, TEMP, lwork, info)
      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      call pzunmqr(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, WORK, lwork, info)

      deallocate (WORK)
      if (present(flop)) flop = flops_zunmqr(side, m, n, k)
   end subroutine pzunmqrf90

   subroutine pcunmqrf90(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, flop)
      implicit none
      character side, trans
      integer m, n, k, ia, ja, ic, jc
      complex(kind=4) MatA(:, :), MatC(:, :), tau(:)
      integer desca(9), descc(9)
      complex(kind=4), allocatable:: WORK(:)
      integer LWORK, INFO
      complex(kind=4):: TEMP(1)
      real(kind=8), optional::flop

      LWORK = -1
      call pcunmqr(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, TEMP, lwork, info)
      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      call pcunmqr(side, trans, m, n, k, MatA, ia, ja, desca, tau, MatC, ic, jc, descc, WORK, lwork, info)

      deallocate (WORK)
      if (present(flop)) flop = flops_zunmqr(side, m, n, k)
   end subroutine pcunmqrf90



   subroutine pgeqrff90(M, N, Matrix, ia, ja, desca, tau, flop)
      implicit none
      integer M, N, ia, ja
      DT Matrix(:, :), tau(:)
      integer rank
      integer desca(9)
      real(kind=8), optional::flop

#if DAT==0
      call pzgeqrff90(M, N, Matrix, ia, ja, desca, tau, flop)
#elif DAT==1
      call pdgeqrff90(M, N, Matrix, ia, ja, desca, tau, flop)
#elif DAT==2
      call pcgeqrff90(M, N, Matrix, ia, ja, desca, tau, flop)
#elif DAT==3
      call psgeqrff90(M, N, Matrix, ia, ja, desca, tau, flop)
#endif

   end subroutine pgeqrff90

   subroutine pzgeqrff90(M, N, Matrix, ia, ja, desca, tau, flop)
      implicit none

      integer M, N, ia, ja
      complex(kind=8):: Matrix(:, :), tau(:)
      integer rank
      integer desca(9)
      integer LWORK, INFO, ierr
      complex(kind=8):: TEMP(1)
      complex(kind=8), allocatable:: WORK(:)
      real(kind=8), optional::flop

      LWORK = -1
      call PZGEQRF(M, N, Matrix, ia, ja, desca, tau, TEMP, lwork, info)
      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      call PZGEQRF(M, N, Matrix, ia, ja, desca, tau, WORK, lwork, info)

      deallocate (WORK)

      if (present(flop)) flop = flops_zgeqpfmod(m, n, min(m, n))
   end subroutine pzgeqrff90

   subroutine pcgeqrff90(M, N, Matrix, ia, ja, desca, tau, flop)
      implicit none

      integer M, N, ia, ja
      complex(kind=4):: Matrix(:, :), tau(:)
      integer rank
      integer desca(9)
      integer LWORK, INFO, ierr
      complex(kind=4):: TEMP(1)
      complex(kind=4), allocatable:: WORK(:)
      real(kind=8), optional::flop

      LWORK = -1
      call PCGEQRF(M, N, Matrix, ia, ja, desca, tau, TEMP, lwork, info)
      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      call PCGEQRF(M, N, Matrix, ia, ja, desca, tau, WORK, lwork, info)

      deallocate (WORK)

      if (present(flop)) flop = flops_zgeqpfmod(m, n, min(m, n))
   end subroutine pcgeqrff90


   subroutine pdgeqrff90(M, N, Matrix, ia, ja, desca, tau, flop)
      implicit none

      integer M, N, ia, ja
      real(kind=8):: Matrix(:, :), tau(:)
      integer rank
      integer desca(9)
      integer LWORK, INFO, ierr
      real(kind=8):: TEMP(1)
      real(kind=8), allocatable:: WORK(:)
      real(kind=8), optional::flop

      LWORK = -1
      call PDGEQRF(M, N, Matrix, 1, 1, desca, tau, TEMP, lwork, info)
      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      call PDGEQRF(M, N, Matrix, 1, 1, desca, tau, WORK, lwork, info)

      deallocate (WORK)

      if (present(flop)) flop = flops_dgeqpfmod(m, n, min(m, n))
   end subroutine pdgeqrff90

   subroutine psgeqrff90(M, N, Matrix, ia, ja, desca, tau, flop)
      implicit none

      integer M, N, ia, ja
      real(kind=4):: Matrix(:, :), tau(:)
      integer rank
      integer desca(9)
      integer LWORK, INFO, ierr
      real(kind=4):: TEMP(1)
      real(kind=4), allocatable:: WORK(:)
      real(kind=8), optional::flop

      LWORK = -1
      call PSGEQRF(M, N, Matrix, 1, 1, desca, tau, TEMP, lwork, info)
      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      call PSGEQRF(M, N, Matrix, 1, 1, desca, tau, WORK, lwork, info)

      deallocate (WORK)

      if (present(flop)) flop = flops_dgeqpfmod(m, n, min(m, n))
   end subroutine psgeqrff90


   subroutine pgeqpfmodf90(M, N, Matrix, ia, ja, desca, ipiv, tau, JPERM, jpiv, rank, rtol, atol, flop)
      implicit none
      integer M, N, ia, ja
      DT Matrix(:, :), tau(:)
      integer ipiv(:), jpiv(:), JPERM(:)
      integer rank
      integer desca(9)
      real(kind=8)::rtol, atol
      real(kind=8), optional::flop

#if DAT==0
      call pzgeqpfmodf90(M, N, Matrix, ia, ja, desca, ipiv, tau, JPERM, jpiv, rank, rtol, atol, flop)
#elif DAT==1
      call pdgeqpfmodf90(M, N, Matrix, ia, ja, desca, ipiv, tau, JPERM, jpiv, rank, rtol, atol, flop)
#elif DAT==2
      call pcgeqpfmodf90(M, N, Matrix, ia, ja, desca, ipiv, tau, JPERM, jpiv, rank, real(rtol), real(atol), flop)
#elif DAT==3
      call psgeqpfmodf90(M, N, Matrix, ia, ja, desca, ipiv, tau, JPERM, jpiv, rank, real(rtol), real(atol), flop)
#endif


   end subroutine pgeqpfmodf90

   subroutine pzgeqpfmodf90(M, N, Matrix, ia, ja, desca, ipiv, tau, JPERM, jpiv, rank, rtol, atol, flop)
      implicit none
      integer M, N, ia, ja
      complex(kind=8) Matrix(:, :), tau(:)
      integer ipiv(:), jpiv(:), JPERM(:)
      integer rank
      integer desca(9)
      real(kind=8)::rtol, atol
      integer LWORK, LRWORK, INFO, ierr
      complex(kind=8):: TEMP(1)
      real(kind=8), allocatable::RWORK(:)
      complex(kind=8), allocatable:: WORK(:)
      real(kind=8):: RTEMP(1)
      real(kind=8), optional::flop

      LWORK = -1
      LRWORK = -1
      call PZGEQPFmod(M, N, Matrix, 1, 1, desca, ipiv, tau, TEMP, lwork, RTEMP, lrwork, info, JPERM, jpiv, rank, rtol, atol)
      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      lrwork = NINT(dble(RTEMP(1)*2.001 + 1))
      allocate (RWORK(lrwork))
      RWORK = 0
      call PZGEQPFmod(M, N, Matrix, 1, 1, desca, ipiv, tau, WORK, lwork, RWORK, lrwork, info, JPERM, jpiv, rank, rtol, atol)

      deallocate (WORK)
      deallocate (RWORK)

      if (present(flop)) flop = flops_zgeqpfmod(m, n, rank)
   end subroutine pzgeqpfmodf90

   subroutine pcgeqpfmodf90(M, N, Matrix, ia, ja, desca, ipiv, tau, JPERM, jpiv, rank, rtol, atol, flop)
      implicit none
      integer M, N, ia, ja
      complex(kind=4) Matrix(:, :), tau(:)
      integer ipiv(:), jpiv(:), JPERM(:)
      integer rank
      integer desca(9)
      real(kind=4)::rtol, atol
      integer LWORK, LRWORK, INFO, ierr
      complex(kind=4):: TEMP(1)
      real(kind=4), allocatable::RWORK(:)
      complex(kind=4), allocatable:: WORK(:)
      real(kind=4):: RTEMP(1)
      real(kind=8), optional::flop

      LWORK = -1
      LRWORK = -1
      call PCGEQPFmod(M, N, Matrix, 1, 1, desca, ipiv, tau, TEMP, lwork, RTEMP, lrwork, info, JPERM, jpiv, rank, rtol, atol)
      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      lrwork = NINT(dble(RTEMP(1)*2.001 + 1))
      allocate (RWORK(lrwork))
      RWORK = 0
      call PCGEQPFmod(M, N, Matrix, 1, 1, desca, ipiv, tau, WORK, lwork, RWORK, lrwork, info, JPERM, jpiv, rank, rtol, atol)

      deallocate (WORK)
      deallocate (RWORK)

      if (present(flop)) flop = flops_zgeqpfmod(m, n, rank)
   end subroutine pcgeqpfmodf90


   subroutine pdgeqpfmodf90(M, N, Matrix, ia, ja, desca, ipiv, tau, JPERM, jpiv, rank, rtol, atol, flop)
      implicit none
      integer M, N, ia, ja
      real(kind=8) Matrix(:, :), tau(:)
      integer ipiv(:), jpiv(:), JPERM(:)
      integer rank
      integer desca(9)
      real(kind=8)::rtol, atol
      integer LWORK, INFO, ierr
      real(kind=8):: TEMP(1)
      real(kind=8), allocatable:: WORK(:)

      real(kind=8), optional::flop
      LWORK = -1

      call PDGEQPFmod(M, N, Matrix, 1, 1, desca, ipiv, tau, TEMP, lwork, info, JPERM, jpiv, rank, rtol, atol)
      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      call PDGEQPFmod(M, N, Matrix, 1, 1, desca, ipiv, tau, WORK, lwork, info, JPERM, jpiv, rank, rtol, atol)

      deallocate (WORK)
      if (present(flop)) flop = flops_dgeqpfmod(m, n, rank)

   end subroutine pdgeqpfmodf90

   subroutine psgeqpfmodf90(M, N, Matrix, ia, ja, desca, ipiv, tau, JPERM, jpiv, rank, rtol, atol, flop)
      implicit none
      integer M, N, ia, ja
      real(kind=4) Matrix(:, :), tau(:)
      integer ipiv(:), jpiv(:), JPERM(:)
      integer rank
      integer desca(9)
      real(kind=4)::rtol, atol
      integer LWORK, INFO, ierr
      real(kind=4):: TEMP(1)
      real(kind=4), allocatable:: WORK(:)

      real(kind=8), optional::flop
      LWORK = -1

      call PSGEQPFmod(M, N, Matrix, 1, 1, desca, ipiv, tau, TEMP, lwork, info, JPERM, jpiv, rank, rtol, atol)
      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      call PSGEQPFmod(M, N, Matrix, 1, 1, desca, ipiv, tau, WORK, lwork, info, JPERM, jpiv, rank, rtol, atol)

      deallocate (WORK)
      if (present(flop)) flop = flops_dgeqpfmod(m, n, rank)

   end subroutine psgeqpfmodf90


   subroutine pgemr2df90(M, N, MatA, ia, ja, desca, MatB, ib, jb, descb, ictxt)
      implicit none
      integer M, N, ia, ja, ib, jb
      DT MatA(:, :), MatB(:, :)
      integer desca(9), descb(9)
      integer ictxt

#if DAT==0
      call pzgemr2d(M, N, MatA, ia, ja, desca, MatB, ib, jb, descb, ictxt)
#elif DAT==1
      call pdgemr2d(M, N, MatA, ia, ja, desca, MatB, ib, jb, descb, ictxt)
#elif DAT==2
      call pcgemr2d(M, N, MatA, ia, ja, desca, MatB, ib, jb, descb, ictxt)
#elif DAT==3
      call psgemr2d(M, N, MatA, ia, ja, desca, MatB, ib, jb, descb, ictxt)
#endif

   end subroutine pgemr2df90

   subroutine pgemmf90(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc, flop)
      implicit none
      character transa, transb
      integer m, n, k, ia, ja, ib, jb, ic, jc
      DT alpha, beta, a(:, :), b(:, :), c(:, :)
      integer desca(9), descb(9), descc(9)
      real(kind=8), optional::flop

#if DAT==0
      call pzgemm(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
      if (present(flop)) flop = flops_gemm(m, n, k)
#elif DAT==1
      call pdgemm(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
      if (present(flop)) flop = flops_gemm(m, n, k)
#elif DAT==2
      call pcgemm(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
      if (present(flop)) flop = flops_gemm(m, n, k)
#elif DAT==3
      call psgemm(transa, transb, m, n, k, alpha, a, ia, ja, desca, b, ib, jb, descb, beta, c, ic, jc, descc)
      if (present(flop)) flop = flops_gemm(m, n, k)
#endif

   end subroutine pgemmf90

   subroutine ptrsmf90(side, uplo, transa, diag, m, n, alpha, a, ia, ja, desca, b, ib, jb, descb, flop)
      implicit none
      character side, uplo, transa, diag
      integer m, n, ia, ja, ib, jb
      integer desca(9), descb(9)
      DT a(:, :), b(:, :), alpha
      real(kind=8), optional::flop

#if DAT==0
      call pztrsm(side, uplo, transa, diag, m, n, alpha, a, ia, ja, desca, b, ib, jb, descb)
      if (present(flop)) flop = flops_ztrsm(side, m, n)
#elif DAT==1
      call pdtrsm(side, uplo, transa, diag, m, n, alpha, a, ia, ja, desca, b, ib, jb, descb)
      if (present(flop)) flop = flops_dtrsm(side, m, n)
#elif DAT==2
      call pctrsm(side, uplo, transa, diag, m, n, alpha, a, ia, ja, desca, b, ib, jb, descb)
      if (present(flop)) flop = flops_ztrsm(side, m, n)
#elif DAT==3
      call pstrsm(side, uplo, transa, diag, m, n, alpha, a, ia, ja, desca, b, ib, jb, descb)
      if (present(flop)) flop = flops_dtrsm(side, m, n)
#endif

   end subroutine ptrsmf90

   subroutine pgetrff90(m, n, a, ia, ja, desca, ipiv, info, flop)
      implicit none
      integer m, n, ia, ja
      DT a(:, :)
      integer desca(9)
      integer ipiv(:)
      integer info
      real(kind=8), optional::flop

#if DAT==0
      call pzgetrf(m, n, a, ia, ja, desca, ipiv, info)
      if (present(flop)) flop = flops_zgetrf(m, n)
#elif DAT==1
      call pdgetrf(m, n, a, ia, ja, desca, ipiv, info)
      if (present(flop)) flop = flops_dgetrf(m, n)
#elif DAT==2
      call pcgetrf(m, n, a, ia, ja, desca, ipiv, info)
      if (present(flop)) flop = flops_zgetrf(m, n)
#elif DAT==3
      call psgetrf(m, n, a, ia, ja, desca, ipiv, info)
      if (present(flop)) flop = flops_dgetrf(m, n)
#endif

   end subroutine pgetrff90

   subroutine pgetrif90(n, a, ia, ja, desca, ipiv, flop)
      implicit none
      integer n, ia, ja
      DT:: a(:, :)
      integer desca(9)
      integer ipiv(:)
      real(kind=8), optional::flop

#if DAT==0
      call pzgetrif90(n, a, ia, ja, desca, ipiv, flop)
#elif DAT==1
      call pdgetrif90(n, a, ia, ja, desca, ipiv, flop)
#elif DAT==2
      call pcgetrif90(n, a, ia, ja, desca, ipiv, flop)
#elif DAT==3
      call psgetrif90(n, a, ia, ja, desca, ipiv, flop)
#endif


   end subroutine pgetrif90

   subroutine pdgetrif90(n, a, ia, ja, desca, ipiv, flop)
      implicit none
      integer n, ia, ja
      real(kind=8):: a(:, :)
      real(kind=8):: TEMP(1)
      integer desca(9)
      integer ipiv(:)
      integer info
      integer TEMPI(1)
      integer lwork, liwork
      integer, allocatable :: iwork(:)
      real(kind=8), allocatable:: work(:)
      real(kind=8), optional::flop

      call pdgetri(n, a, 1, 1, desca, ipiv, TEMP, -1, TEMPI, -1, info)
      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (work(lwork))
      work = 0
      liwork = TEMPI(1)
      allocate (iwork(liwork))
      iwork = 0
      call pdgetri(n, a, 1, 1, desca, ipiv, work, lwork, iwork, liwork, info)
      deallocate (iwork)
      deallocate (work)

      if (present(flop)) flop = flops_dgetri(n)

   end subroutine pdgetrif90

   subroutine psgetrif90(n, a, ia, ja, desca, ipiv, flop)
      implicit none
      integer n, ia, ja
      real(kind=4):: a(:, :)
      real(kind=4):: TEMP(1)
      integer desca(9)
      integer ipiv(:)
      integer info
      integer TEMPI(1)
      integer lwork, liwork
      integer, allocatable :: iwork(:)
      real(kind=4), allocatable:: work(:)
      real(kind=8), optional::flop

      call psgetri(n, a, 1, 1, desca, ipiv, TEMP, -1, TEMPI, -1, info)
      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (work(lwork))
      work = 0
      liwork = TEMPI(1)
      allocate (iwork(liwork))
      iwork = 0
      call psgetri(n, a, 1, 1, desca, ipiv, work, lwork, iwork, liwork, info)
      deallocate (iwork)
      deallocate (work)

      if (present(flop)) flop = flops_dgetri(n)

   end subroutine psgetrif90


   subroutine pzgetrif90(n, a, ia, ja, desca, ipiv, flop)
      implicit none
      integer n, ia, ja
      complex(kind=8):: a(:, :)
      complex(kind=8):: TEMP(1)
      integer desca(9)
      integer ipiv(:)
      integer info
      integer TEMPI(1)
      integer lwork, liwork
      integer, allocatable :: iwork(:)
      complex(kind=8), allocatable:: work(:)
      real(kind=8), optional::flop

      call pzgetri(n, a, 1, 1, desca, ipiv, TEMP, -1, TEMPI, -1, info)
      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (work(lwork))
      work = 0
      liwork = TEMPI(1)
      allocate (iwork(liwork))
      iwork = 0
      call pzgetri(n, a, 1, 1, desca, ipiv, work, lwork, iwork, liwork, info)

      deallocate (iwork)
      deallocate (work)

      if (present(flop)) flop = flops_zgetri(n)

   end subroutine pzgetrif90


   subroutine pcgetrif90(n, a, ia, ja, desca, ipiv, flop)
      implicit none
      integer n, ia, ja
      complex(kind=4):: a(:, :)
      complex(kind=4):: TEMP(1)
      integer desca(9)
      integer ipiv(:)
      integer info
      integer TEMPI(1)
      integer lwork, liwork
      integer, allocatable :: iwork(:)
      complex(kind=4), allocatable:: work(:)
      real(kind=8), optional::flop

      call pcgetri(n, a, 1, 1, desca, ipiv, TEMP, -1, TEMPI, -1, info)
      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (work(lwork))
      work = 0
      liwork = TEMPI(1)
      allocate (iwork(liwork))
      iwork = 0
      call pcgetri(n, a, 1, 1, desca, ipiv, work, lwork, iwork, liwork, info)

      deallocate (iwork)
      deallocate (work)

      if (present(flop)) flop = flops_zgetri(n)

   end subroutine pcgetrif90

   subroutine pgesvdf90(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, flop)
      implicit none

      character jobu, jobvt
      integer m, n, ia, ja, iu, ju, ivt, jvt
      DT:: a(:, :), u(:, :), vt(:, :)
      integer desca(9), descu(9), descvt(9)
      DTR:: s(:)
      real(kind=8), optional::flop

#if DAT==0
      call pzgesvdf90(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, flop)
#elif DAT==1
      call pdgesvdf90(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, flop)
#elif DAT==2
      call pcgesvdf90(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, flop)
#elif DAT==3
      call psgesvdf90(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, flop)
#endif


   end subroutine pgesvdf90

   subroutine pdgesvdf90(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, flop)
      implicit none

      character jobu, jobvt
      integer m, n, ia, ja, iu, ju, ivt, jvt
      real(kind=8):: a(:, :), u(:, :), vt(:, :)
      integer desca(9), descu(9), descvt(9)
      real(kind=8):: s(:)
      real(kind=8):: TEMP(1)
      integer LWORK, mnmax, mnmin
      real(kind=8), allocatable:: WORK(:)
      integer info
      real(kind=8), optional::flop
      mnmax = max(m, n)
      mnmin = min(m, n)

      lwork = -1
      call pdgesvd(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, TEMP, lwork, info)

      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      call pdgesvd(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, WORK, lwork, info)

      deallocate (WORK)
      if (present(flop)) flop = flops_dgesvd(m, n)

   end subroutine pdgesvdf90

   subroutine psgesvdf90(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, flop)
      implicit none

      character jobu, jobvt
      integer m, n, ia, ja, iu, ju, ivt, jvt
      real(kind=4):: a(:, :), u(:, :), vt(:, :)
      integer desca(9), descu(9), descvt(9)
      real(kind=4):: s(:)
      real(kind=4):: TEMP(1)
      integer LWORK, mnmax, mnmin
      real(kind=4), allocatable:: WORK(:)
      integer info
      real(kind=8), optional::flop
      mnmax = max(m, n)
      mnmin = min(m, n)

      lwork = -1
      call psgesvd(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, TEMP, lwork, info)

      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      call psgesvd(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, WORK, lwork, info)

      deallocate (WORK)
      if (present(flop)) flop = flops_dgesvd(m, n)

   end subroutine psgesvdf90


   subroutine pzgesvdf90(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, flop)
      implicit none

      character jobu, jobvt
      integer m, n, ia, ja, iu, ju, ivt, jvt
      complex(kind=8):: a(:, :), u(:, :), vt(:, :)
      integer desca(9), descu(9), descvt(9)
      real(kind=8):: s(:)
      complex(kind=8):: TEMP(1)
      integer LWORK, mnmax, mnmin
      complex(kind=8), allocatable:: WORK(:)
      real(kind=8), allocatable::RWORK(:)
      integer info
      real(kind=8), optional::flop

      mnmax = max(m, n)
      mnmin = min(m, n)

      allocate (rwork(1 + 4*mnmax))
      rwork = 0
      lwork = -1
      call pzgesvd(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, TEMP, lwork, rwork, info)

      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      call pzgesvd(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, WORK, lwork, rwork, info)

      deallocate (WORK, rwork)
      if (present(flop)) flop = flops_zgesvd(m, n)

   end subroutine pzgesvdf90

   subroutine pcgesvdf90(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, flop)
      implicit none

      character jobu, jobvt
      integer m, n, ia, ja, iu, ju, ivt, jvt
      complex(kind=4):: a(:, :), u(:, :), vt(:, :)
      integer desca(9), descu(9), descvt(9)
      real(kind=4):: s(:)
      complex(kind=4):: TEMP(1)
      integer LWORK, mnmax, mnmin
      complex(kind=4), allocatable:: WORK(:)
      real(kind=4), allocatable::RWORK(:)
      integer info
      real(kind=8), optional::flop

      mnmax = max(m, n)
      mnmin = min(m, n)

      allocate (rwork(1 + 4*mnmax))
      rwork = 0
      lwork = -1
      call pcgesvd(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, TEMP, lwork, rwork, info)

      lwork = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(lwork))
      WORK = 0
      call pcgesvd(jobu, jobvt, m, n, a, ia, ja, desca, s, u, iu, ju, descu, vt, ivt, jvt, descvt, WORK, lwork, rwork, info)

      deallocate (WORK, rwork)
      if (present(flop)) flop = flops_zgesvd(m, n)

   end subroutine pcgesvdf90


   subroutine pgeeigf90(Matrix, N, desca, eigval, Q, flops)
      implicit none
      integer N
      DTR norm1,norm0
      DTC eigval(:)
      DT Matrix(:, :),  Q(:,:), VL(1,1)
      integer desca(9),descp(9),desctau0(9),desctau1(9)
      real(kind=8), optional::flops
      real(kind=8)::flop,t1,t2
      integer nprow, npcol, myrow, mycol, myArows, myAcols, myAcolsTau0
      integer ctxt, ii, jj, myi, myj, IWORK(1), iproc, jproc, info
      DT, allocatable:: Matrix0(:,:),tau0(:), tau1(:)
      integer,allocatable:: IPIV(:)
      INTEGER,parameter:: BLOCK_CYCLIC_2D=1, CSRC_=8, CTXT_=2, DLEN_=9, DTYPE_=1,lld_=9, mb_=5, m_=3, nb_=6, n_=4, rsrc_=7

      if (present(flops))flops=0

      ctxt = desca(2)
      call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
      if (myrow /= -1 .and. mycol /= -1) then
      myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
      myAcols = numroc_wp(N, nbslpk, mycol, 0, npcol)
      allocate(Matrix0(max(1,myArows), max(1,myAcols)))
      if(myArows>0 .and. myAcols>0)Matrix0 = Matrix

      ! note the the use of myrow makes sure the vector is duplicated per process row
      CALL descset( desctau1, 1, N, 1, nbslpk, myrow, desca( csrc_ ), ctxt, 1)
      allocate(tau1(max(1,myAcols)))
      tau1=0

      ! note the the use of myrow makes sure the vector is duplicated per process row
      CALL descset( desctau0, 1, N-1, 1, nbslpk, myrow, desca( csrc_ ), ctxt, 1)
      myAcolsTau0 = numroc_wp(N-1, nbslpk, mycol, 0, npcol)
      allocate(tau0(max(1,myAcolsTau0)))
      tau0=0

      ! compute Heisenberg reduction Matrix=Q H Q^h
      t1 = MPI_Wtime()
      call pgehrdf90(N, 1, N, Matrix, 1, 1, desca, tau0, flop)
      if (present(flops))flops=flops+flop
      t2 = MPI_Wtime()
      if(myrow==0 .and. mycol==0)write(*,*)"pgehrdf90 time:",t2-t1
      if(myArows>0 .and. myAcols>0)Q = Matrix

      do jj=1,N
      do ii=jj+2,N
      call g2l(jj, N, npcol, nbslpk, jproc, myj)
      call g2l(ii, N, nprow, nbslpk, iproc, myi)
      if(jproc==mycol .and. iproc==myrow)then
            Matrix(myi,myj)=0
      endif
      enddo
      enddo

      ! norm1 = pfnorm(N, N, Q, 1, 1, desca, '1')
      ! write(*,*)'Q norm0:', norm1
      ! write(*,*)abs(Q(1,:))


      ! shift right Q by 1 column and zero out the upper triangular part to call pun_or_gqrf90
      allocate(IPIV(max(1,myAcols)+nbslpk))
      IPIV=0
      do myj=1,myAcols
         call l2g(myj, mycol, N, npcol, nbslpk, jj)
         if(jj==1)then
            IPIV(myj)=1
         else
            IPIV(myj)=jj-1
         endif
      enddo

      ! note the the use of myrow makes sure the vector is duplicated per process row
      CALL descset( descp, 1, desca( n_ ) + desca( nb_ )*npcol, 1, nbslpk, myrow, desca( csrc_ ), ctxt, 1)

      t1 = MPI_Wtime()
#if DAT==0
      call pzlapiv('B', 'C', 'R', N, N, Q, 1, 1, desca, IPIV, 1, 1, descp, IWORK)
#elif DAT==1
      call pdlapiv('B', 'C', 'R', N, N, Q, 1, 1, desca, IPIV, 1, 1, descp, IWORK)
#elif DAT==2
      call pclapiv('B', 'C', 'R', N, N, Q, 1, 1, desca, IPIV, 1, 1, descp, IWORK)
#elif DAT==3
      call pslapiv('B', 'C', 'R', N, N, Q, 1, 1, desca, IPIV, 1, 1, descp, IWORK)
#endif

      ! norm1 = pfnorm(N, N, Q, 1, 1, desca, '1')
      ! write(*,*)'Q norm1:', norm1
      ! write(*,*) abs(Q(1,:))

      do jj=1,N
      do ii=1,jj
      call g2l(jj, N, npcol, nbslpk, jproc, myj)
      call g2l(ii, N, nprow, nbslpk, iproc, myi)
      if(jproc==mycol .and. iproc==myrow)then
         if(ii==jj)then
            Q(myi,myj)=1d0
         else
            Q(myi,myj)=0
         endif
      endif
      enddo
      enddo


#if DAT==0
      call pzcopy(N-1,tau0,1,1,desctau0,1, tau1, 1, 2, desctau1,1)
#elif DAT==1
      call pdcopy(N-1,tau0,1,1,desctau0,1, tau1, 1, 2, desctau1,1)
#elif DAT==2
      call pccopy(N-1,tau0,1,1,desctau0,1, tau1, 1, 2, desctau1,1)
#elif DAT==3
      call pscopy(N-1,tau0,1,1,desctau0,1, tau1, 1, 2, desctau1,1)
#endif
      t2 = MPI_Wtime()
      if(myrow==0 .and. mycol==0)write(*,*)"pxlapiv+pxcopy time:",t2-t1

      ! write(*,*)myArows, myAcols,desca
      ! norm1 = pfnorm(N, N, Q, 1, 1, desca, '1')
      ! write(*,*)'Q norm2:', norm1, myArows, myAcols, abs(tau1)


      ! form Q from the Heisenberg reduction Matrix=Q H Q^h
      t1 = MPI_Wtime()
      call pun_or_gqrf90(ctxt, Q, tau1, N, N, N, desca, 1, 1, flop)
      if (present(flops))flops=flops+flop
      t2 = MPI_Wtime()
      if(myrow==0 .and. mycol==0)write(*,*)"pun_or_gqrf90 time:",t2-t1


      ! norm1 = pfnorm(N, N, Q, 1, 1, desca, '1')
      ! write(*,*)'Q norm3:', norm1

      ! compute the Schur decomposition H=S T S^h, and return Q=QS
      t1 = MPI_Wtime()
      call plahqrf90(.true., .true., N, 1, N, Matrix, desca, eigval, 1, N, Q, flop)
      if (present(flops))flops=flops+flop
      t2 = MPI_Wtime()
      if(myrow==0 .and. mycol==0)write(*,*)"plahqrf90 time:",t2-t1

      ! norm1 = pfnorm(N, N, Q, 1, 1, desca, '1')
      ! write(*,*)'Q norm4:', norm1

      ! compute the eigen decomposition T=Z D Z^h, and return Z=QZ
      t1 = MPI_Wtime()
      call ptrevcf90('R', N, Matrix, desca, VL, desca, Q, desca,flop)
      if (present(flops))flops=flops+flop
      t2 = MPI_Wtime()
      if(myrow==0 .and. mycol==0)write(*,*)"ptrevcf90 time:",t2-t1
      ! norm1 = pfnorm(N, N, Q, 1, 1, desca, '1')
      ! write(*,*)'Q norm5:', norm1

! check the error of the eigen decomposition
#if 1
      Matrix=Q
      norm0 = pfnorm(N, N, Matrix0, 1, 1, desca, '1')
      do jj=1,N
         call g2l(jj, N, npcol, nbslpk, jproc, myj)
         if(jproc==mycol)then
            Matrix(:,myj) = Matrix(:,myj)*eigval(jj)
         endif
      enddo
      call pgemmf90('N', 'N', N, N, N, BPACK_cone, Matrix0, 1, 1, desca, Q, 1, 1, desca, -BPACK_cone, Matrix, 1, 1, desca)
      norm1 = pfnorm(N, N, Matrix, 1, 1, desca, '1')
      if(0==mycol .and. 0==myrow)then
      write(*,*)'Error of pgeeigf90: ', norm1/norm0
      write(*,*)'Max and min eigenvalues: ', maxval(abs(eigval)), minval(abs(eigval))
      endif
#endif



      Matrix = Matrix0
      deallocate(Matrix0)
      deallocate(IPIV)
      deallocate(tau0)
      deallocate(tau1)
endif
   end subroutine pgeeigf90



   ! input matrix is upper triangular, it's modified but is restored on return
   subroutine ptrevcf90(SIDE, N, Matrix, DESCT, VL, DESCVL,VR, DESCVR,flop)


      implicit none
      character side
      integer N
      DT Matrix(:, :), VL(:,:), VR(:,:)
      integer DESCT(9),DESCVL(9),DESCVR(9)
      real(kind=8), optional::flop

#if DAT==0
      call pztrevcf90(SIDE, N, Matrix, DESCT, VL, DESCVL,VR, DESCVR,flop)
#elif DAT==1
      write(*,*)"pdtrevc hasn't been implemented in scalapack"
#elif DAT==2
      call pctrevcf90(SIDE, N, Matrix, DESCT, VL, DESCVL,VR, DESCVR,flop)
#elif DAT==3
      write(*,*)"pstrevc hasn't been implemented in scalapack"
#endif

   end subroutine ptrevcf90


   subroutine pztrevcf90(SIDE, N, Matrix, DESCT, VL, DESCVL,VR, DESCVR,flop)
      implicit none
      character side, howmny
      integer N,MM,M
      logical select(N)
      complex(kind=8) Matrix(:, :), VL(:,:), VR(:,:)
      complex(kind=8) WORK(2*N)
      real(kind=8) RWORK(N)
      integer DESCT(9),DESCVL(9),DESCVR(9)
      real(kind=8), optional::flop
      integer INFO

      howmny='B'
      select=.true.
      M=N
      MM=N

      call pztrevc( SIDE, HOWMNY, SELECT, N, Matrix, DESCT, VL, DESCVL, VR, DESCVR, MM, M, WORK, RWORK, INFO )
      if (present(flop)) flop = 0 ! LOT typically ignored
   end subroutine pztrevcf90


   subroutine pctrevcf90(SIDE, N, Matrix, DESCT, VL, DESCVL,VR, DESCVR,flop)
      implicit none
      character side, howmny
      integer N,MM,M
      logical select(N)
      complex(kind=4) Matrix(:, :), VL(:,:), VR(:,:)
      complex(kind=4) WORK(2*N)
      real(kind=4) RWORK(N)
      integer DESCT(9),DESCVL(9),DESCVR(9)
      real(kind=8), optional::flop
      integer INFO

      howmny='B'
      select=.true.
      M=N
      MM=N

      call pctrevc( SIDE, HOWMNY, SELECT, N, Matrix, DESCT, VL, DESCVL, VR, DESCVR, MM, M, WORK, RWORK, INFO )
      if (present(flop)) flop = 0 ! LOT typically ignored
   end subroutine pctrevcf90


   subroutine pgehrdf90(N, ilo, ihi, Matrix, ia, ja, desca, tau, flop)
      implicit none
      integer N, ilo, ihi
      DT Matrix(:, :),tau(:)
      integer desca(9)
      integer ia,ja
      real(kind=8), optional::flop

#if DAT==0
      call pzgehrdf90(N, ilo, ihi, Matrix, ia, ja, desca, tau, flop)
#elif DAT==1
      call pdgehrdf90(N, ilo, ihi, Matrix, ia, ja, desca, tau, flop)
#elif DAT==2
      call pcgehrdf90(N, ilo, ihi, Matrix, ia, ja, desca, tau, flop)
#elif DAT==3
      call psgehrdf90(N, ilo, ihi, Matrix, ia, ja, desca, tau, flop)
#endif

   end subroutine pgehrdf90


   subroutine pdgehrdf90(N, ilo, ihi, Matrix, ia, ja, desca, tau, flop)
      implicit none
      integer N, ilo, ihi
      real(kind=8) Matrix(:, :),tau(:)
      integer desca(9)
      integer ia,ja
      real(kind=8), optional::flop

      integer LWORK, INFO
      real(kind=8):: TEMP(1)
      real(kind=8), allocatable:: WORK(:)


      LWORK = -1
      call PDGEHRD(N, ilo, ihi, Matrix, ia, ja, desca, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call PDGEHRD(N, ilo, ihi, Matrix, ia, ja, desca, tau, WORK, LWORK, INFO)

      deallocate(WORK)
      if (present(flop)) flop = 10.0/3.0*N**3d0

   end subroutine pdgehrdf90

   subroutine psgehrdf90(N, ilo, ihi, Matrix, ia, ja, desca, tau, flop)
      implicit none
      integer N, ilo, ihi
      real(kind=4) Matrix(:, :),tau(:)
      integer desca(9)
      integer ia,ja
      real(kind=8), optional::flop

      integer LWORK, INFO
      real(kind=4):: TEMP(1)
      real(kind=4), allocatable:: WORK(:)


      LWORK = -1
      call PSGEHRD(N, ilo, ihi, Matrix, ia, ja, desca, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call PSGEHRD(N, ilo, ihi, Matrix, ia, ja, desca, tau, WORK, LWORK, INFO)

      deallocate(WORK)
      if (present(flop)) flop = 10.0/3.0*N**3d0

   end subroutine psgehrdf90

   subroutine pzgehrdf90(N, ilo, ihi, Matrix, ia, ja, desca, tau, flop)
      implicit none
      integer N, ilo, ihi
      complex(kind=8) Matrix(:, :),tau(:)
      integer desca(9)
      integer ia,ja
      real(kind=8), optional::flop

      integer LWORK, INFO
      complex(kind=8):: TEMP(1)
      complex(kind=8), allocatable:: WORK(:)


      LWORK = -1
      call PZGEHRD(N, ilo, ihi, Matrix, ia, ja, desca, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call PZGEHRD(N, ilo, ihi, Matrix, ia, ja, desca, tau, WORK, LWORK, INFO)

      deallocate(WORK)
      if (present(flop)) flop = 6.0* 10.0/3.0*N**3d0

   end subroutine pzgehrdf90

   subroutine pcgehrdf90(N, ilo, ihi, Matrix, ia, ja, desca, tau, flop)
      implicit none
      integer N, ilo, ihi
      complex(kind=4) Matrix(:, :),tau(:)
      integer desca(9)
      integer ia,ja
      real(kind=8), optional::flop

      integer LWORK, INFO
      complex(kind=4):: TEMP(1)
      complex(kind=4), allocatable:: WORK(:)


      LWORK = -1
      call PCGEHRD(N, ilo, ihi, Matrix, ia, ja, desca, tau, TEMP, LWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call PCGEHRD(N, ilo, ihi, Matrix, ia, ja, desca, tau, WORK, LWORK, INFO)

      deallocate(WORK)
      if (present(flop)) flop = 6.0* 10.0/3.0*N**3d0

   end subroutine pcgehrdf90

   subroutine plahqrf90(wantt, wantz, N, ilo, ihi, Matrix, desca, eigval, iloz, ihiz, Z, flop)
      implicit none
      logical wantt, wantz
      integer N, ilo, ihi, iloz, ihiz
      DT Matrix(:, :), Z(:,:)
      integer desca(9)
      DTC eigval(:)
      real(kind=8), optional::flop

#if DAT==0
      call pzlahqrf90(wantt, wantz, N, ilo, ihi, Matrix, desca, eigval, iloz, ihiz, Z, flop)
#elif DAT==1
      call pdlahqrf90(wantt, wantz, N, ilo, ihi, Matrix, desca, eigval, iloz, ihiz, Z, flop)
#elif DAT==2
      call pclahqrf90(wantt, wantz, N, ilo, ihi, Matrix, desca, eigval, iloz, ihiz, Z, flop)
#elif DAT==3
      call pslahqrf90(wantt, wantz, N, ilo, ihi, Matrix, desca, eigval, iloz, ihiz, Z, flop)
#endif

   end subroutine plahqrf90


   subroutine pslahqrf90(wantt, wantz, N, ilo, ihi, Matrix, desca, eigval, iloz, ihiz, Z, flop)
      implicit none
      logical wantt, wantz
      integer N, ilo, ihi, iloz, ihiz
      real(kind=4) Matrix(:, :), Z(:,:), WR(N), WI(N)
      integer desca(9)
      complex(kind=4) eigval(:)

      integer LWORK, ILWORK, INFO
      real(kind=4):: TEMP(1)
      integer::IWORK(1) ! not referenced, place holder for a future release
      real(kind=4), allocatable:: WORK(:)

      real(kind=8), optional::flop

      LWORK = -1
      ILWORK = 1
      call PSLAHQR( WANTT, WANTZ, N, ILO, IHI, Matrix, DESCA, WR, WI, ILOZ, IHIZ, Z, DESCA, TEMP, LWORK, IWORK, ILWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call PSLAHQR( WANTT, WANTZ, N, ILO, IHI, Matrix, DESCA, WR, WI, ILOZ, IHIZ, Z, DESCA, WORK, LWORK, IWORK, ILWORK, INFO)

      eigval = WR + WI*BPACK_junit

      deallocate (WORK)
      if (present(flop)) flop = N**3d0 ! very rough estimate as the QR iteration convergence is unknown in prior

   end subroutine pslahqrf90


   subroutine pdlahqrf90(wantt, wantz, N, ilo, ihi, Matrix, desca, eigval, iloz, ihiz, Z, flop)
      implicit none
      logical wantt, wantz
      integer N, ilo, ihi, iloz, ihiz
      real(kind=8) Matrix(:, :), Z(:,:), WR(N), WI(N)
      integer desca(9)
      complex(kind=8) eigval(:)

      integer LWORK, ILWORK, INFO
      real(kind=8):: TEMP(1)
      integer::IWORK(1) ! not referenced, place holder for a future release
      real(kind=8), allocatable:: WORK(:)

      real(kind=8), optional::flop

      LWORK = -1
      ILWORK = 1
      call PDLAHQR( WANTT, WANTZ, N, ILO, IHI, Matrix, DESCA, WR, WI, ILOZ, IHIZ, Z, DESCA, TEMP, LWORK, IWORK, ILWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call PDLAHQR( WANTT, WANTZ, N, ILO, IHI, Matrix, DESCA, WR, WI, ILOZ, IHIZ, Z, DESCA, WORK, LWORK, IWORK, ILWORK, INFO)

      eigval = WR + WI*BPACK_junit

      deallocate (WORK)
      if (present(flop)) flop = N**3d0 ! very rough estimate as the QR iteration convergence is unknown in prior

   end subroutine pdlahqrf90



   subroutine pzlahqrf90(wantt, wantz, N, ilo, ihi, Matrix, desca, eigval, iloz, ihiz, Z, flop)
      implicit none
      logical wantt, wantz
      integer N, ilo, ihi, iloz, ihiz
      complex(kind=8) Matrix(:, :), Z(:,:)
      integer desca(9)
      complex(kind=8) eigval(:)

      integer LWORK, ILWORK, INFO
      complex(kind=8):: TEMP(1)
      integer::IWORK(1) ! not referenced, place holder for a future release
      complex(kind=8), allocatable:: WORK(:)

      real(kind=8), optional::flop

      LWORK = -1
      ILWORK = 1
      call PZLAHQR( WANTT, WANTZ, N, ILO, IHI, Matrix, DESCA, eigval, ILOZ, IHIZ, Z, DESCA, TEMP, LWORK, IWORK, ILWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call PZLAHQR( WANTT, WANTZ, N, ILO, IHI, Matrix, DESCA, eigval, ILOZ, IHIZ, Z, DESCA, WORK, LWORK, IWORK, ILWORK, INFO)

      deallocate (WORK)
      if (present(flop)) flop = 4.*N**3d0 ! very rough estimate as the QR iteration convergence is unknown in prior

   end subroutine pzlahqrf90

   subroutine pclahqrf90(wantt, wantz, N, ilo, ihi, Matrix, desca, eigval, iloz, ihiz, Z, flop)
      implicit none
      logical wantt, wantz
      integer N, ilo, ihi, iloz, ihiz
      complex(kind=4) Matrix(:, :), Z(:,:)
      integer desca(9)
      complex(kind=4) eigval(:)

      integer LWORK, ILWORK, INFO
      complex(kind=4):: TEMP(1)
      integer::IWORK(1) ! not referenced, place holder for a future release
      complex(kind=4), allocatable:: WORK(:)

      real(kind=8), optional::flop

      LWORK = -1
      ILWORK = 1
      call PCLAHQR( WANTT, WANTZ, N, ILO, IHI, Matrix, DESCA, eigval, ILOZ, IHIZ, Z, DESCA, TEMP, LWORK, IWORK, ILWORK, INFO)

      LWORK = NINT(dble(TEMP(1)*2.001 + 1))
      allocate (WORK(LWORK))
      WORK = 0
      call PCLAHQR( WANTT, WANTZ, N, ILO, IHI, Matrix, DESCA, eigval, ILOZ, IHIZ, Z, DESCA, WORK, LWORK, IWORK, ILWORK, INFO)

      deallocate (WORK)
      if (present(flop)) flop = 4.*N**3d0 ! very rough estimate as the QR iteration convergence is unknown in prior

   end subroutine pclahqrf90

   real(kind=8) function flops_zgesdd(m, n)
      implicit none
      integer m, n
      if (m > n) then
         flops_zgesdd = 4.*(8.*m*n*n + 4./3.*n*n*n)
      else
         flops_zgesdd = 4.*(8.*n*m*m + 4./3.*m*m*m)
      endif
   end function flops_zgesdd

   real(kind=8) function flops_dgesdd(m, n)
      implicit none
      integer m, n
      if (m > n) then
         flops_dgesdd = 8.*m*n*n + 4./3.*n*n*n
      else
         flops_dgesdd = 8.*n*m*m + 4./3.*m*m*m
      endif
   end function flops_dgesdd

   real(kind=8) function flops_zgesvd(m, n)
      implicit none
      integer m, n
      if (m > n) then
         flops_zgesvd = 4.*(12.*m*n*n + 16./3.*n*n*n)
      else
         flops_zgesvd = 4.*(12.*n*m*m + 16./3.*m*m*m)
      endif
   end function flops_zgesvd

   real(kind=8) function flops_dgesvd(m, n)
      implicit none
      integer m, n
      if (m > n) then
         flops_dgesvd = 12.*m*n*n + 16./3.*n*n*n
      else
         flops_dgesvd = 12.*n*m*m + 16./3.*m*m*m
      endif
   end function flops_dgesvd

   real(kind=8) function flops_dgeqpfmod(m, n, k)
      implicit none
      integer m, n, k
      if (m > n) then
         flops_dgeqpfmod = 2.*m*n*n - 2./3.*n*n*n - (2.*(m - k)*(n - k)*(n - k) - 2./3.*(n - k)*(n - k)*(n - k))
      else
         flops_dgeqpfmod = 2.*n*m*m - 2./3.*m*m*m - (2.*(n - k)*(m - k)*(m - k) - 2./3.*(m - k)*(m - k)*(m - k))
      endif
   end function flops_dgeqpfmod

   real(kind=8) function flops_zgeqpfmod(m, n, k)
      implicit none
      integer m, n, k
      if (m > n) then
         flops_zgeqpfmod = 4.*(2.*m*n*n - 2./3.*n*n*n - (2.*(m - k)*(n - k)*(n - k) - 2./3.*(n - k)*(n - k)*(n - k)))
      else
         flops_zgeqpfmod = 4.*(2.*n*m*m - 2./3.*m*m*m - (2.*(n - k)*(m - k)*(m - k) - 2./3.*(m - k)*(m - k)*(m - k)))
      endif
   end function flops_zgeqpfmod

   real(kind=8) function fmuls_geqrf(m, n)
      implicit none
      integer m, n
      if (m > n) then
         fmuls_geqrf = m*n*n - 1./3.*n*n*n + m*n + 0.5*n*n + 23./6.*n
      else
         fmuls_geqrf = n*m*m - 1./3.*m*m*m + 2*n*m - 0.5*m*m + 23./6.*m
      endif
   end function fmuls_geqrf
   real(kind=8) function fadds_geqrf(m, n)
      implicit none
      integer m, n
      if (m > n) then
         fadds_geqrf = m*n*n - 1./3.*n*n*n + 0.5*n*n + 5./6.*n
      else
         fadds_geqrf = n*m*m - 1./3.*m*m*m + n*m - 0.5*m*m + 5./6.*m
      endif
   end function fadds_geqrf
   real(kind=8) function flops_zgeqrf(m, n)
      implicit none
      integer m, n
      flops_zgeqrf = 6.*fmuls_geqrf(m, n) + 2.*fadds_geqrf(m, n)
   end function flops_zgeqrf
   real(kind=8) function flops_dgeqrf(m, n)
      implicit none
      integer m, n
      flops_dgeqrf = fmuls_geqrf(m, n) + fadds_geqrf(m, n)
   end function flops_dgeqrf

   real(kind=8) function fmuls_ungqr(m, n, k)
      implicit none
      integer m, n, k
      fmuls_ungqr = 2.*m*n*k - (m + n)*k*k + 2./3.*k*k*k + 2.*n*k - k*k - 5./3.*k
   end function fmuls_ungqr
   real(kind=8) function fadds_ungqr(m, n, k)
      implicit none
      integer m, n, k
      fadds_ungqr = 2.*m*n*k - (m + n)*k*k + 2./3.*k*k*k + n*k - m*k + 1./3.*k
   end function fadds_ungqr
   real(kind=8) function flops_zungqr(m, n, k)
      implicit none
      integer m, n, k
      flops_zungqr = 6.*fmuls_ungqr(m, n, k) + 2.*fadds_ungqr(m, n, k)
   end function flops_zungqr
   real(kind=8) function flops_dungqr(m, n, k)
      implicit none
      integer m, n, k
      flops_dungqr = fmuls_ungqr(m, n, k) + fadds_ungqr(m, n, k)
   end function flops_dungqr

   real(kind=8) function fmuls_unmqr(side, m, n, k)
      integer m, n, k
      character side
      if (side == 'L') then
         fmuls_unmqr = 2.*n*m*k - n*k*k + 2.*n*k
      else
         fmuls_unmqr = 2.*n*m*k - m*k*k + m*k + n*k - 0.5*k*k + 0.5*k
      endif
   end function fmuls_unmqr

   real(kind=8) function fadds_unmqr(side, m, n, k)
      integer m, n, k
      character side
      if (side == 'L') then
         fadds_unmqr = 2.*n*m*k - n*k*k + n*k
      else
         fadds_unmqr = 2.*n*m*k - m*k*k + m*k
      endif
   end function fadds_unmqr

   real(kind=8) function flops_zunmqr(side, m, n, k)
      integer m, n, k
      character side
      flops_zunmqr = 6.*fmuls_unmqr(side, m, n, k) + 2.*fadds_unmqr(side, m, n, k)
   end function flops_zunmqr

   real(kind=8) function flops_dunmqr(side, m, n, k)
      integer m, n, k
      character side
      flops_dunmqr = fmuls_unmqr(side, m, n, k) + fadds_unmqr(side, m, n, k)
   end function flops_dunmqr

   real(kind=8) function fmuls_getrf(m, n)
      implicit none
      integer m, n

      if (m > n) then
         fmuls_getrf = 0.5*m*n*n - 1./6.*n*n*n + 0.5*m*n - 0.5*n*n + 2./3.*n
      else
         fmuls_getrf = 0.5*n*m*m - 1./6.*m*m*m + 0.5*n*m - 0.5*m*m + 2./3.*m
      endif
   end function fmuls_getrf
   real(kind=8) function fadds_getrf(m, n)
      implicit none
      integer m, n

      if (m > n) then
         fadds_getrf = 0.5*m*n*n - 1./6.*n*n*n - 0.5*m*n + 1./6.*n
      else
         fadds_getrf = 0.5*n*m*m - 1./6.*m*m*m - 0.5*n*m + 1./6.*m
      endif
   end function fadds_getrf
   real(kind=8) function flops_zgetrf(m, n)
      implicit none
      integer m, n
      flops_zgetrf = 6.*fmuls_getrf(m, n) + 2.*fadds_getrf(m, n)
   end function flops_zgetrf
   real(kind=8) function flops_dgetrf(m, n)
      implicit none
      integer m, n
      flops_dgetrf = fmuls_getrf(m, n) + fadds_getrf(m, n)
   end function flops_dgetrf

   real(kind=8) function fmuls_getrs(n, nrhs)
      implicit none
      integer n, nrhs
      fmuls_getrs = nrhs*n*n
   end function fmuls_getrs
   real(kind=8) function fadds_getrs(n, nrhs)
      implicit none
      integer n, nrhs
      fadds_getrs = nrhs*n*(n - 1)
   end function fadds_getrs
   real(kind=8) function flops_zgetrs(n, nrhs)
      implicit none
      integer n, nrhs
      flops_zgetrs = 6.*fmuls_getrs(n, nrhs) + 2.*fadds_getrs(n, nrhs)
   end function flops_zgetrs
   real(kind=8) function flops_dgetrs(n, nrhs)
      implicit none
      integer n, nrhs
      flops_dgetrs = fmuls_getrs(n, nrhs) + fadds_getrs(n, nrhs)
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
      integer m, n
      character side
      if (side == 'L') then
         fmuls_trsm = 0.5*n*m*(m + 1)
      elseif (side == 'R') then
         fmuls_trsm = 0.5*m*n*(n + 1)
      endif
   end function fmuls_trsm
   real(kind=8) function fadds_trsm(side, m, n)
      integer m, n
      character side
      if (side == 'L') then
         fadds_trsm = 0.5*n*m*(m - 1)
      elseif (side == 'R') then
         fadds_trsm = 0.5*m*n*(n - 1)
      endif
   end function fadds_trsm
   real(kind=8) function flops_ztrsm(side, m, n)
      integer m, n
      character side
      flops_ztrsm = 6.*fmuls_trsm(side, m, n) + 2.*fadds_trsm(side, m, n)
   end function flops_ztrsm
   real(kind=8) function flops_dtrsm(side, m, n)
      integer m, n
      character side
      flops_dtrsm = fmuls_trsm(side, m, n) + fadds_trsm(side, m, n)
   end function flops_dtrsm

   real(kind=8) function fmuls_gemm(m, n, k)
      implicit none
      integer m, n, k
      fmuls_gemm = m*n*k
   end function fmuls_gemm
   real(kind=8) function fadds_gemm(m, n, k)
      implicit none
      integer m, n, k
      fadds_gemm = m*n*k
   end function fadds_gemm
   real(kind=8) function flops_gemm(m, n, k)
      implicit none
      integer m, n, k
#if DAT==0
   flops_gemm = 6.*fmuls_gemm(m, n, k) + 2.*fadds_gemm(m, n, k)
#elif DAT==1
   flops_gemm = fmuls_gemm(m, n, k) + fadds_gemm(m, n, k)
#elif DAT==2
   flops_gemm = 6.*fmuls_gemm(m, n, k) + 2.*fadds_gemm(m, n, k)
#elif DAT==3
   flops_gemm = fmuls_gemm(m, n, k) + fadds_gemm(m, n, k)
#endif

   end function flops_gemm

end module MISC_DenseLA
