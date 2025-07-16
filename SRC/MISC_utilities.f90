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

!> @file MISC_utilities.f90
!> @brief Low-level miscellaneous subroutines and functions

#include "ButterflyPACK_config.fi"

! ! ! ! ! ! ! the following can be commented out to avoid multiple definition for intel compilers, which means for now I'm not using VSL at all.
! #ifdef HAVE_MKL
! #ifndef USEVSL
! #define USEVSL
! include "mkl_vsl.f90"
! #endif
! #endif

module MISC_Utilities
   use BPACK_DEFS
   use MISC_DenseLA
   use BPACK_linkedlist

#ifdef Intel
   USE IFPORT
#endif
#ifdef HAVE_OPENMP
   use omp_lib
#endif


   integer, parameter :: int64 = selected_int_kind(18)

contains

! #ifndef mymacro(x)
! #define mymacro(x) print *, "Now giving information about ", "x" ; \
   ! call mysub( x, size(x,1), size(x,2) ) ; \
   ! print *, "About to do function on ", "x"; \
   ! call dofunction(x) ; \
   ! print *, "x"," is a nice matrix! Huzzah!"
! #endif

   function seq_wtime ( )

   implicit none

   integer count
   integer count_max
   integer count_rate
   real ( kind = 8 ) seq_wtime

   call system_clock ( count, count_rate, count_max )

   seq_wtime = real ( count, kind = 8 ) / real ( count_rate, kind = 8 )

   return
   end function seq_wtime



   subroutine linspaceI(startI, endI, N, array)
      implicit none
      integer startI, endI, N
      integer array(1:N)
      real(kind=8)::rtemp
      integer i

      rtemp = dble(endI - startI + 1)/dble(N)

      array(1) = int(dble(endI - startI + 1)/dble(N)/2.) + startI - 1
      if (array(1) == startI - 1) then
         array(1) = startI
      endif

      do i = 2, N
         array(i) = array(1) + int(dble(i - 1)*rtemp)
      enddo

   end subroutine linspaceI

   subroutine copysubmat_assumeshape(A, B, m, n, trans)
      implicit none
      integer m, n
      DT::A(m, n), B(m, n)
      real(kind=8)::n1, n2
      integer ii, jj, ijind
      character trans

      ! n1 = MPI_Wtime()
      if (trans == 'N') then
         do ii = 1, m
         do jj = 1, n
            B(ii , jj ) = A(ii , jj )
         end do
         end do
      elseif (trans == 'T') then
         do ii = 1, m
         do jj = 1, n
            B(jj , ii ) = A(ii , jj )
         end do
         end do
      endif
      ! n2 = MPI_Wtime()
      ! time_memcopy = time_memcopy + n2-n1

   end subroutine copysubmat_assumeshape

   subroutine copysubmat(A, ais, ajs, B, bis, bjs, m, n, trans)
      implicit none
      integer ais, ajs, bis, bjs, m, n
      DT::A(:, :), B(:, :)
      real(kind=8)::n1, n2
      integer ii, jj, ijind
      character trans

      ! n1 = MPI_Wtime()
      if (trans == 'N') then
         do ii = 1, m
         do jj = 1, n
            B(ii + bis - 1, jj + bjs - 1) = A(ii + ais - 1, jj + ajs - 1)
         end do
         end do
      elseif (trans == 'T') then
         do ii = 1, m
         do jj = 1, n
            B(jj + bis - 1, ii + bjs - 1) = A(ii + ais - 1, jj + ajs - 1)
         end do
         end do
      endif
      ! n2 = MPI_Wtime()
      ! time_memcopy = time_memcopy + n2-n1

   end subroutine copysubmat

   subroutine copymatT(A, B, m, n)
      implicit none
      integer m, n ! dimensions of A
      DT::A(:, :), B(:, :)
      real(kind=8)::n1, n2
      integer ii, jj, ijind

      ! n1 = MPI_Wtime()

      ! !$omp parallel do default(shared) private(ii,jj)
      do ii = 1, m
      do jj = 1, n
         B(jj, ii) = A(ii, jj)
      end do
      end do
      ! !$omp end parallel do

      ! n2 = MPI_Wtime()
      ! time_memcopy = time_memcopy + n2-n1

   end subroutine copymatT


   function myisnan(a)
      implicit none
      DTR a
      logical myisnan
#if __GNUC__ < 5
      myisnan=isnan(a)
#else
      myisnan=ieee_is_nan(a)
#endif
   end function myisnan

   function isnanMat(A, m, n)
      implicit none
      logical:: isnanMat
      DT::A(:, :)
      integer m, n, ii, jj
      isnanMat = .false.
      do ii = 1, m
      do jj = 1, n
         isnanMat = isnanMat .or. myisnan(abs(A(ii, jj)))
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

   subroutine LR_ReCompression(matU, matV, Singular, M, N, rmax, rank, SVD_tolerance, Flops)


      implicit none

      integer N, M, rmax
      real(kind=8), optional:: Flops
      real(kind=8):: flop
      real(kind=8) SVD_tolerance
      integer i, j, ii, jj, indx, rank_1, rank_2
      integer rank, ranknew, ldaU, ldaV

      DT::matU(:, :), matV(:, :)
      DTR::Singular(:)
      DT, allocatable :: QQ1(:, :), RR1(:, :), QQ2(:, :), RR2(:, :), UUsml(:, :), VVsml(:, :), tau_Q(:), mattemp(:, :), matU1(:, :), matV1(:, :)
      DTR, allocatable :: Singularsml(:)
      integer, allocatable :: jpvt1(:), jpvt2(:)

      if (present(Flops)) Flops = 0d0

      ldaU = size(matU, 1)
      ldaV = size(matV, 1)

      call assert(rmax <= min(M, N), 'rmax too big in LR_ReCompression')

      allocate (QQ1(M, rmax))
      ! call copymatN(matU(1:M,1:rmax),QQ1,M,rmax)
      QQ1 = matU(1:M, 1:rmax)
      allocate (tau_Q(rmax))
      ! allocate (jpvt1(rmax))
      ! jpvt1=0
      ! call geqp3f90(QQ1,jpvt1,tau_Q)
      call geqrff90(QQ1, tau_Q, flop=flop)
      if (present(Flops)) Flops = Flops + flop

      allocate (RR1(rmax, rmax))
      RR1 = 0d0
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, rmax
         do i = 1, j
            RR1(i, j) = QQ1(i, j)
         enddo
      enddo
      ! !$omp end parallel do
      call un_or_gqrf90(QQ1, tau_Q, M, rmax, rmax, flop=flop)
      deallocate (tau_Q)
      ! deallocate(jpvt1)
      if (present(Flops)) Flops = Flops + flop

      allocate (QQ2(N, rmax))
      call copymatT(matV(1:rmax, 1:N), QQ2, rmax, N)
      allocate (tau_Q(rmax))
      ! call geqrff90(QQ2,tau_Q)
      ! allocate (jpvt2(rmax))
      ! jpvt2=0
      call geqrff90(QQ2, tau_Q, flop=flop)

      if (present(Flops)) Flops = Flops + flop

      allocate (RR2(rmax, rmax))
      RR2 = 0d0
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, rmax
         do i = 1, j
            RR2(i, j) = QQ2(i, j)
         enddo
      enddo
      ! !$omp end parallel do

      call un_or_gqrf90(QQ2, tau_Q, N, rmax, rmax, flop=flop)
      ! deallocate(jpvt2)
      deallocate (tau_Q)
      if (present(Flops)) Flops = Flops + flop

      allocate (mattemp(rmax, rmax))
      mattemp = 0
      ! call zgemm('N','T',rmax,rmax,rmax, BPACK_cone, RR1, rmax,RR2,rmax,BPACK_czero,mattemp,rmax)
      call gemmf90(RR1, rmax, RR2, rmax, mattemp, rmax, 'N', 'T', rmax, rmax, rmax, BPACK_cone, BPACK_czero, flop=flop)
      if (present(Flops)) Flops = Flops + flop
      allocate (UUsml(rmax, rmax), VVsml(rmax, rmax), Singularsml(rmax))
      call SVD_Truncate(mattemp, rmax, rmax, rmax, UUsml, VVsml, Singularsml, SVD_tolerance, BPACK_SafeUnderflow, ranknew, flop=flop)
      if (present(Flops)) Flops = Flops + flop
      ! do i = 1, ranknew
      !    UUsml(1:rmax, i) = UUsml(1:rmax, i)*Singularsml(i)
      ! enddo
      ! call zgemm('N','N',M,ranknew,rmax, BPACK_cone, QQ1, M,UUsml,rmax,BPACK_czero,matU,ldaU)
      call gemmf90(QQ1, M, UUsml, rmax, matU, ldaU, 'N', 'N', M, ranknew, rmax, BPACK_cone, BPACK_czero, flop=flop)
      if (present(Flops)) Flops = Flops + flop
      ! call zgemm('N','T',ranknew,N,rmax, BPACK_cone, VVsml, rmax,QQ2,N,BPACK_czero,matV,ldaV)
      call gemmf90(VVsml, rmax, QQ2, N, matV, ldaV, 'N', 'T', ranknew, N, rmax, BPACK_cone, BPACK_czero, flop=flop)
      if (present(Flops)) Flops = Flops + flop

      rank = ranknew
      Singular(1:ranknew) = Singularsml(1:ranknew)

      deallocate (mattemp, RR1, QQ1, UUsml, VVsml, Singularsml)
      deallocate (QQ2, RR2)

   end subroutine LR_ReCompression


   subroutine LR_Add(chara,U1,V1,U2,V2,rank1,rank2,ranknew,matU,matV,M,N,SVD_tolerance,flops)
      implicit none
      character:: chara
      integer rank,ranknew,M,N,rank1,rank2,ranknew1,ranknew2
      DT::matU(M,rank1+rank2),matU1(M,rank1+rank2),matV(rank1+rank2,N),matV1(rank1+rank2,N)
      DT::U1(M,rank1),V1(N,rank1),U2(M,rank2),V2(N,rank2)
      integer ii,jj,i,j
      DT :: QQ1(M, rank1+rank2), RR1(min(M,rank1+rank2), rank1+rank2),QQ3(rank1+rank2,N), RR3(min(N,rank1+rank2),N), QQ2(N, rank1+rank2), RR2(min(N,rank1+rank2), rank1+rank2), UUsml(min(M,rank1+rank2),min(min(M,rank1+rank2), min(N,rank1+rank2))), VVsml(min(min(M,rank1+rank2), min(N,rank1+rank2)), min(N,rank1+rank2)), tau_Q(rank1+rank2), mattemp(min(M,rank1+rank2), min(N,rank1+rank2))
      real(kind=8) :: error, signs, Singularsml(min(min(M,rank1+rank2), min(N,rank1+rank2))), SVD_tolerance, flop
      real(kind=8),optional::flops


      if (chara == '+') signs=1d0
      if (chara == '-') signs=-1d0

      if (present(flops)) flops = 0
      rank=rank1+rank2
      do ii=1,M
          do jj=1,rank1
              matU(ii,jj) = U1(ii,jj)
          enddo
      enddo
      do ii=1,M
          do jj=1,rank2
              matU(ii,jj+rank1) = signs*U2(ii,jj)
          enddo
      enddo
      do ii=1,rank1
          do jj=1,N
              matV(ii,jj) = V1(jj,ii)
          enddo
      enddo
      do ii=1,rank2
          do jj=1,N
              matV(ii+rank1,jj) = V2(jj,ii)
          enddo
      enddo

      ranknew = rank1+rank2


!!!!!! ACA+SVD
#if 1
      call LR_Add_ACA(matU, matV, M, N, rank, ranknew, SVD_tolerance, SVD_tolerance,error,flops=flop)
      if (present(flops)) flops = flops + flop
#endif


!!!!!! QR
#if 0
      QQ1 = matU(1:M, 1:rank)
      call RRQR_LQ(QQ1, M, rank, min(M,rank), QQ1, RR1, SVD_tolerance, ranknew1, 'R', flops=flop)
      if (present(flops)) flops = flops + flop

      QQ3 = matV(1:rank, 1:N)
      call RRQR_LQ(QQ3, rank, N, min(N,rank), QQ3, RR3, SVD_tolerance, ranknew2, 'L', flops=flop)
      if (present(flops)) flops = flops + flop

      if(ranknew1<=ranknew2)then
          ranknew = ranknew1
          matU(1:M, 1:ranknew1)=QQ1(1:M, 1:ranknew1)
          matV1=matV
          call gemmf90(RR1, min(M,rank), matV1, rank, matV, rank, 'N', 'N', ranknew1, N, rank, BPACK_cone, BPACK_czero, flop=flop)
          if (present(flops)) flops = flops + flop
      else
          ranknew = ranknew2
          matV(1:ranknew2, 1:N)=RR3(1:ranknew2,1:N)
          matU1=matU
          call gemmf90(matU1, M, QQ3, rank, matU, M, 'N', 'N', M, ranknew2, rank, BPACK_cone, BPACK_czero, flop=flop)
          if (present(flops)) flops = flops + flop
      endif
#endif


!!!!!! QR+SVD
#if 0
      QQ1 = matU(1:M, 1:rank)
      call geqrff90(QQ1, tau_Q, flop=flop)
      if (present(flops)) flops = flops + flop

      RR1 = 0d0
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, rank
         do i = 1, min(j,M)
            RR1(i, j) = QQ1(i, j)
         enddo
      enddo
      ! !$omp end parallel do
      call un_or_gqrf90(QQ1, tau_Q, M, min(M,rank), min(M,rank), flop=flop)
      if (present(flops)) flops = flops + flop

      call copymatT(matV(1:rank, 1:N), QQ2, rank, N)
      call geqrff90(QQ2, tau_Q, flop=flop)
      if (present(flops)) flops = flops + flop

      RR2 = 0d0
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, rank
         do i = 1, min(N,j)
            RR2(i, j) = QQ2(i, j)
         enddo
      enddo
      ! !$omp end parallel do
      call un_or_gqrf90(QQ2, tau_Q, N, min(N,rank), min(N,rank), flop=flop)
      if (present(flops)) flops = flops + flop

      mattemp = 0
      call gemmf90(RR1, min(M,rank), RR2, min(N,rank), mattemp, min(M,rank), 'N', 'T', min(M,rank), min(N,rank), rank, BPACK_cone, BPACK_czero, flop=flop)
      if (present(flops)) flops = flops + flop
      call SVD_Truncate(mattemp, min(M,rank), min(N,rank), min(min(M,rank), min(N,rank)), UUsml, VVsml, Singularsml, SVD_tolerance, BPACK_SafeUnderflow, ranknew, flop=flop)
      if (present(flops)) flops = flops + flop
      call gemmf90(QQ1, M, UUsml, min(M,rank), matU, M, 'N', 'N', M, ranknew, min(M,rank), BPACK_cone, BPACK_czero, flop=flop)
      if (present(flops)) flops = flops + flop
      call gemmf90(VVsml, min(min(M,rank), min(N,rank)), QQ2, N, matV, rank, 'N', 'T', ranknew, N, min(N,rank), BPACK_cone, BPACK_czero, flop=flop)
      if (present(flops)) flops = flops + flop

      do jj=1,ranknew
          matU(:,jj) = matU(:,jj)*Singularsml(jj)
      enddo
#endif

  end subroutine LR_Add




  subroutine LR_Add_ACA(matU, matV, rankmax_r, rankmax_c, rmax, rank, tolerance, SVD_tolerance,error,flops)


      implicit none

      integer i, j, ii, jj, indx, rank_1, rank_2
      integer blocks, index_j, group_nn, rank_blocks
      integer group_m, group_n, size_of_groupm, size_of_groupn
      real(kind=8) norm_Z, norm_U, norm_V, tolerance, SVD_tolerance, dist
      integer edgefine_m, edgefine_n, level_butterfly, level_blocks
      integer edge_m, edge_n, header_m, header_n, Dimn, mn, Navr, itr
      integer rank, ranknew, row, column, rankmax, rankmax_c, rankmax_r, rankmax_min, rmax, idxs_r, idxs_c, frow
      DT value_Z, maxvalue
      DT inner_U, inner_V, ctemp, value_UVs
      real(kind=8) inner_UV, n1, n2, a, error, flop
      integer, allocatable:: select_column(:), select_row(:)
      DT::matU(rankmax_r, rmax), matV(rmax, rankmax_c),matU0(rankmax_r, rmax), matV0(rmax, rankmax_c)
      DT::matr(1, rankmax_c), matc(rankmax_r, 1)
      DTR::Singular(rmax)
      DT, allocatable:: row_R(:), column_R(:), value_UV(:)
      real(kind=8), allocatable:: norm_row_R(:), norm_column_R(:), norm_UVavrbynorm_Z(:)
      DT, allocatable :: QQ1(:, :), RR1(:, :), QQ2(:, :), RR2(:, :), UUsml(:, :), VVsml(:, :), tau_Q(:), mattemp(:, :), matU1(:, :), matV1(:, :)
      DTR, allocatable :: Singularsml(:)
      real(kind=8),optional::flops

      if (present(flops)) flops=0
      matU0 = matU
      matV0 = matV
      frow=1


      Navr = 1 !5 !10
      itr = 1
      allocate (norm_UVavrbynorm_Z(Navr))
      norm_UVavrbynorm_Z = 0

      n1 = MPI_Wtime()

      allocate (select_column(rankmax_c))
      allocate (select_row(rankmax_r))
      allocate (value_UV(max(rankmax_c, rankmax_r)))
      value_UV = 0

      rankmax_min = min(rmax,min(rankmax_r, rankmax_c))
      norm_Z = 0
      select_column = 0
      select_row = 0

      allocate (row_R(rankmax_c), column_R(rankmax_r))
      allocate (norm_row_R(rankmax_c), norm_column_R(rankmax_r))
      row_R = 0
      column_R = 0
      norm_row_R = 0
      norm_column_R = 0

      select_row(1) = frow

      call gemmf77('N', 'N', 1, rankmax_c, rmax, BPACK_cone, matU0(select_row(1), 1), rankmax_r, matV0, rmax, BPACK_czero, matr, 1)

      row_R = matr(1, :)
      norm_row_R = dble(row_R*conjg(cmplx(row_R, kind=8)))

      select_column(1) = maxloc(norm_row_R, 1)
      maxvalue = row_R(select_column(1))

      if (abs(maxvalue) < BPACK_SafeUnderflow) then

         do ii = 1, 100
            a = 0
            call random_number(a)
            select_row(1) = floor_safe(a*(rankmax_r - 1)) + 1

            call gemmf77('N', 'N', 1, rankmax_c, rmax, BPACK_cone, matU0(select_row(1), 1), rankmax_r, matV0, rmax, BPACK_czero, matr, 1)
            row_R = matr(1, :)
            norm_row_R = dble(row_R*conjg(cmplx(row_R, kind=8)))

            select_column(1) = maxloc(norm_row_R, 1)
            maxvalue = row_R(select_column(1))
            if (abs(maxvalue) > BPACK_SafeUnderflow) exit
         end do
         if (abs(maxvalue) < BPACK_SafeUnderflow) then
            rank = 1
            matU(:, 1) = 0
            matV(1, :) = 0
            Singular(1) = 0

            deallocate (select_column)
            deallocate (select_row)
            deallocate (value_UV)
            deallocate (row_R, column_R)
            deallocate (norm_row_R, norm_column_R)
            deallocate (norm_UVavrbynorm_Z)
            return
         endif
      end if

      ! !$omp parallel do default(shared) private(j)
      do j = 1, rankmax_c
         row_R(j) = row_R(j)/maxvalue
      enddo
      ! !$omp end parallel do
      ! !$omp parallel do default(shared) private(j)
      do j = 1, rankmax_c
         matV(1, j) = row_R(j)
      enddo
      ! !$omp end parallel do

      call gemmf77('N', 'N', rankmax_r, 1, rmax, BPACK_cone, matU0, rankmax_r, matV0(1,select_column(1)), rmax, BPACK_czero, matc, rankmax_r)
      column_R = matc(:, 1)
      norm_column_R = dble(column_R*conjg(cmplx(column_R, kind=8)))

      norm_column_R(select_row(1)) = 0

      ! !$omp parallel do default(shared) private(i)
      do i = 1, rankmax_r
         matU(i, 1) = column_R(i)
      enddo
      ! !$omp end parallel do

      norm_U = norm_vector(column_R, rankmax_r)
      norm_V = norm_vector(row_R, rankmax_c)
      norm_Z = norm_Z + norm_U*norm_V
      if (norm_Z > BPACK_SafeUnderflow) then
         norm_UVavrbynorm_Z(itr) = norm_U*norm_V/norm_Z
      else
         norm_UVavrbynorm_Z(itr) = 0
      endif
      itr = mod(itr, Navr) + 1

      ! if(rankmax<2)write(*,*)'rankmax'
      select_row(2) = maxloc(norm_column_R, 1)

      rank = 1
      ! write(*,*)column_R,row_R
      ! write(*,*)norm_Z*ACA_tolerance_forward**2,norm_U*norm_V,'hehe',ACA_tolerance_forward
      do while (tolerance**2 < (sum(norm_UVavrbynorm_Z)/min(Navr, rank)) .and. rank < rankmax_min)

         call gemmf77('N', 'N', 1, rankmax_c, rmax, BPACK_cone, matU0(select_row(rank + 1), 1), rankmax_r, matV0, rmax, BPACK_czero, matr, 1)
         row_R = matr(1, :)
         call gemmf77('N', 'N', 1, rankmax_c, rank, BPACK_cone, matU(select_row(rank + 1), 1), rankmax_r, matV, rmax, BPACK_czero, value_UV, 1)

         row_R = row_R - value_UV(1:rankmax_c)
         norm_row_R = dble(row_R*conjg(cmplx(row_R, kind=8)))
         if (present(flops)) flops = flops + rank*rankmax_c

         do i = 1, rank
            norm_row_R(select_column(i)) = 0
         enddo

         select_column(rank + 1) = maxloc(norm_row_R, 1)
         maxvalue = row_R(select_column(rank + 1))

         if (abs(maxvalue) < BPACK_SafeUnderflow) then
            ! write(*,*)'warning: zero pivot ',maxvalue,' in ACA, exiting with residual', sqrt((sum(norm_UVavrbynorm_Z)/Navr))
            exit
            row_R = 0
         else
            do j = 1, rankmax_c
               row_R(j) = row_R(j)/maxvalue
            enddo
         endif
         ! !$omp parallel do default(shared) private(j)

         ! !$omp end parallel do
         ! !$omp parallel do default(shared) private(j)
         do j = 1, rankmax_c
            matV(rank + 1, j) = row_R(j)
         enddo
         ! !$omp end parallel do


         call gemmf77('N', 'N', rankmax_r, 1, rmax, BPACK_cone, matU0, rankmax_r, matV0(1,select_column(rank+1)), rmax, BPACK_czero, matc, rankmax_r)
         column_R = matc(:, 1)
         call gemmf77('N', 'N', rankmax_r, 1, rank, BPACK_cone, matU, rankmax_r, matV(1, select_column(rank + 1)), rmax, BPACK_czero, value_UV, rankmax_r)
         column_R = column_R - value_UV(1:rankmax_r)
         norm_column_R = dble(column_R*conjg(cmplx(column_R, kind=8)))
         if (present(flops)) flops = flops + rank*rankmax_r

         do i = 1, rank + 1
            norm_column_R(select_row(i)) = 0
         enddo

         ! !$omp parallel do default(shared) private(i)
         do i = 1, rankmax_r
            matU(i, rank + 1) = column_R(i)
         enddo
         ! !$omp end parallel do

         norm_U = norm_vector(column_R, rankmax_r)
         norm_V = norm_vector(row_R, rankmax_c)

         inner_UV = 0
#ifdef HAVE_OPENMP
         !!$omp parallel do default(shared) private(j,i,ctemp,inner_V,inner_U) reduction(+:inner_UV)
#endif
         do j = 1, rank
            inner_U = 0
            inner_V = 0
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_U)
            do i = 1, rankmax_r
               ctemp = matU(i, rank + 1)*conjg(cmplx(matU(i, j), kind=8))
               inner_U = inner_U + ctemp
            enddo
            ! !$omp end parallel do
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_V)
            do i = 1, rankmax_c
               ctemp = matV(rank + 1, i)*conjg(cmplx(matV(j, i), kind=8))
               inner_V = inner_V + ctemp
            enddo
            ! !$omp end parallel do
            inner_UV = inner_UV + 2*dble(inner_U*inner_V)
         enddo
#ifdef HAVE_OPENMP
         !!$omp end parallel do
#endif
         norm_Z = norm_Z + inner_UV + norm_U*norm_V

         ! ! write(*,*)norm_Z,inner_UV,norm_U,norm_V,maxvalue,rank,'gan'
         ! if(myisnan(sqrt(norm_Z)))then
         ! write(*,*)inner_UV,norm_U,norm_V,maxvalue
         ! stop
         ! endif

         if (norm_Z > BPACK_SafeUnderflow) then
            norm_UVavrbynorm_Z(itr) = norm_U*norm_V/norm_Z
         else
            norm_UVavrbynorm_Z(itr) = 0
         endif
         itr = mod(itr, Navr) + 1

         if (present(flops)) flops = flops +  rank*rankmax_c + rank*rankmax_r

         rank = rank + 1
         if (rank > rmax) then
            ! write(*,*)'increase rmax',rank,rmax
            ! stop
            exit
         end if
         if (rank < rankmax_min) then
            select_row(rank + 1) = maxloc(norm_column_R, 1)
         endif

         if (norm_Z < 0) exit

         ! write(*,*)sqrt((sum(norm_UVavrbynorm_Z)/Navr)),sqrt(norm_U*norm_V),sqrt(norm_Z),rank,rankmax_min,tolerance**2<(sum(norm_UVavrbynorm_Z)/Navr)

      enddo
      ! write(*,*)sqrt((sum(norm_UVavrbynorm_Z)/Navr)),sqrt(norm_U*norm_V),sqrt(norm_Z),rank,rankmax_min,tolerance**2<(sum(norm_UVavrbynorm_Z)/Navr)

      error = sqrt((sum(norm_UVavrbynorm_Z)/Navr))

      ! write(*,*)select_row(1:rank),select_column(1:rank)

      deallocate (row_R, column_R)
      deallocate (norm_row_R, norm_column_R)
      deallocate (norm_UVavrbynorm_Z)

      n2 = MPI_Wtime()
      ! time_tmp = time_tmp + n2 - n1



!   ! ACA followed by SVD

!         allocate (QQ1(rankmax_r, rank))
!         QQ1 = matU(1:rankmax_r, 1:rank)
!         ! call copymatN(matU(1:rankmax_r,1:rank),QQ1,rankmax_r,rank)
!         allocate (tau_Q(rank))
!         call geqrff90(QQ1, tau_Q, flop=flop)
!         stats%Flop_Fill = stats%Flop_Fill + flop

!         allocate (RR1(rank, rank))
!         RR1 = 0d0
!         ! !$omp parallel do default(shared) private(i,j)
!         do j = 1, rank
!            do i = 1, j
!               RR1(i, j) = QQ1(i, j)
!            enddo
!         enddo
!         ! !$omp end parallel do
!         call un_or_gqrf90(QQ1, tau_Q, rankmax_r, rank, rank, flop=flop)
!         deallocate (tau_Q)
!         stats%Flop_Fill = stats%Flop_Fill + flop

!         allocate (QQ2(rankmax_c, rank))
!         call copymatT(matV(1:rank, 1:rankmax_c), QQ2, rank, rankmax_c)
!         allocate (tau_Q(rank))
!         call geqrff90(QQ2, tau_Q, flop=flop)
!         stats%Flop_Fill = stats%Flop_Fill + flop

!         allocate (RR2(rank, rank))
!         RR2 = 0d0
!         ! !$omp parallel do default(shared) private(i,j)
!         do j = 1, rank
!            do i = 1, j
!               RR2(i, j) = QQ2(i, j)
!            enddo
!         enddo
!         ! !$omp end parallel do
!         call un_or_gqrf90(QQ2, tau_Q, rankmax_c, rank, rank, flop=flop)
!         deallocate (tau_Q)
!         stats%Flop_Fill = stats%Flop_Fill + flop

!         allocate (mattemp(rank, rank))
!         mattemp = 0
!         call gemmf90(RR1, rank, RR2, rank, mattemp, rank, 'N', 'T', rank, rank, rank, BPACK_cone, BPACK_czero, flop=flop)
!         ! call zgemm('N','T',rank,rank,rank, BPACK_cone, RR1, rank,RR2,rank,BPACK_czero,mattemp,rank)
!         stats%Flop_Fill = stats%Flop_Fill + flop
!         allocate (UUsml(rank, rank), VVsml(rank, rank), Singularsml(rank))
!         call SVD_Truncate(mattemp, rank, rank, rank, UUsml, VVsml, Singularsml, SVD_tolerance, BPACK_SafeUnderflow, ranknew, flop=flop)
!         stats%Flop_Fill = stats%Flop_Fill + flop
!         ! call zgemm('N','N',rankmax_r,ranknew,rank, BPACK_cone, QQ1, rankmax_r,UUsml,rank,BPACK_czero,matU,rankmax_r)
!         call gemmf90(QQ1, rankmax_r, UUsml, rank, matU, rankmax_r, 'N', 'N', rankmax_r, ranknew, rank, BPACK_cone, BPACK_czero, flop=flop)
!         stats%Flop_Fill = stats%Flop_Fill + flop
!         ! call zgemm('N','T',ranknew,rankmax_c,rank, BPACK_cone, VVsml, rank,QQ2,rankmax_c,BPACK_czero,matV,rmax)
!         call gemmf90(VVsml, rank, QQ2, rankmax_c, matV, rmax, 'N', 'T', ranknew, rankmax_c, rank, BPACK_cone, BPACK_czero, flop=flop)
!         stats%Flop_Fill = stats%Flop_Fill + flop

!         rank = ranknew
!         Singular(1:ranknew) = Singularsml(1:ranknew)

!         do jj=1,ranknew
!             matU(:,jj) = matU(:,jj)*Singularsml(jj)
!         enddo

!         deallocate (mattemp, RR1, QQ1, UUsml, VVsml, Singularsml)
!         deallocate (QQ2, RR2)


      deallocate (select_column)
      deallocate (select_row)
      deallocate (value_UV)

      return

   end subroutine LR_Add_ACA



   subroutine LR_FnormUp(matU, matV, M, N, rskip, ruv, rup, ldV, normUV, normUVupdate, tolerance, Flops)


      implicit none

      integer N, M, rskip, ruv, rup, ldV
      real(kind=8) normUVupdate, normUV, inner_UV
      real(kind=8), optional:: Flops
      real(kind=8):: flop
      real(kind=8) tolerance
      integer i, j, k, ii, jj, indx, rank_1, rank_2
      integer rank, ranknew, mn, ranknew1, ranknew2
      DT::ctemp, inner_V, inner_U
      DT::matU(:, :), matV(:, :)
      DT, allocatable :: QQ1(:, :), RR1(:, :), QQ2(:, :), RR2(:, :), UUsml(:, :), VVsml(:, :), tau_Q(:), mattemp(:, :), matU1(:, :), matV1(:, :), UU(:, :), VV(:, :)
      DTR, allocatable :: Singularsml(:), Singular(:)
      integer, allocatable :: jpvt(:), jpvt1(:), jpvt2(:)

      if (present(Flops)) Flops = 0

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



      !!!  rskip means that the first rskip columns/rows of matU/matV are not accounted for in normUV
      inner_UV = 0
      if (ruv-rskip > 0) then
         allocate (UU(ruv-rskip, rup))
         UU = 0
         call gemmf77('T', 'N', ruv-rskip, rup, M, BPACK_cone, matU(1, rskip + 1), M, matU(1, ruv + 1), M, BPACK_czero, UU, ruv-rskip)
         allocate (VV(ruv-rskip, rup))
         VV = 0
         call gemmf77('N', 'T', ruv-rskip, rup, N, BPACK_cone, matV(rskip+1,1), ldV, matV(ruv + 1, 1), ldV, BPACK_czero, VV, ruv-rskip)
         do j = 1, ruv-rskip
            do k = 1, rup
               inner_UV = inner_UV + 2*dble(UU(j, k)*VV(j, k))
            enddo
         enddo
         deallocate (UU)
         deallocate (VV)
      endif

      normUV = sqrt(normUV**2d0 + normUVupdate**2d0 + inner_UV)

      if (present(Flops)) Flops = Flops + (M + N)*ruv*rup

   end subroutine LR_FnormUp

   subroutine LR_Fnorm(matU, matV, M, N, rmax, norm, tolerance, Flops)


      implicit none

      integer N, M, rmax
      real(kind=8) norm
      real(kind=8), optional:: Flops
      real(kind=8):: flop
      real(kind=8) tolerance
      integer i, j, ii, jj, indx, rank_1, rank_2
      integer rank, ranknew, mn, ranknew1, ranknew2

      DT::matU(:, :), matV(:, :)
      DT, allocatable :: QQ1(:, :), RR1(:, :), QQ2(:, :), RR2(:, :), UUsml(:, :), VVsml(:, :), tau_Q(:), mattemp(:, :), matU1(:, :), matV1(:, :), UU(:, :), VV(:, :)
      DTR, allocatable :: Singularsml(:), Singular(:)
      integer, allocatable :: jpvt(:), jpvt1(:), jpvt2(:)

      if (present(Flops)) Flops = 0

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
      ! call zgemm('N','N',mn,N,rmax, BPACK_cone, VV, mn,matV,size(matV,1),BPACK_czero,mattemp,mn)

      ! norm = fnorm(mattemp,mn,N)

      ! deallocate(QQ1)
      ! deallocate(UU)
      ! deallocate(VV)
      ! deallocate(Singular)
      ! deallocate(mattemp)

      allocate (QQ1(M, rmax))
      ! call copymatN(matU(1:M,1:rmax),QQ1,M,rmax)
      QQ1 = matU(1:M, 1:rmax)
      allocate (tau_Q(rmax))
      call geqrff90(QQ1, tau_Q, flop=flop)
      if (present(Flops)) Flops = Flops + flop
      deallocate (tau_Q)

      allocate (RR1(rmax, rmax))
      RR1 = 0d0
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, rmax
         do i = 1, j
            RR1(i, j) = QQ1(i, j)
         enddo
      enddo
      ! !$omp end parallel do

      allocate (QQ2(N, rmax))
      call copymatT(matV(1:rmax, 1:N), QQ2, rmax, N)
      allocate (tau_Q(rmax))
      call geqrff90(QQ2, tau_Q, flop=flop)
      if (present(Flops)) Flops = Flops + flop
      deallocate (tau_Q)

      allocate (RR2(rmax, rmax))
      RR2 = 0d0
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, rmax
         do i = 1, j
            RR2(i, j) = QQ2(i, j)
         enddo
      enddo
      ! !$omp end parallel do

      allocate (mattemp(rmax, rmax))
      mattemp = 0
      ! call zgemm('N','T',rmax,rmax,rmax, BPACK_cone, RR1, rmax,,rmax,BPACK_czero,mattemp,rmax)
      call gemmf90(RR1, rmax, RR2, rmax, mattemp, rmax, 'N', 'T', rmax, rmax, rmax, BPACK_cone, BPACK_czero, flop=flop)
      if (present(Flops)) Flops = Flops + flop
      norm = fnorm(mattemp, rmax, rmax)

      deallocate (mattemp, RR1, QQ1)
      deallocate (QQ2, RR2)

      ! allocate(QQ1(M,rmax))
      ! QQ1 = matU(1:M,1:rmax)
      ! allocate (tau_Q(rmax))
      ! allocate (jpvt1(rmax))
      ! jpvt1=0
      ! call geqp3modf90(QQ1,jpvt1,tau_Q,tolerance,BPACK_SafeUnderflow,ranknew1,flop=flop)
      ! if(present(Flops))Flops= Flops + flop
      ! deallocate(tau_Q)

      ! allocate(QQ2(N,rmax))
      ! call copymatT(matV(1:rmax,1:N),QQ2,rmax,N)
      ! allocate (tau_Q(rmax))
      ! allocate (jpvt2(rmax))
      ! jpvt2=0
      ! call geqp3modf90(QQ2,jpvt2,tau_Q,tolerance,BPACK_SafeUnderflow,ranknew2,flop=flop)
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
      ! ! call zgemm('N','T',ranknew1,ranknew2,rmax, BPACK_cone, RR1, ranknew1,RR2,ranknew2,BPACK_czero,mattemp,ranknew1)
      ! call gemmf90(RR1, ranknew1,RR2,ranknew2,mattemp,ranknew1,'N','T',ranknew1,ranknew2,rmax, BPACK_cone,BPACK_czero,flop=flop)

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
      ! call geqp3modf90(QQ1,jpvt,tau_Q,tolerance,BPACK_SafeUnderflow,ranknew,flop=flop)
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
      ! ! call zgemm('N','N',ranknew,N,rmax, BPACK_cone, RR1, ranknew,QQ2,rmax,BPACK_czero,mattemp,ranknew)
      ! call gemmf90(RR1, ranknew,QQ2,rmax,mattemp,ranknew,'N','N',ranknew,N,rmax, BPACK_cone,BPACK_czero,flop=flop)

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

   subroutine GetRank(M, N, mat, rank, eps, flop)
!
!
      implicit none
      integer M, N, rank, mn, i, mnl, ii, jj
      DT::mat(:, :)
      real(kind=8) eps
      DT, allocatable :: UU(:, :), VV(:, :), Atmp(:, :), A_tmp(:, :), tau(:), mat1(:, :)
      DTR, allocatable :: Singular(:)
      integer, allocatable :: jpvt(:)
      real(kind=8), optional::flop

      allocate (mat1(M, N))

      mat1 = mat
      mn = min(M, N)
      mnl = max(M, N)
      allocate (UU(M, mn))
      allocate (VV(mn, N))
      allocate (Singular(mn))

      if (myisnan(fnorm(mat, M, N))) then
         write (*, *) 'input matrix NAN in GetRank'
         stop
      end if

      call gesvd_robust(mat1, Singular, UU, VV, M, N, mn, flop)
! write(*,*)Singular,'hh'
      if (Singular(1) < BPACK_SafeUnderflow) then
         rank = 1
         deallocate (UU, VV, Singular)
      else
         if (myisnan(sum(Singular))) then
            deallocate (UU, VV, Singular)
            write (*, *) 'gesvd_robust wrong in GetRank, switching to QR'

            allocate (Atmp(mnl, mn))
            allocate (A_tmp(mn, mn))

            if (M >= N) then
               Atmp = mat
            else
               call copymatT(mat, Atmp, M, N)
            end if

            allocate (jpvt(mn))
            allocate (tau(mn))

            ! RRQR
            jpvt = 0
            call geqp3f90(Atmp, jpvt, tau, flop)
            if (myisnan(fnorm(Atmp, mnl, mn))) then
               write (*, *) 'Q or R has NAN in GetRank'
               stop
            end if
            A_tmp = 0
#ifdef HAVE_OPENMP
            !!$omp parallel do default(shared) private(ii,jj)
#endif
            do ii = 1, mn
               do jj = ii, mn
                  A_tmp(ii, jj) = Atmp(ii, jj)
               enddo
            enddo
#ifdef HAVE_OPENMP
            !!$omp end parallel do
#endif
            rank = mn
            do i = 1, mn
               if (abs(A_tmp(i, i))/abs(A_tmp(1, 1))/mnl <= eps) then
                  rank = i
                  if (abs(A_tmp(i, i)) < BPACK_SafeUnderflow) rank = i - 1
                  exit
               end if
            end do

            deallocate (jpvt)
            deallocate (tau)
            deallocate (Atmp)
            deallocate (A_tmp)

         else

            ! write(*,*)Singular,'hh'

            rank = mn
            do i = 1, mn
               if (Singular(i)/Singular(1) <= eps .or. Singular(i) < 1d-60) then
                  rank = i
                  if (Singular(i) < Singular(1)*eps/10 .or. Singular(i) < 1d-60) rank = i - 1
                  exit
               end if
            end do

            deallocate (UU, VV, Singular)

         end if
      endif
      deallocate (mat1)

   end subroutine GetRank

   subroutine PComputeRange(M_p, N, mat, rank, eps, ptree, pgno, Flops, norm_tol)

      implicit none
      integer M, N, M_loc, rank, mn, i, mnl, ii, jj, rrflag
      DT::mat(:, :)
      real(kind=8) eps
      logical::small
      DT, allocatable :: UU(:, :), VV(:, :), Atmp(:, :), A_tmp(:, :), tau(:), mat1D(:, :), mat2D(:, :)
      DTR, allocatable :: Singular(:)
      real(kind=8) :: flop,norm
      integer, allocatable :: jpvt(:)
      integer, allocatable :: ipiv(:), jpiv(:), JPERM(:)
      integer:: M_p(:, :)
      integer, allocatable:: M_p_1D(:, :)
      type(proctree)::ptree
      integer pgno, proc, nproc, nb1Dc, nb1Dr, ctxt1D, ctxt, idxs_o, idxe_o, ierr
      integer myArows, myAcols, info, nprow, npcol, myrow, mycol, taun
      integer::descsMat1D(9), descsMat2D(9)
      real(kind=8), optional:: Flops
      real(kind=8), optional:: norm_tol

      rank = 0

      if (present(Flops)) Flops = 0d0

      if (IOwnPgrp(ptree, pgno)) then

         proc = ptree%MyID - ptree%pgrp(pgno)%head
         nproc = ptree%pgrp(pgno)%nproc
         M_loc = M_p(proc + 1, 2) - M_p(proc + 1, 1) + 1
         M = M_p(nproc, 2)

         if (nproc == 1) then
            call ComputeRange(M, N, mat, rank, 1, eps, Flops=flop,norm_tol=norm_tol)
            if (present(Flops)) Flops = Flops + flop
         else

            mn = min(M, N)


            small=.False.
            if(present(norm_tol))then
               norm = fnorm(mat,M_loc,N)**2d0
               call MPI_ALLREDUCE(MPI_IN_PLACE, norm, 1, MPI_double_precision, MPI_SUM, ptree%pgrp(pgno)%Comm, ierr)
               norm = sqrt(norm)
               if(norm/sqrt(dble(N))<norm_tol)small=.True.
            endif

            if(small .eqv. .True.)then
               rank = 1
               mat(:, 1) = 0d0
            else
               !!!!>**** generate 2D grid blacs quantities
               ctxt = ptree%pgrp(pgno)%ctxt
               call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
               if (myrow /= -1 .and. mycol /= -1) then
                  myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
                  myAcols = numroc_wp(N, nbslpk, mycol, 0, npcol)
                  call descinit_wp(descsMat2D, M, N, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
                  allocate (mat2D(max(1,myArows), max(1,myAcols)))
                  mat2D = 0
               else
                  descsMat2D(2) = -1
                  allocate (mat2D(1, 1))
                  mat2D = 0
               endif

               !!!!>**** redistribution of input matrix
               call Redistribute1Dto2D(mat, M_p, 0, pgno, mat2D, M, 0, pgno, N, ptree)

               ! Compute RRQR
               if (myrow /= -1 .and. mycol /= -1) then
                  allocate (ipiv(myAcols))
                  ipiv = 0
                  taun = numroc_wp(mn, nbslpk, mycol, 0, npcol)
                  allocate (tau(taun))
                  tau = 0
                  allocate (jpiv(N))
                  jpiv = 0
                  allocate (JPERM(N))
                  JPERM = 0
                  call pgeqpfmodf90(M, N, mat2D, 1, 1, descsMat2D, ipiv, tau, JPERM, jpiv, rank, eps, BPACK_SafeUnderflow, flop=flop)
                  if (present(Flops)) Flops = Flops + flop/dble(nprow*npcol)

                  if (rank > 0) then
                     call pun_or_gqrf90(ctxt, mat2D, tau, M, rank, rank, descsMat2D, 1, 1, flop=flop)
                     if (present(Flops)) Flops = Flops + flop/dble(nprow*npcol)
                  else
                     rank = 1
                     if (myArows > 0 .and. myAcols > 0) mat2D = 0d0
                  endif

                  deallocate (ipiv)
                  deallocate (tau)
                  deallocate (jpiv)
                  deallocate (JPERM)

                  if (myisnan(fnorm(mat2D, max(1,myArows), max(1,myAcols)))) then
                     write (*, *) 'Q or R has NAN in PComputeRange'
                     stop
                  end if
               endif

               call MPI_ALLREDUCE(MPI_IN_PLACE, rank, 1, MPI_integer, MPI_MAX, ptree%pgrp(pgno)%Comm, ierr)

               !!!!>**** redistribution of output matrix
               call Redistribute2Dto1D(mat2D, M, 0, pgno, mat, M_p, 0, pgno, N, ptree)

               deallocate (mat2D)

            endif
         endif
      endif

   end subroutine PComputeRange

   subroutine ComputeRange(M, N, mat, rank, rrflag, eps, Flops,norm_tol)

      implicit none
      logical::small
      integer M, N, rank, mn, i, mnl, ii, jj, rrflag
      DT::mat(:, :)
      real(kind=8) eps
      DT, allocatable :: UU(:, :), VV(:, :), Atmp(:, :), A_tmp(:, :), tau(:)
      DTR, allocatable :: Singular(:)
      integer, allocatable :: jpvt(:)
      real(kind=8), optional:: Flops
      real(kind=8), optional:: norm_tol
      real(kind=8):: flop
      DTR:: matnorm
      if (present(Flops)) Flops = 0d0

      mn = min(M, N)
      matnorm=fnorm(mat, M, N)
      if (myisnan(matnorm)) then
         write (*, *) 'input matrix NAN in ComputeRange'
         stop
      end if

      allocate (jpvt(N))
      jpvt = 0
      allocate (tau(mn))
      if (rrflag == 1) then
         ! RRQR
         small=.False.
         if(present(norm_tol))then
            if(matnorm/sqrt(dble(N))<norm_tol)small=.True.
         endif

         if(small .eqv. .True.)then
            rank = 1
            mat(1:M, 1) = 0d0
         else
            call geqp3modf90(mat, jpvt, tau, eps, BPACK_SafeUnderflow, rank, flop=flop)
            if (present(Flops)) Flops = Flops + flop

            if (rank > 0) then
               call un_or_gqrf90(mat, tau, M, rank, rank, flop=flop)
               if (present(Flops)) Flops = Flops + flop
            else
               rank = 1
               mat(1:M, 1) = 0d0
            endif
         endif

      else
         call geqp3f90(mat, jpvt, tau, flop=flop)
         if (present(Flops)) Flops = Flops + flop
         call un_or_gqrf90(mat, tau, M, N, mn, flop=flop)
         if (present(Flops)) Flops = Flops + flop
         rank = mn
      endif
      if (myisnan(fnorm(mat, M, N))) then
         write (*, *) 'Q or R has NAN in ComputeRange'
         stop
      end if

      deallocate (jpvt)
      deallocate (tau)

   end subroutine ComputeRange

   subroutine CheckRandomizedLR(M, N, mat, tolerance)

      implicit none
      integer M, N, rank, mn, i, j, k, rankmax_c, rankmax_r, rankmax_min, flag0, rank_new, rmax
      DT::mat(M, N), ctemp
      DT, allocatable::mat1(:, :), mat2(:, :)
      real(kind=8) tolerance, Smax
      DT, allocatable :: UU(:, :), VV(:, :), matrix_U(:, :), matrix_V(:, :), U_new(:, :), V_new(:, :), U_new1(:, :), V_new1(:, :), test_in(:, :), test_out1(:, :), test_out2(:, :), test_out3(:, :)
      DTR, allocatable :: Singular(:)
      integer, allocatable:: select_row(:), select_column(:)
      DT, allocatable::MatrixSubselection(:, :)

      allocate (mat1(M, N))
      allocate (mat2(M, N))
      call GetRank(M, N, mat, rank, tolerance)
      rankmax_r = rank*3
      write (*, *) rankmax_r, min(m, n)
      if (rankmax_r > min(m, n)) rankmax_r = min(m, n)
      rankmax_c = rankmax_r
      rankmax_min = min(rankmax_r, rankmax_c)

! rankmax_r = M
! rankmax_c = N
! rankmax_min = min(rankmax_r,rankmax_c)

! method 1: SeudoSkeleton

      allocate (select_row(rankmax_r), select_column(rankmax_c))
      call linspaceI(1, M, rankmax_r, select_row)
      call linspaceI(1, N, rankmax_c, select_column)

      allocate (MatrixSubselection(rankmax_r, rankmax_c))
      MatrixSubselection = mat(select_row, select_column)

      allocate (matrix_U(M, rankmax_c), matrix_V(N, rankmax_r))
      matrix_U = mat(1:M, select_column)
      do i = 1, N
      do j = 1, rankmax_r
         matrix_V(i, j) = mat(select_row(j), i)
      end do
      end do

      deallocate (select_column, select_row)

      ! write(*,*)tolerance
      allocate (UU(rankmax_r, rankmax_min), VV(rankmax_min, rankmax_c), Singular(rankmax_min))
      call gesvd_robust(MatrixSubselection, Singular, UU, VV, rankmax_r, rankmax_c, rankmax_min)
      deallocate (MatrixSubselection)

      flag0 = 0; i = 0
      do while (flag0 == 0 .and. i < rankmax_min)
         i = i + 1
         if (Singular(i)/Singular(1) <= tolerance/10) then
            flag0 = 1
         endif
      enddo

      rank_new = i
! if(rank_new>rank)rank_new=rank

      write (*, *) 'rank=', rank, 'rank_new', rank_new

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

      allocate (U_new(M, rank_new))
      allocate (V_new(rank_new, N))

      do j = 1, rank_new
         do i = 1, M
            ctemp = 0
            do k = 1, rankmax_c
               ctemp = ctemp + matrix_U(i, k)*conjg(cmplx(VV(j, k), kind=8))
            enddo
            U_new(i, j) = ctemp
         enddo
      enddo

      do j = 1, rank_new
         do i = 1, N
            ctemp = 0
            do k = 1, rankmax_r
               ctemp = ctemp + conjg(cmplx(UU(k, j), kind=8))*matrix_V(i, k)
            enddo
            V_new(j, i) = ctemp/Singular(j)
         enddo
      enddo

      deallocate (matrix_U, VV)
      deallocate (matrix_V, UU, Singular)

! call gemm_omp(U_new,V_new,mat1,M,N,rank_new)
      call gemmf90(U_new, M, V_new, rank_new, mat1, M, 'N', 'N', M, N, rank_new, BPACK_cone, BPACK_czero)

      write (*, *) M, N
      do j = 1, N
      do i = 1, M
         write (211, *) dble(mat1(i, j))
         write (212, *) aimag(cmplx(mat1(i, j), kind=8))
      end do
      end do

! write(*,*)'F-norm residual:', fnorm(mat-mat1,M,N)/fnorm(mat,M,N),fnorm(mat-mat2,M,N)/fnorm(mat,M,N)
      write (*, *) 'F-norm residual:', fnorm(mat - mat1, M, N)/fnorm(mat, M, N)
      allocate (test_in(N, 1))
      allocate (test_out1(M, 1))
      allocate (test_out2(M, 1))
! allocate(test_out3(M,1))
      do i = 1, N
         call random_dp_number(test_in(i, 1))
      end do

! call gemm_omp(mat,test_in,test_out1,M,1,N)
      call gemmf90(mat, M, test_in, N, test_out1, M, 'N', 'N', M, 1, N, BPACK_cone, BPACK_czero)
! call gemm_omp(mat1,test_in,test_out2,M,1,N)
      call gemmf90(mat1, M, test_in, N, test_out2, M, 'N', 'N', M, 1, N, BPACK_cone, BPACK_czero)

! call gemm_omp(mat2,test_in,test_out3,M,1,N)
! write(*,*)'testing vector error:', fnorm(test_out1-test_out2,M,1)/fnorm(test_out1,M,1),fnorm(test_out1-test_out3,M,1)/fnorm(test_out1,M,1)
      write (*, *) 'testing vector error:', fnorm(test_out1 - test_out2, M, 1)/fnorm(test_out1, M, 1)

! method 2: ACA

      rmax = min(M, N)
      allocate (matrix_U(M, rmax))
      allocate (matrix_V(rmax, N))

      call ACA_CompressionFull(mat, matrix_U, matrix_V, M, N, rmax, rank_new, tolerance*0.3d0, tolerance)

! call gemm_omp(matrix_U(1:M,1:rank_new),matrix_V(1:rank_new,1:N),mat2,M,N,rank_new)
      call gemmf90(matrix_U, M, matrix_V, rmax, mat2, M, 'N', 'N', M, N, rank_new, BPACK_cone, BPACK_czero)

      write (*, *) 'F-norm residual:', fnorm(mat - mat2, M, N)/fnorm(mat, M, N), ' rank:', rank_new

      deallocate (mat1)
      deallocate (mat2)

   end subroutine CheckRandomizedLR

! generate random permutation of 1:N
   subroutine rperm(N, p)

      integer, intent(in):: N
      integer:: p(N)

      integer:: i
      integer:: k, j, ipj, itemp, m
      real(kind=8), dimension(100) :: u

      call assert(N < 1d9, 'In rperm, N too large')
! write(*,*)p(1)
      p = (/(i, i=1, N)/)

! Generate up to 100 U(0,1) numbers at a time.
      do i = 1, N, 100
         m = min(N - i + 1, 100)
         call random_number(u)
         do j = 1, m
            ipj = i + j - 1
            k = int(u(j)*(N - ipj + 1)) + ipj
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

      call random_seed(size=n)
      allocate (seed(n))
! First try if the OS provides a random number generator

! Fallback to XOR:ing the current time and pid. The PID is
! useful in case one launches multiple instances of the same
! program in parallel.
      call system_clock(t)
      if (t == 0) then
         call date_and_time(values=dt)
         t = (dt(1) - 1970)*365_int64*24*60*60*1000 &
             + dt(2)*31_int64*24*60*60*1000 &
             + dt(3)*24_int64*60*60*1000 &
             + dt(5)*60*60*1000 &
             + dt(6)*60*1000 + dt(7)*1000 &
             + dt(8)
      end if
#ifdef CRAY
      pid = 0
#else
      pid = getpid()
#endif
      t = ieor(t, int(pid, kind(t)))
      do i = 1, n
         seed(i) = lcg(t)
      end do
      call random_seed(put=seed)
      deallocate(seed)
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
         s = mod(s*279470273_int64, 4294967291_int64)
         lcg = int(mod(s, int(huge(0), int64)), kind(0))
      end function lcg
   end subroutine init_random_seed

   subroutine random_dp_number(val)
      implicit none
      real(kind=8):: a = 0, b = 0, c = 0, d = 0
      integer seed
      DT val

#if DAT==0
         ! Uniform distribution
         call random_number(a)
         call random_number(b)
         call random_number(c)
         call random_number(d)
         if (c < 0.5d0) then
            a = -a
         endif
         if (d < 0.5d0) then
            b = -b
         endif
         val = a + BPACK_junit*b

         ! ! Normal distribution
         ! call random_number(a)
         ! seed = a*10000000
         ! val =  c8_normal_01 ( seed )
#elif DAT==1
         ! Uniform distribution
         call random_number(a)
         val = a*2d0 - 1d0

         ! ! Normal distribution
         ! call random_number(a)
         ! seed = a*10000000
         ! val =  dble(c8_normal_01 ( seed ))

#elif DAT==2
         ! Uniform distribution
         call random_number(a)
         call random_number(b)
         call random_number(c)
         call random_number(d)
         if (c < 0.5d0) then
            a = -a
         endif
         if (d < 0.5d0) then
            b = -b
         endif
         val = a + BPACK_junit*b

         ! ! Normal distribution
         ! call random_number(a)
         ! seed = a*10000000
         ! val =  c8_normal_01 ( seed )
#elif DAT==3
         ! Uniform distribution
         call random_number(a)
         val = a*2d0 - 1d0

         ! ! Normal distribution
         ! call random_number(a)
         ! seed = a*10000000
         ! val =  dble(c8_normal_01 ( seed ))
#endif

      return
   end subroutine random_dp_number

   complex(kind=8) function c8_normal_01(seed)

!>*****************************************************************************80
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
      integer(kind=4) seed
      real(kind=8) v1
      real(kind=8) v2
      real(kind=8) x_c
      real(kind=8) x_r

      v1 = r8_uniform_01(seed)
      v2 = r8_uniform_01(seed)

      x_r = sqrt(-2.0D+00*log(v1))*cos(2.0D+00*r8_pi*v2)
      x_c = sqrt(-2.0D+00*log(v1))*sin(2.0D+00*r8_pi*v2)

      c8_normal_01 = cmplx(x_r, x_c, kind=8)

      return
   end function c8_normal_01

   real(kind=8) function r8_normal_01(seed)
!>*****************************************************************************80
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
      integer(kind=4) seed
      real(kind=8) x

      r1 = r8_uniform_01(seed)
      r2 = r8_uniform_01(seed)
      x = sqrt(-2.0D+00*log(r1))*cos(2.0D+00*r8_pi*r2)

      r8_normal_01 = x

      return
   end function r8_normal_01

   real(kind=8) function r8_uniform_01(seed)

!>*****************************************************************************80
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

      integer(kind=4) k
      integer(kind=4) seed

      k = seed/127773

      seed = 16807*(seed - k*127773) - k*2836

      if (seed < 0) then
         seed = seed + 2147483647
      end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
      r8_uniform_01 = real(seed, kind=8)*4.656612875D-10

      return
   end function r8_uniform_01

! ****************************************************************************************** !
! ***************                          assert function                               !>*********** !
! ***************                                                                                                      !>*********** !
! ****************************************************************************************** !
! output an error message is the first argument is false
   subroutine assert(statement, msg)
      logical::statement
      character(*)::msg
#ifndef NDEBUG
      if (.not. statement) then
         write (*, *) msg
         stop
      end if
#endif
   end subroutine assert

   function floor_safe(input)
      real(kind=8)::input
      integer floor_safe
      integer input_nint
      input_nint = NINT(input)
      if (abs(input_nint - input) < BPACK_SafeEps) then
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
      if (abs(input_nint - input) < BPACK_SafeEps) then
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
      if (abs(input_nint - input) < BPACK_SafeEps) then
         INT_safe = input_nint
      else
         INT_safe = INT(input)
      end if
   end function INT_safe

   subroutine cscalar(a, b, c)
      implicit none
      real(kind=8) b(3)
      complex(kind=8) a(3), c
      complex(kind=8) ax, ay, az
      real(kind=8) bx, by, bz
      ax = a(1); ay = a(2); az = a(3)
      bx = b(1); by = b(2); bz = b(3)
      c = ax*bx + ay*by + az*bz
      return
   end subroutine cscalar

   subroutine scalar(a, b, c)
      implicit none
      real(kind=8) a(3), b(3), c
      real(kind=8) ax, ay, az
      real(kind=8) bx, by, bz
      ax = a(1); ay = a(2); az = a(3)
      bx = b(1); by = b(2); bz = b(3)
      c = ax*bx + ay*by + az*bz
      return
   end subroutine scalar

   subroutine rrcurl(a, b, c)
      implicit none
      real(kind=8) a(3), b(3), c(3)
      real(kind=8) ax, ay, az
      real(kind=8) bx, by, bz
      real(kind=8) cx, cy, cz
      ax = a(1); ay = a(2); az = a(3)
      bx = b(1); by = b(2); bz = b(3)

      cx = ay*bz - by*az
      cy = -ax*bz + bx*az
      cz = ax*by - bx*ay
      c(1) = cx; c(2) = cy; c(3) = cz
      return
   end subroutine rrcurl

   subroutine rccurl(a, b, c)
      implicit none
      real(kind=8) a(3)
      complex(kind=8) b(3), c(3)
      real(kind=8) ax, ay, az
      complex(kind=8) bx, by, bz
      complex(kind=8) cx, cy, cz
      ax = a(1); ay = a(2); az = a(3)
      bx = b(1); by = b(2); bz = b(3)

      cx = ay*bz - by*az
      cy = -ax*bz + bx*az
      cz = ax*by - bx*ay
      c(1) = cx; c(2) = cy; c(3) = cz
      return
   end subroutine rccurl

   subroutine cccurl(a, b, c)
      implicit none
      complex(kind=8) a(3), b(3), c(3)
      complex(kind=8) ax, ay, az
      complex(kind=8) bx, by, bz
      complex(kind=8) cx, cy, cz
      ax = a(1); ay = a(2); az = a(3)
      bx = b(1); by = b(2); bz = b(3)

      cx = ay*bz - by*az
      cy = -ax*bz + bx*az
      cz = ax*by - bx*ay
      c(1) = cx; c(2) = cy; c(3) = cz
      return
   end subroutine cccurl


   subroutine cscalar_sp(a, b, c)
      implicit none
      real(kind=4) b(3)
      complex(kind=4) a(3), c
      complex(kind=4) ax, ay, az
      real(kind=4) bx, by, bz
      ax = a(1); ay = a(2); az = a(3)
      bx = b(1); by = b(2); bz = b(3)
      c = ax*bx + ay*by + az*bz
      return
    end subroutine cscalar_sp

    subroutine scalar_sp(a, b, c)
      implicit none
      real(kind=4) a(3), b(3), c
      real(kind=4) ax, ay, az
      real(kind=4) bx, by, bz
      ax = a(1); ay = a(2); az = a(3)
      bx = b(1); by = b(2); bz = b(3)
      c = ax*bx + ay*by + az*bz
      return
    end subroutine scalar_sp

    subroutine rrcurl_sp(a, b, c)
      implicit none
      real(kind=4) a(3), b(3), c(3)
      real(kind=4) ax, ay, az
      real(kind=4) bx, by, bz
      real(kind=4) cx, cy, cz
      ax = a(1); ay = a(2); az = a(3)
      bx = b(1); by = b(2); bz = b(3)

      cx = ay*bz - by*az
      cy = -ax*bz + bx*az
      cz = ax*by - bx*ay
      c(1) = cx; c(2) = cy; c(3) = cz
      return
    end subroutine rrcurl_sp

    subroutine rccurl_sp(a, b, c)
      implicit none
      real(kind=4) a(3)
      complex(kind=4) b(3), c(3)
      real(kind=4) ax, ay, az
      complex(kind=4) bx, by, bz
      complex(kind=4) cx, cy, cz
      ax = a(1); ay = a(2); az = a(3)
      bx = b(1); by = b(2); bz = b(3)

      cx = ay*bz - by*az
      cy = -ax*bz + bx*az
      cz = ax*by - bx*ay
      c(1) = cx; c(2) = cy; c(3) = cz
      return
    end subroutine rccurl_sp

    subroutine cccurl_sp(a, b, c)
      implicit none
      complex(kind=4) a(3), b(3), c(3)
      complex(kind=4) ax, ay, az
      complex(kind=4) bx, by, bz
      complex(kind=4) cx, cy, cz
      ax = a(1); ay = a(2); az = a(3)
      bx = b(1); by = b(2); bz = b(3)

      cx = ay*bz - by*az
      cy = -ax*bz + bx*az
      cz = ax*by - bx*ay
      c(1) = cx; c(2) = cy; c(3) = cz
      return
   end subroutine cccurl_sp



   real(kind=8) function norm_vector(vector, n)
      implicit none
      integer n
      DT:: vector(n)
      integer i
      real(kind=8) sum

      sum = 0.0d0
      do i = 1, n
         sum = sum + dble(vector(i)*conjg(cmplx(vector(i), kind=8)))
      enddo

      norm_vector = sum

      return
   end function norm_vector

   subroutine LinearSolve(m, n, k, A, b, x, eps_r, verbose, Flops)
      !
      !
      implicit none

      integer m, n, k, mn_min
      DT:: A(m, n), x(n, k), b(m, k)
      DTR, allocatable::Singular(:)
      DT, allocatable:: Atmp(:, :), xtmp(:, :), tau(:), A_tmp(:, :), UU(:, :), VV(:, :), UU_h(:, :), VV_h(:, :), VV_h2(:, :), z_arb(:, :), x_arb(:, :), matrixtemp(:, :), A_tmp_rank(:, :), xtmp_rank(:, :), xtmp_rank3(:, :), A_tmp_rank2(:, :), xtmp_rank2(:, :)
      real(kind=8):: eps_r
      integer ii, jj, i, j, flag0, rank
      integer, allocatable:: jpvt(:)
      real(kind=8), optional::Flops
      real(kind=8)::flop
      integer::verbose

      DT:: alpha, beta
      alpha = 1d0
      beta = 0d0

      if (present(Flops)) Flops = 0

      allocate (xtmp(n, k))
      allocate (tau(n))
      allocate (jpvt(n))
      allocate (A_tmp(n, n))
      allocate (Atmp(m, n))
      allocate (UU(m, n))
      allocate (VV(n, n))
      allocate (Singular(n))

      ! if (m < n) write (*, *) m, n
      ! call assert(m >= n, 'm should not be less than n for least square')

      if (fnorm(b, m, k) < BPACK_SafeUnderflow) then
!                write(*,*)'warning: RHS zero in least square. |b|= ',fnorm(b,m,k),'size b: ',m,k,'size A',m,n
         x = 0
      else

         Atmp = A

         ! SVD
         call gesvd_robust(Atmp, Singular, UU, VV, m, n, n, flop)
         if (present(Flops)) Flops = Flops + flop

! !!!!!!!  If SVD fails, uncomment the following If statement, but the code might become slow
         ! if(myisnan(sum(Singular)))then

         ! write(*,*)'gesvd wrong in LinearSolve, switching to QR'

         ! ! call GetRank(m,n,Atmp,rank,Rank_detection_factor)
         ! ! write(*,*)rank,'kao kao'

         ! ! stop
         ! Atmp = A

         ! ! RRQR
         ! jpvt = 0
         ! call geqp3f90(Atmp,jpvt,tau)
         ! if(myisnan(fnorm(Atmp,m,n)))then
         ! write(*,*)'Q or R has NAN in LinearSolve'
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
         ! if(abs(A_tmp(i,i))<BPACK_SafeUnderflow)rank = i -1
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

         ! if(myisnan(fnorm(x,n,k)))then
         ! write(*,*)'trisolve has NAN in LinearSolve'
         ! stop
         ! end if

         ! deallocate(A_tmp_rank,xtmp_rank)
         ! else
         if (Singular(1) < BPACK_SafeUnderflow) then
            if(verbose>=2)write (*, *) 'warning: Matrix zero in LinearSolve'
            rank = 1
            x = 0
         else
            rank = n
            do i = 1, n
               if (Singular(i)/Singular(1) <= eps_r) then
                  rank = i
                  if (Singular(i) < Singular(1)*eps_r/10) rank = i - 1
                  exit
               end if
            end do

            allocate (UU_h(rank, m))
            allocate (VV_h(n, rank))
            allocate (matrixtemp(rank, k))

            do ii = 1, rank
            do jj = 1, m
               UU_h(ii, jj) = conjg(cmplx(UU(jj, ii), kind=8))
            end do
            end do

            do ii = 1, n
            do jj = 1, rank
               VV_h(ii, jj) = conjg(cmplx(VV(jj, ii), kind=8))/Singular(jj)
            end do
            end do

            call gemmf90(UU_h, rank, b, m, matrixtemp, rank, 'N', 'N', rank, k, m, alpha, beta, flop)
            if (present(Flops)) Flops = Flops + flop
            call gemmf90(VV_h, n, matrixtemp, rank, x, n, 'N', 'N', n, k, rank, alpha, beta, flop)
            if (present(Flops)) Flops = Flops + flop

            deallocate (UU_h, VV_h, matrixtemp)
         end if
         ! end if

         do ii = 1, k
         if (myisnan(abs(sum(x(:, ii))))) then
            ! do jj =1,rank
            ! write(*,*)jj,'hh',A_tmp_rank(jj,:)
            ! end do
            write (*, *) 'hh', rank, Singular, fnorm(A, m, n)
            stop
         end if
         end do

         ! deallocate(A_tmp_rank,xtmp_rank)

         ! A = Atmp

      end if

      deallocate (Singular)
      deallocate (xtmp)
      deallocate (tau)
      deallocate (jpvt)
      deallocate (A_tmp)
      deallocate (Atmp)
      deallocate (UU)
      deallocate (VV)

   end subroutine LinearSolve

   subroutine GeneralInverse(m, n, A, A_inv, eps_r, Flops)
      !
      !
      implicit none

      integer m, n, mn_min
      DT:: A(:, :), A_inv(:, :)
      DTR, allocatable::Singular(:)
      DT, allocatable:: Atmp(:, :), tau(:), UU(:, :), VV(:, :), UU_h(:, :), VV_h(:, :), matrixtemp(:, :), A_tmp_rank(:, :), xtmp_rank(:, :), xtmp_rank3(:, :), A_tmp_rank2(:, :), xtmp_rank2(:, :)
      real(kind=8):: eps_r
      integer ii, jj, i, j, flag0, rank
      integer, allocatable:: jpvt(:)
      real(kind=8), optional::Flops
      real(kind=8)::flop

      DT:: alpha, beta
      alpha = 1d0
      beta = 0d0

      if (present(Flops)) Flops = 0

      A_inv = 0

      allocate (Atmp(m, n))
      Atmp = 0
      ! SVD
      Atmp = A
      mn_min = min(m, n)
      allocate (Singular(mn_min))
      allocate (UU(m, mn_min))
      allocate (VV(mn_min, n))

      Singular = 0
      UU = 0
      VV = 0

      call gesvd_robust(Atmp, Singular, UU, VV, m, n, mn_min, flop=flop)
      if (present(Flops)) Flops = Flops + flop

      if (Singular(1) < BPACK_SafeUnderflow) then
         write (*, *) 'Warning: A zero in GeneralInverse'
         A_inv = 0
      else
         rank = mn_min
         do i = 1, mn_min
            if (Singular(i)/Singular(1) <= eps_r) then
               rank = i
               if (Singular(i) < Singular(1)*eps_r/10) rank = i - 1
               exit
            end if
         end do

         allocate (UU_h(rank, m))
         allocate (VV_h(n, rank))
         UU_h = 0
         VV_h = 0

         do ii = 1, rank
         do jj = 1, m
            UU_h(ii, jj) = conjg(cmplx(UU(jj, ii), kind=8))
         end do
         end do

         do ii = 1, n
         do jj = 1, rank
            VV_h(ii, jj) = conjg(cmplx(VV(jj, ii), kind=8))/Singular(jj)
         end do
         end do

         call gemmf90(VV_h, n, UU_h, rank, A_inv, n, 'N', 'N', n, m, rank, alpha, beta, flop=flop)
         if (present(Flops)) Flops = Flops + flop

         deallocate (UU_h, VV_h)
      endif

      deallocate (Singular)
      ! deallocate(jpvt)
      deallocate (UU)
      deallocate (VV)
      deallocate (Atmp)

   end subroutine GeneralInverse



   subroutine PGeneralInverse(m, n, A, A_inv, eps_r, ctxt, Flops)
      implicit none
      integer m,n,mn,rank, iproc, myi,ii,jj
      DT::A(:,:), A_inv(:,:)
      DT,allocatable::UU(:,:),VV(:,:)
      real(kind=8)::eps_r
      integer ctxt
      real(kind=8),optional::Flops
      real(kind=8)::flop

      integer:: nprow, npcol, myrow, mycol, myArows, myAcols, info
      integer::descsMatA(9), descsMatAinv(9), descsMatU(9), descsMatV(9)
      DTR, allocatable :: Singular(:)

      if (present(Flops))Flops=0

      call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
      if (myrow /= -1 .and. mycol /= -1) then
         myArows = numroc_wp(m, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(n, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsMatA, m, n, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)

         myArows = numroc_wp(n, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(m, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsMatAinv, n, m, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)

         mn = min(m, n)
         myArows = numroc_wp(m, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(mn, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsMatU, m, mn, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (UU(max(1,myArows), max(1,myAcols)))
         myArows = numroc_wp(mn, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(n, nbslpk, mycol, 0, npcol)
         call descinit_wp(descsMatV, mn, n, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
         allocate (VV(max(1,myArows), max(1,myAcols)))
         allocate (Singular(mn))

         call PSVD_Truncate(m, n, A, descsMatA, UU, VV, descsMatU, descsMatV, Singular, eps_r, rank, ctxt, BPACK_SafeUnderflow, flop=flop)
         if (present(Flops))Flops=Flops+flop

         do ii=1,size(UU,1)
         do jj=1,size(UU,2)
            UU(ii, jj) = conjg(cmplx(UU(ii, jj), kind=8))
         enddo
         enddo
         do ii=1,size(VV,1)
         do jj=1,size(VV,2)
            VV(ii, jj) = conjg(cmplx(VV(ii, jj), kind=8))
         enddo
         enddo
         do ii = 1, rank
            call g2l(ii, rank, nprow, nbslpk, iproc, myi)
            if (iproc == myrow) then
               VV(myi, :) = VV(myi, :)/Singular(ii)
            endif
         enddo
         deallocate (Singular)


         A_inv=0
         call pgemmf90('T', 'T', n, m, rank, BPACK_cone, VV, 1, 1, descsMatV, UU, 1, 1, descsMatU, BPACK_czero, A_inv, 1, 1, descsMatAinv, flop=flop)
         if (present(Flops))Flops=Flops+flop

         deallocate(UU,VV)
      endif

   end subroutine PGeneralInverse



   subroutine RandomizedSVD(matRcol, matZRcol, matRrow, matZcRrow, matU, matV, Singular, rankmax_r, rankmax_c, rmax, rank, tolerance, SVD_tolerance, Flops)
      !
      !

      implicit none

      integer i, j, ii, jj, indx, rank_1, rank_2
      integer blocks, index_j, group_nn, rank_blocks
      integer group_m, group_n, size_of_groupm, size_of_groupn
      real(kind=8) norm_Z, norm_U, norm_V, tolerance, SVD_tolerance
      integer edgefine_m, edgefine_n, level_butterfly, level_blocks
      integer edge_m, edge_n, header_m, header_n
      integer rank, rank1, rank2, rank12, ranknew, row, column, rankmax, rankmax_c, rankmax_r, rankmax_min, rmax, idxs_r, idxs_c
      DT value_Z, value_UV, maxvalue
      DT inner_U, inner_V, ctemp
      real(kind=8) inner_UV
      integer, allocatable:: select_column(:), select_row(:)
      DT::matU(rankmax_r, rmax), matV(rmax, rankmax_c), matRcol(rankmax_c, rmax), matZRcol(rankmax_r, rmax), matRrow(rankmax_r, rmax), matZcRrow(rankmax_c, rmax)
      DT, allocatable::matZRcol1(:, :), matZcRrow1(:, :), tau(:)
      DT, allocatable::QQ1(:, :), RR1(:, :), QQ2(:, :), RR2(:, :), RrowcZRcol(:, :)
      DTR::Singular(rmax)
      DT, allocatable:: row_R(:), column_R(:), matM(:, :), matrixtemp(:, :)
      real(kind=8), allocatable:: norm_row_R(:), norm_column_R(:)

      DT, allocatable :: RrowcQ1(:, :), RrowcQ1inv(:, :), Q2cRcol(:, :), Q2cRcolinv(:, :), QQ2tmp(:, :), RR2tmp(:, :), UUsml(:, :), VVsml(:, :), tau_Q(:), mattemp(:, :), matU1(:, :), matV1(:, :)
      DTR, allocatable :: Singularsml(:)
      integer, allocatable:: jpvt(:)
      real(kind=8), optional::Flops
      real(kind=8)::flop

      if (present(Flops)) Flops = 0

      allocate (matZRcol1(rankmax_r, rmax))
      allocate (matZcRrow1(rankmax_c, rmax))
      allocate (tau(rmax))
      allocate (jpvt(rmax))
      allocate (QQ1(rankmax_r, rmax))
      allocate (RR1(rmax, rmax))
      allocate (QQ2(rankmax_c, rmax))
      allocate (RR2(rmax, rmax))
      allocate (RrowcZRcol(rmax, rmax))

      rankmax_min = min(rankmax_r, rankmax_c)
      matU = 0
      matV = 0
      Singular = 0

      QQ1 = matZRcol
      jpvt = 0
      tau = 0
      call geqp3f90(QQ1, jpvt, tau, flop=flop)
      if (present(Flops)) Flops = Flops + flop
      !  call geqrff90(QQ1,tau)
      RR1 = 0d0
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, rmax
         do i = 1, j
            RR1(i, j) = QQ1(i, j)
         enddo
      enddo
      ! !$omp end parallel do
      call un_or_gqrf90(QQ1, tau, rankmax_r, rmax, rmax, flop=flop)
      if (present(Flops)) Flops = Flops + flop

      ! write(*,*)fnorm(QQ1,rankmax_r,rmax),rankmax_r,rmax,fnorm(RR1,rmax,rmax),rmax,rmax,'really'

      if (abs(RR1(1, 1)) < BPACK_SafeUnderflow) then
         rank1 = 1
      else
         rank1 = rmax
         do i = 1, rmax
            if (abs(RR1(i, i))/abs(RR1(1, 1))/rankmax_r <= tolerance) then   ! be careful with the tolerance here
               rank1 = i
               if (abs(RR1(i, i)) < BPACK_SafeUnderflow) rank1 = i - 1
               exit
            end if
         end do
      endif
! write(*,*)shape(QQ2),shape(matZcRrow)

      QQ2 = matZcRrow
      jpvt = 0
      tau = 0
      call geqp3f90(QQ2, jpvt, tau, flop=flop)
      if (present(Flops)) Flops = Flops + flop
      !  call geqrff90(QQ2,tau)
      RR2 = 0d0
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, rmax
         do i = 1, j
            RR2(i, j) = QQ2(i, j)
         enddo
      enddo
      ! !$omp end parallel do
      call un_or_gqrf90(QQ2, tau, rankmax_c, rmax, rmax, flop=flop)
      if (present(Flops)) Flops = Flops + flop

      if (abs(RR2(1, 1)) < BPACK_SafeUnderflow) then
         rank2 = 1
      else
         rank2 = rmax
         do i = 1, rmax
            if (abs(RR2(i, i))/abs(RR2(1, 1))/rankmax_c <= tolerance) then  ! be careful with the tolerance here
               rank2 = i
               if (abs(RR2(i, i)) < BPACK_SafeUnderflow) rank2 = i - 1
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

      allocate (RrowcQ1(rmax, rank1))
      RrowcQ1 = 0
      allocate (RrowcQ1inv(rank1, rmax))
      RrowcQ1inv = 0
      ! call gemmHN_omp(matRrow,QQ1(1:rankmax_r,1:rank1),RrowcQ1,rmax,rank1,rankmax_r)
      call gemmf90(matRrow, rankmax_r, QQ1, rankmax_r, RrowcQ1, rmax, 'C', 'N', rmax, rank1, rankmax_r, BPACK_cone, BPACK_czero, flop=flop)
      if (present(Flops)) Flops = Flops + flop

      call GeneralInverse(rmax, rank1, RrowcQ1, RrowcQ1inv, tolerance, Flops=flop)
      if (present(Flops)) Flops = Flops + flop
      deallocate (RrowcQ1)

      ! allocate(ipiv(rmax))
      ! call getrff90(RrowcQ1,ipiv)
      ! call getrif90(RrowcQ1,ipiv)
      ! RrowcQ1inv = RrowcQ1
      ! deallocate(ipiv)
      ! deallocate(RrowcQ1)

      allocate (Q2cRcol(rank2, rmax))
      Q2cRcol = 0
      allocate (Q2cRcolinv(rmax, rank2))
      Q2cRcolinv = 0
      ! call gemmHN_omp(QQ2(1:rankmax_c,1:rank2),matRcol,Q2cRcol,rank2,rmax,rankmax_c)
      call gemmf90(QQ2, rankmax_c, matRcol, rankmax_c, Q2cRcol, rank2, 'C', 'N', rank2, rmax, rankmax_c, BPACK_cone, BPACK_czero, flop=flop)
      if (present(Flops)) Flops = Flops + flop
      call GeneralInverse(rank2, rmax, Q2cRcol, Q2cRcolinv, tolerance, Flops=flop)
      if (present(Flops)) Flops = Flops + flop
      deallocate (Q2cRcol)

      ! allocate(ipiv(rmax))
      ! call getrff90(Q2cRcol,ipiv)
      ! call getrif90(Q2cRcol,ipiv)
      ! Q2cRcolinv = Q2cRcol
      ! deallocate(ipiv)
      ! deallocate(Q2cRcol)

      ! call gemmHN_omp(matRrow,matZRcol,RrowcZRcol,rmax,rmax,rankmax_r)
      call gemmf90(matRrow, rankmax_r, matZRcol, rankmax_r, RrowcZRcol, rmax, 'C', 'N', rmax, rmax, rankmax_r, BPACK_cone, BPACK_czero, flop=flop)
      if (present(Flops)) Flops = Flops + flop

      allocate (matrixtemp(rmax, rank2))
      matrixtemp = 0
      allocate (matM(rank1, rank2))
      matM = 0
      ! call gemm_omp(RrowcZRcol,Q2cRcolinv,matrixtemp,rmax,rank2,rmax)
      call gemmf90(RrowcZRcol, rmax, Q2cRcolinv, rmax, matrixtemp, rmax, 'N', 'N', rmax, rank2, rmax, BPACK_cone, BPACK_czero, flop=flop)
      if (present(Flops)) Flops = Flops + flop

      ! call gemm_omp(RrowcQ1inv,matrixtemp,matM,rank1,rank2,rmax)
      call gemmf90(RrowcQ1inv, rank1, matrixtemp, rmax, matM, rank1, 'N', 'N', rank1, rank2, rmax, BPACK_cone, BPACK_czero, flop=flop)
      if (present(Flops)) Flops = Flops + flop

      deallocate (matrixtemp, RrowcQ1inv, Q2cRcolinv)

      ! allocate(matrixtemp(rankmax_r,rank2))
      ! allocate(matM(rank1,rank2))
      ! call gemm_omp(matZRcol,Q2cRcolinv,matrixtemp,rankmax_r,rank2,rmax)
      ! call gemmHN_omp(QQ1(1:rankmax_r,1:rank1),matrixtemp,matM,rank1,rankmax_r,rank2)
      ! deallocate(matrixtemp,RrowcQ1inv,Q2cRcolinv)

      rank12 = min(rank1, rank2)
      allocate (UUsml(rank1, rank12), VVsml(rank12, rank2), Singularsml(rank12))
      UUsml = 0
      VVsml = 0
      Singularsml = 0
      call SVD_Truncate(matM, rank1, rank2, rank12, UUsml, VVsml, Singularsml, SVD_tolerance, BPACK_SafeUnderflow, rank, flop=flop)
      if (present(Flops)) Flops = Flops + flop

      ! write(111,*)UUsml(1:rank1,1:rank)
      ! stop

      ! call gemm_omp(QQ1(1:rankmax_r,1:rank1),UUsml(1:rank1,1:rank),matU(1:rankmax_r,1:rank),rankmax_r,rank,rank1)
      call gemmf90(QQ1, rankmax_r, UUsml, rank1, matU, rankmax_r, 'N', 'N', rankmax_r, rank, rank1, BPACK_cone, BPACK_czero, flop=flop)
      if (present(Flops)) Flops = Flops + flop

      ! call gemmNH_omp(VVsml(1:rank,1:rank2),QQ2(1:rankmax_c,1:rank2),matV(1:rank,1:rankmax_c),rank,rankmax_c,rank2)
      call gemmf90(VVsml, rank12, QQ2, rankmax_c, matV, rmax, 'N', 'H', rank, rankmax_c, rank2, BPACK_cone, BPACK_czero, flop=flop)
      if (present(Flops)) Flops = Flops + flop

      Singular(1:rank) = Singularsml(1:rank)
      deallocate (UUsml, VVsml, Singularsml, matM)

      deallocate (matZRcol1)
      deallocate (matZcRrow1)
      deallocate (tau)
      deallocate (jpvt)
      deallocate (QQ1)
      deallocate (RR1)
      deallocate (QQ2)
      deallocate (RR2)
      deallocate (RrowcZRcol)

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
      ! call SVD_Truncate(matM,rankmax_r,rankmax_c,rank12,UUsml,VVsml,Singularsml,SVD_tolerance,BPACK_SafeUnderflow, rank)
      ! matU(1:rankmax_r,1:rank) = UUsml(1:rankmax_r,1:rank)
      ! matV(1:rank,1:rankmax_c) = VVsml(1:rank,1:rankmax_c)
      ! Singular(1:rank) = Singularsml(1:rank)
      ! deallocate(UUsml,VVsml,Singularsml,matM)

      return

   end subroutine RandomizedSVD

   subroutine RandomSubMat(ms, me, ns, ne, k, A, Oflag)
      implicit none
      integer ms, me, ns, ne, k, m, n
      DT:: A(:, :)
      DT, allocatable::matrix_small(:, :)
      integer:: Oflag

      m = me - ms + 1
      n = ne - ns + 1
      call assert(m > 0, 'm<=0 in RandomSubMat')
      call assert(n > 0, 'n<=0 in RandomSubMat')

      allocate (matrix_small(m, n))
      call RandomMat(m, n, k, matrix_small, Oflag)
      A(ms:me, ns:ne) = matrix_small
      deallocate (matrix_small)

   end subroutine RandomSubMat

   subroutine RandomMat(m, n, k, A, Oflag)

      !
      !
#ifdef USEVSL
      use mkl_vsl_type
      use mkl_vsl
#endif

      implicit none

      integer m, n, k, mn_min, ktmp
      DT:: A(:, :)
      real(kind=8):: c
      DTR, allocatable::Singular(:), Ar(:, :), Ai(:, :)
      complex(kind=8):: ctemp
      complex(kind=8), allocatable:: UU(:, :), VV(:, :)
      integer ii, jj, kk, i, j, flag0, rank
      integer:: Oflag
#ifdef USEVSL
      type(VSL_STREAM_STATE) ::stream(2)
#endif
      integer brng, method, seedd, ierror

      ktmp = k
      if (ktmp > min(m, n)) then
         ! write(*,*)'k is not properly set in RandomMat'
         ktmp = min(m, n)
      end if
      ! call assert(m<=n,'m>n')

      call assert(k <= min(m, n), 'k too large in RandomMat')

#ifdef USEVSL

#if DAT==0
         allocate (Ar(m, n))
         allocate (Ai(m, n))
         brng = VSL_BRNG_MCG31
         method = VSL_RNG_METHOD_UNIFORM_STD
         call random_number(c)
         seedd = NINT(1000*c)
         ierror = vslnewstream(stream(1), brng, seedd)
         call random_number(c)
         seedd = NINT(1000*c)
         ierror = vslnewstream(stream(2), brng, seedd)
         ierror = vdrnguniform(method, stream(1), M*N, Ar, -1d0, 1d0)
         ierror = vdrnguniform(method, stream(2), M*N, Ai, -1d0, 1d0)
         ! ierror=vdrnggaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream, M*N, mati, 1d0, 1d0)
         ierror = vsldeletestream(stream(1))
         ierror = vsldeletestream(stream(2))
         A = Ar + BPACK_junit*Ai
         deallocate (Ar)
         deallocate (Ai)
#elif DAT==1
         allocate (Ar(m, n))
         brng = VSL_BRNG_MCG31
         method = VSL_RNG_METHOD_UNIFORM_STD
         call random_number(c)
         seedd = NINT(1000*c)
         ierror = vslnewstream(stream(1), brng, seedd)
         call random_number(c)
         seedd = NINT(1000*c)
         ierror = vsrnguniform(method, stream(1), M*N, Ar, -1d0, 1d0)
         ierror = vsldeletestream(stream(1))
         A = Ar
         deallocate (Ar)
#elif DAT==2
         allocate (Ar(m, n))
         allocate (Ai(m, n))
         brng = VSL_BRNG_MCG31
         method = VSL_RNG_METHOD_UNIFORM_STD
         call random_number(c)
         seedd = NINT(1000*c)
         ierror = vslnewstream(stream(1), brng, seedd)
         call random_number(c)
         seedd = NINT(1000*c)
         ierror = vslnewstream(stream(2), brng, seedd)
         ierror = vsrnguniform(method, stream(1), M*N, Ar, -1d0, 1d0)
         ierror = vsrnguniform(method, stream(2), M*N, Ai, -1d0, 1d0)
         ! ierror=vdrnggaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream, M*N, mati, 1d0, 1d0)
         ierror = vsldeletestream(stream(1))
         ierror = vsldeletestream(stream(2))
         A = Ar + BPACK_junit*Ai
         deallocate (Ar)
         deallocate (Ai)
#elif DAT==3
         allocate (Ar(m, n))
         brng = VSL_BRNG_MCG31
         method = VSL_RNG_METHOD_UNIFORM_STD
         call random_number(c)
         seedd = NINT(1000*c)
         ierror = vslnewstream(stream(1), brng, seedd)
         call random_number(c)
         seedd = NINT(1000*c)
         ierror = vdrnguniform(method, stream(1), M*N, Ar, -1d0, 1d0)
         ierror = vsldeletestream(stream(1))
         A = Ar
         deallocate (Ar)
#endif
#else
      do ii = 1, m
      do jj = 1, n
         call random_dp_number(A(ii, jj))
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
      ! ! !!!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
      ! do jj=1, n
      ! do ii=1, m
      ! ctemp=(0d0,0d0)
      ! do kk=1, k
      ! ctemp=ctemp+UU(ii,kk)*VV(kk,jj)*Singular(kk)
      ! enddo
      ! A(ii,jj)=ctemp
      ! enddo
      ! enddo
      ! ! !!!$omp end parallel do
      ! ! end if

      ! deallocate(Singular)
      ! deallocate(UU)
      ! deallocate(VV)
! end select

#endif

   end subroutine RandomMat

   subroutine ID_Selection(Mat, select_column, select_row, m, n, rank, tolerance)


      implicit none

      integer i, j, ii, jj, indx, rank_1, rank_2
      integer blocks, index_j, group_nn, rank_blocks
      integer group_m, group_n, size_of_groupm, size_of_groupn
      real(kind=8) tolerance
      integer edgefine_m, edgefine_n, level_butterfly, level_blocks
      integer edge_m, edge_n, header_m, header_n
      integer rank, row, column, rankmax, n, m, rank_c, rank_r
      DT value_Z, value_UV, maxvalue, Mat(m, n)
      DT, allocatable::Mat1(:, :), Mat1T(:, :)
      DT inner_U, inner_V, ctemp
      real(kind=8) inner_UV, flop
      integer select_column(n), select_row(m)

      DT, allocatable:: row_R(:), column_R(:), matrixtemp_V(:, :), matrixtemp_U(:, :), tau(:)
      real(kind=8), allocatable:: norm_row_R(:), norm_column_R(:)
      integer, allocatable :: jpvt(:)

      select_column = 0
      select_row = 0

      allocate (Mat1(m, n))
      Mat1 = Mat
      allocate (Mat1T(n, m))
      call copymatT(Mat, Mat1T, m, n)

      allocate (jpvt(max(m, n)))
      allocate (tau(max(m, n)))

      jpvt = 0
      call geqp3modf90(Mat1T, jpvt, tau, tolerance, BPACK_SafeUnderflow, rank_r)
      select_row(1:rank_r) = jpvt(1:rank_r)
      if (rank_r == 0) then
         rank_r = 1
         select_row(1) = 1
      endif

      jpvt = 0
      call geqp3modf90(Mat1, jpvt, tau, tolerance, BPACK_SafeUnderflow, rank_c)
      select_column(1:rank_c) = jpvt(1:rank_c)
      if (rank_c == 0) then
         rank_c = 1
         select_column(1) = 1
      endif

      rank = min(rank_c, rank_r)

      deallocate (jpvt)
      deallocate (tau)
      deallocate (Mat1)
      deallocate (Mat1T)

      return

   end subroutine ID_Selection

   subroutine ACA_CompressionFull(mat, matU, matV, rankmax_r, rankmax_c, rmax, rank, tolerance, SVD_tolerance)


      implicit none

      integer i, j, ii, jj, indx, rank_1, rank_2
      integer blocks, index_j, group_nn, rank_blocks
      integer group_m, group_n, size_of_groupm, size_of_groupn
      real(kind=8) norm_Z, norm_U, norm_V, tolerance, SVD_tolerance
      integer edgefine_m, edgefine_n, level_butterfly, level_blocks
      integer edge_m, edge_n, header_m, header_n
      integer rank, ranknew, row, column, rankmax, rankmax_c, rankmax_r, rankmax_min, rmax
      DT value_Z, value_UV, maxvalue
      DT inner_U, inner_V, ctemp
      real(kind=8) inner_UV
      integer:: select_column(rmax+1), select_row(rmax+1)
      DT::mat(rankmax_r, rankmax_c), matU(rankmax_r, rmax), matV(rmax, rankmax_c)

      DT, allocatable:: row_R(:), column_R(:)
      real(kind=8), allocatable:: norm_row_R(:), norm_column_R(:)

      DT, allocatable :: QQ1(:, :), RR1(:, :), QQ2(:, :), RR2(:, :), QQ2tmp(:, :), RR2tmp(:, :), UUsml(:, :), VVsml(:, :), tau_Q(:), mattemp(:, :), matU1(:, :), matV1(:, :)
      DTR, allocatable :: Singularsml(:)

      rankmax_min = min(rankmax_r, rankmax_c)
      norm_Z = 0
      select_column = 0
      select_row = 0

      allocate (row_R(rankmax_c), column_R(rankmax_r))
      allocate (norm_row_R(rankmax_c), norm_column_R(rankmax_r))

      select_row(1) = 1

      ! !$omp parallel do default(shared) private(j,value_Z)
      do j = 1, rankmax_c
         value_Z = mat(select_row(1), j)
         row_R(j) = value_Z
         norm_row_R(j) = dble(value_Z*conjg(cmplx(value_Z, kind=8)))
      enddo
      ! !$omp end parallel do

      select_column(1) = maxloc(norm_row_R, 1)
      maxvalue = row_R(select_column(1))

      ! !$omp parallel do default(shared) private(j)
      do j = 1, rankmax_c
         row_R(j) = row_R(j)/maxvalue
      enddo
      ! !$omp end parallel do
      ! !$omp parallel do default(shared) private(j)
      do j = 1, rankmax_c
         matV(1, j) = row_R(j)
      enddo
      ! !$omp end parallel do

      ! !$omp parallel do default(shared) private(i,value_Z)
      do i = 1, rankmax_r
         value_Z = mat(i, select_column(1))
         column_R(i) = value_Z
         norm_column_R(i) = dble(value_Z*conjg(cmplx(value_Z, kind=8)))
      enddo
      ! !$omp end parallel do

      norm_column_R(select_row(1)) = 0

      ! !$omp parallel do default(shared) private(i)
      do i = 1, rankmax_r
         matU(i, 1) = column_R(i)
      enddo
      ! !$omp end parallel do

      norm_U = norm_vector(column_R, rankmax_r)
      norm_V = norm_vector(row_R, rankmax_c)
      norm_Z = norm_Z + norm_U*norm_V

      ! if(rankmax<2)write(*,*)'rankmax'
      select_row(2) = maxloc(norm_column_R, 1)

      rank = 1
      ! write(*,*)column_R,row_R
      ! write(*,*)norm_Z*ACA_tolerance_forward**2,norm_U*norm_V,'hehe',ACA_tolerance_forward
      do while (norm_Z*tolerance**2 < norm_U*norm_V .and. rank < rankmax_min)

         ! !$omp parallel do default(shared) private(j,i,value_Z,value_UV)
         do j = 1, rankmax_c
            value_Z = mat(select_row(rank + 1), j)
            value_UV = 0
            do i = 1, rank
               value_UV = value_UV + matU(select_row(rank + 1), i)*matV(i, j)
            enddo
            row_R(j) = value_Z - value_UV
            norm_row_R(j) = dble(row_R(j)*conjg(cmplx(row_R(j), kind=8)))
         enddo
         ! !$omp end parallel do

         do i = 1, rank
            norm_row_R(select_column(i)) = 0
         enddo

         select_column(rank + 1) = maxloc(norm_row_R, 1)
         maxvalue = row_R(select_column(rank + 1))

         ! !$omp parallel do default(shared) private(j)
         do j = 1, rankmax_c
            row_R(j) = row_R(j)/maxvalue
         enddo
         ! !$omp end parallel do
         ! !$omp parallel do default(shared) private(j)
         do j = 1, rankmax_c
            matV(rank + 1, j) = row_R(j)
         enddo
         ! !$omp end parallel do

         ! !$omp parallel do default(shared) private(i,j,value_Z,value_UV)
         do i = 1, rankmax_r
            value_Z = mat(i, select_column(rank + 1))
            value_UV = 0
            do j = 1, rank
               value_UV = value_UV + matU(i, j)*matV(j, select_column(rank + 1))
            enddo
            column_R(i) = value_Z - value_UV
            norm_column_R(i) = dble(column_R(i)*conjg(cmplx(column_R(i), kind=8)))
         enddo
         ! !$omp end parallel do

         do i = 1, rank + 1
            norm_column_R(select_row(i)) = 0
         enddo

         ! !$omp parallel do default(shared) private(i)
         do i = 1, rankmax_r
            matU(i, rank + 1) = column_R(i)
         enddo
         ! !$omp end parallel do

         norm_U = norm_vector(column_R, rankmax_r)
         norm_V = norm_vector(row_R, rankmax_c)
         inner_UV = 0
         do j = 1, rank
            inner_U = 0
            inner_V = 0
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_U)
            do i = 1, rankmax_r
               ctemp = matU(i, rank + 1)*matU(i, j)
               inner_U = inner_U + ctemp*conjg(cmplx(ctemp, kind=8))
            enddo
            ! !$omp end parallel do
            ! !$omp parallel do default(shared) private(i,ctemp) reduction(+:inner_V)
            do i = 1, rankmax_c
               ctemp = matV(rank + 1, i)*matV(j, i)
               inner_V = inner_V + ctemp*conjg(cmplx(ctemp, kind=8))
            enddo
            ! !$omp end parallel do
            inner_UV = inner_UV + dble(2*sqrt(inner_U*inner_V))
         enddo
         norm_Z = norm_Z + inner_UV + norm_U*norm_V

         rank = rank + 1
         if (rank > rmax) then
            write (*, *) 'increase rmax', rank, rmax
            stop
         end if
         if (rank < rankmax_min) then
            select_row(rank + 1) = maxloc(norm_column_R, 1)
         endif

      enddo

      ! write(*,*)select_row(1:rank),select_column(1:rank)

      deallocate (row_R, column_R)
      deallocate (norm_row_R, norm_column_R)

! ACA followed by SVD

      allocate (QQ1(rankmax_r, rank))
      ! call copymatN(matU(1:rankmax_r,1:rank),QQ1,rankmax_r,rank)
      QQ1 = matU(1:rankmax_r, 1:rank)
      allocate (tau_Q(rank))
      call geqrff90(QQ1, tau_Q)

      allocate (RR1(rank, rank))
      RR1 = 0d0
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, rank
         do i = 1, j
            RR1(i, j) = QQ1(i, j)
         enddo
      enddo
      ! !$omp end parallel do
      call un_or_gqrf90(QQ1, tau_Q, rankmax_r, rank, rank)
      deallocate (tau_Q)

      allocate (QQ2tmp(rankmax_c, rank))
      call copymatT(matV(1:rank, 1:rankmax_c), QQ2tmp, rank, rankmax_c)
      allocate (tau_Q(rank))
      call geqrff90(QQ2tmp, tau_Q)

      allocate (RR2tmp(rank, rank))
      RR2tmp = 0d0
      ! !$omp parallel do default(shared) private(i,j)
      do j = 1, rank
         do i = 1, j
            RR2tmp(i, j) = QQ2tmp(i, j)
         enddo
      enddo
      ! !$omp end parallel do
      call un_or_gqrf90(QQ2tmp, tau_Q, rankmax_c, rank, rank)
      deallocate (tau_Q)

      allocate (QQ2(rank, rankmax_c))
      call copymatT(QQ2tmp, QQ2, rankmax_c, rank)
      allocate (RR2(rank, rank))
      call copymatT(RR2tmp, RR2, rank, rank)

      ! allocate(matU1(rankmax_r,rank))
      ! allocate(matV1(rank,rankmax_c))
      ! call gemm_omp(QQ1,RR1,matU1,rankmax_r,rank,rank)

      ! call gemm_omp(RR2,QQ2,matV1,rank,rankmax_c,rank)

      ! write(*,*)fnorm(matU1-matU(1:rankmax_r,1:rank),rankmax_r,rank),fnorm(matV1-matV(1:rank,1:rankmax_c),rank,rankmax_c)

      deallocate (QQ2tmp, RR2tmp)
      allocate (mattemp(rank, rank))
      ! call gemm_omp(RR1,RR2,mattemp,rank,rank,rank)
      call gemmf90(RR1, rank, RR2, rank, mattemp, rank, 'N', 'N', rank, rank, rank, BPACK_cone, BPACK_czero)

      allocate (UUsml(rank, rank), VVsml(rank, rank), Singularsml(rank))
      call SVD_Truncate(mattemp, rank, rank, rank, UUsml, VVsml, Singularsml, SVD_tolerance, BPACK_SafeUnderflow, ranknew)

      ! call gemm_omp(QQ1,UUsml(1:rank,1:ranknew),matU(1:rankmax_r,1:ranknew),rankmax_r,ranknew,rank)
      call gemmf90(QQ1, rankmax_r, UUsml, rank, matU, rankmax_r, 'N', 'N', rankmax_r, ranknew, rank, BPACK_cone, BPACK_czero)
      ! call gemm_omp(VVsml(1:ranknew,1:rank),QQ2,matV(1:ranknew,1:rankmax_c),ranknew,rankmax_c,rank)
      call gemmf90(VVsml, rank, QQ2, rank, matV, rmax, 'N', 'N', ranknew, rankmax_c, rank, BPACK_cone, BPACK_czero)
      ! write(*,*)'aca rank:',rank,'after svd',ranknew

      rank = ranknew
      do i = 1, ranknew
         matU(1:rankmax_r, i) = matU(1:rankmax_r, i)*Singularsml(i)
      end do
      deallocate (mattemp, RR1, RR2, QQ1, QQ2, UUsml, VVsml, Singularsml)

      return

   end subroutine ACA_CompressionFull

   recursive subroutine RRQR_LQ(mat, mm, nn, mn, UU, VV, tolerance, rank, lr, flops)
!
!
      implicit none
      integer mm, nn, mn, rank, ii, jj, rank_new
      real(kind=8):: tolerance
      DT::mat(mm, nn), UU(mm, mn), VV(mn, nn)
      DT, allocatable::mat0(:, :), UUt(:, :), VVt(:, :), tau(:), RR1(:, :)
      DTR:: Singular(mn)
      integer::i, j, flag
      real(kind=8), optional::flops
      real(kind=8)::flop
      character::lr
      integer, allocatable :: jpiv(:)

      if (lr == 'L') then ! LQ
         allocate (mat0(nn, mm))
         allocate (UUt(nn, mn))
         allocate (VVt(mn, mm))
         call copymatT(mat, mat0, mm, nn)
         call RRQR_LQ(mat0, nn, mm, mn, UUt, VVt, tolerance, rank, 'R', flops)
         call copymatT(UUt, VV, nn, mn)
         call copymatT(VVt, UU, mn, mm)
         deallocate (mat0)
         deallocate (UUt)
         deallocate (VVt)
      elseif (lr == 'R') then ! QR
         if (present(flops)) flops = 0
         allocate (mat0(mm, nn))
         mat0 = mat
         allocate (jpiv(nn))
         jpiv = 0
         allocate (tau(mn))
         tau = 0
         call geqp3modf90(mat0, jpiv, tau, tolerance, BPACK_SafeUnderflow, rank, flop=flop)
         if (present(flops)) flops = flops + flop
         if (rank > 0) then
            allocate (RR1(1:rank, nn))
            RR1 = 0
            do j = 1, nn
               do i = 1, min(j, rank)
                  RR1(i, j) = mat0(i, j)
               enddo
            enddo

            call un_or_gqrf90(mat0, tau, mm, rank, rank, flop=flop)
            if (present(flops)) flops = flops + flop
            UU(1:mm, 1:rank) = mat0(1:mm, 1:rank)
            do j = 1, nn
               VV(1:rank, jpiv(j)) = RR1(1:rank, j)
            enddo
            deallocate (RR1)
         else
            rank = 1
            UU = 0
            VV = 0
         endif

         deallocate (mat0)
         deallocate (jpiv)
         deallocate (tau)
      endif

   end subroutine RRQR_LQ

   subroutine SVD_Truncate(mat, mm, nn, mn, UU, VV, Singular, tolerance_rel, tolerance_abs, rank, flop)
!
!
      implicit none
      integer mm, nn, mn, rank, ii, jj
      real(kind=8):: tolerance_rel, tolerance_abs
      DT::mat(mm, nn), UU(mm, mn), VV(mn, nn)
      DT, allocatable::mat0(:, :)
      DTR:: Singular(mn)
      integer::i, flag
      real(kind=8), optional::flop

      UU = 0
      VV = 0
      Singular = 0

      allocate (mat0(mm, nn))
      mat0 = mat
      call gesvd_robust(mat0, Singular, UU, VV, mm, nn, mn, flop=flop)
      rank = mn

      if (Singular(1) < BPACK_SafeUnderflow) then
         rank = 1
         UU = 0
         VV = 0
         Singular = 0
      else
         rank = mn
         do i = 1, mn
            if (Singular(i) <= max(Singular(1)*tolerance_rel,tolerance_abs)) then
               rank = i
               if (Singular(i) < max(Singular(1)*tolerance_rel,tolerance_abs)/10) rank = i - 1
               exit
            end if
         end do
      endif
      deallocate (mat0)

   end subroutine SVD_Truncate



   subroutine RRQR_SVD(mat, mm, nn, mn, rmax, UU, VV, Singular, tolerance, rank, flop)
         implicit none
         integer mm, nn, mn, mns, rmax, rank, rank1, ii, jj
         real(kind=8):: tolerance
         DT::mat(mm, nn), UU(mm, rmax), VV(rmax, nn), UU1(mm, mn), VV1(mn, nn)
         DT, allocatable::mats(:, :),UUs(:,:),VVs(:,:)
         DTR:: Singular(rmax)
         integer::i, flag
         real(kind=8), optional::flop
         real(kind=8)::flop0

         call RRQR_LQ(mat, mm, nn, mn, UU1, VV1, tolerance, rank1, 'R', flops=flop)
         call assert(rank1<=rmax,'RRQR_SVD requires increasing rmax')
         allocate(mats(rank1,nn))
         mats = VV1(1:rank1,1:nn)
         mns = min(rank1,nn)
         allocate(UUs(1:rank1,1:mns))
         allocate(VVs(1:mns,1:nn))
         call SVD_Truncate(mats, rank1, nn, mns, UUs, VVs, Singular, tolerance, BPACK_SafeUnderflow, rank, flop=flop0)
         if (present(flop)) flop = flop + flop0

         VV(1:rank,1:nn) = VVs(1:rank,1:nn)

         UU=0
         call gemmf90(UU1, mm, UUs, rank1, UU, mm, 'N', 'N', mm, rank, rank1, BPACK_cone, BPACK_czero, flop=flop0)
         if (present(flop)) flop = flop + flop0


         deallocate (mats)
         deallocate (UUs)
         deallocate (VVs)

   end subroutine RRQR_SVD




   subroutine PSVD_Truncate(mm, nn, mat, descMat, UU, VV, descUU, descVV, Singular, tolerance, rank, ctxt, tolerance_abs,flop)
      implicit none
      integer mm, nn, mnmin, rank, ii, jj
      real(kind=8):: tolerance, tolerance_abs
      DT::mat(:, :), UU(:, :), VV(:, :)
      DT, allocatable::mat0(:, :)
      DTR:: Singular(:)
      integer::i, flag
      real(kind=8), optional::flop
      integer::descMat(9), descUU(9), descVV(9)
      integer::ctxt, nprow, npcol, myrow, mycol, iproc, myi

      call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
      mnmin = min(mm, nn)

      call pgesvdf90('V', 'V', mm, nn, mat, 1, 1, descMat, Singular, UU, 1, 1, descUU, VV, 1, 1, descVV, flop=flop)

      rank = mnmin
      if (Singular(1) > BPACK_SafeUnderflow) then
         rank = mnmin
         do i = 1, mnmin
            if (Singular(i)/Singular(1) <= tolerance .or. Singular(i)<=tolerance_abs) then
               rank = i
               if (Singular(i)<=tolerance_abs)exit
               if (Singular(i) < Singular(1)*tolerance/10) rank = i - 1
               exit
            end if
         end do
      else
         rank = 1
         Singular(1) = 0
         VV = 0
         UU = 0
      endif

   end subroutine PSVD_Truncate

! sort the first dimension of the array
   subroutine PIKSRT_DBLE_Multi(N, M, ARR)
      implicit none
      integer j, i, N, M
      real(kind=8) ARR(N, M)
      real(kind=8), allocatable::a(:)
      allocate (a(M))
      do j = 2, N
         a = ARR(j, :)
         do i = j - 1, 1, -1
            if (ARR(i, 1) <= a(1)) goto 10
            ARR(i + 1, :) = ARR(i, :)
         end do
         i = 0
10       ARR(i + 1, :) = a
      end do
      deallocate (a)
      return
   end subroutine PIKSRT_DBLE_Multi

! sort the first dimension of the array
   subroutine PIKSRT_INT_Multi(N, M, ARR)
      implicit none
      integer j, i, N, M
      integer ARR(N, M)
      integer, allocatable::a(:)
      allocate (a(M))
      do j = 2, N
         a = ARR(j, :)
         do i = j - 1, 1, -1
            if (ARR(i, 1) <= a(1)) goto 11
            ARR(i + 1, :) = ARR(i, :)
         end do
         i = 0
11       ARR(i + 1, :) = a
      end do
      deallocate (a)
      return
   end subroutine PIKSRT_INT_Multi

!>*** remove the duplicates in an integer array
   subroutine remove_dup_int(array, nin, nout)
      implicit none
      integer array(nin)
      real(kind=8) array1(nin)
      integer order(nin)
      integer nin, nout, last, ii

      array1(1:nin) = array(1:nin)
      call quick_sort(array1, order, nin)

      nout = 0
      last = 0
      do ii = 1, nin
         if (NINT(array1(ii)) /= last) then
            nout = nout + 1
            last = NINT(array1(ii))
            array(nout) = last
         endif
      enddo
   end subroutine remove_dup_int

   subroutine binary_search(N, x, val, mid)

      implicit none
      integer N
      integer x(N)
      integer :: range, start, finish, mid
      integer :: val

      start = 1
      finish = N
      range = finish - start
      mid = (start + finish)/2

      do while (x(mid) /= val .and. range > 0)
         if (val > x(mid)) then
            start = mid + 1
         else
            finish = mid - 1
         end if
         range = finish - start
         mid = (start + finish)/2
      end do

      if (x(mid) /= val) then
         mid = -1
      end if
   end subroutine binary_search

   subroutine quick_sort(list, order, n)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.

      IMPLICIT NONE
      integer n
      real(kind=8) list(n)
      INTEGER order(n)

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
            DO j = i + 1, right_end
               IF (list(i) > list(j)) THEN
                  temp = list(i); list(i) = list(j); list(j) = temp
                  itemp = order(i); order(i) = order(j); order(j) = itemp
               END IF
            END DO
         END DO

      END subroutine interchange_sort

   END subroutine quick_sort




   subroutine quick_sort_int(list, order, n)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.

      IMPLICIT NONE
      integer n
      integer list(n)
      INTEGER order(n)

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
         integer                :: reference, temp
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
         integer                :: temp

         DO i = left_end, right_end - 1
            DO j = i + 1, right_end
               IF (list(i) > list(j)) THEN
                  temp = list(i); list(i) = list(j); list(j) = temp
                  itemp = order(i); order(i) = order(j); order(j) = itemp
               END IF
            END DO
         END DO

      END subroutine interchange_sort

   END subroutine quick_sort_int




   !>******* convert from a gloal Cartesian coordinate to a local Cartesian coordinate (origin, xunit, yunit, zunit) and then convert to the local spherial coordinate
   subroutine Cart2Sph_Loc(xin, yin, zin, origin, xunit, yunit, zunit, r, theta, phi)
      implicit none
      real(kind=8), intent(in)::xin, yin, zin, origin(3),xunit(3),yunit(3),zunit(3)
      real(kind=8), intent(out)::r, theta, phi
      real(kind=8):: x, y, z, rho(3),origin1(3)

      origin1=0
      rho(1)=xin
      rho(2)=yin
      rho(3)=zin
      rho = rho - origin
      x = dot_product(rho,xunit)
      y = dot_product(rho,yunit)
      z = dot_product(rho,zunit)

      call Cart2Sph(x, y, z, origin1, r, theta, phi)

   end subroutine Cart2Sph_Loc


   subroutine Cart2Sph(xin, yin, zin, origin, r, theta, phi)
      implicit none
      real(kind=8), intent(in)::xin, yin, zin, origin(3)
      real(kind=8), intent(out)::r, theta, phi
      real(kind=8):: x, y, z

      x = xin - origin(1)
      y = yin - origin(2)
      z = zin - origin(3)

      r = sqrt(x**2 + y**2 + z**2)
      theta = acos(z/r)
      if (r == 0) theta = 0
      if (x == 0 .and. y >= 0) then
         phi = BPACK_pi/2
      else if (x == 0 .and. y < 0) then
         phi = 3*BPACK_pi/2

      else
         if (y == 0 .and. x >= 0) then
            phi = 0
         else if (y == 0 .and. x < 0) then
            phi = BPACK_pi
         else

            phi = atan(y/x)
            if (phi > 0 .and. x < 0) phi = phi + BPACK_pi
            if (phi < 0 .and. x > 0) phi = phi + 2*BPACK_pi
            if (phi < 0 .and. x < 0) phi = phi + BPACK_pi
         end if
      end if
   end subroutine Cart2Sph

   complex(kind=8) function Hankel02_Func(x)


      implicit none

      real(kind=8) x
      complex(kind=8) y

      Hankel02_Func = BesselJ0_func(x) - BPACK_junit*BesselY0_func(x)

      return

   end function Hankel02_Func

   real(kind=8) function BesselJ0_func(x)

      implicit none

      real(kind=8) x, z, ax
      real(kind=8) y, rtemp1, rtemp2, xx

      ax = abs(x)
      if (ax < 8d0) then
         y = x*x
         rtemp1 = 57568490574.0d0 + y*(-13362690354.0d0 + y*(651619640.7d0 + y*(-11214424.18d0 + y*(77392.33017d0 + y*(-184.9052456d0)))))
         rtemp2 = 57568490411.0d0 + y*(1029532985.0d0 + y*(9494680.718d0 + y*(59272.64853d0 + y*(267.8532712d0 + y*1.0d0))))
         BesselJ0_func = rtemp1/rtemp2
      else
         z = 8.0d0/ax
         y = z*z

         xx = ax - 0.785398164d0
         rtemp1 = 1.0d0 + y*(-0.1098628627d-2 + y*(0.2734510407d-4 + y*(-0.2073370639d-5 + y*0.2093887211d-6)))
         rtemp2 = -0.1562499995d-1 + y*(0.1430488765d-3 + y*(-0.6911147651d-5 + y*(0.7621095161d-6 - y*0.934935152d-7)))
         BesselJ0_func = sqrt(0.636619772d0/ax)*(cos(xx)*rtemp1 - z*sin(xx)*rtemp2)
      endif

      return

   end function BesselJ0_func

   real(kind=8) function BesselY0_func(x)

      implicit none

      real(kind=8) x, z, ax
      real(kind=8) y, rtemp1, rtemp2, xx

      if (x < 8.0d0) then
         y = x*x
         rtemp1 = -2957821389.0d0 + y*(7062834065.0d0 + y*(-512359803.6d0 + y*(10879881.29d0 + y*(-86327.92757d0 + y*228.4622733d0))))
         rtemp2 = 40076544269.0d0 + y*(745249964.8d0 + y*(7189466.438d0 + y*(47447.26470d0 + y*(226.1030244d0 + y*1.0d0))))
         BesselY0_func = (rtemp1/rtemp2) + 0.636619772d0*BesselJ0_func(x)*LOG(x)
      else
         z = 8.0d0/x
         y = z*z

         xx = x - 0.785398164d0
         rtemp1 = 1.0d0 + y*(-0.1098628627d-2 + y*(0.2734510407d-4 + y*(-0.2073370639d-5 + y*0.2093887211d-6)))
         rtemp2 = -0.1562499995d-1 + y*(0.1430488765d-3 + y*(-0.6911147651d-5 + y*(0.7621095161d-6 - y*0.934935152d-7)))
         BesselY0_func = sqrt(0.636619772d0/x)*(sin(xx)*rtemp1 + z*cos(xx)*rtemp2)
      endif

      return

   end function BesselY0_func






   ! **************************************************************************************** !
   ! *****             bessjyV                                                          ***** !
   ! *****             Nmax: the maximum order                                          ***** !
   ! *****             x: the input x for j_n                                           ***** !
   ! *****             sjv: bessel function                                   ***** !
   ! *****             bf:upward,downward flag (upward is not stable)                   ***** !
   ! **************************************************************************************** !
   ! calculate first kind bessel function j_n(x) with n=0,1,...,Nmax iteratively
   ! Be careful, this routine is not accurate for large argument (x)
   SUBROUTINE bessjyV(Nmax,x,sjv,bf)
      implicit none
      integer Nmax,ii
      real*8 x
      real*8,dimension(Nmax+1)::sjv
      integer bf,kk
      integer,parameter:: kkmax = 20
      real*8:: b0,b1,b_1,b_2,a0,a_1,b_Nmax,b_Nmax_1,tmp1,tmp2,fk_Nmax,fk_Nmax_1,fk_Nmax_array(0:kkmax),fk_Nmax_1_array(0:kkmax)
      real*8:: b_n(Nmax+1)
      integer non_O_idx
      integer underflow
      ! numerical recipe
      if(abs(x)<1d-16)then
         sjv(:)=0d0
         sjv(1)=1d0
      else
         if(bf==1)then
            call bessjy(x,0d0,sjv(1))
            if(Nmax==0)return
            call bessjy(x,1d0,sjv(2))
            do ii=2,Nmax
               sjv(ii+1) = (2*ii-2)*sjv(ii)/x - sjv(ii-1)
            end do
         else
            call bessjy(x,dble(Nmax),sjv(Nmax+1))
            if(Nmax==0)return
            call bessjy(x,dble(Nmax-1),sjv(Nmax))
            do ii=Nmax-2,0,-1
               sjv(ii+1) = (2*ii+2)*sjv(ii+2)/x - sjv(ii+3)
            end do
         end if
      end if

      underflow = 0
      tmp1 = sum(sjv)
      if(tmp1 .NE. tmp1)underflow = 1
      if(sum(abs(sjv))==0d0)then
         underflow = 1
      else
         do ii = Nmax+1,1,-1
            if(sjv(ii)/=0)then
               non_O_idx = ii
               exit
            end if
         end do
         if(abs(sjv(non_O_idx))<1d-300)then
            underflow = 1
         end if
      end if

      if(underflow == 1)then
         write(*,*)"warning: underflow in bessjyV"
      end if

   end SUBROUTINE	bessjyV


   ! **************************************************************************************** !
   ! *****       bessjy: calculate J_v(x), v real, x real                              ****** !
   ! **************************************************************************************** !
      SUBROUTINE bessjy(x,xnu,rj)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: x,xnu
      REAL*8, INTENT(OUT) :: rj
      INTEGER*8, PARAMETER :: MAXIT=1000000 !originally MAXIT = 1000, the bigger MAXIT, the larger supported x
      REAL*8, PARAMETER :: XMIN=2.0d0,EPS=1.0d-10,FPMIN=1.0d-30
      REAL*8, PARAMETER :: PI_D=3.141592653589793238462643383279502884197d0
      INTEGER*4 :: i,isign,l,nl
      REAL*8 :: ry,rjp,ryp
      REAL*8 :: a,b,c,d,del,del1,e,f,fact,fact2,fact3,ff,gam,gam1,gam2,&
         gammi,gampl,h,p,pimu,pimu2,q,r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp,&
         ry1,rymu,rymup,rytemp,sum,sum1,w,x2,xi,xi2,xmu,xmu2
      COMPLEX(kind=8) :: aa,bb,cc,dd,dl,pq
      nl=merge(int(xnu+0.5d0), max(0,int(xnu-x+1.5d0)), x < XMIN)
      xmu=xnu-nl
      xmu2=xmu*xmu
      if(x <=0)stop
      xi=1.0d0/x
      xi2=2.0d0*xi
      w=xi2/PI_D
      isign=1
      h=xnu*xi
      if (h < FPMIN) h=FPMIN
      b=xi2*xnu
      d=0.0
      c=h
      do i=1,MAXIT
         b=b+xi2
         d=b-d
         if (abs(d) < FPMIN) d=FPMIN
         c=b-1.0d0/c
         if (abs(c) < FPMIN) c=FPMIN
         d=1.0d0/d
         del=c*d
         h=del*h
         if (d < 0.0) isign=-isign
         if (abs(del-1.0d0) < EPS) exit
      end do
      call assert(i<=MAXIT,'x too large; try asymptotic expansion')

      rjl=isign*FPMIN
      rjpl=h*rjl
      rjl1=rjl
      rjp1=rjpl
      fact=xnu*xi
      do l=nl,1,-1
         rjtemp=fact*rjl+rjpl
         fact=fact-xi
         rjpl=fact*rjtemp-rjl
         rjl=rjtemp
      end do
      if (rjl == 0.0) rjl=EPS
      f=rjpl/rjl
      if (x < XMIN) then
         x2=0.5d0*x
         pimu=PI_D*xmu
         if (abs(pimu) < EPS) then
            fact=1.0
         else
            fact=pimu/sin(pimu)
         end if
         d=-log(x2)
         e=xmu*d
         if (abs(e) < EPS) then
            fact2=1.0
         else
            fact2=sinh(e)/e
         end if
         call beschb(xmu,gam1,gam2,gampl,gammi)
         ff=2.0d0/PI_D*fact*(gam1*cosh(e)+gam2*fact2*d)
         e=exp(e)
         p=e/(gampl*PI_D)
         q=1.0d0/(e*PI_D*gammi)
         pimu2=0.5d0*pimu
         if (abs(pimu2) < EPS) then
            fact3=1.0
         else
            fact3=sin(pimu2)/pimu2
         end if
         r=PI_D*pimu2*fact3*fact3
         c=1.0
         d=-x2*x2
         sum=ff+r*q
         sum1=p
         do i=1,MAXIT
            ff=(i*ff+p+q)/(i*i-xmu2)
            c=c*d/i
            p=p/(i-xmu)
            q=q/(i+xmu)
            del=c*(ff+r*q)
            sum=sum+del
            del1=c*p-i*del
            sum1=sum1+del1
            if (abs(del) < (1.0d0+abs(sum))*EPS) exit
         end do
         call assert(i<=MAXIT,'bessy series failed to converge')
         rymu=-sum
         ry1=-sum1*xi2
         rymup=xmu*xi*rymu-ry1
         rjmu=w/(rymup-f*rymu)
      else
         a=0.25d0-xmu2
         pq=cmplx(-0.5d0*xi,1.0d0,kind=8)
         aa=cmplx(0.0d0,xi*a,kind=8)
         bb=cmplx(2.0d0*x,2.0d0,kind=8)
         cc=bb+aa/pq
         dd=1.0d0/bb
         pq=cc*dd*pq
         do i=2,MAXIT
            a=a+2*(i-1)
            bb=bb+cmplx(0.0d0,2.0d0,kind=8)
            dd=a*dd+bb
            if (absc(dd) < FPMIN) dd=FPMIN
            cc=bb+a/cc
            if (absc(cc) < FPMIN) cc=FPMIN
            dd=1.0d0/dd
            dl=cc*dd
            pq=pq*dl
            if (absc(dl-1.0d0) < EPS) exit
         end do
         call assert(i<=MAXIT,'cf2 failed in bessjy')
         p=dble(pq)
         q=DIMAG(pq)
         gam=(p-f)/q
         rjmu=sqrt(w/((p-f)*gam+q))
         rjmu=sign(rjmu,rjl)
         rymu=rjmu*gam
         rymup=rymu*(p+q/gam)
         ry1=xmu*xi*rymu-rymup
      end if
      fact=rjmu/rjl
      rj=rjl1*fact
      rjp=rjp1*fact
      do i=1,nl
         rytemp=(xmu+i)*xi2*ry1-rymu
         rymu=ry1
         ry1=rytemp
      end do
      ry=rymu
      ryp=xnu*xi*rymu-ry1
      CONTAINS
   !BL
      FUNCTION absc(z)
      IMPLICIT NONE
      COMPLEX(kind = 8), INTENT(IN) :: z
      REAL*8 :: absc
      absc=abs(dble(z))+abs(DIMAG(z))
      END FUNCTION absc

      END SUBROUTINE bessjy


   ! **************************************************************************************** !
   ! *****       beschb: ......                                                        ****** !
   ! **************************************************************************************** !

      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: x
      REAL*8, INTENT(OUT) :: gam1,gam2,gampl,gammi
      INTEGER*4, PARAMETER :: NUSE1=5,NUSE2=5
      REAL*8 :: xx
      REAL*8, DIMENSION(7) :: c1=(/-1.142022680371168d0,&
         6.5165112670737d-3,3.087090173086d-4,-3.4706269649d-6,&
         6.9437664d-9,3.67795d-11,-1.356d-13/)
      REAL*8, DIMENSION(8) :: c2=(/1.843740587300905d0,&
         -7.68528408447867d-2,1.2719271366546d-3,&
         -4.9717367042d-6, -3.31261198d-8,2.423096d-10,&
         -1.702d-13,-1.49d-15/)
      xx=8.0*x*x-1.0
      gam1=chebev(-1.0d0,1.0d0,c1(1:NUSE1),xx)
      gam2=chebev(-1.0d0,1.0d0,c2(1:NUSE2),xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      END SUBROUTINE beschb




   ! **************************************************************************************** !
   ! *****       chebev: ......                                                        ****** !
   ! **************************************************************************************** !
      FUNCTION chebev(a,b,c,x)
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: a,b,x
      REAL*8, DIMENSION(:), INTENT(IN) :: c
      REAL*8 :: chebev
      INTEGER*4 :: j,m
      REAL*8 :: d,dd,sv,y,y2
      call assert((x-a)*(x-b) <= 0.0,'x not in range in chebev_s')
      m=size(c)
      d=0.0
      dd=0.0
      y=(2.0*x-a-b)/(b-a)
      y2=2.0*y
      do j=m,2,-1
         sv=d
         d=y2*d-dd+c(j)
         dd=sv
      end do
      chebev=y*d-dd+0.5*c(1)
      END FUNCTION chebev




!>**** create a array treeleaf holding size of each leaf box of a tree with nlevel levels (0<=level<=nlevel)
   recursive subroutine CreateLeaf_Natural(nlevel, level, group, idxs, idxe, treeleaf)
      implicit none
      integer nlevel, level, group, idxs, idxe
      integer treeleaf(2**nlevel)

      if (level == nlevel) then
         treeleaf(group - 2**nlevel + 1) = idxe - idxs + 1
      else
         call CreateLeaf_Natural(nlevel, level + 1, 2*group, idxs, int((idxs + idxe)/2d0), treeleaf)
         call CreateLeaf_Natural(nlevel, level + 1, 2*group + 1, int((idxs + idxe)/2d0) + 1, idxe, treeleaf)
      endif
   end subroutine CreateLeaf_Natural

   subroutine NumberingPtree(ptree)
      implicit none
      type(proctree)::ptree
      integer :: level, group

      ptree%pgrp(1)%head = 0; ptree%pgrp(1)%tail = ptree%nproc - 1; ptree%pgrp(1)%nproc = ptree%nproc
      do level = 0, ptree%nlevel - 1
         do group = 2**level, 2**(level + 1) - 1
            if (level < ptree%nlevel - 1) then
               if (ptree%pgrp(group)%nproc == 1) then
                  ptree%pgrp(2*group)%head = ptree%pgrp(group)%head
                  ptree%pgrp(2*group)%tail = ptree%pgrp(group)%tail
                  ptree%pgrp(2*group)%nproc = ptree%pgrp(group)%nproc
                  ptree%pgrp(2*group + 1)%head = ptree%pgrp(group)%head
                  ptree%pgrp(2*group + 1)%tail = ptree%pgrp(group)%tail
                  ptree%pgrp(2*group + 1)%nproc = ptree%pgrp(group)%nproc
               else
                  ptree%pgrp(2*group)%head = ptree%pgrp(group)%head
                  ptree%pgrp(2*group)%tail = int((ptree%pgrp(group)%head + ptree%pgrp(group)%tail)/2)
                  ptree%pgrp(2*group)%nproc = ptree%pgrp(2*group)%tail - ptree%pgrp(2*group)%head + 1

                  ptree%pgrp(2*group + 1)%head = ptree%pgrp(2*group)%tail + 1
                  ptree%pgrp(2*group + 1)%tail = ptree%pgrp(group)%tail
                  ptree%pgrp(2*group + 1)%nproc = ptree%pgrp(2*group + 1)%tail - ptree%pgrp(2*group + 1)%head + 1
               endif
            end if
         end do
      end do
   end subroutine NumberingPtree

!>**** computation of the local butterfly block ranges owned by this MPI rank
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
   subroutine GetLocalBlockRange(ptree, pgno, level, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, dir)
      implicit none
      type(proctree)::ptree
      integer :: level_p, level, ll, group, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, pgno, pgno1, found, nleaf, ith
      character::dir
      if (IOwnPgrp(ptree, pgno)) then
         found = 0
         level_p = 0
         pgno1 = pgno
         ith = 1
         do while (found == 0)
            if (ptree%pgrp(pgno1)%nproc == 1 .or. level_p == level_butterfly) then
               found = 1
               exit
            endif
            if (IOwnPgrp(ptree, pgno1*2)) then
               pgno1 = pgno1*2
               ith = ith*2
            elseif (IOwnPgrp(ptree, pgno1*2 + 1)) then
               pgno1 = pgno1*2 + 1
               ith = ith*2 + 1
            endif
            level_p = level_p + 1
         enddo
         ith = ith - 2**level_p + 1
         call assert(level_p <= level_butterfly, 'too many processes sharing this block')
         nleaf = 2**(level_butterfly - level_p)
         if (dir == 'R') then
            idx_c = (ith - 1)*nleaf + 1
            nc = nleaf
            inc_c = 1
            idx_r = 1
            nr = 1
            inc_r = 1
            do ll = 1, level
               if (ll /= level_butterfly + 1) then
                  idx_r = 2*idx_r - mod(idx_c, 2) ! odd/even column indices mapped to odd/even row indices
                  idx_c = floor((idx_c - 1)/2d0) + 1
                  if (nc > 1) then ! if at previous level column blocks are contiguous, double #row and half #column
                     nc = nc/2
                     nr = nr*2
                  else         ! otherwise double row increments
                     inc_r = inc_r*2
                  endif
               endif
            enddo
         elseif (dir == 'C') then
            idx_r = (ith - 1)*nleaf + 1
            nr = nleaf
            inc_r = 1
            idx_c = 1
            nc = 1
            inc_c = 1
            do ll = level_butterfly, level, -1
               if (ll /= 0) then
                  idx_c = 2*idx_c - mod(idx_r, 2) ! odd/even row indices mapped to odd/even column indices
                  idx_r = floor((idx_r - 1)/2d0) + 1
                  if (nr > 1) then ! if at previous level row blocks are contiguous, double #column and half #row
                     nr = nr/2
                     nc = nc*2
                  else
                     inc_c = inc_c*2 ! otherwise double column increments
                  endif
               endif
            enddo
         endif

         if (ptree%MyID /= ptree%pgrp(pgno1)%head .and. level > 0 .and. level < level_butterfly + 1) then ! for the kernel levels: if several processe own one block, only the head process has it.
            idx_r = 0
            inc_r = 0
            nr = 0
            idx_c = 0
            inc_c = 0
            nc = 0
         endif
      else
         idx_r = 0
         inc_r = 0
         nr = 0
         idx_c = 0
         inc_c = 0
         nc = 0
      endif
   end subroutine GetLocalBlockRange




   !>****  reallocate a(:) to a new size while preserving old data. If new_size<old_size, the old data is truncated.
   subroutine array_resize(a, new_size)
      implicit none
      DT, allocatable, intent(inout) :: a(:)
      integer      :: new_size

      integer :: old_size
      DT, allocatable :: tmp(:)

      call assert(new_size>=0,"array new size incorrect")

      if(.not. allocated(a))then
         if(new_size>0)allocate(a(new_size))
      else
         old_size = size(a)
         allocate(tmp(old_size))
         tmp = a
         deallocate(a)

         allocate(a(new_size))
         a=0d0
         a(1:min(old_size,new_size)) = tmp(1:min(old_size,new_size))
         deallocate(tmp)
      endif
   end subroutine array_resize


   !>**** convert single index to multi-index, assuming first index is the fastest
   subroutine SingleIndexToMultiIndex(Ndim,dims, single_index_in, multi_index)
      implicit none
      integer:: Ndim
      integer:: dims(Ndim)
      integer :: single_index, single_index_in, multi_index(Ndim)
      integer:: i,product,dim_i
      integer*8:: nelem

      nelem=1
      do dim_i=1,Ndim
         nelem = nelem * dims(dim_i)
      enddo
      call assert(nelem<=2.14D9,'integer overflow in SingleIndexToMultiIndex')

      single_index = single_index_in

      ! Initialize multi_index
      multi_index = 0
      product=1

      ! Calculate Multi-Indices
      do i = 1,Ndim
         multi_index(i) = mod((single_index - 1)/product, dims(i)) + 1
         single_index = single_index - (multi_index(i) - 1) *product
         product = product * dims(i)
      end do

   end subroutine SingleIndexToMultiIndex


   !>**** convert multi-index to single-index, assuming first index is the fastest
   subroutine MultiIndexToSingleIndex(Ndim,dims, single_index, multi_index)
      implicit none
      integer:: Ndim
      integer:: dims(Ndim)
      integer :: single_index, multi_index(Ndim)
      integer :: i, product, dim_i
      integer*8 :: nelem

      nelem=1
      do dim_i=1,Ndim
         nelem = nelem * dims(dim_i)
      enddo
      call assert(nelem<=2.14D9,'integer overflow in MultiIndexToSingleIndex')

      ! Initialize single_index
      single_index = 1

      ! Calculate Single Index
      product = 1
      do i = 1, Ndim
         single_index = single_index + (multi_index(i) - 1) * product
         product = product * dims(i)
      end do

   end subroutine MultiIndexToSingleIndex







!>**** computation of the local (multi-dimensinonal) butterfly block ranges owned by this MPI rank
   !ptree: process tree
   !pgno: the process group number that shares this butterfly
   !level_butterfly: number of butterfly levels
   !level: the specified level that requires the computation 0<=level<=level_butterfly+1
   !ndim: number of dimensions
   !dim_s: which dimension to split first at the level 1 of the process tree
   !idx_r(ndim): starting index of local row blocks
   !inc_r(ndim): increments of local row blocks
   !nr(ndim): number of local row blocks
   !idx_c(ndim): starting index of local column blocks
   !inc_c(ndim): increments of local column blocks
   !nc(ndim): number of local column blocks
   !dir: 'R': row-wise ordering (2^level row blocks and 2^(level_butterfly-level) column blocks) or 'C': column-wise ordering (2^(level-1) row blocks and 2^(level_butterfly-level+1) column blocks)
   subroutine GetLocalBlockRange_MD(ptree, pgno, level, level_butterfly, ndim, dim_s, idx_r, inc_r, nr, idx_c, inc_c, nc, dir)
      implicit none
      type(proctree)::ptree
      integer :: ndim, dim_s, dim_i
      integer :: level_p(ndim), level, ll, group, level_butterfly, level_half, idx_r(ndim), inc_r(ndim), nr(ndim), idx_c(ndim), inc_c(ndim), nc(ndim), pgno, pgno1, found, nleaf(ndim), ith(ndim)
      character::dir

      level_half = floor_safe(dble(level_butterfly)/2d0) ! from outer to inner

      if (IOwnPgrp(ptree, pgno)) then
         found = 0
         level_p = 0
         pgno1 = pgno
         ith = 1
         dim_i = dim_s
         do while (found == 0)
            if (ptree%pgrp(pgno1)%nproc == 1 .or. minval(level_p) == level_half) then
               found = 1
               exit
            endif
            if (IOwnPgrp(ptree, pgno1*2)) then
               pgno1 = pgno1*2
               ith(dim_i) = ith(dim_i)*2
            elseif (IOwnPgrp(ptree, pgno1*2 + 1)) then
               pgno1 = pgno1*2 + 1
               ith(dim_i) = ith(dim_i)*2 + 1
            endif
            level_p(dim_i) = level_p(dim_i) + 1
            dim_i = dim_i + 1
            dim_i = mod(dim_i-1,ndim)+1 ! reset dim to 1 if dim=ndim+1
         enddo
         ith = ith - 2**level_p + 1
         call assert(minval(level_p) <= level_half, 'too many processes sharing this block as GetLocalBlockRange_MD parallelization goes at most to the middle level')
         nleaf = 2**(level_butterfly - level_p)

         do dim_i=1,ndim
            if (dir == 'R') then
               idx_c(dim_i) = (ith(dim_i) - 1)*nleaf(dim_i) + 1
               nc(dim_i) = nleaf(dim_i)
               inc_c(dim_i) = 1
               idx_r(dim_i) = 1
               nr(dim_i) = 1
               inc_r(dim_i) = 1
               do ll = 1, level
                  if (ll /= level_butterfly + 1) then
                     idx_r(dim_i) = 2*idx_r(dim_i) - mod(idx_c(dim_i), 2) ! odd/even column indices mapped to odd/even row indices
                     idx_c(dim_i) = floor((idx_c(dim_i) - 1)/2d0) + 1
                     if (nc(dim_i) > 1) then ! if at previous level column blocks are contiguous, double #row and half #column
                        nc(dim_i) = nc(dim_i)/2
                        nr(dim_i) = nr(dim_i)*2
                     else         ! otherwise double row increments
                        write(*,*)'nc(dim_i) == 1 should not happen for GetLocalBlockRange_MD'
                        stop
                     endif
                  endif
               enddo
            elseif (dir == 'C') then
               idx_r(dim_i) = (ith(dim_i) - 1)*nleaf(dim_i) + 1
               nr(dim_i) = nleaf(dim_i)
               inc_r(dim_i) = 1
               idx_c(dim_i) = 1
               nc(dim_i) = 1
               inc_c(dim_i) = 1
               do ll = level_butterfly, level, -1
                  if (ll /= 0) then
                     idx_c(dim_i) = 2*idx_c(dim_i) - mod(idx_r(dim_i), 2) ! odd/even row indices mapped to odd/even column indices
                     idx_r(dim_i) = floor((idx_r(dim_i) - 1)/2d0) + 1
                     if (nr(dim_i) > 1) then ! if at previous level row blocks are contiguous, double #column and half #row
                        nr(dim_i) = nr(dim_i)/2
                        nc(dim_i) = nc(dim_i)*2
                     else
                        write(*,*)'nr(dim_i) == 1 should not happen for GetLocalBlockRange_MD'
                        stop
                     endif
                  endif
               enddo
            endif

            if (ptree%MyID /= ptree%pgrp(pgno1)%head .and. level > 0 .and. level < level_butterfly + 1) then ! for the kernel levels: if several processe own one block, only the head process has it.
               idx_r(dim_i) = 0
               inc_r(dim_i) = 0
               nr(dim_i) = 0
               idx_c(dim_i) = 0
               inc_c(dim_i) = 0
               nc(dim_i) = 0
               write(*,*)'pgno1 should only have 1 rank in GetLocalBlockRange_MD'
               stop
            endif
         enddo

      else
         idx_r = 0
         inc_r = 0
         nr = 0
         idx_c = 0
         inc_c = 0
         nc = 0
      endif
   end subroutine GetLocalBlockRange_MD



subroutine TensorUnfoldingReshape(Ndim,dims_ref_old,dims_ref_new,offsets_ref,ld_old,ld_new,data_in,trans_in,data_out,trans_out,loopnew)
   implicit none
   integer Ndim,loopnew
   integer dims_ref_old(Ndim),dims_ref_new(Ndim),offsets_ref(Ndim), idx_ref(Ndim), dims_new(Ndim-1),dims_new_scalar, idx_new(Ndim-1),idx_new_scalar,dims_old(Ndim-1),dims_old_scalar, idx_old(Ndim-1),idx_old_scalar
   integer ld_old,ld_new, dim_i, i, j, j1
   character::trans_in,trans_out
   DT::data_in(:,:)
   DT,allocatable::data_out(:,:)

   i=0
   do dim_i=1,Ndim
      if(dim_i/=ld_old)then
         i = i+ 1
         dims_old(i) = dims_ref_old(dim_i)
      else
         dims_old_scalar = dims_ref_old(dim_i)
      endif
   enddo

   i=0
   do dim_i=1,Ndim
      if(dim_i/=ld_new)then
         i = i+ 1
         dims_new(i) = dims_ref_new(dim_i)
      else
         dims_new_scalar = dims_ref_new(dim_i)
      endif
   enddo
   if(.not. allocated(data_out))then
      if(trans_out=='N')then
         allocate(data_out(dims_new_scalar,product(dims_new)))
      else
         allocate(data_out(product(dims_new),dims_new_scalar))
      endif
      data_out=0
   endif

   if(loopnew==1)then
      do idx_new_scalar=1,dims_new_scalar
      do j=1,product(dims_new)
         call SingleIndexToMultiIndex(Ndim-1, dims_new, j, idx_new)
         !! get the Ndim dimensional index idx_ref
         i=0
         do dim_i=1,Ndim
         if(dim_i/=ld_new)then
            i=i+1
            idx_ref(dim_i) = idx_new(i) + offsets_ref(dim_i)
         else
            idx_ref(dim_i) = idx_new_scalar + offsets_ref(dim_i)
         endif
         enddo
         !! get the Ndim-1 dimensional index idx_old
         i=0
         do dim_i=1,Ndim
         if(dim_i/=ld_old)then
            i=i+1
            idx_old(i)=idx_ref(dim_i)
         else
            idx_old_scalar=idx_ref(dim_i)
         endif
         enddo

         call MultiIndexToSingleIndex(Ndim-1, dims_old, j1, idx_old)

         if(trans_out=='N' .and. trans_in=='N')then
            data_out(idx_new_scalar,j) = data_in(idx_old_scalar,j1)
         else if(trans_out=='N' .and. trans_in=='T')then
            data_out(idx_new_scalar,j) = data_in(j1,idx_old_scalar)
         else if(trans_out=='T' .and. trans_in=='N')then
            data_out(j,idx_new_scalar) = data_in(idx_old_scalar,j1)
         else if(trans_out=='T' .and. trans_in=='T')then
            data_out(j,idx_new_scalar) = data_in(j1,idx_old_scalar)
         endif
      enddo
      enddo

   else
      do idx_old_scalar=1,dims_old_scalar
      do j=1,product(dims_old)
         call SingleIndexToMultiIndex(Ndim-1, dims_old, j, idx_old)
         !! get the Ndim dimensional index idx_ref
         i=0
         do dim_i=1,Ndim
         if(dim_i/=ld_old)then
            i=i+1
            idx_ref(dim_i) = idx_old(i) + offsets_ref(dim_i)
         else
            idx_ref(dim_i) = idx_old_scalar + offsets_ref(dim_i)
         endif
         enddo
         !! get the Ndim-1 dimensional index idx_new
         i=0
         do dim_i=1,Ndim
         if(dim_i/=ld_new)then
            i=i+1
            idx_new(i)=idx_ref(dim_i)
         else
            idx_new_scalar=idx_ref(dim_i)
         endif
         enddo

         call MultiIndexToSingleIndex(Ndim-1, dims_new, j1, idx_new)

         if(trans_out=='N' .and. trans_in=='N')then
            data_out(idx_new_scalar,j1) = data_in(idx_old_scalar,j)
         else if(trans_out=='N' .and. trans_in=='T')then
            data_out(idx_new_scalar,j1) = data_in(j,idx_old_scalar)
         else if(trans_out=='T' .and. trans_in=='N')then
            data_out(j1,idx_new_scalar) = data_in(idx_old_scalar,j)
         else if(trans_out=='T' .and. trans_in=='T')then
            data_out(j1,idx_new_scalar) = data_in(j,idx_old_scalar)
         endif
      enddo
      enddo

   endif


end subroutine TensorUnfoldingReshape




subroutine TensorPermute(tensor_data, shape_in, perm)
   implicit none
   DT :: tensor_data(:)
   integer  :: shape_in(:)
   integer   :: perm(:)
   DT, allocatable :: tmp_data(:)
   integer,allocatable :: shape_out(:)

   integer :: ndims, totalSize
   integer :: lin,lin_new
   integer, allocatable :: idxs(:), new_idxs(:)
   ndims = size(shape_in)
   allocate(shape_out(ndims))
   shape_out = shape_in(perm)

   totalSize = size(tensor_data)
   allocate(tmp_data(totalSize))

   allocate(idxs(ndims), new_idxs(ndims))

   tmp_data = 0d0

   do lin=1,totalSize
      call SingleIndexToMultiIndex(ndims,shape_in, lin, idxs)
      ! dimension i of new index = dimension perm(i) of old
      new_idxs = idxs(perm)
      call MultiIndexToSingleIndex(ndims,shape_out, lin_new, new_idxs)
      tmp_data(lin_new ) = tensor_data(lin)
   end do

   tensor_data = tmp_data
   deallocate(tmp_data)
   deallocate(shape_out)
   deallocate(idxs, new_idxs)
end subroutine TensorPermute





!>**** computation of the sub process group that handles one block of the outtermost factor of a butterfly. Note: if level_butterfly=0, then pgno_sub=pgno
   !ptree: process tree
   !pgno: the process group number that shares this butterfly
   !level_butterfly: number of butterfly levels
   !pgno_sub: the sub process group number that shares one block of this butterfly
   subroutine GetPgno_Sub(ptree, pgno, level_butterfly, pgno_sub)
      implicit none
      type(proctree)::ptree
      integer :: level_p, level, ll, group, level_butterfly, idx_r, inc_r, nr, idx_c, inc_c, nc, pgno, pgno_sub, found, nleaf, ith
      character::dir
      pgno_sub = -1
      if (IOwnPgrp(ptree, pgno)) then
         found = 0
         level_p = 0
         pgno_sub = pgno
         do while (found == 0)
            if (ptree%pgrp(pgno_sub)%nproc == 1 .or. level_p == level_butterfly) then
               found = 1
               exit
            endif
            if (IOwnPgrp(ptree, pgno_sub*2)) then
               pgno_sub = pgno_sub*2
            elseif (IOwnPgrp(ptree, pgno_sub*2 + 1)) then
               pgno_sub = pgno_sub*2 + 1
            endif
            level_p = level_p + 1
         enddo
      endif
   end subroutine GetPgno_Sub

!>**** computation of the process group number "pgno_sub" that shares the (index_i,index_j,level) block. Note for blocks in the kernels, only the head process in pgno_sub is active; for blocks in the outtermost factors, all processes could be active
   !ptree: process tree
   !pgno: the process group number that shares this butterfly
   !level_butterfly: number of butterfly levels
   !level: the level where the block resides 0<=level<=level_butterfly+1
   !index_i: row index of the block
   !index_j: column index of the block
   !dir: 'R': row-wise ordering (2^level row blocks and 2^(level_butterfly-level) column blocks) or 'C': column-wise ordering (2^(level-1) row blocks and 2^(level_butterfly-level+1) column blocks)
   !pgno_sub: the process group number that shares the (index_i,index_j) block
   subroutine GetBlockPID(ptree, pgno, level, level_butterfly, index_i, index_j, dir, pgno_sub)
      implicit none
      type(proctree)::ptree
      integer :: pid, level_p, num_blocks, level, ll, group, level_butterfly, idx_r, idx_c, idx, index_i, index_j, pgno, pgno_sub, found, nleaf, ith
      character::dir

      level_p = min(level_butterfly, ptree%nlevel - GetTreelevel(pgno))
      idx_r = index_i
      idx_c = index_j

      if (dir == 'R') then
         do ll = level, 1, -1
            if (ll /= level_butterfly + 1) then
               idx_c = 2*idx_c - mod(idx_r, 2) ! odd/even row indices mapped to odd/even column indices
               idx_r = floor_safe((idx_r - 1)/2d0) + 1
            endif
         enddo
         idx = idx_c
      elseif (dir == 'C') then
         do ll = level, level_butterfly
            if (ll /= 0) then
               idx_r = 2*idx_r - mod(idx_c, 2) ! odd/even column indices mapped to odd/even row indices
               idx_c = floor_safe((idx_c - 1)/2d0) + 1
            endif
         enddo
         idx = idx_r
      endif
      idx = idx + 2**level_butterfly - 1 ! index of this leaf-level group in the corresponding level_butterfly-level tree
      do ll = 1, level_butterfly - level_p
         idx = floor_safe(idx/2d0)
      enddo
      idx = idx - 2**level_p + 1 ! index of this group's ancestor group in the corresponding level_p-level tree
      pgno_sub = pgno*2**level_p ! index of the first process group at pgno's level + level_p of the process tree
      pgno_sub = pgno_sub + idx - 1 ! index of the process group that handles group idx in the process tree

      ! pid = ptree%pgrp(pgno_sub)%head
   end subroutine GetBlockPID


!>**** computation of the process group number "pgno_sub" that shares the multi-dimensional (index_i(ndim),index_j(ndim),level) block. Note for blocks in the kernels, only the head process in pgno_sub is active; for blocks in the outtermost factors, all processes could be active
   !ptree: process tree
   !pgno: the process group number that shares this multi-dimensional butterfly
   !ndim: number of dimensions
   !ndim: number of dimensions
   !dim_s: which dimension to split first at the level 1 of the process tree
   !level_butterfly: number of butterfly levels
   !level: the level where the block resides 0<=level<=level_butterfly+1
   !index_i(ndim): row index of the block
   !index_j(ndim): column index of the block
   !dir: 'R': row-wise ordering (2^level row blocks and 2^(level_butterfly-level) column blocks) or 'C': column-wise ordering (2^(level-1) row blocks and 2^(level_butterfly-level+1) column blocks)
   !pgno_sub: the process group number that shares the (index_i,index_j) block
   subroutine GetBlockPID_MD(ptree, pgno, ndim, dim_s, level, level_butterfly, index_i, index_j, dir, pgno_sub)
      implicit none
      type(proctree)::ptree
      integer:: ndim,dim_i,dim_s, idx_tmp
      integer:: idx_r(ndim), idx_c(ndim), idx(ndim), index_i(ndim), index_j(ndim), level_ps(ndim), idx_p(ndim)
      integer :: pid, pp, level_p, num_blocks, level, ll, group, level_butterfly, pgno, pgno_sub, found
      character::dir

      idx_r = index_i
      idx_c = index_j

      do dim_i=1,ndim
         if (dir == 'R') then
            do ll = level, 1, -1
               if (ll /= level_butterfly + 1) then
                  idx_c(dim_i) = 2*idx_c(dim_i) - mod(idx_r(dim_i), 2) ! odd/even row indices mapped to odd/even column indices
                  idx_r(dim_i) = floor_safe((idx_r(dim_i) - 1)/2d0) + 1
               endif
            enddo
            idx(dim_i) = idx_c(dim_i)
         elseif (dir == 'C') then
            do ll = level, level_butterfly
               if (ll /= 0) then
                  idx_r(dim_i) = 2*idx_r(dim_i) - mod(idx_c(dim_i), 2) ! odd/even column indices mapped to odd/even row indices
                  idx_c(dim_i) = floor_safe((idx_c(dim_i) - 1)/2d0) + 1
               endif
            enddo
            idx(dim_i) = idx_r(dim_i)
         endif
         idx(dim_i) = idx(dim_i) + 2**level_butterfly - 1 ! index of this leaf-level group in the corresponding level_butterfly-level tree
      enddo

      level_p = ptree%nlevel - GetTreelevel(pgno)
      pgno_sub = pgno
      dim_i = dim_s
      level_ps=0
      idx_p = 1
      do pp=1,level_p
         idx_tmp = idx(dim_i)
         do ll = 1, level_butterfly - level_ps(dim_i)-1
            idx_tmp = floor_safe(idx_tmp/2d0)
         enddo
         if(idx_tmp==idx_p(dim_i)*2)then
            pgno_sub=pgno_sub*2
         elseif(idx_tmp==idx_p(dim_i)*2+1)then
            pgno_sub=pgno_sub*2+1
         else
            write(*,*)"something wrong for the value of idx_tmp and idx_p", idx_tmp, idx_p(dim_i)
            stop
         endif
         idx_p(dim_i) = idx_tmp
         level_ps(dim_i) = level_ps(dim_i) + 1
         dim_i = dim_i + 1
         dim_i = mod(dim_i-1,ndim)+1 ! reset dim to 1 if dim=ndim+1
      enddo

   end subroutine GetBlockPID_MD






!>**** get the pgno for the multi-dimensinoal group in msh.
   !ptree: process tree
   !ndim: number of dimensions
   !group(ndim): multi-index for the group
   integer function GetMshGroup_Pgno(ptree, ndim, group)
      implicit none
      type(proctree)::ptree
      integer:: ndim,dim_i, idx_tmp, Maxgrp
      integer:: idx_p(ndim)
      integer :: pp, level, ll, group(*), pgno


      level = GetTreelevel(group(1)) - 1
      pgno=1
      Maxgrp = 2**(ptree%nlevel) - 1
      idx_p = 1
      do pp=1,level
         do dim_i=1,Ndim
            idx_tmp = group(dim_i)
            do ll = 1, level - pp
               idx_tmp = floor_safe(idx_tmp/2d0)
            enddo
            if(idx_tmp==idx_p(dim_i)*2)then
               if(pgno*2<=Maxgrp)pgno=pgno*2
            elseif(idx_tmp==idx_p(dim_i)*2+1)then
               if(pgno*2+1<=Maxgrp)pgno=pgno*2+1
            else
               write(*,*)"something wrong for the value of idx_tmp and idx_p", idx_tmp, idx_p(dim_i)
               stop
            endif

            idx_p(dim_i) = idx_tmp
         enddo
      enddo

      GetMshGroup_Pgno= pgno
   end function GetMshGroup_Pgno


   subroutine blacs_exit_wrp(flag)
      integer flag
#ifdef HAVE_MPI
      call blacs_exit(flag)
#endif
   end subroutine blacs_exit_wrp

   subroutine blacs_gridexit_wrp(ctxt)
      integer ctxt
#ifdef HAVE_MPI
      call blacs_gridexit(ctxt)
#endif
   end subroutine blacs_gridexit_wrp


   subroutine blacs_gridmap_wrp(ctxt, pmap, ldu, nprow, npcol)
      integer ctxt,ldu,nprow,npcol
      integer pmap(:,:)
#ifdef HAVE_MPI
      call blacs_gridmap(ctxt, pmap, ldu, nprow, npcol)
#else
      ctxt=1
#endif
   end subroutine blacs_gridmap_wrp

   integer function sys2blacs_handle_wrp(comm)
      integer comm
#ifdef HAVE_MPI
      integer, external :: sys2blacs_handle
      sys2blacs_handle_wrp=sys2blacs_handle(comm)
#else
      sys2blacs_handle_wrp=1
#endif
   end function sys2blacs_handle_wrp

   subroutine CreatePtree(nmpi, groupmembers, MPI_Comm_base, ptree)
      implicit none
      integer nmpi, MPI_Comm_base, MPI_Group_base, MPI_Group_H, MPI_Group_parent, MPI_Comm_parent, groupmembers(nmpi)
      type(proctree)::ptree, ptreecol, ptreerow
      integer :: ierr, Maxgrp, Maxgrpcol, Maxgrprow, MyID_old, level, group, icontxt, ii, jj, kk, pid
      integer :: nprow, npcol, myrow, mycol, nproc, nproc_tmp, nproc_tmp1, nsprow, nsprow1, nspcol, nlevelrow, nlevelcol
      integer, allocatable::pmap(:, :), groupmembers_sml(:), MPI_Group_H_sml(:)



      call MPI_Comm_rank(MPI_Comm_base, MyID_old, ierr)
      call MPI_Comm_group(MPI_Comm_base, MPI_Group_base, ierr)
      call MPI_Group_incl(MPI_Group_base, nmpi, groupmembers, MPI_Group_H, ierr)
      call MPI_Comm_Create(MPI_Comm_base, MPI_Group_H, ptree%Comm, ierr)

      if (ptree%Comm /= MPI_COMM_NULL) then
         call MPI_Comm_size(ptree%Comm, ptree%nproc, ierr)
         call MPI_Comm_rank(ptree%Comm, ptree%MyID, ierr)

         call assert(groupmembers(ptree%MyID + 1) == MyID_old, 'it is assumed the new ID of the Ith proc in groupmembers is I-1')

         ptree%nlevel = ceiling_safe(log(dble(ptree%nproc))/log(2d0)) + 1
         Maxgrp = 2**(ptree%nlevel) - 1
         allocate (ptree%pgrp(Maxgrp))
         call NumberingPtree(ptree)
         allocate (MPI_Group_H_sml(Maxgrp))
         MPI_Group_H_sml = MPI_GROUP_NULL
         do group = 1, Maxgrp
            nproc = ptree%pgrp(group)%nproc
            if (group == 1) then
               MPI_Comm_parent = ptree%Comm
               MPI_Group_parent = MPI_Group_H
            else
               MPI_Comm_parent = ptree%pgrp(floor_safe(group/2d0))%Comm
               MPI_Group_parent = MPI_Group_H_sml(floor_safe(group/2d0))
            endif

            ! ! create the 2D grids as square as possible
            nprow = floor_safe(sqrt(dble(nproc)))
            npcol = floor_safe(nproc/dble(nprow))

            ! ! the following guarantees column dimension is at most one more level than row dimension, this makes parallel ACA implementation easier
            ! nlevelrow = ceiling_safe(log(dble(nprow)) / log(2d0))+1
            ! nlevelcol = ceiling_safe(log(dble(npcol)) / log(2d0))+1
            ! if(nlevelcol>nlevelrow+1)then
            ! npcol = 2**nprow
            ! endif

            ! ! trail to power of 2 grids, the following two lines can be removed
            ! nproc_tmp = 2**floor_safe(log10(dble(nproc))/log10(2d0))
            ! nprow  = 2**floor_safe(log10(sqrt(dble(nproc_tmp)))/log10(2d0))
            ! npcol = floor_safe(nproc/dble(nprow))

            ptree%pgrp(group)%nprow = nprow
            ptree%pgrp(group)%npcol = npcol
            ptree%pgrp(group)%ctxt = -1
            ptree%pgrp(group)%ctxt1D = -1
            ptree%pgrp(group)%ctxt1DCol = -1
            ptree%pgrp(group)%ctxt_head = -1
            ptree%pgrp(group)%Comm = MPI_COMM_NULL

            if (MPI_Comm_parent /= MPI_COMM_NULL) then

               ! create the local communicator for this tree node
               allocate (groupmembers_sml(ptree%pgrp(group)%nproc))
               do ii = 1, ptree%pgrp(group)%nproc
                  groupmembers_sml(ii) = ii - 1
               enddo
               if (mod(group, 2) == 1 .and. group > 1) then
               if (ptree%pgrp(floor_safe(group/2d0))%nproc > 1) then ! this makes sure the two childs are not the same pgrp as the parent
                  groupmembers_sml = groupmembers_sml + ptree%pgrp(group - 1)%nproc
               endif
               endif
               call MPI_Group_incl(MPI_Group_parent, ptree%pgrp(group)%nproc, groupmembers_sml, MPI_Group_H_sml(group), ierr)
               call MPI_Comm_Create(MPI_Comm_parent, MPI_Group_H_sml(group), ptree%pgrp(group)%Comm, ierr)
               deallocate (groupmembers_sml)

               ! write(*,*)ptree%MyID,group,ptree%pgrp(group)%Comm,ptree%pgrp(group)%Comm==MPI_COMM_NULL,ierr

               ! call MPI_barrier(MPI_Comm_parent,ierr)

            endif

            if (ptree%pgrp(group)%Comm /= MPI_COMM_NULL) then

               allocate (pmap(nprow, npcol))
               do jj = 1, npcol
               do ii = 1, nprow   ! 'row major here'
                  kk = npcol*(ii - 1) + jj
                  pmap(ii, jj) = kk - 1
               enddo
               enddo

               ! the context involving 2D grids
               ptree%pgrp(group)%ctxt = sys2blacs_handle_wrp(ptree%pgrp(group)%Comm)
               call blacs_gridmap_wrp(ptree%pgrp(group)%ctxt, pmap, nprow, nprow, npcol)
               deallocate (pmap)

               ! the context involving 1D grids non-cyclic row
               allocate (pmap(nproc, 1))
               do kk = 1, nproc
                  pmap(kk, 1) = kk - 1
               enddo

               ptree%pgrp(group)%ctxt1D = sys2blacs_handle_wrp(ptree%pgrp(group)%Comm)
               call blacs_gridmap_wrp(ptree%pgrp(group)%ctxt1D, pmap, nproc, nproc, 1)
               deallocate (pmap)

               ! the context involving 1D grids non-cyclic column
               allocate (pmap(1, nproc))
               do kk = 1, nproc
                  pmap(1, kk) = kk - 1
               enddo
               ptree%pgrp(group)%ctxt1DCol = sys2blacs_handle_wrp(ptree%pgrp(group)%Comm)
               call blacs_gridmap_wrp(ptree%pgrp(group)%ctxt1DCol, pmap, 1, 1, nproc)
               deallocate (pmap)

               ! the context involving head proc only
               ptree%pgrp(group)%ctxt_head = sys2blacs_handle_wrp(ptree%pgrp(group)%Comm)
               allocate (pmap(1, 1))
               pmap=0
               call blacs_gridmap_wrp(ptree%pgrp(group)%ctxt_head, pmap, 1, 1, 1)
               deallocate (pmap)

               ! call MPI_barrier(ptree%pgrp(group)%Comm,ierr)

            endif

            ! create the hierarchical process grids used for parallel ACA
            nsprow = ptree%pgrp(group)%nprow
            nspcol = ptree%pgrp(group)%npcol
            ptreecol%nproc = nspcol
            ptreecol%nlevel = ceiling_safe(log(dble(ptreecol%nproc))/log(2d0)) + 1
            Maxgrpcol = 2**(ptreecol%nlevel) - 1
            allocate (ptreecol%pgrp(Maxgrpcol))
            call NumberingPtree(ptreecol)
            ptreerow%nproc = nsprow
            ptreerow%nlevel = ceiling_safe(log(dble(ptreerow%nproc))/log(2d0)) + 1
            Maxgrprow = 2**(ptreerow%nlevel) - 1
            allocate (ptreerow%pgrp(Maxgrprow))
            call NumberingPtree(ptreerow)

            ! allocate(ptree%pgrp(group)%gd)
            ! ptree%pgrp(group)%gd%nsprow=nsprow
            ! ptree%pgrp(group)%gd%nspcol=nspcol
            ! ptree%pgrp(group)%gd%hprow=0
            ! ptree%pgrp(group)%gd%hpcol=0
            ! ptree%pgrp(group)%gd%gprow=1
            ! ptree%pgrp(group)%gd%gpcol=1
            ! call CreateNewGrid(ptree%pgrp(group)%gd,0,ptree,ptreerow,ptreecol,nspcol,MPI_Group_H_sml(group),ptree%pgrp(group)%Comm)

            deallocate (ptreecol%pgrp)
            deallocate (ptreerow%pgrp)

            ! call MPI_barrier(ptree%Comm,ierr)
         enddo

         ! call MPI_barrier(ptree%Comm,ierr)
         do group = 1, Maxgrp
            if (MPI_Group_H_sml(group) /= MPI_GROUP_NULL) call MPI_Group_Free(MPI_Group_H_sml(group), ierr)
         enddo
         deallocate (MPI_Group_H_sml)

         ! call MPI_barrier(ptree%Comm,ierr)

         ! if(ptree%MyID==Main_ID)then
         ! do group=1,Maxgrp
         ! write(*,'(A5,I5,A9,I5,A6,I5,A6,I5,A6,I5,A13,I5,A13,I5)')'myid',ptree%MyID,'group no',group,'nproc',ptree%pgrp(group)%nproc,'nprow',ptree%pgrp(group)%nprow,'npcol',ptree%pgrp(group)%npcol,'grid%nsprow',ptree%pgrp(group)%gd%nsprow,'grid%nspcol',ptree%pgrp(group)%gd%nspcol
         ! enddo
         ! endif

      end if

      ! call MPI_barrier(ptree%Comm,ierr)

      call MPI_Group_Free(MPI_Group_base, ierr)
      call MPI_Group_Free(MPI_Group_H, ierr)

      ! stop
   end subroutine CreatePtree

! !>******** Create a new square grid gd.
! ! Note: Ideally, MPI_Comm_parent should involve smaller process counts with increasing recursion level, but I haven't figured out the local numbering in groupmembers_sml and pmap
! recursive subroutine CreateNewGrid(gd,cridx,ptree,ptreerow,ptreecol,nspcol_parent,MPI_Group_parent,MPI_Comm_parent)
   ! implicit none
   ! type(grid)::gd
   ! integer cridx,Iown,nproc
   ! type(proctree)::ptree,ptreerow,ptreecol
   ! integer nspcol_parent,ii,jj,kk,nsproc
   ! integer MPI_Group_parent,MPI_Group_H_sml,ierr,MPI_Comm_parent
   ! integer,allocatable::pmap(:,:),groupmembers_sml(:)
   ! integer,external :: sys2blacs_handle_wrp

   ! gd%ctxt=-1
   ! gd%Comm=MPI_COMM_NULL
   ! MPI_Group_H_sml=MPI_GROUP_NULL
   ! if(MPI_Comm_parent/=MPI_COMM_NULL)then

   ! ! create the local communicator for this grid
   ! nsproc = gd%nsprow*gd%nspcol
   ! Iown=0
   ! allocate(groupmembers_sml(nsproc))
   ! do jj=1,gd%nspcol
   ! do ii=1,gd%nsprow   ! 'row major here'
   ! kk=nspcol_parent*(gd%hprow+ii-1)+jj+gd%hpcol
   ! groupmembers_sml(jj+(ii-1)*gd%nspcol)=kk-1
   ! ! if(kk+ptree%pgrp(group)%head-1==31)Iown=1
   ! enddo
   ! enddo
   ! ! if(ptree%MyID==0 )write(*,*)'size',nsproc,'cridx',cridx
   ! call MPI_Comm_size(MPI_Comm_parent,nproc,ierr)
   ! call MPI_Group_incl(MPI_Group_parent, nsproc, groupmembers_sml, MPI_Group_H_sml, ierr)
   ! call MPI_Comm_Create(MPI_Comm_parent,MPI_Group_H_sml,gd%Comm,ierr)
   ! deallocate(groupmembers_sml)
   ! ! call MPI_Group_Free(MPI_Group_H_sml,ierr)
   ! endif

   ! if(gd%Comm/=MPI_COMM_NULL)then
   ! allocate(pmap(gd%nsprow,gd%nspcol))
   ! do jj=1,gd%nspcol
   ! do ii=1,gd%nsprow   ! 'row major here'
   ! kk=gd%nspcol*(ii-1)+jj
   ! pmap(ii,jj)=kk-1
   ! enddo
   ! enddo

   ! ! the context involving 2D grids
   ! gd%ctxt = sys2blacs_handle_wrp(gd%Comm)
   ! call blacs_gridmap_wrp( gd%ctxt, pmap, gd%nsprow, gd%nsprow, gd%nspcol )
   ! deallocate(pmap)
   ! endif

   ! if(cridx<ptreerow%nlevel+ptreecol%nlevel-2)then
   ! allocate(gd%gdc(2))
   ! if(mod(cridx+1,2)==1)then
   ! do ii=1,2
   ! gd%gdc(ii)%gprow=gd%gprow
   ! gd%gdc(ii)%hprow=0
   ! gd%gdc(ii)%nsprow=gd%nsprow
   ! gd%gdc(ii)%gpcol=gd%gpcol*2+ii-1
   ! gd%gdc(ii)%hpcol=0
   ! if(ii==2 .and. ptreecol%pgrp(gd%gpcol)%nproc>1)then
   ! gd%gdc(ii)%hpcol= ptreecol%pgrp(gd%gpcol*2)%nproc
   ! endif
   ! gd%gdc(ii)%nspcol=ptreecol%pgrp(gd%gdc(ii)%gpcol)%nproc
   ! enddo
   ! else
   ! do ii=1,2
   ! gd%gdc(ii)%gpcol=gd%gpcol
   ! gd%gdc(ii)%hpcol=0
   ! gd%gdc(ii)%nspcol=gd%nspcol
   ! gd%gdc(ii)%gprow=gd%gprow*2+ii-1
   ! gd%gdc(ii)%hprow=0
   ! if(ii==2 .and. ptreerow%pgrp(gd%gprow)%nproc>1)then
   ! gd%gdc(ii)%hprow= ptreerow%pgrp(gd%gprow*2)%nproc
   ! endif
   ! gd%gdc(ii)%nsprow=ptreerow%pgrp(gd%gdc(ii)%gprow)%nproc
   ! enddo
   ! endif

   ! do ii=1,2
   ! call CreateNewGrid(gd%gdc(ii),cridx+1,ptree,ptreerow,ptreecol,gd%nspcol,MPI_Group_H_sml,gd%Comm)
   ! enddo
   ! endif
   ! if(MPI_Group_H_sml/=MPI_GROUP_NULL)call MPI_Group_Free(MPI_Group_H_sml,ierr)

! end subroutine CreateNewGrid

! redistribute array 1D block array dat_i distributed among process group pgno_i to 1D block array dat_o distributed among process group pgno_o, M_p_i/M_p_o denote the starting index of each process, head_i/head_o denote the global index of the first element (among all processes) in the dat_i/dat_o
   subroutine Redistribute1Dto1D(dat_i, ldi, M_p_i, head_i, pgno_i, dat_o, ldo, M_p_o, head_o, pgno_o, N, ptree, addflag)
      implicit none
      integer ldi, ldo
      DT::dat_i(ldi, *), dat_o(ldo, *)
      integer pgno_i, pgno_o, N
      integer M_p_i(:, :), M_p_o(:, :)
      integer nproc_i, nproc_o, idxs_i, idxs_o, idxe_i, idxe_o, ii, jj, iii, jjj
      type(proctree)::ptree
      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      integer tag, Nreqs, Nreqr, recvid, sendid, ierr, head_i, head_o, sizes, sizer, offs, offr
      integer,optional::addflag

      if (pgno_i == pgno_o .and. ptree%pgrp(pgno_i)%nproc == 1) then
         idxs_i = M_p_i(1, 1) + head_i
         idxe_i = M_p_i(1, 2) + head_i
         idxs_o = M_p_o(1, 1) + head_o
         idxe_o = M_p_o(1, 2) + head_o
         if (idxs_o <= idxe_i .and. idxe_o >= idxs_i) then
            offs = max(idxs_i, idxs_o) - idxs_i
            sizes = min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1
            offr = max(idxs_i, idxs_o) - idxs_o
            sizer = min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1
            if(present(addflag))then
               dat_o(offr + 1:offr + sizer, 1:N) = dat_o(offr + 1:offr + sizer, 1:N) + dat_i(offs + 1:offs + sizes, 1:N)
            else
               dat_o(offr + 1:offr + sizer, 1:N) = dat_i(offs + 1:offs + sizes, 1:N)
            endif
         endif
      else

         nproc_i = ptree%pgrp(pgno_i)%nproc
         nproc_o = ptree%pgrp(pgno_o)%nproc
         tag = pgno_o

         allocate (statuss(MPI_status_size, nproc_o))
         allocate (statusr(MPI_status_size, nproc_i))
         allocate (S_req(nproc_o))
         allocate (R_req(nproc_i))

         allocate (sendquant(nproc_o))
         do ii = 1, nproc_o
            sendquant(ii)%size = 0
         enddo

         allocate (recvquant(nproc_i))
         do ii = 1, nproc_i
            recvquant(ii)%size = 0
         enddo

         if (IOwnPgrp(ptree, pgno_i)) then
            ii = ptree%myid - ptree%pgrp(pgno_i)%head + 1
            idxs_i = M_p_i(ii, 1) + head_i
            idxe_i = M_p_i(ii, 2) + head_i

            do jj = 1, nproc_o
               idxs_o = M_p_o(jj, 1) + head_o
               idxe_o = M_p_o(jj, 2) + head_o
               if (idxs_o <= idxe_i .and. idxe_o >= idxs_i) then
                  sendquant(jj)%offset = max(idxs_i, idxs_o) - idxs_i
                  sendquant(jj)%size = min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1
                  if (sendquant(jj)%size > 0) then
                     allocate (sendquant(jj)%dat(sendquant(jj)%size, N))
                     sendquant(jj)%dat = dat_i(sendquant(jj)%offset + 1:sendquant(jj)%offset + sendquant(jj)%size, 1:N)
                  endif
               endif
            enddo
         endif

         if (IOwnPgrp(ptree, pgno_o)) then
            jj = ptree%myid - ptree%pgrp(pgno_o)%head + 1
            idxs_o = M_p_o(jj, 1) + head_o
            idxe_o = M_p_o(jj, 2) + head_o

            do ii = 1, nproc_i
               idxs_i = M_p_i(ii, 1) + head_i
               idxe_i = M_p_i(ii, 2) + head_i
               if (idxs_o <= idxe_i .and. idxe_o >= idxs_i) then
                  recvquant(ii)%offset = max(idxs_i, idxs_o) - idxs_o
                  recvquant(ii)%size = min(idxe_i, idxe_o) - max(idxs_i, idxs_o) + 1
                  allocate (recvquant(ii)%dat(recvquant(ii)%size, N))
                  recvquant(ii)%dat = 0
               endif
            enddo
         endif

         ! post receive
         Nreqr = 0
         do ii = 1, nproc_i
            if (recvquant(ii)%size > 0) then
               jj = ptree%myid - ptree%pgrp(pgno_o)%head + 1
               sendid = ii + ptree%pgrp(pgno_i)%head - 1
               if (ptree%MyID /= sendid) then
                  Nreqr = Nreqr + 1
                  call MPI_Irecv(recvquant(ii)%dat, recvquant(ii)%size*N, MPI_DT, sendid, tag, ptree%Comm, R_req(Nreqr), ierr)
               endif
            endif
         enddo

         ! post send
         Nreqs = 0
         do jj = 1, nproc_o
            if (sendquant(jj)%size > 0) then
               ii = ptree%myid - ptree%pgrp(pgno_i)%head + 1
               recvid = jj + ptree%pgrp(pgno_o)%head - 1
               if (ptree%MyID == recvid) then
                  recvquant(ii)%dat = sendquant(jj)%dat ! make the direct copy if I own both send and receive pieces
               else
                  Nreqs = Nreqs + 1
                  call MPI_Isend(sendquant(jj)%dat, sendquant(jj)%size*N, MPI_DT, recvid, tag, ptree%Comm, S_req(Nreqs), ierr)
               endif
            endif
         enddo

         if (Nreqs > 0) then
            call MPI_waitall(Nreqs, S_req, statuss, ierr)
         endif
         if (Nreqr > 0) then
            call MPI_waitall(Nreqr, R_req, statusr, ierr)
         endif

         ! copy data from receive buffer
         do ii = 1, nproc_i
            if (recvquant(ii)%size > 0) then
               if(present(addflag))then
                  dat_o(recvquant(ii)%offset + 1:recvquant(ii)%offset + recvquant(ii)%size, 1:N) = dat_o(recvquant(ii)%offset + 1:recvquant(ii)%offset + recvquant(ii)%size, 1:N) + recvquant(ii)%dat
               else
                  dat_o(recvquant(ii)%offset + 1:recvquant(ii)%offset + recvquant(ii)%size, 1:N) = recvquant(ii)%dat
               endif
            endif
         enddo

         ! deallocation
         deallocate (S_req)
         deallocate (R_req)
         deallocate (statuss)
         deallocate (statusr)
         do jj = 1, nproc_o
            if (allocated(sendquant(jj)%dat)) deallocate (sendquant(jj)%dat)
         enddo
         deallocate (sendquant)
         do ii = 1, nproc_i
            if (allocated(recvquant(ii)%dat)) deallocate (recvquant(ii)%dat)
         enddo
         deallocate (recvquant)
      endif

   end subroutine Redistribute1Dto1D



! redistribute array 1D block array dat_i distributed among process group pgno_i to two 1D block array dat_o1 and dat_o2 distributed among process group pgno_o1 and pgno_o1 (which should be the same), M_p_i/M_p_o1/M_p_o2 denote the starting index of each process, head_i/head_o1/head_o2 denote the global index of the first element (among all processes) in the dat_i/dat_o1/dat_o2
   subroutine Redistribute1Dto1D_OnetoTwo(dat_i, ldi, M_p_i, head_i, pgno_i, dat_o1, ldo1, M_p_o1, head_o1, pgno_o1, dat_o2, ldo2, M_p_o2, head_o2, pgno_o2, N, ptree)
      implicit none
      integer ldi,ldo1,ldo2
      DT::dat_i(ldi, *), dat_o1(ldo1, *), dat_o2(ldo2, *)
      integer pgno_i, pgno_o,pgno_o1,pgno_o2, N
      integer M_p_i(:, :), M_p_o1(:, :),M_p_o2(:, :)
      integer nproc_i, nproc_o, idxs_i, idxs_o1, idxs_o2, idxe_i, idxe_o1, idxe_o2, ii, jj, iii, jjj,oo
      type(proctree)::ptree
      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      integer tag, Nreqs, Nreqr, recvid, sendid, ierr, head_i, head_o1, head_o2, sizes, sizer, offs, offr, offset1, offset2,cnt,size1,size2

      call assert(pgno_o1==pgno_o2,'pgno_o1 not equal pgno_o2 in Redistribute1Dto1D_OnetoTwo')
      pgno_o = pgno_o1

      if (pgno_i == pgno_o1 .and. pgno_i == pgno_o2 .and. ptree%pgrp(pgno_i)%nproc == 1) then
         idxs_i = M_p_i(1, 1) + head_i
         idxe_i = M_p_i(1, 2) + head_i
         idxs_o1 = M_p_o1(1, 1) + head_o1
         idxe_o1 = M_p_o1(1, 2) + head_o1
         if (idxs_o1 <= idxe_i .and. idxe_o1 >= idxs_i) then
            offs = max(idxs_i, idxs_o1) - idxs_i
            sizes = min(idxe_i, idxe_o1) - max(idxs_i, idxs_o1) + 1
            offr = max(idxs_i, idxs_o1) - idxs_o1
            sizer = min(idxe_i, idxe_o1) - max(idxs_i, idxs_o1) + 1
            dat_o1(offr + 1:offr + sizer, 1:N) = dat_i(offs + 1:offs + sizes, 1:N)
         endif
         idxs_i = M_p_i(1, 1) + head_i
         idxe_i = M_p_i(1, 2) + head_i
         idxs_o2 = M_p_o2(1, 1) + head_o2
         idxe_o2 = M_p_o2(1, 2) + head_o2
         if (idxs_o2 <= idxe_i .and. idxe_o2 >= idxs_i) then
            offs = max(idxs_i, idxs_o2) - idxs_i
            sizes = min(idxe_i, idxe_o2) - max(idxs_i, idxs_o2) + 1
            offr = max(idxs_i, idxs_o2) - idxs_o2
            sizer = min(idxe_i, idxe_o2) - max(idxs_i, idxs_o2) + 1
            dat_o2(offr + 1:offr + sizer, 1:N) = dat_i(offs + 1:offs + sizes, 1:N)
         endif
      else

         nproc_i = ptree%pgrp(pgno_i)%nproc
         nproc_o = ptree%pgrp(pgno_o)%nproc
         tag = pgno_o

         allocate (statuss(MPI_status_size, nproc_o))
         allocate (statusr(MPI_status_size, nproc_i))
         allocate (S_req(nproc_o))
         allocate (R_req(nproc_i))

         allocate (sendquant(nproc_o))
         do ii = 1, nproc_o
            sendquant(ii)%size = 0
         enddo

         allocate (recvquant(nproc_i))
         do ii = 1, nproc_i
            recvquant(ii)%size = 0
         enddo

         if (IOwnPgrp(ptree, pgno_i)) then
            ii = ptree%myid - ptree%pgrp(pgno_i)%head + 1
            idxs_i = M_p_i(ii, 1) + head_i
            idxe_i = M_p_i(ii, 2) + head_i

            do jj = 1, nproc_o
               sendquant(jj)%size=0
               idxs_o1 = M_p_o1(jj, 1) + head_o1
               idxe_o1 = M_p_o1(jj, 2) + head_o1
               if (idxs_o1 <= idxe_i .and. idxe_o1 >= idxs_i) then
                  size1 = (min(idxe_i, idxe_o1) - max(idxs_i, idxs_o1) + 1)
                  sendquant(jj)%size = sendquant(jj)%size + size1+1
               endif
               idxs_o2 = M_p_o2(jj, 1) + head_o2
               idxe_o2 = M_p_o2(jj, 2) + head_o2
               if (idxs_o2 <= idxe_i .and. idxe_o2 >= idxs_i) then
                  size2 = (min(idxe_i, idxe_o2) - max(idxs_i, idxs_o2) + 1)
                  sendquant(jj)%size = sendquant(jj)%size + size2+1
               endif

               if (sendquant(jj)%size > 0) then
                  allocate (sendquant(jj)%dat(sendquant(jj)%size, N))
               endif

               cnt=0
               if (idxs_o1 <= idxe_i .and. idxe_o1 >= idxs_i) then
                  offset1 = max(idxs_i, idxs_o1) - idxs_i
                  sendquant(jj)%dat(cnt+1,1) = 1
                  ! sendquant(jj)%dat(cnt+2,1) = offset1
                  cnt = cnt +1
                  sendquant(jj)%dat(cnt+1:cnt+size1,1:N) = dat_i(offset1 + 1:offset1 + size1, 1:N)
                  cnt = cnt + size1
               endif
               if (idxs_o2 <= idxe_i .and. idxe_o2 >= idxs_i) then
                  offset2 = max(idxs_i, idxs_o2) - idxs_i
                  sendquant(jj)%dat(cnt+1,1) = 2
                  ! sendquant(jj)%dat(cnt+2,1) = offset2
                  cnt = cnt +1
                  sendquant(jj)%dat(cnt+1:cnt+size2,1:N) = dat_i(offset2 + 1:offset2 + size2, 1:N)
                  cnt = cnt + size2
               endif

            enddo
         endif

         if (IOwnPgrp(ptree, pgno_o)) then
            jj = ptree%myid - ptree%pgrp(pgno_o)%head + 1
            do ii = 1, nproc_i
               idxs_i = M_p_i(ii, 1) + head_i
               idxe_i = M_p_i(ii, 2) + head_i
               recvquant(ii)%size=0

               idxs_o1 = M_p_o1(jj, 1) + head_o1
               idxe_o1 = M_p_o1(jj, 2) + head_o1
               if (idxs_o1 <= idxe_i .and. idxe_o1 >= idxs_i) then
                  ! recvquant(ii)%offset = max(idxs_i, idxs_o) - idxs_o
                  recvquant(ii)%size = recvquant(ii)%size + (min(idxe_i, idxe_o1) - max(idxs_i, idxs_o1) + 1)+1
               endif

               idxs_o2 = M_p_o2(jj, 1) + head_o2
               idxe_o2 = M_p_o2(jj, 2) + head_o2
               if (idxs_o2 <= idxe_i .and. idxe_o2 >= idxs_i) then
                  ! recvquant(ii)%offset = max(idxs_i, idxs_o) - idxs_o
                  recvquant(ii)%size = recvquant(ii)%size + (min(idxe_i, idxe_o2) - max(idxs_i, idxs_o2) + 1)+1
               endif

               if(recvquant(ii)%size>0)then
                  allocate (recvquant(ii)%dat(recvquant(ii)%size, N))
                  recvquant(ii)%dat = 0
               endif

            enddo
         endif

         ! post receive
         Nreqr = 0
         do ii = 1, nproc_i
            if (recvquant(ii)%size > 0) then
               jj = ptree%myid - ptree%pgrp(pgno_o)%head + 1
               sendid = ii + ptree%pgrp(pgno_i)%head - 1
               if (ptree%MyID /= sendid) then
                  Nreqr = Nreqr + 1
                  call MPI_Irecv(recvquant(ii)%dat, recvquant(ii)%size*N, MPI_DT, sendid, tag, ptree%Comm, R_req(Nreqr), ierr)
               endif
            endif
         enddo

         ! post send
         Nreqs = 0
         do jj = 1, nproc_o
            if (sendquant(jj)%size > 0) then
               ii = ptree%myid - ptree%pgrp(pgno_i)%head + 1
               recvid = jj + ptree%pgrp(pgno_o)%head - 1
               if (ptree%MyID == recvid) then
                  recvquant(ii)%dat = sendquant(jj)%dat ! make the direct copy if I own both send and receive pieces
               else
                  Nreqs = Nreqs + 1
                  call MPI_Isend(sendquant(jj)%dat, sendquant(jj)%size*N, MPI_DT, recvid, tag, ptree%Comm, S_req(Nreqs), ierr)
               endif
            endif
         enddo

         if (Nreqs > 0) then
            call MPI_waitall(Nreqs, S_req, statuss, ierr)
         endif
         if (Nreqr > 0) then
            call MPI_waitall(Nreqr, R_req, statusr, ierr)
         endif

         ! copy data from receive buffer
         do ii = 1, nproc_i
            if (recvquant(ii)%size > 0) then
               jj = ptree%myid - ptree%pgrp(pgno_o)%head + 1
               idxs_i = M_p_i(ii, 1) + head_i
               idxe_i = M_p_i(ii, 2) + head_i

               cnt=0
               do while(cnt<recvquant(ii)%size)
               cnt = cnt +1
               oo=NINT(dble(recvquant(ii)%dat(cnt,1)))
               if(oo==1)then
                  idxs_o1 = M_p_o1(jj, 1) + head_o1
                  idxe_o1 = M_p_o1(jj, 2) + head_o1
                  size1 = (min(idxe_i, idxe_o1) - max(idxs_i, idxs_o1) + 1)
                  offset1 = max(idxs_i, idxs_o1) - idxs_o1
                  dat_o1(offset1 + 1:offset1 + size1, 1:N) = recvquant(ii)%dat(cnt+1:cnt+size1,1:N)
                  cnt = cnt +size1
               else if(oo==2)then
                  idxs_o2 = M_p_o2(jj, 1) + head_o2
                  idxe_o2 = M_p_o2(jj, 2) + head_o2
                  size2 = (min(idxe_i, idxe_o2) - max(idxs_i, idxs_o2) + 1)
                  offset2 = max(idxs_i, idxs_o2) - idxs_o2
                  dat_o2(offset2 + 1:offset2 + size2, 1:N) = recvquant(ii)%dat(cnt+1:cnt+size2,1:N)
                  cnt = cnt +size2
               endif
               enddo


            endif
         enddo

         ! deallocation
         deallocate (S_req)
         deallocate (R_req)
         deallocate (statuss)
         deallocate (statusr)
         do jj = 1, nproc_o
            if (allocated(sendquant(jj)%dat)) deallocate (sendquant(jj)%dat)
         enddo
         deallocate (sendquant)
         do ii = 1, nproc_i
            if (allocated(recvquant(ii)%dat)) deallocate (recvquant(ii)%dat)
         enddo
         deallocate (recvquant)
      endif

   end subroutine Redistribute1Dto1D_OnetoTwo



! redistribute two 1D block array dat_i1 and dat_i2 distributed among process group pgno_i1 and pgno_i2 (which should be the same) to one 1D block array dat_o distributed among process group pgno_o, M_p_i1/M_p_i2/M_p_o denote the starting index of each process, head_i1/head_i2/head_o denote the global index of the first element (among all processes) in the dat_i1/dat_i2/dat_o
   subroutine Redistribute1Dto1D_TwotoOne(dat_i1, ldi1, M_p_i1, head_i1, pgno_i1,dat_i2, ldi2, M_p_i2, head_i2, pgno_i2, dat_o, ldo, M_p_o, head_o, pgno_o, N, ptree)
      implicit none
      integer ldi1, ldi2, ldo
      DT::dat_o(ldo, *), dat_i1(ldi1, *), dat_i2(ldi2, *)
      integer pgno_i, pgno_o,pgno_i1,pgno_i2, N
      integer M_p_o(:, :), M_p_i1(:, :),M_p_i2(:, :)
      integer nproc_i, nproc_o, idxs_o, idxs_i1, idxs_i2, idxe_o, idxe_i1, idxe_i2, ii, jj, iii, jjj,ss
      type(proctree)::ptree
      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      integer tag, Nreqs, Nreqr, recvid, sendid, ierr, head_o, head_i1, head_i2, sizes, sizer, offs, offr, offset1, offset2,cnt,size1,size2

      call assert(pgno_i1==pgno_i2,'pgno_i1 not equal pgno_i2 in Redistribute1Dto1D_TwotoOne')
      pgno_i = pgno_i1

      if (pgno_o == pgno_i1 .and. pgno_o == pgno_i2 .and. ptree%pgrp(pgno_o)%nproc == 1) then
         idxs_i1 = M_p_i1(1, 1) + head_i1
         idxe_i1 = M_p_i1(1, 2) + head_i1
         idxs_o = M_p_o(1, 1) + head_o
         idxe_o = M_p_o(1, 2) + head_o
         if (idxs_o <= idxe_i1 .and. idxe_o >= idxs_i1) then
            offs = max(idxs_i1, idxs_o) - idxs_i1
            sizes = min(idxe_i1, idxe_o) - max(idxs_i1, idxs_o) + 1
            offr = max(idxs_i1, idxs_o) - idxs_o
            sizer = min(idxe_i1, idxe_o) - max(idxs_i1, idxs_o) + 1
            dat_o(offr + 1:offr + sizer, 1:N) = dat_i1(offs + 1:offs + sizes, 1:N)
         endif
         idxs_i2 = M_p_i2(1, 1) + head_i2
         idxe_i2 = M_p_i2(1, 2) + head_i2
         idxs_o = M_p_o(1, 1) + head_o
         idxe_o = M_p_o(1, 2) + head_o
         if (idxs_o <= idxe_i2 .and. idxe_o >= idxs_i2) then
            offs = max(idxs_i2, idxs_o) - idxs_i2
            sizes = min(idxe_i2, idxe_o) - max(idxs_i2, idxs_o) + 1
            offr = max(idxs_i2, idxs_o) - idxs_o
            sizer = min(idxe_i2, idxe_o) - max(idxs_i2, idxs_o) + 1
            dat_o(offr + 1:offr + sizer, 1:N) = dat_i2(offs + 1:offs + sizes, 1:N)
         endif
      else

         nproc_i = ptree%pgrp(pgno_i)%nproc
         nproc_o = ptree%pgrp(pgno_o)%nproc
         tag = pgno_o

         allocate (statuss(MPI_status_size, nproc_o))
         allocate (statusr(MPI_status_size, nproc_i))
         allocate (S_req(nproc_o))
         allocate (R_req(nproc_i))

         allocate (sendquant(nproc_o))
         do ii = 1, nproc_o
            sendquant(ii)%size = 0
         enddo

         allocate (recvquant(nproc_i))
         do ii = 1, nproc_i
            recvquant(ii)%size = 0
         enddo

         if (IOwnPgrp(ptree, pgno_i)) then
            ii = ptree%myid - ptree%pgrp(pgno_i)%head + 1
            do jj = 1, nproc_o
               sendquant(jj)%size=0
               idxs_o = M_p_o(jj, 1) + head_o
               idxe_o = M_p_o(jj, 2) + head_o

               idxs_i1 = M_p_i1(ii, 1) + head_i1
               idxe_i1 = M_p_i1(ii, 2) + head_i1
               if (idxs_o <= idxe_i1 .and. idxe_o >= idxs_i1) then
                  size1 = (min(idxe_i1, idxe_o) - max(idxs_i1, idxs_o) + 1)
                  sendquant(jj)%size = sendquant(jj)%size + size1+1
               endif
               idxs_i2 = M_p_i2(ii, 1) + head_i2
               idxe_i2 = M_p_i2(ii, 2) + head_i2
               if (idxs_o <= idxe_i2 .and. idxe_o >= idxs_i2) then
                  size2 = (min(idxe_i2, idxe_o) - max(idxs_i2, idxs_o) + 1)
                  sendquant(jj)%size = sendquant(jj)%size + size2+1
               endif

               if (sendquant(jj)%size > 0) then
                  allocate (sendquant(jj)%dat(sendquant(jj)%size, N))
               endif

               cnt=0
               if (idxs_o <= idxe_i1 .and. idxe_o >= idxs_i1) then
                  offset1 = max(idxs_i1, idxs_o) - idxs_i1
                  sendquant(jj)%dat(cnt+1,1) = 1
                  cnt = cnt +1
                  sendquant(jj)%dat(cnt+1:cnt+size1,1:N) = dat_i1(offset1 + 1:offset1 + size1, 1:N)
                  cnt = cnt + size1
               endif
               if (idxs_o <= idxe_i2 .and. idxe_o >= idxs_i2) then
                  offset2 = max(idxs_i2, idxs_o) - idxs_i2
                  sendquant(jj)%dat(cnt+1,1) = 2
                  cnt = cnt +1
                  sendquant(jj)%dat(cnt+1:cnt+size2,1:N) = dat_i2(offset2 + 1:offset2 + size2, 1:N)
                  cnt = cnt + size2
               endif

            enddo
         endif

         if (IOwnPgrp(ptree, pgno_o)) then
            jj = ptree%myid - ptree%pgrp(pgno_o)%head + 1
            do ii = 1, nproc_i
               recvquant(ii)%size=0

               idxs_o = M_p_o(jj, 1) + head_o
               idxe_o = M_p_o(jj, 2) + head_o

               idxs_i1 = M_p_i1(ii, 1) + head_i1
               idxe_i1 = M_p_i1(ii, 2) + head_i1
               if (idxs_o <= idxe_i1 .and. idxe_o >= idxs_i1) then
                  recvquant(ii)%size = recvquant(ii)%size + (min(idxe_i1, idxe_o) - max(idxs_i1, idxs_o) + 1)+1
               endif
               idxs_i2 = M_p_i2(ii, 1) + head_i2
               idxe_i2 = M_p_i2(ii, 2) + head_i2
               if (idxs_o <= idxe_i2 .and. idxe_o >= idxs_i2) then
                  recvquant(ii)%size = recvquant(ii)%size + (min(idxe_i2, idxe_o) - max(idxs_i2, idxs_o) + 1)+1
               endif

               if(recvquant(ii)%size>0)then
                  allocate (recvquant(ii)%dat(recvquant(ii)%size, N))
                  recvquant(ii)%dat = 0
               endif

            enddo
         endif

         ! post receive
         Nreqr = 0
         do ii = 1, nproc_i
            if (recvquant(ii)%size > 0) then
               jj = ptree%myid - ptree%pgrp(pgno_o)%head + 1
               sendid = ii + ptree%pgrp(pgno_i)%head - 1
               if (ptree%MyID /= sendid) then
                  Nreqr = Nreqr + 1
                  call MPI_Irecv(recvquant(ii)%dat, recvquant(ii)%size*N, MPI_DT, sendid, tag, ptree%Comm, R_req(Nreqr), ierr)
               endif
            endif
         enddo

         ! post send
         Nreqs = 0
         do jj = 1, nproc_o
            if (sendquant(jj)%size > 0) then
               ii = ptree%myid - ptree%pgrp(pgno_i)%head + 1
               recvid = jj + ptree%pgrp(pgno_o)%head - 1
               if (ptree%MyID == recvid) then
                  recvquant(ii)%dat = sendquant(jj)%dat ! make the direct copy if I own both send and receive pieces
               else
                  Nreqs = Nreqs + 1
                  call MPI_Isend(sendquant(jj)%dat, sendquant(jj)%size*N, MPI_DT, recvid, tag, ptree%Comm, S_req(Nreqs), ierr)
               endif
            endif
         enddo

         if (Nreqs > 0) then
            call MPI_waitall(Nreqs, S_req, statuss, ierr)
         endif
         if (Nreqr > 0) then
            call MPI_waitall(Nreqr, R_req, statusr, ierr)
         endif

         ! copy data from receive buffer
         do ii = 1, nproc_i
            if (recvquant(ii)%size > 0) then
               jj = ptree%myid - ptree%pgrp(pgno_o)%head + 1
               idxs_o = M_p_o(jj, 1) + head_o
               idxe_o = M_p_o(jj, 2) + head_o

               cnt=0
               do while(cnt<recvquant(ii)%size)
               cnt = cnt +1
               ss=NINT(dble(recvquant(ii)%dat(cnt,1)))
               if(ss==1)then
                  idxs_i1 = M_p_i1(ii, 1) + head_i1
                  idxe_i1 = M_p_i1(ii, 2) + head_i1
                  size1 = (min(idxe_i1, idxe_o) - max(idxs_i1, idxs_o) + 1)
                  offset1 = max(idxs_i1, idxs_o) - idxs_o
                  dat_o(offset1 + 1:offset1 + size1, 1:N) = recvquant(ii)%dat(cnt+1:cnt+size1,1:N)
                  cnt = cnt +size1
               else if(ss==2)then
                  idxs_i2 = M_p_i2(ii, 1) + head_i2
                  idxe_i2 = M_p_i2(ii, 2) + head_i2
                  size2 = (min(idxe_i2, idxe_o) - max(idxs_i2, idxs_o) + 1)
                  offset2 = max(idxs_i2, idxs_o) - idxs_o
                  dat_o(offset2 + 1:offset2 + size2, 1:N) = recvquant(ii)%dat(cnt+1:cnt+size2,1:N)
                  cnt = cnt +size2
               endif
               enddo
            endif
         enddo

         ! deallocation
         deallocate (S_req)
         deallocate (R_req)
         deallocate (statuss)
         deallocate (statusr)
         do jj = 1, nproc_o
            if (allocated(sendquant(jj)%dat)) deallocate (sendquant(jj)%dat)
         enddo
         deallocate (sendquant)
         do ii = 1, nproc_i
            if (allocated(recvquant(ii)%dat)) deallocate (recvquant(ii)%dat)
         enddo
         deallocate (recvquant)
      endif

   end subroutine Redistribute1Dto1D_TwotoOne




! redistribute 1D block array dat_i distributed among process group pgno_i to 2D block array dat_o distributed among process group pgno_o
   subroutine Redistribute1Dto2D(dat_i, M_p_i, head_i, pgno_i, dat_o, M, head_o, pgno_o, N, ptree)
      implicit none
      integer ld
      DT::dat_i(:, :), dat_o(:, :)
      DT, allocatable::dat_1D(:, :)
      integer pgno_i, pgno_o, N, M
      integer M_p_i(:, :)
      integer nproc_i, nproc_o, idxs_i, idxs_o, idxe_i, idxe_o, ii, jj, iii, jjj
      type(proctree)::ptree
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      integer, allocatable:: M_p_1D(:, :)
      integer tag, Nreqs, Nreqr, recvid, sendid, ierr, head_i, head_o, sizes, sizer, offs, offr
      integer ctxt1D, ctxt1DCol, nproc, nprow, npcol, myrow, mycol, nb1Dc, nb1Dr, myArows, myAcols
      integer::desc1D(9), desc2D(9)
      integer::ctxt, info

      ctxt1DCol = ptree%pgrp(pgno_o)%ctxt1DCol
      ctxt1D = ptree%pgrp(pgno_o)%ctxt1D
      nproc = ptree%pgrp(pgno_o)%nproc
      nb1Dc = N
      nb1Dr = ceiling_safe(M/dble(nproc))
      allocate (M_p_1D(nproc, 2))
      do ii = 1, nproc
         M_p_1D(ii, 1) = (ii - 1)*nb1Dr + 1
         M_p_1D(ii, 2) = ii*nb1Dr
      enddo
      M_p_1D(nproc, 2) = M
      jj = ptree%myid - ptree%pgrp(pgno_o)%head + 1
      myArows = M_p_1D(jj, 2) - M_p_1D(jj, 1) + 1
      myAcols = N
      call descinit_wp(desc1D, M, N, nb1Dr, nb1Dc, 0, 0, ctxt1D, max(myArows, 1), info)
! write(*,*)ptree%MyID,M,N,myArows,myAcols,'1D',nproc
      ld=0
      if (myArows > 0 .and. myAcols > 0) then
         allocate (dat_1D(max(1,myArows), max(1,myAcols)))
         ld = max(1,myArows)
         dat_1D = 0
      endif
      ctxt = ptree%pgrp(pgno_o)%ctxt
      call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
      if (myrow /= -1 .and. mycol /= -1) then
         myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(N, nbslpk, mycol, 0, npcol)
         call descinit_wp(desc2D, M, N, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
      else
         desc2D(2) = -1
      endif

      call Redistribute1Dto1D(dat_i, size(dat_i,1), M_p_i, head_i, pgno_i, dat_1D, ld, M_p_1D, head_o, pgno_o, N, ptree)

! write(*,*)ptree%MyID,M,N,myArows,myAcols,size(dat_1D,1),size(dat_1D,2),size(dat_o,1),size(dat_o,2),'2D'
      call pgemr2df90(M, N, dat_1D, 1, 1, desc1D, dat_o, 1, 1, desc2D, ctxt1DCol)

      deallocate (M_p_1D)
      if (allocated(dat_1D)) deallocate (dat_1D)

   end subroutine Redistribute1Dto2D

! redistribute 2D block array dat_i distributed among process group pgno_i to 1D block array dat_o distributed among process group pgno_o
   subroutine Redistribute2Dto1D(dat_i, M, head_i, pgno_i, dat_o, M_p_o, head_o, pgno_o, N, ptree)
      implicit none
      integer ld
      DT::dat_i(:, :), dat_o(:, :)
      DT, allocatable::dat_1D(:, :)
      integer pgno_i, pgno_o, N, M
      integer M_p_o(:, :)
      integer nproc_i, nproc_o, idxs_i, idxs_o, idxe_i, idxe_o, ii, jj, iii, jjj
      type(proctree)::ptree
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      integer, allocatable:: M_p_1D(:, :)
      integer tag, Nreqs, Nreqr, recvid, sendid, ierr, head_i, head_o, sizes, sizer, offs, offr
      integer ctxt1D, ctxt1DCol, nproc, nprow, npcol, myrow, mycol, nb1Dc, nb1Dr, myArows, myAcols
      integer::desc1D(9), desc2D(9)
      integer::ctxt, info

      ctxt1DCol = ptree%pgrp(pgno_i)%ctxt1DCol
      ctxt1D = ptree%pgrp(pgno_i)%ctxt1D
      nproc = ptree%pgrp(pgno_i)%nproc
      nb1Dc = N
      nb1Dr = ceiling_safe(M/dble(nproc))
      allocate (M_p_1D(nproc, 2))
      do ii = 1, nproc
         M_p_1D(ii, 1) = (ii - 1)*nb1Dr + 1
         M_p_1D(ii, 2) = ii*nb1Dr
      enddo
      M_p_1D(nproc, 2) = M
      jj = ptree%myid - ptree%pgrp(pgno_i)%head + 1
      myArows = M_p_1D(jj, 2) - M_p_1D(jj, 1) + 1
      myAcols = N
      call descinit_wp(desc1D, M, N, nb1Dr, nb1Dc, 0, 0, ctxt1D, max(myArows, 1), info)
! write(*,*)ptree%MyID,M,N,myArows,myAcols,'1D',nproc
      ld=0
      if (myArows > 0 .and. myAcols > 0) then
         allocate (dat_1D(max(1,myArows), max(1,myAcols)))
         ld = max(1,myArows)
         dat_1D = 0
      endif

      ctxt = ptree%pgrp(pgno_i)%ctxt
      call blacs_gridinfo_wrp(ctxt, nprow, npcol, myrow, mycol)
      if (myrow /= -1 .and. mycol /= -1) then
         myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
         myAcols = numroc_wp(N, nbslpk, mycol, 0, npcol)
         call descinit_wp(desc2D, M, N, nbslpk, nbslpk, 0, 0, ctxt, max(myArows, 1), info)
      else
         myArows = 0
         myAcols = 0
         desc2D(2) = -1
      endif

! write(*,*)ptree%MyID,M,N,myArows,myAcols,size(dat_i,1),size(dat_i,2),size(dat_1D,1),size(dat_1D,2),'2D2D',isnanMat(dat_i,size(dat_i,1),size(dat_i,2)),isnanMat(dat_1D,size(dat_1D,1),size(dat_1D,2)),myrow,mycol,pgno_i,ctxt
      call pgemr2df90(M, N, dat_i, 1, 1, desc2D, dat_1D, 1, 1, desc1D, ctxt1DCol)
! write(*,*)ptree%MyID,M,N,myArows,myAcols,size(dat_1D,1),size(dat_1D,2),size(dat_o,1),size(dat_o,2),'1D1D',isnanMat(dat_1D,size(dat_1D,1),size(dat_1D,2)),isnanMat(dat_o,size(dat_o,1),size(dat_o,2)),myrow,mycol
      call Redistribute1Dto1D(dat_1D, ld, M_p_1D, head_i, pgno_i, dat_o, size(dat_o,1), M_p_o, head_o, pgno_o, N, ptree)

      deallocate (M_p_1D)
      if (allocated(dat_1D)) deallocate (dat_1D)

   end subroutine Redistribute2Dto1D



! redistribute array 1D block array (the array has shape like n^d x k and is redistributed among n^d) dat_i distributed among process group pgno_i to 1D block array dat_o distributed among process group pgno_o, M_p_i/M_p_o denote the starting index of each process of each dimension, head_i/head_o denote the global index of the first element (among all processes) in the dat_i/dat_o of each dimension
   subroutine Redistribute1Dto1D_MD(Ndim, dat_i, ldi, M_p_i, head_i, pgno_i, dat_o, ldo, M_p_o, head_o, pgno_o, N, ptree, addflag)
      implicit none
      integer Ndim
      integer ldi(Ndim), ldo(Ndim), dim_i, idx_MD(Ndim), idx_s_MD(Ndim),idx_r_MD(Ndim), dims_md(Ndim),idx_r_scalar,idx_s_scalar
      DT::dat_i(product(ldi), *), dat_o(product(ldo), *)
      integer pgno_i, pgno_o, N
      integer M_p_i(:, :, :), M_p_o(:, :, :)
      integer nproc_i, nproc_o, idxs_i(Ndim), idxs_o(Ndim), idxe_i(Ndim), idxe_o(Ndim), ii, jj, iii, jjj
      type(proctree)::ptree
      type(commquant1D), allocatable::sendquant(:), recvquant(:)
      integer, allocatable::S_req(:), R_req(:)
      integer, allocatable:: statuss(:, :), statusr(:, :)
      integer tag, Nreqs, Nreqr, recvid, sendid, ierr, head_i(Ndim), head_o(Ndim), sizes(Ndim), sizer(Ndim), offs(Ndim), offr(Ndim)
      integer,optional::addflag

      if (pgno_i == pgno_o .and. ptree%pgrp(pgno_i)%nproc == 1) then
         idxs_i = M_p_i(1, 1, :) + head_i
         idxe_i = M_p_i(1, 2, :) + head_i
         idxs_o = M_p_o(1, 1, :) + head_o
         idxe_o = M_p_o(1, 2, :) + head_o
         if (ALL(idxs_o <= idxe_i) .and. ALL(idxe_o >= idxs_i)) then
            do dim_i=1,Ndim
               offs(dim_i) = max(idxs_i(dim_i), idxs_o(dim_i)) - idxs_i(dim_i)
               sizes(dim_i) = min(idxe_i(dim_i), idxe_o(dim_i)) - max(idxs_i(dim_i), idxs_o(dim_i)) + 1
               offr(dim_i) = max(idxs_i(dim_i), idxs_o(dim_i)) - idxs_o(dim_i)
               sizer(dim_i) = sizes(dim_i)
            enddo

            do iii=1,product(sizes)
               call SingleIndexToMultiIndex(Ndim, sizes, iii, idx_MD)
               idx_r_MD = idx_MD + offr
               idx_s_MD = idx_MD + offs
               call MultiIndexToSingleIndex(Ndim, ldi, idx_s_scalar, idx_s_MD)
               call MultiIndexToSingleIndex(Ndim, ldo, idx_r_scalar, idx_r_MD)
               if(present(addflag))then
                  dat_o(idx_r_scalar, 1:N) = dat_o(idx_r_scalar, 1:N) + dat_i(idx_s_scalar, 1:N)
               else
                  dat_o(idx_r_scalar, 1:N) = dat_i(idx_s_scalar, 1:N)
               endif
            enddo
         endif
      else

         nproc_i = ptree%pgrp(pgno_i)%nproc
         nproc_o = ptree%pgrp(pgno_o)%nproc
         tag = pgno_o

         allocate (statuss(MPI_status_size, nproc_o))
         allocate (statusr(MPI_status_size, nproc_i))
         allocate (S_req(nproc_o))
         allocate (R_req(nproc_i))

         allocate (sendquant(nproc_o))
         do ii = 1, nproc_o
            allocate(sendquant(ii)%offset_md(Ndim))
            allocate(sendquant(ii)%size_md(Ndim))
            sendquant(ii)%size_md = 0
            sendquant(ii)%offset_md = -1
         enddo

         allocate (recvquant(nproc_i))
         do ii = 1, nproc_i
            allocate(recvquant(ii)%offset_md(Ndim))
            allocate(recvquant(ii)%size_md(Ndim))
            recvquant(ii)%size_md = 0
            recvquant(ii)%offset_md = -1
         enddo

         if (IOwnPgrp(ptree, pgno_i)) then
            ii = ptree%myid - ptree%pgrp(pgno_i)%head + 1
            idxs_i = M_p_i(ii, 1, :) + head_i
            idxe_i = M_p_i(ii, 2, :) + head_i

            do jj = 1, nproc_o
               idxs_o = M_p_o(jj, 1, :) + head_o
               idxe_o = M_p_o(jj, 2, :) + head_o
               if (ALL(idxs_o <= idxe_i) .and. ALL(idxe_o >= idxs_i)) then
                  do dim_i=1,Ndim
                     sendquant(jj)%offset_md(dim_i) = max(idxs_i(dim_i), idxs_o(dim_i)) - idxs_i(dim_i)
                     sendquant(jj)%size_md(dim_i) = min(idxe_i(dim_i), idxe_o(dim_i)) - max(idxs_i(dim_i), idxs_o(dim_i)) + 1
                  enddo
                  if (ALL(sendquant(jj)%size_md > 0)) then
                     allocate (sendquant(jj)%dat(product(sendquant(jj)%size_md), N))
                     do iii=1,product(sendquant(jj)%size_md)
                        call SingleIndexToMultiIndex(Ndim, sendquant(jj)%size_md, iii, idx_MD)
                        idx_s_MD = idx_MD + sendquant(jj)%offset_md
                        call MultiIndexToSingleIndex(Ndim, ldi, idx_s_scalar, idx_s_MD)
                        sendquant(jj)%dat(iii,1:N) = dat_i(idx_s_scalar, 1:N)
                     enddo
                  endif
               endif
            enddo
         endif

         if (IOwnPgrp(ptree, pgno_o)) then
            jj = ptree%myid - ptree%pgrp(pgno_o)%head + 1
            idxs_o = M_p_o(jj, 1, :) + head_o
            idxe_o = M_p_o(jj, 2, :) + head_o

            do ii = 1, nproc_i
               idxs_i = M_p_i(ii, 1, :) + head_i
               idxe_i = M_p_i(ii, 2, :) + head_i
               if (ALL(idxs_o <= idxe_i) .and. ALL(idxe_o >= idxs_i)) then
                  do dim_i=1,Ndim
                     recvquant(ii)%offset_md(dim_i) = max(idxs_i(dim_i), idxs_o(dim_i)) - idxs_o(dim_i)
                     recvquant(ii)%size_md(dim_i) = min(idxe_i(dim_i), idxe_o(dim_i)) - max(idxs_i(dim_i), idxs_o(dim_i)) + 1
                  enddo
                  if (ALL(recvquant(ii)%size_md > 0)) then
                     allocate (recvquant(ii)%dat(product(recvquant(ii)%size_md), N))
                     recvquant(ii)%dat = 0
                  endif
               endif
            enddo
         endif

         ! post receive
         Nreqr = 0
         do ii = 1, nproc_i
            if (ALL(recvquant(ii)%size_md > 0)) then
               jj = ptree%myid - ptree%pgrp(pgno_o)%head + 1
               sendid = ii + ptree%pgrp(pgno_i)%head - 1
               if (ptree%MyID /= sendid) then
                  Nreqr = Nreqr + 1
                  call MPI_Irecv(recvquant(ii)%dat, product(recvquant(ii)%size_md)*N, MPI_DT, sendid, tag, ptree%Comm, R_req(Nreqr), ierr)
               endif
            endif
         enddo

         ! post send
         Nreqs = 0
         do jj = 1, nproc_o
            if (ALL(sendquant(jj)%size_md > 0)) then
               ii = ptree%myid - ptree%pgrp(pgno_i)%head + 1
               recvid = jj + ptree%pgrp(pgno_o)%head - 1
               if (ptree%MyID == recvid) then
                  recvquant(ii)%dat = sendquant(jj)%dat ! make the direct copy if I own both send and receive pieces
               else
                  Nreqs = Nreqs + 1
                  call MPI_Isend(sendquant(jj)%dat, product(sendquant(jj)%size_md)*N, MPI_DT, recvid, tag, ptree%Comm, S_req(Nreqs), ierr)
               endif
            endif
         enddo

         if (Nreqs > 0) then
            call MPI_waitall(Nreqs, S_req, statuss, ierr)
         endif
         if (Nreqr > 0) then
            call MPI_waitall(Nreqr, R_req, statusr, ierr)
         endif

         ! copy data from receive buffer
         do ii = 1, nproc_i
            if (ALL(recvquant(ii)%size_md > 0)) then
               do iii=1,product(recvquant(ii)%size_md)
                  call SingleIndexToMultiIndex(Ndim, recvquant(ii)%size_md, iii, idx_MD)
                  idx_r_MD = idx_MD + recvquant(ii)%offset_md
                  call MultiIndexToSingleIndex(Ndim, ldo, idx_r_scalar, idx_r_MD)
                  if(present(addflag))then
                     dat_o(idx_r_scalar, 1:N) = dat_o(idx_r_scalar, 1:N) + recvquant(ii)%dat(iii,1:N)
                  else
                     dat_o(idx_r_scalar, 1:N) = recvquant(ii)%dat(iii,1:N)
                  endif
               enddo
            endif
         enddo

         ! deallocation
         deallocate (S_req)
         deallocate (R_req)
         deallocate (statuss)
         deallocate (statusr)
         do jj = 1, nproc_o
            if (allocated(sendquant(jj)%dat)) deallocate (sendquant(jj)%dat)
            if (allocated(sendquant(jj)%offset_md)) deallocate (sendquant(jj)%offset_md)
            if (allocated(sendquant(jj)%size_md)) deallocate (sendquant(jj)%size_md)
         enddo
         deallocate (sendquant)
         do ii = 1, nproc_i
            if (allocated(recvquant(ii)%dat)) deallocate (recvquant(ii)%dat)
            if (allocated(recvquant(ii)%offset_md)) deallocate (recvquant(ii)%offset_md)
            if (allocated(recvquant(ii)%size_md)) deallocate (recvquant(ii)%size_md)
         enddo
         deallocate (recvquant)
      endif

   end subroutine Redistribute1Dto1D_MD





! get the level (indexed from 1) of a node in a tree. gno is the node number starting from root (1). Note that the process tree levels are indexed from 1, the basis_group levels are indexed from 0
   integer function GetTreelevel(gno)
      implicit none
      integer gno, ii, level
      ii = gno
      level = 0
      do while (ii /= 0)
         ii = INT(ii/2d0)
         level = level + 1
      enddo
      GetTreelevel = level
   end function GetTreelevel

! check if I share this process group
   logical function IOwnPgrp(ptree, pgno)
      implicit none
      integer pgno
      type(proctree)::ptree
      if(pgno==-1)then
         IOwnPgrp = .false.
      else
         IOwnPgrp = ptree%MyID >= ptree%pgrp(pgno)%head .and. ptree%MyID <= ptree%pgrp(pgno)%tail
      endif
   end function IOwnPgrp


   integer function lcm(a, b)
      integer:: a, b
      lcm = a*b/gcd(a, b)
   end function lcm

   integer function gcd(a, b)
      integer :: a, b, t, as, bs
      as = a
      bs = b
      do while (bs /= 0)
         t = bs
         bs = mod(as, bs)
         as = t
      end do
      gcd = abs(as)
   end function gcd



   integer function blacs_pnum_wp(NPROW, NPCOL, PROW, PCOL)
      integer :: ICONTXT, PROW, PCOL, NPROW, NPCOL
      ! integer :: blacs_pnum ! blacs routine
      ! blacs_pnum_wp = blacs_pnum(ICONTXT, PROW, PCOL)

      ! 'row major here'
      blacs_pnum_wp = NPCOL*PROW + PCOL

   end function blacs_pnum_wp


! get my row and column index in a 2D grid that is column major consisting of contiguous proc ids
   subroutine Gridinfo_2D(pmap, MyID, myrow, mycol)
      implicit none
      integer pmap(3)
      integer myrow, mycol, nprow, npcol, idstart, idend, MyID
      nprow = pmap(1)
      npcol = pmap(2)
      idstart = pmap(3)
      idend = pmap(3) + nprow*npcol - 1
      myrow = -1
      mycol = -1
      if (MyID >= idstart .and. MyID <= idend) then
         mycol = (MyID - idstart)/nprow
         myrow = mod(MyID - idstart, nprow)
      endif
   end subroutine Gridinfo_2D


   subroutine Array1DtoPointer2D(x, p, n1, n2)
      integer::         n1, n2
      DT, target::    x(n1, n2)
      DT, pointer::   p(:,:)
      p => x
   end subroutine Array1DtoPointer2D

   subroutine LogMemory(stats, mem)
      implicit none
      type(Hstat)::stats
      real(kind=8):: mem
      stats%Mem_Current = stats%Mem_Current + mem
      if(mem>0)stats%Mem_Peak = max(stats%Mem_Peak,stats%Mem_Current)
   end subroutine LogMemory




   subroutine get_graph_colors_JP(rows,ia,ja,colors)
      implicit none
      integer ja(:)
      integer rows,ia(rows+1),colors(rows),weights(rows),csp,remaining,i,c,largest,jp,neighbor
      integer,allocatable::forbidden(:)

      colors=0
      call rperm(rows, weights)
      csp = maxval(ia(2:rows+1) - ia(1:rows))
      remaining = rows

      do while(remaining > 0)
          do i = 1,rows
              ! Make sure the node isn't already colored
              if(colors(i) == 0)then

               ! Check if this node has the largest weight amongst its
               ! uncolored neighbors
               largest = 1;
               do jp = ia(i),ia(i+1)-1
                  neighbor = ja(jp)
                  if(colors(neighbor) == 0 .and. weights(i) < weights(neighbor))then
                     largest = 0
                     exit
                  endif
               enddo

               if(largest == 1)then
                  ! Get the smallest color that hasn't been assigned to this
                  ! node's neighbors
                  allocate(forbidden(csp))
                  forbidden=0
                  do jp = ia(i),ia(i+1)-1
                        neighbor = ja(jp)
                        if(colors(neighbor) /= 0)then
                           forbidden(colors(neighbor)) = 1
                        endif
                  enddo
                  do c = 1,csp
                        if (forbidden(c) == 0)then
                           colors(i) = c
                           remaining = remaining - 1
                           exit
                        endif
                  enddo
                  ! This node should have been assigned a color so if it hasn't
                  ! then something went wrong
                  if (colors(i) == 0)then
                        write(*,*)'Error assigning color to node '
                        stop
                  endif
                  deallocate(forbidden)
               endif
            endif
          enddo
      enddo

   end subroutine get_graph_colors_JP

   ! bit reversal operation for an integer n ranging from 1 to 2^bits
integer function bit_reverse(n, bits)
    integer, intent(in) :: n, bits
    integer :: i, result, nwork, bitarray(bits)

	nwork=n-1
	do i=1,bits
		bitarray(bits-i+1)=mod(nwork,2)
		nwork=int(nwork/2)
	enddo

    result = 0
    do i=1,bits
		result = result + 2**(i-1)*bitarray(i)
    end do

    bit_reverse = result+1
end function bit_reverse




#ifdef HAVE_ZFP
   !>**** ZFP compression of fullmat into fullmatZFP%buffer_r/buffer_i into with relative tolerance tol_comp. If the flag already_compressed is set to 1, only the fullmat will be deleted.
   subroutine ZFP_Compress(fullmat, fullmatZFP, M, N, tol_comp, already_compressed)

      implicit none
      DT,pointer::fullmat(:,:)
      integer M,N
      type(zfpquant)::fullmatZFP
      real(kind=8)::tol_comp,tol_abs

      ! input/decompressed arrays
      type(c_ptr) :: array_c_ptr
      DTR, dimension(:, :), allocatable, target :: input_array
      ! zfp_field
      type(zFORp_field) :: field
      ! bitstream
      character, dimension(:), allocatable, target :: buffer
      type(c_ptr) :: buffer_c_ptr
      integer (kind=8) buffer_size_bytes, bitstream_offset_bytes,tmpint8
      type(zFORp_bitstream) :: bitstream
      ! zfp_stream
      type(zFORp_stream) :: stream, dstream
      real (kind=8) :: tol_result, maxcal
      DTR::maxval
      integer :: zfp_type, already_compressed

      if(already_compressed==0)then
         maxval=fnorm(fullmat,M,N,'M')
         tol_abs = tol_comp * maxval

         ! setup zfp_field
         allocate(input_array(M,N))
         input_array=dble(fullmat)
         array_c_ptr = c_loc(input_array)
         zfp_type = DTZFP
         field = zFORp_field_2d(array_c_ptr, zfp_type, M, N)

         ! setup bitstream
         tmpint8 = N*DTRBytes
         buffer_size_bytes = M*tmpint8
         allocate(buffer(buffer_size_bytes*2)) ! adding an extra factor of 2 as the compressed array can be larger than the input data
         buffer_c_ptr = c_loc(buffer)
         bitstream = zFORp_bitstream_stream_open(buffer_c_ptr, buffer_size_bytes)

         ! setup zfp_stream
         fullmatZFP%stream_r = zFORp_stream_open(bitstream)
         tol_result=zFORp_stream_set_accuracy(fullmatZFP%stream_r,tol_abs)

         ! compress
         bitstream_offset_bytes = zFORp_compress(fullmatZFP%stream_r, field)
         allocate(fullmatZFP%buffer_r(bitstream_offset_bytes))
         fullmatZFP%buffer_r=buffer(1:bitstream_offset_bytes)
         call zFORp_field_free(field)
         call zFORp_bitstream_stream_close(bitstream)
         deallocate(input_array)
         deallocate(buffer)


#if DAT==0 || DAT==2
         ! setup zfp_field
         allocate(input_array(M,N))
         input_array=aimag(fullmat)
         array_c_ptr = c_loc(input_array)
         zfp_type = DTZFP
         field = zFORp_field_2d(array_c_ptr, zfp_type, M, N)

         ! setup bitstream
         tmpint8 = N*DTRBytes
         buffer_size_bytes = M*tmpint8
         allocate(buffer(buffer_size_bytes*2)) ! adding an extra factor of 2 as the compressed array can be larger than the input data
         buffer_c_ptr = c_loc(buffer)
         bitstream = zFORp_bitstream_stream_open(buffer_c_ptr, buffer_size_bytes)

         ! setup zfp_stream
         fullmatZFP%stream_i = zFORp_stream_open(bitstream)
         tol_result=zFORp_stream_set_accuracy(fullmatZFP%stream_i,tol_abs)

         ! compress
         bitstream_offset_bytes = zFORp_compress(fullmatZFP%stream_i, field)
         allocate(fullmatZFP%buffer_i(bitstream_offset_bytes))
         fullmatZFP%buffer_i=buffer(1:bitstream_offset_bytes)
         call zFORp_field_free(field)
         call zFORp_bitstream_stream_close(bitstream)
         deallocate(input_array)
         deallocate(buffer)
#endif
      endif
      deallocate(fullmat)
      fullmat=>null()

   end subroutine ZFP_Compress


   !>**** ZFP decompression of fullmatZFP%buffer_r/buffer_i into fullmat. tol_used is the relative tolerance in the compression phase. The flag keep_compressed incidates whether the compressed data will be deleted.
   subroutine ZFP_Decompress(fullmat, fullmatZFP, M, N, tol_used, keep_compressed)

      implicit none
      DT,pointer::fullmat(:,:)
      integer M,N
      type(zfpquant)::fullmatZFP

      ! input/decompressed arrays
      type(c_ptr) :: array_c_ptr
      DTR::maxval
      DTR, dimension(:, :), allocatable, target :: decompressed_array
      ! zfp_field
      type(zFORp_field) :: field
      ! bitstream
      character, dimension(:), allocatable, target :: buffer
      type(c_ptr) :: buffer_c_ptr
      integer (kind=8) buffer_size_bytes, bitstream_offset_bytes
      type(zFORp_bitstream) :: bitstream
      ! zfp_stream
      type(zFORp_stream) :: stream
      real (kind=8) :: tol_result,tol_used
      integer :: zfp_type, keep_compressed

      tol_result = zFORp_stream_accuracy(fullmatZFP%stream_r)

      allocate(fullmat(M,N))

      ! setup zfp_field
      allocate(decompressed_array(M,N))
      array_c_ptr = c_loc(decompressed_array)
      zfp_type = DTZFP
      field = zFORp_field_2d(array_c_ptr, zfp_type, M, N)

      ! setup bitstream
      buffer_size_bytes=size(fullmatZFP%buffer_r,1)
      allocate(buffer(buffer_size_bytes))
      buffer = fullmatZFP%buffer_r
      buffer_c_ptr = c_loc(buffer)
      bitstream = zFORp_bitstream_stream_open(buffer_c_ptr, buffer_size_bytes)
      call zFORp_stream_set_bit_stream(fullmatZFP%stream_r, bitstream);

      ! decompress
      bitstream_offset_bytes = zFORp_decompress(fullmatZFP%stream_r, field)
      fullmat=decompressed_array

      ! deallocations
      call zFORp_bitstream_stream_close(bitstream)
      call zFORp_field_free(field)

      deallocate(decompressed_array)
      deallocate(buffer)
      if(keep_compressed==0)then
         call zFORp_stream_close(fullmatZFP%stream_r)
         deallocate(fullmatZFP%buffer_r)
      endif

#if DAT==0 || DAT==2
      ! setup zfp_field
      allocate(decompressed_array(M,N))
      array_c_ptr = c_loc(decompressed_array)
      zfp_type = DTZFP
      field = zFORp_field_2d(array_c_ptr, zfp_type, M, N)

      ! setup bitstream
      buffer_size_bytes=size(fullmatZFP%buffer_i,1)
      allocate(buffer(buffer_size_bytes))
      buffer = fullmatZFP%buffer_i
      buffer_c_ptr = c_loc(buffer)
      bitstream = zFORp_bitstream_stream_open(buffer_c_ptr, buffer_size_bytes)
      call zFORp_stream_set_bit_stream(fullmatZFP%stream_i, bitstream);

      ! decompress
      bitstream_offset_bytes = zFORp_decompress(fullmatZFP%stream_i, field)
      fullmat=fullmat+BPACK_junit*decompressed_array

      ! deallocations
      call zFORp_bitstream_stream_close(bitstream)
      call zFORp_field_free(field)

      deallocate(decompressed_array)
      deallocate(buffer)

      if(keep_compressed==0)then
         call zFORp_stream_close(fullmatZFP%stream_i)
         deallocate(fullmatZFP%buffer_i)
      endif
#endif

   maxval=fnorm(fullmat,M,N,'M')
   if(maxval<BPACK_SafeEps)then
      tol_used = BPACK_SafeEps
   else
      tol_used = tol_result/maxval
   endif

   end subroutine ZFP_Decompress

#endif



!>**** TT decomposition via SVD of a full array
!b: the d-dimensional input tensor, stored as 1D array
!eps_in: the relative compression tolerance
!t: the TT decomposition. At input, t contains only d, n (and m_n)
!ZFP_flag: whether to further compress the core with zfp
subroutine TT_Compress_SVD(b, eps_in, t, ZFP_flag)
   implicit none
   DT   :: b(:)
   real*8   :: eps_in, tmpsize
   type(TTtype) :: t
   integer,optional :: ZFP_flag
   integer :: d, i, ii,jj, dd
   integer :: pos
   integer :: m, chunk_n,mn_min
   DT, allocatable :: c(:)
   DTR, allocatable :: s(:)
   DT, allocatable :: u(:,:), vt(:,:), mat(:,:)
   integer :: lwork, info, ranki
   integer :: old_size, new_size, final_block, block_size, shape_flat(1)
   integer,allocatable :: shape_hd(:), perm(:)
   real*8  :: ep_s
   DT, pointer:: data_buffer(:,:), data_buffer1(:,:)
   type(zfpquant)::tmpzfp

   integer :: rowdim, coldim
   integer :: lda, ldu, ldvt

   integer :: idx, currentsize

   d = size(t%n)

   allocate(t%r(d+1))
   t%r = 1
   allocate(t%psa(d+1))
   t%psa = 0

   t%d = d

   allocate(c(size(b)))
   c = b


   if(t%mpo==1)then ! MPO
      allocate(shape_hd(2*d))
      shape_hd(1:d) = t%m_n(:,1)
      shape_hd(1+d:2*d) = t%m_n(:,2)
      allocate(perm(2*d))
      perm = 0
      do idx=1,d
        perm(2*(idx-1)+1 ) = idx            ! row bits
        perm(2*(idx-1)+2 ) = d + idx        ! col bits
      enddo
      call TensorPermute(c, shape_hd, perm)
      deallocate(perm)
      deallocate(shape_hd)
   endif




   ! ep = eps_in/sqrt(d-1)  (avoid d=1 corner case)
   if (d>1) then
      ep_s = eps_in / dsqrt(d-1.d0)
   else
      ep_s = eps_in
   end if

   pos = 1

! write(*,*)'c',c(1:512)
   tmpsize=0
   !=== TT-SVD main loop for first d-1 cores
   do i=1,d-1
      m       = t%n(i)*t%r(i)
   !  write(*,*)i,d,t%n(i),t%r(i),'dim'
      chunk_n = size(c)/m


      ! reshape c => (m x chunk_n)
      allocate(mat(m, chunk_n))
      mat = reshape(c, shape(mat))
   !  write(*,*)'c', c
   !  write(*,*)'re c', mat

      ! SVD of mat => U*S*V^T
      mn_min = min(m,chunk_n)
      allocate(u(m,mn_min), vt(mn_min,chunk_n))
      allocate(s(mn_min))


      call SVD_Truncate(mat, m, chunk_n, mn_min, u, vt, s, ep_s, BPACK_SafeUnderflow, ranki)

      t%r(i+1) = ranki

      ! Store U(:,1:ranki) into TT core
      currentsize=0
      if(allocated(t%core))currentsize=size(t%core)
      block_size = t%r(i)*t%n(i)*ranki
      call array_resize(t%core, currentsize + block_size)
      shape_flat(1) = block_size
      t%core(pos : pos+block_size-1) = reshape(u(1:m,1:ranki), shape_flat)

      ! allocate(data_buffer1(block_size,1))
      ! data_buffer1(:,1)=t%core(pos : pos+block_size-1)
      ! call ZFP_Compress(data_buffer1,tmpzfp,block_size,1, ep_s,0)
      ! write(*,*)'core',i,block_size*16, (SIZEOF(tmpzfp%buffer_r)+SIZEOF(tmpzfp%buffer_i)),dble(block_size*16)/(SIZEOF(tmpzfp%buffer_r)+SIZEOF(tmpzfp%buffer_i))
      ! tmpsize = tmpsize + (SIZEOF(tmpzfp%buffer_r)+SIZEOF(tmpzfp%buffer_i))
      ! deallocate(tmpzfp%buffer_r)
      ! deallocate(tmpzfp%buffer_i)

      pos = pos + block_size

      ! Construct c = V(1:ranki, :) * diag(s(1:ranki))
   do ii=1, ranki
      vt(ii,:) = vt(ii,:)*s(ii)
   end do
   deallocate(c)
   allocate(c(ranki*chunk_n))
   c = reshape(vt(1:ranki,1:chunk_n), shape(c))

      deallocate(u, vt, s, mat)
   end do

   !=== Last core
   ! We have c(:) dimension = ( t%r(d)*t%n(d)*r(d+1) ) => but r(d+1) might be unknown until now
   ! r(d+1) should come from size(c)/( t%r(d)* t%n(d) )
   t%r(d+1) = size(c) / ( t%r(d)*t%n(d) )
   final_block = t%r(d)*t%n(d)*t%r(d+1)

   currentsize=0
   if(allocated(t%core))currentsize=size(t%core)
   call array_resize(t%core, currentsize + final_block)
   t%core(pos: pos+final_block-1) = c(:)

      ! allocate(data_buffer1(final_block,1))
      ! data_buffer1(:,1)=t%core(pos : pos+final_block-1)
      ! call ZFP_Compress(data_buffer1,tmpzfp,final_block,1, eps_in,0)
      ! write(*,*)'core L',final_block*16, (SIZEOF(tmpzfp%buffer_r)+SIZEOF(tmpzfp%buffer_i)),dble(final_block*16)/(SIZEOF(tmpzfp%buffer_r)+SIZEOF(tmpzfp%buffer_i))
      ! tmpsize = tmpsize + (SIZEOF(tmpzfp%buffer_r)+SIZEOF(tmpzfp%buffer_i))
      ! deallocate(tmpzfp%buffer_r)
      ! deallocate(tmpzfp%buffer_i)




#if HAVE_ZFP
   if(present(ZFP_flag))then
      if(ZFP_flag==1)then
         allocate(data_buffer(size(t%core),1))
         data_buffer(:,1)=t%core
         call ZFP_Compress(data_buffer,t%coreZFP,size(t%core),1, eps_in/10,0)
         ! write(*,*)SIZEOF(t%core),(SIZEOF(t%coreZFP%buffer_r)+SIZEOF(t%coreZFP%buffer_i)),tmpsize
         deallocate(t%core)
      endif
   endif
#endif




   ! prefix sums
   t%psa(1)=1
   do idx=1,d
      t%psa(idx+1) = t%psa(idx) + t%n(idx)*t%r(idx)*t%r(idx+1)
   end do

end subroutine TT_Compress_SVD


integer function next_power_of_2(n)
  implicit none
  integer, intent(in) :: n
  integer :: power

  if (n <= 1) then
     power = 2 ! power is at least 2 instead of 1 for convenience
  else
     power = 1
     do while (power < n)
        power = power * 2
     end do
  end if
  next_power_of_2 = power
end function next_power_of_2


!>**** QTT decomposition via SVD of a full array (zero padding each mode size to power of 2)
!b: the d_org-dimensional input tensor, stored as 1D array
!tol: the relative compression tolerance
!t: the QTT decomposition. At input, t contains only d_org, n_org (and m_n_org)
!ZFP_flag: whether to further compress the core with zfp
subroutine QTT_Compress_SVD(b, tol, t, ZFP_flag)
   implicit none
   DT   :: b(:)
   real*8   :: tol
   type(TTtype) :: t
   integer,optional::ZFP_flag

   integer :: dd, i, ii,jj,N, d_org, iii,iii_new,s, tot_size,num_twos,idx
   integer, allocatable:: dims_org(:), dims(:), idx_MD_org(:), idx_MD(:),perm(:),shapes(:),shape_hd(:), Ns(:)
   DT, allocatable :: tensor_hd(:)

   allocate(Ns(t%d_org)) ! we enforce that for MPO, each dimension is padded to the same next power of 2
   if(t%mpo==1)then ! MPO
      do dd=1,t%d_org
         Ns(dd) = next_power_of_2(max(t%m_n_org(dd,1),t%m_n_org(dd,2)))
      enddo
      d_org = t%d_org*2
      allocate(tensor_hd(product(Ns)**2))
      tensor_hd = 0
      allocate(dims_org(d_org))
      dims_org(1:t%d_org) = t%m_n_org(1:t%d_org,1)
      dims_org(1+t%d_org:2*t%d_org) = t%m_n_org(1:t%d_org,2)
      allocate(dims(d_org))
      dims(1:t%d_org) = Ns
      dims(1+t%d_org:2*t%d_org) = Ns

   else ! MPS
      do dd=1,t%d_org
         Ns(dd) = next_power_of_2(t%n_org(dd))
      enddo
      d_org = t%d_org
      allocate(tensor_hd(product(Ns)))
      tensor_hd = 0
      allocate(dims_org(d_org))
      dims_org(1:t%d_org) = t%n_org(1:t%d_org)
      allocate(dims(d_org))
      dims = Ns
   endif


   allocate(idx_MD_org(d_org))
   do iii=1,product(dims_org)
      call SingleIndexToMultiIndex(d_org, dims_org, iii, idx_MD_org)
      call MultiIndexToSingleIndex(d_org, dims, iii_new, idx_MD_org)
      tensor_hd(iii_new) = b(iii)
   enddo
   deallocate(dims)
   deallocate(dims_org)
   deallocate(idx_MD_org)

   s = 0
   do dd=1,t%d_org
      s =s + int(log(real(Ns(dd)))/log(2d0))
   enddo
   num_twos = s


   if(t%mpo==1)then ! MPO

      allocate(shapes(num_twos))
      shapes = 4
      ! Build TT decomposition
      t%d = num_twos
      allocate(t%n(num_twos))
      t%n = shapes
      shapes = 2
      allocate(t%m_n(num_twos,2))
      t%m_n(:,1)=shapes
      t%m_n(:,2)=shapes
      call TT_Compress_SVD(tensor_hd, tol, t, ZFP_flag)
      deallocate(shapes)
      deallocate(tensor_hd)

   else ! MPS
      allocate(shapes(num_twos))
      shapes = 2
      ! Build TT decomposition
      t%d = num_twos
      allocate(t%n(num_twos))
      t%n = shapes
      call TT_Compress_SVD(tensor_hd, tol, t, ZFP_flag)
      deallocate(shapes)
      deallocate(tensor_hd)
   endif
   deallocate(Ns)

end subroutine QTT_Compress_SVD







subroutine TT_Decompress(t, b)
   implicit none
   ! Input TT tensor
   type(TTtype)  :: t
   ! Output full array
   DT, allocatable  :: b(:)
   DT, pointer ::data_buffer(:,:)
   integer :: d, i, prod_prev, prod_curr, total_size, idx
   DT, allocatable :: Y_in(:), Y(:)
   integer :: r_left, r_right, n_k
   integer, allocatable:: inv_perm(:),perm(:),shape_hd(:)
   real(kind=8)::tol_used

   d = t%d
   total_size = product(t%n)
   r_left  = t%r(1)         ! should be 1
   n_k     = t%n(1)
   r_right = t%r(2)
   prod_prev = n_k


#if HAVE_ZFP
   if(allocated(t%coreZFP%buffer_r))then
      allocate(data_buffer(t%psa(d+1)-1, 1))
      call ZFP_Decompress(data_buffer,t%coreZFP, t%psa(d+1)-1, 1, tol_used, 1)
      allocate(t%core(t%psa(d+1)-1))
      t%core = data_buffer(:,1)
      deallocate(data_buffer)
   endif
#endif


   allocate(Y(n_k * r_right))
   Y = t%core(t%psa(1) : t%psa(2)-1)

   do i = 2, d
      n_k = t%n(i)
      r_left = t%r(i)
      r_right = t%r(i+1)
      prod_curr = prod_prev * n_k

      allocate(Y_in(prod_prev*r_left))
      Y_in=Y
      deallocate(Y)
      allocate(Y(prod_prev * n_k * r_right))

      call gemmf77('N', 'N', prod_prev, n_k * r_right, r_left, BPACK_cone, Y_in, prod_prev, t%core(t%psa(i):t%psa(i+1)-1), r_left, BPACK_czero, Y, prod_prev)

      prod_prev = prod_curr

      deallocate(Y_in)
   end do


   if(t%mpo==1)then ! MPO
      allocate(shape_hd(2*d))
      shape_hd(1:d) = t%m_n(:,1)
      shape_hd(1+d:2*d) = t%m_n(:,2)
      allocate(perm(2*d))
      perm = 0
      do idx=1,d
        perm(2*(idx-1)+1 ) = idx            ! row bits
        perm(2*(idx-1)+2 ) = d + idx        ! col bits
      enddo

      allocate(inv_perm(2*d))
      inv_perm = 0
      do idx = 1, 2*d
         inv_perm( perm(idx) ) = idx
      end do
      call TensorPermute(Y, shape_hd, inv_perm)
      deallocate(perm)
      deallocate(inv_perm)
      deallocate(shape_hd)
   endif


   allocate(b(total_size))
   b = Y
   deallocate(Y)

#if HAVE_ZFP
   if(allocated(t%coreZFP%buffer_r))then
      deallocate(t%core)
   endif
#endif


end subroutine TT_Decompress



subroutine QTT_Decompress(t, b)
   implicit none
   ! Input QTT tensor
   type(TTtype)  :: t
   ! Output full array
   DT, allocatable  :: b(:), tensor_hd(:)
   integer :: d, i, ii,jj,N, d_org, iii,iii_new,s,dd, tot_size,num_twos,idx
   integer, allocatable:: dims_org(:), dims(:), idx_MD_org(:), idx_MD(:),perm(:),inv_perm(:),shapes(:),shape_hd(:),Ns(:)


   call TT_Decompress(t,tensor_hd)

   allocate(Ns(t%d_org)) ! we enforce that for MPO, each dimension is padded to the same next power of 2

   if(t%mpo==1)then ! MPO
      do dd=1,t%d_org
         Ns(dd) = next_power_of_2(max(t%m_n_org(dd,1),t%m_n_org(dd,2)))
      enddo
      allocate(b(product(t%m_n_org(:,1))*product(t%m_n_org(:,2))))
      b = 0
      d_org = t%d_org*2
      allocate(dims_org(d_org))
      dims_org(1:t%d_org) = t%m_n_org(1:t%d_org,1)
      dims_org(1+t%d_org:2*t%d_org) = t%m_n_org(1:t%d_org,2)
      allocate(dims(d_org))
      dims(1:t%d_org) = Ns
      dims(1+t%d_org:2*t%d_org) = Ns
   else ! MPS
      do dd=1,t%d_org
         Ns(dd) = next_power_of_2(t%n_org(dd))
      enddo
      allocate(b(product(t%n_org)))
      b = 0
      d_org = t%d_org
      allocate(dims_org(d_org))
      dims_org(1:t%d_org) = t%n_org(1:t%d_org)
      allocate(dims(d_org))
      dims = Ns
   endif

   allocate(idx_MD_org(d_org))
   do iii=1,product(dims_org)
      call SingleIndexToMultiIndex(d_org, dims_org, iii, idx_MD_org)
      call MultiIndexToSingleIndex(d_org, dims, iii_new, idx_MD_org)
      b(iii) = tensor_hd(iii_new)
   enddo
   deallocate(dims)
   deallocate(dims_org)
   deallocate(idx_MD_org)
   deallocate(Ns)

end subroutine QTT_Decompress



subroutine TT_Delete(t,keep_zfpstream)
   type(TTtype)  :: t
   integer,optional :: keep_zfpstream

   if(allocated(t%n_org))deallocate(t%n_org)
   if(allocated(t%n))deallocate(t%n)
   if(allocated(t%m_n))deallocate(t%m_n)
   if(allocated(t%m_n_org))deallocate(t%m_n_org)
   if(allocated(t%r))deallocate(t%r)
   if(allocated(t%psa))deallocate(t%psa)
   if(allocated(t%core))deallocate(t%core)


#if HAVE_ZFP
   if (allocated(t%coreZFP%buffer_r)) then
      deallocate (t%coreZFP%buffer_r)
      if(.not. present(keep_zfpstream))call zFORp_stream_close(t%coreZFP%stream_r)
   endif
   if (allocated(t%coreZFP%buffer_i))then
      deallocate (t%coreZFP%buffer_i)
      if(.not. present(keep_zfpstream))call zFORp_stream_close(t%coreZFP%stream_i)
   endif
#endif



end subroutine TT_Delete



subroutine TT_Copy(t,t_new)
   type(TTtype)  :: t,t_new
   integer :: mm

   t_new%mpo=t%mpo
   t_new%d=t%d
   t_new%d_org=t%d_org

   if(allocated(t%n_org))then
      allocate(t_new%n_org(size(t%n_org)))
      t_new%n_org = t%n_org
   endif
   if(allocated(t%n))then
      allocate(t_new%n(size(t%n)))
      t_new%n = t%n
   endif
   if(allocated(t%m_n))then
      allocate(t_new%m_n(size(t%m_n,1),size(t%m_n,2)))
      t_new%m_n = t%m_n
   endif
   if(allocated(t%m_n_org))then
      allocate(t_new%m_n_org(size(t%m_n_org,1),size(t%m_n_org,2)))
      t_new%m_n_org = t%m_n_org
   endif
   if(allocated(t%r))then
      allocate(t_new%r(size(t%r)))
      t_new%r = t%r
   endif
   if(allocated(t%psa))then
      allocate(t_new%psa(size(t%psa)))
      t_new%psa = t%psa
   endif
   if(allocated(t%core))then
      allocate(t_new%core(size(t%core)))
      t_new%core = t%core
   endif


#if HAVE_ZFP
   if(allocated(t%coreZFP%buffer_r))then
      mm = size(t%coreZFP%buffer_r,1)
      allocate(t_new%coreZFP%buffer_r(mm))
      t_new%coreZFP%buffer_r = t%coreZFP%buffer_r
      t_new%coreZFP%stream_r = t%coreZFP%stream_r
   endif
   if(allocated(t%coreZFP%buffer_i))then
      mm = size(t%coreZFP%buffer_i,1)
      allocate(t_new%coreZFP%buffer_i(mm))
      t_new%coreZFP%buffer_i = t%coreZFP%buffer_i
      t_new%coreZFP%stream_i = t%coreZFP%stream_i
   endif
#endif


end subroutine TT_Copy


!>**** TT application/contraction assuming the input vector is a full array
!tmat: the MPO
!b: the input array
!c: the output array
subroutine TT_Apply_Fullvec(tmat, b, c)
   implicit none
   type(TTtype)  :: tmat
   DT        :: b(:)
   DT       :: c(:)

   integer :: d, k
   integer, allocatable :: psa(:), r(:)
   DT, allocatable :: cra(:)
   DT, allocatable :: ctmp(:)
   integer :: rb
   integer :: rk, rk1, n_k, m_k, nvec
   integer :: crSize, restDim
   integer :: irow, icol, i3, shape_tmp(3), perm_tmp(3), shape_tmp4(4), perm_tmp4(4)
   DT, allocatable :: cr(:), cc(:)
   DT, pointer::data_buffer(:,:)
   real(kind=8)::tol_used

   d  = tmat%d

#if HAVE_ZFP
   if(allocated(tmat%coreZFP%buffer_r))then
      allocate(data_buffer(tmat%psa(d+1)-1, 1))
      call ZFP_Decompress(data_buffer,tmat%coreZFP, tmat%psa(d+1)-1, 1, tol_used, 1)
      allocate(tmat%core(tmat%psa(d+1)-1))
      tmat%core = data_buffer(:,1)
      deallocate(data_buffer)
   endif
#endif
   nvec = size(b)/product(tmat%m_n(:,2))
   allocate(cc(size(b)))
   cc = b

   do k=1,d
      rk = tmat%r(k)
      rk1 = tmat%r(k+1)
      m_k  = tmat%m_n(k,1)
      n_k  = tmat%m_n(k,2)

      crSize = tmat%psa(k+1)-tmat%psa(k)
      allocate(cr(crSize))
      cr = tmat%core(tmat%psa(k):tmat%psa(k+1)-1)

      shape_tmp4(1)=rk
      shape_tmp4(2)=m_k
      shape_tmp4(3)=n_k
      shape_tmp4(4)=rk1
      perm_tmp4(1) =2
      perm_tmp4(2) =4
      perm_tmp4(3) =1
      perm_tmp4(4) =3
      ! [rk,m_k,n_k,rk1] -> [m_k,rk1,rk,n_k]
      call TensorPermute(cr, shape_tmp4, perm_tmp4)

      restDim = size(cc)/(rk*n_k)
      allocate(ctmp(m_k*rk1*restDim))
      ctmp=0
      call gemmf77('N', 'N', m_k*rk1, restDim, n_k*rk, BPACK_cone, cr, m_k*rk1, cc, rk*n_k, BPACK_czero, ctmp, m_k*rk1)


      deallocate(cc)
      allocate(cc(m_k*rk1*restDim))

      shape_tmp(1)=m_k
      shape_tmp(2)=rk1*restDim/nvec
      shape_tmp(3)=nvec
      perm_tmp(1) =2
      perm_tmp(2) =1
      perm_tmp(3) =3
      ! [m_k,rk1*restDim/nvec,nvec] -> [rk1*restDim/nvec,m_k,nvec]
      call TensorPermute(ctmp, shape_tmp, perm_tmp)

      cc = ctmp
      deallocate(ctmp, cr)
   end do

   c = cc
   deallocate(cc)

#if HAVE_ZFP
   if(allocated(tmat%coreZFP%buffer_r))then
      deallocate(tmat%core)
   endif
#endif

end subroutine TT_Apply_Fullvec




!>**** QTT application/contraction assuming the input vector is a full array (the input and output will be possibly zero padded)
!tmat: the MPO
!b: the input array
!c: the output array
subroutine QTT_Apply_Fullvec(tmat, b, c)
   implicit none
   type(TTtype)  :: tmat
   DT        :: b(:)
   DT       :: c(:)
   DT,allocatable :: c_tmp(:), b_tmp(:)
   integer :: d,dd,d_org,iii,iii_new,nvec,nelem,nelem_new,vv
   integer,allocatable:: dims_org(:),dims(:),idx_MD_org(:), Ns(:)

   allocate(Ns(tmat%d_org)) ! we enforce that for MPO, each dimension is padded to the same next power of 2

   d_org = tmat%d_org
   do dd=1,tmat%d_org
      Ns(dd) = next_power_of_2(max(tmat%m_n_org(dd,1),tmat%m_n_org(dd,2)))
   enddo

   allocate(dims_org(d_org))
   dims_org(1:d_org) = tmat%m_n_org(1:d_org,2)
   nelem = product(dims_org)
   nvec = size(b)/nelem
   nelem_new = product(Ns)

   allocate(b_tmp(product(Ns)*nvec))
   b_tmp = 0

   allocate(dims(d_org))
   dims = Ns
   allocate(idx_MD_org(d_org))

   do vv=1,nvec
   do iii=1,product(dims_org)
      call SingleIndexToMultiIndex(d_org, dims_org, iii, idx_MD_org)
      call MultiIndexToSingleIndex(d_org, dims, iii_new, idx_MD_org)
      b_tmp(iii_new+(vv-1)*nelem_new) = b(iii+(vv-1)*nelem)
   enddo
   enddo

   allocate(c_tmp(product(Ns)*nvec))
   c_tmp = 0
   call TT_Apply_Fullvec(tmat, b_tmp, c_tmp)

   dims_org(1:d_org) = tmat%m_n_org(1:d_org,1)
   nelem = product(dims_org)
   nelem_new = product(Ns)

   do vv=1,nvec
   do iii=1,product(dims_org)
      call SingleIndexToMultiIndex(d_org, dims_org, iii, idx_MD_org)
      call MultiIndexToSingleIndex(d_org, dims, iii_new, idx_MD_org)
      c(iii+(vv-1)*nelem) = c_tmp(iii_new+(vv-1)*nelem_new)
   enddo
   enddo

   deallocate(dims)
   deallocate(dims_org)
   deallocate(idx_MD_org)
   deallocate(b_tmp)
   deallocate(c_tmp)
   deallocate(Ns)

end subroutine QTT_Apply_Fullvec


!>**** TT application/contraction assuming the input vector is also in TT format
!tmat: the MPO
!b: the input MPS
!c: the output MPS
subroutine TT_Apply_TTvec(tmat, b, c)
   implicit none
   type(TTtype)  :: tmat
   type(TTtype) :: b
   type(TTtype) :: c
   integer :: rbk, rbk1, rak, rak1, n_k, m_k,shape_tmp(2), perm_tmp(2), shape_tmp5(5), perm_tmp5(5), shape_tmp3(3), perm_tmp3(3)
   DT, allocatable :: cr(:), cc(:), ctmp(:)

   integer :: d, i
   integer, allocatable :: psp(:)



   d = tmat%d
   c%d = d
   allocate(c%n(d))
   c%n = tmat%m_n(:,2)

   allocate(c%r(d+1))
   do i=1,d+1
      c%r(i) = tmat%r(i)*b%r(i)
   end do


   allocate(c%psa(d+1))
   c%psa(1) = 1
   do i=1,d
      c%psa(i+1) = c%psa(i) + c%n(i)*c%r(i)*c%r(i+1)
   end do

   allocate(c%core(c%psa(d+1)-1))
   c%core = 0d0

   do i=1,d
      m_k  = tmat%m_n(i,1)
      n_k  = tmat%m_n(i,2)
      allocate(cr(tmat%psa(i+1)-tmat%psa(i)))
      allocate(cc(b%psa(i+1)-b%psa(i)))
      cr=tmat%core(tmat%psa(i):tmat%psa(i+1)-1)
      cc=b%core(b%psa(i):b%psa(i+1)-1)
      rak = tmat%r(i)
      rak1 = tmat%r(i+1)
      rbk = b%r(i)
      rbk1 = b%r(i+1)

      shape_tmp(1)=rak*m_k*n_k
      shape_tmp(2)=rak1
      perm_tmp(1) =2
      perm_tmp(2) =1
      ! [rak,m_k,n_k,rak1] -> [rak1,rak,m_k,n_k]
      call TensorPermute(cr, shape_tmp, perm_tmp)

      shape_tmp3(1)=rbk
      shape_tmp3(2)=n_k
      shape_tmp3(3)=rbk1
      perm_tmp3(1) =2
      perm_tmp3(2) =1
      perm_tmp3(3) =3
      ! [rbk,n_k,rbk1] -> [n_k,rbk,rbk1]
      call TensorPermute(cc, shape_tmp3, perm_tmp3)

      allocate(ctmp(rak1*rak*m_k*rbk*rbk1))
      ctmp=0
      call gemmf77('N', 'N', rak1*rak*m_k, rbk*rbk1, n_k, BPACK_cone, cr, rak1*rak*m_k, cc, n_k, BPACK_czero, ctmp, rak1*rak*m_k)


      shape_tmp5(1)=rak1
      shape_tmp5(2)=rak
      shape_tmp5(3)=m_k
      shape_tmp5(4)=rbk
      shape_tmp5(5)=rbk1
      perm_tmp5(1) =2
      perm_tmp5(2) =4
      perm_tmp5(3) =3
      perm_tmp5(4) =1
      perm_tmp5(5) =5
      ! [rak1,rak,m_k,rbk,rbk1] -> [rak,rbk,m_k,rak1,rbk1]
      call TensorPermute(ctmp, shape_tmp5, perm_tmp5)
      c%core(c%psa(i):c%psa(i+1)-1) = ctmp

      deallocate(ctmp,cc,cr)
   enddo

end subroutine TT_Apply_TTvec




end module MISC_Utilities
