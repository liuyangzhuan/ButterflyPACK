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
!> @file
!> @brief This is an example that solves a 3D EFIE/MFIE/CFIE system for electromagnetics scattering.
!> @details Note that instead of the use of precision dependent subroutine/module/type names "z_", one can also use the following \n
!> #define DAT 0 \n
!> #include "zButterflyPACK_config.fi" \n
!> which will macro replace precision-independent subroutine/module/type names "X" with "z_X" defined in SRC_DOUBLECOMLEX with double-complex precision

#define DAT 0
#include "zButterflyPACK_config.fi"

program test_proxy_selection

   use BPACK_DEFS
	use MISC_Utilities
	use MISC_DenseLA


implicit none


type(TreeNode), pointer :: root => null()

integer :: a, b, m_pick, k, N, nodeID, mid,i,j
integer, allocatable :: indices_pick(:)
integer :: count,ii,jj,ranknew,mn

complex(kind=8),allocatable:: matA(:,:), matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:),matA_recon(:,:), core(:,:), core_row(:,:), mats_interp(:,:), mats_skel(:,:), mattmp1(:,:), mattmp2(:,:), mattmp3(:,:),UU(:,:),VV(:,:)

complex(kind=8),allocatable::mats(:,:)
complex(kind=8),allocatable::Qs(:,:),tau(:),tau_row(:)
real(kind=8),allocatable::locs(:,:),locsp(:,:), dists(:), Singular(:)

integer:: nrow0, ncol0, rankmax_c, rankmax_r, rankmax_r1, nx,ny, nxp,nyp,nrow,ncol,nq
real(kind=8)::x,xp,dx,dxp,aa,tol_comp,tol_Rdetect,tol_running, dy,sample_para,matnew_norm

integer, allocatable::jpvt(:),jpvt_row(:), nns(:,:), order(:), select_row(:)
real(kind=8)::n1, n2, dist, wavenum, norm_tol=1.0d-30, tol_next, maxerror, normval
integer leaf, dims(2), inds(2), ndim, edge_n, jjj, knn, len, H
type(mesh)::msh
type(treesamplequant)::treequant(1)

nrow0=1000000
ncol0=100
m_pick=4
leaf=64
tol_comp=1e-4
tol_Rdetect=tol_comp   !1e-10
ndim=1
knn=10
sample_para=2d0





if(ndim==0)then
   ! ---- Read from file ----
   open(unit=11, file='matrix.dat', form='unformatted', status='old')
   read(11) nrow,ncol
   write(*,*)"matrix read from the file:", nrow,ncol
   allocate(matA(nrow,ncol))
   read(11) matA
   close(11)
else if(ndim==1)then
   nrow=nrow0
   ncol=ncol0
   allocate(matA(nrow,ncol))
   matA = 0
   dx = 1d0/(nrow+1)
   dxp = 1d0/(ncol+1)

   allocate(nns(ncol,knn))
   allocate(dists(nrow))
   allocate(order(nrow))
   nns=0
   do jj=1,ncol
      do ii=1,nrow
         x = dx*(ii-1)
         xp = dx*(jj-1)+1d0
         dists(ii) = abs(x-xp)
      enddo
      call quick_sort(dists, order, nrow)
      nns(jj, 1:knn) = order(1:knn)
   enddo

   do ii=1,nrow
   do jj=1,ncol
      x = dx*(ii-1)
      xp = dx*(jj-1)+1d0
      matA(ii,jj) = 1/abs(x-xp)
      ! matA(ii,jj) = 1
      ! call random_number(aa)
      ! matA(ii,jj)=aa
   enddo
   enddo
elseif(ndim==2)then

   call as_square_as_possible(nrow0,nx,ny)
   call as_square_as_possible(ncol0,nxp,nyp)
   nrow=nx*ny
   ncol=nxp*nyp
   allocate(matA(nrow,ncol))
   matA = 0

   dx = 1d0/(nx+1)
   dy = 1d0/(ny+1)

   wavenum = 2*BPACK_junit/(dx*10)

   allocate(locs(2,nrow))
   dims(1)=nx
   dims(2)=ny
   do k=1,nrow
   call SingleIndexToMultiIndex(2,dims, k, inds)
   ii=inds(1)
   jj=inds(2)
   locs(1,k)=ii*dx-dx
   locs(2,k)=jj*dy
   enddo
   allocate(locsp(2,ncol))
   do n=1,ncol
   call SingleIndexToMultiIndex(2,dims, n, inds)
   ii=inds(1)
   jj=inds(2)
   locsp(1,n)=ii*dx + 1
   locsp(2,n)=jj*dy
   enddo

   msh%Nunk = nrow
   allocate (msh%xyz(2, 1:msh%Nunk))
   msh%xyz = locs
   call KD_tree_ordering(leaf, msh)
   locs(:,msh%old2new) = locs(:,:)


   allocate(nns(ncol,knn))
   allocate(dists(nrow))
   allocate(order(nrow))
   nns=0
   do jj=1,ncol
      do ii=1,nrow
         dists(ii) = sqrt(sum((locs(:,ii) - locsp(:,jj))**2d0))
      enddo
      call quick_sort(dists, order, nrow)
      nns(jj, 1:knn) = order(1:knn)
   enddo

   do ii=1,nrow
   do jj=1,ncol
      dist = sqrt(sum((locs(:,ii) - locsp(:,jj))**2d0))
      ! matA(ii,jj) = 1/dist
      matA(ii,jj) = exp(BPACK_junit*dist*wavenum)/dist
      ! matA(ii,jj) = 1
      ! call random_number(aa)
      ! matA(ii,jj)=aa
   enddo
   enddo

endif





!!!!!!!!!!!!!!! The following is tree-sampling-based


   !  type treesamplequant
   !      integer M,N !< matrix dimensions
   !      integer:: header_m=0, header_n=0 !< matrix offsets
   !      ! integer::leaf !< the level of the AVL tree below which m_pick will be reduced
   !      ! integer:: m_pick !< number of new proxies to pick per intervel
   !      integer:: keeprefine !< whether to generate more samples
   !      type(TreeNode), pointer :: root_node !< root node of the AVL tree
   !      type(frame_t), allocatable :: stack(:) !< the stack that keeps track of the tree traversal
   !      integer:: top=0 !< the current iterator for the stack
   !      integer, allocatable:: skel_rows(:), skel_cols(:) !< rows or columns to be skeletonized
   !      integer, allocatable:: proxies(:) !< rows or columns that have been sampled so far
   !      DT, allocatable:: Qs(:,:), mats(:,:) !< basis matrices and original matrices at the current stage
   !  end type treesamplequant


   treequant(1)%header_m=1
   ! treequant(1)%header_n=1
   treequant(1)%M=nrow
   treequant(1)%N=ncol
   ! treequant(1)%leaf=leaf
   treequant(1)%top=1
   treequant(1)%tail=1

   len = max(0, nrow)
   H=int(log(dble(len))/log(2d0))+1
   allocate(treequant(1)%stack(H))
   ! Push root frame
   treequant(1)%stack(1)%lo      = 1
   treequant(1)%stack(1)%hi      = nrow
   treequant(1)%stack(1)%nodeID   = 1
   treequant(1)%stack(1)%m_pick  = m_pick
   treequant(1)%stack(1)%tol_run = 1e30 ! not used for the root node

   allocate(treequant(1)%skel_cols(ncol))
   do j=1,ncol
      treequant(1)%skel_cols(j) = j
   enddo




n1 = MPI_Wtime()
call refine_interval(matA, treequant, 1, tol_Rdetect,leaf)

! write(*,*)treequant(1)%proxies

n2 = MPI_Wtime()
write(*,*)'time for finding the approxies: ', n2-n1






rankmax_r=size(treequant(1)%mats,1)
rankmax_c=size(treequant(1)%mats,2)


! do ii=1,rankmax_r
!    normval=fnorm(treequant(1)%mats(ii:ii,:),1,rankmax_c)
!    if(normval>1e-30)then
!       treequant(1)%mats(ii:ii,:) = treequant(1)%mats(ii:ii,:)/normval
!    endif
! enddo

! do jj=1,rankmax_c
!    normval=fnorm(treequant(1)%mats(:,jj:jj),rankmax_r,1)
!    if(normval>1e-30)then
!       treequant(1)%mats(:,jj:jj) = treequant(1)%mats(:,jj:jj)/normval
!    endif
! enddo

allocate (core(rankmax_r, rankmax_c))
core = treequant(1)%mats
allocate(mats(rankmax_r,rankmax_c))
mats = treequant(1)%mats

allocate (jpvt(max(rankmax_c, rankmax_r)))
allocate (tau(max(rankmax_c, rankmax_r)))
jpvt = 0



call geqp3modf90(core, jpvt, tau, tol_comp*0.5, BPACK_SafeUnderflow, ranknew)


! ! SVD
! mn=min(rankmax_r,rankmax_c)
! allocate(UU(rankmax_r,mn))
! allocate(VV(mn,rankmax_c))
! allocate(Singular(mn))
! UU = 0
! VV = 0
! Singular = 0
! call SVD_Truncate(core, rankmax_r, rankmax_c, mn, UU, VV, Singular, tol_comp*0.01, norm_tol, ranknew)
! deallocate(UU)
! deallocate(VV)
! deallocate(Singular)
! write(*,*)"rank of treequant(1)%mats",ranknew



! allocate(core_row(rankmax_c,rankmax_r))
! allocate (jpvt_row(max(rankmax_c, rankmax_r)))
! allocate (tau_row(max(rankmax_c, rankmax_r)))
! jpvt_row = 0

! call copymatT(treequant(1)%mats,core_row,rankmax_r,rankmax_c)
! call geqp3modf90(core_row, jpvt_row, tau_row, tol_comp*0.1, BPACK_SafeUnderflow, ranknew)
! deallocate(mats)
! allocate(mats(ranknew,rankmax_c))
! mats = treequant(1)%mats(jpvt_row(1:ranknew),:)

! deallocate(core)
! allocate(core(ranknew,rankmax_c))
! core = mats
! call geqp3f90(core, jpvt, tau)






if (ranknew > 0) then
   call un_or_mqrf90(core, tau, mats, 'L', 'C', size(mats,1), ncol, ranknew)
   call trsmf90(core, mats, 'L', 'U', 'N', 'N', ranknew, ncol)
else
   ranknew = 1
   jpvt(1) = 1
   mats = 0
endif

allocate(mats_interp(ranknew,ncol))
mats_interp =mats(1:ranknew, 1:ncol)
allocate(mats_skel(nrow,ranknew))
mats_skel = matA(:,jpvt(1:ranknew))

allocate(matA_recon(nrow,ncol))
matA_recon=0

call gemmf90(mats_skel, nrow, mats_interp, ranknew, matA_recon, nrow, 'N', 'N', nrow, ncol, ranknew, BPACK_cone,BPACK_czero)


         maxerror=0
         do ii=1,nrow
            maxerror = max(maxerror,fnorm(matA_recon(ii:ii,:)-matA(ii:ii,:),1,ncol)/fnorm(matA(ii:ii,:), 1, ncol))
            ! write(*,*)ii,fnorm(matA_recon(ii:ii,:)-matA(ii:ii,:),1,ncol)/fnorm(matA(ii:ii,:), 1, ncol),fnorm(matA(ii:ii,:), 1, ncol)
         enddo


write(*,*)'size of mats (tree-sampling-based):', size(treequant(1)%mats,1),size(treequant(1)%mats,2)
write(*,*)'reconstruction error (for fullmat):',fnorm(matA_recon-matA,nrow,ncol)/fnorm(matA,nrow,ncol), maxerror,'rank:',ranknew

! write(*,*)jpvt(1:ranknew)


call matrix_resize(mats_skel,rankmax_r,ranknew)
mats_skel = treequant(1)%mats(:,jpvt(1:ranknew))
call matrix_resize(matA_recon,rankmax_r,ncol)
matA_recon=0
call gemmf90(mats_skel, rankmax_r, mats_interp, ranknew, matA_recon, rankmax_r, 'N', 'N', rankmax_r, ncol, ranknew, BPACK_cone,BPACK_czero)
write(*,*)'reconstruction error (for treequant(1)%mats):',fnorm(matA_recon-treequant(1)%mats,rankmax_r,ncol)/fnorm(treequant(1)%mats,rankmax_r,ncol), 'rank:',ranknew





         ! count = nrow
         ! nq = size(treequant(1)%Qs, 2)
         ! allocate(mattmp1(nq, count)); mattmp1 = 0
         ! allocate(mattmp2(ncol, count)); mattmp2 = 0

         ! ! Project matnew onto span(Qs) and compute residual
         ! call gemmf90(treequant(1)%Qs, ncol, matA, count, mattmp1, nq, 'C', 'T', nq, count, ncol, BPACK_cone, BPACK_czero)
         ! call gemmf90(treequant(1)%Qs, ncol, mattmp1, nq, mattmp2, ncol, 'N', 'N', ncol, count, nq, BPACK_cone, BPACK_czero)

         ! allocate(mattmp3(ncol, count))
         ! call copymatT(matA, mattmp3, count, ncol)
         ! matnew_norm = fnorm(mattmp3, ncol, count)

         ! mattmp3 = mattmp3 - mattmp2

         ! ! Update tolerance passed to children
         ! tol_next = min(tol_next, max(matnew_norm * tol_Rdetect, norm_tol))

         ! ! Re-project residual to improve numerical stability
         ! call gemmf90(treequant(1)%Qs, ncol, mattmp3, ncol, mattmp1, nq, 'C', 'N', nq, count, ncol, BPACK_cone, BPACK_czero)
         ! call gemmf90(treequant(1)%Qs, ncol, mattmp1, nq, mattmp2, ncol, 'N', 'N', ncol, count, nq, BPACK_cone, BPACK_czero)
         ! ! write(*,*)'node ', nodeID, 'count', count, 'norm after and before the projection', fnorm(mattmp3- mattmp2,ncol,count),matnew_norm
         ! mattmp3 = mattmp3 - mattmp2
         ! maxerror=0
         ! do ii=1,nrow
         !    maxerror = max(maxerror,fnorm(mattmp3(:,ii:ii),ncol,1)/fnorm(matA(ii:ii,:), 1, ncol))
         ! enddo

         ! write(*,*)'error check ', size(treequant(1)%proxies,1), 'fnorm and max error',fnorm(mattmp3, ncol, count),maxerror, 'tol',tol_next,matnew_norm

         ! deallocate(mattmp1)
         ! deallocate(mattmp2)
         ! deallocate(mattmp3)







deallocate(mats_interp)
deallocate(mats_skel)
deallocate(core)
deallocate(jpvt)
deallocate(tau)
deallocate(treequant(1)%mats)
deallocate(matA_recon)


!!!!!!!!!!!!!!! The following is knn-based

if(ndim>0)then

rankmax_r1 = min(ceiling_safe(sample_para*ranknew*2), nrow)
rankmax_c = ncol
allocate (select_row(rankmax_r1 + ncol*knn))
call linspaceI(1, nrow, rankmax_r1, select_row(1:rankmax_r1))
do j = 1, ncol
   edge_n=j
   do jjj = 1, knn
      rankmax_r1 = rankmax_r1 + 1
      select_row(rankmax_r1) = nns(edge_n, jjj)
   enddo
enddo
call remove_dup_int(select_row, rankmax_r1, rankmax_r)


allocate (core(rankmax_r, rankmax_c))
call matrix_resize(mats,rankmax_r, rankmax_c)
mats = matA(select_row(1:rankmax_r),:)
core = mats


allocate (jpvt(max(rankmax_c, rankmax_r)))
allocate (tau(max(rankmax_c, rankmax_r)))
jpvt = 0
call geqp3modf90(core, jpvt, tau, tol_comp, BPACK_SafeUnderflow, ranknew)

if (ranknew > 0) then
   call un_or_mqrf90(core, tau, mats, 'L', 'C', rankmax_r, ncol, ranknew)
   call trsmf90(core, mats, 'L', 'U', 'N', 'N', ranknew, ncol)
else
   ranknew = 1
   jpvt(1) = 1
   mats = 0
endif

allocate(mats_interp(ranknew,ncol))
mats_interp =mats(1:ranknew, 1:ncol)
allocate(mats_skel(nrow,ranknew))
mats_skel = matA(:,jpvt(1:ranknew))

allocate(matA_recon(nrow,ncol))
matA_recon=0

call gemmf90(mats_skel, nrow, mats_interp, ranknew, matA_recon, nrow, 'N', 'N', nrow, ncol, ranknew, BPACK_cone,BPACK_czero)

write(*,*)'size of mats (knn-based):', size(mats,1),size(mats,2)
write(*,*)'reconstruction error:',fnorm(matA_recon-matA,nrow,ncol)/fnorm(matA,nrow,ncol), 'rank:',ranknew

deallocate(mats_interp)
deallocate(mats_skel)
deallocate(core)
deallocate(jpvt)
deallocate(tau)
deallocate(mats)
deallocate(matA_recon)
endif


!!!!!!!!!!!!!!! The following computes the true rank
call ComputeRange(nrow,ncol, matA, ranknew, 1, tol_comp, norm_tol=1d-13)
write(*,*)'rank of matA:', ranknew


call free_tree(root)

contains


   subroutine refine_interval(matZ, treequant, index_ij, tol_comp, leaf)
   implicit none
   real(kind=8), intent(in) :: tol_comp
   integer, intent(in) :: leaf
   integer, allocatable :: indices_pick(:)
   integer :: count
   integer :: index_ij
   type(treesamplequant)::treequant(:)
   complex(kind=8), intent(in) :: matZ(:,:)

   ! Local work arrays
   complex(kind=8), allocatable :: matnew(:,:), mattmp1(:,:), mattmp2(:,:), mattmp3(:,:)
   integer :: nrow, ncol, nq, ranknew, rank
   integer :: mid, avail, m_adjust
   real(kind=8) :: matnew_norm, norm_tol=1.0d-30
   real(kind=8) :: tol_next, residual
   integer :: nodeID, depth, lo, hi, m_pick, keeprefine, top



   do while (treequant(index_ij)%top <=treequant(index_ij)%tail)
   ! do while (treequant(index_ij)%top >0)

   ! Pop next frame
   nrow = treequant(index_ij)%M
   ncol = treequant(index_ij)%N
   lo        = treequant(index_ij)%stack(treequant(index_ij)%top)%lo
   hi        = treequant(index_ij)%stack(treequant(index_ij)%top)%hi
   m_pick    = treequant(index_ij)%stack(treequant(index_ij)%top)%m_pick
   nodeID     = treequant(index_ij)%stack(treequant(index_ij)%top)%nodeID

   ! Base conditions
   if (lo > hi) cycle
   avail = unused_in_interval(treequant(index_ij)%root_node, 1, nrow, lo, hi)
   ! write(*,*)nodeID,lo,hi,m_pick,avail
   if (avail <= 0) then
      treequant(index_ij)%top = treequant(index_ij)%top + 1
      cycle
   endif
   ! Pick up to m_pick unique numbers in [lo, hi]
   allocate(indices_pick(m_pick))

   call pick_exact_m_interval(treequant(index_ij)%root_node, 1, nrow, lo, hi, m_pick, indices_pick, count)





   ! Sampling step
   allocate(matnew(count, ncol))
   matnew = matZ(indices_pick(1:count), :)



   nrow = treequant(index_ij)%M
   ncol = treequant(index_ij)%N
   lo        = treequant(index_ij)%stack(treequant(index_ij)%top)%lo
   hi        = treequant(index_ij)%stack(treequant(index_ij)%top)%hi
   nodeID     = treequant(index_ij)%stack(treequant(index_ij)%top)%nodeID
   depth     = floor_safe(log(dble(nodeID))/log(2d0))
   m_pick    = treequant(index_ij)%stack(treequant(index_ij)%top)%m_pick
   tol_next  = treequant(index_ij)%stack(treequant(index_ij)%top)%tol_run
   m_adjust = treequant(index_ij)%stack(treequant(index_ij)%top)%m_pick
   ! if (leaf > hi - lo + 1) m_adjust = 1

   treequant(index_ij)%top = treequant(index_ij)%top + 1
   if (depth == 0) then
      allocate(treequant(index_ij)%mats(size(matnew,1), size(matnew,2)))
      treequant(index_ij)%mats = matnew
      allocate(treequant(index_ij)%proxies(count))
      treequant(index_ij)%proxies=indices_pick(1:count)
   else
      call matrix_resize(treequant(index_ij)%mats, size(treequant(index_ij)%mats,1) + size(matnew,1), ncol)
      treequant(index_ij)%mats(size(treequant(index_ij)%mats,1) - size(matnew,1) + 1 : size(treequant(index_ij)%mats,1), :) = matnew
      call array_resize_int(treequant(index_ij)%proxies, size(treequant(index_ij)%proxies,1) + count)
      treequant(index_ij)%proxies(size(treequant(index_ij)%proxies,1) - count + 1 : size(treequant(index_ij)%proxies,1)) = indices_pick(1:count)
   end if
   deallocate(indices_pick)


   keeprefine = 1
   if (count == hi - lo + 1) then
      keeprefine = 0

   else if (depth == 0) then
      keeprefine = 1


      allocate (core(count, ncol))
      core = treequant(index_ij)%mats
      allocate(mats(count,ncol))
      mats = treequant(index_ij)%mats

      allocate (jpvt(max(ncol, count)))
      allocate (tau(max(ncol, count)))
      jpvt = 0
      call geqp3modf90(core, jpvt, tau, tol_comp, BPACK_SafeUnderflow, ranknew)
      if (ranknew > 0) then
         call un_or_mqrf90(core, tau, mats, 'L', 'C', size(mats,1), ncol, ranknew)
         call trsmf90(core, mats, 'L', 'U', 'N', 'N', ranknew, ncol)
      else
         ranknew = 1
         jpvt(1) = 1
         mats = 0
      endif

      allocate(treequant(index_ij)%mats_interp(ranknew,ncol))
      treequant(index_ij)%mats_interp =mats(1:ranknew, 1:ncol)
      allocate(treequant(index_ij)%skel_sofar(ranknew))
      treequant(index_ij)%skel_sofar = jpvt(1:ranknew)
      deallocate(mats)
      deallocate(core)
      deallocate(tau)
      deallocate(jpvt)



   else
      ranknew = size(treequant(index_ij)%skel_sofar,1)
      allocate(mats_skel(count,ranknew))
      mats_skel = matnew(:,treequant(index_ij)%skel_sofar)
      allocate(matA_recon(count,ncol))
      matA_recon=0
      call gemmf90(mats_skel, count, treequant(index_ij)%mats_interp, ranknew, matA_recon, count, 'N', 'N', count, ncol, ranknew, BPACK_cone,BPACK_czero)
      ! tol_next = fnorm(matnew,count,ncol)*tol_comp
      ! tol_next = fnorm(treequant(index_ij)%mats,size(treequant(index_ij)%mats,1),size(treequant(index_ij)%mats,2))*tol_comp/(depth+1)
      tol_next = fnorm(treequant(index_ij)%mats,size(treequant(index_ij)%mats,1),size(treequant(index_ij)%mats,2))*tol_comp
      residual = fnorm(matA_recon-matnew,count,ncol)
      deallocate(mats_skel)
      deallocate(matA_recon)

      write(*,*)'node ', nodeID, 'count', count, 'residual with ID', residual, tol_next

      if (residual < tol_next) then
         keeprefine = 0
      else
         keeprefine = 1



         count = size(treequant(index_ij)%mats,1)
         allocate (core(count, ncol))
         core = treequant(index_ij)%mats
         allocate(mats(count,ncol))
         mats = treequant(index_ij)%mats

         allocate (jpvt(max(ncol, count)))
         allocate (tau(max(ncol, count)))
         jpvt = 0
         call geqp3modf90(core, jpvt, tau, tol_comp, BPACK_SafeUnderflow, ranknew)
         if (ranknew > 0) then
            call un_or_mqrf90(core, tau, mats, 'L', 'C', size(mats,1), ncol, ranknew)
            call trsmf90(core, mats, 'L', 'U', 'N', 'N', ranknew, ncol)
         else
            ranknew = 1
            jpvt(1) = 1
            mats = 0
         endif
         call matrix_resize(treequant(index_ij)%mats_interp,ranknew,ncol)
         treequant(index_ij)%mats_interp =mats(1:ranknew, 1:ncol)
         call array_resize_int(treequant(index_ij)%skel_sofar,ranknew)
         treequant(index_ij)%skel_sofar = jpvt(1:ranknew)
         deallocate(mats)
         deallocate(core)
         deallocate(tau)
         deallocate(jpvt)









         ! call ComputeRange(ncol, count, mattmp3, ranknew, 1, tol_comp, norm_tol=norm_tol)
         ! write(*,*)'after ComputeRange: new basis# ', ranknew, 'node', nodeID

         ! call matrix_resize(treequant(index_ij)%Qs, ncol, nq + ranknew)
         ! treequant(index_ij)%Qs(:, 1 + nq : nq + ranknew) = mattmp3(:, 1:ranknew)
         ! rank = nq + ranknew

         ! ! Re-orthogonalize
         ! call ComputeRange(ncol, nq + ranknew, treequant(index_ij)%Qs, rank, 1, tol_comp, norm_tol=norm_tol)
         ! write(*,*)'after Re-orthogonalization: total basis#', rank
         ! call matrix_resize(treequant(index_ij)%Qs, ncol, rank)

         ! deallocate(mattmp3)
      end if

      ! deallocate(mattmp1)
      ! deallocate(mattmp2)
   end if

   deallocate(matnew)

   ! If we should refine, push children (right, then left) so left is processed next
   if (keeprefine == 1) then
      mid = (lo + hi) / 2

      ! Right child [mid+1, hi]
      if (mid + 1 <= hi) then
         if (treequant(index_ij)%tail == size(treequant(index_ij)%stack)) then
            ! In practice H above is enough; grow if you want absolute safety

            call stack_resize(treequant(index_ij)%stack,size(treequant(index_ij)%stack)*2)
            treequant(index_ij)%stack(1:treequant(index_ij)%tail-treequant(index_ij)%top+1)=treequant(index_ij)%stack(treequant(index_ij)%top:treequant(index_ij)%tail)
            treequant(index_ij)%tail = treequant(index_ij)%tail - treequant(index_ij)%top+1
            treequant(index_ij)%top = 1

         end if
         treequant(index_ij)%tail = treequant(index_ij)%tail + 1
         treequant(index_ij)%stack(treequant(index_ij)%tail)%lo      = mid + 1
         treequant(index_ij)%stack(treequant(index_ij)%tail)%hi      = hi
         treequant(index_ij)%stack(treequant(index_ij)%tail)%nodeID   = 2*nodeID + 1
         treequant(index_ij)%stack(treequant(index_ij)%tail)%m_pick  = m_adjust
         treequant(index_ij)%stack(treequant(index_ij)%tail)%tol_run = tol_next
      end if

      ! Left child [lo, mid]
      if (lo <= mid) then
         if (treequant(index_ij)%tail == size(treequant(index_ij)%stack)) then
            call stack_resize(treequant(index_ij)%stack,size(treequant(index_ij)%stack)*2)
            treequant(index_ij)%stack(1:treequant(index_ij)%tail-treequant(index_ij)%top+1)=treequant(index_ij)%stack(treequant(index_ij)%top:treequant(index_ij)%tail)
            treequant(index_ij)%tail = treequant(index_ij)%tail - treequant(index_ij)%top+1
            treequant(index_ij)%top = 1
         endif
            treequant(index_ij)%tail = treequant(index_ij)%tail + 1
            treequant(index_ij)%stack(treequant(index_ij)%tail)%lo      = lo
            treequant(index_ij)%stack(treequant(index_ij)%tail)%hi      = mid
            treequant(index_ij)%stack(treequant(index_ij)%tail)%nodeID   = 2*nodeID
            treequant(index_ij)%stack(treequant(index_ij)%tail)%m_pick  = m_adjust
            treequant(index_ij)%stack(treequant(index_ij)%tail)%tol_run = tol_next
         end if
      end if
   end do

   ! deallocate(stack)
   end subroutine refine_interval



   subroutine KD_tree_ordering(Nmin_leaf, msh)


      implicit none

      integer Cflag
      integer i, j, ii, jj, iii, jjj, kk
      integer level, edge, node, patch, group, group_m, group_n, col_group, row_group, fidx
      integer blocks
      integer center_edge

      integer index_temp
      integer rows(2), cols(2)
      real(kind=8) a, b, c, d, para, xmax, xmin, ymax, ymin, zmax, zmin, seperator, r, theta, phi, phi_tmp
      real(kind=8) radius, radiusmax, radius2, radiusmax2
      real(kind=8), allocatable:: xyzrange(:), xyzmin(:), xyzmax(:), auxpoint(:), groupcenter(:)
      real(kind=8), allocatable :: distance(:), array(:, :), dist_gram(:, :)
      integer level_c, sortdirec, nrow, phi_end, Ninfo_edge, ind_i, ind_j
      real(kind=8) t1, t2
      integer Maxgroup, nlevel_pre, passflag
      character(len=1024)  :: strings
      integer, allocatable :: order(:), edge_temp(:), map_temp(:)
      integer dimn, groupsize, idxstart, Nsmp
      integer Nmin_leaf
      integer Maxlevel, Maxgrp
      type(mesh)::msh

      integer, allocatable:: perms(:), rows_gram(:), cols_gram(:)
      integer Navr, Bidxs, Bidxe, ierr

      !>*************Initialize permutation vector ********
      allocate (msh%new2old(msh%Nunk))
      do ii = 1, msh%Nunk
         msh%new2old(ii) = ii
      end do

      !>************Compute Maxlevel of BPACK tree*******************
      nlevel_pre = 0
      if (allocated(msh%pretree)) then
         nlevel_pre = ceiling_safe(log(dble(size(msh%pretree, 1)))/log(2d0))
      endif
      level = 0; i = 1
      do while (int(msh%Nunk/i) > Nmin_leaf)
         level = level + 1
         i = 2**level
      enddo
      Maxlevel = level
      if (Maxlevel < nlevel_pre) Maxlevel = nlevel_pre


      !>***************************************************

      Maxgroup = 2**(Maxlevel + 1) - 1
      msh%Maxgroup = Maxgroup
      allocate (msh%basis_group(Maxgroup))

      dimn = 0
      if (allocated(msh%xyz)) Dimn = size(msh%xyz, 1)
      if (dimn > 0) then
         allocate (xyzrange(dimn))
         allocate (xyzmin(dimn))
         allocate (xyzmax(dimn))
         allocate (auxpoint(dimn))
         allocate (groupcenter(dimn))
      endif

      !>**** construct the top few levels whose ordering is provided by the user
      msh%basis_group(1)%head = 1; msh%basis_group(1)%tail = msh%Nunk
      do level = nlevel_pre, 0, -1
         idxstart = 1
         do group = 2**level, 2**(level + 1) - 1
            ! msh%basis_group(group)%level=level

            if (level == nlevel_pre) then
               if (nlevel_pre == 0) then
                  groupsize = msh%Nunk
               else
                  groupsize = msh%pretree(group - 2**nlevel_pre + 1)
               endif
               ! call assert(groupsize > 0, 'zero leafsize may not be handled')
               msh%basis_group(group)%head = idxstart
               msh%basis_group(group)%tail = idxstart + groupsize - 1
               idxstart = idxstart + groupsize
            else
               msh%basis_group(group)%head = msh%basis_group(2*group)%head
               msh%basis_group(group)%tail = msh%basis_group(2*group + 1)%tail
            endif

            !>***** the following is needed for the near_or_far function in H matrix, this needs to be improved
            if (allocated(msh%xyz)) then
               Dimn = size(msh%xyz, 1)
               groupcenter(1:dimn) = 0.0d0
               do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                  do ii = 1, dimn
                     groupcenter(ii) = groupcenter(ii) + msh%xyz(ii, msh%new2old(edge))
                  enddo
               enddo
               do ii = 1, dimn
                  groupcenter(ii) = groupcenter(ii)/(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1)
               enddo

               radiusmax = 0.
               do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                  radius = 0
                  do ii = 1, dimn
                     radius = radius + (msh%xyz(ii, msh%new2old(edge)) - groupcenter(ii))**2
                  enddo
                  radius = sqrt(radius)
                  if (radius > radiusmax) then
                     radiusmax = radius
                     center_edge = edge
                  endif
               enddo

               allocate (msh%basis_group(group)%center(dimn))
               msh%basis_group(group)%center = groupcenter
               msh%basis_group(group)%radius = radiusmax
            endif

         enddo
      enddo



         do level = nlevel_pre, Maxlevel
            do group = 2**level, 2**(level + 1) - 1
               ! msh%basis_group(group)%level=level

               if (allocated(msh%xyz)) then
               if (.not. allocated(msh%basis_group(group)%center)) then
                  groupcenter(1:dimn) = 0.0d0
                  ! !$omp parallel do default(shared) private(edge,ii) reduction(+:groupcenter)
                  do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                     do ii = 1, dimn
                        groupcenter(ii) = groupcenter(ii) + msh%xyz(ii, msh%new2old(edge))
                     enddo
                  enddo
                  ! !$omp end parallel do
                  do ii = 1, dimn
                     groupcenter(ii) = groupcenter(ii)/(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1)
                  enddo

                  radiusmax = 0.
                  do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                     radius = 0
                     do ii = 1, dimn
                        radius = radius + (msh%xyz(ii, msh%new2old(edge)) - groupcenter(ii))**2
                     enddo
                     radius = sqrt(radius)
                     if (radius > radiusmax) then
                        radiusmax = radius
                        center_edge = edge
                     endif
                  enddo

                  allocate (msh%basis_group(group)%center(dimn))
                  msh%basis_group(group)%center = groupcenter
                  msh%basis_group(group)%radius = radiusmax
               endif
               endif

               xyzmin = 1d300
               xyzmax = -1d300
               do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                  do ii = 1, Dimn
                     xyzmax(ii) = max(xyzmax(ii), msh%xyz(ii, msh%new2old(edge)))
                     xyzmin(ii) = min(xyzmin(ii), msh%xyz(ii, msh%new2old(edge)))
                  enddo
               enddo
               xyzrange(1:Dimn) = xyzmax(1:Dimn) - xyzmin(1:Dimn)

               nrow = msh%basis_group(group)%tail - msh%basis_group(group)%head + 1
               allocate (distance(nrow))
               sortdirec = maxloc(xyzrange(1:Dimn), 1)
               ! write(*,*)'gaw',sortdirec,xyzrange(1:Dimn)


#ifdef HAVE_OPENMP
                  !$omp parallel do default(shared) private(i)
#endif
                  do i = msh%basis_group(group)%head, msh%basis_group(group)%tail
                     distance(i - msh%basis_group(group)%head + 1) = msh%xyz(sortdirec, msh%new2old(i))
                  enddo
#ifdef HAVE_OPENMP
                  !$omp end parallel do
#endif

               allocate (order(nrow))
               allocate (map_temp(nrow))

               call quick_sort(distance, order, nrow)
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(ii)
#endif
               do ii = 1, nrow
                  map_temp(ii) = msh%new2old(order(ii) + msh%basis_group(group)%head - 1)
               enddo
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
#ifdef HAVE_OPENMP
               !$omp parallel do default(shared) private(ii)
#endif
               do ii = 1, nrow
                  msh%new2old(ii + msh%basis_group(group)%head - 1) = map_temp(ii)
               enddo
#ifdef HAVE_OPENMP
               !$omp end parallel do
#endif
               deallocate (map_temp)
               deallocate (order)

               deallocate (distance)

               if (level < Maxlevel) then

                  ! call assert(msh%basis_group(group)%tail /= msh%basis_group(group)%head, 'detected zero-sized group, try larger leafsizes or smaller MPI counts')
                  msh%basis_group(2*group)%head = msh%basis_group(group)%head
                  msh%basis_group(2*group)%tail = int((msh%basis_group(group)%head + msh%basis_group(group)%tail)/2)
                  msh%basis_group(2*group + 1)%head = msh%basis_group(2*group)%tail + 1
                  msh%basis_group(2*group + 1)%tail = msh%basis_group(group)%tail
               endif
            enddo
         enddo

      !>**** generate tree structures on other processes
      do level = nlevel_pre, Maxlevel
         do group = 2**level, 2**(level + 1) - 1
            ! msh%basis_group(group)%level=level

            if (allocated(msh%xyz)) then
            if (.not. allocated(msh%basis_group(group)%center)) then
               groupcenter(1:dimn) = 0.0d0
               ! !$omp parallel do default(shared) private(edge,ii) reduction(+:groupcenter)
               do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                  do ii = 1, dimn
                     groupcenter(ii) = groupcenter(ii) + msh%xyz(ii, msh%new2old(edge))
                  enddo
               enddo
               ! !$omp end parallel do
               do ii = 1, dimn
                  groupcenter(ii) = groupcenter(ii)/(msh%basis_group(group)%tail - msh%basis_group(group)%head + 1)
               enddo

               radiusmax = 0.
               do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
                  radius = 0
                  do ii = 1, dimn
                     radius = radius + (msh%xyz(ii, msh%new2old(edge)) - groupcenter(ii))**2
                  enddo
                  radius = sqrt(radius)
                  if (radius > radiusmax) then
                     radiusmax = radius
                  endif
               enddo

               allocate (msh%basis_group(group)%center(dimn))
               msh%basis_group(group)%center = groupcenter
               msh%basis_group(group)%radius = radiusmax
            endif
            endif

            xyzmin = 1d300
            xyzmax = -1d300
            do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
               do ii = 1, Dimn
                  xyzmax(ii) = max(xyzmax(ii), msh%xyz(ii, msh%new2old(edge)))
                  xyzmin(ii) = min(xyzmin(ii), msh%xyz(ii, msh%new2old(edge)))
               enddo
            enddo
            xyzrange(1:Dimn) = xyzmax(1:Dimn) - xyzmin(1:Dimn)
            sortdirec = maxloc(xyzrange(1:Dimn), 1)
            seperator = msh%xyz(sortdirec, msh%new2old(int((msh%basis_group(group)%head + msh%basis_group(group)%tail)/2)))
            msh%basis_group(group)%boundary(1) = sortdirec
            msh%basis_group(group)%boundary(2) = seperator

            if (level < Maxlevel) then
               ! call assert(msh%basis_group(group)%tail /= msh%basis_group(group)%head, 'detected zero-sized group, try larger leafsizes or smaller MPI counts')
               msh%basis_group(2*group)%head = msh%basis_group(group)%head
               msh%basis_group(2*group)%tail = int((msh%basis_group(group)%head + msh%basis_group(group)%tail)/2)
               msh%basis_group(2*group + 1)%head = msh%basis_group(2*group)%tail + 1
               msh%basis_group(2*group + 1)%tail = msh%basis_group(group)%tail
            endif
         enddo
      enddo

      if (dimn > 0) then
         deallocate (xyzrange)
         deallocate (xyzmin)
         deallocate (xyzmax)
         deallocate (auxpoint)
         deallocate (groupcenter)
      endif

      allocate (msh%old2new(msh%Nunk))
      do ii = 1, msh%Nunk
         msh%old2new(msh%new2old(ii)) = ii
      end do

      ! do ii=1,msh%Nunk
      ! write(110,*)msh%old2new(ii)
      ! enddo

      !>**********Dump the ordering into a file********************************

#if        0
      write (strings, *) Dimn
      do level = 0, Maxlevel
         do group = 2**level, 2**(level + 1) - 1
            do edge = msh%basis_group(group)%head, msh%basis_group(group)%tail
               write (113, '(I5,I8,'//TRIM(strings)//'Es16.8)') level, group, msh%xyz(1:Dimn, msh%new2old(edge))
            enddo
         enddo
      enddo
#endif



      return

   end subroutine KD_tree_ordering


end program test_proxy_selection








