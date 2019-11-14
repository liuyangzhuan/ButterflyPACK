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
module Bplus_factor
use Bplus_compress
use Bplus_randomizedop

contains


subroutine Full_LU(blocks,option,stats)

    implicit none

	type(Hoption)::option
	type(Hstat)::stats

    integer size_m, size_n
    integer i,j,k,ii,jj,kk
    real*8  T0,T1
    type(matrixblock) :: blocks
    real(kind=8) flop

    T0=OMP_get_wtime()
    size_m=size(blocks%fullmat,1)
    if(option%ILU==0)then
		! do ii=1,size_m
		! do jj=1,size_m
			! write(777,*)dble(blocks%fullmat(ii,jj)),aimag(blocks%fullmat(ii,jj))
		! enddo
		! enddo
		call getrff90(blocks%fullmat,blocks%ipiv,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
		! do ii=1,size_m
		! do jj=1,size_m
			! write(778,*)dble(blocks%fullmat(ii,jj)),aimag(blocks%fullmat(ii,jj))
		! enddo
		! enddo
	else
		do ii=1,size_m
			blocks%ipiv(ii)=ii
		enddo
	endif
	T1=OMP_get_wtime()
    stats%Time_Direct_LU=stats%Time_Direct_LU+T1-T0

    return

end subroutine Full_LU


subroutine Full_add_multiply(block3,chara,block1,block2,h_mat,option,stats,ptree,msh)

    implicit none

	type(Hoption)::option
	type(Hstat)::stats
	type(Hmat)::h_mat
	type(proctree)::ptree
	type(mesh)::msh

    integer level_butterfly, flag
    integer i, j, k, ii, level, mm, nn, kk, rank, level_blocks, mn, group_k
    integer style(3), data_type(3), id1, id2, id3
    character chara
    DT,allocatable::Vin(:,:),Vin1(:,:),fullmat(:,:),fullmatrix(:,:)
    real*8 T0, T1
    type(matrixblock) :: block1, block2, block3

	stats%Flop_Tmp=0

	T0=OMP_get_wtime()
    style(3)=block3%style
    level_blocks=block3%level

	group_k = block1%col_group
	kk=msh%basis_group(group_k)%tail-msh%basis_group(group_k)%head+1

	call assert(style(3)==1,'block3 supposed to be style 1')

	mm=size(block3%fullmat,1)
    nn=size(block3%fullmat,2)

	allocate(Vin(nn,nn))
	Vin = 0d0
	do ii=1,nn
		Vin(ii,ii)=1d0
	enddo
	allocate(Vin1(kk,nn))
	Vin1 = 0d0
	allocate(fullmatrix(mm,nn))
	fullmatrix=0d0

	call Hmat_block_MVP_dat(block2,'N',msh%basis_group(block2%row_group)%head,msh%basis_group(block2%col_group)%head,nn,Vin,Vin1,cone,ptree,stats)
	call Hmat_block_MVP_dat(block1,'N',msh%basis_group(block1%row_group)%head,msh%basis_group(block1%col_group)%head,nn,Vin1,fullmatrix,cone,ptree,stats)

	if (chara=='-')fullmatrix = -fullmatrix
	block3%fullmat = block3%fullmat + fullmatrix
	deallocate(fullmatrix)
	deallocate(Vin)
	deallocate(Vin1)


    T1=OMP_get_wtime()
    stats%Time_Add_Multiply=stats%Time_Add_Multiply+T1-T0
    stats%Flop_Factor = stats%Flop_Factor+stats%Flop_Tmp

    return

end subroutine Full_add_multiply

subroutine Full_add(block3,chara,block1,ptree,stats)

    implicit none

    integer level_butterfly, flag,group_n,group_m
    integer i, j, k, level, mm, nn, rank, level_blocks, mn, ii, jj
    integer style(3), data_type(3), id1, id2, id3
    character chara
    DT,allocatable:: Vin(:,:)
    !logical isNaN
    real*8 T0, T1
    type(matrixblock) :: block1, block3
	type(Hstat):: stats
	type(proctree):: ptree
    DT,allocatable::fullmatrix(:,:)

	stats%Flop_Tmp=0

    style(3)=block3%style
    style(1)=block1%style
    level_blocks=block3%level


    T0=OMP_get_wtime()
    call assert(style(1)/=1,'block1 not supposed to be full')
    call assert(style(3)==1,'block3 supposed to be full')

	group_m = block3%row_group
	group_n = block3%col_group

	mm=block3%M
	nn=block3%N
	allocate(Vin(nn,nn))

	Vin = 0d0
	do ii=1,nn
		Vin(ii,ii)=1d0
	enddo

	allocate(fullmatrix(mm,nn))
	fullmatrix=0d0

	call BF_block_MVP_dat(block1,'N',mm,nn,nn,Vin,fullmatrix,cone,czero,ptree,stats)

	if (chara=='-')fullmatrix = -fullmatrix

	block3%fullmat = block3%fullmat + fullmatrix
	deallocate(fullmatrix)
	deallocate(Vin)


    T1=OMP_get_wtime()
    stats%Time_Add_Multiply=stats%Time_Add_Multiply+T1-T0
    stats%Flop_Factor = stats%Flop_Factor+stats%Flop_Tmp
    return

end subroutine Full_add


subroutine LR_minusBC(ho_bf1,level_c,rowblock,ptree,stats)

    use BPACK_DEFS

	use MISC_Utilities
    implicit none

	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm_diag
    character chara
    real(kind=8) a,b,c,d
    DT ctemp1, ctemp2
	type(matrixblock),pointer::block_o

    ! type(vectorsblock), pointer :: random1, random2

    real(kind=8),allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	! DT,allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:),myA(:,:),BUold(:,:),BVold(:,:),CUold(:,:),CVold(:,:),BU(:,:),BV(:,:),CU(:,:),CV(:,:),BVCU(:,:),BUBVCU(:,:)

	integer Nsub,Ng,unique_nth,level_left_start,ll
	integer*8 idx_start
    integer level_blocks
    integer groupm_start, groupn_start,dimension_rank,rank1,rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n

	real(kind=8)::n2,n1
	type(hobf)::ho_bf1
	type(matrixblock),pointer::block_off1,block_off2
	type(proctree)::ptree
	integer pgno,pgno1,pgno2
	integer descBUold(9),descBVold(9),descCUold(9),descCVold(9), descBU(9),descBV(9),descCU(9),descCV(9),descBVCU(9),descBUBVCU(9)
	integer ctxt1,ctxt2,ctxt,ctxtall,info,myrow,mycol,myArows,myAcols
	type(Hstat)::stats

	ctemp1=1.0d0 ; ctemp2=0.0d0
	block_off1 => ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2-1)%LL(1)%matrices_block(1)
	block_off2 => ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2)%LL(1)%matrices_block(1)

	block_o =>  ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
	call BF_delete(block_o,1)

	block_o%level_butterfly = 0
	block_o%ButterflyU%nblk_loc=1
	block_o%ButterflyU%inc=1
	block_o%ButterflyU%idx=1
	block_o%ButterflyV%nblk_loc=1
	block_o%ButterflyV%inc=1
	block_o%ButterflyV%idx=1


	allocate(block_o%ButterflyU%blocks(1))
	allocate(block_o%ButterflyV%blocks(1))

	pgno = block_o%pgno
	pgno1 = block_off1%pgno
	pgno2 = block_off2%pgno
	! if(ptree%MyID==2)write(*,*)pgno,pgno1,pgno2,'nana'
	call assert(pgno==pgno1,'block_o and block_off1 should be on the same process group')
	call assert(pgno==pgno2,'block_o and block_off2 should be on the same process group')

	mm = block_off1%M
	nn = block_off1%N

	rank = block_off2%rankmax
	block_o%rankmax=rank
	block_o%M_loc = block_off1%M_loc
	block_o%N_loc = block_o%M_loc
	allocate(block_o%M_p(size(block_off1%M_p,1),2))
	block_o%M_p = block_off1%M_p
	allocate(block_o%N_p(size(block_off1%M_p,1),2))
	block_o%N_p = block_off1%M_p

	allocate(block_o%ButterflyU%blocks(1)%matrix(block_o%M_loc,rank))
	block_o%ButterflyU%blocks(1)%matrix =0
	allocate(block_o%ButterflyV%blocks(1)%matrix(block_o%M_loc,rank))
	block_o%ButterflyV%blocks(1)%matrix =block_off2%ButterflyV%blocks(1)%matrix

	stats%Flop_Tmp=0
	call BF_block_MVP_dat(block_off1,'N',block_off1%M_loc,block_off1%N_loc,rank,block_off2%ButterflyU%blocks(1)%matrix,block_o%ButterflyU%blocks(1)%matrix,ctemp1,ctemp2,ptree,stats)
	block_o%ButterflyU%blocks(1)%matrix = -block_o%ButterflyU%blocks(1)%matrix
	! if(ptree%MyID==2)write(*,*)pgno,pgno1,pgno2,'neeeeana'
	stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp
	stats%Flop_Tmp=0

    return

end subroutine LR_minusBC


subroutine LR_SMW(block_o,Memory,ptree,stats,pgno)

    use BPACK_DEFS


    implicit none

    integer level_c,rowblock,kover,rank,kk1,kk2,nprow,npcol
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn
    integer dimension_rank, dimension_m, dimension_n, blocks, groupm, groupn,index_j,index_i
    real(kind=8) a,b,c,d,Memory,flop
    DT ctemp,TEMP(1)
	type(matrixblock)::block_o
	DT, allocatable::matrixtemp(:,:),matrixtemp1(:,:),matrixtemp2(:,:),matrixtemp3(:,:),UU(:,:),VV(:,:),matrix_small(:,:),vin(:,:),vout1(:,:),vout2(:,:),vout3(:,:),matU(:,:)
	real(kind=8), allocatable:: Singular(:)
    integer, allocatable :: ipiv(:),iwork(:)
	type(proctree)::ptree
	integer pgno,ctxt,ctxt_head,myrow,mycol,myArows,myAcols,iproc,myi,jproc,myj,info
	integer descUV(9),descsmall(9),desctemp(9),TEMPI(1)

	integer lwork,liwork,lcmrc,ierr
	DT,allocatable:: work(:)
	type(Hstat)::stats

	ctxt = ptree%pgrp(pgno)%ctxt
	ctxt_head = ptree%pgrp(pgno)%ctxt_head

	rank = size(block_o%ButterflyU%blocks(1)%matrix,2)
	allocate(matrixtemp(rank,rank))
	matrixtemp=0
	allocate(matrixtemp1(rank,rank))
	matrixtemp1=0
	allocate(matU(block_o%M_loc,rank))
	matU = block_o%ButterflyU%blocks(1)%matrix

	! write(*,*)fnorm(block_o%ButterflyV%blocks(1)%matrix,size(block_o%ButterflyV%blocks(1)%matrix,1),size(block_o%ButterflyV%blocks(1)%matrix,2)),fnorm(block_o%ButterflyU%blocks(1)%matrix,size(block_o%ButterflyU%blocks(1)%matrix,1),size(block_o%ButterflyU%blocks(1)%matrix,2)),ptree%MyID,'re',shape(block_o%ButterflyV%blocks(1)%matrix),shape(block_o%ButterflyU%blocks(1)%matrix),shape(matrixtemp),isnanMat(block_o%ButterflyV%blocks(1)%matrix,size(block_o%ButterflyV%blocks(1)%matrix,1),size(block_o%ButterflyV%blocks(1)%matrix,2)),isnanMat(block_o%ButterflyU%blocks(1)%matrix,size(block_o%ButterflyU%blocks(1)%matrix,1),size(block_o%ButterflyU%blocks(1)%matrix,2))

	call gemmf90(block_o%ButterflyV%blocks(1)%matrix,block_o%M_loc,block_o%ButterflyU%blocks(1)%matrix,block_o%M_loc,matrixtemp,rank,'T','N',rank,rank,block_o%M_loc,cone,czero,flop=flop)
	stats%Flop_Factor = stats%Flop_Factor + flop

	! write(*,*)'goog1'
	call assert(MPI_COMM_NULL/=ptree%pgrp(pgno)%Comm,'communicator should not be null 1')
	call MPI_ALLREDUCE(matrixtemp,matrixtemp1,rank*rank,MPI_DT,MPI_SUM,ptree%pgrp(pgno)%Comm,ierr)
	! write(*,*)'goog2'
	do ii=1,rank
		matrixtemp1(ii,ii) = matrixtemp1(ii,ii)+1
	enddo

	! write(*,*)abs(matrixtemp1),rank,'gggddd'

	if(rank<=nbslpk)then

#if 0
		allocate(ipiv(rank))
		ipiv=0
		call getrff90(matrixtemp1,ipiv,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
		call getrif90(matrixtemp1,ipiv,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
		deallocate(ipiv)
#else
		matrixtemp = matrixtemp1
		call GeneralInverse(rank,rank,matrixtemp,matrixtemp1,SafeEps,Flops=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
#endif

	else

		!!!!!! the SVD-based pseudo inverse needs to be implemented later

		call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
		if(myrow/=-1 .and. mycol/=-1)then
			if(ptree%MyID==ptree%pgrp(pgno)%head)then
				call blacs_gridinfo(ctxt_head, nprow, npcol, myrow, mycol)
				myArows = numroc_wp(rank, nbslpk, myrow, 0, nprow)
				call descinit( desctemp, rank, rank, nbslpk, nbslpk, 0, 0, ctxt_head, max(myArows,1), info )
				call assert(info==0,'descinit fail for desctemp')
			else
				desctemp(2)=-1
			endif

			call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
			myArows = numroc_wp(rank, nbslpk, myrow, 0, nprow)
			myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)

			allocate(matrix_small(myArows,myAcols))
			matrix_small=0

			call descinit( descsmall, rank, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
			! if(info/=0)then
				! write(*,*)'nneref',rank,nbslpk,myArows,myAcols,max(myArows,1),ptree%pgrp(pgno)%nproc,ptree%MyID,ptree%pgrp(pgno)%head,pgno,ptree%pgrp(pgno)%nprow,ptree%pgrp(pgno)%npcol,info
			! endif
			call assert(info==0,'descinit fail for descsmall')

			call pgemr2df90(rank, rank, matrixtemp1, 1, 1, desctemp, matrix_small, 1, 1, descsmall, ctxt)

			allocate(ipiv(myArows+nbslpk))
			ipiv=0
			call pgetrff90(rank,rank,matrix_small,1,1,descsmall,ipiv,info,flop=flop)
			stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)

			call pgetrif90(rank,matrix_small,1,1,descsmall,ipiv,flop=flop)
			stats%Flop_Factor = stats%Flop_Factor + flop/dble(nprow*npcol)

			deallocate(ipiv)


			call pgemr2df90(rank, rank, matrix_small, 1, 1, descsmall,matrixtemp1, 1, 1, desctemp, ctxt)
			deallocate(matrix_small)
		endif

		call MPI_Bcast(matrixtemp1,rank*rank,MPI_DT,0,ptree%pgrp(pgno)%Comm,ierr)

	endif


	call gemmf90(matU,block_o%M_loc,matrixtemp1,rank,block_o%ButterflyU%blocks(1)%matrix,block_o%M_loc,'N','N',block_o%M_loc,rank,rank,cone,czero,flop=flop)
	block_o%ButterflyU%blocks(1)%matrix = -block_o%ButterflyU%blocks(1)%matrix
	stats%Flop_Factor = stats%Flop_Factor + flop

	deallocate(matrixtemp,matrixtemp1,matU)



	Memory = 0
	Memory = Memory + SIZEOF(block_o%ButterflyV%blocks(1)%matrix)/1024.0d3
	Memory = Memory + SIZEOF(block_o%ButterflyU%blocks(1)%matrix)/1024.0d3

    return

end subroutine LR_SMW



subroutine LR_Sblock(ho_bf1,level_c,rowblock,ptree,stats)

    use BPACK_DEFS

	use MISC_Utilities
    implicit none

	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,pp,qq
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real(kind=8) a,b,c,d
	type(matrixblock),pointer::block_o,blocks

    type(vectorsblock), pointer :: random1, random2

    real(kind=8),allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_end_diag,idx_start_loc,idx_end_loc
	DT,allocatable::vec_old(:,:),vec_new(:,:)

	integer Nsub,Ng,unique_nth,level_left_start
	integer*8 idx_start
    integer level_blocks,head,tail
    integer groupm_start, groupn_start,dimension_rank
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n

	integer nth_s,nth_e,num_vect_sub,nth
	real(kind=8)::n2,n1
	type(hobf)::ho_bf1
	type(proctree)::ptree
	type(Hstat)::stats



	block_o =>  ho_bf1%levels(level_c)%BP_inverse_update(rowblock)%LL(1)%matrices_block(1)

	! write(*,*)block_o%row_group,block_o%col_group,isnanMat(block_o%ButterflyU%blocks(1)%matrix,size(block_o%ButterflyU%blocks(1)%matrix,1),size(block_o%ButterflyU%blocks(1)%matrix,2)),'dfdU1',ptree%MyID
	! write(*,*)block_o%row_group,block_o%col_group,isnanMat(block_o%ButterflyV%blocks(1)%matrix,size(block_o%ButterflyV%blocks(1)%matrix,1),size(block_o%ButterflyV%blocks(1)%matrix,2)),'dfdV1',ptree%MyID

    level_butterfly=block_o%level_butterfly
    call assert(level_butterfly==0,'Butterfly_Sblock_LowRank only works with LowRank blocks')



	num_blocks=2**level_butterfly


	num_vect_sub = size(block_o%ButterflyU%blocks(1)%matrix,2)
    ! groupm=block_o%row_group  ! Note: row_group and col_group interchanged here

	! get the right multiplied vectors
	pp = ptree%myid-ptree%pgrp(block_o%pgno)%head+1
	idx_start_glo = block_o%headm + block_o%M_p(pp,1) -1



	! mm=block_o%M
	mm=block_o%M_loc
	allocate(vec_old(mm,num_vect_sub))
	allocate(vec_new(mm,num_vect_sub))
	vec_old = block_o%ButterflyU%blocks(1)%matrix
	stats%Flop_Tmp=0
	do level = ho_bf1%Maxlevel+1,level_c+1,-1
		N_diag = 2**(level-level_c-1)
		idx_start_diag = max((rowblock-1)*N_diag+1,ho_bf1%levels(level)%Bidxs)
		idx_end_diag = min(rowblock*N_diag,ho_bf1%levels(level)%Bidxe)
		vec_new = 0

		n1 = OMP_get_wtime()
		do ii = idx_start_diag,idx_end_diag

			if(associated(ho_bf1%levels(level)%BP_inverse(ii)%LL))then
			blocks=>ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)
			if(IOwnPgrp(ptree,blocks%pgno))then

			qq = ptree%myid-ptree%pgrp(blocks%pgno)%head+1
			head = blocks%headm + blocks%M_p(qq,1) -1
			tail = head + blocks%M_loc - 1
			idx_start_loc = head-idx_start_glo+1
			idx_end_loc = tail-idx_start_glo+1
			if(level==ho_bf1%Maxlevel+1)then
				call Full_block_MVP_dat(blocks,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,&
				&vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),cone,czero)
			else
				call BF_block_MVP_inverse_dat(ho_bf1,level,ii,'N',idx_end_loc-idx_start_loc+1,num_vect_sub,vec_old(idx_start_loc:idx_end_loc,1:num_vect_sub),vec_new(idx_start_loc:idx_end_loc,1:num_vect_sub),ptree,stats)
			endif

			endif
			endif
		end do
		n2 = OMP_get_wtime()
		! time_tmp = time_tmp + n2 - n1


		vec_old = vec_new
	end do
	! ! write(*,*)vec_new(1,1),RandomVectors_InOutput(2)%vector(1,1)
	block_o%ButterflyU%blocks(1)%matrix = vec_new
	! write(*,*)block_o%row_group,block_o%col_group,isnanMat(block_o%ButterflyU%blocks(1)%matrix,size(block_o%ButterflyU%blocks(1)%matrix,1),size(block_o%ButterflyU%blocks(1)%matrix,2)),'dfdU',ptree%MyID
	! write(*,*)block_o%row_group,block_o%col_group,isnanMat(block_o%ButterflyV%blocks(1)%matrix,size(block_o%ButterflyV%blocks(1)%matrix,1),size(block_o%ButterflyV%blocks(1)%matrix,2)),'dfdV',ptree%MyID
	deallocate(vec_old)
	deallocate(vec_new)

	stats%Flop_Factor = stats%Flop_Factor + stats%Flop_Tmp

    return

end subroutine LR_Sblock



subroutine LR_A_minusBDinvC(partitioned_block,ptree,option,stats)
   use BPACK_DEFS


   implicit none
   integer level, ii, num_vect_sub,mv,nv
   DT,allocatable :: V_tmp(:,:),V_tmp2(:,:),Vin_tmp(:,:),Vinter_tmp(:,:),Vout_tmp(:,:),U_tmp(:,:)
   DT :: ctemp1,ctemp2,a,b
   type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
   integer groupn,groupm,mm,nn,rank,rank1,rank2
   type(matrixblock)::partitioned_block
	type(proctree)::ptree
	type(Hstat)::stats
	type(Hoption)::option
	integer ctxt, nprow, npcol, myrow, mycol,pgno,myArows,myAcols,iproc,myi,jproc,myj,info
	DT,allocatable:: matU(:,:),matV(:,:),matU2D(:,:),matV2D(:,:),matU2Dnew(:,:),matV2Dnew(:,:)
	integer::descsMatU2D(9),descsMatV2D(9),descsMatSml(9),descsMatSmlRR1(9),descsMatSmlRR2(9),descsMatU2Dnew(9),descsMatV2Dnew(9),descsUUSml(9),descsVVSml(9),descsUU_u(9),descsVV_u(9),descsUU_v(9),descsVV_v(9)
	DT, allocatable :: QQ1(:,:), RR1(:,:),QQ2(:,:), RR2(:,:), UUsml(:,:), VVsml(:,:),UUu(:,:),UUv(:,:), VVu(:,:),VVv(:,:),tau_Q(:),mattemp(:,:),matU1(:,:),matV1(:,:)
	real(kind=8), allocatable :: Singularsml(:),Singularuv(:)
	real(kind=8)::flop
	integer ierr,mn1,mn2,jj,ranknew


		blocks_A => partitioned_block%sons(1,1)
		blocks_B => partitioned_block%sons(1,2)
		blocks_C => partitioned_block%sons(2,1)
		blocks_D => partitioned_block%sons(2,2)


		groupn=blocks_B%col_group    ! Note: row_group and col_group interchanged here
		nn=blocks_B%N
		groupm=blocks_B%row_group    ! Note: row_group and col_group interchanged here
		mm=blocks_B%M
		num_vect_sub = blocks_B%rankmax
		allocate(Vin_tmp(blocks_B%N_loc,num_vect_sub))
		Vin_tmp = blocks_B%ButterflyV%blocks(1)%matrix
		allocate(Vout_tmp(blocks_B%M_loc,num_vect_sub))
		Vout_tmp=0
		allocate(Vinter_tmp(blocks_B%N_loc,num_vect_sub))
		Vinter_tmp=0


		ctemp1=1.0d0 ; ctemp2=0.0d0
		call BF_block_MVP_dat(blocks_D,'T',nn,nn,num_vect_sub,Vin_tmp,Vinter_tmp,ctemp1,ctemp2,ptree,stats)
		Vinter_tmp = Vinter_tmp + Vin_tmp
		ctemp1=1.0d0 ; ctemp2=0.0d0
		call BF_block_MVP_dat(blocks_C,'T',nn,mm,num_vect_sub,Vinter_tmp,Vout_tmp,ctemp1,ctemp2,ptree,stats)
		rank = blocks_B%rankmax + blocks_A%rankmax
		allocate(matU(blocks_A%M_loc,rank))
		matU(:,1:blocks_A%rankmax)=blocks_A%ButterflyU%blocks(1)%matrix
		matU(:,1+blocks_A%rankmax:rank)=blocks_B%ButterflyU%blocks(1)%matrix
		allocate(matV(blocks_A%N_loc,rank))
		matV(:,1:blocks_A%rankmax)=blocks_A%ButterflyV%blocks(1)%matrix
		matV(:,1+blocks_A%rankmax:rank)=-Vout_tmp

		! deallocate(blocks_A%ButterflyU%blocks(1)%matrix)
		! deallocate(blocks_A%ButterflyV%blocks(1)%matrix)
		! blocks_A%rankmax=rank
		! blocks_A%rankmin=rank
		! allocate(blocks_A%ButterflyU%blocks(1)%matrix(blocks_A%M_loc,rank))
		! allocate(blocks_A%ButterflyV%blocks(1)%matrix(blocks_A%N_loc,rank))
		! blocks_A%ButterflyU%blocks(1)%matrix=matU
		! blocks_A%ButterflyV%blocks(1)%matrix=matV


! SVD recompression

	!!!!**** generate 2D grid blacs quantities for matU
	pgno = blocks_A%pgno
	ctxt = ptree%pgrp(pgno)%ctxt
	call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
	if(myrow/=-1 .and. mycol/=-1)then
		myArows = numroc_wp(blocks_A%M, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
		call descinit( descsMatU2D, blocks_A%M, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		allocate(matU2D(myArows,myAcols))
		matU2D=0

	else
		descsMatU2D(2)=-1
		allocate(matU2D(1,1))
		matU2D=0
	endif
	!!!!**** redistribution of input matrix
	call Redistribute1Dto2D(matU,blocks_A%M_p,0,pgno,matU2D,blocks_A%M,0,pgno,rank,ptree)

	!!!!**** generate 2D grid blacs quantities for matV transpose
	! ctxt = ptree%pgrp(pgno)%ctxt
	call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
	if(myrow/=-1 .and. mycol/=-1)then
		myArows = numroc_wp(blocks_A%N, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
		call descinit( descsMatV2D, blocks_A%N, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		allocate(matV2D(myArows,myAcols))
		matV2D=0
	else
		descsMatV2D(2)=-1
		allocate(matV2D(1,1))
		matV2D=0
	endif
	!!!!**** redistribution of input matrix


	call Redistribute1Dto2D(matV,blocks_A%N_p,0,pgno,matV2D,blocks_A%N,0,pgno,rank,ptree)


	if(myrow/=-1 .and. mycol/=-1)then


 ! GCC 9 was segfaulting when using PXGEQRF when it calls PXGEQR2 ->PXLARFG->PXSCAL
#if __GNUC__ < 9

		mn1=min(blocks_A%M,rank)
		myArows = numroc_wp(blocks_A%M, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
		allocate (tau_Q(myAcols))
		call pgeqrff90(blocks_A%M,rank,matU2D,1,1,descsMatU2D,tau_Q,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop

		myArows = numroc_wp(mn1, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
		call descinit( descsMatSmlRR1, mn1, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		allocate (RR1(myArows,myAcols))
		RR1=0d0
		do myj=1,myAcols
			call l2g(myj,mycol,mn1,npcol,nbslpk,jj)
			do myi=1,myArows
				call l2g(myi,myrow,rank,nprow,nbslpk,ii)
				if(ii<=jj)RR1(myi,myj)=matU2D(myi,myj)
			enddo
		enddo
		call pun_or_gqrf90(ctxt,matU2D,tau_Q,blocks_A%M,mn1,mn1,descsMatU2D,1,1,flop=flop)
		deallocate(tau_Q)
		stats%Flop_Factor = stats%Flop_Factor + flop

		mn2=min(blocks_A%N,rank)
		myArows = numroc_wp(blocks_A%N, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
		allocate (tau_Q(myAcols))
		call pgeqrff90(blocks_A%N,rank,matV2D,1,1,descsMatV2D,tau_Q,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop

		myArows = numroc_wp(mn2, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
		call descinit( descsMatSmlRR2, mn2, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		allocate (RR2(myArows,myAcols))
		RR2=0d0
		do myj=1,myAcols
			call l2g(myj,mycol,mn2,npcol,nbslpk,jj)
			do myi=1,myArows
				call l2g(myi,myrow,rank,nprow,nbslpk,ii)
				if(ii<=jj)RR2(myi,myj)=matV2D(myi,myj)
			enddo
		enddo
		call pun_or_gqrf90(ctxt,matV2D,tau_Q,blocks_A%N,mn2,mn2,descsMatV2D,1,1,flop=flop)
		deallocate(tau_Q)
		stats%Flop_Factor = stats%Flop_Factor + flop

		myArows = numroc_wp(mn1, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(mn2, nbslpk, mycol, 0, npcol)
		allocate(mattemp(myArows,myAcols))
		mattemp=0
		call descinit( descsMatSml, mn1, mn2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		call pgemmf90('N','T',mn1,mn2,rank,cone,RR1,1,1,descsMatSmlRR1,RR2,1,1,descsMatSmlRR2,czero,mattemp,1,1,descsMatSml,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
		myArows = numroc_wp(mn1, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(min(mn1,mn2), nbslpk, mycol, 0, npcol)
		call descinit( descsUUSml, mn1, min(mn1,mn2), nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		allocate(UUsml(myArows,myAcols))
		myArows = numroc_wp(min(mn1,mn2), nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(mn2, nbslpk, mycol, 0, npcol)
		call descinit( descsVVSml, min(mn1,mn2), mn2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		allocate(VVsml(myArows,myAcols))

		allocate(Singularsml(min(mn1,mn2)))
		call PSVD_Truncate(mn1,mn2,mattemp,descsMatSml,UUsml,VVsml,descsUUSml,descsVVSml,Singularsml,option%tol_rand,ranknew,ctxt,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
		myArows = numroc_wp(blocks_A%M, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
		allocate(matU2Dnew(myArows,myAcols))
		matU2Dnew=0
		call descinit( descsMatU2Dnew, blocks_A%M, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		call pgemmf90('N','N',blocks_A%M,ranknew,mn1,cone,matU2D,1,1,descsMatU2D,UUsml,1,1,descsUUSml,czero,matU2Dnew,1,1,descsMatU2Dnew,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
		do myj=1,myAcols
			call l2g(myj,mycol,ranknew,npcol,nbslpk,jj)
			matU2Dnew(:,myj)=matU2Dnew(:,myj)*Singularsml(jj)
		enddo


		myArows = numroc_wp(blocks_A%N, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
		allocate(matV2Dnew(myArows,myAcols))
		matV2Dnew=0
		call descinit( descsMatV2Dnew, blocks_A%N, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		call pgemmf90('N','T',blocks_A%N,ranknew,mn2,cone,matV2D,1,1,descsMatV2D,VVsml,1,1,descsVVSml,czero,matV2Dnew,1,1,descsMatV2Dnew,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
		rank = ranknew

		deallocate(mattemp,RR1,UUsml,VVsml,Singularsml)
		deallocate(RR2)

#else
		mn1=min(blocks_A%M,rank)
		myArows = numroc_wp(blocks_A%M, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(mn1, nbslpk, mycol, 0, npcol)
		call descinit( descsUU_u, blocks_A%M, mn1, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		allocate(UUu(myArows,myAcols))
		myArows = numroc_wp(mn1, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
		call descinit( descsVV_u, mn1, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		allocate(VVu(myArows,myAcols))
		allocate(Singularuv(mn1))
		call PSVD_Truncate(blocks_A%M,rank,matU2D,descsMatU2D,UUu,VVu,descsUU_u,descsVV_u,Singularuv,option%tol_rand,rank1,ctxt,flop=flop)
		do ii=1,rank1
			call g2l(ii,rank1,nprow,nbslpk,iproc,myi)
			if(iproc==myrow)then
				VVu(myi,:) = VVu(myi,:)*Singularuv(ii)
			endif
		enddo
		deallocate(Singularuv)


		mn2=min(blocks_A%N,rank)
		myArows = numroc_wp(blocks_A%N, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(mn2, nbslpk, mycol, 0, npcol)
		call descinit( descsUU_v, blocks_A%N, mn2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		allocate(UUv(myArows,myAcols))
		myArows = numroc_wp(mn2, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
		call descinit( descsVV_v, mn2, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		allocate(VVv(myArows,myAcols))
		allocate(Singularuv(mn1))
		call PSVD_Truncate(blocks_A%N,rank,matV2D,descsMatV2D,UUv,VVv,descsUU_v,descsVV_v,Singularuv,option%tol_rand,rank2,ctxt,flop=flop)
		do ii=1,rank2
			call g2l(ii,rank2,nprow,nbslpk,iproc,myi)
			if(iproc==myrow)then
				VVv(myi,:) = VVv(myi,:)*Singularuv(ii)
			endif
		enddo
		deallocate(Singularuv)


		myArows = numroc_wp(rank1, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rank2, nbslpk, mycol, 0, npcol)
		allocate(mattemp(myArows,myAcols))
		mattemp=0
		call descinit( descsMatSml, rank1, rank2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		call pgemmf90('N','T',rank1,rank2,rank,cone,VVu,1,1,descsVV_u,VVv,1,1,descsVV_v,czero,mattemp,1,1,descsMatSml,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
		myArows = numroc_wp(rank1, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(min(rank1,rank2), nbslpk, mycol, 0, npcol)
		call descinit( descsUUSml, rank1, min(rank1,rank2), nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		allocate(UUsml(myArows,myAcols))
		myArows = numroc_wp(min(rank1,rank2), nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(rank2, nbslpk, mycol, 0, npcol)
		call descinit( descsVVSml, min(rank1,rank2), rank2, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		allocate(VVsml(myArows,myAcols))

		allocate(Singularsml(min(rank1,rank2)))
		call PSVD_Truncate(rank1,rank2,mattemp,descsMatSml,UUsml,VVsml,descsUUSml,descsVVSml,Singularsml,option%tol_rand,ranknew,ctxt,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
		myArows = numroc_wp(blocks_A%M, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
		allocate(matU2Dnew(myArows,myAcols))
		matU2Dnew=0
		call descinit( descsMatU2Dnew, blocks_A%M, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		call pgemmf90('N','N',blocks_A%M,ranknew,rank1,cone,UUu,1,1,descsUU_u,UUsml,1,1,descsUUSml,czero,matU2Dnew,1,1,descsMatU2Dnew,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
		do myj=1,myAcols
			call l2g(myj,mycol,ranknew,npcol,nbslpk,jj)
			matU2Dnew(:,myj)=matU2Dnew(:,myj)*Singularsml(jj)
		enddo


		myArows = numroc_wp(blocks_A%N, nbslpk, myrow, 0, nprow)
		myAcols = numroc_wp(ranknew, nbslpk, mycol, 0, npcol)
		allocate(matV2Dnew(myArows,myAcols))
		matV2Dnew=0
		call descinit( descsMatV2Dnew, blocks_A%N, ranknew, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
		call pgemmf90('N','T',blocks_A%N,ranknew,rank2,cone,UUv,1,1,descsUU_v,VVsml,1,1,descsVVSml,czero,matV2Dnew,1,1,descsMatV2Dnew,flop=flop)
		stats%Flop_Factor = stats%Flop_Factor + flop
		rank = ranknew

		deallocate(mattemp,UUsml,VVsml,Singularsml)
		deallocate(UUu,VVu,UUv,VVv)


#endif


	else
		allocate(matU2Dnew(1,1))
		allocate(matV2Dnew(1,1))
	endif
	call MPI_Bcast(rank,1,MPI_INTEGER,Main_ID,ptree%pgrp(pgno)%Comm,ierr)

	deallocate(blocks_A%ButterflyU%blocks(1)%matrix)
	deallocate(blocks_A%ButterflyV%blocks(1)%matrix)
	blocks_A%rankmax=rank
	blocks_A%rankmin=rank
	allocate(blocks_A%ButterflyU%blocks(1)%matrix(blocks_A%M_loc,rank))
	allocate(blocks_A%ButterflyV%blocks(1)%matrix(blocks_A%N_loc,rank))

	call Redistribute2Dto1D(matU2Dnew,blocks_A%M,0,pgno,blocks_A%ButterflyU%blocks(1)%matrix,blocks_A%M_p,0,pgno,rank,ptree)
	call Redistribute2Dto1D(matV2Dnew,blocks_A%N,0,pgno,blocks_A%ButterflyV%blocks(1)%matrix,blocks_A%N_p,0,pgno,rank,ptree)

	deallocate(matU2D)
	deallocate(matV2D)
	deallocate(Vin_tmp,Vout_tmp,Vinter_tmp)
	deallocate(matU,matV)
	deallocate(matU2Dnew,matV2Dnew)

end subroutine LR_A_minusBDinvC

subroutine BF_inverse_schur_partitionedinverse(ho_bf1,level_c,rowblock,error_inout,option,stats,ptree,msh)

    use BPACK_DEFS
	use MISC_Utilities


    use omp_lib

    implicit none

	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt,ll
    character chara
    real(kind=8) T0
    type(matrixblock),pointer::block_o,block_off1,block_off2
    integer rank_new_max,rank0
	real(kind=8):: rank_new_avr,error
	integer niter
	real(kind=8):: error_inout,rate,err_avr
	integer itermax,ntry
	real(kind=8):: n1,n2,Memory
	type(Hoption)::option
	type(Hstat)::stats
	type(hobf)::ho_bf1
	type(proctree)::ptree
	type(mesh)::msh
	integer pgno

	error_inout=0

	block_off1 => ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2-1)%LL(1)%matrices_block(1)
	block_off2 => ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2)%LL(1)%matrices_block(1)


	block_o => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
	block_o%level_butterfly = block_off1%level_butterfly
	level_butterfly=block_o%level_butterfly

	Memory = 0

	if(block_off1%level_butterfly==0 .or. block_off2%level_butterfly==0)then
		call LR_minusBC(ho_bf1,level_c,rowblock,ptree,stats)
	else
		ho_bf1%ind_lv=level_c
		ho_bf1%ind_bk=rowblock
		rank0 = max(block_off1%rankmax,block_off2%rankmax)
		rate=1.2d0
		call BF_randomized(block_o%pgno,level_butterfly,rank0,rate,block_o,ho_bf1,BF_block_MVP_inverse_minusBC_dat,error,'minusBC',option,stats,ptree,msh)
		stats%Flop_Factor=stats%Flop_Factor+stats%Flop_Tmp
		error_inout = max(error_inout, error)
	endif

	pgno = block_o%pgno

	n1 = OMP_get_wtime()
	! if(block_o%level==3)then
	if(level_butterfly>=option%schulzlevel)then
		call BF_inverse_schulziteration_IplusButter(block_o,error,option,stats,ptree,msh)
	else
		call BF_inverse_partitionedinverse_IplusButter(block_o,level_butterfly,0,option,error,stats,ptree,msh,pgno)
	endif

	error_inout = max(error_inout, error)

	n2 = OMP_get_wtime()
	stats%Time_SMW=stats%Time_SMW + n2-n1
	! write(*,*)'I+B Inversion Time:',n2-n1

	if(ptree%MyID==Main_ID .and. option%verbosity>=1)write(*,'(A10,I5,A6,I3,A8,I3,A11,Es14.7)')'OneL No. ',rowblock,' rank:',block_o%rankmax,' L_butt:',block_o%level_butterfly,' error:',error_inout


    return

end subroutine BF_inverse_schur_partitionedinverse

subroutine BF_inverse_schulziteration_IplusButter(block_o,error_inout,option,stats,ptree,msh)

    use BPACK_DEFS
	use MISC_Utilities


    use omp_lib

    implicit none

	integer level_c,rowblock
    integer groupm,blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt,ll
    character chara
    real(kind=8) T0
    type(matrixblock)::block_o,block_Xn
    integer rank_new_max,rank0,num_vect
	real(kind=8):: rank_new_avr,error
	integer niter
	real(kind=8):: error_inout,rate,err_avr
	integer itermax,ntry,converged
	real(kind=8):: n1,n2,Memory,memory_temp,norm1,norm2
	type(Hoption)::option
	type(Hstat)::stats
	type(proctree)::ptree
	type(mesh)::msh
	type(schulz_operand)::schulz_op
	DT,allocatable::VecIn(:,:),VecOut(:,:),VecBuff(:,:)
	DT::ctemp1,ctemp2
	character(len=10)::iternumber
	integer ierr

	error_inout=0
	level_butterfly=block_o%level_butterfly

	Memory = 0

	mm = block_o%M_loc

	num_vect=1
	allocate(VecIn(mm,num_vect))
	VecIn=0
	allocate(VecOut(mm,num_vect))
	VecOut=0
	allocate(VecBuff(mm,num_vect))
	VecBuff=0

	call BF_copy('N',block_o,schulz_op%matrices_block,memory_temp)
	call BF_copy('N',block_o,block_Xn,memory_temp)


	call BF_compute_schulz_init(schulz_op,option,ptree,stats)


	itermax=100
	converged=0
	! n1 = OMP_get_wtime()
	do ii=1,itermax

		write(iternumber ,  "(I4)") ii

		rank0 = block_Xn%rankmax

		rate=1.2d0
		call BF_randomized(block_Xn%pgno,level_butterfly,rank0,rate,block_Xn,schulz_op,BF_block_MVP_schulz_dat,error,'schulz iter'//TRIM(iternumber),option,stats,ptree,msh,ii)
		stats%Flop_Factor=stats%Flop_Factor+stats%Flop_Tmp

		if(schulz_op%order==2)schulz_op%scale=schulz_op%scale*(2-schulz_op%scale)
		if(schulz_op%order==3)schulz_op%scale=schulz_op%scale*(3 - 3*schulz_op%scale + schulz_op%scale**2d0)

		! test error

		ctemp1=1.0d0 ; ctemp2=0.0d0
		call RandomMat(mm,num_vect,min(mm,num_vect),VecIn,1)
		! XnR
		call BF_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,'N',mm,mm,num_vect,VecIn,VecBuff,ctemp1,ctemp2,ptree,stats,ii+1)


		! AXnR
		call BF_block_MVP_dat(schulz_op%matrices_block,'N',mm,mm,num_vect,VecBuff,VecOut,ctemp1,ctemp2,ptree,stats)
		VecOut = 	VecBuff+VecOut

		norm1 = fnorm(VecOut-VecIn,mm,num_vect)**2d0
		norm2 = fnorm(VecIn,mm,num_vect)**2d0
		call MPI_ALLREDUCE(MPI_IN_PLACE, norm1, 1,MPI_double_precision, MPI_SUM, ptree%pgrp(schulz_op%matrices_block%pgno)%Comm,ierr)
		call MPI_ALLREDUCE(MPI_IN_PLACE, norm2, 1,MPI_double_precision, MPI_SUM, ptree%pgrp(schulz_op%matrices_block%pgno)%Comm,ierr)
		error_inout = sqrt(norm1)/sqrt(norm2)


		if(ptree%MyID==Main_ID .and. option%verbosity>=1)write(*,'(A22,A6,I3,A8,I2,A8,I3,A7,Es14.7)')' Schultz ',' rank:',block_Xn%rankmax,' Iter:',ii,' L_butt:',block_Xn%level_butterfly,' error:',error_inout


		if(error_inout<option%tol_rand)then
			converged=1
			exit
		endif

		if(isnan(error_inout))then
			converged=0
			exit
		endif


	enddo
	! n2 = OMP_get_wtime()


	if(converged==0)then
		write(*,*)'Schulz Iteration does not converge'
		stop
	else
		! write(*,*)'Schulz Iteration Time:',n2-n1
		call BF_delete(block_o,1)
		call BF_get_rank(block_Xn,ptree)
		rank_new_max = block_Xn%rankmax
		call BF_copy_delete(block_Xn,block_o,Memory)
		call BF_delete(schulz_op%matrices_block,1)
		if(allocated(schulz_op%diags))deallocate(schulz_op%diags)
	endif

	deallocate(VecIn)
	deallocate(VecOut)
	deallocate(VecBuff)


    return

end subroutine BF_inverse_schulziteration_IplusButter







subroutine BF_compute_schulz_init(schulz_op,option,ptree,stats)

    use BPACK_DEFS
	use MISC_Utilities


    use omp_lib

    implicit none

    integer level_butterfly
    integer mm, nn, mn,ii
    real(kind=8) T0

	real(kind=8):: error
	integer niter,groupm,groupn
	real(kind=8):: error_inout
	integer num_vect,rank,ranktmp,q,qq
	real(kind=8):: n1,n2,memory_temp,flop
	type(Hoption)::option
	type(schulz_operand)::schulz_op
	real(kind=8), allocatable:: Singular(:)
	DT, allocatable::UU(:,:),VV(:,:),RandVectIn(:,:),RandVectOut(:,:),matrixtmp(:,:),matrixtmp1(:,:)
	type(proctree)::ptree
	type(Hstat)::stats

	stats%Flop_tmp=0

	schulz_op%order=option%schulzorder

	error_inout=0

	level_butterfly=schulz_op%matrices_block%level_butterfly

	mm = schulz_op%matrices_block%M_loc
	nn=mm
	num_vect=min(10,nn)

	allocate(RandVectIn(nn,num_vect))
	allocate(RandVectOut(mm,num_vect))
	RandVectOut=0
	call RandomMat(nn,num_vect,min(nn,num_vect),RandVectIn,1)

	! computation of AR
	call BF_block_MVP_dat(schulz_op%matrices_block,'N',mm,nn,num_vect,RandVectIn,RandVectOut,cone,czero,ptree,stats)
	RandVectOut = RandVectIn+RandVectOut


	! power iteration of order q, the following is prone to roundoff error, see algorithm 4.4 Halko 2010
	q=6
	do qq=1,q
		RandVectOut=conjg(cmplx(RandVectOut,kind=8))

		call BF_block_MVP_dat(schulz_op%matrices_block,'T',mm,nn,num_vect,RandVectOut,RandVectIn,cone,czero,ptree,stats)
		RandVectIn = RandVectOut+RandVectIn

		RandVectIn=conjg(cmplx(RandVectIn,kind=8))

		call BF_block_MVP_dat(schulz_op%matrices_block,'N',mm,nn,num_vect,RandVectIn,RandVectOut,cone,czero,ptree,stats)
		RandVectOut = RandVectIn+RandVectOut

	enddo



	! computation of range Q of AR
	call PComputeRange(schulz_op%matrices_block%M_p,num_vect,RandVectOut,ranktmp,option%tol_Rdetect,ptree,schulz_op%matrices_block%pgno,flop)
	stats%Flop_Tmp = stats%Flop_Tmp + flop


	! computation of B = Q^c*A
	RandVectOut=conjg(cmplx(RandVectOut,kind=8))
	call BF_block_MVP_dat(schulz_op%matrices_block,'T',mm,nn,num_vect,RandVectOut,RandVectIn,cone,czero,ptree,stats)
	RandVectIn =RandVectOut+RandVectIn
	RandVectOut=conjg(cmplx(RandVectOut,kind=8))

	! computation of singular values of B
	mn=min(schulz_op%matrices_block%M,ranktmp)
	allocate(Singular(mn))
	Singular=0
	call PSVDTruncateSigma(schulz_op%matrices_block,RandVectIn,ranktmp,rank,Singular,option,stats,ptree,flop)
	stats%Flop_Tmp = stats%Flop_Tmp + flop
	schulz_op%A2norm=Singular(1)
	deallocate(Singular)


	deallocate(RandVectIn)
	deallocate(RandVectOut)



	! allocate(matrixtmp1(nn,nn))
	! matrixtmp1=0
	! do ii=1,nn
		! matrixtmp1(ii,ii)=1d0
	! enddo
	! allocate(matrixtmp(nn,nn))
	! matrixtmp=0
	! call BF_block_MVP_dat(schulz_op%matrices_block,'N',mm,nn,nn,matrixtmp1,matrixtmp,cone,czero)
	! matrixtmp = matrixtmp+matrixtmp1
	! allocate (UU(nn,nn),VV(nn,nn),Singular(nn))
	! call SVD_Truncate(matrixtmp,nn,nn,nn,UU,VV,Singular,option%tol_comp,rank)
	! write(*,*)Singular(1),schulz_op%A2norm,'nimade'
	! schulz_op%A2norm=Singular(1)
	! deallocate(UU,VV,Singular)
	! deallocate(matrixtmp)
	! deallocate(matrixtmp1)


	stats%Flop_factor = stats%Flop_tmp

end subroutine BF_compute_schulz_init


recursive subroutine BF_inverse_partitionedinverse_IplusButter(blocks_io,level_butterfly_target,recurlevel,option,error_inout,stats,ptree,msh,pgno)

    use BPACK_DEFS
	use MISC_Utilities


    use omp_lib

    implicit none

	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, recurlevel, mm, nn, ii, jj,tt,kk1,kk2,rank,err_cnt
    character chara
    real(kind=8) T0,err_avr
    type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
    type(matrixblock)::blocks_io
    type(matrixblock)::blocks_schur
    integer rank_new_max,rank0
	real(kind=8):: rank_new_avr,error,rate
	integer niter
	real(kind=8):: error_inout
	integer itermax,ntry
	real(kind=8):: n1,n2,Memory
	DT, allocatable::matrix_small(:,:)
	type(Hoption)::option
	type(Hstat)::stats
	integer level_butterfly_target,pgno,pgno1
	type(proctree)::ptree
	type(mesh)::msh
	integer ierr
	type(matrixblock)::partitioned_block

	error_inout=0

	if(blocks_io%level_butterfly==0)then
		call LR_SMW(blocks_io,Memory,ptree,stats,pgno)
		return
    else
		allocate(partitioned_block%sons(2,2))

		blocks_A => partitioned_block%sons(1,1)
		blocks_B => partitioned_block%sons(1,2)
		blocks_C => partitioned_block%sons(2,1)
		blocks_D => partitioned_block%sons(2,2)

		! split into four smaller butterflies
		n1 = OMP_get_wtime()
		call BF_split(blocks_io, partitioned_block,ptree,stats,msh)
		n2 = OMP_get_wtime()
		stats%Time_split = stats%Time_split + n2-n1




		if(IOwnPgrp(ptree,blocks_D%pgno))then

			! partitioned inverse of D
			! level_butterfly=level_butterfly_target-1
			level_butterfly=blocks_D%level_butterfly
			pgno1 = blocks_D%pgno
			call BF_inverse_partitionedinverse_IplusButter(blocks_D,level_butterfly,recurlevel+1,option,error,stats,ptree,msh,pgno1)
			error_inout = max(error_inout, error)

			! construct the schur complement A-BD^-1C
			! level_butterfly = level_butterfly_target-1
			level_butterfly = blocks_A%level_butterfly

			! write(*,*)'A-BDC',level_butterfly,level
			call BF_get_rank_ABCD(partitioned_block,rank0)
			if(level_butterfly==0)then
				call LR_A_minusBDinvC(partitioned_block,ptree,option,stats)
			else
				rate=1.2d0
				call BF_randomized(blocks_A%pgno,level_butterfly,rank0,rate,blocks_A,partitioned_block,BF_block_MVP_inverse_A_minusBDinvC_dat,error,'A-BD^-1C',option,stats,ptree,msh)
				stats%Flop_Factor=stats%Flop_Factor+stats%Flop_Tmp
				error_inout = max(error_inout, error)
			endif

			! write(*,*)'ddd1'
			! partitioned inverse of the schur complement
			! level_butterfly=level_butterfly_target-1
			level_butterfly=blocks_A%level_butterfly
			pgno1 = blocks_D%pgno
			call BF_inverse_partitionedinverse_IplusButter(blocks_A,level_butterfly,recurlevel+1,option,error,stats,ptree,msh,pgno1)
			error_inout = max(error_inout, error)
			call BF_get_rank_ABCD(partitioned_block,rank0)
		else
			rank0=0
		endif
		call MPI_ALLREDUCE(MPI_IN_PLACE, rank0, 1,MPI_integer, MPI_MAX, ptree%pgrp(blocks_io%pgno)%Comm,ierr)
		call MPI_ALLREDUCE(MPI_IN_PLACE, error_inout, 1,MPI_double_precision, MPI_MAX, ptree%pgrp(blocks_io%pgno)%Comm,ierr)


		if(level_butterfly_target<=1)then ! try to use deterministic algorithms for merging four LRs into a bigger LR
			pgno1 = blocks_io%pgno  ! recording pgno as blocks_io%pgno will be set to partitioned_block%sons(1,1)%pgno in LR_ABCDinverse
			blocks_io%level_butterfly=0
			call LR_ABCDinverse(partitioned_block,blocks_io,ptree,stats,option,msh)
			call BF_ReDistribute_Inplace(blocks_io,pgno1,stats,ptree,msh)
		else
			level_butterfly = level_butterfly_target
			rate=1.2d0
			call BF_randomized(blocks_io%pgno,level_butterfly,rank0,rate,blocks_io,partitioned_block,BF_block_MVP_inverse_ABCD_dat,error,'ABCDinverse',option,stats,ptree,msh)
			stats%Flop_Factor=stats%Flop_Factor+stats%Flop_Tmp
			error_inout = max(error_inout, error)
		endif



		! stop

		if(option%verbosity>=2 .and. recurlevel==0 .and. ptree%MyID==Main_ID)write(*,'(A23,A6,I3,A8,I3,A11,Es14.7)')' RecursiveI ',' rank:',blocks_io%rankmax,' L_butt:',blocks_io%level_butterfly,' error:',error_inout

		do ii=1,2
		do jj=1,2
			call BF_delete(partitioned_block%sons(ii,jj),1)
		enddo
		enddo
		deallocate(partitioned_block%sons)

		return

	end if
end subroutine BF_inverse_partitionedinverse_IplusButter




recursive subroutine LR_ABCDMerge(blocks,rank,option,msh,stats,ptree,pgno,gd,cridx)
use BPACK_DEFS
implicit none
    integer rank,ranktmp
    integer header_m, header_n
    integer N,M,i,j,ii,jj,myi,myj,iproc,jproc,rmax,mn
	type(mesh)::msh
	type(Hoption)::option
	type(matrixblock)::blocks,blockc(2)
	type(proctree)::ptree
	integer pgno
	type(grid),pointer::gd
	type(grid),pointer::gdc1,gdc2
	integer:: cridx,info
	DT,allocatable:: UU(:,:), VV(:,:),matU(:,:),matV(:,:),matU1(:,:),matV1(:,:),matU2(:,:),matV2(:,:),tmp(:,:),matU1D(:,:),matV1D(:,:),Vin(:,:),Vout1(:,:),Vout2(:,:),Vinter(:,:),Fullmat(:,:),QQ1(:,:),matU2D(:,:),matV2D(:,:)
	real(kind=8),allocatable::Singular(:)
	integer nsproc1,nsproc2,nprow,npcol,nprow1D,npcol1D,myrow,mycol,nprow1,npcol1,myrow1,mycol1,nprow2,npcol2,myrow2,mycol2,myArows,myAcols,M1,N1,M2,N2,rank1,rank2,ierr
	integer::descsmatU(9),descsmatV(9),descsmatU1(9),descsmatV1(9),descsmatU2(9),descsmatV2(9),descUU(9),descVV(9),descsmatU1c(9),descsmatU2c(9),descsmatV1c(9),descsmatV2c(9),descButterflyV(9),descButterflyU(9),descButterU1D(9),descButterV1D(9),descVin(9),descVout(9),descVinter(9),descFull(9)
	integer dims(6),dims_tmp(6) ! M1,N1,rank1,M2,N2,rank2
	DT:: TEMP(1)
	integer LWORK,mnmax,mnmin,rank_new
	DT,allocatable:: WORK(:)
	real(kind=8),allocatable::RWORK(:),center(:)
	real(kind=8):: rtemp,dist,error,rtemp1,rtemp0,fnorm1,fnorm0,flop
	integer :: nb1Dc, nb1Dr,frow,Dimn,edge_n,edge_m,MyID,Ntest,rmaxc,rmaxr
	integer,allocatable::M_p(:,:),N_p(:,:),mrange(:),nrange(:)
	type(Hstat)::stats
	integer::passflag=0

	rank=0
	blocks%ButterflyU%idx=1
	blocks%ButterflyU%inc=1
	blocks%ButterflyU%nblk_loc=1
	blocks%ButterflyV%idx=1
	blocks%ButterflyV%inc=1
	blocks%ButterflyV%nblk_loc=1

	call blacs_gridinfo(gd%ctxt, nprow, npcol, myrow, mycol)
	if(.not. (associated(blocks%sons)))then !  reach bottom level
		! !!!!!!! check error
	else
		if(allocated(blocks%ButterflyU%blocks(1)%matrix))then  ! no need to do merge as LR is already built in parallel
			rank=blocks%rankmax
			goto 101
		endif
100		if(.not. associated(gd%gdc))then
			gdc1=>gd
			gdc2=>gd
		else
			gdc1=>gd%gdc(1)
			gdc2=>gd%gdc(2)
		endif

		if(myrow/=-1 .and. mycol/=-1)then
			nsproc1 = gdc1%nsprow*gdc1%nspcol
			nsproc2 = gdc2%nsprow*gdc2%nspcol

			dims_tmp(1:3)=0
			call blacs_gridinfo(gdc1%ctxt, nprow1, npcol1, myrow1, mycol1)

			! if(ptree%MyID==31)write(*,*)ptree%MyID,nprow1,npcol1,myrow1,mycol1,'dddd'
			if(nprow1/=-1 .and. npcol1/=-1)then
				call LR_ABCDMerge(blocks%sons(1,1),rank,option,msh,stats,ptree,pgno,gdc1,cridx+1)
				dims_tmp(1)=blocks%sons(1,1)%M
				dims_tmp(2)=blocks%sons(1,1)%N
				dims_tmp(3)=blocks%sons(1,1)%rankmax
			endif

			dims_tmp(4:6)=0
			call blacs_gridinfo(gdc2%ctxt, nprow2, npcol2, myrow2, mycol2)
			! if(ptree%MyID==31)write(*,*)ptree%MyID,nprow2,npcol2,myrow2,mycol2,'dddd2'
			if(nprow2/=-1 .and. npcol2/=-1)then
				call LR_ABCDMerge(blocks%sons(2,1),rank,option,msh,stats,ptree,pgno,gdc2,cridx+1)
				dims_tmp(4)=blocks%sons(2,1)%M
				dims_tmp(5)=blocks%sons(2,1)%N
				dims_tmp(6)=blocks%sons(2,1)%rankmax
			endif
			! write(*,*)ptree%MyID,cridx+1
			call MPI_ALLREDUCE(dims_tmp,dims,6,MPI_INTEGER,&
			 MPI_MAX,gd%Comm,ierr)

			M1=dims(1)
			N1=dims(2)
			rank1=dims(3)
			M2=dims(4)
			N2=dims(5)
			rank2=dims(6)

			if(mod(cridx+1,2)==1)then  ! merge along column dimension

				call assert(M1==M2,'M1/=M2 in column merge')

				myArows = numroc_wp(M1, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank1+rank2, nbslpk, mycol, 0, npcol)
				allocate(matU(myArows,myAcols))
				call descinit( descsmatU, M1, rank1+rank2, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descsmatU')

				myArows = numroc_wp(N1, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank1, nbslpk, mycol, 0, npcol)
				allocate(matV1(myArows,myAcols))
				call descinit( descsmatV1, N1, rank1, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descsmatV1')

				myArows = numroc_wp(N2, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank2, nbslpk, mycol, 0, npcol)
				allocate(matV2(myArows,myAcols))
				call descinit( descsmatV2, N2, rank2, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descsmatV2')



				if(nprow1/=-1 .and. npcol1/=-1)then
					! redistribute U1
					myArows = numroc_wp(M1, nbslpk, myrow1, 0, nprow1)
					myAcols = numroc_wp(rank1, nbslpk, mycol1, 0, npcol1)
					call descinit( descsmatU1c, M1, rank1, nbslpk, nbslpk, 0, 0, gdc1%ctxt, max(myArows,1), info )
					call assert(info==0,'descinit fail for descsmatU1c')

					call pgemr2df90(M1, rank1, blocks%sons(1,1)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU1c, matU, 1, 1, descsmatU, gd%ctxt)
					deallocate(blocks%sons(1,1)%ButterflyU%blocks(1)%matrix)
					deallocate(blocks%sons(1,1)%ButterflyU%blocks)

					! redistribute V1
					myArows = numroc_wp(N1, nbslpk, myrow1, 0, nprow1)
					myAcols = numroc_wp(rank1, nbslpk, mycol1, 0, npcol1)
					call descinit( descsmatV1c, N1, rank1, nbslpk, nbslpk, 0, 0, gdc1%ctxt, max(myArows,1), info )
					call assert(info==0,'descinit fail for descsmatV1c')
					call pgemr2df90(N1, rank1, blocks%sons(1,1)%ButterflyV%blocks(1)%matrix, 1, 1, descsmatV1c, matV1, 1, 1, descsmatV1, gd%ctxt)
					deallocate(blocks%sons(1,1)%ButterflyV%blocks(1)%matrix)
					deallocate(blocks%sons(1,1)%ButterflyV%blocks)
				else
					descsmatU1c(2)=-1
					call pgemr2df90(M1, rank1, tmp, 1, 1, descsmatU1c, matU, 1, 1, descsmatU, gd%ctxt)
					descsmatV1c(2)=-1
					call pgemr2df90(N1, rank1, tmp, 1, 1, descsmatV1c, matV1, 1, 1, descsmatV1, gd%ctxt)
				endif

				if(nprow2/=-1 .and. npcol2/=-1)then
					! redistribute U2
					myArows = numroc_wp(M2, nbslpk, myrow2, 0, nprow2)
					myAcols = numroc_wp(rank2, nbslpk, mycol2, 0, npcol2)
					call descinit( descsmatU2c, M2, rank2, nbslpk, nbslpk, 0, 0, gdc2%ctxt, max(myArows,1), info )
					call assert(info==0,'descinit fail for descsmatU2c')
					call pgemr2df90(M2, rank2, blocks%sons(2,1)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU2c, matU, 1, 1+rank1, descsmatU, gd%ctxt)
					deallocate(blocks%sons(2,1)%ButterflyU%blocks(1)%matrix)
					deallocate(blocks%sons(2,1)%ButterflyU%blocks)

					! redistribute V2
					myArows = numroc_wp(N2, nbslpk, myrow2, 0, nprow2)
					myAcols = numroc_wp(rank2, nbslpk, mycol2, 0, npcol2)
					call descinit( descsmatV2c, N2, rank2, nbslpk, nbslpk, 0, 0, gdc2%ctxt, max(myArows,1), info )
					call assert(info==0,'descinit fail for descsmatV2c')
					call pgemr2df90(N2, rank2, blocks%sons(2,1)%ButterflyV%blocks(1)%matrix, 1, 1, descsmatV2c, matV2, 1, 1, descsmatV2, gd%ctxt)
					deallocate(blocks%sons(2,1)%ButterflyV%blocks(1)%matrix)
					deallocate(blocks%sons(2,1)%ButterflyV%blocks)
				else
					descsmatU2c(2)=-1
					call pgemr2df90(M2, rank2, tmp, 1, 1, descsmatU2c, matU, 1, 1+rank1, descsmatU, gd%ctxt)
					descsmatV2c(2)=-1
					call pgemr2df90(N2, rank2, tmp, 1, 1, descsmatV2c, matV2, 1, 1, descsmatV2, gd%ctxt)
				endif

				! compute truncated SVD on matU
				mnmin=min(M1,rank1+rank2)

				myArows = numroc_wp(M1, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(mnmin, nbslpk, mycol, 0, npcol)
				allocate(UU(myArows,myAcols))
				call descinit( descUU, M1, mnmin, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descUU')
				UU=0
				myArows = numroc_wp(mnmin, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank1+rank2, nbslpk, mycol, 0, npcol)
				allocate(VV(myArows,myAcols))
				call descinit( descVV, mnmin, rank1+rank2, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descVV')
				VV=0

				allocate(Singular(mnmin))
				Singular=0


				call PSVD_Truncate(M1, rank1+rank2,matU,descsmatU,UU,VV,descUU,descVV,Singular,option%tol_comp,rank,gd%ctxt,flop=flop)
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

				do ii=1,rank
					call g2l(ii,rank,nprow,nbslpk,iproc,myi)
					if(iproc==myrow)then
						VV(myi,:) = VV(myi,:)*Singular(ii)
					endif
				enddo



				! compute butterfly U and V
				blocks%rankmax = rank
				blocks%rankmin = rank

				myArows = numroc_wp(blocks%N, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
				allocate(blocks%ButterflyV%blocks(1)%matrix(myArows,myAcols))
				blocks%ButterflyV%blocks(1)%matrix=0

				call descinit( descButterflyV, blocks%N, rank, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descButterflyV')

				myArows = numroc_wp(blocks%M, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
				allocate(blocks%ButterflyU%blocks(1)%matrix(myArows,myAcols))
				call descinit( descButterflyU, blocks%M, rank, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descButterflyU')
				blocks%ButterflyU%blocks(1)%matrix=UU(1:myArows,1:myAcols)



				call pgemmf90('N','T',N1,rank,rank1,cone, matV1,1,1,descsmatV1,VV,1,1,descVV,czero,blocks%ButterflyV%blocks(1)%matrix,1,1,descButterflyV,flop=flop)
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)
				call pgemmf90('N','T',N2,rank,rank2,cone, matV2,1,1,descsmatV2,VV,1,1+rank1,descVV,czero,blocks%ButterflyV%blocks(1)%matrix,1+N1,1,descButterflyV,flop=flop)
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)
				deallocate(UU,VV,Singular,matV1,matV2,matU)
				deallocate(blocks%sons)

			else
				call assert(N1==N2,'N1/=N2 in row merge')
				myArows = numroc_wp(N1, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank1+rank2, nbslpk, mycol, 0, npcol)
				allocate(matV(myArows,myAcols))
				call descinit( descsmatV, N1, rank1+rank2, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descsmatV')

				myArows = numroc_wp(M1, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank1, nbslpk, mycol, 0, npcol)
				allocate(matU1(myArows,myAcols))
				call descinit( descsmatU1, M1, rank1, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descsmatU1')

				myArows = numroc_wp(M2, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank2, nbslpk, mycol, 0, npcol)
				allocate(matU2(myArows,myAcols))
				call descinit( descsmatU2, M2, rank2, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descsmatU2')



				if(nprow1/=-1 .and. npcol1/=-1)then
					! redistribute U1
					myArows = numroc_wp(M1, nbslpk, myrow1, 0, nprow1)
					myAcols = numroc_wp(rank1, nbslpk, mycol1, 0, npcol1)
					call descinit( descsmatU1c, M1, rank1, nbslpk, nbslpk, 0, 0, gdc1%ctxt, max(myArows,1), info )
					call assert(info==0,'descinit fail for descsmatU1c')
					! write(*,*)shape(blocks%sons(1,1)%ButterflyU%blocks(1)%matrix),shape(matU1),rank1,M1,blocks%sons(1,1)%rankmax
					call pgemr2df90(M1, rank1, blocks%sons(1,1)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU1c, matU1, 1, 1, descsmatU1, gd%ctxt)
					deallocate(blocks%sons(1,1)%ButterflyU%blocks(1)%matrix)
					deallocate(blocks%sons(1,1)%ButterflyU%blocks)

					! redistribute V1
					myArows = numroc_wp(N1, nbslpk, myrow1, 0, nprow1)
					myAcols = numroc_wp(rank1, nbslpk, mycol1, 0, npcol1)
					call descinit( descsmatV1c, N1, rank1, nbslpk, nbslpk, 0, 0, gdc1%ctxt, max(myArows,1), info )
					call assert(info==0,'descinit fail for descsmatV1c')
					call pgemr2df90(N1, rank1, blocks%sons(1,1)%ButterflyV%blocks(1)%matrix, 1, 1, descsmatV1c, matV, 1, 1, descsmatV, gd%ctxt)
					deallocate(blocks%sons(1,1)%ButterflyV%blocks(1)%matrix)
					deallocate(blocks%sons(1,1)%ButterflyV%blocks)
				else
					descsmatU1c(2)=-1
					call pgemr2df90(M1, rank1, tmp, 1, 1, descsmatU1c, matU1, 1, 1, descsmatU1, gd%ctxt)
					descsmatV1c(2)=-1
					call pgemr2df90(N1, rank1, tmp, 1, 1, descsmatV1c, matV, 1, 1, descsmatV, gd%ctxt)
				endif

				if(nprow2/=-1 .and. npcol2/=-1)then
					! redistribute U2
					myArows = numroc_wp(M2, nbslpk, myrow2, 0, nprow2)
					myAcols = numroc_wp(rank2, nbslpk, mycol2, 0, npcol2)
					call descinit( descsmatU2c, M2, rank2, nbslpk, nbslpk, 0, 0, gdc2%ctxt, max(myArows,1), info )
					call assert(info==0,'descinit fail for descsmatU2c')
					call pgemr2df90(M2, rank2, blocks%sons(2,1)%ButterflyU%blocks(1)%matrix, 1, 1, descsmatU2c, matU2, 1, 1, descsmatU2, gd%ctxt)
					deallocate(blocks%sons(2,1)%ButterflyU%blocks(1)%matrix)
					deallocate(blocks%sons(2,1)%ButterflyU%blocks)

					! redistribute V2
					myArows = numroc_wp(N2, nbslpk, myrow2, 0, nprow2)
					myAcols = numroc_wp(rank2, nbslpk, mycol2, 0, npcol2)
					call descinit( descsmatV2c, N2, rank2, nbslpk, nbslpk, 0, 0, gdc2%ctxt, max(myArows,1), info )
					call assert(info==0,'descinit fail for descsmatV2c')
					call pgemr2df90(N2, rank2, blocks%sons(2,1)%ButterflyV%blocks(1)%matrix, 1, 1, descsmatV2c, matV, 1, 1+rank1, descsmatV, gd%ctxt)
					deallocate(blocks%sons(2,1)%ButterflyV%blocks(1)%matrix)
					deallocate(blocks%sons(2,1)%ButterflyV%blocks)
				else
					descsmatU2c(2)=-1
					call pgemr2df90(M2, rank2, tmp, 1, 1, descsmatU2c, matU2, 1, 1, descsmatU2, gd%ctxt)
					descsmatV2c(2)=-1
					call pgemr2df90(N2, rank2, tmp, 1, 1, descsmatV2c, matV, 1, 1+rank1, descsmatV, gd%ctxt)
				endif

				! compute truncated SVD on matV
				mnmax=max(N1,rank1+rank2)
				mnmin=min(N1,rank1+rank2)

				myArows = numroc_wp(N1, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(mnmin, nbslpk, mycol, 0, npcol)
				allocate(UU(myArows,myAcols))
				call descinit( descUU, N1, mnmin, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descUU')
				UU=0
				myArows = numroc_wp(mnmin, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank1+rank2, nbslpk, mycol, 0, npcol)
				allocate(VV(myArows,myAcols))
				call descinit( descVV, mnmin, rank1+rank2, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descVV')
				VV=0

				allocate(Singular(mnmin))
				Singular=0

				call PSVD_Truncate(N1, rank1+rank2,matV,descsmatV,UU,VV,descUU,descVV,Singular,option%tol_comp,rank,gd%ctxt,flop=flop)
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)
				do ii=1,rank
					call g2l(ii,rank,nprow,nbslpk,iproc,myi)
					if(iproc==myrow)then
						VV(myi,:) = VV(myi,:)*Singular(ii)
					endif
				enddo

				! compute butterfly U and V
				blocks%rankmax = rank
				blocks%rankmin = rank

				myArows = numroc_wp(blocks%N, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
				allocate(blocks%ButterflyV%blocks(1)%matrix(myArows,myAcols))
				blocks%ButterflyV%blocks(1)%matrix=UU(1:myArows,1:myAcols)

				call descinit( descButterflyV, blocks%N, rank, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descButterflyV')

				myArows = numroc_wp(blocks%M, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
				allocate(blocks%ButterflyU%blocks(1)%matrix(myArows,myAcols))
				call descinit( descButterflyU, blocks%M, rank, nbslpk, nbslpk, 0, 0, gd%ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descButterflyU')
				blocks%ButterflyU%blocks(1)%matrix=0



				call pgemmf90('N','T',M1,rank,rank1,cone, matU1,1,1,descsmatU1,VV,1,1,descVV,czero,blocks%ButterflyU%blocks(1)%matrix,1,1,descButterflyU,flop=flop)
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

				call pgemmf90('N','T',M2,rank,rank2,cone, matU2,1,1,descsmatU2,VV,1,1+rank1,descVV,czero,blocks%ButterflyU%blocks(1)%matrix,1+M1,1,descButterflyU,flop=flop)
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)

				deallocate(blocks%sons)
				deallocate(UU,VV,Singular,matU1,matU2,matV)
			endif
		endif

		! write(*,*)ptree%MyID,cridx,rank,'hei',blocks%M,blocks%N

101			if(cridx==0)then

			!! the following is needed when there is idle procs in the process grids
			ranktmp=rank
			call MPI_ALLREDUCE(ranktmp,rank,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)
			blocks%rankmax = rank
			blocks%rankmin = rank


			! distribute UV factor into 1D grid
			! write(*,*)rank,'wocanide',ptree%MyID
			call blacs_gridinfo(gd%ctxt, nprow, npcol, myrow, mycol)
			if(myrow/=-1 .and. mycol/=-1)then
				myArows = numroc_wp(blocks%M, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
				! if(myArows>0 .and. myAcols>0)then
					allocate(matU2D(myArows,myAcols))
					matU2D = blocks%ButterflyU%blocks(1)%matrix
				! endif
				myArows = numroc_wp(blocks%N, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
				! if(myArows>0 .and. myAcols>0)then
					allocate(matV2D(myArows,myAcols))
					matV2D = blocks%ButterflyV%blocks(1)%matrix
				! endif
			else
				allocate(matU2D(1,1))  ! required for Redistribute2Dto1D
				matU2D=0
				allocate(matV2D(1,1))  ! required for Redistribute2Dto1D
				matV2D=0
			endif

			if(allocated(blocks%ButterflyU%blocks(1)%matrix))deallocate(blocks%ButterflyU%blocks(1)%matrix)
			allocate(blocks%ButterflyU%blocks(1)%matrix(blocks%M_loc,rank))
			blocks%ButterflyU%blocks(1)%matrix=0

			if(allocated(blocks%ButterflyV%blocks(1)%matrix))deallocate(blocks%ButterflyV%blocks(1)%matrix)
			allocate(blocks%ButterflyV%blocks(1)%matrix(blocks%N_loc,rank))
			blocks%ButterflyV%blocks(1)%matrix=0
! write(*,*) blocks%row_group,blocks%col_group,isnanMat(matU2D,size(matU2D,1),size(matU2D,2)),isnanMat(matV2D,size(matV2D,1),size(matV2D,2)),'nima222bi',ptree%MyID
			! write(*,*)size(matU2D,1),size(matU2D,2),ptree%MyID,'dddddddd',myrow,mycol,myArows,myAcols
			call Redistribute2Dto1D(matU2D,blocks%M,0,pgno,blocks%ButterflyU%blocks(1)%matrix,blocks%M_p,0,pgno,rank,ptree)
			call Redistribute2Dto1D(matV2D,blocks%N,0,pgno,blocks%ButterflyV%blocks(1)%matrix,blocks%N_p,0,pgno,rank,ptree)

			! write(*,*) fnorm(blocks%ButterflyU%blocks(1)%matrix,blocks%M_loc,rank),fnorm(blocks%ButterflyV%blocks(1)%matrix,blocks%N_loc,rank),isnanMat(blocks%ButterflyU%blocks(1)%matrix,blocks%M_loc,rank),'nimabi'
			! write(*,*) blocks%row_group,blocks%col_group,isnanMat(blocks%ButterflyU%blocks(1)%matrix,blocks%M_loc,rank),'nimabi',ptree%MyID


			if(allocated(matU2D))deallocate(matU2D)
			if(allocated(matV2D))deallocate(matV2D)

		endif

	endif


end subroutine LR_ABCDMerge



!**** Merge four child LR into a bigger one
	!partitioned_block: (input) partitioned_block%sons store the four smaller BF
	!ptree: (input) process tree
	!stats: (input) statistics
	!option: (input) containing compression options
	!msh: (input) containing mesh
	!blocks_o: (inout) the parent BF
subroutine LR_ABCDinverse(partitioned_block,blocks_o,ptree,stats,option,msh)
use BPACK_DEFS
implicit none
	type(matrixblock)::partitioned_block
	type(matrixblock)::blocks_o
    integer rank,ranktmp,leafsize
    integer header_m, header_n
    integer N,M,i,j,ii,jj,myi,myj,iproc,jproc,mn
	type(mesh)::msh
	type(Hoption)::option
	type(proctree)::ptree
	integer pgno
	type(grid),pointer::gd
	integer:: cridx,info
	integer::mrange_dummy(1),nrange_dummy(1)
	type(Hstat)::stats
	integer::passflag=0
	integer::frow,rmax
	real(kind=8)::error
	DT:: mat_dummy(1,1)


	call BF_delete(blocks_o,1)

	blocks_o%pgno=partitioned_block%sons(1,1)%pgno
	call ComputeParallelIndices(blocks_o,blocks_o%pgno,ptree,msh)

	call LR_BuildABCD(blocks_o,partitioned_block,option,msh,stats,ptree,partitioned_block%sons(1,1)%pgno,ptree%pgrp(partitioned_block%sons(1,1)%pgno)%gd,0)

	call LR_ABCDMerge(blocks_o,rank,option,msh,stats,ptree,partitioned_block%sons(1,1)%pgno,ptree%pgrp(partitioned_block%sons(1,1)%pgno)%gd,0)







end subroutine LR_ABCDinverse




!**** Use low-rank arithmetic to form the four child LRs in inverse_ABCD
	!partitioned_block: (input) partitioned_block%sons store the four smaller LR
	!ptree: (input) process tree
	!stats: (input) statistics
	!option: (input) containing compression options
	!msh: (input) containing mesh
	!blocks: (inout) the current LR in the recursion
	!pgno: (in) the process group used for the four children
	!gd: (in) the process grid from process group pgno
recursive subroutine LR_BuildABCD(blocks,partitioned_block,option,msh,stats,ptree,pgno,gd,cridx)
use BPACK_DEFS
implicit none
	type(matrixblock)::partitioned_block
    integer rank,ranktmp
    integer header_m, header_n
    integer N,M,i,j,ii,jj,myi,myj,iproc,jproc,rmax,mn
	type(mesh)::msh
	type(Hoption)::option
	type(matrixblock)::blocks
	type(proctree)::ptree
	integer pgno
	type(grid),pointer::gd
	type(grid),pointer::gdc1,gdc2
	integer:: cridx,info
	DT,allocatable:: UU(:,:),UU1(:,:),UU2(:,:), VV(:,:),VV1(:,:),VV2(:,:),SS1(:,:),TT1(:,:),matU(:,:),matV(:,:),matU1(:,:),matV1(:,:),matU2(:,:),matV2(:,:),tmp(:,:),matU1D(:,:),matV1D(:,:),Vin(:,:),Vout1(:,:),Vout2(:,:),Vinter(:,:),Fullmat(:,:),QQ1(:,:),matU2D(:,:),matV2D(:,:)
	real(kind=8),allocatable::Singular(:)
	integer nsproc1,nsproc2,nprow,npcol,nprow1D,npcol1D,myrow,mycol,nprow1,npcol1,myrow1,mycol1,nprow2,npcol2,myrow2,mycol2,myArows,myAcols,M1,N1,M2,N2,rank1,rank2,ierr
	integer::descsmatU(9),descsmatV(9),descsmatU1(9),descsmatV1(9),descsmatU2(9),descsmatV2(9),descUU(9),descVV(9),descsmatU1c(9),descsmatU2c(9),descsmatV1c(9),descsmatV2c(9),descButterflyV(9),descButterflyU(9),descButterU1D(9),descButterV1D(9),descVin(9),descVout(9),descVinter(9),descFull(9)
	integer dims(6),dims_tmp(6) ! M1,N1,rank1,M2,N2,rank2
	DT:: TEMP(1)
	integer LWORK,mnmax,mnmin,rank_new
	DT,allocatable:: WORK(:)
	real(kind=8),allocatable::RWORK(:),center(:)
	real(kind=8):: rtemp,dist,error,rtemp1,rtemp0,fnorm1,fnorm0,flop
	integer :: nb1Dc, nb1Dr,frow,Dimn,edge_n,edge_m,MyID,Ntest,rmaxc,rmaxr
	integer,allocatable::M_p(:,:),N_p(:,:),mrange(:),nrange(:)
	type(Hstat)::stats
	integer::passflag=0
	type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D
	integer ctxt


		blocks_A => partitioned_block%sons(1,1)
		blocks_B => partitioned_block%sons(1,2)
		blocks_C => partitioned_block%sons(2,1)
		blocks_D => partitioned_block%sons(2,2)

		rank=0
		blocks%ButterflyU%idx=1
		blocks%ButterflyU%inc=1
		blocks%ButterflyU%nblk_loc=1
		blocks%ButterflyV%idx=1
		blocks%ButterflyV%inc=1
		blocks%ButterflyV%nblk_loc=1
		allocate(blocks%ButterflyU%blocks(1))
		allocate(blocks%ButterflyV%blocks(1))


		if(cridx==2)then ! reach bottom level
			blocks%level = blocks_A%level

			if(blocks%row_group== blocks_A%row_group.and. blocks%col_group==blocks_A%col_group)then
				rank=blocks_A%rankmax
				allocate(UU(blocks_A%M_loc,rank))
				UU=blocks_A%ButterflyU%blocks(1)%matrix
				allocate(VV(blocks_A%N_loc,rank))
				VV=blocks_A%ButterflyV%blocks(1)%matrix
			elseif(blocks%row_group== blocks_B%row_group.and. blocks%col_group==blocks_B%col_group)then
				rank=blocks_B%rankmax
				allocate(UU(blocks_B%M_loc,rank))
				UU=0
				call BF_block_MVP_dat(blocks_A,'N',blocks_A%M_loc,blocks_A%N_loc,rank,blocks_B%ButterflyU%blocks(1)%matrix,UU,cone,czero,ptree,stats)
				UU = -UU - blocks_B%ButterflyU%blocks(1)%matrix
				allocate(VV(blocks_B%N_loc,rank))
				VV=0
				call BF_block_MVP_dat(blocks_D,'T',blocks_D%M_loc,blocks_D%N_loc,rank,blocks_B%ButterflyV%blocks(1)%matrix,VV,cone,czero,ptree,stats)
				VV = VV + blocks_B%ButterflyV%blocks(1)%matrix

			elseif(blocks%row_group== blocks_C%row_group.and. blocks%col_group==blocks_C%col_group)then
				rank=blocks_C%rankmax
				allocate(UU(blocks_C%M_loc,rank))
				UU=0
				call BF_block_MVP_dat(blocks_D,'N',blocks_D%M_loc,blocks_D%N_loc,rank,blocks_C%ButterflyU%blocks(1)%matrix,UU,cone,czero,ptree,stats)
				UU = -UU - blocks_C%ButterflyU%blocks(1)%matrix
				allocate(VV(blocks_C%N_loc,rank))
				VV=0
				call BF_block_MVP_dat(blocks_A,'T',blocks_A%M_loc,blocks_A%N_loc,rank,blocks_C%ButterflyV%blocks(1)%matrix,VV,cone,czero,ptree,stats)
				VV = VV + blocks_C%ButterflyV%blocks(1)%matrix

			elseif(blocks%row_group== blocks_D%row_group.and. blocks%col_group==blocks_D%col_group)then
				rank1=blocks_C%rankmax
				allocate(UU1(blocks_C%M_loc,rank1))
				UU1=0
				call BF_block_MVP_dat(blocks_D,'N',blocks_D%M_loc,blocks_D%N_loc,rank1,blocks_C%ButterflyU%blocks(1)%matrix,UU1,cone,czero,ptree,stats)
				UU1 = UU1 + blocks_C%ButterflyU%blocks(1)%matrix
				allocate(TT1(blocks_A%N_loc,rank1))
				TT1=0
				call BF_block_MVP_dat(blocks_A,'T',blocks_A%M_loc,blocks_A%N_loc,rank1,blocks_C%ButterflyV%blocks(1)%matrix,TT1,cone,czero,ptree,stats)
				TT1 = TT1 + blocks_C%ButterflyV%blocks(1)%matrix
				allocate(SS1(blocks_B%N_loc,rank1))
				SS1=0
				call BF_block_MVP_dat(blocks_B,'T',blocks_B%M_loc,blocks_B%N_loc,rank1,TT1,SS1,cone,czero,ptree,stats)
				allocate(VV1(blocks_B%N_loc,rank1))
				VV1=0
				call BF_block_MVP_dat(blocks_D,'T',blocks_D%M_loc,blocks_D%N_loc,rank1,SS1,VV1,cone,czero,ptree,stats)
				VV1 = VV1 + SS1

				rank2=blocks_D%rankmax
				allocate(UU2(blocks_D%M_loc,rank2))
				UU2 = blocks_D%ButterflyU%blocks(1)%matrix
				allocate(VV2(blocks_D%N_loc,rank2))
				VV2 = blocks_D%ButterflyV%blocks(1)%matrix
				rank=rank1+rank2

				allocate(UU(blocks_D%M_loc,rank))
				UU(:,1:rank1)=UU1
				UU(:,1+rank1:rank)=UU2
				allocate(VV(blocks_D%N_loc,rank))
				VV(:,1:rank1)=VV1
				VV(:,1+rank1:rank)=VV2

				deallocate(UU1,UU2,VV1,VV2,TT1,SS1)

			endif

			blocks%rankmax=rank
			call ComputeParallelIndices(blocks,pgno,ptree,msh)



			!!!!**** generate 2D grid blacs quantities
			ctxt = gd%ctxt
			call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)
			if(myrow/=-1 .and. mycol/=-1)then
				myArows = numroc_wp(blocks%M, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
				call descinit( descsmatU, blocks%M, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
				allocate(blocks%ButterflyU%blocks(1)%matrix(myArows,myAcols))
				blocks%ButterflyU%blocks(1)%matrix=0

				myArows = numroc_wp(blocks%N, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)
				call descinit( descsmatV, blocks%N, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
				allocate(blocks%ButterflyV%blocks(1)%matrix(myArows,myAcols))
				blocks%ButterflyV%blocks(1)%matrix=0
			else
				descsmatU(2)=-1
				descsmatV(2)=-1
				allocate(blocks%ButterflyU%blocks(1)%matrix(1,1))
				allocate(blocks%ButterflyV%blocks(1)%matrix(1,1))
			endif


			!!!!**** redistribution of UV factors
			call Redistribute1Dto2D(UU,blocks%M_p,0,pgno,blocks%ButterflyU%blocks(1)%matrix,blocks%M,0,pgno,rank,ptree)
			call Redistribute1Dto2D(VV,blocks%N_p,0,pgno,blocks%ButterflyV%blocks(1)%matrix,blocks%N,0,pgno,rank,ptree)


			deallocate(UU,VV)

		else
			if(.not. associated(gd%gdc))then
				gdc1=>gd
				gdc2=>gd
			else
				gdc1=>gd%gdc(1)
				gdc2=>gd%gdc(2)
			endif

			call blacs_gridinfo(gd%ctxt, nprow, npcol, myrow, mycol)
			if(myrow/=-1 .and. mycol/=-1)then
				nsproc1 = gdc1%nsprow*gdc1%nspcol
				nsproc2 = gdc2%nsprow*gdc2%nspcol
				allocate(blocks%sons(2,1))
				call blacs_gridinfo(gdc1%ctxt, nprow1, npcol1, myrow1, mycol1)
				if(nprow1/=-1 .and. npcol1/=-1)then

					! proportional mapping along row or column dimensions
					if(mod(cridx+1,2)==1)then  ! split along column dimension
						blocks%sons(1,1)%row_group=blocks%row_group
						blocks%sons(1,1)%headm = blocks%headm
						blocks%sons(1,1)%col_group=blocks%col_group*2
						blocks%sons(1,1)%headn = blocks%headn
					else  ! split along row dimension
						blocks%sons(1,1)%col_group=blocks%col_group
						blocks%sons(1,1)%headn = blocks%headn
						blocks%sons(1,1)%row_group=blocks%row_group*2
						blocks%sons(1,1)%headm = blocks%headm
					endif

					blocks%sons(1,1)%M = msh%basis_group(blocks%sons(1,1)%row_group)%tail-msh%basis_group(blocks%sons(1,1)%row_group)%head+1
					blocks%sons(1,1)%N = msh%basis_group(blocks%sons(1,1)%col_group)%tail-msh%basis_group(blocks%sons(1,1)%col_group)%head+1


					call LR_BuildABCD(blocks%sons(1,1),partitioned_block,option,msh,stats,ptree,pgno,gdc1,cridx+1)
				endif

				call blacs_gridinfo(gdc2%ctxt, nprow2, npcol2, myrow2, mycol2)
				if(nprow2/=-1 .and. npcol2/=-1)then

					! proportional mapping along row or column dimensions
					if(mod(cridx+1,2)==1)then  ! split along column dimension
						blocks%sons(2,1)%row_group=blocks%row_group
						blocks%sons(2,1)%headm = blocks%headm
						blocks%sons(2,1)%col_group=blocks%col_group*2+1
						blocks%sons(2,1)%headn = blocks%headn + blocks%sons(1,1)%N
					else  ! split along row dimension
						blocks%sons(2,1)%col_group=blocks%col_group
						blocks%sons(2,1)%headn = blocks%headn
						blocks%sons(2,1)%row_group=blocks%row_group*2+1
						blocks%sons(2,1)%headm = blocks%headm + blocks%sons(1,1)%M
					endif
					blocks%sons(2,1)%M = msh%basis_group(blocks%sons(2,1)%row_group)%tail-msh%basis_group(blocks%sons(2,1)%row_group)%head+1
					blocks%sons(2,1)%N = msh%basis_group(blocks%sons(2,1)%col_group)%tail-msh%basis_group(blocks%sons(2,1)%col_group)%head+1


					call LR_BuildABCD(blocks%sons(2,1),partitioned_block,option,msh,stats,ptree,pgno,gdc2,cridx+1)
				endif
			endif
		endif

end subroutine LR_BuildABCD








!**** Merge four child butterflies (of the same level) into a bigger one
	!partitioned_block: (input) partitioned_block%sons store the four smaller BF
	!ptree: (input) process tree
	!stats: (input) statistics
	!option: (input) containing compression options
	!msh: (input) containing mesh
	!blocks_o: (inout) the parent BF
subroutine BF_Aggregate(partitioned_block,blocks_o,ptree,stats,option,msh)
    use BPACK_DEFS
	use MISC_Utilities

    implicit none
	integer level_p,ADflag,iii,jjj
	integer mm1,mm2,nn1,nn2,M1,M2,N1,N2,ii,jj,ss,kk,j,i,mm,nn
	integer level_butterfly, num_blocks, level_butterfly_c, level_butterfly_o, level_butterfly_dummy,level_final,num_blocks_c,level,num_col,num_row,num_rowson,num_colson

    type(matrixblock)::partitioned_block
    type(matrixblock)::blocks_o,blocks_dummyL(2),blocks_dummyR(2)
    type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D,blocks
	DT,allocatable:: matrixtemp1(:,:),matrixtemp2(:,:),vin(:,:),vout1(:,:),vout2(:,:)
	DT::ctemp1,ctemp2
	type(proctree)::ptree
	type(Hstat)::stats
	type(Hoption)::option
	type(mesh)::msh
	integer pgno_i,pgno_o,pgno
	real(kind=8), allocatable:: Singular(:)
	DT, allocatable::UU(:,:),VV(:,:)
	integer idx_r,inc_r,nr,idx_c,inc_c,nc,idx,inc,index_i,index_j,mn,nblk_loc,num_blk,rank,row_group,col_group,level_o
	integer ierr
	character mode

	blocks_A=>partitioned_block%sons(1,1)
	blocks_B=>partitioned_block%sons(1,2)
	blocks_C=>partitioned_block%sons(2,1)
	blocks_D=>partitioned_block%sons(2,2)

	pgno_i = blocks_A%pgno
	pgno_o = blocks_o%pgno
	pgno = min(pgno_i,pgno_o)


	if(IOwnPgrp(ptree,pgno_o))then
		level_butterfly_o = blocks_o%level_butterfly
	else
		level_butterfly_o=-1
	endif
	call MPI_ALLREDUCE(MPI_IN_PLACE,level_butterfly_o,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)

	if(IOwnPgrp(ptree,pgno_i))then
		level_butterfly_c = blocks_A%level_butterfly
	else
		level_butterfly_c=-1
	endif
	call MPI_ALLREDUCE(MPI_IN_PLACE,level_butterfly_c,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(pgno)%Comm,ierr)

	call assert(level_butterfly_o==2+level_butterfly_c,'BF_ABCD only supports merging L-2-level BFs into a L-level BF')



	call BF_delete(blocks_o,0)
	if(.not. allocated(blocks_o%ButterflyKerl))then
		allocate(blocks_o%ButterflyKerl(level_butterfly_o))
	endif
	blocks_o%pgno=pgno_i
	blocks_o%level_half = floor_safe(dble(level_butterfly_o)/2d0) ! from outer to inner


	num_blocks=2**level_butterfly_o
	!****** row-wise ordering from right side
	do level=0, blocks_o%level_half

		call GetLocalBlockRange(ptree,blocks_o%pgno,level,level_butterfly_o,idx_r,inc_r,nr,idx_c,inc_c,nc,'R')
		if(level==0)then
			blocks_o%ButterflyV%idx=idx_c
			blocks_o%ButterflyV%inc=inc_c
			blocks_o%ButterflyV%nblk_loc=nc
			blocks_o%ButterflyV%num_blk=num_blocks
		elseif(level==level_butterfly_o+1)then
			blocks_o%ButterflyU%idx=idx_r
			blocks_o%ButterflyU%inc=inc_r
			blocks_o%ButterflyU%nblk_loc=nr
			blocks_o%ButterflyU%num_blk=num_blocks
		else
			blocks_o%ButterflyKerl(level)%num_row=2**level
			blocks_o%ButterflyKerl(level)%num_col=2**(level_butterfly_o-level+1)
			blocks_o%ButterflyKerl(level)%idx_r=idx_r
			blocks_o%ButterflyKerl(level)%inc_r=inc_r
			blocks_o%ButterflyKerl(level)%nr=nr
			blocks_o%ButterflyKerl(level)%idx_c=idx_c*2-1
			blocks_o%ButterflyKerl(level)%inc_c=inc_c
			blocks_o%ButterflyKerl(level)%nc=nc*2
		endif
	enddo

	!****** column-wise ordering from left side
	level_final=blocks_o%level_half+1
	do level=level_butterfly_o+1,level_final, -1

		call GetLocalBlockRange(ptree,blocks_o%pgno,level,level_butterfly_o,idx_r,inc_r,nr,idx_c,inc_c,nc,'C')

		if(level==0)then
			blocks_o%ButterflyV%idx=idx_c
			blocks_o%ButterflyV%inc=inc_c
			blocks_o%ButterflyV%nblk_loc=nc
			blocks_o%ButterflyV%num_blk=num_blocks
		elseif(level==level_butterfly_o+1)then
			blocks_o%ButterflyU%idx=idx_r
			blocks_o%ButterflyU%inc=inc_r
			blocks_o%ButterflyU%nblk_loc=nr
			blocks_o%ButterflyU%num_blk=num_blocks
		else
			blocks_o%ButterflyKerl(level)%num_row=2**level
			blocks_o%ButterflyKerl(level)%num_col=2**(level_butterfly_o-level+1)
			blocks_o%ButterflyKerl(level)%idx_r=idx_r*2-1
			blocks_o%ButterflyKerl(level)%inc_r=inc_r
			blocks_o%ButterflyKerl(level)%nr=nr*2
			blocks_o%ButterflyKerl(level)%idx_c=idx_c
			blocks_o%ButterflyKerl(level)%inc_c=inc_c
			blocks_o%ButterflyKerl(level)%nc=nc
		endif
	enddo



	!!!! level level_butterfly_o+1 and level_butterfly_o only requires flop operations and communication
	if(IOwnPgrp(ptree,pgno_i))then
	level_butterfly_dummy = max(level_butterfly_o-1,0)
	write(*,*)1,'s',ptree%MyID,ptree%pgrp(pgno_i)%nproc,ptree%pgrp(pgno_o)%nproc
	do iii=1,2
		blocks_dummyL(iii)%level_butterfly=level_butterfly_dummy
		blocks_dummyL(iii)%level_half=level_butterfly_dummy-1 ! the value of level_half is only used to guarantee column-wise all2all
		mode='C'
		call GetLocalBlockRange(ptree,pgno_i,level_butterfly_dummy+1,level_butterfly_dummy,idx_r,inc_r,nr,idx_c,inc_c,nc,mode)
		num_blk=2**level_butterfly_dummy
		idx=idx_r
		inc=inc_r
		nblk_loc=nr
		if(nblk_loc>0)then
		! blocks_dummyL(iii)%ButterflyU%idx=idx+(iii-1)*num_blk/2
		blocks_dummyL(iii)%ButterflyU%idx=idx
		blocks_dummyL(iii)%ButterflyU%inc=inc
		blocks_dummyL(iii)%ButterflyU%nblk_loc=nblk_loc
		blocks_dummyL(iii)%ButterflyU%num_blk=num_blk
		blocks_dummyL(iii)%pgno=pgno_i
		allocate(blocks_dummyL(iii)%ButterflyU%blocks(blocks_dummyL(iii)%ButterflyU%nblk_loc))
		endif


		allocate(blocks_dummyL(iii)%ButterflyKerl(1:level_butterfly_dummy)) ! only need one kernel level
		! level_o = level_butterfly_o
		! num_row=2**level_o
		! num_col=2**(level_butterfly_o-level_o+1)
		num_row=2**level_butterfly_dummy
		num_col=2

		call GetLocalBlockRange(ptree,pgno_i,level_butterfly_dummy,level_butterfly_dummy,idx_r,inc_r,nr,idx_c,inc_c,nc,mode)
		if(mode=='R')then
			idx_c=idx_c*2-1
			inc_c=inc_c
			nc=nc*2
		elseif(mode=='C')then
			idx_r=idx_r*2-1
			inc_r=inc_r
			nr=nr*2
		endif


		if(nr>0 .and. nc>0)then
		blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%idx_r=idx_r
		! blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%idx_r=idx_r+(iii-1)*num_row/2
		blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%idx_c=idx_c
		blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%inc_r=inc_r
		blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%inc_c=inc_c
		blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%nr=nr
		blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%nc=nc
		blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%num_row=num_row
		blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%num_col=num_col
		allocate(blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%nr,blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%nc))
		endif

		if(nblk_loc>0 .and. nr>0 .and. nc>0)then
			do ii=1,partitioned_block%sons(iii,1)%ButterflyU%nblk_loc
				nn1=size(partitioned_block%sons(iii,1)%ButterflyU%blocks(ii)%matrix,2)
				nn2=size(partitioned_block%sons(iii,2)%ButterflyU%blocks(ii)%matrix,2)
				mm=size(partitioned_block%sons(iii,1)%ButterflyU%blocks(ii)%matrix,1)
				allocate(matrixtemp1(mm,nn1+nn2))
				index_i = (ii-1)*blocks_dummyL(iii)%ButterflyU%inc+blocks_dummyL(iii)%ButterflyU%idx
				row_group=blocks_o%row_group*2**level_butterfly_o+(index_i*2-1)-1
				mm1=msh%basis_group(row_group)%tail-msh%basis_group(row_group)%head+1
				mm2=mm-mm1
				matrixtemp1(:,1:nn1)=partitioned_block%sons(iii,1)%ButterflyU%blocks(ii)%matrix
				matrixtemp1(:,1+nn1:nn1+nn2)=partitioned_block%sons(iii,2)%ButterflyU%blocks(ii)%matrix

				!!!! low rank for the first half rows
				allocate(matrixtemp2(mm1,nn1+nn2))
				matrixtemp2=matrixtemp1(1:mm1,:)
				mn=min(mm1,nn1+nn2)
				allocate (UU(mm1,mn),VV(mn,nn1+nn2),Singular(mn))
				call SVD_Truncate(matrixtemp2,mm1,nn1+nn2,mn,UU,VV,Singular,option%tol_comp,rank)
				do ss=1,rank
				UU(:,ss)=UU(:,ss)*Singular(ss)
				enddo
				allocate(blocks_dummyL(iii)%ButterflyU%blocks(ii*2-1)%matrix(mm1,rank))
				blocks_dummyL(iii)%ButterflyU%blocks(ii*2-1)%matrix=UU(:,1:rank)

				allocate(blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2-1,1)%matrix(rank,nn1))
				blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2-1,1)%matrix = VV(1:rank,1:nn1)
				allocate(blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2-1,2)%matrix(rank,nn2))
				blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2-1,2)%matrix = VV(1:rank,1+nn1:nn1+nn2)
				deallocate(UU,VV,Singular,matrixtemp2)

				!!!! low rank for the second half rows
				allocate(matrixtemp2(mm2,nn1+nn2))
				matrixtemp2=matrixtemp1(1+mm1:mm,:)
				mn=min(mm2,nn1+nn2)
				allocate (UU(mm2,mn),VV(mn,nn1+nn2),Singular(mn))
				call SVD_Truncate(matrixtemp2,mm2,nn1+nn2,mn,UU,VV,Singular,option%tol_comp,rank)
				do ss=1,rank
				UU(:,ss)=UU(:,ss)*Singular(ss)
				enddo
				allocate(blocks_dummyL(iii)%ButterflyU%blocks(ii*2)%matrix(mm2,rank))
				blocks_dummyL(iii)%ButterflyU%blocks(ii*2)%matrix=UU(:,1:rank)

				allocate(blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2,1)%matrix(rank,nn1))
				blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2,1)%matrix = VV(1:rank,1:nn1)
				allocate(blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2,2)%matrix(rank,nn2))
				blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy)%blocks(ii*2,2)%matrix = VV(1:rank,1+nn1:nn1+nn2)
				deallocate(UU,VV,Singular,matrixtemp2)
				deallocate(matrixtemp1)
			enddo
		endif
		num_blk=2**level_butterfly_dummy
		call BF_all2all_UV(blocks_dummyL(iii),pgno_i,blocks_dummyL(iii)%ButterflyU,level_butterfly_dummy+1,(iii-1)*num_blk,blocks_o,pgno_i,blocks_o%ButterflyU,level_butterfly_o+1,stats,ptree)
		call BF_all2all_ker(blocks_dummyL(iii),pgno_i,blocks_dummyL(iii)%ButterflyKerl(level_butterfly_dummy),level_butterfly_dummy,(iii-1)*num_blk,0,blocks_o,pgno_i,blocks_o%ButterflyKerl(level_butterfly_o),level_butterfly_o,stats,ptree)

	enddo


	do jjj=1,2
		blocks_dummyR(jjj)%level_butterfly=level_butterfly_dummy
		blocks_dummyR(jjj)%level_half=level_butterfly_dummy+1  ! the value of level_half is only used to guarantee row-wise all2all
		mode='R'
		call GetLocalBlockRange(ptree,pgno_i,0,level_butterfly_dummy,idx_r,inc_r,nr,idx_c,inc_c,nc,mode)
		num_blk=2**level_butterfly_dummy
		! num_blk=2**level_butterfly_o
		idx=idx_c
		inc=inc_c
		nblk_loc=nc
		if(nblk_loc>0)then
		! blocks_dummyR(jjj)%ButterflyV%idx=idx+(jjj-1)*num_blk/2
		blocks_dummyR(jjj)%ButterflyV%idx=idx
		blocks_dummyR(jjj)%ButterflyV%inc=inc
		blocks_dummyR(jjj)%ButterflyV%nblk_loc=nblk_loc
		blocks_dummyR(jjj)%ButterflyV%num_blk=num_blk
		blocks_dummyR(jjj)%pgno=pgno_i
		allocate(blocks_dummyR(jjj)%ButterflyV%blocks(blocks_dummyR(jjj)%ButterflyV%nblk_loc))
		endif


		allocate(blocks_dummyR(jjj)%ButterflyKerl(1:level_butterfly_dummy)) ! only need one kernel level
		! level_o = 1
		! num_row=2**level_o
		! num_col=2**(level_butterfly_o-level_o+1)
		num_row=2
		num_col=2**level_butterfly_dummy
		call GetLocalBlockRange(ptree,pgno_i,1,level_butterfly_dummy,idx_r,inc_r,nr,idx_c,inc_c,nc,mode)
		if(mode=='R')then
			idx_c=idx_c*2-1
			inc_c=inc_c
			nc=nc*2
		elseif(mode=='C')then
			idx_r=idx_r*2-1
			inc_r=inc_r
			nr=nr*2
		endif


		if(nr>0 .and. nc>0)then
		blocks_dummyR(jjj)%ButterflyKerl(1)%idx_r=idx_r
		! blocks_dummyR(jjj)%ButterflyKerl(1)%idx_c=idx_c+(jjj-1)*num_row/2
		blocks_dummyR(jjj)%ButterflyKerl(1)%idx_c=idx_c
		blocks_dummyR(jjj)%ButterflyKerl(1)%inc_r=inc_r
		blocks_dummyR(jjj)%ButterflyKerl(1)%inc_c=inc_c
		blocks_dummyR(jjj)%ButterflyKerl(1)%nr=nr
		blocks_dummyR(jjj)%ButterflyKerl(1)%nc=nc
		blocks_dummyR(jjj)%ButterflyKerl(1)%num_row=num_row
		blocks_dummyR(jjj)%ButterflyKerl(1)%num_col=num_col
		allocate(blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(blocks_dummyR(jjj)%ButterflyKerl(1)%nr,blocks_dummyR(jjj)%ButterflyKerl(1)%nc))
		endif

		if(nblk_loc>0 .and. nr>0 .and. nc>0)then
			do ii=1,partitioned_block%sons(1,jjj)%ButterflyV%nblk_loc
				mm1=size(partitioned_block%sons(1,jjj)%ButterflyV%blocks(ii)%matrix,2)
				mm2=size(partitioned_block%sons(2,jjj)%ButterflyV%blocks(ii)%matrix,2)
				nn=size(partitioned_block%sons(1,jjj)%ButterflyV%blocks(ii)%matrix,1)
				allocate(matrixtemp1(mm1+mm2,nn))



				index_j = (ii-1)*blocks_dummyR(jjj)%ButterflyV%inc+blocks_dummyR(jjj)%ButterflyV%idx
				col_group=blocks_o%col_group*2**level_butterfly_o+(index_j*2-1)-1
				nn1=msh%basis_group(col_group)%tail-msh%basis_group(col_group)%head+1
				nn2=nn-nn1
				call copymatT(partitioned_block%sons(1,jjj)%ButterflyV%blocks(ii)%matrix,matrixtemp1(1:mm1,:),nn,mm1)
				call copymatT(partitioned_block%sons(2,jjj)%ButterflyV%blocks(ii)%matrix,matrixtemp1(1+mm1:mm1+mm2,:),nn,mm2)

				!!!! low rank for the first half rows
				allocate(matrixtemp2(mm1+mm2,nn1))
				matrixtemp2=matrixtemp1(:,1:nn1)
				mn=min(mm1+mm2,nn1)
				allocate (UU(mm1+mm2,mn),VV(mn,nn1),Singular(mn))
				call SVD_Truncate(matrixtemp2,mm1+mm2,nn1,mn,UU,VV,Singular,option%tol_comp,rank)
				do ss=1,rank
				VV(ss,:)=VV(ss,:)*Singular(ss)
				enddo
				allocate(blocks_dummyR(jjj)%ButterflyV%blocks(ii*2-1)%matrix(nn1,rank))
				call copymatT(VV(1:rank,:),blocks_dummyR(jjj)%ButterflyV%blocks(ii*2-1)%matrix,rank,nn1)
				allocate(blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(1,ii*2-1)%matrix(mm1,rank))
				blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(1,ii*2-1)%matrix = UU(1:mm1,1:rank)
				allocate(blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(2,ii*2-1)%matrix(mm2,rank))
				blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(2,ii*2-1)%matrix = UU(1+mm1:mm1+mm2,1:rank)
				deallocate(UU,VV,Singular,matrixtemp2)

				!!!! low rank for the second half rows
				allocate(matrixtemp2(mm1+mm2,nn2))
				matrixtemp2=matrixtemp1(:,1+nn1:nn1+nn2)
				mn=min(mm1+mm2,nn2)
				allocate (UU(mm1+mm2,mn),VV(mn,nn2),Singular(mn))
				call SVD_Truncate(matrixtemp2,mm1+mm2,nn2,mn,UU,VV,Singular,option%tol_comp,rank)
				do ss=1,rank
				VV(ss,:)=VV(ss,:)*Singular(ss)
				enddo
				allocate(blocks_dummyR(jjj)%ButterflyV%blocks(ii*2)%matrix(nn2,rank))
				call copymatT(VV(1:rank,:),blocks_dummyR(jjj)%ButterflyV%blocks(ii*2)%matrix,rank,nn2)
				allocate(blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(1,ii*2)%matrix(mm1,rank))
				blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(1,ii*2)%matrix = UU(1:mm1,1:rank)
				allocate(blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(2,ii*2)%matrix(mm2,rank))
				blocks_dummyR(jjj)%ButterflyKerl(1)%blocks(2,ii*2)%matrix = UU(1+mm1:mm1+mm2,1:rank)
				deallocate(UU,VV,Singular,matrixtemp2)
				deallocate(matrixtemp1)
			enddo
		endif
		num_blk=2**level_butterfly_dummy
		call BF_all2all_UV(blocks_dummyR(jjj),pgno_i,blocks_dummyR(jjj)%ButterflyV,0,(jjj-1)*num_blk,blocks_o,pgno_i,blocks_o%ButterflyV,0,stats,ptree)
		call BF_all2all_ker(blocks_dummyR(jjj),pgno_i,blocks_dummyR(jjj)%ButterflyKerl(1),1,0,(jjj-1)*num_blk,blocks_o,pgno_i,blocks_o%ButterflyKerl(1),1,stats,ptree)
	enddo

	do iii=1,2
	do jjj=1,2
	do level=2,blocks_o%level_butterfly-1
		num_row=2**level/2
		num_col=2**(level_butterfly_o-level+1)/2
		call BF_all2all_ker(partitioned_block%sons(iii,jjj),pgno_i,partitioned_block%sons(iii,jjj)%ButterflyKerl(level-1),level-1,(iii-1)*num_row,(jjj-1)*num_col,blocks_o,pgno_i,blocks_o%ButterflyKerl(level),level,stats,ptree)
	enddo
	enddo
	enddo

	call BF_delete(blocks_dummyL(1),1)
	call BF_delete(blocks_dummyL(2),1)
	call BF_delete(blocks_dummyR(1),1)
	call BF_delete(blocks_dummyR(2),1)
	endif

	!!!! redistribute the parent BF to its desired process group
	if(IOwnPgrp(ptree,pgno))then
		call BF_ReDistribute_Inplace(blocks_o,pgno_o,stats,ptree,msh)
	endif

	if(IOwnPgrp(ptree,pgno_o))then
	call BF_get_rank(blocks_o,ptree)
	endif


	! mm1=msh%basis_group(blocks_A%row_group)%tail-msh%basis_group(blocks_A%row_group)%head+1
	! nn1=msh%basis_group(blocks_A%col_group)%tail-msh%basis_group(blocks_A%col_group)%head+1
	! mm2=msh%basis_group(blocks_D%row_group)%tail-msh%basis_group(blocks_D%row_group)%head+1
	! nn2=msh%basis_group(blocks_D%col_group)%tail-msh%basis_group(blocks_D%col_group)%head+1


	! allocate(vin(nn1+nn2,1))
	! vin = 1
	! allocate(vout1(mm1+mm2,1))
	! vout1 = 0
	! allocate(vout2(mm1+mm2,1))
	! vout2 = 0

	! ctemp1 = 1d0; ctemp2 = 0d0
	! call BF_block_MVP_dat(blocks_i,'N',mm1+mm2,nn1+nn2,1,vin,vout1,ctemp1,ctemp2)

	! ctemp1 = 1d0; ctemp2 = 0d0
	! call BF_block_MVP_dat(blocks_A,'N',mm1,nn1,1,vin(1:nn1,:),vout2(1:mm1,:),ctemp1,ctemp2)
	! ctemp1 = 1d0; ctemp2 = 1d0
	! call BF_block_MVP_dat(blocks_B,'N',mm1,nn2,1,vin(1+nn1:nn1+nn2,:),vout2(1:mm1,:),ctemp1,ctemp2)

	! ctemp1 = 1d0; ctemp2 = 0d0
	! call BF_block_MVP_dat(blocks_C,'N',mm2,nn1,1,vin(1:nn1,:),vout2(1+mm1:mm1+mm2,:),ctemp1,ctemp2)
	! ctemp1 = 1d0; ctemp2 = 1d0
	! call BF_block_MVP_dat(blocks_D,'N',mm2,nn2,1,vin(1+nn1:nn1+nn2,:),vout2(1+mm1:mm1+mm2,:),ctemp1,ctemp2)

	! write(*,*)'spliting error:',fnorm(vout1-vout2,mm1+mm2,1)/fnorm(vout1,mm1+mm2,1)
	! deallocate(vin,vout1,vout2)


end subroutine BF_Aggregate
subroutine BF_split(blocks_i,blocks_o,ptree,stats,msh)
    use BPACK_DEFS
	use MISC_Utilities


    use omp_lib

    implicit none
	integer level_p,ADflag,iii,jjj
	integer mm1,mm2,nn1,nn2,M1,M2,N1,N2,ii,jj,kk,j,i,mm,nn
	integer level_butterfly, num_blocks, level_butterfly_c, num_blocks_c,level,num_col,num_row,num_rowson,num_colson

    type(matrixblock),target::blocks_i
    type(matrixblock)::blocks_o,blocks_dummy
    type(matrixblock),pointer::blocks_A,blocks_B,blocks_C,blocks_D,blocks
	DT,allocatable:: matrixtemp1(:,:),matrixtemp2(:,:),vin(:,:),vout1(:,:),vout2(:,:)
	DT::ctemp1,ctemp2
	type(mesh)::msh
	type(proctree)::ptree
	type(Hstat)::stats
	integer pgno


	blocks_A=>blocks_o%sons(1,1)
	blocks_B=>blocks_o%sons(1,2)
	blocks_C=>blocks_o%sons(2,1)
	blocks_D=>blocks_o%sons(2,2)

	if(blocks_i%level_butterfly==0)then
		level_butterfly=0
	else
		level_butterfly=max(blocks_i%level_butterfly-2,0)
	endif


	!*** try to use the same process group as blocks_i
	pgno = blocks_i%pgno
	do while(level_butterfly<ptree%nlevel-GetTreelevel(pgno))
		pgno = pgno*2
	enddo


	do iii=1,2
	do jjj=1,2
		blocks=>blocks_o%sons(iii,jjj)
		blocks%level = blocks_i%level+1
		blocks%row_group = blocks_i%row_group*2+iii-1
		blocks%col_group = blocks_i%col_group*2+jjj-1
		blocks%style = blocks_i%style
		blocks%headm = msh%basis_group(blocks%row_group)%head
		blocks%M = msh%basis_group(blocks%row_group)%tail-msh%basis_group(blocks%row_group)%head+1
		blocks%headn = msh%basis_group(blocks%col_group)%head
		blocks%N = msh%basis_group(blocks%col_group)%tail-msh%basis_group(blocks%col_group)%head+1
		blocks%pgno = pgno
		call ComputeParallelIndices(blocks,pgno,ptree,msh)
	enddo
	enddo

	if(blocks_i%level_butterfly==0)then

		kk = size(blocks_i%ButterflyU%blocks(1)%matrix,2)
		do iii=1,2
		do jjj=1,2
			blocks=>blocks_o%sons(iii,jjj)
			blocks%level_butterfly=0
			blocks%level_half = 0
			allocate(blocks%ButterflyU%blocks(1))
			allocate(blocks%ButterflyV%blocks(1))
			if(IOwnPgrp(ptree,blocks%pgno))then
				allocate(blocks%ButterflyU%blocks(1)%matrix(blocks%M_loc,kk))
				allocate(blocks%ButterflyV%blocks(1)%matrix(blocks%N_loc,kk))
				blocks%rankmax = kk
				blocks%rankmin = kk
				blocks%ButterflyU%nblk_loc=1
				blocks%ButterflyU%inc=1
				blocks%ButterflyU%idx=1
				blocks%ButterflyV%nblk_loc=1
				blocks%ButterflyV%inc=1
				blocks%ButterflyV%idx=1
			endif
			call Redistribute1Dto1D(blocks_i%ButterflyU%blocks(1)%matrix,blocks_i%M_p,blocks_i%headm,blocks_i%pgno,blocks%ButterflyU%blocks(1)%matrix,blocks%M_p,blocks%headm,blocks%pgno,kk,ptree)
			call Redistribute1Dto1D(blocks_i%ButterflyV%blocks(1)%matrix,blocks_i%N_p,blocks_i%headn,blocks_i%pgno,blocks%ButterflyV%blocks(1)%matrix,blocks%N_p,blocks%headn,blocks%pgno,kk,ptree)
		enddo
		enddo

	else

	   !**** first redistribute blocks_i into blocks_dummy%sons of the same butterfly levels
		blocks_dummy%level_butterfly = blocks_i%level_butterfly
		blocks_dummy%level_half = blocks_i%level_half
		allocate(blocks_dummy%sons(2,2))
		do ii=1,2
		do jj=1,2
			allocate(blocks_dummy%sons(ii,jj)%ButterflyKerl(blocks_dummy%level_butterfly))
			blocks_dummy%sons(ii,jj)%level_butterfly=blocks_dummy%level_butterfly
			blocks_dummy%sons(ii,jj)%level_half=blocks_dummy%level_half
		enddo
		enddo


		do level=0,blocks_dummy%level_butterfly+1
			if(level==0)then
				call BF_all2all_V_split(blocks_i,blocks_i%pgno,level,blocks_dummy,blocks_o%sons(1,1)%pgno,level,stats,ptree)
			elseif(level==blocks_i%level_butterfly+1)then
				call BF_all2all_U_split(blocks_i,blocks_i%pgno,level,blocks_dummy,blocks_o%sons(1,1)%pgno,level,stats,ptree)
			else
				call BF_all2all_ker_split(blocks_i,blocks_i%pgno,level,blocks_dummy,blocks_o%sons(1,1)%pgno,level,stats,ptree)
			endif
		enddo
		!**** next convert blocks_dummy%sons into  blocks_o%sons
		call BF_convert_to_smallBF(blocks_dummy,blocks_o,stats,ptree)
		do ii=1,2
		do jj=1,2
			call BF_delete(blocks_dummy%sons(ii,jj),1)
		enddo
		enddo
		deallocate(blocks_dummy%sons)
	endif

	do ii=1,2
	do jj=1,2
		call BF_get_rank(blocks_o%sons(ii,jj),ptree)
	enddo
	enddo

	! mm1=msh%basis_group(blocks_A%row_group)%tail-msh%basis_group(blocks_A%row_group)%head+1
	! nn1=msh%basis_group(blocks_A%col_group)%tail-msh%basis_group(blocks_A%col_group)%head+1
	! mm2=msh%basis_group(blocks_D%row_group)%tail-msh%basis_group(blocks_D%row_group)%head+1
	! nn2=msh%basis_group(blocks_D%col_group)%tail-msh%basis_group(blocks_D%col_group)%head+1


	! allocate(vin(nn1+nn2,1))
	! vin = 1
	! allocate(vout1(mm1+mm2,1))
	! vout1 = 0
	! allocate(vout2(mm1+mm2,1))
	! vout2 = 0

	! ctemp1 = 1d0; ctemp2 = 0d0
	! call BF_block_MVP_dat(blocks_i,'N',mm1+mm2,nn1+nn2,1,vin,vout1,ctemp1,ctemp2)

	! ctemp1 = 1d0; ctemp2 = 0d0
	! call BF_block_MVP_dat(blocks_A,'N',mm1,nn1,1,vin(1:nn1,:),vout2(1:mm1,:),ctemp1,ctemp2)
	! ctemp1 = 1d0; ctemp2 = 1d0
	! call BF_block_MVP_dat(blocks_B,'N',mm1,nn2,1,vin(1+nn1:nn1+nn2,:),vout2(1:mm1,:),ctemp1,ctemp2)

	! ctemp1 = 1d0; ctemp2 = 0d0
	! call BF_block_MVP_dat(blocks_C,'N',mm2,nn1,1,vin(1:nn1,:),vout2(1+mm1:mm1+mm2,:),ctemp1,ctemp2)
	! ctemp1 = 1d0; ctemp2 = 1d0
	! call BF_block_MVP_dat(blocks_D,'N',mm2,nn2,1,vin(1+nn1:nn1+nn2,:),vout2(1+mm1:mm1+mm2,:),ctemp1,ctemp2)

	! write(*,*)'spliting error:',fnorm(vout1-vout2,mm1+mm2,1)/fnorm(vout1,mm1+mm2,1)
	! deallocate(vin,vout1,vout2)


end subroutine BF_split





subroutine BF_get_rank_ABCD(partitioned_block,rankmax)

    use BPACK_DEFS
    implicit none

	integer rankmax,ii,jj
    type(matrixblock)::partitioned_block

	rankmax = -1000
	do ii=1,2
	do jj=1,2
	rankmax = max(rankmax,partitioned_block%sons(ii,jj)%rankmax)
	enddo
	enddo
end subroutine BF_get_rank_ABCD



!**** Update one off-diagonal block in HODLR compressed as
! Bplus/Butterfly/LR by multiplying on it left the inverse of diagonal block
! If LR, call LR_Sblock; if butterfly, call BF_randomized; if Bplus, call Bplus_randomized_constr
	!ho_bf1: working HODLR
	!level_c: level# of the block in HODLR
	!rowblock: block# of the block at this level in HODLR
	!option: containing compression options
	!stats: statistics
	!ptree: process tree
subroutine Bplus_Sblock_randomized_memfree(ho_bf1,level_c,rowblock,option,stats,ptree,msh)

    use BPACK_DEFS
	use omp_lib
    implicit none

	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt
    character chara
    real(kind=8) T0
    type(blockplus),pointer::bplus
	type(matrixblock)::block_old
	type(matrixblock),pointer::block_o
    integer::rank_new_max
	real(kind=8)::rank_new_avr,error,rate,rankrate_inner,rankrate_outter
	integer niter,rank,ntry,rank0,rank0_inner,rank0_outter
	real(kind=8):: error_inout
	real(kind=8):: n1,n2
	type(Hoption)::option
	type(Hstat)::stats
	type(hobf)::ho_bf1
	type(proctree)::ptree
	type(mesh)::msh

	error_inout=0

	call Bplus_copy(ho_bf1%levels(level_c)%BP(rowblock),ho_bf1%levels(level_c)%BP_inverse_update(rowblock))
	!!!!!!! the forward block BP can be deleted if not used in solution phase


    bplus =>  ho_bf1%levels(level_c)%BP_inverse_update(rowblock)
	if(bplus%Lplus==1)then

		block_o =>  ho_bf1%levels(level_c)%BP_inverse_update(rowblock)%LL(1)%matrices_block(1)
		level_butterfly=block_o%level_butterfly

		if(level_butterfly==0)then
			call LR_Sblock(ho_bf1,level_c,rowblock,ptree,stats)
		else
			ho_bf1%ind_lv=level_c
			ho_bf1%ind_bk=rowblock
			rank0 = block_o%rankmax
			rate = 1.2d0
			call BF_randomized(block_o%pgno,level_butterfly,rank0,rate,block_o,ho_bf1,BF_block_MVP_Sblock_dat,error_inout,'Sblock',option,stats,ptree,msh,msh)
			stats%Flop_Factor=stats%Flop_Factor+stats%Flop_Tmp
		end if

		if(ptree%MyID==Main_ID .and. option%verbosity>=1)write(*,'(A10,I5,A6,I3,A8,I3,A11,Es14.7)')'OneL No. ',rowblock,' rank:',block_o%rankmax,' L_butt:',block_o%level_butterfly,' error:',error_inout


	else

		ho_bf1%ind_lv=level_c
		ho_bf1%ind_bk=rowblock
		Bplus =>  ho_bf1%levels(level_c)%BP_inverse_update(rowblock)
		block_o =>  ho_bf1%levels(level_c)%BP_inverse_update(rowblock)%LL(1)%matrices_block(1)

		rank0_inner = Bplus%LL(2)%rankmax
		rankrate_inner = 2.0d0

		rank0_outter = block_o%rankmax
		rankrate_outter=1.2d0
		level_butterfly = block_o%level_butterfly
		call Bplus_randomized_constr(level_butterfly,Bplus,ho_bf1,rank0_inner,rankrate_inner,Bplus_block_MVP_Sblock_dat,rank0_outter,rankrate_outter,Bplus_block_MVP_Outter_Sblock_dat,error_inout,'Sblock+',option,stats,ptree,msh)

		block_o =>  ho_bf1%levels(level_c)%BP_inverse_update(rowblock)%LL(1)%matrices_block(1)

		if(option%verbosity>=1 .and. ptree%myid==ptree%pgrp(Bplus%LL(1)%matrices_block(1)%pgno)%head)write(*,'(A10,I5,A6,I3,A8,I3,A11,Es14.7)')'Mult No. ',rowblock,' rank:',block_o%rankmax,' L_butt:',block_o%level_butterfly,' error:',error_inout


	end if

    return

end subroutine Bplus_Sblock_randomized_memfree



subroutine Bplus_inverse_schur_partitionedinverse(ho_bf1,level_c,rowblock,option,stats,ptree,msh)

    use BPACK_DEFS
	use MISC_Utilities


    use omp_lib
    use Bplus_compress

    implicit none

	integer level_c,rowblock,ierr
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks,level_butterfly_loc
    integer num_col, num_row, level, mm, nn, ii, jj,tt,ll,llplus,bb,mmb
    character chara
    real(kind=8) T0
    type(matrixblock),pointer::block_o,block_off1,block_off2
    type(matrixblock),pointer::blocks_o_D
    type(matrixblock)::block_tmp
	type(blockplus),pointer::Bplus,Bplus_schur
    integer rank_new_max
	real(kind=8):: rank_new_avr,error,err_avr,err_max
	integer niter
	real(kind=8):: error_inout,rate,rankrate_inner,rankrate_outter
	integer itermax,ntry,cnt,cnt_partial
	real(kind=8):: n1,n2,n3,n4,Memory
	integer rank0,rank0_inner,rank0_outter,Lplus,level_BP,levelm,groupm_start,ij_loc,edge_s,edge_e,edge_first,idx_end_m_ref,idx_start_m_ref,idx_start_b,idx_end_b
	DT,allocatable:: matin(:,:),matout(:,:),matin_tmp(:,:),matout_tmp(:,:)
	DT:: ctemp1,ctemp2
	integer, allocatable :: ipiv(:)
	type(Hoption)::option
	type(Hstat)::stats
	type(hobf)::ho_bf1
	type(matrixblock):: agent_block
	type(blockplus):: agent_bplus
	type(proctree) :: ptree
	type(mesh) :: msh


    bplus => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)

	if(bplus%Lplus==1)then
		call BF_inverse_schur_partitionedinverse(ho_bf1,level_c,rowblock,error_inout,option,stats,ptree,msh)
	else
		ctemp1 = 1d0
		ctemp2 = 0d0

		block_off1 => ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2-1)%LL(1)%matrices_block(1)
		block_off2 => ho_bf1%levels(level_c)%BP_inverse_update(rowblock*2)%LL(1)%matrices_block(1)
		block_o => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
	! write(*,*)block_o%row_group,block_o%col_group,level_c,rowblock,block_o%level,'diao'
		block_o%level_butterfly = block_off1%level_butterfly

		Memory = 0

		error_inout=0

		ho_bf1%ind_lv=level_c
		ho_bf1%ind_bk=rowblock

		rank0_inner = ho_bf1%levels(level_c)%BP_inverse_update(2*rowblock-1)%LL(2)%rankmax
		rankrate_inner = 2.0d0

		rank0_outter = max(block_off1%rankmax,block_off2%rankmax)
		rankrate_outter=1.2d0

		level_butterfly = block_o%level_butterfly

		call Bplus_randomized_constr(level_butterfly,Bplus,ho_bf1,rank0_inner,rankrate_inner,Bplus_block_MVP_minusBC_dat,rank0_outter,rankrate_outter,Bplus_block_MVP_Outter_minusBC_dat,error,'mBC+',option,stats,ptree,msh)
		error_inout = max(error_inout, error)


		! write(*,*)'good!!!!'
		! stop
		n1 = OMP_get_wtime()
		Bplus => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)
		Lplus = ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)%Lplus
		do llplus =Lplus,1,-1
			do bb=1,Bplus%LL(llplus)%Nbound
				block_o => Bplus%LL(llplus)%matrices_block(bb)
				if(IOwnPgrp(ptree,block_o%pgno))then
					!!!!! partial update butterflies at level llplus from left B1 = D^-1xB
					if(llplus/=Lplus)then
						rank0 = block_o%rankmax
						rate=1.2d0
						level_butterfly = block_o%level_butterfly
						call BF_randomized(block_o%pgno,level_butterfly,rank0,rate,block_o,Bplus,Bplus_block_MVP_diagBinvB_dat,error,'L update',option,stats,ptree,msh,msh)
						stats%Flop_Factor=stats%Flop_Factor+stats%Flop_Tmp
						error_inout = max(error_inout, error)
					endif

					!!!!! invert I+B1 to be I+B2
					level_butterfly=block_o%level_butterfly
					call BF_inverse_partitionedinverse_IplusButter(block_o,level_butterfly,0,option,error,stats,ptree,msh,block_o%pgno)
					error_inout = max(error_inout, error)


					if(llplus/=Lplus)then
						rank0 = block_o%rankmax
						rate=1.2d0
						level_butterfly = block_o%level_butterfly
						call BF_randomized(block_o%pgno,level_butterfly,rank0,rate,block_o,Bplus,Bplus_block_MVP_BdiagBinv_dat,error,'R update',option,stats,ptree,msh,msh)
						stats%Flop_Factor=stats%Flop_Factor+stats%Flop_Tmp
						error_inout = max(error_inout, error)
					endif
				endif
			end do
		end do
		n2 = OMP_get_wtime()
		stats%Time_SMW=stats%Time_SMW + n2-n1


		do ll=1,Bplus%Lplus
		Bplus%LL(ll)%rankmax=0
		do bb=1,Bplus%LL(ll)%Nbound
			Bplus%LL(ll)%rankmax=max(Bplus%LL(ll)%rankmax,Bplus%LL(ll)%matrices_block(bb)%rankmax)
		enddo
		call MPI_ALLREDUCE(MPI_IN_PLACE,Bplus%LL(ll)%rankmax,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(Bplus%LL(1)%matrices_block(1)%pgno)%Comm,ierr)
		end do

		rank_new_max = 0
		do ll=1,Lplus
			rank_new_max = max(rank_new_max,Bplus%LL(ll)%rankmax)
		end do

		if(option%verbosity>=1 .and. ptree%myid==ptree%pgrp(Bplus%LL(1)%matrices_block(1)%pgno)%head)write(*,'(A10,I5,A6,I3,A8,I3,A11,Es14.7)')'Mult No. ',rowblock,' rank:',rank_new_max,' L_butt:',Bplus%LL(1)%matrices_block(1)%level_butterfly,' error:',error_inout

	endif

    return

end subroutine Bplus_inverse_schur_partitionedinverse



end module Bplus_factor
