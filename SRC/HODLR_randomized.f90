#include "HODLR_config.fi"
module HODLR_randomMVP
use Bplus_randomized 
use HODLR_Solve_Mul
use Bplus_compress


contains

subroutine matvec_user(trans,M,N,num_vect,Vin,Vout,ker)
	
	class(*),pointer :: quant
	integer, INTENT(IN):: M,N,num_vect
	DT::Vin(:,:),Vout(:,:) 
	! type(mesh)::msh
	! type(proctree)::ptree
	type(kernelquant)::ker
	! type(Hstat)::stats
	character trans
	
	procedure(F_MatVec), POINTER :: proc
	proc => ker%FuncMatvec
	call proc(trans,M,N,num_vect,Vin,Vout,ker%QuantApp)

	return
	
end subroutine matvec_user	


subroutine HODLR_randomized(ho_bf1,blackbox_HODLR_MVP,Nloc,Memory,error,option,stats,ker,ptree,msh)
	
	
    use HODLR_DEFS
    implicit none
	real(kind=8):: n1,n2,n3,n4,Memory,error_inout,error_lastiter,Memtmp,tmpfact,error,tmp1,tmp2,norm1,norm2
	integer level_c,level_butterfly,bb,rank_new_max,ii,groupm,groupn,Nloc,rank_max_lastiter,rank_max_lastlevel,rank_pre_max,converged
	type(matrixblock),pointer::block_o,block_ref
	DT,allocatable::Vin(:,:),Vout1(:,:),Vout2(:,:)
	type(matrixblock),allocatable::block_rand(:)
	type(hobf)::ho_bf1
	type(Hoption)::option
	type(Hstat)::stats
	real(kind=8):: time_gemm1
	type(kernelquant)::ker
	procedure(MatVec)::blackbox_HODLR_MVP
	type(proctree)::ptree
	type(mesh)::msh
	integer Bidxs,Bidxe,ierr,tt
	
	if(.not. allocated(stats%rankmax_of_level))allocate (stats%rankmax_of_level(ho_bf1%Maxlevel))
	stats%rankmax_of_level = 0

	rank_max_lastlevel = option%rank0
	
	
	n3 = OMP_get_wtime()
	Memory = 0
	do level_c = 1,ho_bf1%Maxlevel+1
		if(level_c==ho_bf1%Maxlevel+1)then
			call HODLR_randomized_OneL_Fullmat(ho_bf1,blackbox_HODLR_MVP,Nloc,level_c,Memtmp,ker,ptree,stats,msh)
		else
			if(level_c>option%LRlevel)then
				level_butterfly = 0   
			else 
				level_butterfly=int((ho_bf1%Maxlevel-level_c)/2)*2
			endif


			if(level_c/=ho_bf1%Maxlevel+1)then
				Bidxs = ho_bf1%levels(level_c)%Bidxs*2-1
				Bidxe = ho_bf1%levels(level_c)%Bidxe*2
			else
				Bidxs = ho_bf1%levels(level_c)%Bidxs
				Bidxe = ho_bf1%levels(level_c)%Bidxe		
			endif			
			
			converged=0
			rank_max_lastiter = rank_max_lastlevel
			error_lastiter=Bigvalue
			do tt = 1,option%itermax
				rank_pre_max = ceiling_safe(rank_max_lastlevel*option%rankrate**(tt-1))+3
							
				if(level_butterfly==0)then
					n1 = OMP_get_wtime()
					allocate (block_rand(Bidxe-Bidxs+1))
					do bb = Bidxs,Bidxe
					! if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
						groupm=ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)%row_group         ! Note: row_group and col_group interchanged here   
						groupn=ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)%col_group         ! Note: row_group and col_group interchanged here  
						call BF_Init_randomized(level_butterfly,rank_pre_max,groupm,groupn,ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1),block_rand(bb-Bidxs+1),msh,1)						
					! endif
					enddo
					n2 = OMP_get_wtime()
					stats%Time_random(1) = stats%Time_random(1) + n2-n1				

					call HODLR_randomized_OneL_Lowrank(ho_bf1,block_rand,blackbox_HODLR_MVP,Nloc,level_c,rank_pre_max,option,ker,ptree,stats,msh)
				else 

					n1 = OMP_get_wtime()
					allocate (block_rand(Bidxe-Bidxs+1))
					do bb = Bidxs,Bidxe
					! if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
						groupm=ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)%row_group         ! Note: row_group and col_group interchanged here   
						groupn=ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)%col_group         ! Note: row_group and col_group interchanged here  
						call BF_Init_randomized(level_butterfly,rank_pre_max,groupm,groupn,ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1),block_rand(bb-Bidxs+1),msh,0)						
					! endif
					enddo
					n2 = OMP_get_wtime()
					stats%Time_random(1) = stats%Time_random(1) + n2-n1				
				
					n1 = OMP_get_wtime()
					call HODLR_Reconstruction_LL(ho_bf1,block_rand,blackbox_HODLR_MVP,Nloc,level_c,level_butterfly,option,stats,ker,ptree,msh)
					n2 = OMP_get_wtime()
					write(*,*)'reconstructLL: ', n2-n1, 'vecCNT',vecCNT

					n1 = OMP_get_wtime()
					call HODLR_Reconstruction_RR(ho_bf1,block_rand,blackbox_HODLR_MVP,Nloc,level_c,level_butterfly,option,stats,ker,ptree,msh)
					n2 = OMP_get_wtime()
					write(*,*)'reconstructRR: ', n2-n1, 'vecCNT',vecCNT
				end if
				
				
				call HODLR_Test_Error_RR(ho_bf1,block_rand,blackbox_HODLR_MVP,Nloc,level_c,error_inout,ker,ptree,stats,msh)
				
				rank_new_max = 0
				do bb = Bidxs,Bidxe
				if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
					call get_butterfly_minmaxrank(block_rand(bb-Bidxs+1))
					rank_new_max = max(rank_new_max,block_rand(bb-Bidxs+1)%rankmax)
				endif	
				end do
				call MPI_ALLREDUCE(MPI_IN_PLACE ,rank_new_max,1,MPI_INTEGER,MPI_MAX,ptree%Comm,ierr)	
				
				
				if(ptree%MyID==Main_ID)write(*,'(A10,I5,A6,I3,A8,I3, A8,I3,A7,Es14.7,A9,I5)')' Level ',level_c,' rank:',rank_new_max,' Ntrial:',tt,' L_butt:',level_butterfly,' error:',error_inout,' #sample:',rank_pre_max
					
				! !!!!*** terminate if 1. error small enough or 2. error not decreasing or 3. rank not increasing	
				! if(error_inout>option%tol_rand .and. error_inout<error_lastiter .and. ((rank_new_max>rank_max_lastiter .and. tt>1).or.tt==1))then
				
				! !!!!*** terminate if 1. error small enough or 2. error not decreasing or 3. rank smaller than num_vec					
				! if(error_inout>option%tol_rand .and. error_inout<error_lastiter .and. rank_new_max==rank_pre_max)then
				
				!!!!*** terminate if 1. error small enough or 2. rank smaller than num_vec					
				if(error_inout>option%tol_rand .and. rank_new_max==rank_pre_max)then				
					do bb = Bidxs,Bidxe
						call delete_blocks(block_rand(bb-Bidxs+1),1)
					end do
					deallocate(block_rand)
					error_lastiter=	error_inout
					rank_max_lastiter = rank_new_max					
				else 
					do bb = Bidxs,Bidxe
					if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
						block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
						call copy_butterfly('N',block_rand(bb-Bidxs+1),block_o,Memtmp)
						Memory = Memory + Memtmp
					endif	
					end do

					do bb = Bidxs,Bidxe
						call delete_blocks(block_rand(bb-Bidxs+1),1)
					end do					
					deallocate (block_rand)
					stats%rankmax_of_level(level_c) = rank_new_max
					rank_max_lastlevel = rank_new_max
					converged=1
					exit
				endif
				
			end do
			if(converged==0)then
				write(*,*)'randomized scheme not converged. level: ',level_c,' rank:',rank_new_max,' L_butt:',level_butterfly,' error:',error_inout
				stop
			end if
		end if
	end do

	stats%Mem_Comp_for=stats%Mem_Comp_for+Memory
	n4 = OMP_get_wtime()
	stats%Time_Fill = stats%Time_Fill + n4-n3 		
	if(ptree%MyID==Main_ID)write(*,*)  'rankmax_of_level:',stats%rankmax_of_level
	
	allocate(Vin(Nloc,1))
	allocate(Vout1(Nloc,1))
	allocate(Vout2(Nloc,1))
	do ii=1,Nloc
		call random_dp_number(Vin(ii,1))
	end do

	call blackbox_HODLR_MVP('N',Nloc,Nloc,1,Vin,Vout1,ker)
	call MVM_Z_forward('N',Nloc,1,1,ho_bf1%Maxlevel+1,Vin,Vout2,ho_bf1,ptree,stats)
	
	tmp1 = fnorm(Vout2-Vout1,Nloc,1)**2d0
	call MPI_ALLREDUCE(tmp1, norm1, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
	tmp2 = fnorm(Vout1,Nloc,1)**2d0
	call MPI_ALLREDUCE(tmp2, norm2, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
	error = sqrt(norm1)/sqrt(norm2)	
	
	deallocate(Vin,Vout1,Vout2)
	
end subroutine HODLR_randomized



subroutine HODLR_randomized_OneL_Lowrank(ho_bf1,block_rand,blackbox_HODLR_MVP,Nloc,level_c,rmax,option,ker,ptree,stats,msh)
	
	
    use HODLR_DEFS
    implicit none
	real(kind=8):: n1,n2,Memory,error_inout,Memtmp
	integer mn,rankref,level_c,rmax,rmaxloc,level_butterfly,bb,bb1,bb_inv,rank_new_max,rank,num_vect,groupn,groupm,header_n,header_m,tailer_m,tailer_n,ii,jj,k,mm,nn
	type(matrixblock),pointer::block_o,block_ref,block_inv
	DT,allocatable::RandVectTmp(:,:)
	DT, allocatable :: matRcol(:,:),matZRcol(:,:),matRrow(:,:),matZcRrow(:,:),mattmp1(:,:),mattmp2(:,:),matrixtemp(:,:),matrixtemp1(:,:),matrixtempQ(:,:),matrixtempin(:,:),matrixtempout(:,:)
	real(kind=8), allocatable:: Singular(:)
	DT::ctemp1,ctemp2
	DT, allocatable::UU(:,:),VV(:,:)
	integer q,qq,Nloc,pp
	integer,allocatable::perms(:), ranks(:) 
	type(hobf)::ho_bf1
	type(matrixblock)::block_rand(:)
	type(Hoption)::option
	DT,allocatable:: RandVectInR(:,:),RandVectOutR(:,:),RandVectInL(:,:),RandVectOutL(:,:)
	type(kernelquant)::ker
	type(proctree)::ptree
	type(Hstat)::stats
	type(mesh)::msh
	procedure(MatVec)::blackbox_HODLR_MVP
	integer Bidxs,Bidxe,head,tail,idx_start_loc,idx_end_loc,ierr

	
	
	Bidxs = ho_bf1%levels(level_c)%Bidxs*2-1
	Bidxe = ho_bf1%levels(level_c)%Bidxe*2
	
	
	level_butterfly = 0
	
	
	ctemp1=1.0d0 ; ctemp2=0.0d0	
	
	rank_new_max = 0
	
	num_vect = rmax
	allocate(RandVectInR(Nloc,num_vect))
	RandVectInR=0
	allocate(RandVectOutR(Nloc,num_vect))
		
	allocate(ranks(Bidxe-Bidxs+1))
	ranks=rmax
	
	! write(*,*)Bidxs,Bidxe,ptree%MyID,'wocao'
	
	do bb = Bidxs,Bidxe
		if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then	
			block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
   
			mm=block_o%M_loc
			
			pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
			header_m = block_o%M_p(pp,1) + block_o%headm -1			
			
			k=header_m-msh%idxs	
			ranks(bb-Bidxs+1)=min(min(block_o%M,block_o%N),num_vect)
			
			allocate(matrixtemp(mm,num_vect))
			call RandomMat(mm,num_vect,min(mm,num_vect),matrixtemp,1)
			
			! write(*,*)1+k,mm+k,header_m,msh%idxs,msh%idxe,block_o%row_group,block_o%col_group,ptree%MyID
			
			RandVectInR(1+k:mm+k,1:num_vect) = matrixtemp
			deallocate(matrixtemp)
		endif	
	end do
	
	
	call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP,'N', RandVectInR, RandVectOutR, Nloc,level_c,num_vect,ker,ptree,stats,msh)
	
	! power iteration of order q, the following is prone to roundoff error, see algorithm 4.4 Halko 2010
	q=0
	do qq=1,q
		RandVectOutR=conjg(cmplx(RandVectOutR,kind=8))
		call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP,'T', RandVectOutR, RandVectInR, Nloc,level_c,num_vect,ker,ptree,stats,msh)
		RandVectInR=conjg(cmplx(RandVectInR,kind=8))
		call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP,'N', RandVectInR, RandVectOutR, Nloc,level_c,num_vect,ker,ptree,stats,msh)
	enddo
	
	
	! computation of range Q
	do bb_inv = ho_bf1%levels(level_c)%Bidxs,ho_bf1%levels(level_c)%Bidxe
		
		pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(bb_inv)%pgno)%head + 1
		head = ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_p(pp,1) + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%headn -1
		tail = head + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_loc -1
		idx_start_loc = head-msh%idxs+1
		idx_end_loc = tail-msh%idxs+1			
		
		call PComputeRange_twoforward(ho_bf1,level_c,Bidxs,bb_inv,ranks,RandVectOutR(idx_start_loc:idx_end_loc,1:num_vect),option%tol_comp*1D-1,ptree,stats)	
	

		! call MPI_ALLREDUCE(MPI_IN_PLACE ,ranks(bb_inv*2-1-Bidxs+1),1,MPI_INTEGER,MPI_MIN,ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%pgno)%Comm,ierr)
		! call MPI_ALLREDUCE(MPI_IN_PLACE ,ranks(bb_inv*2-Bidxs+1),1,MPI_INTEGER,MPI_MIN,ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%pgno)%Comm,ierr)		
		
		! do bb=bb_inv*2-1,bb_inv*2 	
			! if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then	
				! block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1) 
				! mm=block_o%M_loc

				! pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
				! header_m = block_o%M_p(pp,1) + block_o%headm -1				
				! ! header_m=block_o%headm
				! k=header_m-msh%idxs
				
				! allocate(matrixtemp(mm,ranks(bb-Bidxs+1)))
				! matrixtemp = RandVectOutR(1+k:mm+k,1:ranks(bb-Bidxs+1))
				
				! ! write(*,*)block_o%row_group,block_o%col_group,fnorm(matrixtemp(1:mm,1:ranks(bb-Bidxs+1)),mm,ranks(bb-Bidxs+1)),'AR'				
				
				! ! call ComputeRange(mm,ranks(bb),matrixtemp,rank,1,option%tol_comp*1D-1)	
				! call PComputeRange(block_o%M_p,ranks(bb-Bidxs+1),matrixtemp,rank,option%tol_comp*1D-1,ptree,ho_bf1%levels(level_c)%BP(bb)%pgno)
				
				! ranks(bb-Bidxs+1) = rank
				! RandVectOutR(1+k:mm+k,1:ranks(bb-Bidxs+1)) = matrixtemp(1:mm,1:ranks(bb-Bidxs+1))
				
				! ! write(*,*)block_o%row_group,block_o%col_group,fnorm(matrixtemp(1:mm,1:ranks(bb-Bidxs+1)),mm,ranks(bb-Bidxs+1)),'range',1+k,mm+k,ranks(bb-Bidxs+1)
				
				! deallocate(matrixtemp)
			! endif	
		! end do	
		
		! call MPI_ALLREDUCE(MPI_IN_PLACE ,ranks(bb_inv*2-1-Bidxs+1),1,MPI_INTEGER,MPI_MIN,ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%pgno)%Comm,ierr)
		! call MPI_ALLREDUCE(MPI_IN_PLACE ,ranks(bb_inv*2-Bidxs+1),1,MPI_INTEGER,MPI_MIN,ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%pgno)%Comm,ierr)			
		
		
	end do

	!!!!!!!!!!!!!!!!!!!! need add a redistrubtion here 
	
	! computation of B = Q^c*A
	RandVectOutR=conjg(cmplx(RandVectOutR,kind=8))
	call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP,'T', RandVectOutR, RandVectInR, Nloc,level_c,num_vect,ker,ptree,stats,msh)
	RandVectOutR=conjg(cmplx(RandVectOutR,kind=8))
	
	! computation of SVD of B and LR of A
	
	do bb_inv = ho_bf1%levels(level_c)%Bidxs,ho_bf1%levels(level_c)%Bidxe
		
		pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(bb_inv)%pgno)%head + 1
		head = ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_p(pp,1) + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%headn -1
		tail = head + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_loc -1
		idx_start_loc = head-msh%idxs+1
		idx_end_loc = tail-msh%idxs+1			
		
		call PSVDTruncate_twoforward(ho_bf1,level_c,Bidxs,bb_inv,ranks,RandVectOutR(idx_start_loc:idx_end_loc,1:num_vect),RandVectInR(idx_start_loc:idx_end_loc,1:num_vect),block_rand,option%tol_comp,ptree,stats)
			
		! do bb=bb_inv*2-1,bb_inv*2 	
			! block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
			! block_inv => ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)
			! mm=block_o%M
			! nn=block_o%N

! if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
			! pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
			! header_m = block_o%M_p(pp,1) + block_o%headm -1		
			! allocate(matrixtemp1(nn,ranks(bb-Bidxs+1)))
! endif

! allocate(matrixtempin(idx_end_loc-idx_start_loc+1,1:ranks(bb-Bidxs+1)))
! matrixtempin=RandVectInR(idx_start_loc:idx_end_loc,1:ranks(bb-Bidxs+1))
! call Redistribute1Dto1D(matrixtempin,block_inv%M_p,block_inv%headm,block_inv%pgno,matrixtemp1,block_o%N_p,block_o%headn,block_o%pgno,ranks(bb-Bidxs+1),ptree)  
! deallocate(matrixtempin)
			
! if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
			! allocate(matrixtemp(ranks(bb-Bidxs+1),nn))		
			! call copymatT(matrixtemp1,matrixtemp,nn,ranks(bb-Bidxs+1))
			! mn=min(ranks(bb-Bidxs+1),nn)
			! allocate (UU(ranks(bb-Bidxs+1),mn),VV(mn,nn),Singular(mn))
			! call SVD_Truncate(matrixtemp,ranks(bb-Bidxs+1),nn,mn,UU,VV,Singular,option%tol_comp,rank)				
			! do ii=1,rank
				! UU(:,ii) = UU(:,ii)* Singular(ii)
			! end do	

			! rank_new_max = max(rank_new_max,rank)					
			
			! call delete_blocks(block_rand(bb-Bidxs+1),1)
			
			! block_rand(bb-Bidxs+1)%style = 2
			! block_rand(bb-Bidxs+1)%level_butterfly = 0
			! block_rand(bb-Bidxs+1)%rankmax = rank
			! block_rand(bb-Bidxs+1)%rankmin = rank

			! block_rand(bb-Bidxs+1)%row_group=-1
			! block_rand(bb-Bidxs+1)%col_group=-1		
			
			
			! allocate(block_rand(bb-Bidxs+1)%ButterflyU%blocks(1))
			! allocate(block_rand(bb-Bidxs+1)%ButterflyV%blocks(1))		
			
			! allocate(block_rand(bb-Bidxs+1)%ButterflyV%blocks(1)%matrix(nn,rank))
			! call copymatT(VV(1:rank,1:nn),block_rand(bb-Bidxs+1)%ButterflyV%blocks(1)%matrix,rank,nn)
			! allocate(block_rand(bb-Bidxs+1)%ButterflyU%blocks(1)%matrix(mm,rank))
			! allocate(matrixtempQ(mm,ranks(bb-Bidxs+1)))
			! matrixtempQ=0
! endif
			
! allocate(matrixtempin(idx_end_loc-idx_start_loc+1,1:ranks(bb-Bidxs+1)))
! matrixtempin=0
! matrixtempin = RandVectOutR(idx_start_loc:idx_end_loc,1:ranks(bb-Bidxs+1))

! write(*,*)block_o%row_group,block_o%col_group,fnorm(matrixtempin,idx_end_loc-idx_start_loc+1,ranks(bb-Bidxs+1)),'matrixtempin',idx_start_loc,idx_end_loc,ranks(bb-Bidxs+1)


! call Redistribute1Dto1D(matrixtempin,block_inv%M_p,block_inv%headm,block_inv%pgno,matrixtempQ,block_o%M_p,block_o%headm,block_o%pgno,ranks(bb-Bidxs+1),ptree)  
! deallocate(matrixtempin)			
			
	

! if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then

! write(*,*)block_o%row_group,block_o%col_group,fnorm(matrixtempQ,mm,ranks(bb-Bidxs+1)),'matrixtempQ'

	
			! ! call gemm_omp(matrixtempQ,UU(1:ranks(bb),1:rank),block_rand(bb)%ButterflyU%blocks(1)%matrix,mm,rank,ranks(bb))
			! call gemmf90(matrixtempQ,mm,UU,ranks(bb-Bidxs+1),block_rand(bb-Bidxs+1)%ButterflyU%blocks(1)%matrix,mm,'N','N',mm,rank,ranks(bb-Bidxs+1),cone,czero)
			
			
			! deallocate(matrixtemp,matrixtemp1,matrixtempQ,UU,VV,Singular)
			
! endif			
		! end do
		
	end do		

	deallocate(RandVectOutR,RandVectInR,ranks)


end subroutine HODLR_randomized_OneL_Lowrank


subroutine HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP,trans, VectIn, VectOut, Nloc,level_c,num_vect,ker,ptree,stats,msh)
	
	
    use HODLR_DEFS
    implicit none
	real(kind=8):: n1,n2,Memory,error_inout,Memtmp
	integer Nloc,mn,rankref,level_c,rmax,rmaxloc,bb,rank_new_max,rank,num_vect,groupn,groupm,header_n,header_m,tailer_m,tailer_n,ii,jj,k,mm,nn
	type(matrixblock),pointer::block_o,block_ref
	DT,allocatable::RandVectTmp(:,:)
	DT, allocatable :: matRcol(:,:),matZRcol(:,:),matRrow(:,:),matZcRrow(:,:),mattmp1(:,:),mattmp2(:,:),matrixtemp(:,:)
	real(kind=8), allocatable:: Singular(:)
	DT::ctemp1,ctemp2
	complex (kind=8), allocatable::UU(:,:),VV(:,:)
	integer,allocatable::perms(:) 
	character trans
	DT::VectIn(:,:),VectOut(:,:)
	DT,allocatable:: RandVectIn(:,:),RandVectOut(:,:)
	type(hobf)::ho_bf1
	type(kernelquant)::ker
	type(proctree)::ptree
	type(Hstat)::stats
	type(mesh)::msh
	integer Bidxs,Bidxe,pp
	
	procedure(MatVec)::blackbox_HODLR_MVP
	
	ctemp1=1.0d0 ; ctemp2=0.0d0	
	Memory = 0
	rank_new_max = 0
	
		if(level_c/=ho_bf1%Maxlevel+1)then
			Bidxs = ho_bf1%levels(level_c)%Bidxs*2-1
			Bidxe = ho_bf1%levels(level_c)%Bidxe*2
			
			if(trans=='N')then

				VectOut = 0
				
				allocate(RandVectIn(Nloc,num_vect))
				allocate(RandVectOut(Nloc,num_vect))		
				allocate(RandVectTmp(Nloc,num_vect))
				
		
				! Compute the odd block MVP first
				RandVectIn=0
				do bb =Bidxs,Bidxe
					if(mod(bb,2)==1)then
					if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb+1)%pgno))then
						block_o =>  ho_bf1%levels(level_c)%BP(bb+1)%LL(1)%matrices_block(1)
						mm=block_o%M_loc
						
						pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
						header_m = block_o%M_p(pp,1) + block_o%headm -1	
						k=header_m-msh%idxs	
						RandVectIn(1+k:mm+k,1:num_vect) = VectIn(1+k:mm+k,1:num_vect)
					endif
					endif				
				end do	
				
				call blackbox_HODLR_MVP('N',Nloc,Nloc,num_vect,RandVectIn,RandVectTmp,ker)
				call MVM_Z_forward('N',Nloc,num_vect,1,level_c-1,RandVectIn,RandVectOut,ho_bf1,ptree,stats)
				RandVectOut = RandVectTmp-RandVectOut			
				stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp

				do bb =Bidxs,Bidxe
					if(mod(bb,2)==1)then
					if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then			
						block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
						groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   

						mm=block_o%M_loc
						pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
						header_m = block_o%M_p(pp,1) + block_o%headm -1
						k=header_m-msh%idxs	
						VectOut(1+k:mm+k,1:num_vect) = RandVectOut(1+k:mm+k,1:num_vect)
					endif
					endif
				end do	

				! Compute the even block MVP next
				RandVectIn=0
				do bb =Bidxs,Bidxe
					if(mod(bb,2)==0)then
					if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb-1)%pgno))then
						block_o =>  ho_bf1%levels(level_c)%BP(bb-1)%LL(1)%matrices_block(1)   
						mm=block_o%M_loc
						pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
						header_m = block_o%M_p(pp,1) + block_o%headm -1
						k=header_m-msh%idxs	
						RandVectIn(1+k:mm+k,1:num_vect) = VectIn(1+k:mm+k,1:num_vect)
					endif
					endif				
				end do
				
				call blackbox_HODLR_MVP('N',Nloc,Nloc,num_vect,RandVectIn,RandVectTmp,ker)
				
				call MVM_Z_forward('N',Nloc,num_vect,1,level_c-1,RandVectIn,RandVectOut,ho_bf1,ptree,stats)
				RandVectOut = RandVectTmp-RandVectOut			
				stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp
				
				do bb =Bidxs,Bidxe
					if(mod(bb,2)==0)then
					if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then			
						block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
						groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   

						mm=block_o%M_loc
						pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
						header_m = block_o%M_p(pp,1) + block_o%headm -1
						k=header_m-msh%idxs	
						VectOut(1+k:mm+k,1:num_vect) = RandVectOut(1+k:mm+k,1:num_vect)
					endif
					endif
				end do	
					
				deallocate(RandVectIn)	
				deallocate(RandVectOut)	
				deallocate(RandVectTmp)	
				
			else if(trans=='T')then
				VectOut = 0
				
				allocate(RandVectIn(Nloc,num_vect))
				allocate(RandVectOut(Nloc,num_vect))
				
				allocate(RandVectTmp(Nloc,num_vect))
				
				! Compute the odd block MVP first
				RandVectIn=0
				
				do bb =Bidxs,Bidxe
					if(mod(bb,2)==1)then
					if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
						block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
						groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
						mm=block_o%M_loc
						pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
						header_m = block_o%M_p(pp,1) + block_o%headm -1
						k=header_m-msh%idxs	
						RandVectIn(1+k:mm+k,1:num_vect) = VectIn(1+k:mm+k,1:num_vect)
					endif
					endif
				end do	
			
				call blackbox_HODLR_MVP('T',Nloc,Nloc,num_vect,RandVectIn,RandVectTmp,ker)
				
				call MVM_Z_forward('T',Nloc,num_vect,1,level_c-1,RandVectIn,RandVectOut,ho_bf1,ptree,stats)
				RandVectOut = RandVectTmp-RandVectOut		
				stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp

				do bb =Bidxs,Bidxe
					if(mod(bb,2)==1)then
					if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb+1)%pgno))then						
						block_o =>  ho_bf1%levels(level_c)%BP(bb+1)%LL(1)%matrices_block(1)
						mm=block_o%M_loc
						pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
						header_m = block_o%M_p(pp,1) + block_o%headm -1
						k=header_m-msh%idxs	
						VectOut(1+k:mm+k,1:num_vect)=RandVectOut(1+k:mm+k,1:num_vect)				
					endif
					endif
				end do			

				! Compute the even block MVP next
				RandVectIn=0
				do bb =Bidxs,Bidxe
					if(mod(bb,2)==0)then
					if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
						block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
						groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
						mm=block_o%M_loc
						pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
						header_m = block_o%M_p(pp,1) + block_o%headm -1
						k=header_m-msh%idxs	
						RandVectIn(1+k:mm+k,1:num_vect) = VectIn(1+k:mm+k,1:num_vect)
					endif
					endif
				end do	
			
				call blackbox_HODLR_MVP('T',Nloc,Nloc,num_vect,RandVectIn,RandVectTmp,ker)
				call MVM_Z_forward('T',Nloc,num_vect,1,level_c-1,RandVectIn,RandVectOut,ho_bf1,ptree,stats)
				stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp
				RandVectOut = RandVectTmp-RandVectOut		
			
				do bb =Bidxs,Bidxe
					if(mod(bb,2)==0)then
					if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb-1)%pgno))then						
						block_o =>  ho_bf1%levels(level_c)%BP(bb-1)%LL(1)%matrices_block(1)
						mm=block_o%M_loc
						pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
						header_m = block_o%M_p(pp,1) + block_o%headm -1
						k=header_m-msh%idxs	
						VectOut(1+k:mm+k,1:num_vect)=RandVectOut(1+k:mm+k,1:num_vect)				
					endif
					endif
				end do	
				
				deallocate(RandVectIn)	
				deallocate(RandVectOut)	
				deallocate(RandVectTmp)	
				
			endif			
			
			
			
		else
			Bidxs = ho_bf1%levels(level_c)%Bidxs
			Bidxe = ho_bf1%levels(level_c)%Bidxe
			
			VectOut=0
			allocate(RandVectTmp(Nloc,num_vect))
			! Compute the odd block MVP first
			call blackbox_HODLR_MVP('N',Nloc,Nloc,num_vect,VectIn,RandVectTmp,ker)
			call MVM_Z_forward('N',Nloc,num_vect,1,level_c-1,VectIn,VectOut,ho_bf1,ptree,stats)
			VectOut = RandVectTmp-VectOut		
			stats%Flop_Fill = stats%Flop_Fill + stats%Flop_Tmp
			deallocate(RandVectTmp)
		endif
		

		
end subroutine HODLR_MVP_randomized_OneL



subroutine HODLR_randomized_OneL_Fullmat(ho_bf1,blackbox_HODLR_MVP,N,level_c,Memory,ker,ptree,stats,msh)
	
	
    use HODLR_DEFS
    implicit none
	real(kind=8):: n1,n2,Memory,error_inout,Memtmp
	integer N,rankref,level_c,rmaxloc,level_butterfly,bb,rank_new_max,rank,num_vect,groupn,groupm,header_n,header_m,tailer_m,tailer_n,ii,jj,k,mm,nn
	type(matrixblock),pointer::block_o,block_ref
	DT,allocatable::RandVectTmp(:,:)
	DT, allocatable :: matRcol(:,:),matZRcol(:,:),matRrow(:,:),matZcRrow(:,:),mattmp1(:,:),mattmp2(:,:)
	real(kind=8), allocatable:: Singular(:)
	DT::ctemp1,ctemp2
	type(hobf)::ho_bf1
	DT,allocatable:: RandVectInR(:,:),RandVectOutR(:,:),RandVectInL(:,:),RandVectOutL(:,:)
	type(kernelquant)::ker
	type(proctree)::ptree
	type(Hstat)::stats
	type(mesh)::msh
	integer ierr,tempi
	integer Bidxs,Bidxe,N_unk_loc
	
	procedure(MatVec)::blackbox_HODLR_MVP
	
	
	if(level_c/=ho_bf1%Maxlevel+1)then
		Bidxs = ho_bf1%levels(level_c)%Bidxs*2-1
		Bidxe = ho_bf1%levels(level_c)%Bidxe*2
	else
		Bidxs = ho_bf1%levels(level_c)%Bidxs
		Bidxe = ho_bf1%levels(level_c)%Bidxe		
	endif	
	
	
	ctemp1=1.0d0 ; ctemp2=0.0d0	
	
	Memory = 0
	rank_new_max = 0
	
	num_vect = 0
	
	do bb =Bidxs,Bidxe
		if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
			block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
			nn=block_o%N
			num_vect = max(num_vect,nn)
		endif
	end do	
	call MPI_ALLREDUCE(MPI_IN_PLACE ,num_vect,1,MPI_INTEGER,MPI_MAX,ptree%Comm,ierr)
	
	N_unk_loc = msh%idxe-msh%idxs+1
	
	allocate(RandVectInR(N_unk_loc,num_vect))
	RandVectInR=0	
	allocate(RandVectTmp(N_unk_loc,num_vect))
	allocate(RandVectOutR(N_unk_loc,num_vect))
	
	do bb=ho_bf1%levels(level_c)%Bidxs, ho_bf1%levels(level_c)%Bidxe
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		mm=block_o%M_loc
		header_m=block_o%headm
		k=header_m-msh%idxs	
		do ii=1, mm
			RandVectInR(ii+k,ii)=1d0
		enddo
	end do		
	
	call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP,'N', RandVectInR, RandVectOutR, N_unk_loc,level_c,num_vect,ker,ptree,stats,msh)
	
	
	do bb =Bidxs,Bidxe
		if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
			block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
			groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
			mm=block_o%M_loc
			header_m=block_o%headm
			k=header_m-msh%idxs
			
			nn=block_o%N_loc
			
			block_o%style = 1
			allocate(block_o%fullmat(mm,nn))
			! call copymatN(RandVectOutR(k+1:k+mm,1:nn),block_o%fullmat,mm,nn)
			block_o%fullmat = RandVectOutR(k+1:k+mm,1:nn)
			Memory = Memory + SIZEOF(block_o%fullmat)/1024.0d3
		endif
	end do		
	
	deallocate(RandVectInR,RandVectOutR,RandVectTmp)
	
	
	if(ptree%MyID==Main_ID)write(*,'(A10,I5,A13)')' Level ',level_c,' fullmat done'
	
	
end subroutine HODLR_randomized_OneL_Fullmat
		
	


subroutine HODLR_Reconstruction_LL(ho_bf1,block_rand,blackbox_HODLR_MVP,N,level_c,level_butterfly,option,stats,ker,ptree,msh)
    
    use HODLR_DEFS
    implicit none
	
    integer level_c,rowblock,N
    integer group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj,bb, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub,num_vect_subsub
    DT ctemp, a, b
    character chara
	integer level_right_start,num_col,num_row
    
    ! type(matricesblock), pointer :: blocks
    ! type(RandomBlock), pointer :: random
    type(RandomBlock),allocatable :: vec_rand(:)
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
    real(kind=8)::n1,n2

    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
    type(matrixblock),pointer::block_o

	
    integer::rank_new_max,dimension_rank
	real(kind=8)::rank_new_avr,error 
	DT,allocatable::matrixtmp(:,:)
	integer niter,unique_nth  ! level# where each block is touched only once  
	real(kind=8):: error_inout
	integer,allocatable::perms(:)
	
    type(matrixblock)::block_rand(:)
	type(hobf)::ho_bf1	
	type(Hoption)::option
	type(Hstat)::stats
	type(kernelquant)::ker
	type(proctree)::ptree
	type(mesh)::msh
	
	procedure(MatVec)::blackbox_HODLR_MVP
	
	
	
	
	vecCNT = 0
	
    allocate (vec_rand(ho_bf1%levels(level_c)%N_block_forward))
    

    dimension_rank = block_rand(1)%dimension_rank   ! be careful here
	num_blocks=2**level_butterfly	
	num_vect_subsub= dimension_rank+5 ! be careful with the oversampling factor here
	Nbind = 1	
	num_vect_sub = num_vect_subsub*Nbind
	level_right_start = floor_safe(level_butterfly/2d0)
	
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		allocate(vec_rand(bb)%RandomVectorLL(0:level_butterfly+2)) 
		call BF_Init_RandVect_Empty('T',vec_rand(bb),num_vect_sub,block_rand(bb),stats)	
	end do	
	
	do unique_nth = 0,level_right_start
		Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later	
		Ng = 2**level_butterfly/Nsub
		
		do ii = 1,Nsub/Nbind	
			nth_s = (ii-1)*Nbind+1
			nth_e = ii*Nbind
			
			n1 = OMP_get_wtime()
			call HODLR_Randomized_Vectors_LL(ho_bf1,block_rand,vec_rand,blackbox_HODLR_MVP,N,level_c,level_butterfly,nth_s,nth_e,num_vect_sub,unique_nth,ker,ptree,stats,msh)
			vecCNT = vecCNT + num_vect_sub*2
			
			n2 = OMP_get_wtime()
			! time_getvec = time_getvec + n2-n1
			stats%Time_Random(2) = stats%Time_Random(2) + n2-n1
			
			n1 = OMP_get_wtime()
			do bb=1,ho_bf1%levels(level_c)%N_block_forward
				call BF_Resolving_Butterfly_LL_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth,block_rand(bb),vec_rand(bb),option,stats)
			end do
			n2 = OMP_get_wtime()
			stats%Time_Random(3) = stats%Time_Random(3) + n2-n1
		end do
	end do
	
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		call BF_Delete_RandVect('T',vec_rand(bb),level_butterfly)
	end do
	deallocate(vec_rand)

    return
    
end subroutine HODLR_Reconstruction_LL







subroutine HODLR_Reconstruction_RR(ho_bf1,block_rand,blackbox_HODLR_MVP,N,level_c,level_butterfly,option,stats,ker,ptree,msh)
    
    use HODLR_DEFS
    implicit none
	
    integer level_c,rowblock
    integer N, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj,bb, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub,num_vect_subsub
    DT ctemp, a, b
    character chara
	integer level_right_start,num_col,num_row
    
    ! type(matricesblock), pointer :: blocks
    ! type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
    real(kind=8)::n1,n2

    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
    type(matrixblock),pointer::block_o
    integer::rank_new_max,dimension_rank,level_left_start
	real(kind=8)::rank_new_avr 
	DT,allocatable::matrixtmp(:,:)
	integer niter,unique_nth  ! level# where each block is touched only once  
	real(kind=8):: error_inout,rtemp
	integer,allocatable::perms(:)
	type(RandomBlock),allocatable :: vec_rand(:)
	
	type(hobf)::ho_bf1
	type(matrixblock)::block_rand(:)
	type(Hoption)::option
	type(Hstat)::stats
	type(kernelquant)::ker
	type(proctree)::ptree
	type(mesh)::msh
	
	procedure(MatVec)::blackbox_HODLR_MVP
	
	
	vecCNT = 0
	

    dimension_rank = block_rand(1)%dimension_rank   ! be careful here
	num_blocks=2**level_butterfly	
	num_vect_subsub= dimension_rank+5 ! be careful with the oversampling factor here
	Nbind = 1	
	num_vect_sub = num_vect_subsub*Nbind
	
	allocate (vec_rand(ho_bf1%levels(level_c)%N_block_forward))
	
	level_left_start= floor_safe(level_butterfly/2d0)+1    !  check here later
	
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		allocate (vec_rand(bb)%RandomVectorRR(0:level_butterfly+2)) 
		call BF_Init_RandVect_Empty('N',vec_rand(bb),num_vect_sub,block_rand(bb),stats)
	end do
		
	do unique_nth=level_butterfly+1,level_left_start,-1
		if(mod(level_butterfly,2)==0)then
			Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start)))    !  check here later
		else 
			Nsub = NINT(2*2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start)))
		end if			
		Ng = 2**level_butterfly/Nsub
		
		do ii = 1,Nsub/Nbind	
			nth_s = (ii-1)*Nbind+1
			nth_e = ii*Nbind
			
			n1 = OMP_get_wtime()

			call HODLR_Randomized_Vectors_RR(ho_bf1,block_rand,vec_rand,blackbox_HODLR_MVP,N,level_c,level_butterfly,nth_s,nth_e,num_vect_sub,unique_nth,ker,ptree,stats,msh)

			vecCNT = vecCNT + num_vect_sub*2
			
			n2 = OMP_get_wtime()
			! time_getvec = time_getvec + n2-n1
			stats%Time_Random(2) = stats%Time_Random(2) + n2-n1
			
			n1 = OMP_get_wtime()
			do bb=1,ho_bf1%levels(level_c)%N_block_forward				
				call BF_Resolving_Butterfly_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth,block_rand(bb),vec_rand(bb),option,stats)
			end do
			n2 = OMP_get_wtime()
			stats%Time_Random(3) = stats%Time_Random(3) + n2-n1	
		end do
	end do
	

	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		call BF_Delete_RandVect('N',vec_rand(bb),level_butterfly)
	end do
	deallocate(vec_rand)
	
	! write(*,*)'more cool'	
    return
    
end subroutine HODLR_Reconstruction_RR









subroutine HODLR_Test_Error_RR(ho_bf1,block_rand,blackbox_HODLR_MVP,Nloc,level_c,error,ker,ptree,stats,msh)

    use HODLR_DEFS
    implicit none
    
	integer nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,groupm,groupn,bb,pp
    integer Nloc,mm,nn
    real(kind=8) a,b,c,d, condition_number,norm1_R,norm2_R,norm3_R,norm4_R
    DT ctemp, ctemp1, ctemp2
    
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
	integer Nsub,Ng,num_vect,nth_s,nth_e,level_butterfly,ierr
	integer*8 idx_start
	real(kind=8)::error,tmp1,tmp2,norm1,norm2
	integer level_c,rowblock,dimension_m,header_m,tailer_m,header_n,tailer_n
	DT,allocatable::RandomVectors_Output_ref(:,:)
	type(matrixblock),pointer::block_o
	
	type(hobf)::ho_bf1
	type(matrixblock)::block_rand(:)
	type(vectorsblock),pointer:: RandomVectors_InOutput(:)
	type(kernelquant)::ker
	type(proctree)::ptree
	procedure(MatVec)::blackbox_HODLR_MVP
	type(Hstat)::stats
	type(mesh)::msh
	integer Bidxs,Bidxe,bb_inv,head,tail,idx_start_loc,idx_end_loc
	
	if(level_c/=ho_bf1%Maxlevel+1)then
		Bidxs = ho_bf1%levels(level_c)%Bidxs*2-1
		Bidxe = ho_bf1%levels(level_c)%Bidxe*2
	else
		Bidxs = ho_bf1%levels(level_c)%Bidxs
		Bidxe = ho_bf1%levels(level_c)%Bidxe		
	endif

	
	num_vect=1
	allocate (RandomVectors_InOutput(3))
	do ii=1,3
		allocate (RandomVectors_InOutput(ii)%vector(Nloc,num_vect))
		RandomVectors_InOutput(ii)%vector = 0d0
	end do
	
	do bb =Bidxs,Bidxe
		if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
			block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
			mm=block_o%M_loc
			pp = ptree%MyID - ptree%pgrp(block_o%pgno)%head + 1
			header_m = block_o%M_p(pp,1) + block_o%headm -1	
			k=header_m-msh%idxs	
			do jj=1,num_vect
				do ii=1, mm
					call random_dp_number(RandomVectors_InOutput(1)%vector(ii+k,jj))	! matrixtemp1(jj,ii) ! 
				enddo
			enddo
			
			! !!! this is needed to use MVM_Z_forward, otherwise I need to write a seperate subroutine for matvec of block_rand at one level
			! write(*,*)block_rand(bb-Bidxs+1)%row_group,block_rand(bb-Bidxs+1)%col_group,'b copy',bb-Bidxs+1,ptree%MyID
			! call copy_butterfly('N',block_rand(bb-Bidxs+1),block_o)
			! write(*,*)block_o%row_group,block_o%col_group,'a copy'
		endif
	end do	
	

	call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP,'N', RandomVectors_InOutput(1)%vector, RandomVectors_InOutput(3)%vector, Nloc,level_c,num_vect,ker,ptree,stats,msh)
	
	
	do bb_inv = ho_bf1%levels(level_c)%Bidxs,ho_bf1%levels(level_c)%Bidxe	
		pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level_c)%BP_inverse(bb_inv)%pgno)%head + 1
		head = ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_p(pp,1) + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%headn -1
		tail = head + ho_bf1%levels(level_c)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)%N_loc -1
		idx_start_loc = head-msh%idxs+1
		idx_end_loc = tail-msh%idxs+1					
		if(level_c==ho_bf1%Maxlevel+1)then	
			call fullmat_block_MVP_dat(block_rand(bb_inv-Bidxs+1),'N',idx_end_loc-idx_start_loc+1,num_vect,&
			&RandomVectors_InOutput(1)%vector(idx_start_loc:idx_end_loc,1:num_vect),RandomVectors_InOutput(2)%vector(idx_start_loc:idx_end_loc,1:num_vect),cone,czero)			
		else
			call BF_block_MVP_twoforward_dat(ho_bf1,level_c,bb_inv,block_rand,'N',idx_end_loc-idx_start_loc+1,num_vect,RandomVectors_InOutput(1)%vector(idx_start_loc:idx_end_loc,1:num_vect),RandomVectors_InOutput(2)%vector(idx_start_loc:idx_end_loc,1:num_vect),cone,czero,ptree,stats)
		endif	
	end do		
	

	! call MVM_Z_forward('N',Nloc,num_vect,level_c,level_c,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,ho_bf1,ptree,stats)
	! do bb =Bidxs,Bidxe
		! if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(bb)%pgno))then
			! block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)			
			! !!! this is needed to use MVM_Z_forward, otherwise I need to write a seperate subroutine for matvec of block_rand at one level
			! ! write(*,*)block_o%row_group,block_o%col_group,'b delete'
			! call delete_blocks(block_o,1)
			! ! write(*,*)block_o%row_group,block_o%col_group,'a delete'
		! endif
	! end do	
	
	
	
	tmp1 = fnorm(RandomVectors_InOutput(2)%vector-RandomVectors_InOutput(3)%vector,Nloc,1)**2d0
	call MPI_ALLREDUCE(tmp1, norm1, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
	tmp2 = fnorm(RandomVectors_InOutput(3)%vector,Nloc,1)**2d0
	call MPI_ALLREDUCE(tmp2, norm2, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
	error = sqrt(norm1)/sqrt(norm2)		
	! error=0
	do ii=1,3
		deallocate (RandomVectors_InOutput(ii)%vector)
	end do		
	deallocate(RandomVectors_InOutput)
		

    return                

end subroutine HODLR_Test_Error_RR



subroutine HODLR_Randomized_Vectors_LL(ho_bf1,block_rand,vec_rand,blackbox_HODLR_MVP,N,level_c,level_butterfly,nth_s,nth_e,num_vect_sub,unique_nth,ker,ptree,stats,msh)

    use HODLR_DEFS
    
	use misc
    implicit none
    
	integer level_c,rowblock,unique_nth
    integer N,i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,bb
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real(kind=8) a,b,c,d
    ! DT ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o
	
    type(vectorsblock), pointer :: random1, random2
    
    real(kind=8),allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	DT,allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:)
	
	integer Nsub,Ng
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,level_right_start
	! type(RandomBlock), pointer :: random
	real(kind=8)::n1,n2
	
	type(hobf)::ho_bf1
	type(matrixblock)::block_rand(:)
	type(RandomBlock):: vec_rand(:)
	type(vectorsblock),pointer:: RandomVectors_InOutput(:)
	type(kernelquant)::ker
	type(proctree)::ptree
	type(Hstat)::stats
	type(mesh)::msh
	
	procedure(MatVec)::blackbox_HODLR_MVP
	
	
	
    ! ctemp1=1.0d0 ; ctemp2=0.0d0	
	! block_o =>  ho_bf1%levels(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 

	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	

    num_blocks=2**level_butterfly
	level_right_start = floor_safe(level_butterfly/2d0)
	Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later		
    Ng = 2**level_butterfly/Nsub
	! dimension_rank =block_rand(1)%dimension_rank     
	
	
	allocate (RandomVectors_InOutput(3))
	do ii=1,3
		allocate (RandomVectors_InOutput(ii)%vector(N,num_vect_sub))
		RandomVectors_InOutput(ii)%vector = 0d0
	end do
	
	call assert(ho_bf1%levels(level_c)%N_block_forward>1,'N_block_forward not correct')
	! call assert(oddeven==0 .or. oddeven==1,'oddeven not correct')
	
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		! mm=msh%basis_group(groupm)%tail-msh%basis_group(groupm)%head+1 
	
		groupm_start=groupm*2**(level_butterfly)
		
		do nth= nth_s,nth_e
			do i=1, num_blocks
				if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
					header_m=msh%basis_group(groupm_start+i-1)%head
					tailer_m=msh%basis_group(groupm_start+i-1)%tail
					mm=tailer_m-header_m+1
					k=header_m-1	

					allocate(matrixtemp1(num_vect_subsub,mm))
					call RandomMat(num_vect_subsub,mm,min(mm,num_vect_subsub),matrixtemp1,0)
					
					! !$omp parallel do default(shared) private(ii,jj)
					 do jj=1,num_vect_subsub
						 do ii=1, mm
							 call random_dp_number(RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj))	! matrixtemp1(jj,ii) ! 
						 enddo
					 enddo
					 ! !$omp end parallel do
					 deallocate(matrixtemp1)
					 
				 end if
			end do
		end do		
	
	end do
	
	
	call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP,'T', RandomVectors_InOutput(1)%vector, RandomVectors_InOutput(3)%vector, N,level_c,num_vect_sub,ker,ptree,stats,msh)
	
	
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		! mm=msh%basis_group(groupm)%tail-msh%basis_group(groupm)%head+1 
		groupm_start=groupm*2**(level_butterfly)
		! random=>random_Block(bb)
		do i=1, num_blocks
			header_m=msh%basis_group(groupm_start+i-1)%head
			k = header_m - 1
			mm=size(block_rand(bb)%ButterflyU%blocks(i)%matrix,1)
			! !$omp parallel do default(shared) private(ii,jj)
			do ii=1, mm
				do jj=1, num_vect_sub
					vec_rand(bb)%RandomVectorLL(0)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
				enddo
			enddo
			! !$omp end parallel do
			! k=k+mm
		enddo 
		
		
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		! nn=msh%basis_group(groupn)%tail-msh%basis_group(groupn)%head+1 
		groupn_start=groupn*2**(level_butterfly)
		do i=1, num_blocks
			header_n=msh%basis_group(groupn_start+i-1)%head
			k = header_n - 1
			nn=size(block_rand(bb)%ButterflyV%blocks(i)%matrix,1)
			! !$omp parallel do default(shared) private(ii,jj)
			do ii=1, nn
				do jj=1, num_vect_sub
					vec_rand(bb)%RandomVectorLL(level_butterfly+2)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
				enddo
			enddo
			! !$omp end parallel do
			! k=k+nn
		enddo 			
	end do


    ! !$omp parallel do default(shared) private(i)
    do i=1, 3
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    ! !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine HODLR_Randomized_Vectors_LL





subroutine HODLR_Randomized_Vectors_RR(ho_bf1,block_rand,vec_rand,blackbox_HODLR_MVP,N,level_c,level_butterfly,nth_s,nth_e,num_vect_sub,unique_nth,ker,ptree,stats,msh)

    use HODLR_DEFS
    
	use misc
    implicit none
    
	integer level_c,rowblock,unique_nth
    integer N,i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,bb
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real(kind=8) a,b,c,d
    ! DT ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o
	
    ! type(vectorsblock), pointer :: random1, random2
    
    real(kind=8),allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	DT,allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:)
	
	integer Nsub,Ng
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,level_left_start
	! type(RandomBlock), pointer :: random
	real(kind=8)::n1,n2

	type(hobf)::ho_bf1
	type(matrixblock)::block_rand(:)
	type(RandomBlock):: vec_rand(:)
	type(vectorsblock),pointer:: RandomVectors_InOutput(:)
	type(kernelquant)::ker
	procedure(MatVec)::blackbox_HODLR_MVP
	type(proctree)::ptree	
	type(Hstat)::stats
	type(mesh)::msh
	
    ! ctemp1=1.0d0 ; ctemp2=0.0d0	
	! block_o =>  ho_bf1%levels(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 

	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
    num_blocks=2**level_butterfly
	level_left_start= floor_safe(level_butterfly/2d0)+1
	
	if(mod(level_butterfly,2)==0)then
		Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start)))    !  check here later
	else 
		Nsub = NINT(2*2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start)))
	end if	
	Ng = 2**level_butterfly/Nsub    
	
	allocate (RandomVectors_InOutput(3))
	do ii=1,3
		allocate (RandomVectors_InOutput(ii)%vector(N,num_vect_sub))
		RandomVectors_InOutput(ii)%vector = 0d0
	end do
	
	call assert(ho_bf1%levels(level_c)%N_block_forward>1,'N_block_forward not correct')
	! call assert(oddeven==0 .or. oddeven==1,'oddeven not correct')
	
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		! nn=msh%basis_group(groupn)%tail-msh%basis_group(groupn)%head+1 
	
		groupn_start=groupn*2**(level_butterfly)
		
		do nth= nth_s,nth_e
			do i=1, num_blocks
				if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
					header_n=msh%basis_group(groupn_start+i-1)%head
					tailer_n=msh%basis_group(groupn_start+i-1)%tail
					nn=tailer_n-header_n+1
					k=header_n-1	

					allocate(matrixtemp1(num_vect_subsub,nn))
					call RandomMat(num_vect_subsub,nn,min(num_vect_subsub,nn),matrixtemp1,0)
					
					! !$omp parallel do default(shared) private(ii,jj)
					 do jj=1,num_vect_subsub
						 do ii=1, nn
							 call random_dp_number(RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj))	! matrixtemp1(jj,ii) ! 
						 enddo
					 enddo
					 ! !$omp end parallel do
					 deallocate(matrixtemp1)
					 
				 end if
			end do
		end do		
	
	end do
	
	

	call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP,'N', RandomVectors_InOutput(1)%vector, RandomVectors_InOutput(3)%vector, N,level_c,num_vect_sub,ker,ptree,stats,msh)
	
	
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		! nn=msh%basis_group(groupn)%tail-msh%basis_group(groupn)%head+1 
		groupn_start=groupn*2**(level_butterfly)
		! random=>random_Block(bb)
		do i=1, num_blocks
			header_n=msh%basis_group(groupn_start+i-1)%head
			k = header_n - 1		
			nn=size(block_rand(bb)%ButterflyV%blocks(i)%matrix,1)
			! !$omp parallel do default(shared) private(ii,jj)
			do ii=1, nn
				do jj=1, num_vect_sub
					vec_rand(bb)%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
				enddo
			enddo
			! !$omp end parallel do
			! k=k+nn
		enddo 
		
		
		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		! mm=msh%basis_group(groupm)%tail-msh%basis_group(groupm)%head+1 
		groupm_start=groupm*2**(level_butterfly)
		do i=1, num_blocks
			header_m=msh%basis_group(groupm_start+i-1)%head
			k = header_m - 1		
			mm=size(block_rand(bb)%ButterflyU%blocks(i)%matrix,1)
			! !$omp parallel do default(shared) private(ii,jj)
			do ii=1, mm
				do jj=1, num_vect_sub
					vec_rand(bb)%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
				enddo
			enddo
			! !$omp end parallel do
			! k=k+mm
		enddo 			
	end do


    ! !$omp parallel do default(shared) private(i)
    do i=1, 3
        deallocate (RandomVectors_InOutput(i)%vector)
    enddo
    ! !$omp end parallel do
    deallocate (RandomVectors_InOutput)		
	
	
    return                

end subroutine HODLR_Randomized_Vectors_RR

subroutine PComputeRange_twoforward(ho_bf1,level,Bidxs,ii,ranks,AR,eps,ptree,stats)
   use HODLR_DEFS
   implicit none
   integer ranks(:)
   integer level, ii, bb
   DT :: AR(:,:)
   DT,pointer :: matrixtemp1(:,:),matrixtemp2(:,:),matrixtemp(:,:)
   type(matrixblock),pointer::block_o,block_inv,block_schur,block_off1,block_off2,block_off
   integer tempi,groupn,groupm,mm(2),mm1,nn1,mm2,nn2,ierr,nin1,nout1,nin2,nout2,offout(2),offout1,offout2,rank
   type(hobf)::ho_bf1
   type(proctree)::ptree
   type(Hstat)::stats
   real(kind=8)::n1,n2,eps,flop
   integer,pointer::M_p(:,:)
   integer Bidxs
	
	block_off1 => ho_bf1%levels(level)%BP(ii*2-1)%LL(1)%matrices_block(1)	
	block_off2 => ho_bf1%levels(level)%BP(ii*2)%LL(1)%matrices_block(1)
	block_inv => ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)	

	mm(1)=block_off1%M_loc  
	mm(2)=block_off2%M_loc		
	
	offout(1) = 0
	offout(2) = block_off1%M
	
	!!!*** redistribute AR from process layout of hodlr to the process layout of block_off1 and block_off2
	call MPI_ALLREDUCE(MPI_IN_PLACE ,ranks(ii*2-1-Bidxs+1),1,MPI_INTEGER,MPI_MIN,ptree%pgrp(block_inv%pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE ,ranks(ii*2-Bidxs+1),1,MPI_INTEGER,MPI_MIN,ptree%pgrp(block_inv%pgno)%Comm,ierr)
	
	do bb=1,2
	block_off => ho_bf1%levels(level)%BP(ii*2-1+bb-1)%LL(1)%matrices_block(1)
	M_p => block_off%M_p
	if(mm(bb)>0)then
		if(bb==1)then
		allocate(matrixtemp1(mm(bb),ranks(ii*2-1+bb-1-Bidxs+1)))
		matrixtemp=>matrixtemp1
		endif
		if(bb==2)then
		allocate(matrixtemp2(mm(bb),ranks(ii*2-1+bb-1-Bidxs+1)))
		matrixtemp=>matrixtemp2
		endif		
	endif	
	n1 = OMP_get_wtime()
	call Redistribute1Dto1D(AR,block_inv%M_p,0,block_inv%pgno,matrixtemp,M_p,offout(bb),block_off%pgno,ranks(ii*2-1+bb-1-Bidxs+1),ptree)
    n2 = OMP_get_wtime()
    stats%Time_RedistV = stats%Time_RedistV + n2-n1		
	enddo
	
	!!!*** compute range of AR from QR for block_off1 and block_off2
	do bb=1,2
		block_off => ho_bf1%levels(level)%BP(ii*2-1+bb-1)%LL(1)%matrices_block(1)
		M_p => block_off%M_p
		if(bb==1)matrixtemp=>matrixtemp1
		if(bb==2)matrixtemp=>matrixtemp2
			
		if(mm(bb)>0)then
			call PComputeRange(M_p,ranks(ii*2-1+bb-1-Bidxs+1),matrixtemp,rank,eps,ptree,block_off%pgno,flop)
			! write(*,*)ptree%MyID,ranks(ii*2-1+bb-1-Bidxs+1),rank,'aha'
			ranks(ii*2-1+bb-1-Bidxs+1) = rank
			stats%Flop_Fill = stats%Flop_Fill + flop
		endif
	enddo
	
	call MPI_ALLREDUCE(MPI_IN_PLACE ,ranks(ii*2-1-Bidxs+1),1,MPI_INTEGER,MPI_MIN,ptree%pgrp(block_inv%pgno)%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE ,ranks(ii*2-Bidxs+1),1,MPI_INTEGER,MPI_MIN,ptree%pgrp(block_inv%pgno)%Comm,ierr)
	! write(*,*)'wonima',ii*2-1-Bidxs+1,ranks(ii*2-1-Bidxs+1),ii*2-Bidxs+1,ranks(ii*2-Bidxs+1),ptree%MyID,mm(1)
	
	
	
	!!!*** redistribute AR from process layout of block_off1 and block_off2  to the process layout of hodlr
	do bb=1,2
		block_off => ho_bf1%levels(level)%BP(ii*2-1+bb-1)%LL(1)%matrices_block(1)
		M_p => block_off%M_p
		if(bb==1)matrixtemp=>matrixtemp1
		if(bb==2)matrixtemp=>matrixtemp2
			
		n1 = OMP_get_wtime()
		call Redistribute1Dto1D(matrixtemp,M_p,offout(bb),block_off%pgno,AR,block_inv%M_p,0,block_inv%pgno,ranks(ii*2-1+bb-1-Bidxs+1),ptree)	
		n2 = OMP_get_wtime()
		stats%Time_RedistV = stats%Time_RedistV + n2-n1	
		if(mm(bb)>0)then
			deallocate(matrixtemp)
		endif	
	enddo	
 
end subroutine PComputeRange_twoforward



!!!!!***** this subroutine is part of the randomized SVD. 
! Given B^T = (Q^cA)^T (N_loc x ranks(bb)) and Q (M_loc x ranks(bb)) in the process layout of hodlr, it computes SVD B=USV and output A = (QU)*(SV)
subroutine PSVDTruncate_twoforward(ho_bf1,level,Bidxs,bb_inv,ranks,Q,QcA_trans,block_rand,eps,ptree,stats)
   use HODLR_DEFS
   implicit none
   integer ranks(:)
   integer level, ii, bb,bb_inv
   DT :: Q(:,:),QcA_trans(:,:)
   DT,pointer :: mat1(:,:),mat2(:,:),mat(:,:),matQ1(:,:),matQ2(:,:),matQ(:,:),matQ2D(:,:),matQcA_trans1(:,:),matQcA_trans2(:,:),matQcA_trans(:,:),matQcA_trans2D(:,:),matQUt2D(:,:),UU(:,:),VV(:,:),mattemp(:,:)
   type(matrixblock),pointer::block_o,block_inv,block_schur,block_off1,block_off2,block_off
   integer groupn,groupm,mm(2),nn(2),ierr,nin1,nout1,nin2,nout2,offM(2),offN(2),offout1,offout2,rank
   type(hobf)::ho_bf1
   type(proctree)::ptree
   type(Hstat)::stats
   real(kind=8)::n1,n2,eps,flop
   integer,pointer::M_p(:,:),N_p(:,:)
   real(kind=8),pointer::Singular(:)
   integer descQ2D(9), descQcA_trans2D(9), descUU(9), descVV(9),descQUt2D(9)
   integer tempi,ctxt,info,iproc,jproc,myi,myj,myArows,myAcols,myrow,mycol,nprow,npcol,M,N,mnmin
   type(matrixblock)::block_rand(:)
   integer Bidxs
   stats%Flop_Tmp=0
	block_off1 => ho_bf1%levels(level)%BP(bb_inv*2-1)%LL(1)%matrices_block(1)	
	block_off2 => ho_bf1%levels(level)%BP(bb_inv*2)%LL(1)%matrices_block(1)
	block_inv => ho_bf1%levels(level)%BP_inverse(bb_inv)%LL(1)%matrices_block(1)	

	mm(1)=block_off1%M_loc  
	mm(2)=block_off2%M_loc		
	nn(1)=block_off1%N_loc  
	nn(2)=block_off2%N_loc	
	
	offM(1) = 0
	offM(2) = block_off1%M
	offN(1) = block_off1%M
	offN(2) = 0
	
	!!!*** redistribute Q and B^T from process layout of hodlr to the process layout of block_off1 and block_off2
	n1 = OMP_get_wtime()
	do bb=1,2
		block_off => ho_bf1%levels(level)%BP(bb_inv*2-1+bb-1)%LL(1)%matrices_block(1)
		M_p => block_off%M_p
		N_p => block_off%N_p
		if(mm(bb)>0)then
			if(bb==1)then
			allocate(matQ1(mm(bb),ranks(bb_inv*2-1+bb-1-Bidxs+1)))
			matQ=>matQ1
			endif
			if(bb==2)then
			allocate(matQ2(mm(bb),ranks(bb_inv*2-1+bb-1-Bidxs+1)))
			matQ=>matQ2
			endif		
		endif	
		call Redistribute1Dto1D(Q,block_inv%M_p,0,block_inv%pgno,matQ,M_p,offM(bb),block_off%pgno,ranks(bb_inv*2-1+bb-1-Bidxs+1),ptree)
		
		
		if(nn(bb)>0)then
			if(bb==1)then
			allocate(matQcA_trans1(nn(bb),ranks(bb_inv*2-1+bb-1-Bidxs+1)))
			matQcA_trans=>matQcA_trans1
			endif
			if(bb==2)then
			allocate(matQcA_trans2(nn(bb),ranks(bb_inv*2-1+bb-1-Bidxs+1)))
			matQcA_trans=>matQcA_trans2
			endif		
		endif	
		call Redistribute1Dto1D(QcA_trans,block_inv%N_p,0,block_inv%pgno,matQcA_trans,N_p,offN(bb),block_off%pgno,ranks(bb_inv*2-1+bb-1-Bidxs+1),ptree)	
	enddo
	n2 = OMP_get_wtime()
	stats%Time_RedistV = stats%Time_RedistV + n2-n1	
		
	!!!*** compute B^T=V^TS^TU^T and A = (QU)*(SV) 
	do bb=1,2
		block_off => ho_bf1%levels(level)%BP(bb_inv*2-1+bb-1)%LL(1)%matrices_block(1)
		M_p => block_off%M_p
		N_p => block_off%N_p
		M = block_off%M
		N = block_off%N
		if(bb==1)matQcA_trans=>matQcA_trans1
		if(bb==2)matQcA_trans=>matQcA_trans2
		if(bb==1)matQ=>matQ1
		if(bb==2)matQ=>matQ2
		if(mm(bb)>0)then
		
			!!!!**** generate 2D grid blacs quantities	
			ctxt = ptree%pgrp(block_off%pgno)%ctxt		
			call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)	
			if(myrow/=-1 .and. mycol/=-1)then
				myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(ranks(bb_inv*2-1+bb-1-Bidxs+1), nbslpk, mycol, 0, npcol)
				! write(*,*)ptree%MyID,'descQ2D',M, ranks(bb_inv*2-1+bb-1-Bidxs+1)
				call descinit( descQ2D, M, ranks(bb_inv*2-1+bb-1-Bidxs+1), nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descQ2D')
				allocate(matQ2D(myArows,myAcols))
				matQ2D=0			
				
				myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(ranks(bb_inv*2-1+bb-1-Bidxs+1), nbslpk, mycol, 0, npcol)
				! write(*,*)ptree%MyID,'descQcA_trans2D',N, ranks(bb_inv*2-1+bb-1-Bidxs+1)
				call descinit( descQcA_trans2D, N, ranks(bb_inv*2-1+bb-1-Bidxs+1), nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
				call assert(info==0,'descinit fail for descQcA_trans2D')
				allocate(MatQcA_trans2D(myArows,myAcols))
				MatQcA_trans2D=0				
				
				mnmin=min(N,ranks(bb_inv*2-1+bb-1-Bidxs+1))

				myArows = numroc_wp(N, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(mnmin, nbslpk, mycol, 0, npcol)		
				allocate(UU(myArows,myAcols))
				! write(*,*)ptree%MyID,'descUU',N, mnmin
				call descinit( descUU, N, mnmin, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descUU')
				UU=0
				
				myArows = numroc_wp(mnmin, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(ranks(bb_inv*2-1+bb-1-Bidxs+1), nbslpk, mycol, 0, npcol)		
				allocate(VV(myArows,myAcols))
				! write(*,*)ptree%MyID,'descVV', mnmin, ranks(bb_inv*2-1+bb-1-Bidxs+1)
				call descinit( descVV, mnmin, ranks(bb_inv*2-1+bb-1-Bidxs+1), nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descVV')
				VV=0
				
				allocate(Singular(mnmin))
				Singular=0

			else 
				descQ2D(2)=-1
				descQcA_trans2D(2)=-1
				descUU(2)=-1
				descVV(2)=-1
				allocate(matQ2D(1,1))   ! required for Redistribute1Dto2D
				matQ2D=0
				allocate(matQcA_trans2D(1,1)) ! required for Redistribute1Dto2D
				matQcA_trans2D=0
				allocate(UU(1,1))  ! required for Redistribute2Dto1D
				UU=0
				allocate(VV(1,1))  
				VV=0
			endif
			
			!!!!**** redistribution into 2D grid
			call Redistribute1Dto2D(matQ,M_p,0,block_off%pgno,matQ2D,M,0,block_off%pgno,ranks(bb_inv*2-1+bb-1-Bidxs+1),ptree)	
			call Redistribute1Dto2D(matQcA_trans,N_p,0,block_off%pgno,matQcA_trans2D,N,0,block_off%pgno,ranks(bb_inv*2-1+bb-1-Bidxs+1),ptree)	
			
			! write(*,*)block_off%row_group,block_off%col_group,fnorm(matQcA_trans,N,ranks(bb_inv*2-1+bb-1-Bidxs+1)),'matQcA_trans'
			! write(*,*)block_off%row_group,block_off%col_group,fnorm(matQ,M,ranks(bb_inv*2-1+bb-1-Bidxs+1)),'matQ'
			
			
			!!!!**** compute B^T=V^TS^TU^T
			rank=0
			if(myrow/=-1 .and. mycol/=-1)then
				call PSVD_Truncate(N, ranks(bb_inv*2-1+bb-1-Bidxs+1),matQcA_trans2D,descQcA_trans2D,UU,VV,descUU,descVV,Singular,eps,rank,ctxt,flop=flop)
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol) 
				do ii=1,rank
					call g2l(ii,rank,npcol,nbslpk,jproc,myj)
					if(jproc==mycol)then
						UU(:,myj) = UU(:,myj)*Singular(ii) 		
					endif
				enddo


				myArows = numroc_wp(M, nbslpk, myrow, 0, nprow)
				myAcols = numroc_wp(rank, nbslpk, mycol, 0, npcol)		
				allocate(matQUt2D(myArows,myAcols))
				! write(*,*)'descQUt2D', M, rank
				call descinit( descQUt2D, M, rank, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )	
				call assert(info==0,'descinit fail for descQUt2D')
				matQUt2D=0				
				
				call pgemmf90('N','T',M,rank,ranks(bb_inv*2-1+bb-1-Bidxs+1),cone, matQ2D,1,1,descQ2D,VV,1,1,descVV,czero,matQUt2D,1,1,descQUt2D,flop=flop)
				stats%Flop_Fill = stats%Flop_Fill + flop/dble(nprow*npcol)				
			else 
				allocate(matQUt2D(1,1)) ! required for Redistribute2Dto1D
			endif
			
			!!!!**** copy results to block_rand
			call MPI_ALLREDUCE(MPI_IN_PLACE ,rank,1,MPI_INTEGER,MPI_MAX,ptree%pgrp(block_off%pgno)%Comm,ierr)
			call delete_blocks(block_rand(bb_inv*2-1+bb-1-Bidxs+1),0)
			! block_rand(bb_inv*2-1+bb-1-Bidxs+1)%style = 2
			! block_rand(bb_inv*2-1+bb-1-Bidxs+1)%level_butterfly = 0
			block_rand(bb_inv*2-1+bb-1-Bidxs+1)%rankmax = rank
			block_rand(bb_inv*2-1+bb-1-Bidxs+1)%rankmin = rank

			! block_rand(bb_inv*2-1+bb-1-Bidxs+1)%row_group=block_off%row_group
			! block_rand(bb_inv*2-1+bb-1-Bidxs+1)%col_group=block_off%col_group			
			
			allocate(block_rand(bb_inv*2-1+bb-1-Bidxs+1)%ButterflyU%blocks(1))
			allocate(block_rand(bb_inv*2-1+bb-1-Bidxs+1)%ButterflyV%blocks(1))		
			allocate(block_rand(bb_inv*2-1+bb-1-Bidxs+1)%ButterflyU%blocks(1)%matrix(block_off%M_loc,rank))
			allocate(block_rand(bb_inv*2-1+bb-1-Bidxs+1)%ButterflyV%blocks(1)%matrix(block_off%N_loc,rank))
						
			
			!!!!**** redistribution into 1D grid conformal to leaf sizes
			call Redistribute2Dto1D(matQUt2D,M,0,block_off%pgno,block_rand(bb_inv*2-1+bb-1-Bidxs+1)%ButterflyU%blocks(1)%matrix,block_off%M_p,0,block_off%pgno,rank,ptree)	
			call Redistribute2Dto1D(UU,N,0,block_off%pgno,block_rand(bb_inv*2-1+bb-1-Bidxs+1)%ButterflyV%blocks(1)%matrix,block_off%N_p,0,block_off%pgno,rank,ptree)	

			! allocate(mattemp(block_off%M_loc,block_off%N_loc))
			! mattemp=0
			! call gemmf90(block_rand(bb_inv*2-1+bb-1-Bidxs+1)%ButterflyU%blocks(1)%matrix,block_off%M_loc,block_rand(bb_inv*2-1+bb-1-Bidxs+1)%ButterflyV%blocks(1)%matrix,block_off%N_loc,mattemp,block_off%M_loc,'N','T',block_off%M_loc,block_off%N_loc,rank,cone,czero)
			
			! write(*,*)block_off%row_group,block_off%col_group,fnorm(block_rand(bb_inv*2-1+bb-1-Bidxs+1)%ButterflyU%blocks(1)%matrix,block_off%M_loc,rank),fnorm(block_rand(bb_inv*2-1+bb-1-Bidxs+1)%ButterflyV%blocks(1)%matrix,block_off%N_loc,rank),'aha'
			! ! write(*,*)block_off%row_group,block_off%col_group,fnorm(mattemp,block_off%M_loc,rank),'aha'
			! deallocate(mattemp)
			
			if(myrow/=-1 .and. mycol/=-1)then
				! deallocate(matQ2D)
				! deallocate(MatQcA_trans2D)
				! deallocate(UU)
				! deallocate(VV)
				deallocate(Singular)
				
			endif
				deallocate(matQ2D)
				deallocate(MatQcA_trans2D)
				deallocate(UU)
				deallocate(VV)
				deallocate(matQUt2D)
				! deallocate(Singular)
				! deallocate(matQUt2D)			
			deallocate(matQcA_trans)
			deallocate(matQ)
		endif
	enddo
 
end subroutine PSVDTruncate_twoforward

end module HODLR_randomMVP
