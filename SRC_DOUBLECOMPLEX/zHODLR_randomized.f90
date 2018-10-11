#include "HODLR_config.fi"
module z_HODLR_randomMVP
use z_Bplus_randomized 
use z_HODLR_Solve_Mul
use z_Bplus_compress


contains

subroutine HODLR_randomized(ho_bf1,blackbox_HODLR_MVP_randomized_dat,N,rankmax,Memory,error,option,stats,operand,ptree,msh)
	
	
    use z_HODLR_DEFS
    implicit none
	real(kind=8):: n1,n2,Memory,error_inout,Memtmp,tmpfact,error
	integer level_c,rankmax,level_butterfly,bb,rank_new_max,ii,groupm,groupn,N
	type(matrixblock),pointer::block_o,block_ref
	DT,allocatable::Vin(:,:),Vout1(:,:),Vout2(:,:)
	type(matrixblock),allocatable::block_rand(:)
	type(hobf)::ho_bf1
	type(Hoption)::option
	type(Hstat)::stats
	real(kind=8):: time_gemm1
	class(*)::operand
	procedure(HOBF_MVP_blk)::blackbox_HODLR_MVP_randomized_dat
	type(proctree)::ptree
	type(mesh)::msh
	
	if(.not. allocated(stats%rankmax_of_level))allocate (stats%rankmax_of_level(ho_bf1%Maxlevel))
	stats%rankmax_of_level = 0

	Memory = 0
	do level_c = 1,ho_bf1%Maxlevel+1
		if(level_c==ho_bf1%Maxlevel+1)then
			call HODLR_randomized_OneL_Fullmat(ho_bf1,blackbox_HODLR_MVP_randomized_dat,N,level_c,Memtmp,operand,ptree,stats,msh)
		else
			if(level_c>option%LRlevel)then
				level_butterfly = 0   !!! uncomment out this to enable low-rank only reconstruction
			else 
				level_butterfly=int((ho_bf1%Maxlevel-level_c)/2)*2
			endif

			! write(*,*)'before init'
			n1 = OMP_get_wtime()
			allocate (block_rand(ho_bf1%levels(level_c)%N_block_forward))
			do bb = 1, ho_bf1%levels(level_c)%N_block_forward
				groupm=ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)%row_group         ! Note: row_group and col_group interchanged here   
				groupn=ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)%col_group         ! Note: row_group and col_group interchanged here   
				call BF_Init_randomized(level_butterfly,rankmax,groupm,groupn,ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1),block_rand(bb))						
			enddo
			n2 = OMP_get_wtime()
			stats%Time_random(1) = stats%Time_random(1) + n2-n1
			
			! if(level_butterfly==0)then
				! call HODLR_MVP_randomized_OneL_Lowrank(ho_bf1,block_rand,blackbox_HODLR_MVP_randomized_dat,N,level_c,rankmax,option,operand,ptree,stats,msh)
			! else 

				n1 = OMP_get_wtime()
				call HODLR_Reconstruction_LL(ho_bf1,block_rand,blackbox_HODLR_MVP_randomized_dat,N,level_c,level_butterfly,option,stats,operand,ptree,msh)
				n2 = OMP_get_wtime()
				write(*,*)'reconstructLL: ', n2-n1, 'vecCNT',vecCNT

				n1 = OMP_get_wtime()
				call HODLR_Reconstruction_RR(ho_bf1,block_rand,blackbox_HODLR_MVP_randomized_dat,N,level_c,level_butterfly,option,stats,operand,ptree,msh)
				n2 = OMP_get_wtime()
				write(*,*)'reconstructRR: ', n2-n1, 'vecCNT',vecCNT
			! end if
			
			
			call HODLR_Test_Error_RR(ho_bf1,block_rand,blackbox_HODLR_MVP_randomized_dat,N,level_c,error_inout,operand,ptree,stats,msh)
			
			rank_new_max = 0
			do bb = 1, ho_bf1%levels(level_c)%N_block_forward
				block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				call get_butterfly_minmaxrank(block_rand(bb))
				rank_new_max = max(rank_new_max,block_rand(bb)%rankmax)
				call copy_butterfly('N',block_rand(bb),block_o,Memtmp)
				Memory = Memory + Memtmp
				call delete_blocks(block_rand(bb))
			end do
			write(*,'(A10,I5,A6,I3,A8,I3,A7,Es14.7)')' Level ',level_c,' rank:',rank_new_max,' L_butt:',level_butterfly,' error:',error_inout
			deallocate (block_rand)
			stats%rankmax_of_level(level_c) = rank_new_max

			
		end if
	end do

	write(*,*)  'rankmax_of_level:',stats%rankmax_of_level
	
	
	
	allocate(Vin(N,1))
	allocate(Vout1(N,1))
	allocate(Vout2(N,1))
	do ii=1,N
		call random_dp_number(Vin(ii,1))
	end do

	call blackbox_HODLR_MVP_randomized_dat('N',N,1,Vin,Vout1,msh,operand)
	call MVM_Z_forward('N',N,1,1,ho_bf1%Maxlevel+1,Vin,Vout2,ho_bf1,ptree,stats)
	
	error = fnorm(Vout2-Vout1,N,1)/fnorm(Vout1,N,1)
	deallocate(Vin,Vout1,Vout2)
	
end subroutine HODLR_randomized



subroutine HODLR_MVP_randomized_OneL_Lowrank(ho_bf1,block_rand,blackbox_HODLR_MVP_randomized_dat,N,level_c,rmax,option,operand,ptree,stats,msh)
	
	
    use z_HODLR_DEFS
    implicit none
	real(kind=8):: n1,n2,Memory,error_inout,Memtmp
	integer mn,rankref,level_c,rmax,rmaxloc,level_butterfly,bb,rank_new_max,rank,num_vect,groupn,groupm,header_n,header_m,tailer_m,tailer_n,ii,jj,k,mm,nn
	type(matrixblock),pointer::block_o,block_ref
	DT,allocatable::RandVectTmp(:,:)
	DT, allocatable :: matRcol(:,:),matZRcol(:,:),matRrow(:,:),matZcRrow(:,:),mattmp1(:,:),mattmp2(:,:),matrixtemp(:,:),matrixtempQ(:,:)
	real(kind=8), allocatable:: Singular(:)
	DT::ctemp1,ctemp2
	DT, allocatable::UU(:,:),VV(:,:)
	integer q,qq,N
	integer,allocatable::perms(:), ranks(:) 
	type(hobf)::ho_bf1
	type(matrixblock)::block_rand(:)
	type(Hoption)::option
	DT,allocatable:: RandVectInR(:,:),RandVectOutR(:,:),RandVectInL(:,:),RandVectOutL(:,:)
	class(*)::operand
	type(proctree)::ptree
	type(Hstat)::stats
	type(mesh)::msh
	
	procedure(HOBF_MVP_blk)::blackbox_HODLR_MVP_randomized_dat
	
	
	level_butterfly = 0
	
	
	ctemp1=1.0d0 ; ctemp2=0.0d0	
	
	rank_new_max = 0
	
	num_vect = rmax
	allocate(RandVectInR(N,num_vect))
	RandVectInR=0
	allocate(RandVectOutR(N,num_vect))
	
	allocate(ranks(ho_bf1%levels(level_c)%N_block_forward))
	
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=block_o%N
		header_n=block_o%headn
		tailer_n=nn+header_n-1
		k=header_n-1	
		
		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=block_o%M
		
		ranks(bb)=min(min(mm,nn),num_vect)
		
		allocate(matrixtemp(nn,num_vect))
		call RandomMat(nn,num_vect,min(nn,num_vect),matrixtemp,1)
		RandVectInR(1+k:nn+k,1:num_vect) = matrixtemp
		deallocate(matrixtemp)
	end do	
	
	call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP_randomized_dat,'N', RandVectInR, RandVectOutR, N,level_c,num_vect,operand,ptree,stats,msh)
	
	! power iteration of order q, the following is prone to roundoff error, see algorithm 4.4 Halko 2010
	q=0
	do qq=1,q
		RandVectOutR=conjg(cmplx(RandVectOutR,kind=8))
		call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP_randomized_dat,'T', RandVectOutR, RandVectInR, N,level_c,num_vect,operand,ptree,stats,msh)
		RandVectInR=conjg(cmplx(RandVectInR,kind=8))
		call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP_randomized_dat,'N', RandVectInR, RandVectOutR, N,level_c,num_vect,operand,ptree,stats,msh)
	enddo
	
	
	! computation of range Q
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=block_o%M
		header_m=block_o%headm
		tailer_m=mm+header_m-1
		
		k=header_m-1
		
		allocate(matrixtemp(mm,ranks(bb)))
		matrixtemp = RandVectOutR(header_m:header_m+mm-1,1:ranks(bb))
		
		call ComputeRange(mm,ranks(bb),matrixtemp,rank,0,option%tol_comp)			
		ranks(bb) = rank
		RandVectOutR(header_m:header_m+mm-1,1:ranks(bb)) = matrixtemp(1:mm,1:ranks(bb))
		deallocate(matrixtemp)
	end do

	! computation of B = Q^c*A
	RandVectOutR=conjg(cmplx(RandVectOutR,kind=8))
	call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP_randomized_dat,'T', RandVectOutR, RandVectInR, N,level_c,num_vect,operand,ptree,stats,msh)
	RandVectOutR=conjg(cmplx(RandVectOutR,kind=8))
	
	! computation of SVD of B and LR of A
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=block_o%M
		header_m=block_o%headm
		tailer_m=mm+header_m-1
		

		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=block_o%N
		header_n=block_o%headn
		tailer_n=nn+header_n-1		
		
		allocate(matrixtemp(ranks(bb),nn))
		call copymatT(RandVectInR(header_n:header_n+nn-1,1:ranks(bb)),matrixtemp,nn,ranks(bb))
		
				
		mn=min(ranks(bb),nn)
		allocate (UU(ranks(bb),mn),VV(mn,nn),Singular(mn))
		call SVD_Truncate(matrixtemp,ranks(bb),nn,mn,UU,VV,Singular,option%tol_comp,rank)				
		do ii=1,rank
			UU(:,ii) = UU(:,ii)* Singular(ii)
		end do	

		rank_new_max = max(rank_new_max,rank)					
		
		call delete_blocks(block_rand(bb))
		
		block_rand(bb)%style = 2
		block_rand(bb)%level_butterfly = 0
		block_rand(bb)%rankmax = rank
		block_rand(bb)%rankmin = rank

		block_rand(bb)%row_group=-1
		block_rand(bb)%col_group=-1		
		
		
		allocate(block_rand(bb)%ButterflyU%blocks(1))
		allocate(block_rand(bb)%ButterflyV%blocks(1))		
		
		allocate(block_rand(bb)%ButterflyV%blocks(1)%matrix(nn,rank))
		call copymatT(VV(1:rank,1:nn),block_rand(bb)%ButterflyV%blocks(1)%matrix,rank,nn)
		allocate(block_rand(bb)%ButterflyU%blocks(1)%matrix(mm,rank))
		allocate(matrixtempQ(mm,ranks(bb)))
		matrixtempQ = RandVectOutR(header_m:header_m+mm-1,1:ranks(bb))
		
		! call gemm_omp(matrixtempQ,UU(1:ranks(bb),1:rank),block_rand(bb)%ButterflyU%blocks(1)%matrix,mm,rank,ranks(bb))
		call gemmf90(matrixtempQ,mm,UU,ranks(bb),block_rand(bb)%ButterflyU%blocks(1)%matrix,mm,'N','N',mm,rank,ranks(bb),cone,czero)
		
		
		deallocate(matrixtemp,matrixtempQ,UU,VV,Singular)

	end do		

	deallocate(RandVectOutR,RandVectInR,ranks)


end subroutine HODLR_MVP_randomized_OneL_Lowrank


subroutine HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP_randomized_dat,trans, VectIn, VectOut, N,level_c,num_vect,operand,ptree,stats,msh)
	
	
    use z_HODLR_DEFS
    implicit none
	real(kind=8):: n1,n2,Memory,error_inout,Memtmp
	integer N,mn,rankref,level_c,rmax,rmaxloc,level_butterfly,bb,rank_new_max,rank,num_vect,groupn,groupm,header_n,header_m,tailer_m,tailer_n,ii,jj,k,mm,nn
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
	class(*)::operand
	type(proctree)::ptree
	type(Hstat)::stats
	type(mesh)::msh
	
	procedure(HOBF_MVP_blk)::blackbox_HODLR_MVP_randomized_dat

	
	level_butterfly = 0
	
		! write(*,*)'haha1'
		ctemp1=1.0d0 ; ctemp2=0.0d0	
		
		Memory = 0
		rank_new_max = 0
		
		! num_vect = rmax
		
		if(trans=='N')then

			VectOut = 0
			
			allocate(RandVectIn(N,num_vect))
			allocate(RandVectOut(N,num_vect))
			
			allocate(RandVectTmp(N,num_vect))
			
			! Compute the odd block MVP first
			RandVectIn=0
			do bb=1,ho_bf1%levels(level_c)%N_block_forward-1,2
				block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
				nn=block_o%N
				header_n=block_o%headn
				tailer_n=nn+header_n-1
				! nn=tailer_n-header_n+1
				k=header_n-1	
				RandVectIn(1+k:nn+k,1:num_vect) = VectIn(1+k:nn+k,1:num_vect)
			end do	
			
			call blackbox_HODLR_MVP_randomized_dat('N',N,num_vect,RandVectIn,RandVectTmp,msh,operand)
			call MVM_Z_forward('N',N,num_vect,1,level_c-1,RandVectIn,RandVectOut,ho_bf1,ptree,stats)
			RandVectOut = RandVectTmp-RandVectOut			
				
			do bb=1,ho_bf1%levels(level_c)%N_block_forward-1,2
				block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   

				mm=block_o%M
				header_m=block_o%headm
				tailer_m = mm + header_m-1
				VectOut(header_m:header_m+mm-1,1:num_vect) = RandVectOut(header_m:header_m+mm-1,1:num_vect)
			end do	

			! Compute the even block MVP next
			RandVectIn=0
			do bb=2,ho_bf1%levels(level_c)%N_block_forward,2
				block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
				nn=block_o%N
				header_n=block_o%headn
				tailer_n=nn+header_n-1
				k=header_n-1	
				RandVectIn(1+k:nn+k,1:num_vect) = VectIn(1+k:nn+k,1:num_vect)
			end do	
			
			call blackbox_HODLR_MVP_randomized_dat('N',N,num_vect,RandVectIn,RandVectTmp,msh,operand)
			
			call MVM_Z_forward('N',N,num_vect,1,level_c-1,RandVectIn,RandVectOut,ho_bf1,ptree,stats)
			RandVectOut = RandVectTmp-RandVectOut			
				
			do bb=2,ho_bf1%levels(level_c)%N_block_forward,2
				block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
				mm=block_o%M
				header_m=block_o%headm
				tailer_m = mm + header_m-1
				VectOut(header_m:header_m+mm-1,1:num_vect) = RandVectOut(header_m:header_m+mm-1,1:num_vect)
			end do
				
			deallocate(RandVectIn)	
			deallocate(RandVectOut)	
			deallocate(RandVectTmp)	
			
		else if(trans=='T')then
			VectOut = 0
			
			allocate(RandVectIn(N,num_vect))
			allocate(RandVectOut(N,num_vect))
			
			allocate(RandVectTmp(N,num_vect))
			
			! Compute the odd block MVP first
			RandVectIn=0
			do bb=1,ho_bf1%levels(level_c)%N_block_forward-1,2
				block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
				mm=block_o%M
				header_m=block_o%headm
				tailer_m = mm + header_m-1
				k=header_m-1
				RandVectIn(1+k:mm+k,1:num_vect) = VectIn(1+k:mm+k,1:num_vect)
			end do	
		
			call blackbox_HODLR_MVP_randomized_dat('T',N,num_vect,RandVectIn,RandVectTmp,msh,operand)
			
			call MVM_Z_forward('T',N,num_vect,1,level_c-1,RandVectIn,RandVectOut,ho_bf1,ptree,stats)
			RandVectOut = RandVectTmp-RandVectOut		
		
			do bb=1,ho_bf1%levels(level_c)%N_block_forward-1,2
				block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
				nn=block_o%N
				header_n=block_o%headn
				tailer_n=nn+header_n-1
				VectOut(header_n:header_n+nn-1,1:num_vect)=RandVectOut(header_n:header_n+nn-1,1:num_vect)
			end do			

			! Compute the even block MVP next
			RandVectIn=0
			do bb=2,ho_bf1%levels(level_c)%N_block_forward,2
				block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
				mm=block_o%M
				header_m=block_o%headm
				tailer_m = mm + header_m-1
				k=header_m-1
				RandVectIn(1+k:mm+k,1:num_vect) = VectIn(1+k:mm+k,1:num_vect)
			end do	
		
			call blackbox_HODLR_MVP_randomized_dat('T',N,num_vect,RandVectIn,RandVectTmp,msh,operand)
			call MVM_Z_forward('T',N,num_vect,1,level_c-1,RandVectIn,RandVectOut,ho_bf1,ptree,stats)
			RandVectOut = RandVectTmp-RandVectOut		
		
			do bb=2,ho_bf1%levels(level_c)%N_block_forward,2
				block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
				nn=block_o%N
				header_n=block_o%headn
				tailer_n=nn+header_n-1
				VectOut(header_n:header_n+nn-1,1:num_vect)=RandVectOut(header_n:header_n+nn-1,1:num_vect)
			end do
			
			deallocate(RandVectIn)	
			deallocate(RandVectOut)	
			deallocate(RandVectTmp)	
			
		endif
		
end subroutine HODLR_MVP_randomized_OneL



subroutine HODLR_randomized_OneL_Fullmat(ho_bf1,blackbox_HODLR_MVP_randomized_dat,N,level_c,Memory,operand,ptree,stats,msh)
	
	
    use z_HODLR_DEFS
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
	class(*)::operand
	type(proctree)::ptree
	type(Hstat)::stats
	type(mesh)::msh
	
	procedure(HOBF_MVP_blk)::blackbox_HODLR_MVP_randomized_dat
	
	
	ctemp1=1.0d0 ; ctemp2=0.0d0	
	
	Memory = 0
	rank_new_max = 0
	
	
	num_vect = 0
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		nn=block_o%N
		header_n=block_o%headn
		tailer_n=nn+header_n-1
		num_vect = max(num_vect,nn)
	end do	
	
	allocate(RandVectInR(N,num_vect))
	RandVectInR=0	
	allocate(RandVectTmp(N,num_vect))
	allocate(RandVectOutR(N,num_vect))
	
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=block_o%N
		header_n=block_o%headn
		tailer_n=nn+header_n-1
		k=header_n-1	
		do ii=1, nn
			RandVectInR(ii+k,ii)=1d0
		enddo
	end do		
	
	call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP_randomized_dat,'N', RandVectInR, RandVectOutR, N,level_c,num_vect,operand,ptree,stats,msh)
	
	
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=block_o%M
		header_m=block_o%headm
		tailer_m = mm + header_m-1
		k=header_m-1
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		
		nn=block_o%N
		header_n=block_o%headn
		tailer_n=nn+header_n-1
		
		block_o%style = 1
		allocate(block_o%fullmat(mm,nn))
		! call copymatN(RandVectOutR(k+1:k+mm,1:nn),block_o%fullmat,mm,nn)
		block_o%fullmat = RandVectOutR(k+1:k+mm,1:nn)
		Memory = Memory + SIZEOF(block_o%fullmat)/1024.0d3
	end do		
	
	deallocate(RandVectInR,RandVectOutR,RandVectTmp)
	
	
	write(*,'(A10,I5,A13)')' Level ',level_c,' fullmat done'
	
	
end subroutine HODLR_randomized_OneL_Fullmat
		
	


subroutine HODLR_Reconstruction_LL(ho_bf1,block_rand,blackbox_HODLR_MVP_randomized_dat,N,level_c,level_butterfly,option,stats,operand,ptree,msh)
    
    use z_HODLR_DEFS
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
	class(*)::operand
	type(proctree)::ptree
	type(mesh)::msh
	
	procedure(HOBF_MVP_blk)::blackbox_HODLR_MVP_randomized_dat
	
	
	
	
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
			call HODLR_Randomized_Vectors_LL(ho_bf1,block_rand,vec_rand,blackbox_HODLR_MVP_randomized_dat,N,level_c,level_butterfly,nth_s,nth_e,num_vect_sub,unique_nth,operand,ptree,stats,msh)
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







subroutine HODLR_Reconstruction_RR(ho_bf1,block_rand,blackbox_HODLR_MVP_randomized_dat,N,level_c,level_butterfly,option,stats,operand,ptree,msh)
    
    use z_HODLR_DEFS
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
	class(*)::operand
	type(proctree)::ptree
	type(mesh)::msh
	
	procedure(HOBF_MVP_blk)::blackbox_HODLR_MVP_randomized_dat
	
	
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

			call HODLR_Randomized_Vectors_RR(ho_bf1,block_rand,vec_rand,blackbox_HODLR_MVP_randomized_dat,N,level_c,level_butterfly,nth_s,nth_e,num_vect_sub,unique_nth,operand,ptree,stats,msh)

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









subroutine HODLR_Test_Error_RR(ho_bf1,block_rand,blackbox_HODLR_MVP_randomized_dat,N,level_c,error,operand,ptree,stats,msh)

    use z_HODLR_DEFS
    implicit none
    
	integer nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,groupm,groupn,bb
    integer N,mm,nn
    real(kind=8) a,b,c,d, condition_number,norm1_R,norm2_R,norm3_R,norm4_R
    DT ctemp, ctemp1, ctemp2
    
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
	integer Nsub,Ng,num_vect,nth_s,nth_e,level_butterfly
	integer*8 idx_start
	real(kind=8)::error
	integer level_c,rowblock,dimension_m,header_m,tailer_m,header_n,tailer_n
	DT,allocatable::RandomVectors_Output_ref(:,:)
	type(matrixblock),pointer::block_o
	
	type(hobf)::ho_bf1
	type(matrixblock)::block_rand(:)
	type(vectorsblock),pointer:: RandomVectors_InOutput(:)
	class(*)::operand
	type(proctree)::ptree
	procedure(HOBF_MVP_blk)::blackbox_HODLR_MVP_randomized_dat
	type(Hstat)::stats
	type(mesh)::msh
	
	
	num_vect=1
	allocate (RandomVectors_InOutput(3))
	do ii=1,3
		allocate (RandomVectors_InOutput(ii)%vector(N,num_vect))
		RandomVectors_InOutput(ii)%vector = 0d0
	end do
	

	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=block_o%N
		header_n=block_o%headn
		tailer_n=nn+header_n-1

		k=header_n-1	
		do jj=1,num_vect
			do ii=1, nn
				call random_dp_number(RandomVectors_InOutput(1)%vector(ii+k,jj))	! matrixtemp1(jj,ii) ! 
			enddo
		enddo
	end do	
	

	call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP_randomized_dat,'N', RandomVectors_InOutput(1)%vector, RandomVectors_InOutput(3)%vector, N,level_c,num_vect,operand,ptree,stats,msh)
		
		
	ctemp1=1.0d0 ; ctemp2=0.0d0	
	RandomVectors_InOutput(2)%vector = 0d0	
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		! write(*,*)bb,ho_bf1%levels(level_c)%N_block_forward-1,'ha?'
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=block_o%M
		header_m=block_o%headm
		tailer_m = mm + header_m-1
		
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=block_o%N
		header_n=block_o%headn
		tailer_n=nn+header_n-1
	
		call BF_block_MVP_dat(block_rand(bb),'N',mm,nn,num_vect,RandomVectors_InOutput(1)%vector(header_n:header_n+nn-1,:),RandomVectors_InOutput(2)%vector(header_m:header_m+mm-1,:),ctemp1,ctemp2,ptree,stats)	
		
	end do		

	 
	error = fnorm(RandomVectors_InOutput(2)%vector-RandomVectors_InOutput(3)%vector,mm,num_vect)/fnorm(RandomVectors_InOutput(3)%vector,mm,num_vect)
	
	do ii=1,3
		deallocate (RandomVectors_InOutput(ii)%vector)
	end do		
	deallocate(RandomVectors_InOutput)
		

    return                

end subroutine HODLR_Test_Error_RR


subroutine HODLR_MVP_randomized_Fullmat(trans,N,num_vect,Vin,Vout,msh,operand)
	implicit none 
	character trans
	DT Vin(:,:),Vout(:,:)
	DT,allocatable:: Vin_tmp(:,:),Vout_tmp(:,:)
	DT ctemp,a,b
	integer ii,jj,num_vect,nn,fl_transpose,kk,black_step,N
	real(kind=8) n1,n2,tmp(2)
	type(mesh)::msh
	class(*)::operand
	
	select TYPE(operand)   
    type is (kernelquant)
	
		n1 = OMP_get_wtime()
		
		allocate(Vin_tmp(N,num_vect))
		allocate(Vout_tmp(N,num_vect))
		do ii=1,N
		Vin_tmp(msh%new2old(ii),:)=Vin(ii,:)
		enddo	
		a=1d0
		b=0d0
		call gemmf90(operand%matZ_glo,N,Vin_tmp,N,Vout_tmp,N,trans,'N',N,num_vect,N,a,b)	

		do ii=1,N
		Vout(ii,:)=Vout_tmp(msh%new2old(ii),:)
		enddo
		
		deallocate(Vout_tmp)
		deallocate(Vin_tmp)

		n2 = OMP_get_wtime()
		! time_gemm1 = time_gemm1 + n2 - n1
	end select
	
end subroutine HODLR_MVP_randomized_Fullmat



subroutine HODLR_Randomized_Vectors_LL(ho_bf1,block_rand,vec_rand,blackbox_HODLR_MVP_randomized_dat,N,level_c,level_butterfly,nth_s,nth_e,num_vect_sub,unique_nth,operand,ptree,stats,msh)

    use z_HODLR_DEFS
    
	use z_misc
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
	class(*)::operand
	type(proctree)::ptree
	type(Hstat)::stats
	type(mesh)::msh
	
	procedure(HOBF_MVP_blk)::blackbox_HODLR_MVP_randomized_dat
	
	
	
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
		! mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	
		groupm_start=groupm*2**(level_butterfly)
		
		do nth= nth_s,nth_e
			do i=1, num_blocks
				if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
					header_m=basis_group(groupm_start+i-1)%head
					tailer_m=basis_group(groupm_start+i-1)%tail
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
	
	
	call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP_randomized_dat,'T', RandomVectors_InOutput(1)%vector, RandomVectors_InOutput(3)%vector, N,level_c,num_vect_sub,operand,ptree,stats,msh)
	
	
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		! mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		groupm_start=groupm*2**(level_butterfly)
		! random=>random_Block(bb)
		do i=1, num_blocks
			header_m=basis_group(groupm_start+i-1)%head
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
		! nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		groupn_start=groupn*2**(level_butterfly)
		do i=1, num_blocks
			header_n=basis_group(groupn_start+i-1)%head
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





subroutine HODLR_Randomized_Vectors_RR(ho_bf1,block_rand,vec_rand,blackbox_HODLR_MVP_randomized_dat,N,level_c,level_butterfly,nth_s,nth_e,num_vect_sub,unique_nth,operand,ptree,stats,msh)

    use z_HODLR_DEFS
    
	use z_misc
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
	class(*)::operand
	procedure(HOBF_MVP_blk)::blackbox_HODLR_MVP_randomized_dat
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
		! nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
		groupn_start=groupn*2**(level_butterfly)
		
		do nth= nth_s,nth_e
			do i=1, num_blocks
				if(i>=(nth-1)*Ng+1 .and. i<=nth*Ng)then	
					header_n=basis_group(groupn_start+i-1)%head
					tailer_n=basis_group(groupn_start+i-1)%tail
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
	
	

	call HODLR_MVP_randomized_OneL(ho_bf1,blackbox_HODLR_MVP_randomized_dat,'N', RandomVectors_InOutput(1)%vector, RandomVectors_InOutput(3)%vector, N,level_c,num_vect_sub,operand,ptree,stats,msh)
	
	
	do bb=1,ho_bf1%levels(level_c)%N_block_forward
		block_o =>  ho_bf1%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		! nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		groupn_start=groupn*2**(level_butterfly)
		! random=>random_Block(bb)
		do i=1, num_blocks
			header_n=basis_group(groupn_start+i-1)%head
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
		! mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		groupm_start=groupm*2**(level_butterfly)
		do i=1, num_blocks
			header_m=basis_group(groupm_start+i-1)%head
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


end module z_HODLR_randomMVP
