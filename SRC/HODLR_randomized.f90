module HODLR_randomMVP
use Utilites_randomized 
! use H_structure
use Butterfly_inversion
contains



subroutine HODLR_MVP(rankmax,Memory,error)
	! use lapack95
	! use blas95
    use MODULE_FILE
    implicit none
	real*8:: n1,n2,Memory,error_inout,Memtmp,tmpfact,error
	integer level_c,rankmax,level_butterfly,bb,rank_new_max,ii
	type(matrixblock),pointer::block_o,block_ref
	complex(kind=8),allocatable::Vin(:,:),Vout1(:,:),Vout2(:,:)
	
	if(.not. allocated(rankmax_of_level))allocate (rankmax_of_level(Maxlevel_for_blocks))
	rankmax_of_level = 0
	
	tmpfact = Rank_detection_factor
	Rank_detection_factor = Rank_detection_factor !* 0.001
	
	Memory = 0
	do level_c = 1,Maxlevel_for_blocks+1
		if(level_c==Maxlevel_for_blocks+1)then
			call HODLR_MVP_OneL_Fullmat(level_c,Memtmp)
		else 
			level_butterfly=int((maxlevel_for_blocks-level_c)/2)*2
		
			level_butterfly = 0   !!! uncomment out this to enable low-rank only reconstruction
		
			! if(level_butterfly==0)then
				! call HODLR_MVP_OneL_Lowrank(level_c,rankmax,Memtmp)
				! ! write(*,*)'sorry man'
				! ! stop
				! Memory = Memory + Memtmp
			! else 
				time_gemm1 = 0
				! write(*,*)'before init'
				call Initialize_Butterfly_HODLR_MVP(rankmax,level_c,level_butterfly)
				! write(*,*)'before reconstructLL'
				n1 = OMP_get_wtime()
				call Reconstruction_LL_HODLR_MVP(level_c,level_butterfly)
				n2 = OMP_get_wtime()
				write(*,*)'reconstructLL: ', n2-n1, 'time_gemm1',time_gemm1,'vecCNT',vecCNT
				time_gemm1 = 0
				! write(*,*)'before reconstructRR'
				n1 = OMP_get_wtime()
				call Reconstruction_RR_HODLR_MVP(level_c,level_butterfly,error_inout)
				n2 = OMP_get_wtime()
				write(*,*)'reconstructRR: ', n2-n1, 'time_gemm1',time_gemm1,'vecCNT',vecCNT
				
				rank_new_max = 0
				do bb = 1, ho_bf%levels(level_c)%N_block_forward
					block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
					call get_randomizedbutterfly_minmaxrank(butterfly_block_randomized(bb))
					rank_new_max = max(rank_new_max,butterfly_block_randomized(bb)%rankmax)
					call copy_randomizedbutterfly(butterfly_block_randomized(bb),block_o,Memtmp)
					Memory = Memory + Memtmp
					call Delete_randomized_butterflys(bb)
				end do
				write(*,'(A10,I5,A6,I3,A8,I3,A7,Es14.7)')' Level ',level_c,' rank:',rank_new_max,' L_butt:',level_butterfly,' error:',error_inout
				deallocate (butterfly_block_randomized)
				rankmax_of_level(level_c) = rank_new_max
			! end if
		end if
	end do

	Rank_detection_factor = tmpfact
	write(*,*)  'rankmax_of_level:',rankmax_of_level
	
	
	
	allocate(Vin(Maxedge,1))
	allocate(Vout1(Maxedge,1))
	allocate(Vout2(Maxedge,1))
	do ii=1,Maxedge
		Vin(ii,1) = random_complex_number()
	end do
	call GetOutputs_BlackBox('N',Vin,Vout1,1)
	call MVM_Z_forward_partial('N',Maxedge,1,Maxlevel_for_blocks+1,Vin,Vout2,ho_bf)	
	error = fnorm(Vout2-Vout1,Maxedge,1)/fnorm(Vout1,Maxedge,1)
	deallocate(Vin,Vout1,Vout2)
	
end subroutine HODLR_MVP



subroutine HODLR_MVP_OneL_Lowrank(level_c,rmax,Memory)
	! use lapack95
	! use blas95
    use MODULE_FILE
    implicit none
	real*8:: n1,n2,Memory,error_inout,Memtmp
	integer mn,rankref,level_c,rmax,rmaxloc,level_butterfly,bb,rank_new_max,rank,num_vect,groupn,groupm,header_n,header_m,tailer_m,tailer_n,ii,jj,k,mm,nn
	type(matrixblock),pointer::block_o,block_ref
	complex(kind=8),allocatable::RandVectTmp(:,:)
	complex(kind=8), allocatable :: matRcol(:,:),matZRcol(:,:),matRrow(:,:),matZcRrow(:,:),mattmp1(:,:),mattmp2(:,:),matrixtemp(:,:),matrixtempQ(:,:)
	real*8, allocatable:: Singular(:)
	complex(kind=8)::ctemp1,ctemp2
	complex (kind=8), allocatable::UU(:,:),VV(:,:)
	integer q,qq
	integer,allocatable::perms(:), ranks(:) 
	
	level_butterfly = 0
	
	if(level_c>0)then ! This MVP-based randomized scheme is stable only when the random vectors are "good enough"
		
		ctemp1=1.0d0 ; ctemp2=0.0d0	
		
		Memory = 0
		rank_new_max = 0
		
		num_vect = rmax
		allocate(RandVectInR(Maxedge,num_vect))
		RandVectInR=0
		allocate(RandVectOutR(Maxedge,num_vect))
		
		allocate(ranks(ho_bf%levels(level_c)%N_block_forward))
		
		do bb=1,ho_bf%levels(level_c)%N_block_forward
			block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
			groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
			nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
			header_n=basis_group(groupn)%head
			tailer_n=basis_group(groupn)%tail
			nn=tailer_n-header_n+1
			k=header_n-1	
			
			groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
			mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
			
			ranks(bb)=min(min(mm,nn),num_vect)
			
			allocate(matrixtemp(nn,num_vect))
			call RandomMat(nn,num_vect,min(nn,num_vect),matrixtemp,1)
			RandVectInR(1+k:nn+k,1:num_vect) = matrixtemp
			deallocate(matrixtemp)
		end do	
		
		call HODLR_MVP_OneL('N', RandVectInR, RandVectOutR, level_c,num_vect)
		
		! power iteration of order q, the following is prone to roundoff error, see algorithm 4.4 Halko 2010
		q=0
		do qq=1,q
			RandVectOutR=conjg(RandVectOutR)
			call HODLR_MVP_OneL('T', RandVectOutR, RandVectInR, level_c,num_vect)
			RandVectInR=conjg(RandVectInR)
			call HODLR_MVP_OneL('N', RandVectInR, RandVectOutR, level_c,num_vect)
		enddo
		
		
		! computation of range Q
		do bb=1,ho_bf%levels(level_c)%N_block_forward
			block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
			groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
			header_m=basis_group(groupm)%head
			tailer_m=basis_group(groupm)%tail
			mm=tailer_m-header_m+1 
			k=header_m-1
			
			allocate(matrixtemp(mm,ranks(bb)))
			matrixtemp = RandVectOutR(header_m:header_m+mm-1,1:ranks(bb))
			
			call ComputeRange(mm,ranks(bb),matrixtemp,rank,0,SVD_tolerance_factor)			
			ranks(bb) = rank
			RandVectOutR(header_m:header_m+mm-1,1:ranks(bb)) = matrixtemp(1:mm,1:ranks(bb))
			deallocate(matrixtemp)
		end do

		! computation of B = Q^c*A
		RandVectOutR=conjg(RandVectOutR)
		call HODLR_MVP_OneL('T', RandVectOutR, RandVectInR, level_c,num_vect)
		RandVectOutR=conjg(RandVectOutR)
		
		! computation of SVD of B and LR of A
		do bb=1,ho_bf%levels(level_c)%N_block_forward
			block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
			groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
			header_m=basis_group(groupm)%head
			tailer_m=basis_group(groupm)%tail
			mm=tailer_m-header_m+1 
	
			groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
			nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
			header_n=basis_group(groupn)%head
			tailer_n=basis_group(groupn)%tail			
			
			allocate(matrixtemp(ranks(bb),nn))
			call copymatT_omp(RandVectInR(header_n:header_n+nn-1,1:ranks(bb)),matrixtemp,nn,ranks(bb))
			
					
			mn=min(ranks(bb),nn)
			allocate (UU(ranks(bb),mn),VV(mn,nn),Singular(mn))
			call SVD_Truncate(matrixtemp,ranks(bb),nn,mn,UU,VV,Singular,SVD_tolerance_factor,rank)				
			do ii=1,rank
				UU(:,ii) = UU(:,ii)* Singular(ii)
			end do	

			rank_new_max = max(rank_new_max,rank)					
			block_o%style = 2
			block_o%level_butterfly = 0
			block_o%rankmax = rank
			block_o%rankmin = rank
			allocate(block_o%ButterflyU(1))
			allocate(block_o%ButterflyV(1))		
			
			allocate(block_o%ButterflyV(1)%matrix(nn,rank))
			call copymatT_omp(VV(1:rank,1:nn),block_o%ButterflyV(1)%matrix,rank,nn)
			allocate(block_o%ButterflyU(1)%matrix(mm,rank))
			allocate(matrixtempQ(mm,ranks(bb)))
			matrixtempQ = RandVectOutR(header_m:header_m+mm-1,1:ranks(bb))
			
			call gemm_omp(matrixtempQ,UU(1:ranks(bb),1:rank),block_o%ButterflyU(1)%matrix,mm,ranks(bb),rank)
			
			deallocate(matrixtemp,matrixtempQ,UU,VV,Singular)
			Memory = Memory + SIZEOF(block_o%ButterflyU(1)%matrix)/1024.0d3
			Memory = Memory + SIZEOF(block_o%ButterflyV(1)%matrix)/1024.0d3
		end do		
	
		deallocate(RandVectOutR,RandVectInR,ranks)

	else  ! This MVP-based scheme contructs the full matrix first and then compress. Very stable. 
	
	
		ctemp1=1.0d0 ; ctemp2=0.0d0	
		
		Memory = 0
		rank_new_max = 0
		
		
		num_vect = 0
		do bb=1,ho_bf%levels(level_c)%N_block_forward
			block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
			groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
			nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
			num_vect = max(num_vect,nn)
		end do	
		
		allocate(RandVectInR(Maxedge,num_vect))
		RandVectInR=0	
		allocate(RandVectTmp(Maxedge,num_vect))
		allocate(RandVectOutR(Maxedge,num_vect))
		
		do bb=1,ho_bf%levels(level_c)%N_block_forward
			block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
			groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
			nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
			header_n=basis_group(groupn)%head
			tailer_n=basis_group(groupn)%tail
			k=header_n-1	
			do ii=1, nn
				RandVectInR(ii+k,ii)=1d0
			enddo
		end do		
		
		call HODLR_MVP_OneL('N', RandVectInR, RandVectOutR, level_c,num_vect)
		
		do bb=1,ho_bf%levels(level_c)%N_block_forward
			block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
			groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
			mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
			header_m=basis_group(groupm)%head
			tailer_m=basis_group(groupm)%tail
			k=header_m-1
			groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
			nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
			
			
			allocate(matrixtemp(mm,nn))
			call copymatN_omp(RandVectOutR(k+1:k+mm,1:nn),matrixtemp,mm,nn)			
			mn=min(mm,nn)
			allocate (UU(mm,mn),VV(mn,nn),Singular(mn))
			call SVD_Truncate(matrixtemp,mm,nn,mn,UU,VV,Singular,SVD_tolerance_factor,rank)				
			do ii=1,rank
				UU(:,ii) = UU(:,ii)* Singular(ii)
			end do			
						
			rank_new_max = max(rank_new_max,rank)					
			block_o%style = 2
			block_o%level_butterfly = 0
			block_o%rankmax = rank
			block_o%rankmin = rank
			allocate(block_o%ButterflyU(1))
			allocate(block_o%ButterflyV(1))		
			
			allocate(block_o%ButterflyV(1)%matrix(nn,rank))
			call copymatT_omp(VV(1:rank,1:nn),block_o%ButterflyV(1)%matrix,rank,nn)
			allocate(block_o%ButterflyU(1)%matrix(mm,rank))
			call copymatN_omp(UU(1:mm,1:rank),block_o%ButterflyU(1)%matrix,mm,rank)
			deallocate(matrixtemp,UU,VV,Singular)
			Memory = Memory + SIZEOF(block_o%ButterflyU(1)%matrix)/1024.0d3
			Memory = Memory + SIZEOF(block_o%ButterflyV(1)%matrix)/1024.0d3

		end do		
		
		deallocate(RandVectInR,RandVectOutR,RandVectTmp)	
	
	end if
	
! write(*,*)'haha6'
	call Test_Error_RR_HODLR_MVP_Lowrank(level_c,error_inout)

	write(*,'(A10,I5,A6,I3,A8,I3,A7,Es14.7)')' Level ',level_c,' rank:',rank_new_max,' L_butt:',level_butterfly,' error:',error_inout
	rankmax_of_level(level_c) = rank_new_max
	
	
end subroutine HODLR_MVP_OneL_Lowrank


subroutine HODLR_MVP_OneL(trans, VectIn, VectOut, level_c,num_vect)
	! use lapack95
	! use blas95
    use MODULE_FILE
    implicit none
	real*8:: n1,n2,Memory,error_inout,Memtmp
	integer mn,rankref,level_c,rmax,rmaxloc,level_butterfly,bb,rank_new_max,rank,num_vect,groupn,groupm,header_n,header_m,tailer_m,tailer_n,ii,jj,k,mm,nn
	type(matrixblock),pointer::block_o,block_ref
	complex(kind=8),allocatable::RandVectTmp(:,:)
	complex(kind=8), allocatable :: matRcol(:,:),matZRcol(:,:),matRrow(:,:),matZcRrow(:,:),mattmp1(:,:),mattmp2(:,:),matrixtemp(:,:)
	real*8, allocatable:: Singular(:)
	complex(kind=8)::ctemp1,ctemp2
	complex (kind=8), allocatable::UU(:,:),VV(:,:)
	integer,allocatable::perms(:) 
	character trans
	complex(kind=8)::VectIn(:,:),VectOut(:,:)
	complex(kind=8),allocatable:: RandVectIn(:,:),RandVectOut(:,:)
	level_butterfly = 0
	
		! write(*,*)'haha1'
		ctemp1=1.0d0 ; ctemp2=0.0d0	
		
		Memory = 0
		rank_new_max = 0
		
		! num_vect = rmax
		
		if(trans=='N')then

			VectOut = 0
			
			allocate(RandVectIn(Maxedge,num_vect))
			allocate(RandVectOut(Maxedge,num_vect))
			
			allocate(RandVectTmp(Maxedge,num_vect))
			
			! Compute the odd block MVP first
			RandVectIn=0
			do bb=1,ho_bf%levels(level_c)%N_block_forward-1,2
				block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
				nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
				header_n=basis_group(groupn)%head
				tailer_n=basis_group(groupn)%tail
				nn=tailer_n-header_n+1
				k=header_n-1	
				RandVectIn(1+k:nn+k,1:num_vect) = VectIn(1+k:nn+k,1:num_vect)
			end do	
			
			call GetOutputs_BlackBox('N',RandVectIn,RandVectTmp,num_vect)
			call MVM_Z_forward_partial('N',Maxedge,num_vect,level_c-1,RandVectIn,RandVectOut,ho_bf)
			RandVectOut = RandVectTmp-RandVectOut			
				
			do bb=1,ho_bf%levels(level_c)%N_block_forward-1,2
				block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
				mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
				header_m=basis_group(groupm)%head
				tailer_m=basis_group(groupm)%tail
				VectOut(header_m:header_m+mm-1,1:num_vect) = RandVectOut(header_m:header_m+mm-1,1:num_vect)
			end do	

			! Compute the even block MVP next
			RandVectIn=0
			do bb=2,ho_bf%levels(level_c)%N_block_forward,2
				block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
				nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
				header_n=basis_group(groupn)%head
				tailer_n=basis_group(groupn)%tail
				nn=tailer_n-header_n+1
				k=header_n-1	
				RandVectIn(1+k:nn+k,1:num_vect) = VectIn(1+k:nn+k,1:num_vect)
			end do	
			
			call GetOutputs_BlackBox('N',RandVectIn,RandVectTmp,num_vect)
			call MVM_Z_forward_partial('N',Maxedge,num_vect,level_c-1,RandVectIn,RandVectOut,ho_bf)
			RandVectOut = RandVectTmp-RandVectOut			
				
			do bb=2,ho_bf%levels(level_c)%N_block_forward,2
				block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
				mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
				header_m=basis_group(groupm)%head
				tailer_m=basis_group(groupm)%tail
				VectOut(header_m:header_m+mm-1,1:num_vect) = RandVectOut(header_m:header_m+mm-1,1:num_vect)
			end do
				
			deallocate(RandVectIn)	
			deallocate(RandVectOut)	
			deallocate(RandVectTmp)	
			
		else if(trans=='T')then
			VectOut = 0
			
			allocate(RandVectIn(Maxedge,num_vect))
			allocate(RandVectOut(Maxedge,num_vect))
			
			allocate(RandVectTmp(Maxedge,num_vect))
			
			! Compute the odd block MVP first
			RandVectIn=0
			do bb=1,ho_bf%levels(level_c)%N_block_forward-1,2
				block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
				mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
				header_m=basis_group(groupm)%head
				tailer_m=basis_group(groupm)%tail
				mm=tailer_m-header_m+1
				k=header_m-1
				RandVectIn(1+k:mm+k,1:num_vect) = VectIn(1+k:mm+k,1:num_vect)
			end do	
		
			call GetOutputs_BlackBox('T',RandVectIn,RandVectTmp,num_vect)
			call MVM_Z_forward_partial('T',Maxedge,num_vect,level_c-1,RandVectIn,RandVectOut,ho_bf)
			RandVectOut = RandVectTmp-RandVectOut		
		
			do bb=1,ho_bf%levels(level_c)%N_block_forward-1,2
				block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
				nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
				header_n=basis_group(groupn)%head
				tailer_n=basis_group(groupn)%tail	
				VectOut(header_n:header_n+nn-1,1:num_vect)=RandVectOut(header_n:header_n+nn-1,1:num_vect)
			end do			

			! Compute the even block MVP next
			RandVectIn=0
			do bb=2,ho_bf%levels(level_c)%N_block_forward,2
				block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
				mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
				header_m=basis_group(groupm)%head
				tailer_m=basis_group(groupm)%tail
				mm=tailer_m-header_m+1
				k=header_m-1
				RandVectIn(1+k:mm+k,1:num_vect) = VectIn(1+k:mm+k,1:num_vect)
			end do	
		
			call GetOutputs_BlackBox('T',RandVectIn,RandVectTmp,num_vect)
			call MVM_Z_forward_partial('T',Maxedge,num_vect,level_c-1,RandVectIn,RandVectOut,ho_bf)
			RandVectOut = RandVectTmp-RandVectOut		
		
			do bb=2,ho_bf%levels(level_c)%N_block_forward,2
				block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
				groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
				nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
				header_n=basis_group(groupn)%head
				tailer_n=basis_group(groupn)%tail	
				VectOut(header_n:header_n+nn-1,1:num_vect)=RandVectOut(header_n:header_n+nn-1,1:num_vect)
			end do
			
			deallocate(RandVectIn)	
			deallocate(RandVectOut)	
			deallocate(RandVectTmp)	
			
		endif
		
end subroutine HODLR_MVP_OneL















subroutine HODLR_MVP_OneL_Fullmat(level_c,Memory)
	! use lapack95
	! use blas95
    use MODULE_FILE
    implicit none
	real*8:: n1,n2,Memory,error_inout,Memtmp
	integer rankref,level_c,rmaxloc,level_butterfly,bb,rank_new_max,rank,num_vect,groupn,groupm,header_n,header_m,tailer_m,tailer_n,ii,jj,k,mm,nn
	type(matrixblock),pointer::block_o,block_ref
	complex(kind=8),allocatable::RandVectTmp(:,:)
	complex(kind=8), allocatable :: matRcol(:,:),matZRcol(:,:),matRrow(:,:),matZcRrow(:,:),mattmp1(:,:),mattmp2(:,:)
	real*8, allocatable:: Singular(:)
	complex(kind=8)::ctemp1,ctemp2
	
	ctemp1=1.0d0 ; ctemp2=0.0d0	
	
	Memory = 0
	rank_new_max = 0
	
	
	num_vect = 0
	do bb=1,ho_bf%levels(level_c)%N_block_forward
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		num_vect = max(num_vect,nn)
	end do	
	
	allocate(RandVectInR(Maxedge,num_vect))
	RandVectInR=0	
	allocate(RandVectTmp(Maxedge,num_vect))
	allocate(RandVectOutR(Maxedge,num_vect))
	
	do bb=1,ho_bf%levels(level_c)%N_block_forward
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		header_n=basis_group(groupn)%head
		tailer_n=basis_group(groupn)%tail
		k=header_n-1	
		do ii=1, nn
			RandVectInR(ii+k,ii)=1d0
		enddo
	end do		
	
	call GetOutputs_BlackBox('N',RandVectInR,RandVectTmp,num_vect)
	call MVM_Z_forward_partial('N',Maxedge,num_vect,level_c-1,RandVectInR,RandVectOutR,ho_bf)
	RandVectOutR = RandVectTmp-RandVectOutR	
	
	do bb=1,ho_bf%levels(level_c)%N_block_forward
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		header_m=basis_group(groupm)%head
		tailer_m=basis_group(groupm)%tail
		k=header_m-1
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		
		block_o%style = 1
		allocate(block_o%fullmat(mm,nn))
		call copymatN_omp(RandVectOutR(k+1:k+mm,1:nn),block_o%fullmat,mm,nn)
		Memory = Memory + SIZEOF(block_o%fullmat)/1024.0d3
	end do		
	
	deallocate(RandVectInR,RandVectOutR,RandVectTmp)
	
	
	write(*,'(A10,I5,A13)')' Level ',level_c,' fullmat done'
	
	
end subroutine HODLR_MVP_OneL_Fullmat
		
		



subroutine Initialize_Butterfly_HODLR_MVP(rankmax,level_c,level_butterfly)
	use misc
    use MODULE_FILE
	! use lapack95
	! use blas95
    implicit none
    
    integer level_c,rowblock,kover
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn, bb
    integer dimension_rank, dimension_m, dimension_n, blocks, groupm, groupn,tmpi,tmpj,dimension_max
    real*8 a,b,c,d
    complex (kind=8) ctemp
	complex (kind=8), allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real*8, allocatable:: Singular(:)
	integer mn_min,index_i,index_j,rankmax,groupm_start,groupn_start
	! type(matrixblock),pointer::block
	
	
    allocate (butterfly_block_randomized(ho_bf%levels(level_c)%N_block_forward))
    ! level_butterfly=int((maxlevel_for_blocks-level_c)/2)*2
    dimension_rank = rankmax
	num_blocks=2**level_butterfly
	
	
	
	
	do bb = 1, ho_bf%levels(level_c)%N_block_forward
		butterfly_block_randomized(bb)%level_butterfly=level_butterfly
		
		groupm=ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)%row_group         ! Note: row_group and col_group interchanged here   
		groupn=ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)%col_group         ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1
		
		butterfly_block_randomized(bb)%dimension_rank=dimension_rank
		!write (*,*) dimension_rank, int(real(kk)/num_blocks/2.0), mm
		
		allocate (butterfly_block_randomized(bb)%ButterflyU(2**level_butterfly))
		allocate (butterfly_block_randomized(bb)%ButterflyV(2**level_butterfly))
		
		groupm_start=groupm*2**(level_butterfly)
		groupn_start=groupn*2**(level_butterfly)

		
		dimension_max = 2*dimension_rank
		do blocks=1, num_blocks	
			dimension_m=basis_group(groupm_start+blocks-1)%tail-basis_group(groupm_start+blocks-1)%head+1
			dimension_n=basis_group(groupn_start+blocks-1)%tail-basis_group(groupn_start+blocks-1)%head+1
			dimension_max = max(dimension_max,dimension_m)	
			dimension_max = max(dimension_max,dimension_n)	
		end do	
		allocate(butterfly_block_randomized(bb)%KerInv(dimension_max,dimension_max))
		call RandomMat(dimension_max,dimension_max,dimension_max,butterfly_block_randomized(bb)%KerInv,3)
		
		do blocks=1, num_blocks
				
			dimension_m=basis_group(groupm_start+blocks-1)%tail-basis_group(groupm_start+blocks-1)%head+1
			allocate (butterfly_block_randomized(bb)%ButterflyU(blocks)%matrix(dimension_m,dimension_rank))

			
			allocate(matrixtemp1(dimension_rank,dimension_m))
			call RandomMat(dimension_rank,dimension_m,min(dimension_rank,dimension_m),matrixtemp1,0)
			do j=1, dimension_rank
				do i=1, dimension_m
					butterfly_block_randomized(bb)%ButterflyU(blocks)%matrix(i,j) = matrixtemp1(j,i)
				end do
			end do	
			deallocate(matrixtemp1)

			
			dimension_n=basis_group(groupn_start+blocks-1)%tail-basis_group(groupn_start+blocks-1)%head+1
			allocate (butterfly_block_randomized(bb)%ButterflyV(blocks)%matrix(dimension_n,dimension_rank))
			! allocate (butterfly_block_randomized(bb)%ButterflyVInv(blocks)%matrix(dimension_rank,dimension_n))		

			allocate(matrixtemp1(dimension_rank,dimension_n))
			call RandomMat(dimension_rank,dimension_n,min(dimension_rank,dimension_n),matrixtemp1,0)
			do j=1, dimension_rank
				do i=1, dimension_n
					butterfly_block_randomized(bb)%ButterflyV(blocks)%matrix(i,j) = matrixtemp1(j,i)
				end do
			end do	
			deallocate(matrixtemp1)

		enddo

		
		if (level_butterfly/=0) then
			allocate (matrixtemp1(2*dimension_rank,2*dimension_rank))
			allocate (butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly))

			do level=1, level_butterfly
				num_row=2**level
				num_col=2**(level_butterfly-level+1)
				butterfly_block_randomized(bb)%ButterflyKerl(level)%num_row=num_row
				butterfly_block_randomized(bb)%ButterflyKerl(level)%num_col=num_col
				allocate (butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(num_row,num_col))

				do j=1, num_col, 2
					index_j=int((j+1)/2)
					do i=1, num_row, 2
						index_i=int((i+1)/2)
						allocate (butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix(dimension_rank,dimension_rank))
						allocate (butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i+1,j)%matrix(dimension_rank,dimension_rank))
						allocate (butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j+1)%matrix(dimension_rank,dimension_rank))
						allocate (butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i+1,j+1)%matrix(dimension_rank,dimension_rank))

						call RandomMat(2*dimension_rank,2*dimension_rank,2*dimension_rank,matrixtemp1,0)	
						
						! !$omp parallel do default(shared) private(ii,jj)
						do jj=1, dimension_rank
							do ii=1, dimension_rank
								butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix(ii,jj)=matrixtemp1(ii,jj)
								butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i+1,j)%matrix(ii,jj)=matrixtemp1(ii+dimension_rank,jj)
								butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j+1)%matrix(ii,jj)=matrixtemp1(ii,jj+dimension_rank)
								butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i+1,j+1)%matrix(ii,jj)=matrixtemp1(ii+dimension_rank,jj+dimension_rank)
							enddo
						enddo
						! !$omp end parallel do
						
					enddo
				enddo
			enddo
			deallocate (matrixtemp1)
		endif	
		
		
	end do
	
	
   


  
    return

end subroutine Initialize_Butterfly_HODLR_MVP





subroutine Reconstruction_LL_HODLR_MVP(level_c,level_butterfly)
    
    use MODULE_FILE
    implicit none
	
    integer level_c,rowblock
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj,bb, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub,num_vect_subsub
    complex(kind=8) ctemp, a, b
    character chara
	integer level_right_start,num_col,num_row
    
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
    real*8::n1,n2

    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
    type(matrixblock),pointer::block_o
    integer::rank_new_max,dimension_rank
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth  ! level# where each block is touched only once  
	real*8:: error_inout
	integer,allocatable::perms(:)
	
	vecCNT = 0
	
    allocate (Random_Block(ho_bf%levels(level_c)%N_block_forward))
    ! level_butterfly=int((maxlevel_for_blocks-level_c)/2)*2
    dimension_rank = butterfly_block_randomized(1)%dimension_rank   ! be careful here
	num_blocks=2**level_butterfly	
	num_vect_subsub= dimension_rank+5 ! be careful with the oversampling factor here
	Nbind = 1	
	num_vect_sub = num_vect_subsub*Nbind
	
	do bb=1,ho_bf%levels(level_c)%N_block_forward
		allocate (Random_Block(bb)%RandomVectorLL(0:level_butterfly+2)) 
		random=>Random_Block(bb)
		call Init_RandVects('T',random,num_vect_sub,bb)
		level_right_start = floor_safe(level_butterfly/2d0) !  check here later
		call Zero_Butterflys(0,level_right_start,bb)
	end do
		
	do unique_nth = 0,level_right_start
		Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later	
		Ng = 2**level_butterfly/Nsub
		
		do ii = 1,Nsub/Nbind	
			nth_s = (ii-1)*Nbind+1
			nth_e = ii*Nbind
			
			n1 = OMP_get_wtime()
			! First: The odd numbered blocks at level_c
			call Get_Randomized_Vectors_LL_HODLR_MVP(level_c,level_butterfly,nth_s,nth_e,num_vect_sub,unique_nth,1)
			! Second: The even numbered blocks at level_c
			call Get_Randomized_Vectors_LL_HODLR_MVP(level_c,level_butterfly,nth_s,nth_e,num_vect_sub,unique_nth,0)
			vecCNT = vecCNT + num_vect_sub*2
			
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1
			Time_Vector_forward = Time_Vector_forward + n2-n1
			
			n1 = OMP_get_wtime()
			do bb=1,ho_bf%levels(level_c)%N_block_forward
				call Resolving_Butterflys_LL_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth,bb)
			end do
			n2 = OMP_get_wtime()
			time_resolve = time_resolve + n2-n1		
		end do
	end do
	

	do bb=1,ho_bf%levels(level_c)%N_block_forward
		random=>Random_Block(bb)
		do level=0, level_butterfly+2
			num_row=random%RandomVectorLL(level)%num_row
			num_col=random%RandomVectorLL(level)%num_col
			do j=1, num_col
				do i=1, num_row
					deallocate (random%RandomVectorLL(level)%blocks(i,j)%matrix)
					if(level/=0 .and. level/=level_butterfly+2)deallocate (random%RandomVectorLL(level)%blocks(i,j)%null)
				enddo
			enddo
			deallocate (random%RandomVectorLL(level)%blocks)
		enddo
		deallocate (random%RandomVectorLL)
	end do

    return
    
end subroutine Reconstruction_LL_HODLR_MVP







subroutine Reconstruction_RR_HODLR_MVP(level_c,level_butterfly,error)
    
    use MODULE_FILE
    implicit none
	
    integer level_c,rowblock
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj,bb, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub,num_vect_subsub
    complex(kind=8) ctemp, a, b
    character chara
	integer level_right_start,num_col,num_row
    
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    integer Nsub,Ng,nth,nth_s,nth_e
	integer Nbind
    real*8::n1,n2

    integer blocks1, blocks2, blocks3, level_butterfly
    integer tt
    type(matrixblock),pointer::block_o
    integer::rank_new_max,dimension_rank,level_left_start
	real*8::rank_new_avr,error 
	complex(kind=8),allocatable::matrixtmp(:,:)
	integer niter,unique_nth  ! level# where each block is touched only once  
	real*8:: error_inout,rtemp
	integer,allocatable::perms(:)
	
	vecCNT = 0
	
    ! allocate (Random_Block(ho_bf%levels(level_c)%N_block_forward))
    ! level_butterfly=int((maxlevel_for_blocks-level_c)/2)*2
    dimension_rank = butterfly_block_randomized(1)%dimension_rank   ! be careful here
	num_blocks=2**level_butterfly	
	num_vect_subsub= dimension_rank+5 ! be careful with the oversampling factor here
	Nbind = 1	
	num_vect_sub = num_vect_subsub*Nbind
	
	do bb=1,ho_bf%levels(level_c)%N_block_forward
		allocate (Random_Block(bb)%RandomVectorRR(0:level_butterfly+2)) 
		random=>Random_Block(bb)
		call Init_RandVects('N',random,num_vect_sub,bb)
		level_left_start= floor_safe(level_butterfly/2d0)+1   !  check here later
		call Zero_Butterflys(level_left_start,level_butterfly+1,bb)
	end do
		
	do unique_nth=level_butterfly+1,level_left_start,-1
		if(mod(level_butterfly,2)==0)then
			Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))    !  check here later
		else 
			Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))
		end if			
		Ng = 2**level_butterfly/Nsub
		
		do ii = 1,Nsub/Nbind	
			nth_s = (ii-1)*Nbind+1
			nth_e = ii*Nbind
			
			n1 = OMP_get_wtime()
			! First: The odd numbered blocks at level_c
			call Get_Randomized_Vectors_RR_HODLR_MVP(level_c,level_butterfly,nth_s,nth_e,num_vect_sub,unique_nth,1)
			! Second: The even numbered blocks at level_c
			call Get_Randomized_Vectors_RR_HODLR_MVP(level_c,level_butterfly,nth_s,nth_e,num_vect_sub,unique_nth,0)
			vecCNT = vecCNT + num_vect_sub*2
			
			n2 = OMP_get_wtime()
			time_getvec = time_getvec + n2-n1
			Time_Vector_forward = Time_Vector_forward + n2-n1
			
			n1 = OMP_get_wtime()
			do bb=1,ho_bf%levels(level_c)%N_block_forward
				
				
				rtemp = 0
				random=>Random_Block(bb)
				do j=1,num_blocks
					nn=size(random%RandomVectorRR(0)%blocks(1,j)%matrix,1)
					rtemp = rtemp + fnorm(random%RandomVectorRR(0)%blocks(1,j)%matrix,nn,num_vect_subsub)
				end do
				! write(*,*)bb,ho_bf%levels(level_c)%N_block_forward,'woca',unique_nth,level_left_start,rtemp
				
				call Resolving_Butterflys_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth,bb)
			end do
			n2 = OMP_get_wtime()
			time_resolve = time_resolve + n2-n1		
		end do
	end do
	

	do bb=1,ho_bf%levels(level_c)%N_block_forward
		random=>Random_Block(bb)
		do level=0, level_butterfly+2
			num_row=random%RandomVectorRR(level)%num_row
			num_col=random%RandomVectorRR(level)%num_col
			do j=1, num_col
				do i=1, num_row
					deallocate (random%RandomVectorRR(level)%blocks(i,j)%matrix)
					if(level/=0 .and. level/=level_butterfly+2)deallocate (random%RandomVectorRR(level)%blocks(i,j)%null)
				enddo
			enddo
			deallocate (random%RandomVectorRR(level)%blocks)
		enddo
		deallocate (random%RandomVectorRR)
	end do
	deallocate(Random_Block)
! write(*,*)'cool'
	call Test_Error_RR_HODLR_MVP(level_c,error)
! write(*,*)'more cool'	
    return
    
end subroutine Reconstruction_RR_HODLR_MVP









subroutine Test_Error_RR_HODLR_MVP(level_c,error)

    use MODULE_FILE
    implicit none
    
	integer nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,groupm,groupn,bb
    integer mm,nn
    real*8 a,b,c,d, condition_number,norm1_R,norm2_R,norm3_R,norm4_R
    complex(kind=8) ctemp, ctemp1, ctemp2
    
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
	integer Nsub,Ng,num_vect,nth_s,nth_e,level_butterfly
	integer*8 idx_start
	real*8::error
	integer level_c,rowblock,dimension_m,header_m,tailer_m,header_n,tailer_n
	complex(kind=8),allocatable::RandomVectors_Output_ref(:,:)
	type(matrixblock),pointer::block_o
	
	num_vect=1
	allocate (RandomVectors_InOutput(3))
	do ii=1,3
		allocate (RandomVectors_InOutput(ii)%vector(Maxedge,num_vect))
		RandomVectors_InOutput(ii)%vector = 0d0
	end do
	
	! Test the odd block accuracy first
	do bb=1,ho_bf%levels(level_c)%N_block_forward-1,2
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		header_n=basis_group(groupn)%head
		tailer_n=basis_group(groupn)%tail
		nn=tailer_n-header_n+1
		k=header_n-1	
		do jj=1,num_vect
			do ii=1, nn
				RandomVectors_InOutput(1)%vector(ii+k,jj)=random_complex_number()	! matrixtemp1(jj,ii) ! 
			enddo
		enddo
	end do	
	
	call GetOutputs_BlackBox('N',RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,num_vect)
	call MVM_Z_forward_partial('N',Maxedge,num_vect,level_c-1,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector,ho_bf)
	RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(2)%vector-RandomVectors_InOutput(3)%vector
	
	do bb=2,ho_bf%levels(level_c)%N_block_forward,2
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		header_m=basis_group(groupm)%head
		tailer_m=basis_group(groupm)%tail
		mm=tailer_m-header_m+1
		k=header_m-1	
		do jj=1,num_vect
			do ii=1, mm
				RandomVectors_InOutput(3)%vector(ii+k,jj)=0	  
			enddo
		enddo
	end do		
		
		
		! do jj=1,num_vect
			! do ii=1, Maxedge
				! RandomVectors_InOutput(1)%vector(ii,jj)=random_complex_number()	! matrixtemp1(jj,ii) ! 
			! enddo
		! enddo	
		
	ctemp1=1.0d0 ; ctemp2=0.0d0	
	RandomVectors_InOutput(2)%vector = 0d0	
	do bb=1,ho_bf%levels(level_c)%N_block_forward-1,2
		! write(*,*)bb,ho_bf%levels(level_c)%N_block_forward-1,'ha?'
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		header_m=basis_group(groupm)%head
		
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		header_n=basis_group(groupn)%head
		
		call Butterflys_block_MVP_dat('N',mm,nn,num_vect,RandomVectors_InOutput(1)%vector(header_n:header_n+nn-1,:),RandomVectors_InOutput(2)%vector(header_m:header_m+mm-1,:),ctemp1,ctemp2,bb)		
	end do		
		! write(*,*)fnorm(RandomVectors_InOutput(1)%vector,Maxedge,num_vect),fnorm(RandomVectors_InOutput(2)%vector,Maxedge,num_vect),'gan'
	 norm1_R=0 ; norm2_R=0
	 do ii=1, mm
		do jj =1,num_vect
			 norm1_R=norm1_R+abs(RandomVectors_InOutput(3)%vector(ii,jj))**2
			 norm2_R=norm2_R+abs(RandomVectors_InOutput(2)%vector(ii,jj)-RandomVectors_InOutput(3)%vector(ii,jj))**2
		enddo
	 enddo
	 
	 ! Test the even block accuracy next
	do ii=1,3
		RandomVectors_InOutput(ii)%vector = 0d0
	end do	 
	do bb=2,ho_bf%levels(level_c)%N_block_forward,2
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		header_n=basis_group(groupn)%head
		tailer_n=basis_group(groupn)%tail
		nn=tailer_n-header_n+1
		k=header_n-1	
		do jj=1,num_vect
			do ii=1, nn
				RandomVectors_InOutput(1)%vector(ii+k,jj)=random_complex_number()	! matrixtemp1(jj,ii) ! 
			enddo
		enddo
	end do	
	
	call GetOutputs_BlackBox('N',RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,num_vect)
	call MVM_Z_forward_partial('N',Maxedge,num_vect,level_c-1,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector,ho_bf)
	RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(2)%vector-RandomVectors_InOutput(3)%vector
	
	do bb=1,ho_bf%levels(level_c)%N_block_forward-1,2
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		header_m=basis_group(groupm)%head
		tailer_m=basis_group(groupm)%tail
		mm=tailer_m-header_m+1
		k=header_m-1	
		do jj=1,num_vect
			do ii=1, mm
				RandomVectors_InOutput(3)%vector(ii+k,jj)=0	  
			enddo
		enddo
	end do		
		
		
	ctemp1=1.0d0 ; ctemp2=0.0d0	
	RandomVectors_InOutput(2)%vector = 0d0	
	do bb=2,ho_bf%levels(level_c)%N_block_forward,2
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		header_m=basis_group(groupm)%head
		
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		header_n=basis_group(groupn)%head
		
		call Butterflys_block_MVP_dat('N',mm,nn,num_vect,RandomVectors_InOutput(1)%vector(header_n:header_n+nn-1,:),RandomVectors_InOutput(2)%vector(header_m:header_m+mm-1,:),ctemp1,ctemp2,bb)		
	end do		
	

	 do ii=1, mm
		do jj =1,num_vect
			 norm1_R=norm1_R+abs(RandomVectors_InOutput(3)%vector(ii,jj))**2
			 norm2_R=norm2_R+abs(RandomVectors_InOutput(2)%vector(ii,jj)-RandomVectors_InOutput(3)%vector(ii,jj))**2
		enddo
	 enddo	 
	 
	error = sqrt(norm2_R/norm1_R)	
		! write(*,*)norm2_R,norm1_R,'dd'
	
	do ii=1,3
		deallocate (RandomVectors_InOutput(ii)%vector)
	end do		
	deallocate(RandomVectors_InOutput)
		

    return                

end subroutine Test_Error_RR_HODLR_MVP





subroutine Test_Error_RR_HODLR_MVP_Lowrank(level_c,error)

    use MODULE_FILE
    implicit none
    
	integer nth
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,groupm,groupn,bb
    integer mm,nn
    real*8 a,b,c,d, condition_number,norm1_R,norm2_R,norm3_R,norm4_R
    complex(kind=8) ctemp, ctemp1, ctemp2
    
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
	integer Nsub,Ng,num_vect,nth_s,nth_e,level_butterfly
	integer*8 idx_start
	real*8::error
	integer level_c,rowblock,dimension_m,header_m,tailer_m,header_n,tailer_n
	complex(kind=8),allocatable::RandomVectors_Output_ref(:,:)
	type(matrixblock),pointer::block_o
	
	num_vect=1
	allocate (RandomVectors_InOutput(3))
	do ii=1,3
		allocate (RandomVectors_InOutput(ii)%vector(Maxedge,num_vect))
		RandomVectors_InOutput(ii)%vector = 0d0
	end do
	
	! Test the odd block accuracy first
	do bb=1,ho_bf%levels(level_c)%N_block_forward-1,2
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		header_n=basis_group(groupn)%head
		tailer_n=basis_group(groupn)%tail
		nn=tailer_n-header_n+1
		k=header_n-1	
		do jj=1,num_vect
			do ii=1, nn
				RandomVectors_InOutput(1)%vector(ii+k,jj)=random_complex_number()	! matrixtemp1(jj,ii) ! 
			enddo
		enddo
	end do	
	
	call GetOutputs_BlackBox('N',RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,num_vect)
	call MVM_Z_forward_partial('N',Maxedge,num_vect,level_c-1,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector,ho_bf)
	RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(2)%vector-RandomVectors_InOutput(3)%vector
	
	do bb=2,ho_bf%levels(level_c)%N_block_forward,2
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		header_m=basis_group(groupm)%head
		tailer_m=basis_group(groupm)%tail
		mm=tailer_m-header_m+1
		k=header_m-1	
		do jj=1,num_vect
			do ii=1, mm
				RandomVectors_InOutput(3)%vector(ii+k,jj)=0	  
			enddo
		enddo
	end do		
		
		
		! do jj=1,num_vect
			! do ii=1, Maxedge
				! RandomVectors_InOutput(1)%vector(ii,jj)=random_complex_number()	! matrixtemp1(jj,ii) ! 
			! enddo
		! enddo	
		
	ctemp1=1.0d0 ; ctemp2=0.0d0	
	RandomVectors_InOutput(2)%vector = 0d0	
	do bb=1,ho_bf%levels(level_c)%N_block_forward-1,2
		! write(*,*)bb,ho_bf%levels(level_c)%N_block_forward-1,'ha?'
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		header_m=basis_group(groupm)%head
		
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		header_n=basis_group(groupn)%head
		
		call butterfly_block_MVP_randomized_dat(block_o,'N',mm,nn,num_vect,RandomVectors_InOutput(1)%vector(header_n:header_n+nn-1,:),RandomVectors_InOutput(2)%vector(header_m:header_m+mm-1,:),ctemp1,ctemp2)			
	end do	

		! write(*,*)fnorm(RandomVectors_InOutput(1)%vector,Maxedge,num_vect),fnorm(RandomVectors_InOutput(2)%vector,Maxedge,num_vect),'gan'
	 norm1_R=0 ; norm2_R=0
	 do ii=1, mm
		do jj =1,num_vect
			 norm1_R=norm1_R+abs(RandomVectors_InOutput(3)%vector(ii,jj))**2
			 norm2_R=norm2_R+abs(RandomVectors_InOutput(2)%vector(ii,jj)-RandomVectors_InOutput(3)%vector(ii,jj))**2
		enddo
	 enddo
	 
	 ! Test the even block accuracy next
	do ii=1,3
		RandomVectors_InOutput(ii)%vector = 0d0
	end do	 
	do bb=2,ho_bf%levels(level_c)%N_block_forward,2
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		header_n=basis_group(groupn)%head
		tailer_n=basis_group(groupn)%tail
		nn=tailer_n-header_n+1
		k=header_n-1	
		do jj=1,num_vect
			do ii=1, nn
				RandomVectors_InOutput(1)%vector(ii+k,jj)=random_complex_number()	! matrixtemp1(jj,ii) ! 
			enddo
		enddo
	end do	
	
	call GetOutputs_BlackBox('N',RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,num_vect)
	call MVM_Z_forward_partial('N',Maxedge,num_vect,level_c-1,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector,ho_bf)
	RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(2)%vector-RandomVectors_InOutput(3)%vector
	
	do bb=1,ho_bf%levels(level_c)%N_block_forward-1,2
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)
		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		header_m=basis_group(groupm)%head
		tailer_m=basis_group(groupm)%tail
		mm=tailer_m-header_m+1
		k=header_m-1	
		do jj=1,num_vect
			do ii=1, mm
				RandomVectors_InOutput(3)%vector(ii+k,jj)=0	  
			enddo
		enddo
	end do		
		
		
	ctemp1=1.0d0 ; ctemp2=0.0d0	
	RandomVectors_InOutput(2)%vector = 0d0	
	do bb=2,ho_bf%levels(level_c)%N_block_forward,2
		! write(*,*)bb,ho_bf%levels(level_c)%N_block_forward-1,'ha?'
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		header_m=basis_group(groupm)%head
		
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		header_n=basis_group(groupn)%head
		
		call butterfly_block_MVP_randomized_dat(block_o,'N',mm,nn,num_vect,RandomVectors_InOutput(1)%vector(header_n:header_n+nn-1,:),RandomVectors_InOutput(2)%vector(header_m:header_m+mm-1,:),ctemp1,ctemp2)			
	end do		
	

	 do ii=1, mm
		do jj =1,num_vect
			 norm1_R=norm1_R+abs(RandomVectors_InOutput(3)%vector(ii,jj))**2
			 norm2_R=norm2_R+abs(RandomVectors_InOutput(2)%vector(ii,jj)-RandomVectors_InOutput(3)%vector(ii,jj))**2
		enddo
	 enddo	 
	 
	error = sqrt(norm2_R/norm1_R)	
		! write(*,*)norm2_R,norm1_R,'dd'
	
	do ii=1,3
		deallocate (RandomVectors_InOutput(ii)%vector)
	end do		
	deallocate(RandomVectors_InOutput)
		

    return                

end subroutine Test_Error_RR_HODLR_MVP_Lowrank





subroutine GetOutputs_BlackBox(trans,Vin,Vout,num_vect)
	implicit none 
	character trans
	complex(kind=8) Vin(:,:),Vout(:,:)
	complex(kind=8),allocatable:: Vin_tmp(:,:),Vout_tmp(:,:)
	complex(kind=8) ctemp
	integer ii,jj,num_vect,nn,fl_transpose,kk,black_step
	real*8 n1,n2,tmp(2)
	
	n1 = OMP_get_wtime()
	
	if(fullmatflag==1)then
		if(trans=='N')then
			do nn=1,num_vect
				do ii=1,Maxedge
					ctemp = 0d0
					do jj=1,Maxedge
						ctemp = ctemp + matZ_glo(new2old(ii),new2old(jj))*Vin(jj,nn)
					end do
					Vout(ii,nn) = ctemp
				end do
			end do
		else if(trans=='T')then
			do nn=1,num_vect
				do jj=1,Maxedge
					ctemp = 0d0
					do ii=1,Maxedge
						ctemp = ctemp + matZ_glo(new2old(ii),new2old(jj))*Vin(ii,nn)
					end do
					Vout(jj,nn) = ctemp
				end do	
			end do
		end if
		
	else if(fullmatflag==-1)then
		write(*,*)'this option has been removed'
		stop
	else 
		write(*,*)'this option has been removed'
		stop
	end if
	n2 = OMP_get_wtime()
	time_gemm1 = time_gemm1 + n2 - n1

end subroutine GetOutputs_BlackBox




subroutine MVM_Z_forward_partial(trans,Ns,num_vectors,level_end,Vin,Vout,ho_bf1)

    use MODULE_FILE
    ! use lapack95
    implicit none
    
	character trans
	integer Ns, level_end
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	! type(matrixblock),pointer::block_o
	type(blockplus),pointer::bplus_o
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_m,idx_end_m,idx_start_n,idx_end_n
	
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	complex(kind=8)::Vin(:,:),Vout(:,:)
	! complex(kind=8)::Vin(Ns,num_vectors),Vout(Ns,num_vectors)
	type(hobf)::ho_bf1
 
	! num_vectors = 1   
	
	
	! get the right multiplied vectors
	ctemp1=1.0d0 ; ctemp2=1.0d0
	allocate(vec_old(Ns,num_vectors))
	allocate(vec_new(Ns,num_vectors))
	vec_old(1:Ns,1:num_vectors) = Vin
	vec_new = 0


	if(trans=='N')then
		do level = 1,level_end !Maxlevel_for_blocks+1
			do ii =1, ho_bf1%levels(level)%N_block_forward
				bplus_o =>  ho_bf1%levels(level)%BP(ii)
				groupm = bplus_o%row_group
				groupn = bplus_o%col_group
				idx_start_m = basis_group(groupm)%head
				idx_end_m = basis_group(groupm)%tail				
				idx_start_n = basis_group(groupn)%head
				idx_end_n = basis_group(groupn)%tail
				
				if(level==Maxlevel_for_blocks+1)then
					call fullmat_block_MVP_randomized_dat(bplus_o%LL(1)%matrices_block(1),'N',idx_end_m-idx_start_m+1,num_vectors,&
					&vec_old(idx_start_n:idx_end_n,1:num_vectors),vec_new(idx_start_m:idx_end_m,1:num_vectors),ctemp1,ctemp2)
				else 
					call Bplus_block_MVP_randomized_dat(bplus_o,'N',idx_end_m-idx_start_m+1,idx_end_n-idx_start_n+1,num_vectors,&
					&vec_old(idx_start_n:idx_end_n,1:num_vectors),vec_new(idx_start_m:idx_end_m,1:num_vectors),ctemp1,ctemp2)
				end if
			end do				
		end do
	else if(trans=='T')then
		do level = 1,level_end !Maxlevel_for_blocks+1
			do ii =1, ho_bf1%levels(level)%N_block_forward
				bplus_o =>  ho_bf1%levels(level)%BP(ii)
				groupm = bplus_o%row_group
				groupn = bplus_o%col_group
				idx_start_m = basis_group(groupm)%head
				idx_end_m = basis_group(groupm)%tail				
				idx_start_n = basis_group(groupn)%head
				idx_end_n = basis_group(groupn)%tail
				
				if(level==Maxlevel_for_blocks+1)then
					call fullmat_block_MVP_randomized_dat(bplus_o%LL(1)%matrices_block(1),'T',idx_end_m-idx_start_m+1,num_vectors,&
					&vec_old(idx_start_m:idx_end_m,1:num_vectors),vec_new(idx_start_n:idx_end_n,1:num_vectors),ctemp1,ctemp2)
				else 
					! write(*,*)bplus_o%LL(1)%matrices_block(1)%col_group,bplus_o%LL(1)%matrices_block(1)%row_group,bplus_o%Lplus
					call Bplus_block_MVP_randomized_dat(bplus_o,'T',idx_end_m-idx_start_m+1,idx_end_n-idx_start_n+1,num_vectors,&
					&vec_old(idx_start_m:idx_end_m,1:num_vectors),vec_new(idx_start_n:idx_end_n,1:num_vectors),ctemp1,ctemp2)
				end if
			end do				
		end do		
	end if

	Vout = vec_new(1:Ns,1:num_vectors)
	deallocate(vec_old)
	deallocate(vec_new)	
	 
    return                

end subroutine MVM_Z_forward_partial 



subroutine Zero_Butterflys(level_start,level_end,bb)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    ! type(matricesblock), pointer :: blocks
	integer level_start,level_end,level_butterfly,bb
	
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,rank,iii,jjj,index_i,index_j,mm_start,nn_start
    real*8 a,b,c,d
    complex(kind=8) ctemp
    
    complex(kind=8), allocatable :: matrixtemp(:,:), matrixU(:,:), matrixV(:,:)
    real*8, allocatable :: Singular(:)
    
    level_butterfly=butterfly_block_randomized(bb)%level_butterfly
    
    num_blocks=2**level_butterfly
	
	do level = level_start,level_end
		if(level==0)then
			do ii=1, num_blocks
				butterfly_block_randomized(bb)%ButterflyV(ii)%matrix = 0
			end do
		else if(level==level_butterfly+1)then
			do ii=1, num_blocks
				butterfly_block_randomized(bb)%ButterflyU(ii)%matrix = 0
			end do
		else 		
            num_row=2**level
            num_col=2**(level_butterfly-level+1)
			do j=1, num_col, 2
                do i=1, num_row, 2
                    butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix = 0
                    butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i+1,j)%matrix = 0
                    butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j+1)%matrix = 0
                    butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i+1,j+1)%matrix = 0
                enddo
            enddo		
		end if		
	end do
 
    return

end subroutine Zero_Butterflys





subroutine Init_RandVects(chara,random,num_vect_sub,bb)
    
    use MODULE_FILE
    implicit none
    
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, bb, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub,level_butterfly
    complex(kind=8) ctemp, a, b
    character chara
	integer num_col,num_row
    real*8:: mem_vec
	
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock), pointer :: random
    
    level_butterfly=butterfly_block_randomized(bb)%level_butterfly
	
    if (chara=='N') then

        num_blocks=2**level_butterfly
            
        allocate (random%RandomVectorRR(0)%blocks(1,num_blocks))
        random%RandomVectorRR(0)%num_row=1
        random%RandomVectorRR(0)%num_col=num_blocks

        do i=1, num_blocks
            nn=size(butterfly_block_randomized(bb)%butterflyV(i)%matrix,1)
            allocate (random%RandomVectorRR(0)%blocks(1,i)%matrix(nn,num_vect_sub))
			random%RandomVectorRR(0)%blocks(1,i)%matrix = 0
        enddo 

        do level=0, level_butterfly
            if (level==0) then
                num_groupn=num_blocks
                allocate (random%RandomVectorRR(1)%blocks(1,num_groupn))
                random%RandomVectorRR(1)%num_row=1
                random%RandomVectorRR(1)%num_col=num_groupn
                do j=1, num_groupn
                    rank=size(butterfly_block_randomized(bb)%ButterflyV(j)%matrix,2)
                    nn=size(butterfly_block_randomized(bb)%ButterflyV(j)%matrix,1)
                    allocate (random%RandomVectorRR(1)%blocks(1,j)%matrix(rank,num_vect_sub))
					random%RandomVectorRR(1)%blocks(1,j)%matrix = 0
					random%RandomVectorRR(1)%blocks(1,j)%nulldim = rank
					allocate (random%RandomVectorRR(1)%blocks(1,j)%null(rank,rank))
					random%RandomVectorRR(1)%blocks(1,j)%null = 0
					do ii=1,rank
						random%RandomVectorRR(1)%blocks(1,j)%null(ii,ii)=1
					end do
				enddo                    
            else
                num_groupm=butterfly_block_randomized(bb)%ButterflyKerl(level)%num_row
                num_groupn=butterfly_block_randomized(bb)%ButterflyKerl(level)%num_col               
				allocate (random%RandomVectorRR(level+1)%blocks(num_groupm,int(num_groupn/2)))
				random%RandomVectorRR(level+1)%num_row=num_groupm
				random%RandomVectorRR(level+1)%num_col=int(num_groupn/2)                    
				do i=1, num_groupm
					index_i=int((i+1)/2)
					do j=1, num_groupn, 2
						index_j=int((j+1)/2)
						nn1=size(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix,2)
						nn2=size(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j+1)%matrix,2)
						mm=size(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix,1)
						allocate (random%RandomVectorRR(level+1)%blocks(i,index_j)%matrix(mm,num_vect_sub))
						random%RandomVectorRR(level+1)%blocks(i,index_j)%matrix = 0
						random%RandomVectorRR(level+1)%blocks(i,index_j)%nulldim = mm
						allocate (random%RandomVectorRR(level+1)%blocks(i,index_j)%null(mm,mm))
						random%RandomVectorRR(level+1)%blocks(i,index_j)%null = 0
						do ii=1,mm
							random%RandomVectorRR(level+1)%blocks(i,index_j)%null(ii,ii)=1
						end do						
					enddo
				enddo               
            endif
            if (level==level_butterfly) then
                allocate (random%RandomVectorRR(level+2)%blocks(num_blocks,1))
                random%RandomVectorRR(level+2)%num_row=num_blocks
                random%RandomVectorRR(level+2)%num_col=1
                do i=1, num_blocks
                    rank=size(butterfly_block_randomized(bb)%ButterflyU(i)%matrix,2)
                    mm=size(butterfly_block_randomized(bb)%ButterflyU(i)%matrix,1)
                    allocate (random%RandomVectorRR(level+2)%blocks(i,1)%matrix(mm,num_vect_sub))
					random%RandomVectorRR(level+2)%blocks(i,1)%matrix = 0
				enddo
            endif
        enddo      

		mem_vec = 0
		do level=0, level_butterfly+2
			num_row=random%RandomVectorRR(level)%num_row
			num_col=random%RandomVectorRR(level)%num_col
			do j=1, num_col
				do i=1, num_row
					mem_vec =mem_vec + SIZEOF(random%RandomVectorRR(level)%blocks(i,j)%matrix)/1024.0d3
					if(level/=0 .and. level/=level_butterfly+2)mem_vec =mem_vec + SIZEOF(random%RandomVectorRR(level)%blocks(i,j)%null)/1024.0d3
				enddo
			enddo
		enddo
		Memory_int_vec = max(Memory_int_vec,mem_vec) 
        
    elseif (chara=='T') then
    
        num_blocks=2**level_butterfly
        allocate (random%RandomVectorLL(0)%blocks(num_blocks,1))
        random%RandomVectorLL(0)%num_row=num_blocks
        random%RandomVectorLL(0)%num_col=1
        
        do i=1, num_blocks
			mm=size(butterfly_block_randomized(bb)%butterflyU(i)%matrix,1)			
            allocate (random%RandomVectorLL(0)%blocks(i,1)%matrix(mm,num_vect_sub))
			random%RandomVectorLL(0)%blocks(i,1)%matrix = 0
        enddo 
        
        do level=0, level_butterfly
            if (level==0) then
                num_groupm=num_blocks
                allocate (random%RandomVectorLL(1)%blocks(num_groupm,1))
                random%RandomVectorLL(1)%num_row=num_groupm
                random%RandomVectorLL(1)%num_col=1
                do i=1, num_groupm
                    rank=size(butterfly_block_randomized(bb)%ButterflyU(i)%matrix,2)
                    mm=size(butterfly_block_randomized(bb)%ButterflyU(i)%matrix,1)
                    allocate (random%RandomVectorLL(1)%blocks(i,1)%matrix(rank,num_vect_sub))
					random%RandomVectorLL(1)%blocks(i,1)%matrix = 0
					random%RandomVectorLL(1)%blocks(i,1)%nulldim = rank
					allocate (random%RandomVectorLL(1)%blocks(i,1)%null(rank,rank))
					random%RandomVectorLL(1)%blocks(i,1)%null = 0
					do ii=1,rank
						random%RandomVectorLL(1)%blocks(i,1)%null(ii,ii)=1
					end do
				enddo                    
            else
                num_groupm=butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%num_row
                num_groupn=butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%num_col
                
				allocate (random%RandomVectorLL(level+1)%blocks(int(num_groupm/2),num_groupn))
				random%RandomVectorLL(level+1)%num_row=int(num_groupm/2)
				random%RandomVectorLL(level+1)%num_col=num_groupn                    
				do j=1, num_groupn
					index_j=int((j+1)/2)
					do i=1, num_groupm, 2
						index_i=int((i+1)/2)
						mm1=size(butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,1)
						mm2=size(butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j)%matrix,1)
						nn=size(butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,2)
						allocate (random%RandomVectorLL(level+1)%blocks(index_i,j)%matrix(nn,num_vect_sub))
						random%RandomVectorLL(level+1)%blocks(index_i,j)%matrix = 0
						random%RandomVectorLL(level+1)%blocks(index_i,j)%nulldim = nn
						allocate (random%RandomVectorLL(level+1)%blocks(index_i,j)%null(nn,nn))
						random%RandomVectorLL(level+1)%blocks(index_i,j)%null = 0
						do ii =1,nn
							random%RandomVectorLL(level+1)%blocks(index_i,j)%null(ii,ii) = 1
						end do
					enddo
				enddo
                
            endif
            if (level==level_butterfly) then
                allocate (random%RandomVectorLL(level+2)%blocks(1,num_blocks))
                random%RandomVectorLL(level+2)%num_row=1
                random%RandomVectorLL(level+2)%num_col=num_blocks
                do j=1, num_blocks
                    nn=size(butterfly_block_randomized(bb)%ButterflyV(j)%matrix,1)
                    rank=size(butterfly_block_randomized(bb)%ButterflyV(j)%matrix,2)
                    allocate (random%RandomVectorLL(level+2)%blocks(1,j)%matrix(nn,num_vect_sub))
					random%RandomVectorLL(level+2)%blocks(1,j)%matrix = 0
                enddo
            endif
        enddo
    
		mem_vec = 0
		do level=0, level_butterfly+2
			num_row=random%RandomVectorLL(level)%num_row
			num_col=random%RandomVectorLL(level)%num_col
			do j=1, num_col
				do i=1, num_row
					mem_vec =mem_vec + SIZEOF(random%RandomVectorLL(level)%blocks(i,j)%matrix)/1024.0d3
					if(level/=0 .and. level/=level_butterfly+2)mem_vec =mem_vec + SIZEOF(random%RandomVectorLL(level)%blocks(i,j)%null)/1024.0d3
				enddo
			enddo
		enddo
		Memory_int_vec = max(Memory_int_vec,mem_vec) 
	
	
    endif
    

    return
    
end subroutine Init_RandVects




subroutine Resolving_Butterflys_LL_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth,bb)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   
   integer nth_s,nth_e,unique_nth,bb
   integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,mm,kk,level_left,level_right, rs,re,rank,level_right_start,level_left_start
   integer index_i, index_j, iter, vector1, vector2, direction, round, flag
   real*8 a,b,c,d,norm1,norm2,norm3,norm4,norm1L,norm2L,norm3L,norm4L,norm1R,norm2R,norm3R,norm4R,error,errorL,errorR,rtemp,error0,error1,error2
   complex(kind=8) ctemp
   integer kmax
   
   type(RandomBlock), pointer :: random
   
   integer, allocatable :: ipiv(:), kernel_selection(:)
   real*8, allocatable :: Singular(:)
   complex(kind=8), allocatable :: matrixtemp(:,:), vectortemp(:), vectorstemp(:,:), tau(:), vectorstemp1(:,:),vectorstemp2(:,:)
   complex(kind=8), allocatable :: matrixtemp1(:,:),matA(:,:),matB(:,:),matC(:,:),matinv(:,:),matinv1(:,:),matinv2(:,:)
   integer num_vect_sub,num_vect_subsub,nth,ind_r,noe,Ng,dimension_nn,nn1,nn2,ieo,level_butterfly
   real*8::n1,n2
   type(butterflyblock_randomized), pointer :: blocks
   
   call assert(nth_e==nth_e,'Nbind/=1')
   
   blocks => butterfly_block_randomized(bb)   
   level_butterfly=blocks%level_butterfly 
   
   num_blocks=2**level_butterfly
   num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
	
   ! if (level_butterfly/=0) then
       mm=num_vect_subsub !rank

	   
	   level_right_start = floor_safe(level_butterfly/2d0)	!  check here later		   
	   ! ! level_right_start = level_butterfly+1
	   
	   
	   do level_right=0,unique_nth !level_right_start
		   ! kmax = ceiling_safe(rank/dble(2**(level_right_start-level_right)))+1
		   ! if(level_butterfly==9)write(*,*)level_right,kmax
			! write(*,*)level_right,'haha'
		   if (level_right==0) then 
			   do nth = nth_s,nth_e
				   !$omp parallel do default(shared) private(j)
				   do j=1, num_blocks
						call OneVs_LL(j,level_right,unique_nth,num_vect_sub,mm,nth,nth_s,bb)	   
					end do
					!$omp end parallel do
				end do
	
		   elseif (level_right==level_butterfly+1) then
				write(*,*)'the right half scheme should not touch leftmost matrix'
				stop
		   else

			   num_row=blocks%ButterflyKerl(level_right)%num_row
			   num_col=blocks%ButterflyKerl(level_right)%num_col
			   
			   do nth = nth_s,nth_e
				   noe = ceiling_safe(nth*Ng*2d0/num_col)
				   index_i = ceiling_safe(noe/2d0)
				   !$omp parallel do default(shared) private(j,index_j)
				   do j=1, num_col, 2
						index_j=int((j+1)/2)
						call OneKernels_LL(index_i, index_j,noe,level_right,unique_nth,num_vect_sub,mm,nth,nth_s,bb)
				   enddo
				   !$omp end parallel do
			   enddo
		   endif	   
	   end do
	   
   ! endif
   
   
   return

end subroutine Resolving_Butterflys_LL_new

subroutine OneVs_LL(j,level_right,unique_nth,num_vect_sub,mm,nth,nth_s,bb)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none 
   type(butterflyblock_randomized), pointer :: blocks
   complex(kind=8), allocatable :: matA(:,:),matB(:,:),matC(:,:),matinv(:,:)
   integer j,level_right,unique_nth,dimension_nn,mm,rank,num_vect_sub,nth,nth_s,level_butterfly,bb
	
   blocks => butterfly_block_randomized(bb)   
   level_butterfly=blocks%level_butterfly 
   
   if(level_right==unique_nth)then
	   dimension_nn=size(blocks%butterflyV(j)%matrix,1)
	   allocate(matB(mm,dimension_nn))
	   call copymatT_omp(random_Block(bb)%RandomVectorLL(level_butterfly+2)%blocks(1,j)%matrix(1:dimension_nn,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB,dimension_nn,mm)
	   call GetRank(mm,dimension_nn,matB,rank,Rank_detection_factor)
	   if(rank>blocks%dimension_rank)rank = blocks%dimension_rank
							   
	   if(allocated(blocks%butterflyV(j)%matrix))deallocate(blocks%butterflyV(j)%matrix)
	   ! if(allocated(blocks%ButterflyVInv(j)%matrix))deallocate(blocks%ButterflyVInv(j)%matrix)
	   if(allocated(random_Block(bb)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%matrix))deallocate(random_Block(bb)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%matrix)
	   allocate(blocks%butterflyV(j)%matrix(dimension_nn,rank))
	   ! allocate(blocks%ButterflyVInv(j)%matrix(rank,dimension_nn))
	   allocate(random_Block(bb)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%matrix(rank,num_vect_sub))
	   ! call RandomMat(rank,dimension_nn,min(rank,dimension_nn),blocks%ButterflyVInv(j)%matrix,0)
	   
	   allocate(matC(rank,dimension_nn),matA(mm,rank),matinv(dimension_nn,rank))
	   ! call copymatT_omp(blocks%ButterflyVInv(j)%matrix,matinv,rank,dimension_nn)
	   call copymatT_omp(blocks%KerInv(1:rank,1:dimension_nn),matinv,rank,dimension_nn)	
	   if(isnan(fnorm(matB,mm,dimension_nn)) .or. isnan(fnorm(matinv,dimension_nn,rank)))then
		write(*,*)shape(matB),fnorm(matB,mm,dimension_nn),fnorm(matinv,dimension_nn,rank),j,'hei',fnorm(random_Block(bb)%RandomVectorLL(level_butterfly+2)%blocks(1,j)%matrix,dimension_nn,mm)
		stop
	   end if
	   call gemm_omp(matB,matinv,matA,mm,dimension_nn,rank)
	   call LeastSquare(mm,rank,dimension_nn,matA,matB,matC,LS_tolerance)
	   ! write(*,*)fnorm(matC,rank,dimension_nn),'V',level_right,level_butterfly
	   call copymatT_omp(matC,blocks%ButterflyV(j)%matrix,rank,dimension_nn)						   
	   deallocate(matB,matC,matA,matinv)						   
   else 
	   rank=size(blocks%butterflyV(j)%matrix,2)
	   dimension_nn=size(blocks%butterflyV(j)%matrix,1)									
	   allocate(matB(mm,dimension_nn),matA(mm,rank),matinv(dimension_nn,rank))
	   call copymatT_omp(random_Block(bb)%RandomVectorLL(level_butterfly+2)%blocks(1,j)%matrix(1:dimension_nn,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB,dimension_nn,mm)
	   ! call copymatT_omp(blocks%ButterflyVInv(j)%matrix,matinv,rank,dimension_nn)
	   call copymatT_omp(blocks%KerInv(1:rank,1:dimension_nn),matinv,rank,dimension_nn)			
	   if(isnan(fnorm(matB,mm,dimension_nn)) .or. isnan(fnorm(matinv,dimension_nn,rank)))then
		write(*,*)fnorm(matB,mm,dimension_nn),fnorm(matinv,dimension_nn,rank),j,'hei1'
		stop
	   end if
	   call gemm_omp(matB,matinv,matA,mm,dimension_nn,rank)
	   
	   call copymatT_omp(matA,random_Block(bb)%RandomVectorLL(level_butterfly+1)%blocks(1,j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),mm,rank)					   
	   deallocate(matB,matA,matinv)	
   end if   
   
end subroutine OneVs_LL 	


subroutine OneKernels_LL(index_i, index_j,noe,level_right,unique_nth,num_vect_sub,mm,nth,nth_s,bb)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none 
   type(butterflyblock_randomized), pointer :: blocks
   complex(kind=8), allocatable :: matA(:,:),matB(:,:),matC(:,:),matinv(:,:),matinv1(:,:),matinv2(:,:)
   integer bb,index_i,index_j,i,j,level_right,unique_nth,dimension_nn,mm,rank,num_vect_sub,nth,nth_s,nn1,nn2,ieo,noe,rs,re,level_butterfly

   blocks => butterfly_block_randomized(bb)   
   level_butterfly=blocks%level_butterfly 
   
   i = index_i*2-1
   j = index_j*2-1
   ieo = i + 1 - mod(noe,2)

	nn1 = size(random_Block(bb)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j)%matrix,1)
	nn2 = size(random_Block(bb)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j+1)%matrix,1)

	if(level_right==unique_nth)then
		allocate (matB(mm,nn1+nn2))
		call copymatT_omp(random_Block(bb)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j)%matrix(1:nn1,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1:nn1),nn1,mm)
		call copymatT_omp(random_Block(bb)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j+1)%matrix(1:nn2,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1+nn1:nn2+nn1),nn2,mm)
		if(mod(noe,2)==1)then
			call GetRank(mm,nn1+nn2,matB,rank,Rank_detection_factor)
			if(rank>blocks%dimension_rank)rank = blocks%dimension_rank
			
			rs = 1
			re = rank
		else 
			call GetRank(mm,nn1+nn2,matB,rank,Rank_detection_factor)
			if(rank>blocks%dimension_rank)rank = blocks%dimension_rank
									   
			rs = size(blocks%ButterflyKerl(level_right)%blocks(i,j)%matrix,1)+1
			re = rs+rank-1		
		end if


		if(allocated(blocks%ButterflyKerl(level_right)%blocks(ieo,j)%matrix))deallocate(blocks%ButterflyKerl(level_right)%blocks(ieo,j)%matrix)
		if(allocated(blocks%ButterflyKerl(level_right)%blocks(ieo,j+1)%matrix))deallocate(blocks%ButterflyKerl(level_right)%blocks(ieo,j+1)%matrix)
		allocate(blocks%ButterflyKerl(level_right)%blocks(ieo,j)%matrix(rank,nn1))
		allocate(blocks%ButterflyKerl(level_right)%blocks(ieo,j+1)%matrix(rank,nn2))
		if(allocated(random_Block(bb)%RandomVectorLL(level_butterfly-level_right+1)%blocks(ieo,index_j)%matrix))deallocate(random_Block(bb)%RandomVectorLL(level_butterfly-level_right+1)%blocks(ieo,index_j)%matrix)
		allocate(random_Block(bb)%RandomVectorLL(level_butterfly-level_right+1)%blocks(ieo,index_j)%matrix(rank,mm))
		
		
		allocate (matC(rank,nn1+nn2),matA(mm,rank),matinv(nn1+nn2,rank))
		call copymatN_omp(blocks%KerInv(1:nn1+nn2,rs:re),matinv,nn1+nn2,rank)			
		if(isnan(fnorm(matB,mm,nn1+nn2)) .or. isnan(fnorm(matinv,nn1+nn2,rank)))then
		 write(*,*)fnorm(matB,mm,nn1+nn2),fnorm(matinv,nn1+nn2,rank),i,j,'ho'
		 stop
	    end if
		call gemm_omp(matB,matinv,matA,mm,nn1+nn2,rank)
		call LeastSquare(mm,rank,nn1+nn2,matA,matB,matC,LS_tolerance)
		! write(*,*)fnorm(matC,rank,nn1+nn2),'KernelL',level_right,level_butterfly
		call copymatN_omp(matC(1:rank,1:nn1),blocks%ButterflyKerl(level_right)%blocks(ieo,j)%matrix,rank,nn1)
		call copymatN_omp(matC(1:rank,nn1+1:nn1+nn2),blocks%ButterflyKerl(level_right)%blocks(ieo,j+1)%matrix,rank,nn2)							
		deallocate(matB,matC,matA,matinv)
	else 
		if(mod(noe,2)==1)then
			rank = size(blocks%ButterflyKerl(level_right)%blocks(i,j)%matrix,1)
			rs = 1
			re = rank
		else 
			rank = size(blocks%ButterflyKerl(level_right)%blocks(i+1,j)%matrix,1)
			rs = size(blocks%ButterflyKerl(level_right)%blocks(i,j)%matrix,1)+1
			re = rs+rank-1
		end if
		allocate (matB(mm,nn1+nn2),matC(rank,nn1+nn2),matA(mm,rank),matinv(nn1+nn2,rank))
		call copymatT_omp(random_Block(bb)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j)%matrix(1:nn1,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1:nn1),nn1,mm)
		call copymatT_omp(random_Block(bb)%RandomVectorLL(level_butterfly-level_right+2)%blocks(index_i,j+1)%matrix(1:nn2,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1+nn1:nn2+nn1),nn2,mm)
										
		call copymatN_omp(blocks%KerInv(1:nn1+nn2,rs:re),matinv,nn1+nn2,rank)	
		if(isnan(fnorm(matB,mm,nn1+nn2)) .or. isnan(fnorm(matinv,nn1+nn2,rank)))then
		 write(*,*)fnorm(matB,mm,nn1+nn2),fnorm(matinv,nn1+nn2,rank),i,j,'ho1'
		 stop
	    end if		
		call gemm_omp(matB,matinv,matA,mm,nn1+nn2,rank)
		call copymatT_omp(matA,random_Block(bb)%RandomVectorLL(level_butterfly-level_right+1)%blocks(ieo,index_j)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),mm,rank)
		deallocate(matB,matC,matA,matinv)									
	end if
		! write(*,'(I5,I5,I5,I5,I5,Es16.7E3,A2)')unique_nth,level_right,nth,i,j,error0,'L' 


   
end subroutine OneKernels_LL   







subroutine Resolving_Butterflys_RR_new(num_vect_sub,nth_s,nth_e,Ng,unique_nth,bb)

   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none
   
   integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,level_left,level_right, rank,level_right_start,level_left_start,nn1,nn2,rs,re
   integer index_i, index_j, mm, nn, iter, vector1, vector2, direction, round, flag, bb
   real*8 a,b,c,d,norm1,norm2,norm3,norm4,norm1L,norm2L,norm3L,norm4L,norm1R,norm2R,norm3R,norm4R,error,errorL,errorR,rtemp,error0,error1,error2
   complex(kind=8) ctemp
   
   type(butterflyblock_randomized), pointer :: blocks
   type(RandomBlock), pointer :: random
   
   integer, allocatable :: ipiv(:), kernel_selection(:)
   real*8, allocatable :: Singular(:)
   complex(kind=8), allocatable :: matrixtemp(:,:), vectortemp(:), vectorstemp(:,:), tau(:), vectorstemp1(:,:),vectorstemp2(:,:)
   complex(kind=8), allocatable :: matrixtemp1(:,:),matA(:,:),matB(:,:),matC(:,:),matinv(:,:),matinv1(:,:),matinv2(:,:)
   integer num_vect_sub,num_vect_subsub,nth,ind_r,ind_c,noe,Ng,nth_s,nth_e,dimension_mm,dimension_n,jeo,level_butterfly
   real*8::n1,n2
   integer::kmax,unique_nth
   
   call assert(nth_e==nth_e,'Nbind/=1')
   
   blocks => butterfly_block_randomized(bb)   
   level_butterfly=blocks%level_butterfly 

   num_blocks=2**level_butterfly
   num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
   
   ! if (level_butterfly/=0) then
       mm=num_vect_subsub !rank

	   level_left_start= floor_safe(level_butterfly/2d0)+1    !  check here later		   
	   ! ! level_left_start = 0
	   
	   random=>random_Block(bb)
	   if(level_left_start>0 .and. level_left_start==unique_nth)then
			n1 = OMP_get_wtime()
			call Butterflys_partial_MVP_Half('N',0,level_left_start-1,random,num_vect_sub,nth_s,nth_e,Ng,bb)
			! call Butterfly_partial_MVP(blocks,'N',0,level_left_start-1,random,num_vect_sub)
			n2 = OMP_get_wtime()
			time_halfbuttermul = time_halfbuttermul + n2-n1		 
		endif 
	   
	   
	   do level_left = level_butterfly+1,unique_nth,-1
			if (level_left==level_butterfly+1) then
				do nth=nth_s,nth_e
					!$omp parallel do default(shared) private(i)
					do i=1, num_blocks   
						call OneUs_RR(i,level_left,unique_nth,num_vect_sub,mm,nth,nth_s,bb)
					end do
					!$omp end parallel do					
				end do
			elseif (level_left==0) then
				write(*,*)'the left half scheme should not touch rightmost matrix'
				stop
			else 
				! write(*,*)'good'
			   num_row=blocks%ButterflyKerl(level_left)%num_row
			   num_col=blocks%ButterflyKerl(level_left)%num_col
			   do nth = nth_s,nth_e
				   noe = ceiling_safe(nth*Ng*2d0/num_row)
				   index_j = ceiling_safe(noe/2d0)
				   !$omp parallel do default(shared) private(i,index_i)
				   do i=1, num_row, 2
					   index_i=int((i+1)/2)
					   call OneKernels_RR(index_i, index_j,noe,level_left,level_left_start,unique_nth,num_vect_sub,mm,nth,nth_s,bb)									
				   enddo
				   !$omp end parallel do	
			   enddo
				! write(*,*)'good1'			   
			end if
	   end do
   ! endif
   
   
   return

end subroutine Resolving_Butterflys_RR_new


subroutine OneUs_RR(i,level_left,unique_nth,num_vect_sub,mm,nth,nth_s,bb)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none 
   type(butterflyblock_randomized), pointer :: blocks
   complex(kind=8), allocatable :: matA(:,:),matB(:,:),matC(:,:),matinv(:,:)
   integer i,level_left,unique_nth,dimension_mm,mm,rank,num_vect_sub,nth,nth_s,level_butterfly,bb
	
   blocks => butterfly_block_randomized(bb)   
   level_butterfly=blocks%level_butterfly 
   
	if(level_left==unique_nth)then

		if(level_butterfly>0)then
			dimension_mm=size(blocks%butterflyU(i)%matrix,1)	
			allocate(matB(mm,dimension_mm))
			call copymatT_omp(random_Block(bb)%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(1:dimension_mm,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB,dimension_mm,mm)							
			call GetRank(mm,dimension_mm,matB,rank,Rank_detection_factor)
			if(rank>blocks%dimension_rank)rank = blocks%dimension_rank
			
			if(allocated(blocks%butterflyU(i)%matrix))deallocate(blocks%butterflyU(i)%matrix)
			if(allocated(random_Block(bb)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix))deallocate(random_Block(bb)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix)
			allocate(blocks%butterflyU(i)%matrix(dimension_mm,rank))
			allocate(random_Block(bb)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix(rank,num_vect_sub))
			allocate(matC(rank,dimension_mm),matA(mm,rank),matinv(dimension_mm,rank))	
			call copymatT_omp(blocks%KerInv(1:rank,1:dimension_mm),matinv,rank,dimension_mm)
			if(isnan(fnorm(matB,mm,dimension_mm)) .or. isnan(fnorm(matinv,dimension_mm,rank)))then
			 write(*,*)fnorm(matB,mm,dimension_mm),fnorm(matinv,dimension_mm,rank),i,'heee'
			 stop
			end if		
			call gemm_omp(matB,matinv,matA,mm,dimension_mm,rank)							
			call LeastSquare(mm,rank,dimension_mm,matA,matB,matC,LS_tolerance)
			! write(*,*)fnorm(matC,rank,dimension_mm),'U',level_left,level_butterfly
			call copymatT_omp(matC,blocks%ButterflyU(i)%matrix,rank,dimension_mm)							
			deallocate(matB,matC,matA,matinv)
		else 
			dimension_mm=size(blocks%butterflyU(i)%matrix,1)	
			allocate(matB(mm,dimension_mm))
			call copymatT_omp(random_Block(bb)%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(1:dimension_mm,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB,dimension_mm,mm)									
			rank = size(random_Block(bb)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix,1)
			if(allocated(blocks%butterflyU(i)%matrix))deallocate(blocks%butterflyU(i)%matrix)
			allocate(blocks%butterflyU(i)%matrix(dimension_mm,rank))
			allocate(matC(rank,dimension_mm),matA(mm,rank))	
			call copymatT_omp(random_Block(bb)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matA,rank,mm)
			call LeastSquare(mm,rank,dimension_mm,matA,matB,matC,LS_tolerance)
			! write(*,*)fnorm(matC,rank,dimension_mm),'U',level_left,level_butterfly
			call copymatT_omp(matC,blocks%ButterflyU(i)%matrix,rank,dimension_mm)							
			deallocate(matB,matC,matA)			
		endif		
	else 
		dimension_mm=size(blocks%butterflyU(i)%matrix,1)						
		rank=size(blocks%butterflyU(i)%matrix,2)						
		allocate(matB(mm,dimension_mm),matA(mm,rank),matinv(dimension_mm,rank))	
		call copymatT_omp(random_Block(bb)%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(1:dimension_mm,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB,dimension_mm,mm)
		call copymatT_omp(blocks%KerInv(1:rank,1:dimension_mm),matinv,rank,dimension_mm)	
		if(isnan(fnorm(matB,mm,dimension_mm)) .or. isnan(fnorm(matinv,dimension_mm,rank)))then
		 write(*,*)fnorm(matB,mm,dimension_mm),fnorm(matinv,dimension_mm,rank),i,'heee1'
		 stop
	    end if			
		call gemm_omp(matB,matinv,matA,mm,dimension_mm,rank)
		call copymatT_omp(matA,random_Block(bb)%RandomVectorRR(level_butterfly+1)%blocks(i,1)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),mm,rank)
		deallocate(matB,matA,matinv)					
	end if	   
   

end subroutine OneUs_RR  



subroutine OneKernels_RR(index_i, index_j,noe,level_left,level_left_start,unique_nth,num_vect_sub,mm,nth,nth_s,bb)
   use MODULE_FILE
   ! use lapack95
   ! use blas95
   implicit none 
   type(butterflyblock_randomized), pointer :: blocks
   complex(kind=8), allocatable :: matA(:,:),matB(:,:),matC(:,:),matinv(:,:),matinv1(:,:),matinv2(:,:)
   integer index_i,index_j,i,j,bb,level_left,unique_nth,dimension_nn,mm,rank,num_vect_sub,nth,nth_s,nn1,nn2,jeo,noe,rs,re,level_left_start,level_butterfly

   blocks => butterfly_block_randomized(bb)   
   level_butterfly=blocks%level_butterfly 
	
	i = index_i*2-1
	j = index_j*2-1
	jeo = j + 1 - mod(noe,2)					   
	
	nn1 = size(random_Block(bb)%RandomVectorRR(level_left+1)%blocks(i,index_j)%matrix,1)
	nn2 = size(random_Block(bb)%RandomVectorRR(level_left+1)%blocks(i+1,index_j)%matrix,1)

	if(level_left==unique_nth)then
		if(level_left==level_left_start)then
			allocate (matB(mm,nn1+nn2))
			call copymatT_omp(random_Block(bb)%RandomVectorRR(level_left+1)%blocks(i,index_j)%matrix(1:nn1,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1:nn1),nn1,mm)
			call copymatT_omp(random_Block(bb)%RandomVectorRR(level_left+1)%blocks(i+1,index_j)%matrix(1:nn2,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1+nn1:nn2+nn1),nn2,mm)								
			
			if(mod(noe,2)==1)then
				rank = size(random_Block(bb)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix,1)
				rs = 1
				re = rank
			else 
				rank = size(random_Block(bb)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix,1)
				rs = size(blocks%ButterflyKerl(level_left)%blocks(i,j)%matrix,2)+1
				re = rs+rank-1		
			end if																		

			if(allocated(blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix))deallocate(blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix)
			if(allocated(blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix))deallocate(blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix)
			allocate(blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix(nn1,rank))
			allocate(blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix(nn2,rank))
			allocate(matC(rank,nn1+nn2),matA(mm,rank))
			call copymatT_omp(random_Block(bb)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matA,rank,mm)
			call LeastSquare(mm,rank,nn1+nn2,matA,matB,matC,LS_tolerance)
			! write(*,*)fnorm(matC,rank,nn1+nn2),'KernelR',level_left,level_butterfly
			call copymatT_omp(matC(1:rank,1:nn1),blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix,rank,nn1)
			call copymatT_omp(matC(1:rank,nn1+1:nn1+nn2),blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix,rank,nn2)										
			deallocate(matB,matC,matA)

		else 
			allocate (matB(mm,nn1+nn2))
			call copymatT_omp(random_Block(bb)%RandomVectorRR(level_left+1)%blocks(i,index_j)%matrix(1:nn1,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1:nn1),nn1,mm)
			call copymatT_omp(random_Block(bb)%RandomVectorRR(level_left+1)%blocks(i+1,index_j)%matrix(1:nn2,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1+nn1:nn2+nn1),nn2,mm)								
			if(mod(noe,2)==1)then
				call GetRank(mm,nn1+nn2,matB,rank,Rank_detection_factor)
				if(rank>blocks%dimension_rank)rank = blocks%dimension_rank
				
				rs = 1
				re = rank
			else 
				call GetRank(mm,nn1+nn2,matB,rank,Rank_detection_factor)
				if(rank>blocks%dimension_rank)rank = blocks%dimension_rank
				
				rs = size(blocks%ButterflyKerl(level_left)%blocks(i,j)%matrix,2)+1
				re = rs+rank-1		
			end if																		
			if(allocated(blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix))deallocate(blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix)
			if(allocated(blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix))deallocate(blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix)
			allocate(blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix(nn1,rank))
			allocate(blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix(nn2,rank))
			if(allocated(random_Block(bb)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix))deallocate(random_Block(bb)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix)
			allocate(random_Block(bb)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix(rank,mm))
			
			allocate(matC(rank,nn1+nn2),matA(mm,rank),matinv(nn1+nn2,rank))
			call copymatT_OMP(blocks%KerInv(rs:re,1:nn1+nn2),matinv,rank,nn1+nn2)	
			if(isnan(fnorm(matB,mm,nn1+nn2)) .or. isnan(fnorm(matinv,nn1+nn2,rank)))then
			 write(*,*)fnorm(matB,mm,nn1+nn2),fnorm(matinv,nn1+nn2,rank),i,j,'helo'
			 stop
			end if	
			call gemm_omp(matB,matinv,matA,mm,nn1+nn2,rank)
			
			call LeastSquare(mm,rank,nn1+nn2,matA,matB,matC,LS_tolerance)
			! write(*,*)fnorm(matC,rank,nn1+nn2),'KernelR',level_left,level_butterfly
			call copymatT_omp(matC(1:rank,1:nn1),blocks%ButterflyKerl(level_left)%blocks(i,jeo)%matrix,rank,nn1)
			call copymatT_omp(matC(1:rank,nn1+1:nn1+nn2),blocks%ButterflyKerl(level_left)%blocks(i+1,jeo)%matrix,rank,nn2)										
			deallocate(matB,matC,matA,matinv)
			
		end if
	else 
		if(mod(noe,2)==1)then
			rank = size(blocks%ButterflyKerl(level_left)%blocks(i,j)%matrix,2)
			rs = 1
			re = rank
		else 
			rank = size(blocks%ButterflyKerl(level_left)%blocks(i,j+1)%matrix,2)
			rs = size(blocks%ButterflyKerl(level_left)%blocks(i,j)%matrix,2)+1
			re = rs+rank-1
		end if							
		allocate (matB(mm,nn1+nn2),matA(mm,rank),matinv(nn1+nn2,rank))
		call copymatT_omp(random_Block(bb)%RandomVectorRR(level_left+1)%blocks(i,index_j)%matrix(1:nn1,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1:nn1),nn1,mm)
		call copymatT_omp(random_Block(bb)%RandomVectorRR(level_left+1)%blocks(i+1,index_j)%matrix(1:nn2,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),matB(1:mm,1+nn1:nn2+nn1),nn2,mm)								
		
		call copymatT_omp(blocks%KerInv(rs:re,1:nn1+nn2),matinv,rank,nn1+nn2)			
		if(isnan(fnorm(matB,mm,nn1+nn2)) .or. isnan(fnorm(matinv,nn1+nn2,rank)))then
		 write(*,*)fnorm(matB,mm,nn1+nn2),fnorm(matinv,nn1+nn2,rank),i,j,'helo'
		 stop
		end if			
		call gemm_omp(matB,matinv,matA,mm,nn1+nn2,rank)
		call copymatT_omp(matA,random_Block(bb)%RandomVectorRR(level_left)%blocks(index_i,jeo)%matrix(1:rank,(nth-nth_s)*mm+1:(nth-nth_s+1)*mm),mm,rank)
		deallocate(matB,matA,matinv)	
	end if

	! write(*,'(I5,I5,I5,I5,I5,Es16.7E3,A2)')unique_nth,level_left,nth,i,j,error0,'R'
end subroutine OneKernels_RR  




subroutine Butterflys_Partial_MVP_Half(chara,level_start,level_end,random,num_vect_sub,nth_s,nth_e,Ng,bb)
    
    use MODULE_FILE
    implicit none
    
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, ij, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2, bb
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2, level_end, level_start
    complex(kind=8) ctemp, a, b
    character chara
	integer num_vect_sub,num_vect_subsub,nth_s,nth_e,Ng,nth,dimension_rank,level_butterfly
    
    type(RandomBlock) :: random
    
    type(butterfly_Kerl),allocatable :: ButterflyVector(:)
   
	level_butterfly=butterfly_block_randomized(bb)%level_butterfly
	dimension_rank = butterfly_block_randomized(bb)%dimension_rank
    num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)
   
    if (chara=='N') then

        num_blocks=2**level_butterfly
        
        do level=level_start, level_end
            if (level==0) then
                num_groupn=num_blocks
                do nth=nth_s, nth_e
					!$omp parallel do default(shared) private(j,rank,nn,ii,jj,ctemp,kk)
					do j = (nth-1)*Ng+1,nth*Ng
						rank=size(butterfly_block_randomized(bb)%ButterflyV(j)%matrix,2)
						nn=size(butterfly_block_randomized(bb)%ButterflyV(j)%matrix,1)
						random%RandomVectorRR(1)%blocks(1,j)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
						! !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
						do ii=1, rank
							do jj=1, num_vect_subsub
								ctemp=(0.,0.)
								do kk=1, nn
									ctemp=ctemp+butterfly_block_randomized(bb)%ButterflyV(j)%matrix(kk,ii)*random%RandomVectorRR(0)%blocks(1,j)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
								enddo
								random%RandomVectorRR(1)%blocks(1,j)%matrix(ii,jj+(nth-nth_s)*num_vect_subsub)=ctemp
							enddo
						enddo
						! !$omp end parallel do
						
						! write(*,*)level,fnorm(random%RandomVectorRR(0)%blocks(1,j)%matrix,nn,num_vect_subsub),fnorm(random%RandomVectorRR(1)%blocks(1,j)%matrix,rank,num_vect_subsub),fnorm(butterfly_block_randomized(bb)%ButterflyV(j)%matrix,nn,rank),'hi'
							
						
						
					enddo
					!$omp end parallel do					
                enddo
            elseif (level==level_butterfly+1) then
                                 
            else
                num_groupm=butterfly_block_randomized(bb)%ButterflyKerl(level)%num_row
                num_groupn=butterfly_block_randomized(bb)%ButterflyKerl(level)%num_col
					
				!$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,nn1,nn2,mm,nth)
				do ij=1,num_groupm*(num_groupn/2)
					i = (ij-1)/(num_groupn/2)+1
					j = (mod(ij-1,(num_groupn/2)) + 1)*2-1
					index_i=int((i+1)/2)
					index_j=int((j+1)/2)
					
					nn1=size(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix,2)
					nn2=size(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j+1)%matrix,2)
					mm=size(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix,1)
					do nth = nth_s,nth_e
						random%RandomVectorRR(level+1)%blocks(i,index_j)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
						if((j>=(nth-1)*Ng/2**(level-1)+1 .and. j<=nth*Ng/2**(level-1)) .or. &
						& (j+1>=(nth-1)*Ng/2**(level-1)+1 .and. j+1<=nth*Ng/2**(level-1)))then						
							! !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
							do ii=1, mm
								do jj=1, num_vect_subsub
									ctemp=(0.,0.)
									do kk=1, nn1
										ctemp=ctemp+butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix(ii,kk)*random%RandomVectorRR(level)%blocks(index_i,j)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
									enddo
									do kk=1, nn2
										ctemp=ctemp+butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j+1)%matrix(ii,kk)*random%RandomVectorRR(level)%blocks(index_i,j+1)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
									enddo
									random%RandomVectorRR(level+1)%blocks(i,index_j)%matrix(ii,jj+(nth-nth_s)*num_vect_subsub)=ctemp
								enddo
							enddo
							! !$omp end parallel do
							
							! write(*,*)level,fnorm(random%RandomVectorRR(level+1)%blocks(i,index_j)%matrix,mm,num_vect_subsub),fnorm(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix,mm,nn1),fnorm(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j+1)%matrix,mm,nn2),'hi'
							
							
						end if	
					end do	
				enddo
				!$omp end parallel do
            endif
        enddo      
        
                    
    elseif (chara=='T') then
    
        num_blocks=2**level_butterfly

        do level=level_start, level_end
            if (level==0) then
                num_groupm=num_blocks
                do nth=nth_s, nth_e
					!$omp parallel do default(shared) private(i,rank,mm,ii,jj,ctemp,kk)
					do i = (nth-1)*Ng+1,nth*Ng
						rank=size(butterfly_block_randomized(bb)%ButterflyU(i)%matrix,2)
						mm=size(butterfly_block_randomized(bb)%ButterflyU(i)%matrix,1)
						random%RandomVectorLL(1)%blocks(i,1)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
						! !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
						do ii=1, rank
							do jj=1, num_vect_subsub
								ctemp=(0.,0.)
								do kk=1, mm
									ctemp=ctemp+butterfly_block_randomized(bb)%ButterflyU(i)%matrix(kk,ii)*random%RandomVectorLL(0)%blocks(i,1)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
								enddo
								random%RandomVectorLL(1)%blocks(i,1)%matrix(ii,jj+(nth-nth_s)*num_vect_subsub)=ctemp
							enddo
						enddo
						! !$omp end parallel do
					end do
					!$omp end parallel do
                enddo
            elseif (level==level_butterfly+1) then               
            else
                num_groupm=butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%num_row
                num_groupn=butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%num_col             

				!$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,mm1,mm2,nn,nth)
				do ij=1,num_groupn*(num_groupm/2)
					j = (ij-1)/(num_groupm/2)+1
					i = (mod(ij-1,(num_groupm/2)) + 1)*2-1	
					index_j=int((j+1)/2)
					index_i=int((i+1)/2)						
					
					mm1=size(butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,1)
					mm2=size(butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j)%matrix,1)
					nn=size(butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,2)
					do nth = nth_s,nth_e
						random%RandomVectorLL(level+1)%blocks(index_i,j)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
						if((i>=(nth-1)*Ng/2**(level-1)+1 .and. i<=nth*Ng/2**(level-1)) .or. &
						& (i+1>=(nth-1)*Ng/2**(level-1)+1 .and. i+1<=nth*Ng/2**(level-1)))then
							! !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
							do jj=1, nn
								do ii=1, num_vect_subsub
									ctemp=(0.,0.)
									do kk=1, mm1
										ctemp=ctemp+random%RandomVectorLL(level)%blocks(i,index_j)%matrix(kk,ii+(nth-nth_s)*num_vect_subsub)*butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix(kk,jj)
									enddo
									do kk=1, mm2
										ctemp=ctemp+random%RandomVectorLL(level)%blocks(i+1,index_j)%matrix(kk,ii+(nth-nth_s)*num_vect_subsub)*butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j)%matrix(kk,jj)
									enddo
									random%RandomVectorLL(level+1)%blocks(index_i,j)%matrix(jj,ii+(nth-nth_s)*num_vect_subsub)=ctemp
								enddo
							enddo
							! !$omp end parallel do
						end if
					end do
				enddo
				!$omp end parallel do

            endif
        enddo
    
    endif
    
    return
    
end subroutine Butterflys_Partial_MVP_Half




subroutine Butterflys_block_MVP_dat(chara,M,N,Nrnd,random1,random2,a,b,bb)
    
    use MODULE_FILE
	use misc
    implicit none
    
    integer M,N, Nrnd,index_i, index_j, na, nb, index_start, num_vectors
    integer i, j, ii, jj, ij, level, level_butterfly, index_iijj, index_ij, k, k1, k2, kk, intemp1, intemp2
    integer vector_inuse, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, butterflyB_inuse, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, mm1, mm2,levelm
    complex(kind=8) ctemp, a, b
    character chara
	! type(matrixblock)::blocks
    integer:: middleflag,bb
	
    type(butterfly_Kerl),allocatable :: ButterflyVector(:)
    !  complex(kind=8) :: random1(N,Nrnd), random2(M,Nrnd)
        complex(kind=8) :: random1(:,:), random2(:,:)
	complex(kind=8),allocatable::matrixtemp(:,:),matrixtemp1(:,:)
	!  write(*,*)'nima-1'
	integer,allocatable:: arr_acc_m(:),arr_acc_n(:)
	
	middleflag = 0
	if(allocated(butterfly_block_randomized(bb)%ButterflyMiddle))middleflag=1
	
    level_butterfly=butterfly_block_randomized(bb)%level_butterfly
    num_blocks=2**level_butterfly	
	allocate(arr_acc_m(num_blocks))
	allocate(arr_acc_n(num_blocks))
	
	k1=0
	k2=0
	do i=1, num_blocks
		arr_acc_n(i) = k1
		arr_acc_m(i) = k2
		nn=size(butterfly_block_randomized(bb)%ButterflyV(i)%matrix,1)
		k1 =k1 +nn
		mm=size(butterfly_block_randomized(bb)%ButterflyU(i)%matrix,1)
		k2 =k2 +mm		
	enddo 
	
    num_vectors=Nrnd
    ! write(*,*)num_vectors
	! stop
	
	! if(CheckNAN_Butterfly(blocks))then
			! write(*,*)'NAN in 0 butterfly_block_MVP_randomized_dat'
			! stop		
	! end if
	
    if (chara=='N') then
    
		if(isnan(sum(abs(random1(:,1))**2)))then
			write(*,*)'NAN in 1 butterfly_block_MVP_randomized_dat'
			stop
		end if	
	
        level_butterfly=butterfly_block_randomized(bb)%level_butterfly
        num_blocks=2**level_butterfly
	    levelm = ceiling_safe(dble(level_butterfly)/2d0)        
      
        allocate (ButterflyVector(0:level_butterfly+2))
        allocate (ButterflyVector(0)%blocks(1,num_blocks))
        ButterflyVector(0)%num_row=1
        ButterflyVector(0)%num_col=num_blocks
                !  write(*,*)'nima0'
		!$omp parallel do default(shared) private(i,nn,ii,jj)
		do i=1, num_blocks
            nn=size(butterfly_block_randomized(bb)%ButterflyV(i)%matrix,1)
            allocate (ButterflyVector(0)%blocks(1,i)%matrix(nn,num_vectors))
            do ii=1, nn
                do jj=1, num_vectors
                    ButterflyVector(0)%blocks(1,i)%matrix(ii,jj)=random1(ii+arr_acc_n(i),jj)
                enddo
            enddo
		enddo 
		!$omp end parallel do
		
        !  write(*,*)'nima1'
        do level=0, level_butterfly
	        !  write(*,*)'nima1',level
            if (level==0) then
                num_groupn=num_blocks
                allocate (ButterflyVector(1)%blocks(1,num_groupn))
                ButterflyVector(1)%num_row=1
                ButterflyVector(1)%num_col=num_groupn
				
                !$omp parallel do default(shared) private(j,rank,nn,ii,jj,ctemp,kk)
				do j=1, num_groupn
				! write(*,*)num_groupn
                    rank=size(butterfly_block_randomized(bb)%ButterflyV(j)%matrix,2)
                    nn=size(butterfly_block_randomized(bb)%ButterflyV(j)%matrix,1)
                    allocate (ButterflyVector(1)%blocks(1,j)%matrix(rank,num_vectors))
                    ! !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
                    do ii=1, rank
                        do jj=1, num_vectors
                            ctemp=0
                            do kk=1, nn
                                ctemp=ctemp+butterfly_block_randomized(bb)%ButterflyV(j)%matrix(kk,ii)*ButterflyVector(0)%blocks(1,j)%matrix(kk,jj)
                            enddo
                            ButterflyVector(1)%blocks(1,j)%matrix(ii,jj)=ctemp
                        enddo
                    enddo
                    ! !$omp end parallel do
					
					! write(*,*)shape(butterfly_block_randomized(bb)%ButterflyV(j)%matrix),'da'
					! ! if(is_nan_mat_c(butterfly_block_randomized(bb)%ButterflyV(j)%matrix,nn,rank))then
						! ! write(*,*)'V'
						! ! stop
					! ! end if					
                enddo  
				!$omp end parallel do
				
            else
                num_groupm=butterfly_block_randomized(bb)%ButterflyKerl(level)%num_row
                num_groupn=butterfly_block_randomized(bb)%ButterflyKerl(level)%num_col
                if (num_groupn/=1) then
                    allocate (ButterflyVector(level+1)%blocks(num_groupm,int(num_groupn/2)))
                    ButterflyVector(level+1)%num_row=num_groupm
                    ButterflyVector(level+1)%num_col=int(num_groupn/2)                    
                    
					
                    !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,nn1,nn2,mm)
					do ij=1,num_groupm*(num_groupn/2)
						i = (ij-1)/(num_groupn/2)+1
						j = (mod(ij-1,(num_groupn/2)) + 1)*2-1
						index_i=int((i+1)/2)
						index_j=int((j+1)/2)					

						nn1=size(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix,2)
						nn2=size(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j+1)%matrix,2)
						mm=size(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix,1)
						allocate (ButterflyVector(level+1)%blocks(i,index_j)%matrix(mm,num_vectors))
						if(size(ButterflyVector(level)%blocks(index_i,j)%matrix,1)<nn1)then
							write(*,*)butterfly_block_randomized(bb)%level_butterfly,level,'nimade'
							stop
						end if
						! !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
						do ii=1, mm
							do jj=1, num_vectors
								ctemp=0
								do kk=1, nn1
									ctemp=ctemp+butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix(ii,kk)*ButterflyVector(level)%blocks(index_i,j)%matrix(kk,jj)
								enddo
								do kk=1, nn2
									ctemp=ctemp+butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j+1)%matrix(ii,kk)*ButterflyVector(level)%blocks(index_i,j+1)%matrix(kk,jj)
								enddo
								ButterflyVector(level+1)%blocks(i,index_j)%matrix(ii,jj)=ctemp
							enddo
						enddo
						! !$omp end parallel do
				! ! if(is_nan_mat_c(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix,mm,nn1) .or. is_nan_mat_c(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j+1)%matrix,mm,nn2))then
					! ! write(*,*)'kernel'
					! ! stop
				! ! end if
						if(level==levelm .and. middleflag==1 .and. level_butterfly>=2)then	
							! allocate(matrixtemp(mm,num_vectors))
							call gemm_omp(butterfly_block_randomized(bb)%ButterflyMiddle(i,index_j)%matrix,ButterflyVector(level+1)%blocks(i,index_j)%matrix,ButterflyVector(level+1)%blocks(i,index_j)%matrix,mm,mm,num_vectors)
							! ButterflyVector(level+1)%blocks(i,index_j)%matrix = matrixtemp
							! deallocate(matrixtemp)
						end if					
                    enddo
					!$omp end parallel do
                else
                    allocate (ButterflyVector(level+1)%blocks(num_groupm,1))
                    ButterflyVector(level+1)%num_row=num_groupm
                    ButterflyVector(level+1)%num_col=1
                    do i=1, num_groupm
                        index_i=int((i+1)/2)
                        nn=size(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,1)%matrix,2)
                        mm=size(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,1)%matrix,1)
                        allocate (ButterflyVector(level+1)%blocks(i,1)%matrix(mm,num_vectors))
                        !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
                        do ii=1, mm
                            do jj=1, num_vectors
                                ctemp=0                       
                                do kk=1, nn
                                    ctemp=ctemp+butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,1)%matrix(ii,kk)*ButterflyVector(level)%blocks(index_i,1)%matrix(kk,jj)
                                enddo
                                ButterflyVector(level+1)%blocks(i,1)%matrix(ii,jj)=ctemp
                            enddo
                        enddo
                        !$omp end parallel do
					! ! if(is_nan_mat_c(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,1)%matrix,mm,nn))then
						! ! write(*,*)'kernel2'
						! ! stop
					! ! end if
													
						
						
                    enddo
                endif
            endif
            if (level==level_butterfly) then
                allocate (ButterflyVector(level+2)%blocks(num_blocks,1))
                ButterflyVector(level+2)%num_row=num_blocks
                ButterflyVector(level+2)%num_col=1
                !$omp parallel do default(shared) private(i,rank,mm,ii,jj,ctemp,kk)
				do i=1, num_blocks
                    rank=size(butterfly_block_randomized(bb)%ButterflyU(i)%matrix,2)
                    mm=size(butterfly_block_randomized(bb)%ButterflyU(i)%matrix,1)
                    allocate (ButterflyVector(level+2)%blocks(i,1)%matrix(mm,num_vectors))
                    ! !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
                    do ii=1, mm
                        do jj=1, num_vectors
                            ctemp=0
                            do kk=1, rank
                                ctemp=ctemp+butterfly_block_randomized(bb)%ButterflyU(i)%matrix(ii,kk)*ButterflyVector(level+1)%blocks(i,1)%matrix(kk,jj)
                            enddo
                            ButterflyVector(level+2)%blocks(i,1)%matrix(ii,jj)=ctemp
                        enddo
                    enddo
                    ! !$omp end parallel do
					
					! ! if(is_nan_mat_c(butterfly_block_randomized(bb)%ButterflyU(i)%matrix,mm,rank))then
						! ! write(*,*)'U'
						! ! stop
					! ! end if
					
					
                enddo
				!$omp end parallel do
            endif
        enddo
                                  

        !$omp parallel do default(shared) private(index_i,mm,ii,jj)
		do index_i=1, num_blocks
			mm=size(butterfly_block_randomized(bb)%ButterflyU(index_i)%matrix,1)
            do ii=1, mm
                do jj=1, num_vectors
                    random2(ii+arr_acc_m(index_i),jj)=b*random2(ii+arr_acc_m(index_i),jj)+a*ButterflyVector(level_butterfly+2)%blocks(index_i,1)%matrix(ii,jj)
					! if(isnan(abs(b*random2(ii+k,jj)+a*ButterflyVector(level_butterfly+2)%blocks(index_i,1)%matrix(ii,jj))))write(*,*)index_i,ii,k,jj,ButterflyVector(level_butterfly+2)%blocks(index_i,1)%matrix(ii,jj),random2(ii+k,jj),a,b
                enddo
            enddo
        enddo    
		!$omp end parallel do		
 
		if(isnan(sum(abs(random2(:,1))**2)))then
			write(*,*)'NAN in 2 butterfly_block_MVP_randomized_dat',butterfly_block_randomized(bb)%level_butterfly
			stop
		end if
        !deallocate (butterflyvector)
                    
    elseif (chara=='T') then
    
        level_butterfly=butterfly_block_randomized(bb)%level_butterfly
        num_blocks=2**level_butterfly
        levelm = ceiling_safe(dble(level_butterfly)/2d0)        
        
        allocate (ButterflyVector(0:level_butterfly+2))
        allocate (ButterflyVector(0)%blocks(num_blocks,1))
        ButterflyVector(0)%num_row=num_blocks
        ButterflyVector(0)%num_col=1
        
		!$omp parallel do default(shared) private(i,mm,ii,jj)
		do i=1, num_blocks
			mm=size(butterfly_block_randomized(bb)%ButterflyU(i)%matrix,1)
            allocate (ButterflyVector(0)%blocks(i,1)%matrix(mm,num_vectors))
            do ii=1, mm
                do jj=1, num_vectors
                    ButterflyVector(0)%blocks(i,1)%matrix(ii,jj)=random1(ii+arr_acc_m(i),jj)
                enddo
            enddo
        enddo 
        !$omp end parallel do		
        
        do level=0, level_butterfly
            if (level==0) then
                num_groupm=num_blocks
                allocate (ButterflyVector(1)%blocks(num_groupm,1))
                ButterflyVector(1)%num_row=num_groupm
                ButterflyVector(1)%num_col=1
                !$omp parallel do default(shared) private(i,rank,mm,ii,jj,ctemp,kk)
				do i=1, num_groupm
                    rank=size(butterfly_block_randomized(bb)%ButterflyU(i)%matrix,2)
                    mm=size(butterfly_block_randomized(bb)%ButterflyU(i)%matrix,1)
                    allocate (ButterflyVector(1)%blocks(i,1)%matrix(rank,num_vectors))
                    ! !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
                    do ii=1, rank
                        do jj=1, num_vectors
                            ctemp=0
                            do kk=1, mm
                                ctemp=ctemp+butterfly_block_randomized(bb)%ButterflyU(i)%matrix(kk,ii)*ButterflyVector(0)%blocks(i,1)%matrix(kk,jj)
                            enddo
                            ButterflyVector(1)%blocks(i,1)%matrix(ii,jj)=ctemp
                        enddo
                    enddo
                    ! !$omp end parallel do
                enddo                   
				!$omp end parallel do				
            else
                num_groupm=butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%num_row
                num_groupn=butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%num_col
                if (num_groupm/=1) then
                    allocate (ButterflyVector(level+1)%blocks(int(num_groupm/2),num_groupn))
                    ButterflyVector(level+1)%num_row=int(num_groupm/2)
                    ButterflyVector(level+1)%num_col=num_groupn                    
                    
					
                    !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,mm1,mm2,nn)
					do ij=1,num_groupn*(num_groupm/2)
						j = (ij-1)/(num_groupm/2)+1
						i = (mod(ij-1,(num_groupm/2)) + 1)*2-1	
						index_j=int((j+1)/2)
						index_i=int((i+1)/2)					

						mm1=size(butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,1)
						mm2=size(butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j)%matrix,1)
						nn=size(butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,2)
						allocate (ButterflyVector(level+1)%blocks(index_i,j)%matrix(nn,num_vectors))
						! !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
						do jj=1, nn
							do ii=1, num_vectors
								ctemp=0
								do kk=1, mm1
									ctemp=ctemp+ButterflyVector(level)%blocks(i,index_j)%matrix(kk,ii)*butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix(kk,jj)
								enddo
								do kk=1, mm2
									ctemp=ctemp+ButterflyVector(level)%blocks(i+1,index_j)%matrix(kk,ii)*butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j)%matrix(kk,jj)
								enddo
								ButterflyVector(level+1)%blocks(index_i,j)%matrix(jj,ii)=ctemp
							enddo
						enddo
						! !$omp end parallel do
						
						if(level_butterfly-level==levelm .and. middleflag==1 .and. level_butterfly>=2)then	
							
							! allocate(matrixtemp(nn,num_vectors))
							! allocate(matrixtemp1(nn,nn))
							! call copymatT_omp(butterfly_block_randomized(bb)%ButterflyMiddle(index_i,j)%matrix,matrixtemp1,nn,nn)
							! call gemm_omp(matrixtemp1,ButterflyVector(level+1)%blocks(index_i,j)%matrix,matrixtemp,nn,nn,num_vectors)
							! ButterflyVector(level+1)%blocks(index_i,j)%matrix = matrixtemp
							! deallocate(matrixtemp)
							! deallocate(matrixtemp1)	
							
							call gemmTN_omp(butterfly_block_randomized(bb)%ButterflyMiddle(index_i,j)%matrix,ButterflyVector(level+1)%blocks(index_i,j)%matrix,ButterflyVector(level+1)%blocks(index_i,j)%matrix,nn,nn,num_vectors)
							
						end if
                    enddo
					!$omp end parallel do
                else
                    allocate (ButterflyVector(level+1)%blocks(1,num_groupn))
                    ButterflyVector(level+1)%num_row=1
                    ButterflyVector(level+1)%num_col=num_groupn
                    do j=1, num_groupn
                        index_j=int((j+1)/2)
                        nn=size(butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(1,j)%matrix,2)
                        mm=size(butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(1,j)%matrix,1)
                        allocate (ButterflyVector(level+1)%blocks(1,j)%matrix(nn,num_vectors))
                        !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
                        do ii=1, nn
                            do jj=1, num_vectors
                                ctemp=0                           
                                do kk=1, mm
                                    ctemp=ctemp+ButterflyVector(level)%blocks(1,index_j)%matrix(kk,jj)*butterfly_block_randomized(bb)%ButterflyKerl(level_butterfly-level+1)%blocks(1,j)%matrix(kk,ii)
                                enddo
                                ButterflyVector(level+1)%blocks(1,j)%matrix(ii,jj)=ctemp
                            enddo
                        enddo
                        !$omp end parallel do
                    enddo
                endif
            endif
            if (level==level_butterfly) then
                allocate (ButterflyVector(level+2)%blocks(1,num_blocks))
                ButterflyVector(level+2)%num_row=1
                ButterflyVector(level+2)%num_col=num_blocks
                !$omp parallel do default(shared) private(j,rank,nn,ii,jj,ctemp,kk)
				do j=1, num_blocks
                    nn=size(butterfly_block_randomized(bb)%ButterflyV(j)%matrix,1)
                    rank=size(butterfly_block_randomized(bb)%ButterflyV(j)%matrix,2)
                    allocate (ButterflyVector(level+2)%blocks(1,j)%matrix(nn,num_vectors))
                    ! !$omp parallel do default(shared) private(ii,jj,kk,ctemp)
                    do ii=1, nn
                        do jj=1, num_vectors
                            ctemp=0
                            do kk=1, rank
                                ctemp=ctemp+ButterflyVector(level+1)%blocks(1,j)%matrix(kk,jj)*butterfly_block_randomized(bb)%ButterflyV(j)%matrix(ii,kk)
                            enddo
                            ButterflyVector(level+2)%blocks(1,j)%matrix(ii,jj)=ctemp
                        enddo
                    enddo
                    ! !$omp end parallel do
                enddo
				!$omp end parallel do
            endif
        enddo
        
		!$omp parallel do default(shared) private(index_j,nn,ii,jj)
		do index_j=1, num_blocks
			nn=size(butterfly_block_randomized(bb)%ButterflyV(index_j)%matrix,1)
            do ii=1, nn
                do jj=1, num_vectors
                    random2(ii+arr_acc_n(index_j),jj)=b*random2(ii+arr_acc_n(index_j),jj)+a*ButterflyVector(level_butterfly+2)%blocks(1,index_j)%matrix(ii,jj)
                enddo
            enddo
        enddo 
		!$omp end parallel do
    
    endif
    
    ! !$omp parallel do default(shared) private(level,i,j)
    do level=0, level_butterfly+2
        do j=1, ButterflyVector(level)%num_col
            do i=1, ButterflyVector(level)%num_row
                deallocate (ButterflyVector(level)%blocks(i,j)%matrix)
            enddo
        enddo
        deallocate (ButterflyVector(level)%blocks)
    enddo
    ! !$omp end parallel do
    deallocate (ButterflyVector)   
    deallocate(arr_acc_m,arr_acc_n)
	
    return
    
end subroutine Butterflys_block_MVP_dat





subroutine Delete_randomized_butterflys(bb)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj
    character chara
    real*8 T0
    type(matrixblock),pointer::block_o
    integer rank_new_max
	real*8:: rank_new_avr,error
	integer niter,bb
	real*8:: error_inout
	integer itermax,levelm,index_i_m,index_j_m
	! real*8:: n1,n2
		

    level_butterfly=butterfly_block_randomized(bb)%level_butterfly
    num_blocks=2**level_butterfly
	
	
	if(allocated(butterfly_block_randomized(bb)%ButterflyMiddle))then
		levelm = ceiling_safe(dble(level_butterfly)/2d0)
		do index_i_m=1, 2**levelm
			do index_j_m=1, 2**(level_butterfly-levelm)	
				deallocate(butterfly_block_randomized(bb)%ButterflyMiddle(index_i_m,index_j_m)%matrix)
			end do
		end do
		deallocate(butterfly_block_randomized(bb)%ButterflyMiddle)
	end if	

	if(allocated(butterfly_block_randomized(bb)%KerInv))deallocate(butterfly_block_randomized(bb)%KerInv)	
		
    do i=1, num_blocks
		deallocate (butterfly_block_randomized(bb)%ButterflyU(i)%matrix)
        deallocate (butterfly_block_randomized(bb)%ButterflyV(i)%matrix)
    enddo
    deallocate (butterfly_block_randomized(bb)%ButterflyU, butterfly_block_randomized(bb)%ButterflyV)
	
	if (level_butterfly/=0) then
        do level=1, level_butterfly
            num_row=butterfly_block_randomized(bb)%ButterflyKerl(level)%num_row
            num_col=butterfly_block_randomized(bb)%ButterflyKerl(level)%num_col
            do j=1, num_col
                do i=1, num_row
                    mm=size(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix,1)
                    nn=size(butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix,2)
                    deallocate (butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks(i,j)%matrix)
                enddo
            enddo
            deallocate (butterfly_block_randomized(bb)%ButterflyKerl(level)%blocks)
        enddo
        deallocate (butterfly_block_randomized(bb)%ButterflyKerl)
    endif
    
    ! deallocate (butterfly_block_randomized)
    
    return

end subroutine Delete_randomized_butterflys





subroutine Get_Randomized_Vectors_LL_HODLR_MVP(level_c,level_butterfly,nth_s,nth_e,num_vect_sub,unique_nth,oddeven)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock,unique_nth,oddeven
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,bb
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    ! complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:)
	
	integer Nsub,Ng
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,level_right_start
	type(RandomBlock), pointer :: random
	real*8::n1,n2
	
    ! ctemp1=1.0d0 ; ctemp2=0.0d0	
	! block_o =>  ho_bf%levels(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 

	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
    ! level_butterfly=int((maxlevel_for_blocks-level_c)/2)*2
    num_blocks=2**level_butterfly
	level_right_start = floor_safe(level_butterfly/2d0)
	Nsub = NINT(2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(level_right_start-unique_nth)))   !  check here later		
    Ng = 2**level_butterfly/Nsub
	! dimension_rank =butterfly_block_randomized(1)%dimension_rank     
	
	
	allocate (RandomVectors_InOutput(3))
	do ii=1,3
		allocate (RandomVectors_InOutput(ii)%vector(Maxedge,num_vect_sub))
		RandomVectors_InOutput(ii)%vector = 0d0
	end do
	
	call assert(ho_bf%levels(level_c)%N_block_forward>1,'N_block_forward not correct')
	call assert(oddeven==0 .or. oddeven==1,'oddeven not correct')
	
	do bb=2-oddeven,ho_bf%levels(level_c)%N_block_forward-oddeven,2
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
	
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
							 RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=random_complex_number()	! matrixtemp1(jj,ii) ! 
						 enddo
					 enddo
					 ! !$omp end parallel do
					 deallocate(matrixtemp1)
					 
				 end if
			end do
		end do		
	
	end do
	
	

	call GetOutputs_BlackBox('T',RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,num_vect_sub)
	call MVM_Z_forward_partial('T',Maxedge,num_vect_sub,level_c-1,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector,ho_bf)
	RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(2)%vector-RandomVectors_InOutput(3)%vector
	
	
	
	do bb=2-oddeven,ho_bf%levels(level_c)%N_block_forward-oddeven,2
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		groupm_start=groupm*2**(level_butterfly)
		random=>random_Block(bb)
		do i=1, num_blocks
			header_m=basis_group(groupm_start+i-1)%head
			k = header_m - 1
			mm=size(butterfly_block_randomized(bb)%butterflyU(i)%matrix,1)
			! !$omp parallel do default(shared) private(ii,jj)
			do ii=1, mm
				do jj=1, num_vect_sub
					random%RandomVectorLL(0)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
				enddo
			enddo
			! !$omp end parallel do
			! k=k+mm
		enddo 
		
		
		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		groupn_start=groupn*2**(level_butterfly)
		do i=1, num_blocks
			header_n=basis_group(groupn_start+i-1)%head
			k = header_n - 1
			nn=size(butterfly_block_randomized(bb)%butterflyV(i)%matrix,1)
			! !$omp parallel do default(shared) private(ii,jj)
			do ii=1, nn
				do jj=1, num_vect_sub
					random%RandomVectorLL(level_butterfly+2)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
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

end subroutine Get_Randomized_Vectors_LL_HODLR_MVP





subroutine Get_Randomized_Vectors_RR_HODLR_MVP(level_c,level_butterfly,nth_s,nth_e,num_vect_sub,unique_nth,oddeven)

    use MODULE_FILE
    ! use lapack95
	use misc
    implicit none
    
	integer level_c,rowblock,unique_nth,oddeven
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test,bb
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    ! complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:),matrixtemp1(:,:)
	
	integer Nsub,Ng
	integer*8 idx_start   
    integer level_blocks
    integer groupm_start, groupn_start
    integer header_mm, header_nn
	integer header_m, header_n, tailer_m, tailer_n
	
	integer nth_s,nth_e,num_vect_sub,nth,num_vect_subsub,level_left_start
	type(RandomBlock), pointer :: random
	real*8::n1,n2
	
    ! ctemp1=1.0d0 ; ctemp2=0.0d0	
	! block_o =>  ho_bf%levels(level_c)%BP(rowblock)%LL(1)%matrices_block(1) 

	num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)	
    ! level_butterfly=int((maxlevel_for_blocks-level_c)/2)*2
    num_blocks=2**level_butterfly
	level_left_start= floor_safe(level_butterfly/2d0)+1
	
	if(mod(level_butterfly,2)==0)then
		Nsub = 2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))    !  check here later
	else 
		Nsub = 2*2**ceiling_safe((level_butterfly-1)/2d0)/dble(2**(unique_nth-level_left_start))
	end if	
	Ng = 2**level_butterfly/Nsub    
	
	allocate (RandomVectors_InOutput(3))
	do ii=1,3
		allocate (RandomVectors_InOutput(ii)%vector(Maxedge,num_vect_sub))
		RandomVectors_InOutput(ii)%vector = 0d0
	end do
	
	call assert(ho_bf%levels(level_c)%N_block_forward>1,'N_block_forward not correct')
	call assert(oddeven==0 .or. oddeven==1,'oddeven not correct')
	
	do bb=2-oddeven,ho_bf%levels(level_c)%N_block_forward-oddeven,2
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
	
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
							 RandomVectors_InOutput(1)%vector(ii+k,(nth-nth_s)*num_vect_subsub+jj)=random_complex_number()	! matrixtemp1(jj,ii) ! 
						 enddo
					 enddo
					 ! !$omp end parallel do
					 deallocate(matrixtemp1)
					 
				 end if
			end do
		end do		
	
	end do
	
	

	call GetOutputs_BlackBox('N',RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(2)%vector,num_vect_sub)
	call MVM_Z_forward_partial('N',Maxedge,num_vect_sub,level_c-1,RandomVectors_InOutput(1)%vector,RandomVectors_InOutput(3)%vector,ho_bf)
	RandomVectors_InOutput(3)%vector = RandomVectors_InOutput(2)%vector-RandomVectors_InOutput(3)%vector
	
	
	
	do bb=2-oddeven,ho_bf%levels(level_c)%N_block_forward-oddeven,2
		block_o =>  ho_bf%levels(level_c)%BP(bb)%LL(1)%matrices_block(1)

		groupn=block_o%col_group  ! Note: row_group and col_group interchanged here   
		nn=basis_group(groupn)%tail-basis_group(groupn)%head+1 
		groupn_start=groupn*2**(level_butterfly)
		random=>random_Block(bb)
		do i=1, num_blocks
			header_n=basis_group(groupn_start+i-1)%head
			k = header_n - 1		
			nn=size(butterfly_block_randomized(bb)%butterflyV(i)%matrix,1)
			! !$omp parallel do default(shared) private(ii,jj)
			do ii=1, nn
				do jj=1, num_vect_sub
					random%RandomVectorRR(0)%blocks(1,i)%matrix(ii,jj)=RandomVectors_InOutput(1)%vector(ii+k,jj)
				enddo
			enddo
			! !$omp end parallel do
			! k=k+nn
		enddo 
		
		
		groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
		mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 
		groupm_start=groupm*2**(level_butterfly)
		do i=1, num_blocks
			header_m=basis_group(groupm_start+i-1)%head
			k = header_m - 1		
			mm=size(butterfly_block_randomized(bb)%butterflyU(i)%matrix,1)
			! !$omp parallel do default(shared) private(ii,jj)
			do ii=1, mm
				do jj=1, num_vect_sub
					random%RandomVectorRR(level_butterfly+2)%blocks(i,1)%matrix(ii,jj)=RandomVectors_InOutput(3)%vector(ii+k,jj)
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

end subroutine Get_Randomized_Vectors_RR_HODLR_MVP






end module HODLR_randomMVP
