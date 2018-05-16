module Bplus_inversion_schur_partition
use Butterfly_inversion_schur_partition
use Randomized_reconstruction

integer rankthusfarBC
contains 


subroutine Bplus_inverse_schur_partitionedinverse(ho_bf1,level_c,rowblock,option,Memory)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock
	type(blockplus),pointer::bplus
	real*8:: n1,n2,Memory,error_inout
	type(Hoption)::option
	type(hobf)::ho_bf1
	
    bplus => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)	
	
	if(bplus%Lplus==1)then
		call OneL_inverse_schur_partitionedinverse(ho_bf1,level_c,rowblock,error_inout,option,Memory)
	else 
		call MultiL_inverse_schur_partitionedinverse(ho_bf1,level_c,rowblock,option,Memory)
	end if	

end subroutine Bplus_inverse_schur_partitionedinverse




subroutine OneL_inverse_schur_partitionedinverse(ho_bf,level_c,rowblock,error_inout,option,Memory)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt,ll
    character chara
    real*8 T0
    type(matrixblock),pointer::block_o,block_off1,block_off2
    integer rank_new_max,rank0
	real*8:: rank_new_avr,error
	integer niter
	real*8:: error_inout,rate,err_avr
	integer itermax,ntry
	real*8:: n1,n2,Memory
	type(Hoption)::option

	error_inout=0
	
	block_off1 => ho_bf1%levels(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)	
	block_off2 => ho_bf1%levels(level_c)%BP(rowblock*2)%LL(1)%matrices_block(1)	

	
	block_o => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)	
	block_o%level_butterfly = block_off1%level_butterfly	
	level_butterfly=block_o%level_butterfly
	
	Memory = 0
	
	
	ho_bf1%ind_lv=level_c
	ho_bf1%ind_bk=rowblock
	
	rank0 = max(block_off1%rankmax,block_off2%rankmax)
	rate=1.2d0
	call Butterfly_randomized(level_butterfly,rank0,rate,block_o,ho_bf1,butterfly_block_MVP_inverse_minusBC_dat,error,'minusBC',option) 	
	error_inout = max(error_inout, error)
	
	n1 = OMP_get_wtime()
	! if(block_o%level==3)then
	if(level_butterfly>=option%schulzlevel)then
		call Butterfly_inverse_schulziteration_IplusButter(block_o,error,option,Memory)
	else
		call Butterfly_inverse_partitionedinverse_IplusButter(block_o,option,error)
	endif

	error_inout = max(error_inout, error)

	n2 = OMP_get_wtime()
	write(*,*)'I+B Inversion Time:',n2-n1	
	
#if PRNTlevel >= 1
	write(*,'(A10,I5,A6,I3,A8,I3,A11,Es14.7)')'OneL No. ',rowblock,' rank:',blocks_minusBC%rankmax,' L_butt:',blocks_minusBC%level_butterfly,' error:',error_inout	
#endif

	
    return

end subroutine OneL_inverse_schur_partitionedinverse




subroutine Butterfly_inverse_schulziteration_IplusButter(block_o,error_inout,option,Memory)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock
    integer groupm,blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks
    integer num_col, num_row, level, mm, nn, ii, jj,tt,ll
    character chara
    real*8 T0
    type(matrixblock)::block_o,block_Xn
    integer rank_new_max,rank0,num_vect
	real*8:: rank_new_avr,error
	integer niter
	real*8:: error_inout,rate,err_avr
	integer itermax,ntry,converged
	real*8:: n1,n2,Memory,memory_temp
	type(Hoption)::option
	type(schulz_operand)::schulz_op
	complex(kind=8),allocatable::VecIn(:,:),VecOut(:,:),VecBuff(:,:)
	complex(kind=8)::ctemp1,ctemp2
	character(len=10)::iternumber 

	error_inout=0
	level_butterfly=block_o%level_butterfly
	
	Memory = 0

	groupm=block_o%row_group  ! Note: row_group and col_group interchanged here   
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	
	num_vect=1
	allocate(VecIn(mm,num_vect))
	allocate(VecOut(mm,num_vect))
	VecOut=0
	allocate(VecBuff(mm,num_vect))
	VecBuff=0
	
	call copy_butterfly(block_o,schulz_op%matrices_block,memory_temp)
	call copy_butterfly(block_o,block_Xn,memory_temp)
	
	call compute_schulz_init(schulz_op,option)
	
	itermax=100
	converged=0
	! n1 = OMP_get_wtime()
	do ii=1,itermax

		write(iternumber ,  "(I4)") ii
		
		rank0 = block_Xn%rankmax
		
		rate=1.2d0
		call Butterfly_randomized(level_butterfly,rank0,rate,block_Xn,schulz_op,butterfly_block_MVP_schulz_dat,error,'schulz iter'//TRIM(iternumber),option,ii) 	
		
		if(schulz_op%order==2)schulz_op%scale=schulz_op%scale*(2-schulz_op%scale)
		if(schulz_op%order==3)schulz_op%scale=schulz_op%scale*(3 - 3*schulz_op%scale + schulz_op%scale**2d0)
		
		! test error
		
		ctemp1=1.0d0 ; ctemp2=0.0d0
		call RandomMat(mm,num_vect,min(mm,num_vect),VecIn,1)
		! XnR
		call butterfly_block_MVP_schulz_Xn_dat(schulz_op,block_Xn,'N',mm,mm,num_vect,VecIn,VecBuff,ctemp1,ctemp2,ii+1)
		
		
		! AXnR
		call butterfly_block_MVP_randomized_dat(schulz_op%matrices_block,'N',mm,mm,num_vect,VecBuff,VecOut,ctemp1,ctemp2)
		VecOut = 	VecBuff+VecOut			
		error_inout = fnorm(VecOut-VecIn,mm,num_vect)/fnorm(VecIn,mm,num_vect)
		

#if PRNTlevel >= 1			
		write(*,'(A22,A6,I3,A8,I2,A8,I3,A7,Es14.7)')' Schultz ',' rank:',block_Xn%rankmax,' Iter:',ii,' L_butt:',block_Xn%level_butterfly,' error:',error_inout
#endif			
		
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
		call delete_blocks(block_o)
		call get_butterfly_minmaxrank(block_Xn)
		rank_new_max = block_Xn%rankmax
		call copy_delete_randomizedbutterfly(block_Xn,block_o,Memory) 
		call delete_blocks(schulz_op%matrices_block)	
		if(allocated(schulz_op%diags))deallocate(schulz_op%diags)
	endif
	
	deallocate(VecIn)
	deallocate(VecOut)
	deallocate(VecBuff)
	
	
    return

end subroutine Butterfly_inverse_schulziteration_IplusButter







subroutine compute_schulz_init(schulz_op,option)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
    integer level_butterfly
    integer mm, nn, mn,ii
    real*8 T0

	real*8:: error
	integer niter,groupm,groupn
	real*8:: error_inout
	integer num_vect,rank,ranktmp,q,qq
	real*8:: n1,n2,memory_temp
	type(Hoption)::option
	type(schulz_operand)::schulz_op
	real*8, allocatable:: Singular(:)
	complex (kind=8), allocatable::UU(:,:),VV(:,:),RandVectIn(:,:),RandVectOut(:,:),matrixtmp(:,:),matrixtmp1(:,:)
	
	schulz_op%order=option%schulzorder
	
	error_inout=0
	
	level_butterfly=schulz_op%matrices_block%level_butterfly
	
	groupm=schulz_op%matrices_block%row_group  ! Note: row_group and col_group interchanged here   
	mm=basis_group(groupm)%tail-basis_group(groupm)%head+1 	
	nn=mm
	num_vect=min(10,nn)

	allocate(RandVectIn(nn,num_vect))
	allocate(RandVectOut(mm,num_vect))
	RandVectOut=0
	call RandomMat(nn,num_vect,min(nn,num_vect),RandVectIn,1)
	
	! computation of AR
	call butterfly_block_MVP_randomized_dat(schulz_op%matrices_block,'N',mm,nn,num_vect,RandVectIn,RandVectOut,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8))
	RandVectOut = RandVectIn+RandVectOut
	
	
	! power iteration of order q, the following is prone to roundoff error, see algorithm 4.4 Halko 2010
	q=6
	do qq=1,q
		RandVectOut=conjg(RandVectOut)
		
		call butterfly_block_MVP_randomized_dat(schulz_op%matrices_block,'T',mm,nn,num_vect,RandVectOut,RandVectIn,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8))
		RandVectIn = RandVectOut+RandVectIn		
		
		RandVectIn=conjg(RandVectIn)
		
		call butterfly_block_MVP_randomized_dat(schulz_op%matrices_block,'N',mm,nn,num_vect,RandVectIn,RandVectOut,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8))
		RandVectOut = RandVectIn+RandVectOut		
		
	enddo	
	

	
	! computation of range Q of AR
	call ComputeRange(mm,num_vect,RandVectOut,ranktmp,0,option%tol_SVD)	
	
	
	! computation of B = Q^c*A
	RandVectOut=conjg(RandVectOut)
	call butterfly_block_MVP_randomized_dat(schulz_op%matrices_block,'T',mm,nn,num_vect,RandVectOut,RandVectIn,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8))
	RandVectIn =RandVectOut+RandVectIn 
	RandVectOut=conjg(RandVectOut)	
	
	! computation of SVD of B and LR of A
	mn=min(nn,ranktmp)
	allocate (UU(nn,mn),VV(mn,ranktmp),Singular(mn))
	call SVD_Truncate(RandVectIn(1:nn,1:ranktmp),nn,ranktmp,mn,UU,VV,Singular,option%tol_SVD,rank)	
	schulz_op%A2norm=Singular(1)
	
	deallocate(UU,VV,Singular)
	
	deallocate(RandVectIn)
	deallocate(RandVectOut)


	
	
	
	! allocate(matrixtmp1(nn,nn))
	! matrixtmp1=0	
	! do ii=1,nn
		! matrixtmp1(ii,ii)=1d0
	! enddo	
	! allocate(matrixtmp(nn,nn))
	! matrixtmp=0
	! call butterfly_block_MVP_randomized_dat(schulz_op%matrices_block,'N',mm,nn,nn,matrixtmp1,matrixtmp,CMPLX(1d0,0d0,kind=8),CMPLX(0d0,0d0,kind=8))
	! matrixtmp = matrixtmp+matrixtmp1
	! allocate (UU(nn,nn),VV(nn,nn),Singular(nn))
	! call SVD_Truncate(matrixtmp,nn,nn,nn,UU,VV,Singular,option%tol_SVD,rank)	
	! write(*,*)Singular(1),schulz_op%A2norm,'nimade'
	! schulz_op%A2norm=Singular(1)
	! deallocate(UU,VV,Singular)	
	! deallocate(matrixtmp)
	! deallocate(matrixtmp1)
	
	
end subroutine compute_schulz_init




subroutine MultiL_inverse_schur_partitionedinverse(ho_bf1,level_c,rowblock,option,Memory)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
    use Butterfly_compress_forward
	
    implicit none
    
	integer level_c,rowblock
    integer blocks1, blocks2, blocks3, level_butterfly, i, j, k, num_blocks,level_butterfly_loc
    integer num_col, num_row, level, mm, nn, ii, jj,tt,ll,llplus,bb,mmb
    character chara
    real*8 T0
    type(matrixblock),pointer::block_o,block_off1,block_off2
    type(matrixblock),pointer::blocks_o_D
    type(matrixblock)::block_tmp
	type(blockplus),pointer::Bplus,Bplus_schur
    integer rank_new_max
	real*8:: rank_new_avr,error,err_avr,err_max
	integer niter
	real*8:: error_inout,rate,rankrate_inner,rankrate_outter
	integer itermax,ntry,cnt,cnt_partial
	real*8:: n1,n2,Memory
	integer rank0,rank0_inner,rank0_outter,Lplus,level_BP,levelm,groupm_start,ij_loc,edge_s,edge_e,edge_first,idx_end_m_ref,idx_start_m_ref,idx_start_b,idx_end_b
	complex(kind=8),allocatable:: matin(:,:),matout(:,:),matin_tmp(:,:),matout_tmp(:,:)
	complex(kind=8):: ctemp1,ctemp2 
	integer, allocatable :: ipiv(:)
	type(Hoption)::option
	type(hobf)::ho_bf1
	
	ctemp1 = 1d0
	ctemp2 = 0d0
	
	block_off1 => ho_bf1%levels(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)			
	block_off2 => ho_bf1%levels(level_c)%BP(rowblock*2)%LL(1)%matrices_block(1)			
	block_o => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
! write(*,*)block_o%row_group,block_o%col_group,level_c,rowblock,block_o%level,'diao'	
	block_o%level_butterfly = block_off1%level_butterfly	
		
	Memory = 0
	
	error_inout=0

	! call MultiL_diagonal_minusBC(level_c,rowblock) 
	
	
	ho_bf1%ind_lv=level_c
	ho_bf1%ind_bk=rowblock
	Bplus =>  ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)
	
	rank0_inner = ho_bf1%levels(level_c)%BP(2*rowblock-1)%LL(2)%rankmax
	rankrate_inner = 2.0d0
	
	rank0_outter = max(block_off1%rankmax,block_off2%rankmax)
	rankrate_outter=1.2d0	
	
	call Buplus_randomized(Bplus,ho_bf1,rank0_inner,rankrate_inner,Bplus_block_MVP_minusBC_dat,rank0_outter,rankrate_outter,Bplus_block_MVP_Outter_minusBC_dat,error,'mBC+',option)
	error_inout = max(error_inout, error)


	
	
	! write(*,*)'good!!!!'
	! stop

	Bplus => ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)
	Lplus = ho_bf1%levels(level_c)%BP_inverse_schur(rowblock)%Lplus
	do llplus =Lplus,1,-1
		do bb=1,Bplus%LL(llplus)%Nbound
			block_o => Bplus%LL(llplus)%matrices_block(bb)
		
			err_max=0
			!!!!! partial update butterflies at level llplus from left B1 = D^-1xB
			if(llplus/=Lplus)then
			
				n1 = OMP_get_wtime()
				call Butterfly_MoveSingulartoLeft(block_o)
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1
				
				level_butterfly = int((Maxlevel_for_blocks - Bplus%LL(llplus)%matrices_block(1)%level)/2)*2 
				level_BP = Bplus%level
				levelm = ceiling_safe(dble(level_butterfly)/2d0)						
				level_butterfly_loc = levelm
				groupm_start=block_o%row_group*2**levelm
				edge_s =basis_group(block_o%row_group)%head 
				edge_e =basis_group(block_o%row_group)%tail 

				
				
				
				do ii=1,Bplus%LL(llplus+1)%Nbound
					edge_first = basis_group(Bplus%LL(llplus+1)%matrices_block(ii)%row_group)%head
					if(edge_first>=edge_s .and. edge_first<=edge_e)then
						ij_loc = Bplus%LL(llplus+1)%matrices_block(ii)%row_group - groupm_start + 1					
						if(level_butterfly_loc==0)then
							write(*,*)'level_butterfly_loc==0 not done'
							stop
						else 

							! allocate(agent_block(1))

							call Extract_partial_butterfly(block_o,level_butterfly_loc,ij_loc,'L')
							call Extract_partial_Bplus(Bplus,llplus+1,Bplus%LL(llplus+1)%matrices_block(ii)%row_group)							
															

															
							rank0 = agent_block(1)%rankmax
							rate=1.2d0
							level_butterfly = agent_block(1)%level_butterfly
							call Butterfly_randomized(level_butterfly,rank0,rate,agent_block(1),agent_bplus(1),Bplus_block_MVP_BplusB_dat,error,'L small',option,agent_block(1)) 	
							call Copy_butterfly_partial(agent_block(1),block_o,level_butterfly_loc,ij_loc,'L',Memory)						
								
							err_max = max(err_max, error)						
								
							call delete_blocks(agent_block(1))
							deallocate(agent_block)	
							call delete_bplus(agent_bplus(1))
							deallocate(agent_bplus)	

						end if
					end if
				end do
				

#if PRNTlevel >= 2 
				write(*,'(A30,I7,A6,I3,A11,Es14.7)')' L partial: ll ',llplus,' bb:',bb,' error:',err_max
#endif	
			end if
			
			error_inout = max(error_inout, err_max)			
   
			! write(*,*)block_o%level_butterfly,'ahaha'
			
			!!!!! invert I+B1 to be I+B2		
				
			call Butterfly_inverse_partitionedinverse_IplusButter(block_o,option,error)		
			error_inout = max(error_inout, error)		
	
			
			err_max=0
			!!!!! partial update butterflies at level llplus from right B2xD^-1
			if(llplus/=Lplus)then
				! write(*,*)'hhe'
				n1 = OMP_get_wtime()	
				call Butterfly_MoveSingulartoRight(block_o)
				n2 = OMP_get_wtime()
				! time_tmp = time_tmp + n2 - n1		
				! write(*,*)'hhe1'
				level_butterfly = int((Maxlevel_for_blocks - Bplus%LL(llplus)%matrices_block(1)%level)/2)*2 
				level_BP = Bplus%level
				levelm = ceiling_safe(dble(level_butterfly)/2d0)						
				level_butterfly_loc = levelm
				groupm_start=block_o%row_group*2**levelm
				edge_s =basis_group(block_o%row_group)%head 
				edge_e =basis_group(block_o%row_group)%tail 

				
				do ii=1,Bplus%LL(llplus+1)%Nbound
					edge_first = basis_group(Bplus%LL(llplus+1)%matrices_block(ii)%row_group)%head
					if(edge_first>=edge_s .and. edge_first<=edge_e)then
						ij_loc = Bplus%LL(llplus+1)%matrices_block(ii)%row_group - groupm_start + 1					
						if(level_butterfly_loc==0)then
							write(*,*)'level_butterfly_loc==0 not done'
							stop
						else 
! write(*,*)'hhe2',level_butterfly_loc,ij_loc,block_o%level_butterfly
							cnt_partial = cnt_partial + 1
							! allocate(agent_block(1))
						   
							call Extract_partial_butterfly(block_o,level_butterfly_loc,ij_loc,'R')
							call Extract_partial_Bplus(Bplus,llplus+1,Bplus%LL(llplus+1)%matrices_block(ii)%row_group)							
							
							
							
							
							rank0 = agent_block(1)%rankmax
							rate=1.2d0
							level_butterfly = agent_block(1)%level_butterfly
							call Butterfly_randomized(level_butterfly,rank0,rate,agent_block(1),agent_bplus(1),Bplus_block_MVP_BBplus_dat,error,'R small',option,agent_block(1)) 	
							call Copy_butterfly_partial(agent_block(1),block_o,level_butterfly_loc,ij_loc,'R',Memory)						
								
							err_max = max(err_max, error)															
								
							call delete_blocks(agent_block(1))
							deallocate(agent_block)	
							call delete_bplus(agent_bplus(1))
							deallocate(agent_bplus)								

						end if
					end if
				end do
				
				! call Butterfly_sym2asym(block_o)
					
				error_inout = max(error_inout, err_max)	
				
#if PRNTlevel >= 2
				write(*,'(A30,I7,A6,I3,A11,Es14.7)')' R partial: ll ',llplus,' bb:',bb,' error:',err_max
#endif				
			end if
		end do
	end do
	
	
					
	call ComputeMemory_Bplus(Bplus,Memory)	
	
	rank_new_max = 0
	do ll=1,Lplus
		rank_new_max = max(rank_new_max,Bplus%LL(ll)%rankmax)
	end do

#if PRNTlevel >= 1	
	write(*,'(A10,I5,A6,I3,A8,I3,A11,Es14.7)')'Mult No. ',rowblock,' rank:',rank_new_max,' L_butt:',Bplus%LL(1)%matrices_block(1)%level_butterfly,' error:',error_inout		
#endif	
	
    return

end subroutine MultiL_inverse_schur_partitionedinverse


subroutine Extract_partial_Bplus(bplus_i,ll_s,row_group)
use MODULE_FILE
use misc
implicit none 
type(matrixblock),pointer::block_i,block_o
type(blockplus)::bplus_i

integer i, j, ii, jj, iii, jjj,index_ij,mm,nn,rank,index_i,index_j,levelm,index_i_m,index_j_m,ll,bb,bb_o
integer level, blocks, edge, patch, node, group,level_c
integer::block_num,block_num_new,num_blocks,level_butterfly,Nboundall	
real*8::rtemp
integer row_group,ll_s,idx_s,idx_e

call assert(bplus_i%row_group==bplus_i%col_group,'only works for square matrix')

idx_s = basis_group(row_group)%head
idx_e = basis_group(row_group)%tail

allocate(agent_bplus(1))
allocate(agent_bplus(1)%LL(LplusMax))
do ll=1,LplusMax
agent_bplus(1)%LL(ll)%Nbound = 0
end do

agent_bplus(1)%Lplus = bplus_i%Lplus - ll_s + 1
agent_bplus(1)%row_group = 	row_group
agent_bplus(1)%col_group = 	row_group
agent_bplus(1)%level = basis_group(row_group)%level

do ll=1,agent_bplus(1)%Lplus
	agent_bplus(1)%LL(ll)%Nbound = 0
	agent_bplus(1)%LL(ll)%rankmax=bplus_i%LL(ll+ll_s-1)%rankmax
	do bb=1,bplus_i%LL(ll+ll_s-1)%Nbound
		if(basis_group(bplus_i%LL(ll+ll_s-1)%matrices_block(bb)%row_group)%head>=idx_s .and. basis_group(bplus_i%LL(ll+ll_s-1)%matrices_block(bb)%row_group)%tail<=idx_e)then
			agent_bplus(1)%LL(ll)%Nbound = agent_bplus(1)%LL(ll)%Nbound + 1
		end if
	end do
	if(agent_bplus(1)%LL(ll)%Nbound>0)then
		allocate(agent_bplus(1)%LL(ll)%matrices_block(agent_bplus(1)%LL(ll)%Nbound))
	end if
end do


do ll=1,agent_bplus(1)%Lplus
	bb_o = 0
	do bb=1,bplus_i%LL(ll+ll_s-1)%Nbound
		if(basis_group(bplus_i%LL(ll+ll_s-1)%matrices_block(bb)%row_group)%head>=idx_s .and. basis_group(bplus_i%LL(ll+ll_s-1)%matrices_block(bb)%row_group)%tail<=idx_e)then
			bb_o = bb_o + 1
			call copy_butterfly(bplus_i%LL(ll+ll_s-1)%matrices_block(bb),agent_bplus(1)%LL(ll)%matrices_block(bb_o))
		end if
	end do
end do



end subroutine Extract_partial_Bplus





end module Bplus_inversion_schur_partition
