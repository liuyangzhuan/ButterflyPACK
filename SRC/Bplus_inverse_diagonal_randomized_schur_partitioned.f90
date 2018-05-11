module Bplus_inversion_schur_partition
use Butterfly_inversion_schur_partition
use Randomized_reconstruction

integer rankthusfarBC
contains 


subroutine Bplus_inverse_schur_partitionedinverse(level_c,rowblock,Memory)

    use MODULE_FILE
	use misc
	! use lapack95
	! use blas95
    use omp_lib
	
    implicit none
    
	integer level_c,rowblock
	type(blockplus),pointer::bplus
	real*8:: n1,n2,Memory

    bplus => ho_bf%levels(level_c)%BP_inverse_schur(rowblock)	
	
	if(bplus%Lplus==1)then
		call OneL_inverse_schur_partitionedinverse(level_c,rowblock,Memory)
	else 
		call MultiL_inverse_schur_partitionedinverse(level_c,rowblock,Memory)
	end if	
	
	

end subroutine Bplus_inverse_schur_partitionedinverse




subroutine OneL_inverse_schur_partitionedinverse(level_c,rowblock,Memory)

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
    type(matrixblock),pointer::blocks_minusBC
    integer rank_new_max,rank0
	real*8:: rank_new_avr,error
	integer niter
	real*8:: error_inout,rate
	integer itermax,ntry
	real*8:: n1,n2,Memory


	block_off1 => ho_bf%levels(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)	
	block_off2 => ho_bf%levels(level_c)%BP(rowblock*2)%LL(1)%matrices_block(1)	

	
	block_o => ho_bf%levels(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)	
	block_o%level_butterfly = block_off1%level_butterfly	
	level_butterfly=block_o%level_butterfly
	
	Memory = 0
	
	allocate(partitioned_blocks(0:block_o%level_butterfly+1))
	do ll=0,block_o%level_butterfly+1
		allocate(partitioned_blocks(ll)%blocks_A)
		allocate(partitioned_blocks(ll)%blocks_B)
		allocate(partitioned_blocks(ll)%blocks_C)
		allocate(partitioned_blocks(ll)%blocks_D)
	end do
	blocks_minusBC => partitioned_blocks(0)%blocks_D
	
	blocks_minusBC%col_group=block_o%col_group
	blocks_minusBC%row_group=block_o%row_group
	blocks_minusBC%style=block_o%style
	blocks_minusBC%level = block_o%level
	
	
	
	error_avr_glo = 0
	error_cnt = 0	
	
	ho_bf%ind_lv=level_c
	ho_bf%ind_bk=rowblock
	
	rank0 = max(block_off1%rankmax,block_off2%rankmax)
	rate=1.2d0
	call Butterfly_randomized(level_butterfly,rank0,rate,blocks_minusBC,ho_bf,butterfly_block_MVP_inverse_minusBC_dat,error_inout,'minusBC') 	
	

	call Butterfly_inverse_partitionedinverse_IplusButter(0,1)	

	write(*,'(A10,I5,A6,I3,A8,I3,A11,Es14.7)')'OneL No. ',rowblock,' rank:',blocks_minusBC%rankmax,' L_butt:',blocks_minusBC%level_butterfly,' error_avr:',error_inout	

	call copy_delete_butterfly(blocks_minusBC,block_o,Memory)
	! call delete_blocks(blocks_minusBC)	


	do ll=0,block_o%level_butterfly+1
		deallocate(partitioned_blocks(ll)%blocks_A)
		deallocate(partitioned_blocks(ll)%blocks_B)
		deallocate(partitioned_blocks(ll)%blocks_C)
		deallocate(partitioned_blocks(ll)%blocks_D)
	end do
	deallocate(partitioned_blocks)

	
    return

end subroutine OneL_inverse_schur_partitionedinverse



subroutine MultiL_inverse_schur_partitionedinverse(level_c,rowblock,Memory)

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
	real*8:: rank_new_avr,error,error_avr
	integer niter
	real*8:: error_inout,rate,rankrate_inner,rankrate_outter
	integer itermax,ntry,cnt,cnt_partial
	real*8:: n1,n2,Memory
	integer rank0,rank0_inner,rank0_outter,Lplus,level_BP,levelm,groupm_start,ij_loc,edge_s,edge_e,edge_first,idx_end_m_ref,idx_start_m_ref,idx_start_b,idx_end_b
	complex(kind=8),allocatable:: matin(:,:),matout(:,:),matin_tmp(:,:),matout_tmp(:,:)
	complex(kind=8):: ctemp1,ctemp2 
	integer, allocatable :: ipiv(:)
	
	ctemp1 = 1d0
	ctemp2 = 0d0
	
	block_off1 => ho_bf%levels(level_c)%BP(rowblock*2-1)%LL(1)%matrices_block(1)			
	block_off2 => ho_bf%levels(level_c)%BP(rowblock*2)%LL(1)%matrices_block(1)			
	block_o => ho_bf%levels(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
! write(*,*)block_o%row_group,block_o%col_group,level_c,rowblock,block_o%level,'diao'	
	block_o%level_butterfly = block_off1%level_butterfly	
		
	Memory = 0
	

	! call MultiL_diagonal_minusBC(level_c,rowblock) 
	
	
	ho_bf%ind_lv=level_c
	ho_bf%ind_bk=rowblock
	Bplus =>  ho_bf%levels(level_c)%BP_inverse_schur(rowblock)
	
	rank0_inner = ho_bf%levels(level_c)%BP(2*rowblock-1)%LL(2)%rankmax
	rankrate_inner = 2.0d0
	
	rank0_outter = max(block_off1%rankmax,block_off2%rankmax)
	rankrate_outter=1.2d0	
	
	call Buplus_randomized(Bplus,ho_bf,rank0_inner,rankrate_inner,Bplus_block_MVP_minusBC_dat,rank0_outter,rankrate_outter,Bplus_block_MVP_Outter_minusBC_dat,error_inout,'mBC')
	

	
	! ! Bplus_schur => ho_bf%levels(level_c)%BP_inverse_schur(rowblock)
	! ! idx_start_m_ref = basis_group(Bplus_schur%row_group)%head
	! ! idx_end_m_ref = basis_group(Bplus_schur%row_group)%tail
	! ! mm = idx_end_m_ref - idx_start_m_ref + 1
		! ! ! write(*,*)'e',mm
	! ! allocate(matin(mm,mm))
	! ! allocate(matsub_glo(mm,mm))
		! ! ! write(*,*)'f'
	! ! matin = 0
	! ! do ii=1,mm
		! ! matin(ii,ii) = 1
	! ! end do
	! ! ! write(*,*)'g'
	! ! ! call Bplus_block_MVP_randomized_dat(Bplus_schur,'N',mm,mm,mm,matin,matsub_glo,ctemp1,ctemp2,1,Bplus_schur%Lplus)
	! ! call Bplus_block_MVP_randomized_dat(Bplus_schur,'N',mm,mm,mm,matin,matsub_glo,ctemp1,ctemp2,1,1)
	! ! do ii=1,mm
		! ! matsub_glo(ii,ii) = 1 + matsub_glo(ii,ii)
	! ! end do	
		! ! ! write(*,*)'h'
	
	! ! allocate(ipiv(mm))
	! ! call getrf(matsub_glo,ipiv)
	! ! ! write(*,*)shape(matrixtemp1)
	! ! call getri(matsub_glo,ipiv)	
	! ! deallocate(ipiv)	
	
	! ! do ii=1,mm
		! ! matsub_glo(ii,ii) =  matsub_glo(ii,ii) - 1
	! ! end do		
	
	! ! block_o => ho_bf%levels(level_c)%BP_inverse_schur(rowblock)%LL(1)%matrices_block(1)
	! ! block_tmp%row_group = block_o%row_group
	! ! block_tmp%col_group = block_o%col_group
	! ! block_tmp%level = block_o%level
	! ! block_tmp%style = block_o%style
	
	! ! ! write(*,*)'wocaonidddma',block_tmp%level,block_o%level
	
	! ! call Butterfly_compress_N15_givenfullmat(block_tmp,idx_start_m_ref,idx_start_m_ref)
	
	! ! write(*,*)'wocaonima',block_tmp%rankmax,block_tmp%level_butterfly,block_tmp%level,block_o%level
	! stop
	
	error_avr_glo = 0
	error_cnt = 0
	
	! write(*,*)'good!!!!'
	! stop
	Bplus => ho_bf%levels(level_c)%BP_inverse_schur(rowblock)
	Lplus = ho_bf%levels(level_c)%BP_inverse_schur(rowblock)%Lplus
	do llplus =Lplus,1,-1
		do bb=1,Bplus%LL(llplus)%Nbound
			block_o => Bplus%LL(llplus)%matrices_block(bb)
		
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
				error_avr = 0
				cnt = 0
				cnt_partial = 0
				
				
				
				do ii=1,Bplus%LL(llplus+1)%Nbound
					edge_first = basis_group(Bplus%LL(llplus+1)%matrices_block(ii)%row_group)%head
					if(edge_first>=edge_s .and. edge_first<=edge_e)then
						ij_loc = Bplus%LL(llplus+1)%matrices_block(ii)%row_group - groupm_start + 1					
						if(level_butterfly_loc==0)then
							write(*,*)'level_butterfly_loc==0 not done'
							stop
						else 

							cnt_partial = cnt_partial + 1
							! allocate(agent_block(1))

							call Extract_partial_butterfly(block_o,level_butterfly_loc,ij_loc,'L')
							call Extract_partial_Bplus(Bplus,llplus+1,Bplus%LL(llplus+1)%matrices_block(ii)%row_group)							
															

															
							rank0 = agent_block(1)%rankmax
							rate=1.2d0
							level_butterfly = agent_block(1)%level_butterfly
							call Butterfly_randomized(level_butterfly,rank0,rate,agent_block(1),agent_bplus(1),Bplus_block_MVP_BplusB_dat,error_inout,'L small',agent_block(1)) 	
							call Copy_butterfly_partial(agent_block(1),block_o,level_butterfly_loc,ij_loc,'L',Memory)						
								
							error_avr = error_avr + error_inout
							cnt = cnt + 1
							error_avr_glo = error_avr_glo + error_inout
							error_cnt = error_cnt + 1															
								
							call delete_blocks(agent_block(1))
							deallocate(agent_block)	
							call delete_bplus(agent_bplus(1))
							deallocate(agent_bplus)	

						end if
					end if
				end do
				

				
				if(cnt>0)error_avr = error_avr/cnt
				if(verboselevel>=1)write(*,'(A30,I7,A6,I3,A11,Es14.7)')' L partial: ll ',llplus,' bb:',bb,' error_avr:',error_avr
			end if
			
			
   
			! write(*,*)block_o%level_butterfly,'ahaha'
			
			!!!!! invert I+B1 to be I+B2		
			allocate(partitioned_blocks(0:block_o%level_butterfly+1))
			do ll=0,block_o%level_butterfly+1
				allocate(partitioned_blocks(ll)%blocks_A)
				allocate(partitioned_blocks(ll)%blocks_B)
				allocate(partitioned_blocks(ll)%blocks_C)
				allocate(partitioned_blocks(ll)%blocks_D)
			end do
			blocks_o_D => partitioned_blocks(0)%blocks_D
			call copy_delete_butterfly(block_o,blocks_o_D)
			! call delete_blocks(block_o)
			
! write(*,*)blocks_o_D%level_butterfly,blocks_o_D%rankmax,'heihei'
! stop			
			call Butterfly_inverse_partitionedinverse_IplusButter(0,1)		
							   
			call copy_delete_butterfly(blocks_o_D,block_o,Memory)
			! call delete_blocks(blocks_o_D)	
			do ll=0,block_o%level_butterfly+1
				deallocate(partitioned_blocks(ll)%blocks_A)
				deallocate(partitioned_blocks(ll)%blocks_B)
				deallocate(partitioned_blocks(ll)%blocks_C)
				deallocate(partitioned_blocks(ll)%blocks_D)
			end do
			deallocate(partitioned_blocks)			
			
			
		
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
				error_avr = 0
				cnt = 0
				cnt_partial = 0
				
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
							call Butterfly_randomized(level_butterfly,rank0,rate,agent_block(1),agent_bplus(1),Bplus_block_MVP_BBplus_dat,error_inout,'R small',agent_block(1)) 	
							call Copy_butterfly_partial(agent_block(1),block_o,level_butterfly_loc,ij_loc,'R',Memory)						
								
							error_avr = error_avr + error_inout
							cnt = cnt + 1
							error_avr_glo = error_avr_glo + error_inout
							error_cnt = error_cnt + 1															
								
							call delete_blocks(agent_block(1))
							deallocate(agent_block)	
							call delete_bplus(agent_bplus(1))
							deallocate(agent_bplus)								

						end if
					end if
				end do
				
				! call Butterfly_sym2asym(block_o)
				
				
				if(cnt>0)error_avr = error_avr/cnt
				if(verboselevel>=1)write(*,'(A30,I7,A6,I3,A11,Es14.7)')' R partial: ll ',llplus,' bb:',bb,' error_avr:',error_avr
			end if
		end do
	end do
	
	call ComputeMemory_Bplus(Bplus,Memory)	
	
	rank_new_max = 0
	do ll=1,Lplus
		rank_new_max = max(rank_new_max,Bplus%LL(ll)%rankmax)
	end do
	
	write(*,'(A10,I5,A6,I3,A8,I3,A11,Es14.7)')'Mult No. ',rowblock,' rank:',rank_new_max,' L_butt:',Bplus%LL(1)%matrices_block(1)%level_butterfly,' error_avr:',error_avr_glo/error_cnt		
	
	
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
