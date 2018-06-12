module Butterfly_inversion
use Utilites_randomized
! use Butterfly_compression_givenfullmat
use Bplus_rightmultiply
use Randomized_reconstruction

integer,parameter::dp=kind(0.0d0),sp=kind(0.0)
complex(kind=dp),allocatable::r0_initial(:)

contains 



subroutine MVM_Z_factorized(Ns,num_vectors,Vin,Vout,ho_bf1,ptree,stats)

    use MODULE_FILE
    ! use lapack95
    implicit none
    
	integer Ns
	integer level_c,rowblock,head,tail
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors,pp
    integer mm,nn,mn,blocks1,blocks2,blocks3,level_butterfly,groupm,groupn,groupm_diag
    character chara
    real*8 a,b,c,d
    complex(kind=8) ctemp, ctemp1, ctemp2
	type(matrixblock),pointer::block_o
	
    type(vectorsblock), pointer :: random1, random2
    
    real*8,allocatable :: Singular(:)
	integer idx_start_glo,N_diag,idx_start_diag,idx_start_loc,idx_end_loc
	
	complex(kind=8),allocatable::vec_old(:,:),vec_new(:,:)
	! complex(kind=8)::Vin(:),Vout(:)
	complex(kind=8)::Vin(Ns,num_vectors),Vout(Ns,num_vectors)
	type(hobf)::ho_bf1
    type(proctree)::ptree
	type(Hstat)::stats
		
	idx_start_glo = ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)	 
	   		
	! get the right multiplied vectors
	ctemp1=1.0d0 ; ctemp2=0.0d0
	! allocate(vec_old(Ns,num_vectors))
	allocate(vec_new(Ns,num_vectors))
	Vout = Vin
	! write(*,*)'ddddd',Ns,num_vectors
	! write(*,*)'begin'
	
	do level = ho_bf1%Maxlevel+1,1,-1
		vec_new = 0
		do ii = ho_bf1%levels(level)%Bidxs,ho_bf1%levels(level)%Bidxe
			pp = ptree%MyID - ptree%pgrp(ho_bf1%levels(level)%BP_inverse(ii)%pgno)%head + 1
			head = ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%N_p(pp,1) + ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%headn -1
			tail = head + ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1)%N_loc -1
			idx_start_loc = head-idx_start_glo+1
			idx_end_loc = tail-idx_start_glo+1					
		
			if(level==ho_bf1%Maxlevel+1)then	
				call fullmat_block_MVP_randomized_dat(ho_bf1%levels(level)%BP_inverse(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vectors,&
				&Vout(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2)							
			else 
				call Bplus_block_MVP_inverse_dat(ho_bf1,level,ii,'N',idx_end_loc-idx_start_loc+1,num_vectors,Vout(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ptree,stats)
				
				! call OneL_block_MVP_inverse_dat(ho_bf1,level,ii,'N',idx_end_loc-idx_start_loc+1,num_vectors,Vout(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ptree,stats)
				
			endif
		end do
		Vout = vec_new
	end do
	! Vout = vec_new(1:Ns,1)
	! deallocate(vec_old)
	deallocate(vec_new)	
	 
    return                

end subroutine MVM_Z_factorized
 

end module Butterfly_inversion
