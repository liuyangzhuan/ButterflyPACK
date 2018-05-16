module Butterfly_inversion
use Utilites_randomized
! use Butterfly_compression_givenfullmat
use Bplus_rightmultiply
use Randomized_reconstruction

integer,parameter::dp=kind(0.0d0),sp=kind(0.0)
complex(kind=dp),allocatable::r0_initial(:)

contains 



subroutine MVM_Z_factorized(Ns,num_vectors,Vin,Vout,ho_bf1)

    use MODULE_FILE
    ! use lapack95
    implicit none
    
	integer Ns
	integer level_c,rowblock
    integer i,j,k,level,num_blocks,num_row,num_col,ii,jj,kk,test, num_vectors
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
    
	level_c = 1
	rowblock = 1
	
	idx_start_glo = 1	 
	   

	! if(num_vectors==1)then
		
		
		! num_vectors = 1
		! get the right multiplied vectors
		ctemp1=1.0d0 ; ctemp2=0.0d0
		! allocate(vec_old(Ns,num_vectors))
		allocate(vec_new(Ns,num_vectors))
		Vout = Vin

		! write(*,*)'begin'
		
		do level = Maxlevel_for_blocks+1,1,-1
			N_diag = 2**(level-1)
			idx_start_diag = (rowblock-1)*N_diag+1
			vec_new = 0
			do ii = idx_start_diag,idx_start_diag+N_diag-1

				
				! write(*,*)level,ii,idx_start_loc,idx_end_loc,mm,N_diag
				
				! write(*,*)idx_start_loc,idx_end_loc,idx_start_glo,basis_group(groupm_diag)%head,num_vectors,mm !,block_o%col_group,basis_group(block_o%col_group)%head
				if(level==Maxlevel_for_blocks+1)then
					! write(*,*)level,ii
					groupm_diag = ho_bf1%levels(level)%BP(ii)%row_group ! Note: row_group and col_group interchanged here   
					idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
					idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1				
				
					call fullmat_block_MVP_randomized_dat(ho_bf1%levels(level)%BP(ii)%LL(1)%matrices_block(1),'N',idx_end_loc-idx_start_loc+1,num_vectors,&
					&Vout(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors),ctemp1,ctemp2)
					! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)			
				else 
					! write(*,*)level,ii
					groupm_diag = ho_bf1%levels(level)%BP_inverse_schur(ii)%row_group/2 ! Note: row_group and col_group interchanged here   
					idx_start_loc = basis_group(groupm_diag)%head-idx_start_glo+1
					idx_end_loc = basis_group(groupm_diag)%tail-idx_start_glo+1				
					call Bplus_block_MVP_inverse_dat(ho_bf1,level,ii,'N',idx_end_loc-idx_start_loc+1,num_vectors,Vout(idx_start_loc:idx_end_loc,1:num_vectors),vec_new(idx_start_loc:idx_end_loc,1:num_vectors))
					! ! vec_new(idx_start_loc:idx_end_loc,1:num_vectors) = 	vec_old(idx_start_loc:idx_end_loc,1:num_vectors)	
				end if
			end do
			Vout = vec_new
		end do
		! Vout = vec_new(1:Ns,1)
		! deallocate(vec_old)
		deallocate(vec_new)	
	! else 
		! write(*,*)'multiple RHS not implemented yet'
	! end if	 
    return                

end subroutine MVM_Z_factorized
 

end module Butterfly_inversion
