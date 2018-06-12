module Utilites_randomized
use Utilities
use omp_lib 
contains


subroutine Butterfly_Partial_MVP_Half(block_rand,chara,level_start,level_end,random,num_vect_sub,nth_s,nth_e,Ng)
    
    use MODULE_FILE
    implicit none
    
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, ij, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2, level_end, level_start
    complex(kind=8) ctemp, a, b
    character chara
	integer num_vect_sub,num_vect_subsub,nth_s,nth_e,Ng,nth,dimension_rank,level_butterfly
    
    type(RandomBlock) :: random
    
    type(butterfly_Kerl),allocatable :: ButterflyVector(:)
    type(matrixblock)::block_rand
   ! write(*,*)'in '
   
	level_butterfly=block_rand%level_butterfly
	dimension_rank = block_rand%dimension_rank
    num_vect_subsub = num_vect_sub/(nth_e-nth_s+1)
   
    if (chara=='N') then

        num_blocks=2**level_butterfly
        
        do level=level_start, level_end
            if (level==0) then
                num_groupn=num_blocks
                do nth=nth_s, nth_e
					! !$omp parallel do default(shared) private(j,rank,nn,ii,jj,ctemp,kk)
					do j = (nth-1)*Ng+1,nth*Ng
						rank=size(block_rand%ButterflyV%blocks(j)%matrix,2)
						nn=size(block_rand%ButterflyV%blocks(j)%matrix,1)
						if(.not. allocated(random%RandomVectorRR(1)%blocks(1,j)%matrix))allocate(random%RandomVectorRR(1)%blocks(1,j)%matrix(rank,num_vect_subsub))
						random%RandomVectorRR(1)%blocks(1,j)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
						!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
						do jj=1, num_vect_subsub						
							do ii=1, rank
								ctemp=(0.,0.)
								do kk=1, nn
									ctemp=ctemp+block_rand%ButterflyV%blocks(j)%matrix(kk,ii)*random%RandomVectorRR(0)%blocks(1,j)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
								enddo
								random%RandomVectorRR(1)%blocks(1,j)%matrix(ii,jj+(nth-nth_s)*num_vect_subsub)=ctemp
							enddo
						enddo
						!$omp end parallel do
					enddo
					! !$omp end parallel do					
                enddo
            elseif (level==level_butterfly+1) then
                                 
            else
                num_groupm=block_rand%ButterflyKerl(level)%num_row
                num_groupn=block_rand%ButterflyKerl(level)%num_col
					
				! !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,nn1,nn2,mm,nth)
				do ij=1,(num_groupm/2)*(num_groupn/2)
					index_i = (ij-1)/(num_groupn/2)+1
					index_j = mod(ij-1,(num_groupn/2)) + 1
					i = index_i*2-1
					j = index_j*2-1
					
																					  
																						
																					 
					do nth = nth_s,nth_e
																																
						if((j>=(nth-1)*Ng/2**(level-1)+1 .and. j<=nth*Ng/2**(level-1)) .or. &
						& (j+1>=(nth-1)*Ng/2**(level-1)+1 .and. j+1<=nth*Ng/2**(level-1)))then						
							nn1=size(block_rand%ButterflyKerl(level)%blocks(i,j)%matrix,2)
							nn2=size(block_rand%ButterflyKerl(level)%blocks(i,j+1)%matrix,2)
							mm=size(block_rand%ButterflyKerl(level)%blocks(i,j)%matrix,1)
							! write(*,*)ij,i,j,level,'ha',index_i
							if(.not. allocated(random%RandomVectorRR(level+1)%blocks(i,index_j)%matrix))allocate(random%RandomVectorRR(level+1)%blocks(i,index_j)%matrix(mm,num_vect_subsub))
							random%RandomVectorRR(level+1)%blocks(i,index_j)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
							!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
							do jj=1, num_vect_subsub							
								do ii=1, mm
									ctemp=(0.,0.)
									do kk=1, nn1
										ctemp=ctemp+block_rand%ButterflyKerl(level)%blocks(i,j)%matrix(ii,kk)*random%RandomVectorRR(level)%blocks(index_i,j)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
									enddo
									do kk=1, nn2
										ctemp=ctemp+block_rand%ButterflyKerl(level)%blocks(i,j+1)%matrix(ii,kk)*random%RandomVectorRR(level)%blocks(index_i,j+1)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
									enddo
									random%RandomVectorRR(level+1)%blocks(i,index_j)%matrix(ii,jj+(nth-nth_s)*num_vect_subsub)=ctemp
								enddo
							enddo
							!$omp end parallel do
							
							nn1=size(block_rand%ButterflyKerl(level)%blocks(i+1,j)%matrix,2)
							nn2=size(block_rand%ButterflyKerl(level)%blocks(i+1,j+1)%matrix,2)
							mm=size(block_rand%ButterflyKerl(level)%blocks(i+1,j)%matrix,1)
							! write(*,*)ij,i,j,level,'ha',index_i
							if(.not. allocated(random%RandomVectorRR(level+1)%blocks(i+1,index_j)%matrix))allocate(random%RandomVectorRR(level+1)%blocks(i+1,index_j)%matrix(mm,num_vect_subsub))
							random%RandomVectorRR(level+1)%blocks(i+1,index_j)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
							!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
							do jj=1, num_vect_subsub							
								do ii=1, mm
									ctemp=(0.,0.)
									do kk=1, nn1
										ctemp=ctemp+block_rand%ButterflyKerl(level)%blocks(i+1,j)%matrix(ii,kk)*random%RandomVectorRR(level)%blocks(index_i,j)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
									enddo
									do kk=1, nn2
										ctemp=ctemp+block_rand%ButterflyKerl(level)%blocks(i+1,j+1)%matrix(ii,kk)*random%RandomVectorRR(level)%blocks(index_i,j+1)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
									enddo
									random%RandomVectorRR(level+1)%blocks(i+1,index_j)%matrix(ii,jj+(nth-nth_s)*num_vect_subsub)=ctemp
								enddo
							enddo
							!$omp end parallel do							
							
							! write(*,*)ij,i,j,level,'ha done0',index_i
							deallocate(random%RandomVectorRR(level)%blocks(index_i,j)%matrix)
							deallocate(random%RandomVectorRR(level)%blocks(index_i,j+1)%matrix)
							! write(*,*)ij,i,j,level,'ha done',index_i
						end if	
					end do	
				enddo
				! !$omp end parallel do
            endif
        enddo      
        
                    
    elseif (chara=='T') then
    
        num_blocks=2**level_butterfly

        do level=level_start, level_end
            if (level==0) then
                num_groupm=num_blocks
                do nth=nth_s, nth_e
					! !$omp parallel do default(shared) private(i,rank,mm,ii,jj,ctemp,kk)
					do i = (nth-1)*Ng+1,nth*Ng
						rank=size(block_rand%ButterflyU%blocks(i)%matrix,2)
						mm=size(block_rand%ButterflyU%blocks(i)%matrix,1)
						if(.not. allocated(random%RandomVectorLL(1)%blocks(i,1)%matrix))allocate(random%RandomVectorLL(1)%blocks(i,1)%matrix(rank,num_vect_subsub))
						random%RandomVectorLL(1)%blocks(i,1)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
						!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
						do jj=1, num_vect_subsub						
							do ii=1, rank
								ctemp=(0.,0.)
								do kk=1, mm
									ctemp=ctemp+block_rand%ButterflyU%blocks(i)%matrix(kk,ii)*random%RandomVectorLL(0)%blocks(i,1)%matrix(kk,jj+(nth-nth_s)*num_vect_subsub)
								enddo
								random%RandomVectorLL(1)%blocks(i,1)%matrix(ii,jj+(nth-nth_s)*num_vect_subsub)=ctemp
							enddo
						enddo
						!$omp end parallel do
					end do
					! !$omp end parallel do
                enddo
            elseif (level==level_butterfly+1) then               
            else
                num_groupm=block_rand%ButterflyKerl(level_butterfly-level+1)%num_row
                num_groupn=block_rand%ButterflyKerl(level_butterfly-level+1)%num_col             

				! !$omp parallel do default(shared) private(ij,ii,jj,kk,ctemp,i,j,index_i,index_j,mm1,mm2,nn,nth)
				do ij=1,(num_groupn/2)*(num_groupm/2)
					index_j = (ij-1)/(num_groupm/2)+1
					index_i = mod(ij-1,(num_groupm/2)) + 1	
					j = 2*index_j-1
					i = 2*index_i-1						
					
																										
																										  
																									   
					do nth = nth_s,nth_e
																																
						if((i>=(nth-1)*Ng/2**(level-1)+1 .and. i<=nth*Ng/2**(level-1)) .or. &
						& (i+1>=(nth-1)*Ng/2**(level-1)+1 .and. i+1<=nth*Ng/2**(level-1)))then
							mm1=size(block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,1)
							mm2=size(block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j)%matrix,1)
							nn=size(block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix,2)
							if(.not. allocated(random%RandomVectorLL(level+1)%blocks(index_i,j)%matrix))allocate(random%RandomVectorLL(level+1)%blocks(index_i,j)%matrix(nn,num_vect_subsub))
							random%RandomVectorLL(level+1)%blocks(index_i,j)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
							!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
							do ii=1, num_vect_subsub							
								do jj=1, nn
									ctemp=(0.,0.)
									do kk=1, mm1
										ctemp=ctemp+random%RandomVectorLL(level)%blocks(i,index_j)%matrix(kk,ii+(nth-nth_s)*num_vect_subsub)*block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i,j)%matrix(kk,jj)
									enddo
									do kk=1, mm2
										ctemp=ctemp+random%RandomVectorLL(level)%blocks(i+1,index_j)%matrix(kk,ii+(nth-nth_s)*num_vect_subsub)*block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j)%matrix(kk,jj)
									enddo
									random%RandomVectorLL(level+1)%blocks(index_i,j)%matrix(jj,ii+(nth-nth_s)*num_vect_subsub)=ctemp
								enddo
							enddo
							!$omp end parallel do
							mm1=size(block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i,j+1)%matrix,1)
							mm2=size(block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j+1)%matrix,1)
							nn=size(block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i,j+1)%matrix,2)
							if(.not. allocated(random%RandomVectorLL(level+1)%blocks(index_i,j+1)%matrix))allocate(random%RandomVectorLL(level+1)%blocks(index_i,j+1)%matrix(nn,num_vect_subsub))
							random%RandomVectorLL(level+1)%blocks(index_i,j+1)%matrix(:,(nth-nth_s)*num_vect_subsub+1:(nth-nth_s+1)*num_vect_subsub) = 0
							!$omp parallel do default(shared) private(ii,jj,kk,ctemp)
							do ii=1, num_vect_subsub							
								do jj=1, nn
									ctemp=(0.,0.)
									do kk=1, mm1
										ctemp=ctemp+random%RandomVectorLL(level)%blocks(i,index_j)%matrix(kk,ii+(nth-nth_s)*num_vect_subsub)*block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i,j+1)%matrix(kk,jj)
									enddo
									do kk=1, mm2
										ctemp=ctemp+random%RandomVectorLL(level)%blocks(i+1,index_j)%matrix(kk,ii+(nth-nth_s)*num_vect_subsub)*block_rand%ButterflyKerl(level_butterfly-level+1)%blocks(i+1,j+1)%matrix(kk,jj)
									enddo
									random%RandomVectorLL(level+1)%blocks(index_i,j+1)%matrix(jj,ii+(nth-nth_s)*num_vect_subsub)=ctemp
								enddo
							enddo
							!$omp end parallel do							
							deallocate(random%RandomVectorLL(level)%blocks(i,index_j)%matrix)
							deallocate(random%RandomVectorLL(level)%blocks(i+1,index_j)%matrix)
						end if
					end do
				enddo
				! !$omp end parallel do

            endif
        enddo
    
    endif
       ! write(*,*)'out '
    return
    
end subroutine Butterfly_Partial_MVP_Half

    



subroutine condition_internalmatrices(random,chara,level,dimension_rank,level_butterfly,num_vectors,condition_number)

    use MODULE_FILE
    ! use lapack95
    ! use blas95
    implicit none
    
    integer i,j,i_old,j_old,iijj,k,level,num_blocks,num_row,num_col,ii,jj,kk,level_left,level_right, rank
    integer index_i, index_j, mm, nn, iter, vector1, vector2, direction, round, mi, nj, min_mn
    real*8 a,b,c,d,norm1,norm2,norm3,norm4,error,condition_number
    complex(kind=8) ctemp
    character chara
	integer dimension_rank,level_butterfly,num_vectors
    
    type(RandomBlock) :: random
    
    complex(kind=8), allocatable :: matrixtemp1(:,:), UU(:,:), VV(:,:)
    real*8, allocatable :: Singular (:)
    
    num_blocks=2**level_butterfly
    condition_number = 0
	
    if (chara=='R') then
        
        ! if (level==0) then
            ! write(*,*)'level cannot be zero'
			! stop
        ! elseif (level==level_butterfly+2) then
            ! write(*,*)'level cannot be level_butterfly+2'
			! stop
        ! else
            ! rank=dimension_rank
        ! endif
        ! allocate (matrixtemp1(num_blocks*rank,num_vectors))
        ! min_mn=min(num_blocks*rank,num_vectors)
        ! allocate (Singular(min_mn))
        ! num_row=random%RandomVectorRR(level)%num_row
        ! num_col=random%RandomVectorRR(level)%num_col
        ! do i=1, num_row
            ! do j=1, num_col
                ! mm=((i-1)*num_col+j-1)*rank
                ! !$omp parallel do default(shared) private(ii,jj)
                ! do jj=1, num_vectors
                    ! do ii=1, rank
                        ! matrixtemp1(mm+ii,jj)=random%RandomVectorRR(level)%blocks(i,j)%matrix(ii,jj)
                    ! enddo
                ! enddo
                ! !$omp end parallel do
            ! enddo
        ! enddo
        ! call gesvd(matrixtemp1,Singular)
        ! condition_number=Singular(1)/Singular(min_mn)
        ! write(*,*)'haddh',Singular
        ! !open (30,file='Singular_distribution_R.txt')
        ! !do i=1, min_mn
        ! !    write (30,*) i, Singular(i)
        ! !enddo
        ! !close (30)
        ! !pause
        
        ! deallocate (matrixtemp1,Singular)
		
		
		
		
		
        if (level==0) then
            write(*,*)'level cannot be zero'
			stop
        elseif (level==level_butterfly+2) then
            write(*,*)'level cannot be level_butterfly+2'
			stop
        else
            rank=dimension_rank
        endif
        num_row=random%RandomVectorRR(level)%num_row
        num_col=random%RandomVectorRR(level)%num_col

        allocate (matrixtemp1(2*rank,num_vectors))
        min_mn=min(2*rank,num_vectors)
        allocate (Singular(min_mn))
		allocate (UU(2*rank,min_mn))
		allocate (VV(min_mn,num_vectors))
		
		
		iijj = 0
        do i=1, num_row
            do j=1, num_col
				iijj  = iijj + 1
			    if(mod(iijj,2)==0)then
					!$omp parallel do default(shared) private(ii,jj)
					do jj=1, num_vectors
						do ii=1, rank
							matrixtemp1(ii,jj)=random%RandomVectorRR(level)%blocks(i_old,j_old)%matrix(ii,jj)
							matrixtemp1(ii+rank,jj)=random%RandomVectorRR(level)%blocks(i,j)%matrix(ii,jj)
						enddo
					enddo
					!$omp end parallel do
					call gesvd_robust(matrixtemp1,Singular,UU,VV,2*rank,num_vectors,min_mn)
					condition_number=max(condition_number,Singular(1)/Singular(min_mn))
					! write(*,*)'hah',Singular
					! write(*,*)condition_number
				end if
				i_old = i
				j_old = j
			enddo
        enddo
        deallocate (matrixtemp1,Singular,UU,VV)	
		

		
        ! if (level==0) then
            ! write(*,*)'level cannot be zero'
			! stop
        ! elseif (level==level_butterfly+2) then
            ! write(*,*)'level cannot be level_butterfly+2'
			! stop
        ! else
            ! rank=dimension_rank
        ! endif
        ! num_row=random%RandomVectorRR(level)%num_row
        ! num_col=random%RandomVectorRR(level)%num_col

        ! allocate (matrixtemp1(rank,num_vectors))
        ! min_mn=min(rank,num_vectors)
        ! allocate (Singular(min_mn))
		
        ! do i=1, num_row
            ! do j=1, num_col
				! !$omp parallel do default(shared) private(ii,jj)
				! do jj=1, num_vectors
					! do ii=1, rank
						! matrixtemp1(ii,jj)=random%RandomVectorRR(level)%blocks(i,j)%matrix(ii,jj)
					! enddo
				! enddo
				! !$omp end parallel do
				! call gesvd(matrixtemp1,Singular)
				! condition_number=max(condition_number,Singular(1)/Singular(min_mn))
				! ! write(*,*)condition_number
				! write(*,*)Singular(1),Singular(min_mn)
			! enddo
        ! enddo
        ! deallocate (matrixtemp1,Singular)			
		
		
		
		
    elseif (chara=='L') then
    
        ! if (level==0) then
            ! write(*,*)'level cannot be zero'
			! stop
        ! elseif (level==level_butterfly+2) then
            ! write(*,*)'level cannot be level_butterfly+2'
			! stop
        ! else
            ! rank=dimension_rank
        ! endif
        ! allocate (matrixtemp1(num_blocks*rank,num_vectors))
        ! min_mn=min(num_blocks*rank,num_vectors)
        ! allocate (Singular(min_mn))
        ! num_row=random%RandomVectorLL(level)%num_row
        ! num_col=random%RandomVectorLL(level)%num_col
        ! do i=1, num_row
            ! do j=1, num_col
                ! mm=((i-1)*num_col+j-1)*rank
                ! !$omp parallel do default(shared) private(ii,jj)
                ! do jj=1, num_vectors
                    ! do ii=1, rank
                        ! matrixtemp1(ii+mm,jj)=random%RandomVectorLL(level)%blocks(i,j)%matrix(ii,jj)
                    ! enddo
                ! enddo
                ! !$omp end parallel do
            ! enddo
        ! enddo
        ! call gesvd(matrixtemp1,Singular)
		! ! write(*,*)min_mn,num_blocks*rank,num_vectors
        ! condition_number=Singular(1)/Singular(min_mn)
        
        ! deallocate (matrixtemp1,Singular)
        
		
		
		
		
        if (level==0) then
            write(*,*)'level cannot be zero'
			stop
        elseif (level==level_butterfly+2) then
            write(*,*)'level cannot be level_butterfly+2'
			stop
        else
            rank=dimension_rank
        endif
        allocate (matrixtemp1(2*rank,num_vectors))
        min_mn=min(2*rank,num_vectors)
        allocate (Singular(min_mn))
		allocate (UU(2*rank,min_mn))
		allocate (VV(min_mn,num_vectors))
        num_row=random%RandomVectorLL(level)%num_row
        num_col=random%RandomVectorLL(level)%num_col
        iijj = 0
		do j=1, num_col
			do i=1, num_row
                iijj = iijj +1
                if(mod(iijj,2)==0)then
					!$omp parallel do default(shared) private(ii,jj)
					do jj=1, num_vectors
						do ii=1, rank
							matrixtemp1(ii,jj)=random%RandomVectorLL(level)%blocks(i,j)%matrix(ii,jj)
							matrixtemp1(ii+rank,jj)=random%RandomVectorLL(level)%blocks(i_old,j_old)%matrix(ii,jj)
						enddo
					enddo
					!$omp end parallel do
					call gesvd_robust(matrixtemp1,Singular,UU,VV,2*rank,num_vectors,min_mn)
					! write(*,*)min_mn,num_blocks*rank,num_vectors
					condition_number=max(condition_number,Singular(1)/Singular(min_mn))	
				end if	
				i_old = i
				j_old = j
            enddo
        enddo
      
        
        deallocate (matrixtemp1,Singular,UU,VV)		
		
		
    endif
    
    return

end subroutine condition_internalmatrices


subroutine Delete_RandVect(chara,random,level_butterfly)
	use MODULE_FILE
	implicit none 
    character chara
	integer num_col,num_row,level,i,j	
	type(RandomBlock):: random
	integer level_butterfly

	if(chara=='T')then
		do level=0, level_butterfly+2
			num_row=random%RandomVectorLL(level)%num_row
			num_col=random%RandomVectorLL(level)%num_col
			do j=1, num_col
				do i=1, num_row
					if(allocated(random%RandomVectorLL(level)%blocks(i,j)%matrix))deallocate (random%RandomVectorLL(level)%blocks(i,j)%matrix)
					if(level/=0 .and. level/=level_butterfly+2)then
						if(allocated(random%RandomVectorLL(level)%blocks(i,j)%null))deallocate (random%RandomVectorLL(level)%blocks(i,j)%null)
					end if
				enddo
			enddo
			deallocate (random%RandomVectorLL(level)%blocks)
		enddo
		deallocate (random%RandomVectorLL)
	else if(chara=='N')then
		do level=0, level_butterfly+2
			num_row=random%RandomVectorRR(level)%num_row
			num_col=random%RandomVectorRR(level)%num_col
			do j=1, num_col
				do i=1, num_row
					if(allocated(random%RandomVectorRR(level)%blocks(i,j)%matrix))deallocate (random%RandomVectorRR(level)%blocks(i,j)%matrix)
					if(level/=0 .and. level/=level_butterfly+2)then
						if(allocated(random%RandomVectorRR(level)%blocks(i,j)%null))deallocate (random%RandomVectorRR(level)%blocks(i,j)%null)
					end if
				enddo
			enddo
			deallocate (random%RandomVectorRR(level)%blocks)
		enddo
		deallocate (random%RandomVectorRR)		
	end if
end subroutine Delete_RandVect


subroutine Init_RandVect_Empty(chara,random,num_vect_sub,block_rand,stats)
    
    use MODULE_FILE
    implicit none
        real*8:: mem_vec
    integer n, group_m, group_n, group_mm, group_nn, index_i, index_j, na, nb, index_start
    integer i, j, ii, jj, level, groupm_start, groupn_start, index_iijj, index_ij, k, kk, intemp1, intemp2
    integer header_m, header_n, tailer_m, tailer_n, mm, nn, num_blocks, level_define, col_vector
    integer rank1, rank2, rank, num_groupm, num_groupn, header_nn, header_mm, ma, mb
    integer vector_a, vector_b, nn1, nn2, level_blocks, mm1, mm2,num_vect_sub,level_butterfly
    complex(kind=8) ctemp, a, b
    character chara
	integer num_col,num_row
	type(matrixblock):: block_rand
    type(Hstat)::stats
    ! type(matricesblock), pointer :: blocks
    type(RandomBlock) :: random
    
    level_butterfly=block_rand%level_butterfly
	
    if (chara=='N') then

        num_blocks=2**level_butterfly
            
        allocate (random%RandomVectorRR(0)%blocks(1,num_blocks))
        random%RandomVectorRR(0)%num_row=1
        random%RandomVectorRR(0)%num_col=num_blocks
        do i=1, num_blocks
            nn=size(block_rand%ButterflyV%blocks(i)%matrix,1)
            allocate (random%RandomVectorRR(0)%blocks(1,i)%matrix(nn,num_vect_sub))
			random%RandomVectorRR(0)%blocks(1,i)%matrix = 0
        enddo 
        do level=0, level_butterfly
            if (level==0) then
                num_groupn=num_blocks
                allocate (random%RandomVectorRR(1)%blocks(1,num_groupn))
                random%RandomVectorRR(1)%num_row=1
                random%RandomVectorRR(1)%num_col=num_groupn                   
            else
                num_groupm=block_rand%ButterflyKerl(level)%num_row
                num_groupn=block_rand%ButterflyKerl(level)%num_col               
				allocate (random%RandomVectorRR(level+1)%blocks(num_groupm,int(num_groupn/2)))
				random%RandomVectorRR(level+1)%num_row=num_groupm
				random%RandomVectorRR(level+1)%num_col=int(num_groupn/2)                    
				           
            endif
            if (level==level_butterfly) then
                allocate (random%RandomVectorRR(level+2)%blocks(num_blocks,1))
                random%RandomVectorRR(level+2)%num_row=num_blocks
                random%RandomVectorRR(level+2)%num_col=1
                do i=1, num_blocks
                    rank=size(block_rand%ButterflyU%blocks(i)%matrix,2)
                    mm=size(block_rand%ButterflyU%blocks(i)%matrix,1)
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
		stats%Mem_int_vec = max(stats%Mem_int_vec,mem_vec) 
        
    elseif (chara=='T') then
    
        num_blocks=2**level_butterfly
        allocate (random%RandomVectorLL(0)%blocks(num_blocks,1))
        random%RandomVectorLL(0)%num_row=num_blocks
        random%RandomVectorLL(0)%num_col=1
        do i=1, num_blocks
			mm=size(block_rand%ButterflyU%blocks(i)%matrix,1)			
            allocate (random%RandomVectorLL(0)%blocks(i,1)%matrix(mm,num_vect_sub))
			random%RandomVectorLL(0)%blocks(i,1)%matrix = 0
        enddo         
        do level=0, level_butterfly
            if (level==0) then
                num_groupm=num_blocks
                allocate (random%RandomVectorLL(1)%blocks(num_groupm,1))
                random%RandomVectorLL(1)%num_row=num_groupm
                random%RandomVectorLL(1)%num_col=1                
            else
                num_groupm=block_rand%ButterflyKerl(level_butterfly-level+1)%num_row
                num_groupn=block_rand%ButterflyKerl(level_butterfly-level+1)%num_col
                
				allocate (random%RandomVectorLL(level+1)%blocks(int(num_groupm/2),num_groupn))
				random%RandomVectorLL(level+1)%num_row=int(num_groupm/2)
				random%RandomVectorLL(level+1)%num_col=num_groupn                    
            endif
            if (level==level_butterfly) then
                allocate (random%RandomVectorLL(level+2)%blocks(1,num_blocks))
                random%RandomVectorLL(level+2)%num_row=1
                random%RandomVectorLL(level+2)%num_col=num_blocks
                do j=1, num_blocks
                    nn=size(block_rand%ButterflyV%blocks(j)%matrix,1)
                    rank=size(block_rand%ButterflyV%blocks(j)%matrix,2)
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
		stats%Mem_int_vec = max(stats%Mem_int_vec,mem_vec) 
	
	
    endif
    

    return
    
end subroutine Init_RandVect_Empty






subroutine Initialize_Bplus_FromInput(Bplus,Bplus_randomized)
	use misc
    use MODULE_FILE
	! use lapack95
	! use blas95
    implicit none
    
    integer level_c,rowblock
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn
    integer dimension_rank, dimension_m, dimension_n, blocks, groupm, groupn,tmpi,tmpj
    real*8 a,b,c,d
    complex (kind=8) ctemp
	complex (kind=8), allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real*8, allocatable:: Singular(:)
	integer mn_min,index_i,index_j,blocks_idx
	type(matrixblock)::block
	type(blockplus)::Bplus,Bplus_randomized
	integer seed_myid(2)
	integer time(8)
	integer ll,bb
	integer:: level_BP,levelm,groupm_start,Nboundall	
	
	
	! allocate(Bplus_randomized)
	
	Bplus_randomized%level = Bplus%level
	Bplus_randomized%col_group = Bplus%col_group
	Bplus_randomized%row_group = Bplus%row_group
	Bplus_randomized%Lplus = Bplus%Lplus
	Bplus_randomized%boundary = Bplus%boundary
	
	allocate(Bplus_randomized%LL(LplusMax))
	

	! Bplus_randomized%LL(1)%Nbound = 1
	! allocate(Bplus_randomized%LL(1)%matrices_block(1))
	! Bplus_randomized%LL(1)%matrices_block(1)%level = Bplus_randomized%level
	! Bplus_randomized%LL(1)%matrices_block(1)%col_group = Bplus_randomized%col_group
	! Bplus_randomized%LL(1)%matrices_block(1)%row_group = Bplus_randomized%row_group			
	! Bplus_randomized%LL(1)%matrices_block(1)%style = Bplus%LL(1)%matrices_block(1)%style
	allocate(Bplus_randomized%LL(1)%boundary_map(1))
	Bplus_randomized%LL(1)%boundary_map(1) = Bplus%LL(1)%boundary_map(1)
	
	do ll=1,LplusMax
	Bplus_randomized%LL(ll)%Nbound = 0
	end do	
	
	do ll=1,LplusMax-1
		if(Bplus%LL(ll)%Nbound>0)then

			Bplus_randomized%LL(ll)%rankmax = Bplus%LL(ll)%rankmax
			Bplus_randomized%LL(ll)%Nbound = Bplus%LL(ll)%Nbound
			allocate(Bplus_randomized%LL(ll)%matrices_block(Bplus_randomized%LL(ll)%Nbound))
			
			do bb =1,Bplus_randomized%LL(ll)%Nbound
				Bplus_randomized%LL(ll)%matrices_block(bb)%row_group = Bplus%LL(ll)%matrices_block(bb)%row_group
				Bplus_randomized%LL(ll)%matrices_block(bb)%col_group = Bplus%LL(ll)%matrices_block(bb)%col_group
				Bplus_randomized%LL(ll)%matrices_block(bb)%style = Bplus%LL(ll)%matrices_block(bb)%style
				Bplus_randomized%LL(ll)%matrices_block(bb)%level = Bplus%LL(ll)%matrices_block(bb)%level
			end do
			
			
			if(Bplus%LL(ll+1)%Nbound==0)then		
				Bplus_randomized%LL(ll+1)%Nbound=0
			else 
				level_butterfly = Bplus%LL(ll)%matrices_block(1)%level_butterfly 
				! level_butterfly = 0 !!! only lowrank																			   
				level_BP = Bplus%level
				levelm = ceiling_safe(dble(level_butterfly)/2d0)						
				groupm_start=Bplus%LL(ll)%matrices_block(1)%row_group*2**levelm		
				! Nboundall = 2**(Bplus%LL(ll)%matrices_block(1)%level+levelm-level_BP)				
				Nboundall=size(Bplus%LL(ll+1)%boundary_map)
				
				allocate(Bplus_randomized%LL(ll+1)%boundary_map(Nboundall))
				! write(*,*)shape(Bplus%LL(ll+1)%boundary_map),shape(Bplus_randomized%LL(ll+1)%boundary_map),'didi',ll
				
				! write(*,*)'gali',ll,Nboundall,shape(Bplus%LL(ll+1)%boundary_map)
				Bplus_randomized%LL(ll+1)%boundary_map = Bplus%LL(ll+1)%boundary_map
			end if
		else 
			exit
		end if

	end do
	
    return

end subroutine Initialize_Bplus_FromInput


subroutine Initialize_Butterfly_randomized(level_butterfly,rankmax,groupm,groupn,block,block_rand)

    use MODULE_FILE
    implicit none
    
    integer level_c,rowblock,kover
	integer i,j,k,level,num_blocks,blocks3,num_row,num_col,ii,jj,kk,level_butterfly, mm, nn
    integer dimension_max,dimension_rank, dimension_m, dimension_n, blocks, groupm, groupm_start,groupn_start,groupn,index_j,index_i
    real*8 a,b,c,d
    complex (kind=8) ctemp
	type(matrixblock)::block,block_rand
	complex (kind=8), allocatable::matrixtemp1(:,:),UU(:,:),VV(:,:)
	real*8, allocatable:: Singular(:)
	integer rankmax
    type(partitionedblocks)::partitioned_block
	
	! allocate (butterfly_block_randomized(1))

    block_rand%level_butterfly=level_butterfly
    num_blocks=2**level_butterfly
	dimension_rank= rankmax 
	
	! write(*,*)dimension_rank
	 

    block_rand%dimension_rank=dimension_rank
	
    block_rand%style=2
    block_rand%row_group=-1
    block_rand%col_group=-1

	block_rand%M = block%M
	block_rand%N = block%N
	block_rand%headm = block%headm
	block_rand%headn = block%headn
	
	
	block_rand%M_loc = block%M_loc
	block_rand%M_loc_db = block%M_loc_db
	block_rand%N_loc = block%N_loc
	block_rand%N_loc_db = block%N_loc_db
	block_rand%pgno = block%pgno
	block_rand%pgno_db = block%pgno_db

	if(associated(block%N_p))then
		allocate(block_rand%N_p(size(block%N_p,1),2))
		block_rand%N_p = block%N_p
	endif
	if(associated(block%M_p))then
		allocate(block_rand%M_p(size(block%M_p,1),2))
		block_rand%M_p = block%M_p
	endif
	if(associated(block%N_p_db))then
		allocate(block_rand%N_p_db(size(block%N_p_db,1),2))
		block_rand%N_p_db = block%N_p_db
	endif
	if(associated(block%M_p_db))then
		allocate(block_rand%M_p_db(size(block%M_p_db,1),2))
		block_rand%M_p_db = block%M_p_db
	endif	
	
	
	
	groupm_start=groupm*2**level_butterfly
	groupn_start=groupn*2**level_butterfly

	
    allocate (block_rand%ButterflyU%blocks(2**level_butterfly))
    allocate (block_rand%ButterflyV%blocks(2**level_butterfly))

	dimension_max = 2*dimension_rank
	do blocks=1, num_blocks	
		
		if(allocated(block%ButterflyU%blocks) .and. size(block%ButterflyU%blocks)==num_blocks)then
			! write(*,*)size(block%ButterflyU),num_blocks,allocated(block%ButterflyU)
			dimension_m=size(block%ButterflyU%blocks(blocks)%matrix,1)
			dimension_n=size(block%ButterflyV%blocks(blocks)%matrix,1)
		else 	
			dimension_m=basis_group(groupm_start+blocks-1)%tail-basis_group(groupm_start+blocks-1)%head+1
			dimension_n=basis_group(groupn_start+blocks-1)%tail-basis_group(groupn_start+blocks-1)%head+1
		endif

		dimension_max = max(dimension_max,dimension_m)	
		dimension_max = max(dimension_max,dimension_n)
		block_rand%ButterflyU%blocks(blocks)%mdim=dimension_m
		block_rand%ButterflyV%blocks(blocks)%mdim=dimension_n
	end do	
	allocate(block_rand%KerInv(dimension_max,dimension_max))
	call RandomMat(dimension_max,dimension_max,dimension_max,block_rand%KerInv,3)	
    
    do blocks=1, num_blocks
		
		if(allocated(block%ButterflyU%blocks) .and. size(block%ButterflyU%blocks)==num_blocks)then
			dimension_m=size(block%ButterflyU%blocks(blocks)%matrix,1)
		else
			dimension_m= basis_group(groupm_start+blocks-1)%tail-basis_group(groupm_start+blocks-1)%head+1
        endif
		
		allocate (block_rand%ButterflyU%blocks(blocks)%matrix(dimension_m,dimension_rank))

		allocate(matrixtemp1(dimension_rank,dimension_m))
		call RandomMat(dimension_rank,dimension_m,min(dimension_m,dimension_rank),matrixtemp1,0)
        
		! !$omp parallel do default(shared) private(i,j)		
		do j=1, dimension_rank
            do i=1, dimension_m
				block_rand%ButterflyU%blocks(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do
		! !$omp end parallel do
		
		deallocate(matrixtemp1)		
		
		if(allocated(block%ButterflyU%blocks) .and. size(block%ButterflyU%blocks)==num_blocks)then
			dimension_n=size(block%ButterflyV%blocks(blocks)%matrix,1)
		else
			dimension_n=basis_group(groupn_start+blocks-1)%tail-basis_group(groupn_start+blocks-1)%head+1
        endif
		
		allocate (block_rand%ButterflyV%blocks(blocks)%matrix(dimension_n,dimension_rank))

		allocate(matrixtemp1(dimension_rank,dimension_n))
		call RandomMat(dimension_rank,dimension_n,min(dimension_n,dimension_rank),matrixtemp1,0)
		
		! !$omp parallel do default(shared) private(i,j)		
        do j=1, dimension_rank
            do i=1, dimension_n
				block_rand%ButterflyV%blocks(blocks)%matrix(i,j) = matrixtemp1(j,i)
			end do
		end do
		! !$omp end parallel do
		
		deallocate(matrixtemp1)		
		
    enddo

    if (level_butterfly/=0) then
        allocate (matrixtemp1(2*dimension_rank,2*dimension_rank))
        allocate (block_rand%ButterflyKerl(level_butterfly))

        do level=1, level_butterfly
            num_row=2**level
            num_col=2**(level_butterfly-level+1)
            block_rand%ButterflyKerl(level)%num_row=num_row
            block_rand%ButterflyKerl(level)%num_col=num_col
            allocate (block_rand%ButterflyKerl(level)%blocks(num_row,num_col))

        enddo
        deallocate (matrixtemp1)
    endif	
    
    return

end subroutine Initialize_Butterfly_randomized
 
 

end module Utilites_randomized
