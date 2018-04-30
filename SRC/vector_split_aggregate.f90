! subroutine split_Butterfly(blocks,flag)
! 
!     use MODULE_FILE
!     implicit none
! 
!     integer blocks, flag, level
!     integer style, blocks_son(4), style_son(4), size_m, size_n
!     integer mm(2), nn(2), rank, num_blocks, index_mm, index_nn, rank
!     integer i, j, k, level_butterfly, son, index_i, index_j, index_ij, index_iijj
!     complex ctemp
! 
!     style=matrices_block(blocks,flag)%style
! !    if (style/=1 .and. style/=2) then
! !        write (*,*) 'split_Butterfly error!'
! !        pause
! !    endif
!     blocks_son(1)=4*blocks+1
!     blocks_son(2)=blocks_son(1)+1
!     blocks_son(3)=blocks_son(2)+1
!     blocks_son(4)=blocks_son(3)+1
!     
!     if (style==1) then
!     
!         size_m=size(matrices_block(blocks,flag)%fullmat,1)
!         size_n=size(matrices_block(blocks,flag)%fullmat,2)
!         mm(1)=int((1+size_m)/2)
!         mm(2)=size_m-mm(1)
!         nn(1)=int((1+size_n)/2)
!         nn(2)=size_n-nn(1)
!         matrices_block(blocks_son(1:4),flag)%style=1
!         allocate (matrices_block(blocks_son(1),flag)%fullmat(mm(1),nn(1)))
!         allocate (matrices_block(blocks_son(2),flag)%fullmat(mm(2),nn(1)))
!         allocate (matrices_block(blocks_son(3),flag)%fullmat(mm(1),nn(2)))
!         allocate (matrices_block(blocks_son(4),flag)%fullmat(mm(2),nn(2)))
!         call copy_matrices('F',blocks_son(1),flag,1,mm(1),1,nn(1),blocks,flag,1,mm(1),1,nn(1))
!         call copy_matrices('F',blocks_son(2),flag,1,mm(2),1,nn(1),blocks,flag,mm(1)+1,size_m,1,nn(1))
!         call copy_matrices('F',blocks_son(3),flag,1,mm(1),1,nn(2),blocks,flag,1,mm(1),nn(1)+1,size_n)
!         call copy_matrices('F',blocks_son(4),flag,1,mm(2),1,nn(2),blocks,flag,mm(1)+1,size_m,nn(1)+1,size_n)
!         deallocate (matrices_block(blocks,flag)%fullmat)
!         matrices_block(blocks,flag)%style=3
!         
!     elseif (style==2) then
!     
!         level_butterfly=matrices_block(blocks,flag)%level_butterfly-1
!         matrices_block(blocks_son(1:4),flag)%style=2
!         allocate (matrices_block(blocks_son(1),flag)%butterflymatrixP(num_blocks,0:level_butterfly)
!         allocate (matrices_block(blocks_son(1),flag)%butterflymatrixB(num_blocks,2)
!         allocate (matrices_block(blocks_son(2),flag)%butterflymatrixP(num_blocks,0:level_butterfly)
!         allocate (matrices_block(blocks_son(2),flag)%butterflymatrixB(num_blocks,2)
!         allocate (matrices_block(blocks_son(3),flag)%butterflymatrixP(num_blocks,0:level_butterfly)
!         allocate (matrices_block(blocks_son(3),flag)%butterflymatrixB(num_blocks,2)
!         allocate (matrices_block(blocks_son(4),flag)%butterflymatrixP(num_blocks,0:level_butterfly)
!         allocate (matrices_block(blocks_son(4),flag)%butterflymatrixB(num_blocks,2)
!         do level=0, level_butterfly
!             index_ij=0
!             index_mm=2**level
!             index_nn=2**(level_butterfly-level)
!             do index_i=1, index_mm
!                 do index_j=1, 2*index_nn
!                     index_ij=index_ij+1
!                     if (index_i<=int((index_mm+1)/2)) then
!                         if(index_j<=index_nn)
!                             index_iijj=(index_i-1)*index_nn+index_j
!                             call copy_matrixP_butterfly(blocks_son(1),level,index_iijj,blocks,level,index_ij)
!                         else
!                             index_iijj=(index_i-1)*index_nn+index_j-index_nn
!                             call copy_matrixP_butterfly(blocks_son(3),level,index_iijj,blocks,level,index_ij)
!                         endif
!                     endif
!                     if (index_i>=int((index_mm/2)+1) then
!                         if (index_j<=index_nn)
!                             index_iijj=(index_i-int(index_mm/2)-1)*index_nn+index_j
!                             call copy_matrixP_butterfly(blocks_son(2),level,index_iijj,blocks,level,index_ij)
!                         else
!                             index_iijj=(index_i-int(index_mm/2)-1)*index_nn+index_j-index_nn
!                             call copy_matrixP_butterfly(blocks_son(2),level,index_iijj,blocks,level,index_ij)
!                         endif
!                     endif
!                 enddo
!             enddo
!             if (level<level_butterfly) then
!                 deallocate (matrices_block(blocks,flag)%butterflymatrixP(index_ij,level)%matrixP)
!             endif
!         enddo
!         
!         num_blocks=2**level_butterfly
!         butterflyB_inuse=matrices_block(blocks,flag)%butterflyB_inuse
!         matrices_block(blocks_son(1:4),flag)%butterflyB_inuse=1
!         do index_ij=1, 2*num_blocks
!             mm(1)=size(matrices_block(blocks,flag)%butterflymatrixB(index_ij,butterflyB_inuse)%matrixB,1)
!             rank=size(matrices_block(blocks,flag)%butterflymatrixB(index_ij,butterflyB_inuse)%matrixB,2)
!             if (mod(index_ij,2)==1) then
!                 nn(1)=size(matrices_block(blocks,flag)%butterflymatrixP(index_ij,level_butterfly)%matrixP,2)
!                 nn(2)=size(matrices_block(blocks,flag)%butterflymatrixP(index_ij+1,level_butterfly)%matrixP,2)
!             else
!                 nn(1)=size(matrices_block(blocks,flag)%butterflymatrixP(index_ij-1,level_butterfly)%matrixP,2)
!                 nn(2)=size(matrices_block(blocks,flag)%butterflymatrixP(index_ij,level_butterfly)%matrixP,2)
!                 deallocate (matrices_block(blocks,flag)%butterflymatrixP(index_ij-1,level_butterfly)%matrixP)
!                 deallocate (matrices_block(blocks,flag)%butterflymatrixP(index_ij,level_butterfly)%matrixP)
!             endif
!             if (index_ij<=num_blocks) then
!                 allocate (matrices_block(blocks_son(1),flag)%butterflymatrixB(index_ij,1)%matrixB(mm(1),nn(1)))
!                 allocate (matrices_block(blocks_son(3),flag)%butterflymatrixB(index_ij,1)%matrixB(mm(1),nn(2)))
!                 !$omp parallel do default(shared) private(i,j,k,ctemp)
!                 do i=1, mm(1)
!                     do j=1, nn(1)
!                         ctemp=(0.,0.)
!                         do k=1, rank
!                             ctemp=ctemp+matrices_block(blocks,flag)%butterflymatrixB(index_ij,butterflyB_inuse)%matrixB(i,k)*matrices_block(blocks,flag)%butterflymatrixP(index_ij,level_butterfly+1)%matrixP(j,k)
!                         enddo
!                         matrices_block(blocks_son(1),flag)%butterflymatrixB(index_ij,1)%matrixB(i,j)=ctemp
!                     enddo
!                     do j=1, nn(2)
!                         ctemp=(0.,0.)
!                         do k=1, rank
!                             ctemp=ctemp+matrices_block(blocks,flag)%butterflymatrixB(index_ij,butterflyB_inuse)%matrixB(i,k)*matrices_block(blocks,flag)%butterflymatrixP(index_ij,level_butterfly+1)%matrixP(j+nn(1),k)
!                         enddo
!                         matrices_block(blocks_son(3),flag)%butterflymatrixB(index_ij,1)%matrixB(i,j)=ctemp
!                     enddo
!                 enddo
!                 !$omp end parallel do
!             else
!                 allocate (matrices_block(blocks_son(2),flag)%butterflymatrixB(index_ij-num_blocks,1)%matrixB(mm(1),nn(1)))
!                 allocate (matrices_block(blocks_son(4),flag)%butterflymatrixB(index_ij-num_blocks,1)%matrixB(mm(1),nn(2)))
!                 !$omp parallel do default(shared) private(i,j,k,ctemp)
!                 do i=1, mm(1)
!                     do j=1, nn(1)
!                         ctemp=(0.,0.)
!                         do k=1, rank
!                             ctemp=ctemp+matrices_block(blocks,flag)%butterflymatrixB(index_ij,butterflyB_inuse)%matrixB(i,k)*matrices_block(blocks,flag)%butterflymatrixP(index_ij,level_butterfly+1)%matrixP(j,k)
!                         enddo
!                         matrices_block(blocks_son(2),flag)%butterflymatrixB(index_ij-num_blocks,1)%matrixB(i,j)=ctemp
!                     enddo
!                     do j=1, nn(2)
!                         ctemp=(0.,0.)
!                         do k=1, rank
!                             ctemp=ctemp+matrices_block(blocks,flag)%butterflymatrixB(index_ij,butterflyB_inuse)%matrixB(i,k)*matrices_block(blocks,flag)%butterflymatrixP(index_ij,level_butterfly+1)%matrixP(j+nn(1),k)
!                         enddo
!                         matrices_block(blocks_son(4),flag)%butterflymatrixB(index_ij-num_blocks,1)%matrixB(i,j)=ctemp
!                     enddo
!                 enddo
!                 !$omp end parallel do
!             endif
!             deallocate (matrices_block(blocks,flag)%butterflymatrixB(index_ij,butterflyB_inuse)%matrixB)
!             deallocate (matrices_block(blocks,flag)%butterflymatrixP(index_ij,level_butterfly+1)%matrixP)
!         enddo
!  
!         deallocate (matrices_block(blocks,flag)%butterflymatrixB, matrices_block(blocks,flag)%butterflymatrixP)
!         matrices_block(blocks,flag)%style=3
!         
!     endif
! 
!     return
! 
! end subroutine split_Butterfly

subroutine split_vectors(vectors)

    use MODULE_FILE
    implicit none

    integer vectors
    integer style, blocks_son(2), style_son(2), size_m, size_n
    integer mm(2), nn
    integer i, j, k
    complex ctemp

    blocks_son(1)=2*vectors
    blocks_son(2)=blocks_son(1)+1
    size_m=size(vectors_block(vectors)%vector,1)
    size_n=size(vectors_block(vectors)%vector,2)
    mm(1)=int((1+size_m)/2)
    mm(2)=size_m-mm(1)
    allocate (vectors_block(blocks_son(1))%vector(mm(1),size_n))
    allocate (vectors_block(blocks_son(2))%vector(mm(2),size_n))
    vectors_block(blocks_son(1:2))%style=1
    call copy_vectors(blocks_son(1),1,mm(1),1,size_n,vectors,1,mm(1),1,size_n)
    call copy_vectors(blocks_son(2),1,mm(2),1,size_n,vectors,mm(1)+1,size_m,1,size_n)
    deallocate (vectors_block(vectors)%vector)
    vectors_block(vectors)%style=3

    return

end subroutine split_vectors

recursive subroutine aggregate_vectors(vectors)

    use MODULE_FILE
    implicit none

    integer vectors
    integer style, blocks_son(2), style_son(2)
    integer mm(2), nn
    integer i, j, k
    complex ctemp

    blocks_son(1)=2*vectors
    blocks_son(2)=blocks_son(1)+1

    if (vectors_block(blocks_son(1))%style==3) then
        call aggregate_vectors(blocks_son(1))
    endif
    if (vectors_block(blocks_son(2))%style==3) then
        call aggregate_vectors(blocks_son(2))
    endif

    mm(1)=size(vectors_block(blocks_son(1))%vector,1)
    mm(2)=size(vectors_block(blocks_son(2))%vector,1)
    nn=size(vectors_block(blocks_son(2))%vector,2)

    allocate (vectors_block(vectors)%vector(mm(1)+mm(2),nn))
    call copy_vectors(vectors,1,mm(1),1,nn,blocks_son(1),1,mm(1),1,nn)
    call copy_vectors(vectors,mm(1)+1,mm(1)+mm(2),1,nn,blocks_son(2),1,mm(2),1,nn)
    vectors_block(vectors)%style=1
    vectors_block(blocks_son(1:2))%style=4
    deallocate (vectors_block(blocks_son(1))%vector,vectors_block(blocks_son(2))%vector)

    return

end subroutine aggregate_vectors

subroutine copy_vectors (vectors_1,headm_1,tailm_1,headn_1,tailn_1,vectors_2,headm_2,tailm_2,headn_2,tailn_2)

    use MODULE_FILE
    implicit none

    integer vectors_1,vectors_2,headm_1,headn_1,tailm_1,tailn_1,headm_2,headn_2,tailm_2,tailn_2
    integer i,j,mm,nn

    mm=tailm_1-headm_1+1
    nn=tailn_1-headn_1+1

     !$omp parallel do default(shared) private(i,j)
    do j=1,nn
        do i=1,mm
            vectors_block(vectors_1)%vector(headm_1+i-1,headn_1+j-1)=vectors_block(vectors_2)%vector(headm_2+i-1,headn_2+j-1)
        enddo
    enddo
    !$omp end parallel do

    return

end subroutine copy_vectors

! subroutine copy_matrixP_butterfly(blocks_son,level_son,index_iijj,blocks,level,index_ij) 
! 
!     use MODULE_FILE
!     implicit none
!     
!     integer blocks_son, level_son, blocks, level, index_ij, index_iijj
!     integer i, j, k, ii, jj, kk, mm, nn
!     
!     mm=size(matrices_block(blocks,0)%butterflymatrixP(index_ij,level)%matrixP,1)
!     nn=size(matrices_block(blocks,0)%butterflymatrixP(index_ij,level)%matrixP,2)
!     
!     allocate (matrices_block(blocks_son,0)%butterflymatrixP(index_iijj,level_son)%matrixP(mm,nn))
!     
!     !$omp parallel do default(shared) private(i,j)
!     do j=1, nn
!         do i=1, mm
!             matrices_block(blocks_son,0)%butterflymatrixP(index_iijj,level_son)%matrixP(i,j)=matrices_block(blocks,0)%butterflymatrixP(index_ij,level)%matrixP(i,j)
!         enddo
!     enddo
!     !$omp end parallel do
!     
!     return    
! 
! end subroutine copy_matrixP_butterfly

