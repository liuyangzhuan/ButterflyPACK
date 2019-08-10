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
PROGRAM nastran2geo

    implicit none
    
    integer i,j,ii,jj,iii,jjj,k,kk,kkk
    integer node, patch, edge, flag
    integer num_nodes, num_patches
    integer intemp1, intemp2, intemp
    real T0
    character (50) chartemp, chartemp1, chartemp2
	character(len=300)::filename, DATA_DIR    
    real, allocatable :: xyz_nodes(:,:)
    integer, allocatable :: node_patches(:,:)
    
    real a(3),b(3),c(3),r0
    
    T0=secnds(0.0)
    CALL getarg(1, filename)
	open(unit=10,file=trim(filename),status = 'unknown')

!	read(unit=8,fmt=*) filename
	
!    open(unit=10,file=filename,status='unknown')
    
    do i = 1, 1
        read (10,*) 
    enddo
    
    flag=0
    k=0
    do while (flag==0)
         read (10,*) chartemp1
         read (10,*) 
         k=k+1
         if (chartemp1/='GRID*') then
             flag=1
         endif
    enddo    
    num_nodes=k-1
    rewind (10)
    allocate (xyz_nodes(3,num_nodes))
    
    do i = 1, 1
        read (10,*) chartemp
    enddo
    do i= 1, num_nodes
        !read (10,'(A, I, 16X, ES16.9E2, ES16.9E2)') chartemp1, intemp, xyz_nodes(1,i), xyz_nodes(2,i)
        read (10,'(40X,2G16.0)') xyz_nodes(1,i), xyz_nodes(2,i)
        read (10,'(8X,G16.0)') xyz_nodes(3,i)
    enddo
    
    flag=0
    k=0
    do while (flag==0)
         read (10,*) chartemp1
         k=k+1
         if (chartemp1=='ENDDATA') then
             flag=1
         endif
    enddo    
    num_patches=k-1
    rewind (10)
    allocate (node_patches(3,num_patches))
    
     do i = 1, 1
        read (10,*) 
     enddo
     do i= 1, 2*num_nodes
        read (10,*)
    enddo
    
    do i= 1, num_patches
        read (10,'(A8,5I8)') chartemp1, intemp1, intemp2, node_patches(1,i), node_patches(2,i), node_patches(3,i)  
    enddo
     
    close (10)
    
    open(15,file='node.geo')
    open(20,file='elem.geo')
    
    write (15,*) num_nodes
    write (20,*) num_patches
    
    do i=1, num_nodes
        write (15,*) i, xyz_nodes(1:3,i)
    enddo
    
    do i=1, num_patches
        write (20,*) i, node_patches(1:3,i)
    enddo
    
    close (15)
    close (20)
    
    deallocate (xyz_nodes,node_patches)
    
    stop

end PROGRAM nastran2geo
