module geometry_model
contains 

subroutine geo_modeling()

    use MODULE_FILE
	use misc
    implicit none
    
    integer i,j,ii,jj,iii,jjj
    integer intemp
    integer node, patch, edge, flag
    integer node1, node2
    integer node_temp(2)
	integer Dimn
    real*8 a(3),b(3),c(3),r0
    Dimn=3
    
    
    open(11,file=trim(DATA_DIR)//'/node.geo')
    open(111,file=trim(DATA_DIR)//'/elem.geo')
    
    read(11,*)Maxnode
    read(111,*)Maxpatch
    Maxedge=Maxpatch*3/2
    
    allocate(xyz(3,maxnode+Maxedge))
    allocate(node_of_patch(0:3,maxpatch),node_patch_of_edge(0:6,maxedge+1000))
    allocate(normal_of_patch(3,maxpatch))
    
    
    !************xyz****************
    i=1
    do while(i<=maxnode)
        read(11,*)intemp,xyz(1:3,i)
        xyz(1:3,i)=xyz(1:3,i)/Scale
        i=i+1
    enddo
    close(11)
    
    i=1
    if (mesh_normal==1) then
        do while(i<=maxpatch)
            read(111,*)intemp,node_of_patch(1:3,i)
            i=i+1 
        enddo
    elseif (mesh_normal==-1) then
        do while(i<=maxpatch)
            read(111,*)intemp,node_of_patch(3,i),node_of_patch(2,i),node_of_patch(1,i)
            i=i+1 
        enddo
    endif
    close(111)
    
    !************normal_of_patch****************
    
    !$omp parallel do default(shared) private(patch,i,a,b,c,r0)
    do patch=1,maxpatch
        do i=1,3
            a(i)=(xyz(i,node_of_patch(2,patch))-xyz(i,node_of_patch(1,patch)))
            b(i)=(xyz(i,node_of_patch(3,patch))-xyz(i,node_of_patch(1,patch)))
        enddo
        call curl(a,b,c)
        r0=sqrt(c(1)**2+c(2)**2+c(3)**2)
        c(1)=c(1)/r0
        c(2)=c(2)/r0
        c(3)=c(3)/r0
        normal_of_patch(1:3,patch)=c(1:3)	    
    enddo
    !$omp end parallel do
    
    !************node_patch_of_edge****************

    edge=0
    do i=1,maxpatch-1
        do j=i+1,maxpatch
            flag=0;node1=0;node2=0;iii=1
            do ii=1,3
                do jj=1,3
	     	         if(node_of_patch(ii,i)==node_of_patch(jj,j))then
                        flag=flag+1
                        node_temp(iii)=node_of_patch(ii,i)
                        iii=iii+1
                    endif
                enddo
            enddo
            if(flag==2)then
                edge=edge+1
                if(node_temp(1)<node_temp(2)) then
                    node_patch_of_edge(1,edge)=node_temp(1)
                    node_patch_of_edge(2,edge)=node_temp(2)
                else
                    node_patch_of_edge(1,edge)=node_temp(2)
                    node_patch_of_edge(2,edge)=node_temp(1)
                endif
                node_patch_of_edge(3,edge)=i
                node_patch_of_edge(4,edge)=j       ! notice that : i<j  
                node_patch_of_edge(0,edge)=0
            endif
        enddo
    enddo
    
    Maxedge=edge
    
    !$omp parallel do default(shared) private(edge,node_temp,jj,iii,jjj)
    do edge=1,maxedge
	    node_temp(1)=0
	    node_temp(2)=0	    
	    do jj=3,4
             do iii=1,3
                 do jjj=1,2
        	            if(node_of_patch(iii,node_patch_of_edge(jj,edge))==node_patch_of_edge(jjj,edge)) then
                           node_temp(jj-2)=node_temp(jj-2)+iii               
        	            endif
                 enddo
             enddo
         enddo
         node_temp(1)=6-node_temp(1)
         node_temp(2)=6-node_temp(2)
         node_patch_of_edge(5,edge)=node_of_patch(node_temp(1),node_patch_of_edge(3,edge))
         node_patch_of_edge(6,edge)=node_of_patch(node_temp(2),node_patch_of_edge(4,edge))
    enddo
    !$omp end parallel do
    
    node=Maxnode
    do edge=1, Maxedge
        node=node+1
        node_patch_of_edge(0,edge)=node
        do i=1,3
            xyz(i,node)=1./2.*(xyz(i,node_patch_of_edge(1,edge))+xyz(i,node_patch_of_edge(2,edge)))
        enddo
    enddo
	
	
	
	
	
	maxedgelength = 0
	do edge=1,Maxedge
		maxedgelength = max(maxedgelength,sqrt(sum(abs(xyz(:,node_patch_of_edge(1,edge))-xyz(:,node_patch_of_edge(2,edge)))**2)))
	end do	

	minedgelength = 10000000
	do edge=1,Maxedge
		minedgelength = min(minedgelength,sqrt(sum(abs(xyz(:,node_patch_of_edge(1,edge))-xyz(:,node_patch_of_edge(2,edge)))**2)))
	end do	
	
	! write(*,*)	node_xy(1,1:100),sum(node_xy(1,:))
	! stop

	
		
    write (*,*) ''
    write (*,*) 'Maxedge:',Maxedge
	write (*,*) 'minedgelength:',minedgelength
	write (*,*) 'wavelength/minedgelength:',wavelength/minedgelength
	write (*,*) 'maxedgelength:',maxedgelength
	write (*,*) 'wavelength/maxedgelength:',wavelength/maxedgelength
    open (256,file='Info.txt',position='append')
    write (256,*) 'Maxedge:',Maxedge
    close (256)
    write (*,*) '' 
    
    return
    
end subroutine geo_modeling

end module geometry_model