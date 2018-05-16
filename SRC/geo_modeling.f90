module geometry_model
contains 


subroutine geo_modeling()

    use MODULE_FILE
    implicit none
	if(Kernel==EMCURV)then
		call geo_modeling_CURV()
	elseif(Kernel==EMSURF)then
		call geo_modeling_SURF()
	else
		write(*,*)'unknown Kernel for geo_modeling'
		stop
	endif
end subroutine geo_modeling	

subroutine geo_modeling_CURV()

    use MODULE_FILE
    implicit none
    
    integer i,j,ii,jj,iii,jjj
    integer intemp
    integer node, patch, edge, flag
    integer node1, node2, num_node
    integer node_temp(2)
    real(kind=8) dx, xx, yy, rr, theta,L,M,Am,tt,L1,L2,L3
    
    real(kind=8),allocatable :: node_xy_original(:,:)
    integer,allocatable :: num_edge_of_node(:)
    
    real(kind=8) a(3),b(3),c(3),r0, phi_start
    

	Ncorner = 0
	
    if (geo_model==1) then  !*****single strip*****
        
        Delta_ll=2d0/Maxedge
        Maxnode=2*Maxedge+1
        allocate (xyz(2,0:Maxnode-1), node_patch_of_edge(0:2,Maxedge))
        xyz(1,0)=-1d0 ; xyz(2,0)=0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxnode-1
            dx=node*Delta_ll/2
            xyz(1,node)=dx-1d0
            xyz(2,node)=0
        enddo
        !$omp end parallel do
        node_patch_of_edge(1,1)=0 ; node_patch_of_edge(2,1)=2 ; node_patch_of_edge(0,1)=1   ! 1 start 0 center 2 end
        do edge=2, Maxedge
            node_patch_of_edge(1,edge)=node_patch_of_edge(2,edge-1)
            node_patch_of_edge(2,edge)=node_patch_of_edge(2,edge-1)+2
            node_patch_of_edge(0,edge)=node_patch_of_edge(0,edge-1)+2
        enddo           
        
    elseif (geo_model==2) then  !*****corner reflector*****
    
        Delta_ll=2d0/Maxedge
        Maxnode=2*Maxedge+1
        allocate (xyz(2,0:Maxnode-1), node_patch_of_edge(0:2,Maxedge))
        
        xyz(1,0)=-1d0/sqrt(2d0) ; xyz(2,0)=1d0/sqrt(2d0)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*Delta_ll/2
            xyz(1,node)=(dx-1d0)/sqrt(2d0)
            xyz(2,node)=(1d0-dx)/sqrt(2d0)
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*Delta_ll/2
            xyz(1,node+Maxedge)=dx/sqrt(2d0)
            xyz(2,node+Maxedge)=dx/sqrt(2d0)
        enddo
        !$omp end parallel do
        node_patch_of_edge(1,1)=0 ; node_patch_of_edge(2,1)=2 ; node_patch_of_edge(0,1)=1
        do edge=2, Maxedge
            node_patch_of_edge(1,edge)=node_patch_of_edge(2,edge-1)
            node_patch_of_edge(2,edge)=node_patch_of_edge(2,edge-1)+2
            node_patch_of_edge(0,edge)=node_patch_of_edge(0,edge-1)+2
        enddo
		
		Ncorner = 1
		allocate(corner_points(2,Ncorner))
		corner_points(1,1) = 0
		corner_points(2,1) = 0        
    
    elseif (geo_model==3) then  !*****two opposite strips*****
    
        Delta_ll=2d0/Maxedge
        Maxnode=2*Maxedge+2
        allocate (xyz(2,0:Maxnode-1), node_patch_of_edge(0:2,Maxedge))
        
        xyz(1,0)=-0.5d0 ; xyz(2,0)=1.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*Delta_ll/2
            xyz(1,node)=dx-0.5d0
            xyz(2,node)=1.0d0
        enddo
        !$omp end parallel do
        xyz(1,Maxedge+1)=-0.5d0 ; xyz(2,Maxedge+1)=-1.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*Delta_ll/2
            xyz(1,node+Maxedge+1)=dx-0.5d0
            xyz(2,node+Maxedge+1)=-1.0d0
        enddo
        !$omp end parallel do
        
        intemp=int(Maxedge/2)
        
        node_patch_of_edge(1,1)=0 ; node_patch_of_edge(2,1)=2 ; node_patch_of_edge(0,1)=1
        do edge=2, intemp
            node_patch_of_edge(1,edge)=node_patch_of_edge(2,edge-1)
            node_patch_of_edge(2,edge)=node_patch_of_edge(2,edge-1)+2
            node_patch_of_edge(0,edge)=node_patch_of_edge(0,edge-1)+2
        enddo
        
        node_patch_of_edge(1,intemp+1)=Maxnode/2 ; node_patch_of_edge(2,intemp+1)=Maxnode/2+2 ; node_patch_of_edge(0,intemp+1)=Maxnode/2+1
        do edge=intemp+2, Maxedge
            node_patch_of_edge(1,edge)=node_patch_of_edge(2,edge-1)
            node_patch_of_edge(2,edge)=node_patch_of_edge(2,edge-1)+2
            node_patch_of_edge(0,edge)=node_patch_of_edge(0,edge-1)+2
        enddo
    
    
    elseif (geo_model==5) then !************cylinder*****************
        
        Delta_ll=  2.0d0*pi/Maxedge
        Maxnode=2*Maxedge
        allocate (xyz(2,0:Maxnode-1), node_patch_of_edge(0:2,Maxedge))
        xyz(1,0)=1.0d0 ; xyz(2,0)=0.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxnode-1
            dx=node*Delta_ll/2
            xyz(1,node)=cos(dx)
            xyz(2,node)=sin(dx)
        enddo
        !$omp end parallel do
        node_patch_of_edge(1,1)=0 ; node_patch_of_edge(2,1)=2 ; node_patch_of_edge(0,1)=1
        do edge=2, Maxedge-1
            node_patch_of_edge(1,edge)=node_patch_of_edge(2,edge-1)
            node_patch_of_edge(2,edge)=node_patch_of_edge(2,edge-1)+2
            node_patch_of_edge(0,edge)=node_patch_of_edge(0,edge-1)+2
        enddo     
        node_patch_of_edge(1,Maxedge)=node_patch_of_edge(2,Maxedge-1)
        node_patch_of_edge(2,Maxedge)=0
        node_patch_of_edge(0,Maxedge)=node_patch_of_edge(0,Maxedge-1)+2      
        
    elseif (geo_model==6) then !***********cavity******************
    
        Delta_ll=11d0/3d0/Maxedge
        Maxnode=2*Maxedge+1
        allocate (xyz(2,0:Maxnode-1), node_patch_of_edge(0:2,Maxedge))
        dx=Delta_ll/2d0
        xyz(1,0)=-1d0/6. ; xyz(2,0)=1.0d0
        node=0
        do while (flag==0) 
            node=node+1
            xyz(1,node)=xyz(1,node-1)-dx
            xyz(2,node)=xyz(2,node-1)
            if (xyz(1,node)<=-0.5d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            xyz(1,node)=xyz(1,node-1)
            xyz(2,node)=xyz(2,node-1)-dx
            if (xyz(2,node)<=0.) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            xyz(1,node)=xyz(1,node-1)+dx
            xyz(2,node)=xyz(2,node-1)
            if (xyz(1,node)>=0.5d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            xyz(1,node)=xyz(1,node-1)
            xyz(2,node)=xyz(2,node-1)+dx
            if (xyz(2,node)>=1.0d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            xyz(1,node)=xyz(1,node-1)-dx
            xyz(2,node)=xyz(2,node-1)
            if (node>=Maxnode-1) then
                flag=1
            endif
        enddo
        
        node_patch_of_edge(1,1)=0 ; node_patch_of_edge(2,1)=2 ; node_patch_of_edge(0,1)=1
        do edge=2, Maxedge
            node_patch_of_edge(1,edge)=node_patch_of_edge(2,edge-1)
            node_patch_of_edge(2,edge)=node_patch_of_edge(2,edge-1)+2
            node_patch_of_edge(0,edge)=node_patch_of_edge(0,edge-1)+2
        enddo

		Ncorner = 4
		allocate(corner_points(2,Ncorner))
		corner_points(1,1) = -0.5
		corner_points(2,1) = 0     
		corner_points(1,2) = 0.5
		corner_points(2,2) = 0
		corner_points(1,3) = 0.5
		corner_points(2,3) = 1
		corner_points(1,4) = -0.5
		corner_points(2,4) = 1		
		
    elseif (geo_model==7) then   !************open cylinder*****************
    
	
	
        Delta_ll=1d0*pi/Maxedge !2.0d0*pi*5d0/6d0/Maxedge
        Maxnode=2*Maxedge+1
        allocate (xyz(2,0:Maxnode-1), node_patch_of_edge(0:2,Maxedge))
        
        xyz(1,0)=cos(0*pi) ; xyz(2,0)=sin(0*pi)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*Delta_ll
            xyz(1,node*2)=cos(0*pi+dx)
            xyz(2,node*2)=sin(0*pi+dx)
        enddo
        !$omp end parallel do
		
		
        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            xyz(1,node*2-1)=(xyz(1,node*2-2)+xyz(1,node*2))/2d0
            xyz(2,node*2-1)=(xyz(2,node*2-2)+xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do		
		
		
        node_patch_of_edge(1,1)=0 ; node_patch_of_edge(2,1)=2 ; node_patch_of_edge(0,1)=1
        do edge=2, Maxedge
            node_patch_of_edge(1,edge)=node_patch_of_edge(2,edge-1)
            node_patch_of_edge(2,edge)=node_patch_of_edge(2,edge-1)+2
            node_patch_of_edge(0,edge)=node_patch_of_edge(0,edge-1)+2
        enddo           

    elseif (geo_model==8) then   !************corrugated open cylinder*****************
		M = 0.2d0*wavelength
		L = 1.5d0*wavelength
        phi_start = 3d0/2d0*pi
		
		Delta_ll=1d0*pi/Maxedge !2.0d0*pi*5d0/6d0/Maxedge
        Maxnode=2*Maxedge+1
        allocate (xyz(2,0:Maxnode-1), node_patch_of_edge(0:2,Maxedge))
        
        xyz(1,0)=cos(phi_start) ; xyz(2,0)=sin(phi_start)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*Delta_ll
            xyz(1,node*2)=(1+M*sin(2*pi*dx/L))*cos(phi_start+dx)
            xyz(2,node*2)=(1+M*sin(2*pi*dx/L))*sin(phi_start+dx)
        enddo
        !$omp end parallel do
		
		
        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            xyz(1,node*2-1)=(xyz(1,node*2-2)+xyz(1,node*2))/2d0
            xyz(2,node*2-1)=(xyz(2,node*2-2)+xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do		
		
		
        node_patch_of_edge(1,1)=0 ; node_patch_of_edge(2,1)=2 ; node_patch_of_edge(0,1)=1
        do edge=2, Maxedge
            node_patch_of_edge(1,edge)=node_patch_of_edge(2,edge-1)
            node_patch_of_edge(2,edge)=node_patch_of_edge(2,edge-1)+2
            node_patch_of_edge(0,edge)=node_patch_of_edge(0,edge-1)+2
        enddo    

		
    elseif (geo_model==9) then  !*****corrugated corner reflector*****
		M = 0.2d0*wavelength
		L = 1.5d0*wavelength
		
        Delta_ll=2d0/Maxedge
        Maxnode=2*Maxedge+1
        allocate (xyz(2,0:Maxnode-1), node_patch_of_edge(0:2,Maxedge))
        

		Am = M*sin(pi/2-2*pi/L)-M		
		xyz(1,0)=-1d0/sqrt(2d0)+Am/sqrt(2d0); xyz(2,0)=1d0/sqrt(2d0)+Am/sqrt(2d0)
		
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge/2
            dx=node*Delta_ll
			Am = M*sin(2*pi*dx/L+pi/2-2*pi/L)-M
			
            xyz(1,node*2)=(dx-1d0)/sqrt(2d0)+Am/sqrt(2d0)
            xyz(2,node*2)=(1d0-dx)/sqrt(2d0)+Am/sqrt(2d0)
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge/2
            dx=node*Delta_ll
			Am = M*sin(2*pi*(dx+1)/L+pi/2-2*pi/L)-M
			
            xyz(1,node*2+Maxedge)=dx/sqrt(2d0)-Am/sqrt(2d0)
            xyz(2,node*2+Maxedge)=dx/sqrt(2d0)+Am/sqrt(2d0)
        enddo
        !$omp end parallel do

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            xyz(1,node*2-1)=(xyz(1,node*2-2)+xyz(1,node*2))/2d0
            xyz(2,node*2-1)=(xyz(2,node*2-2)+xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        node_patch_of_edge(1,1)=0 ; node_patch_of_edge(2,1)=2 ; node_patch_of_edge(0,1)=1
        do edge=2, Maxedge
            node_patch_of_edge(1,edge)=node_patch_of_edge(2,edge-1)
            node_patch_of_edge(2,edge)=node_patch_of_edge(2,edge-1)+2
            node_patch_of_edge(0,edge)=node_patch_of_edge(0,edge-1)+2
        enddo  		

		Ncorner = 1
		allocate(corner_points(2,Ncorner))
		corner_points(1,1) = 0
		corner_points(2,1) = 0
		
    elseif (geo_model==10) then  !*****cup*****

		L1 = 0
		L2 = sqrt(2d0)
		L3 = sqrt(2d0)
		L = 3
        Delta_ll=(2*L1+2*L2+2*L3+L)/Maxedge
        Maxnode=2*Maxedge+1
        allocate (xyz(2,0:Maxnode-1), node_patch_of_edge(0:2,Maxedge))
        
		
        !$omp parallel do default(shared) private(node,dx,tt)
        do node=0, Maxedge
            dx=node*Delta_ll
			if(dx<=L1)then
				tt = dx
				xyz(1,node*2)= -(L/2-(L1-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				xyz(2,node*2)= (L1+L2+L3)/sqrt(2d0)-tt/sqrt(2d0)			
			else if(dx<=L1+L2)then
				tt = dx - L1
				xyz(1,node*2)= -(L/2-(-L2+L3)/sqrt(2d0))+tt/sqrt(2d0)
				xyz(2,node*2)= (L2+L3)/sqrt(2d0)-tt/sqrt(2d0)
			else if(dx<=L1+L2+L3)then
				tt = dx - (L1+L2)
				xyz(1,node*2)= -(L/2-L3/sqrt(2d0))-tt/sqrt(2d0)
				xyz(2,node*2)= (L3)/sqrt(2d0)-tt/sqrt(2d0)				
			else if(dx<=L1+L2+L3+L)then
				tt = dx - (L1+L2+L3)
				xyz(1,node*2)= -L/2+tt
				xyz(2,node*2)= 0
			else if(dx<=L1+L2+L3+L+L3)then
				tt = dx - (L1+L2+L3+L)
				xyz(1,node*2)= L/2-tt/sqrt(2d0)
				xyz(2,node*2)= tt/sqrt(2d0)		
			else if(dx<=L1+L2+L3+L+L3+L2)then
				tt = dx - (L1+L2+L3+L+L3)
				xyz(1,node*2)= (L/2-L3/sqrt(2d0))+tt/sqrt(2d0)
				xyz(2,node*2)= (L3)/sqrt(2d0)+tt/sqrt(2d0)				
			else
				tt = dx - (L1+L2+L3+L+L3+L2)
				xyz(1,node*2)= (L/2-(-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				xyz(2,node*2)= (L3+L2)/sqrt(2d0)+tt/sqrt(2d0)							
			end if 
        enddo
        !$omp end parallel do			

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            xyz(1,node*2-1)=(xyz(1,node*2-2)+xyz(1,node*2))/2d0
            xyz(2,node*2-1)=(xyz(2,node*2-2)+xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        node_patch_of_edge(1,1)=0 ; node_patch_of_edge(2,1)=2 ; node_patch_of_edge(0,1)=1
        do edge=2, Maxedge
            node_patch_of_edge(1,edge)=node_patch_of_edge(2,edge-1)
            node_patch_of_edge(2,edge)=node_patch_of_edge(2,edge-1)+2
            node_patch_of_edge(0,edge)=node_patch_of_edge(0,edge-1)+2
        enddo  		

		! Ncorner = 4
		! allocate(corner_points(2,Ncorner))
		! corner_points(1,1) = -L/2
		! corner_points(2,1) = 0
		! corner_points(1,2) = L/2
		! corner_points(2,2) = 0
		! corner_points(1,3) = -(L/2-L3/sqrt(2d0))
		! corner_points(2,3) = (L3)/sqrt(2d0)		
		! corner_points(1,4) = (L/2-L3/sqrt(2d0))
		! corner_points(2,4) = (L3)/sqrt(2d0)			

    elseif (geo_model==11) then  !*****longer cup*****
		L1 = 1
		L2 = sqrt(2d0)
		L3 = sqrt(2d0)
		L = 3.5
        Delta_ll=(2*L1+2*L2+2*L3+L)/Maxedge
        Maxnode=2*Maxedge+1
        allocate (xyz(2,0:Maxnode-1), node_patch_of_edge(0:2,Maxedge))
        
		
        !$omp parallel do default(shared) private(node,dx,tt)
        do node=0, Maxedge
            dx=node*Delta_ll
			if(dx<=L1)then
				tt = dx
				xyz(1,node*2)= -(L/2-(L1-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				xyz(2,node*2)= (L1+L2+L3)/sqrt(2d0)-tt/sqrt(2d0)			
			else if(dx<=L1+L2)then
				tt = dx - L1
				xyz(1,node*2)= -(L/2-(-L2+L3)/sqrt(2d0))+tt/sqrt(2d0)
				xyz(2,node*2)= (L2+L3)/sqrt(2d0)-tt/sqrt(2d0)
			else if(dx<=L1+L2+L3)then
				tt = dx - (L1+L2)
				xyz(1,node*2)= -(L/2-L3/sqrt(2d0))-tt/sqrt(2d0)
				xyz(2,node*2)= (L3)/sqrt(2d0)-tt/sqrt(2d0)				
			else if(dx<=L1+L2+L3+L)then
				tt = dx - (L1+L2+L3)
				xyz(1,node*2)= -L/2+tt
				xyz(2,node*2)= 0
			else if(dx<=L1+L2+L3+L+L3)then
				tt = dx - (L1+L2+L3+L)
				xyz(1,node*2)= L/2-tt/sqrt(2d0)
				xyz(2,node*2)= tt/sqrt(2d0)		
			else if(dx<=L1+L2+L3+L+L3+L2)then
				tt = dx - (L1+L2+L3+L+L3)
				xyz(1,node*2)= (L/2-L3/sqrt(2d0))+tt/sqrt(2d0)
				xyz(2,node*2)= (L3)/sqrt(2d0)+tt/sqrt(2d0)				
			else
				tt = dx - (L1+L2+L3+L+L3+L2)
				xyz(1,node*2)= (L/2-(-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				xyz(2,node*2)= (L3+L2)/sqrt(2d0)+tt/sqrt(2d0)							
			end if 
        enddo
        !$omp end parallel do			

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            xyz(1,node*2-1)=(xyz(1,node*2-2)+xyz(1,node*2))/2d0
            xyz(2,node*2-1)=(xyz(2,node*2-2)+xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        node_patch_of_edge(1,1)=0 ; node_patch_of_edge(2,1)=2 ; node_patch_of_edge(0,1)=1
        do edge=2, Maxedge
            node_patch_of_edge(1,edge)=node_patch_of_edge(2,edge-1)
            node_patch_of_edge(2,edge)=node_patch_of_edge(2,edge-1)+2
            node_patch_of_edge(0,edge)=node_patch_of_edge(0,edge-1)+2
        enddo  		

		Ncorner = 6
		allocate(corner_points(2,Ncorner))
		corner_points(1,1) = -(L/2-(-L2+L3)/sqrt(2d0))
		corner_points(2,1) = (L2+L3)/sqrt(2d0)
		corner_points(1,2) = (L/2-(-L2+L3)/sqrt(2d0))
		corner_points(2,2) = (L2+L3)/sqrt(2d0)	
		corner_points(1,3) = -L/2
		corner_points(2,3) = 0
		corner_points(1,4) = L/2
		corner_points(2,4) = 0
		corner_points(1,5) = -(L/2-L3/sqrt(2d0))
		corner_points(2,5) = (L3)/sqrt(2d0)		
		corner_points(1,6) = (L/2-L3/sqrt(2d0))
		corner_points(2,6) = (L3)/sqrt(2d0)		
    endif
    
	maxedgelength = 0
	do edge=1,Maxedge
		maxedgelength = max(maxedgelength,sqrt(sum(abs(xyz(:,node_patch_of_edge(1,edge))-xyz(:,node_patch_of_edge(2,edge))))**2))
	end do	
	
	minedgelength = 10000000
	do edge=1,Maxedge
		minedgelength = min(minedgelength,sqrt(sum(abs(xyz(:,node_patch_of_edge(1,edge))-xyz(:,node_patch_of_edge(2,edge)))**2)))
	end do		
	
	corner_radius = 2*maxedgelength
	
	! write(*,*)	xyz(1,1:100),sum(xyz(1,:))
	! stop

	
		
    write (*,*) ''
    write (*,*) 'Maxedge:',Maxedge
	write (*,*) 'dx:',Delta_ll/2
	write (*,*) 'wavelength/maxedgelength:',wavelength/maxedgelength
    open (256,file='Info.txt',position='append')
    write (256,*) 'Maxedge:',Maxedge
    close (256)
    write (*,*) '' 
    
    return
    
end subroutine geo_modeling_CURV


subroutine geo_modeling_SURF()

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
	
	! write(*,*)	xyz(1,1:100),sum(xyz(1,:))
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
    
end subroutine geo_modeling_SURF

end module geometry_model