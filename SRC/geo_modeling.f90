module geometry_model
contains 


subroutine geo_modeling_CURV(msh,ker,ptree)

    use MODULE_FILE
    implicit none
    type(mesh)::msh
	type(kernelquant)::ker
    integer i,j,ii,jj,iii,jjj
    integer intemp
    integer node, patch, edge, flag
    integer node1, node2, num_node,Maxedge
    integer node_temp(2)
    real(kind=8) dx, xx, yy, rr, theta,L,M,Am,tt,L1,L2,L3
    
    real(kind=8),allocatable :: node_xy_original(:,:)
    integer,allocatable :: num_edge_of_node(:)
    
    real(kind=8) a(3),b(3),c(3),r0, phi_start
	type(proctree)::ptree
	
	Maxedge	= msh%Nunk

	msh%Ncorner = 0
	
    if (msh%model2d==1) then  !*****single strip*****
        
        msh%Delta_ll=2d0/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        msh%xyz(1,0)=-1d0 ; msh%xyz(2,0)=0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, msh%maxnode-1
            dx=node*msh%Delta_ll/2
            msh%xyz(1,node)=dx-1d0
            msh%xyz(2,node)=0
        enddo
        !$omp end parallel do
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1   ! 1 start 0 center 2 end
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo           
        
    elseif (msh%model2d==2) then  !*****corner reflector*****
    
        msh%Delta_ll=2d0/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        
        msh%xyz(1,0)=-1d0/sqrt(2d0) ; msh%xyz(2,0)=1d0/sqrt(2d0)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*msh%Delta_ll/2
            msh%xyz(1,node)=(dx-1d0)/sqrt(2d0)
            msh%xyz(2,node)=(1d0-dx)/sqrt(2d0)
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*msh%Delta_ll/2
            msh%xyz(1,node+Maxedge)=dx/sqrt(2d0)
            msh%xyz(2,node+Maxedge)=dx/sqrt(2d0)
        enddo
        !$omp end parallel do
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo
		
		msh%Ncorner = 1
		allocate(msh%corner_points(2,msh%Ncorner))
		msh%corner_points(1,1) = 0
		msh%corner_points(2,1) = 0        
    
    elseif (msh%model2d==3) then  !*****two opposite strips*****
    
        msh%Delta_ll=2d0/Maxedge
        msh%maxnode=2*Maxedge+2
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        
        msh%xyz(1,0)=-0.5d0 ; msh%xyz(2,0)=1.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*msh%Delta_ll/2
            msh%xyz(1,node)=dx-0.5d0
            msh%xyz(2,node)=1.0d0
        enddo
        !$omp end parallel do
        msh%xyz(1,Maxedge+1)=-0.5d0 ; msh%xyz(2,Maxedge+1)=-1.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*msh%Delta_ll/2
            msh%xyz(1,node+Maxedge+1)=dx-0.5d0
            msh%xyz(2,node+Maxedge+1)=-1.0d0
        enddo
        !$omp end parallel do
        
        intemp=int(Maxedge/2)
        
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, intemp
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo
        
        msh%info_unk(1,intemp+1)=msh%maxnode/2 ; msh%info_unk(2,intemp+1)=msh%maxnode/2+2 ; msh%info_unk(0,intemp+1)=msh%maxnode/2+1
        do edge=intemp+2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo
    
    
    elseif (msh%model2d==5) then !************cylinder*****************
        
        msh%Delta_ll=  2.0d0*pi/Maxedge
        msh%maxnode=2*Maxedge
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        msh%xyz(1,0)=1.0d0 ; msh%xyz(2,0)=0.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, msh%maxnode-1
            dx=node*msh%Delta_ll/2
            msh%xyz(1,node)=cos(dx)
            msh%xyz(2,node)=sin(dx)
        enddo
        !$omp end parallel do
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge-1
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo     
        msh%info_unk(1,Maxedge)=msh%info_unk(2,Maxedge-1)
        msh%info_unk(2,Maxedge)=0
        msh%info_unk(0,Maxedge)=msh%info_unk(0,Maxedge-1)+2      
        
    elseif (msh%model2d==6) then !***********cavity******************
    
        msh%Delta_ll=11d0/3d0/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        dx=msh%Delta_ll/2d0
        msh%xyz(1,0)=-1d0/6. ; msh%xyz(2,0)=1.0d0
        node=0
        do while (flag==0) 
            node=node+1
            msh%xyz(1,node)=msh%xyz(1,node-1)-dx
            msh%xyz(2,node)=msh%xyz(2,node-1)
            if (msh%xyz(1,node)<=-0.5d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            msh%xyz(1,node)=msh%xyz(1,node-1)
            msh%xyz(2,node)=msh%xyz(2,node-1)-dx
            if (msh%xyz(2,node)<=0.) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            msh%xyz(1,node)=msh%xyz(1,node-1)+dx
            msh%xyz(2,node)=msh%xyz(2,node-1)
            if (msh%xyz(1,node)>=0.5d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            msh%xyz(1,node)=msh%xyz(1,node-1)
            msh%xyz(2,node)=msh%xyz(2,node-1)+dx
            if (msh%xyz(2,node)>=1.0d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            msh%xyz(1,node)=msh%xyz(1,node-1)-dx
            msh%xyz(2,node)=msh%xyz(2,node-1)
            if (node>=msh%maxnode-1) then
                flag=1
            endif
        enddo
        
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo

		msh%Ncorner = 4
		allocate(msh%corner_points(2,msh%Ncorner))
		msh%corner_points(1,1) = -0.5
		msh%corner_points(2,1) = 0     
		msh%corner_points(1,2) = 0.5
		msh%corner_points(2,2) = 0
		msh%corner_points(1,3) = 0.5
		msh%corner_points(2,3) = 1
		msh%corner_points(1,4) = -0.5
		msh%corner_points(2,4) = 1		
		
    elseif (msh%model2d==7) then   !************open cylinder*****************
    
	
	
        msh%Delta_ll=1d0*pi/Maxedge !2.0d0*pi*5d0/6d0/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        
        msh%xyz(1,0)=cos(0*pi) ; msh%xyz(2,0)=sin(0*pi)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*msh%Delta_ll
            msh%xyz(1,node*2)=cos(0*pi+dx)
            msh%xyz(2,node*2)=sin(0*pi+dx)
        enddo
        !$omp end parallel do
		
		
        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            msh%xyz(1,node*2-1)=(msh%xyz(1,node*2-2)+msh%xyz(1,node*2))/2d0
            msh%xyz(2,node*2-1)=(msh%xyz(2,node*2-2)+msh%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do		
		
		
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo           

    elseif (msh%model2d==8) then   !************corrugated open cylinder*****************
		M = 0.2d0*ker%wavelength
		L = 1.5d0*ker%wavelength
        phi_start = 3d0/2d0*pi
		
		msh%Delta_ll=1d0*pi/Maxedge !2.0d0*pi*5d0/6d0/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        
        msh%xyz(1,0)=cos(phi_start) ; msh%xyz(2,0)=sin(phi_start)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*msh%Delta_ll
            msh%xyz(1,node*2)=(1+M*sin(2*pi*dx/L))*cos(phi_start+dx)
            msh%xyz(2,node*2)=(1+M*sin(2*pi*dx/L))*sin(phi_start+dx)
        enddo
        !$omp end parallel do
		
		
        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            msh%xyz(1,node*2-1)=(msh%xyz(1,node*2-2)+msh%xyz(1,node*2))/2d0
            msh%xyz(2,node*2-1)=(msh%xyz(2,node*2-2)+msh%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do		
		
		
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo    

		
    elseif (msh%model2d==9) then  !*****corrugated corner reflector*****
		M = 0.2d0*ker%wavelength
		L = 1.5d0*ker%wavelength
		
        msh%Delta_ll=2d0/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        

		Am = M*sin(pi/2-2*pi/L)-M		
		msh%xyz(1,0)=-1d0/sqrt(2d0)+Am/sqrt(2d0); msh%xyz(2,0)=1d0/sqrt(2d0)+Am/sqrt(2d0)
		
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge/2
            dx=node*msh%Delta_ll
			Am = M*sin(2*pi*dx/L+pi/2-2*pi/L)-M
			
            msh%xyz(1,node*2)=(dx-1d0)/sqrt(2d0)+Am/sqrt(2d0)
            msh%xyz(2,node*2)=(1d0-dx)/sqrt(2d0)+Am/sqrt(2d0)
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge/2
            dx=node*msh%Delta_ll
			Am = M*sin(2*pi*(dx+1)/L+pi/2-2*pi/L)-M
			
            msh%xyz(1,node*2+Maxedge)=dx/sqrt(2d0)-Am/sqrt(2d0)
            msh%xyz(2,node*2+Maxedge)=dx/sqrt(2d0)+Am/sqrt(2d0)
        enddo
        !$omp end parallel do

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            msh%xyz(1,node*2-1)=(msh%xyz(1,node*2-2)+msh%xyz(1,node*2))/2d0
            msh%xyz(2,node*2-1)=(msh%xyz(2,node*2-2)+msh%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo  		

		msh%Ncorner = 1
		allocate(msh%corner_points(2,msh%Ncorner))
		msh%corner_points(1,1) = 0
		msh%corner_points(2,1) = 0
		
    elseif (msh%model2d==10) then  !*****cup*****

		L1 = 0
		L2 = sqrt(2d0)
		L3 = sqrt(2d0)
		L = 3
        msh%Delta_ll=(2*L1+2*L2+2*L3+L)/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        
		
        !$omp parallel do default(shared) private(node,dx,tt)
        do node=0, Maxedge
            dx=node*msh%Delta_ll
			if(dx<=L1)then
				tt = dx
				msh%xyz(1,node*2)= -(L/2-(L1-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L1+L2+L3)/sqrt(2d0)-tt/sqrt(2d0)			
			else if(dx<=L1+L2)then
				tt = dx - L1
				msh%xyz(1,node*2)= -(L/2-(-L2+L3)/sqrt(2d0))+tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L2+L3)/sqrt(2d0)-tt/sqrt(2d0)
			else if(dx<=L1+L2+L3)then
				tt = dx - (L1+L2)
				msh%xyz(1,node*2)= -(L/2-L3/sqrt(2d0))-tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L3)/sqrt(2d0)-tt/sqrt(2d0)				
			else if(dx<=L1+L2+L3+L)then
				tt = dx - (L1+L2+L3)
				msh%xyz(1,node*2)= -L/2+tt
				msh%xyz(2,node*2)= 0
			else if(dx<=L1+L2+L3+L+L3)then
				tt = dx - (L1+L2+L3+L)
				msh%xyz(1,node*2)= L/2-tt/sqrt(2d0)
				msh%xyz(2,node*2)= tt/sqrt(2d0)		
			else if(dx<=L1+L2+L3+L+L3+L2)then
				tt = dx - (L1+L2+L3+L+L3)
				msh%xyz(1,node*2)= (L/2-L3/sqrt(2d0))+tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L3)/sqrt(2d0)+tt/sqrt(2d0)				
			else
				tt = dx - (L1+L2+L3+L+L3+L2)
				msh%xyz(1,node*2)= (L/2-(-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L3+L2)/sqrt(2d0)+tt/sqrt(2d0)							
			end if 
        enddo
        !$omp end parallel do			

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            msh%xyz(1,node*2-1)=(msh%xyz(1,node*2-2)+msh%xyz(1,node*2))/2d0
            msh%xyz(2,node*2-1)=(msh%xyz(2,node*2-2)+msh%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo  		

		! msh%Ncorner = 4
		! allocate(msh%corner_points(2,msh%Ncorner))
		! msh%corner_points(1,1) = -L/2
		! msh%corner_points(2,1) = 0
		! msh%corner_points(1,2) = L/2
		! msh%corner_points(2,2) = 0
		! msh%corner_points(1,3) = -(L/2-L3/sqrt(2d0))
		! msh%corner_points(2,3) = (L3)/sqrt(2d0)		
		! msh%corner_points(1,4) = (L/2-L3/sqrt(2d0))
		! msh%corner_points(2,4) = (L3)/sqrt(2d0)			

    elseif (msh%model2d==11) then  !*****longer cup*****
		L1 = 1
		L2 = sqrt(2d0)
		L3 = sqrt(2d0)
		L = 3.5
        msh%Delta_ll=(2*L1+2*L2+2*L3+L)/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        
		
        !$omp parallel do default(shared) private(node,dx,tt)
        do node=0, Maxedge
            dx=node*msh%Delta_ll
			if(dx<=L1)then
				tt = dx
				msh%xyz(1,node*2)= -(L/2-(L1-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L1+L2+L3)/sqrt(2d0)-tt/sqrt(2d0)			
			else if(dx<=L1+L2)then
				tt = dx - L1
				msh%xyz(1,node*2)= -(L/2-(-L2+L3)/sqrt(2d0))+tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L2+L3)/sqrt(2d0)-tt/sqrt(2d0)
			else if(dx<=L1+L2+L3)then
				tt = dx - (L1+L2)
				msh%xyz(1,node*2)= -(L/2-L3/sqrt(2d0))-tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L3)/sqrt(2d0)-tt/sqrt(2d0)				
			else if(dx<=L1+L2+L3+L)then
				tt = dx - (L1+L2+L3)
				msh%xyz(1,node*2)= -L/2+tt
				msh%xyz(2,node*2)= 0
			else if(dx<=L1+L2+L3+L+L3)then
				tt = dx - (L1+L2+L3+L)
				msh%xyz(1,node*2)= L/2-tt/sqrt(2d0)
				msh%xyz(2,node*2)= tt/sqrt(2d0)		
			else if(dx<=L1+L2+L3+L+L3+L2)then
				tt = dx - (L1+L2+L3+L+L3)
				msh%xyz(1,node*2)= (L/2-L3/sqrt(2d0))+tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L3)/sqrt(2d0)+tt/sqrt(2d0)				
			else
				tt = dx - (L1+L2+L3+L+L3+L2)
				msh%xyz(1,node*2)= (L/2-(-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L3+L2)/sqrt(2d0)+tt/sqrt(2d0)							
			end if 
        enddo
        !$omp end parallel do			

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            msh%xyz(1,node*2-1)=(msh%xyz(1,node*2-2)+msh%xyz(1,node*2))/2d0
            msh%xyz(2,node*2-1)=(msh%xyz(2,node*2-2)+msh%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo  		

		msh%Ncorner = 6
		allocate(msh%corner_points(2,msh%Ncorner))
		msh%corner_points(1,1) = -(L/2-(-L2+L3)/sqrt(2d0))
		msh%corner_points(2,1) = (L2+L3)/sqrt(2d0)
		msh%corner_points(1,2) = (L/2-(-L2+L3)/sqrt(2d0))
		msh%corner_points(2,2) = (L2+L3)/sqrt(2d0)	
		msh%corner_points(1,3) = -L/2
		msh%corner_points(2,3) = 0
		msh%corner_points(1,4) = L/2
		msh%corner_points(2,4) = 0
		msh%corner_points(1,5) = -(L/2-L3/sqrt(2d0))
		msh%corner_points(2,5) = (L3)/sqrt(2d0)		
		msh%corner_points(1,6) = (L/2-L3/sqrt(2d0))
		msh%corner_points(2,6) = (L3)/sqrt(2d0)		
    endif
    
	msh%maxedgelength = 0
	do edge=1,Maxedge
		msh%maxedgelength = max(msh%maxedgelength,sqrt(sum(abs(msh%xyz(:,msh%info_unk(1,edge))-msh%xyz(:,msh%info_unk(2,edge))))**2))
	end do	
	
	msh%minedgelength = 10000000
	do edge=1,Maxedge
		msh%minedgelength = min(msh%minedgelength,sqrt(sum(abs(msh%xyz(:,msh%info_unk(1,edge))-msh%xyz(:,msh%info_unk(2,edge)))**2)))
	end do		
	
	msh%corner_radius = 2*msh%maxedgelength
	
	! write(*,*)	msh%xyz(1,1:100),sum(msh%xyz(1,:))
	! stop

	
		
    if(ptree%MyID==Main_ID)write (*,*) ''
    if(ptree%MyID==Main_ID)write (*,*) 'Maxedge:',Maxedge
	if(ptree%MyID==Main_ID)write (*,*) 'dx:',msh%Delta_ll/2
	if(ptree%MyID==Main_ID)write (*,*) 'wavelength/maxedgelength:',ker%wavelength/msh%maxedgelength

    
    return
    
end subroutine geo_modeling_CURV


subroutine geo_modeling_SURF(msh,ker,ptree)

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
	type(mesh)::msh
	type(kernelquant)::ker
	type(proctree)::ptree
	
	integer Maxedge
	
    Dimn=3
    
    open(11,file=trim(DATA_DIR)//'/node.geo')
    open(111,file=trim(DATA_DIR)//'/elem.geo')
    
    read(11,*)msh%maxnode
    read(111,*)msh%maxpatch
    Maxedge=msh%maxpatch*3/2
	
    
	
    allocate(msh%xyz(3,msh%maxnode+Maxedge))
    allocate(msh%node_of_patch(0:3,msh%maxpatch),msh%info_unk(0:6,maxedge+1000))
    allocate(msh%normal_of_patch(3,msh%maxpatch))
    
    
    !************msh%xyz****************
    i=1
    do while(i<=msh%maxnode)
        read(11,*)intemp,msh%xyz(1:3,i)
        msh%xyz(1:3,i)=msh%xyz(1:3,i)/msh%scaling
        i=i+1
    enddo
    close(11)
    
    i=1
    if (msh%mesh_normal==1) then
        do while(i<=msh%maxpatch)
            read(111,*)intemp,msh%node_of_patch(1:3,i)
            i=i+1 
        enddo
    elseif (msh%mesh_normal==-1) then
        do while(i<=msh%maxpatch)
            read(111,*)intemp,msh%node_of_patch(3,i),msh%node_of_patch(2,i),msh%node_of_patch(1,i)
            i=i+1 
        enddo
    endif
    close(111)
    
    !************msh%normal_of_patch****************
    
    !$omp parallel do default(shared) private(patch,i,a,b,c,r0)
    do patch=1,msh%maxpatch
        do i=1,3
            a(i)=(msh%xyz(i,msh%node_of_patch(2,patch))-msh%xyz(i,msh%node_of_patch(1,patch)))
            b(i)=(msh%xyz(i,msh%node_of_patch(3,patch))-msh%xyz(i,msh%node_of_patch(1,patch)))
        enddo
        call curl(a,b,c)
        r0=sqrt(c(1)**2+c(2)**2+c(3)**2)
        c(1)=c(1)/r0
        c(2)=c(2)/r0
        c(3)=c(3)/r0
        msh%normal_of_patch(1:3,patch)=c(1:3)	    
    enddo
    !$omp end parallel do
    
    !************msh%info_unk****************

    edge=0
    do i=1,msh%maxpatch-1
        do j=i+1,msh%maxpatch
            flag=0;node1=0;node2=0;iii=1
            do ii=1,3
                do jj=1,3
	     	         if(msh%node_of_patch(ii,i)==msh%node_of_patch(jj,j))then
                        flag=flag+1
                        node_temp(iii)=msh%node_of_patch(ii,i)
                        iii=iii+1
                    endif
                enddo
            enddo
            if(flag==2)then
                edge=edge+1
                if(node_temp(1)<node_temp(2)) then
                    msh%info_unk(1,edge)=node_temp(1)
                    msh%info_unk(2,edge)=node_temp(2)
                else
                    msh%info_unk(1,edge)=node_temp(2)
                    msh%info_unk(2,edge)=node_temp(1)
                endif
                msh%info_unk(3,edge)=i
                msh%info_unk(4,edge)=j       ! notice that : i<j  
                msh%info_unk(0,edge)=0
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
        	            if(msh%node_of_patch(iii,msh%info_unk(jj,edge))==msh%info_unk(jjj,edge)) then
                           node_temp(jj-2)=node_temp(jj-2)+iii               
        	            endif
                 enddo
             enddo
         enddo
         node_temp(1)=6-node_temp(1)
         node_temp(2)=6-node_temp(2)
         msh%info_unk(5,edge)=msh%node_of_patch(node_temp(1),msh%info_unk(3,edge))
         msh%info_unk(6,edge)=msh%node_of_patch(node_temp(2),msh%info_unk(4,edge))
    enddo
    !$omp end parallel do
    
    node=msh%maxnode
    do edge=1, Maxedge
        node=node+1
        msh%info_unk(0,edge)=node
        do i=1,3
            msh%xyz(i,node)=1./2.*(msh%xyz(i,msh%info_unk(1,edge))+msh%xyz(i,msh%info_unk(2,edge)))
        enddo
    enddo
	
	msh%maxedgelength = 0
	do edge=1,Maxedge
		msh%maxedgelength = max(msh%maxedgelength,sqrt(sum(abs(msh%xyz(:,msh%info_unk(1,edge))-msh%xyz(:,msh%info_unk(2,edge)))**2)))
	end do	

	msh%minedgelength = BigValue
	do edge=1,Maxedge
		msh%minedgelength = min(msh%minedgelength,sqrt(sum(abs(msh%xyz(:,msh%info_unk(1,edge))-msh%xyz(:,msh%info_unk(2,edge)))**2)))
	end do	
	
	! write(*,*)	msh%xyz(1,1:100),sum(msh%xyz(1,:))
	! stop
	msh%Nunk = Maxedge

    write (*,*) ''
    write (*,*) 'Maxedge:',Maxedge
	write (*,*) 'minedgelength:',msh%minedgelength
	write (*,*) 'wavelength/minedgelength:',ker%wavelength/msh%minedgelength
	write (*,*) 'maxedgelength:',msh%maxedgelength
	write (*,*) 'wavelength/maxedgelength:',ker%wavelength/msh%maxedgelength

    write (*,*) '' 
    
    return
    
end subroutine geo_modeling_SURF

end module geometry_model