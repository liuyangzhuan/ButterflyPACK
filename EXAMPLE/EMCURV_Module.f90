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

! Developers: Yang Liu, Xiaoye S. Li.
!             (Lawrence Berkeley National Lab, Computational Research Division).

module EMCURV_MODULE
use z_BPACK_DEFS
use z_misc
implicit none

	!**** define your application-related variables here   
	
	!**** quantities related to geometries, meshes, unknowns and points 

	type quant_EMCURV
		real(kind=8) wavenum    ! CEM: wave number  
		real(kind=8) wavelength  ! CEM: wave length
		real(kind=8) omiga       ! CEM: angular frequency
		real(kind=8) rank_approximate_para1, rank_approximate_para2, rank_approximate_para3 ! CEM: rank estimation parameter  
		integer RCS_static  ! CEM: 1: monostatic or 2: bistatic RCS
		integer RCS_Nsample ! CEM: number of RCS samples
		real(kind=8):: CFIE_alpha ! CEM: combination parameter in CFIE
		
		integer Nunk ! size of the matrix 
		real(kind=8),allocatable:: xyz(:,:)   ! coordinates of the points
		integer,allocatable:: info_unk(:,:)  
		! for 2D mesh: 0 point to coordinates of each edge center (unknown x), 1-2 point to coordinates of each edge vertice  

		! 2D mesh
		integer maxnode ! # of vertices in a mesh
		integer maxedge ! # of edges in a mesh 
		real(kind=8) maxedgelength,minedgelength ! maximum and minimum edge length for 2D and 3D meshes	
		integer model2d ! # shape of 2-D curves: (1=strip; 2=corner reflector; 3=two opposite strips; 4=CR with RRS; 5=cylinder; 6=Rectangle Cavity); 7=half cylinder; 8=corrugated half cylinder; 9=corrugated corner reflector; 10=open polygon; 11=taller open polygon
		real(kind=8) Delta_ll	! edge length of each element in 2D curves	
		integer::Ncorner ! # of corners in 2D curves to facilate partitioning
		real(kind=8),allocatable::corner_points(:,:) ! coordinates of corner points 
		real(kind=8)::corner_radius ! radius of the corner points within which no partitioning is performed 
	end type quant_EMCURV

contains

	subroutine delete_quant_EMCURV(quant)
		implicit none
		type(quant_EMCURV):: quant
		if(allocated(quant%xyz))deallocate(quant%xyz)
		if(allocated(quant%info_unk))deallocate(quant%info_unk)
		if(allocated(quant%corner_points))deallocate(quant%corner_points)
	end subroutine delete_quant_EMCURV


	!**** user-defined subroutine to sample Z_mn
	subroutine Zelem_EMCURV(m,n,value_e,quant)
		
		use z_BPACK_DEFS
		implicit none
		
		integer edge_m, edge_n, i, j, flag
		integer, INTENT(IN):: m,n
		real(kind=8) r_mn, rtemp1, rtemp2
		complex(kind=8) value_e

		class(*),pointer :: quant
		
		select TYPE(quant)
			type is (quant_EMCURV)
				! convert to new indices because msh%info_unk has been reordered

				edge_m = m
				edge_n = n

				! if(mod(m,2)==0)edge_m=m-1

				
				if (edge_m/=edge_n) then
				
					flag=0
					do j=1,2
						do i=1,2
							if (quant%info_unk(i,edge_m)==quant%info_unk(j,edge_n)) then
								flag=1
							endif
						enddo
					enddo
					
					if (flag==1) then
						
						value_e=quant%wavenum*impedence0/4.0*quant%Delta_ll*(1-junit/pi*(3*LOG(3*gamma*quant%wavenum*quant%Delta_ll/4.0)-LOG(gamma*quant%wavenum*quant%Delta_ll/4.0)-2))
					
					else
				
						r_mn=(quant%xyz(1,quant%info_unk(0,edge_m))-quant%xyz(1,quant%info_unk(0,edge_n)))**2+(quant%xyz(2,quant%info_unk(0,edge_m))-quant%xyz(2,quant%info_unk(0,edge_n)))**2
						r_mn=sqrt(r_mn)
						value_e=quant%wavenum*impedence0/4.0*quant%Delta_ll*z_Hankel02_Func(quant%wavenum*r_mn)
					endif
				else
					value_e=quant%wavenum*impedence0/4.0*quant%Delta_ll*(1.0-junit*2.0/pi*(LOG(gamma*quant%wavenum*quant%Delta_ll/4.0)-1.0))
				endif

			class default
				write(*,*)"unexpected type"
				stop
			end select			
			
		return
		
	end subroutine Zelem_EMCURV	
	

	subroutine VV_polar_CURV(dphi,edge,ctemp_1,curr,msh,quant)
		
		use z_BPACK_DEFS
		! use EMCURV_MODULE
		implicit none
		complex(kind=8)::curr
		complex(kind=8) ctemp,phase,ctemp_1
		real(kind=8) dsita,dphi    
		integer edge
		type(z_mesh)::msh
		type(quant_EMCURV)::quant
		
		phase=junit*quant%wavenum*(quant%xyz(1,quant%info_unk(0,edge))*cos(dphi*pi/180.)+quant%xyz(2,quant%info_unk(0,edge))*sin(dphi*pi/180.))
		ctemp_1=curr*quant%Delta_ll*exp(phase)

		return
		
	end subroutine VV_polar_CURV

	subroutine RCS_bistatic_CURV(curr,msh,quant,ptree)
		!integer flag
		use z_BPACK_DEFS
		! use EMCURV_MODULE
		implicit none
		complex(kind=8)::curr(:)
		real(kind=8) rcs
		complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1,ctemp_loc
		real(kind=8) ddphi,dphi
		
		integer i,j,ii,jj,iii,jjj,patch,flag
		real(kind=8) l_edge,l_edgefine
		type(z_mesh)::msh
		integer edge,edge_m,edge_n,ierr
		type(quant_EMCURV)::quant 
		type(z_proctree)::ptree
		
		if(ptree%MyID==Main_ID)open (100, file='VV_bistatic.txt')
	  
		ddphi=180./quant%RCS_Nsample
		
		do i=0, quant%RCS_Nsample   !phi=0
			dphi=i*ddphi
			ctemp_loc=0
			rcs=0
			!$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)   
			do edge=msh%idxs,msh%idxe
				call VV_polar_CURV(dphi,msh%new2old(edge),ctemp_1,curr(edge-msh%idxs+1),msh,quant)
				ctemp_loc=ctemp_loc+ctemp_1    
			enddo
			!$omp end parallel do
			
			call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
			
			rcs=(abs(impedence0*ctemp))**2/4d0*quant%wavenum
			!rcs=rcs/quant%wavelength
			rcs=10*log10(rcs)
			if(ptree%MyID==Main_ID)write(100,*)dphi,rcs
		enddo
		if(ptree%MyID==Main_ID) close(100)
		
		return
		
	end subroutine RCS_bistatic_CURV



	subroutine RCS_monostatic_VV_CURV(dphi,rcs,curr,msh,quant,ptree)

		use z_BPACK_DEFS
		! use EMCURV_MODULE
		implicit none
		
		real(kind=8) rcs
		complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1,ctemp_2
		real(kind=8) dsita,dphi
		integer edge,edge_m,edge_n,ierr
		complex(kind=8):: curr(:)
		type(z_mesh)::msh
		type(quant_EMCURV)::quant
		type(z_proctree)::ptree
		
		ctemp_loc=0
			rcs=0
			!$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
			do edge=msh%idxs,msh%idxe
				call VV_polar_CURV(dphi,msh%new2old(edge),ctemp_1,curr(edge-msh%idxs+1),msh,quant)
				ctemp_loc=ctemp_loc+ctemp_1
			enddo
			!$omp end parallel do
			
			call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
			
			rcs=(abs(impedence0*ctemp))**2/4d0*quant%wavenum
			!rcs=rcs/quant%wavelength
			rcs=10*log10(rcs)
			
		return
		
	end subroutine RCS_monostatic_VV_CURV

	subroutine element_Vinc_VV_CURV(phi,edge,value,msh,quant)

		use z_BPACK_DEFS
		! use EMCURV_MODULE
		implicit none
		
		integer edge
		complex(kind=8) value
		real(kind=8) theta, phi
		complex(kind=8)  phase
		type(z_mesh)::msh
		type(quant_EMCURV)::quant
		
		phase=junit*quant%wavenum*(quant%xyz(1,quant%info_unk(0,edge))*cos(phi*pi/180.)+quant%xyz(2,quant%info_unk(0,edge))*sin(phi*pi/180.))
		value=exp(phase)
		
		return
		
	end subroutine element_Vinc_VV_CURV	
	
	

subroutine geo_modeling_CURV(quant,MPIcomm)

    use z_BPACK_DEFS
    implicit none
    ! type(z_mesh)::msh
	type(quant_EMCURV)::quant
    integer i,j,ii,jj,iii,jjj
    integer intemp
    integer node, patch, edge, flag
    integer node1, node2, num_node,Maxedge
    integer node_temp(2)
    real(kind=8) dx, xx, yy, rr, theta,L,M,Am,tt,L1,L2,L3
    
	
    real(kind=8),allocatable :: node_xy_original(:,:)
    integer,allocatable :: num_edge_of_node(:)
    
    real(kind=8) a(3),b(3),c(3),r0, phi_start
	! type(z_proctree)::ptree
	integer MPIcomm,ierr,MyID
	
	call MPI_Comm_rank(MPIcomm,MyID,ierr)
	
	Maxedge	= quant%Nunk

	quant%Ncorner = 0
	
    if (quant%model2d==1) then  !*****single strip*****
        
        quant%Delta_ll=2d0/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        quant%xyz(1,0)=-1d0 ; quant%xyz(2,0)=0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, quant%maxnode-1
            dx=node*quant%Delta_ll/2
            quant%xyz(1,node)=dx-1d0
            quant%xyz(2,node)=0
        enddo
        !$omp end parallel do
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1   ! 1 start 0 center 2 end
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo           
        
    elseif (quant%model2d==2) then  !*****corner reflector*****
    
        quant%Delta_ll=2d0/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        
        quant%xyz(1,0)=-1d0/sqrt(2d0) ; quant%xyz(2,0)=1d0/sqrt(2d0)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*quant%Delta_ll/2
            quant%xyz(1,node)=(dx-1d0)/sqrt(2d0)
            quant%xyz(2,node)=(1d0-dx)/sqrt(2d0)
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*quant%Delta_ll/2
            quant%xyz(1,node+Maxedge)=dx/sqrt(2d0)
            quant%xyz(2,node+Maxedge)=dx/sqrt(2d0)
        enddo
        !$omp end parallel do
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo
		
		quant%Ncorner = 1
		allocate(quant%corner_points(2,quant%Ncorner))
		quant%corner_points(1,1) = 0
		quant%corner_points(2,1) = 0        
    
    elseif (quant%model2d==3) then  !*****two opposite strips*****
    
        quant%Delta_ll=2d0/Maxedge
        quant%maxnode=2*Maxedge+2
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        
        quant%xyz(1,0)=-0.5d0 ; quant%xyz(2,0)=1.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*quant%Delta_ll/2
            quant%xyz(1,node)=dx-0.5d0
            quant%xyz(2,node)=1.0d0
        enddo
        !$omp end parallel do
        quant%xyz(1,Maxedge+1)=-0.5d0 ; quant%xyz(2,Maxedge+1)=-1.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*quant%Delta_ll/2
            quant%xyz(1,node+Maxedge+1)=dx-0.5d0
            quant%xyz(2,node+Maxedge+1)=-1.0d0
        enddo
        !$omp end parallel do
        
        intemp=int(Maxedge/2)
        
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, intemp
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo
        
        quant%info_unk(1,intemp+1)=quant%maxnode/2 ; quant%info_unk(2,intemp+1)=quant%maxnode/2+2 ; quant%info_unk(0,intemp+1)=quant%maxnode/2+1
        do edge=intemp+2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo
    
    
    elseif (quant%model2d==5) then !************cylinder*****************
        
        quant%Delta_ll=  2.0d0*pi/Maxedge
        quant%maxnode=2*Maxedge
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        quant%xyz(1,0)=1.0d0 ; quant%xyz(2,0)=0.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, quant%maxnode-1
            dx=node*quant%Delta_ll/2
            quant%xyz(1,node)=cos(dx)
            quant%xyz(2,node)=sin(dx)
        enddo
        !$omp end parallel do
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge-1
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo     
        quant%info_unk(1,Maxedge)=quant%info_unk(2,Maxedge-1)
        quant%info_unk(2,Maxedge)=0
        quant%info_unk(0,Maxedge)=quant%info_unk(0,Maxedge-1)+2      
        
    elseif (quant%model2d==6) then !***********cavity******************
    
        quant%Delta_ll=11d0/3d0/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        dx=quant%Delta_ll/2d0
        quant%xyz(1,0)=-1d0/6. ; quant%xyz(2,0)=1.0d0
        node=0
        do while (flag==0) 
            node=node+1
            quant%xyz(1,node)=quant%xyz(1,node-1)-dx
            quant%xyz(2,node)=quant%xyz(2,node-1)
            if (quant%xyz(1,node)<=-0.5d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            quant%xyz(1,node)=quant%xyz(1,node-1)
            quant%xyz(2,node)=quant%xyz(2,node-1)-dx
            if (quant%xyz(2,node)<=0.) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            quant%xyz(1,node)=quant%xyz(1,node-1)+dx
            quant%xyz(2,node)=quant%xyz(2,node-1)
            if (quant%xyz(1,node)>=0.5d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            quant%xyz(1,node)=quant%xyz(1,node-1)
            quant%xyz(2,node)=quant%xyz(2,node-1)+dx
            if (quant%xyz(2,node)>=1.0d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            quant%xyz(1,node)=quant%xyz(1,node-1)-dx
            quant%xyz(2,node)=quant%xyz(2,node-1)
            if (node>=quant%maxnode-1) then
                flag=1
            endif
        enddo
        
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo

		quant%Ncorner = 4
		allocate(quant%corner_points(2,quant%Ncorner))
		quant%corner_points(1,1) = -0.5
		quant%corner_points(2,1) = 0     
		quant%corner_points(1,2) = 0.5
		quant%corner_points(2,2) = 0
		quant%corner_points(1,3) = 0.5
		quant%corner_points(2,3) = 1
		quant%corner_points(1,4) = -0.5
		quant%corner_points(2,4) = 1		
		
    elseif (quant%model2d==7) then   !************open cylinder*****************
    
	
	
        quant%Delta_ll=1d0*pi/Maxedge !2.0d0*pi*5d0/6d0/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        
        quant%xyz(1,0)=cos(0*pi) ; quant%xyz(2,0)=sin(0*pi)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*quant%Delta_ll
            quant%xyz(1,node*2)=cos(0*pi+dx)
            quant%xyz(2,node*2)=sin(0*pi+dx)
        enddo
        !$omp end parallel do
		
		
        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            quant%xyz(1,node*2-1)=(quant%xyz(1,node*2-2)+quant%xyz(1,node*2))/2d0
            quant%xyz(2,node*2-1)=(quant%xyz(2,node*2-2)+quant%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do		
		
		
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo           

    elseif (quant%model2d==8) then   !************corrugated open cylinder*****************
		M = 0.2d0*quant%wavelength
		L = 1.5d0*quant%wavelength
        phi_start = 3d0/2d0*pi
		
		quant%Delta_ll=1d0*pi/Maxedge !2.0d0*pi*5d0/6d0/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        
        quant%xyz(1,0)=cos(phi_start) ; quant%xyz(2,0)=sin(phi_start)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*quant%Delta_ll
            quant%xyz(1,node*2)=(1+M*sin(2*pi*dx/L))*cos(phi_start+dx)
            quant%xyz(2,node*2)=(1+M*sin(2*pi*dx/L))*sin(phi_start+dx)
        enddo
        !$omp end parallel do
		
		
        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            quant%xyz(1,node*2-1)=(quant%xyz(1,node*2-2)+quant%xyz(1,node*2))/2d0
            quant%xyz(2,node*2-1)=(quant%xyz(2,node*2-2)+quant%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do		
		
		
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo    

		
    elseif (quant%model2d==9) then  !*****corrugated corner reflector*****
		M = 0.2d0*quant%wavelength
		L = 1.5d0*quant%wavelength
		
        quant%Delta_ll=2d0/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        

		Am = M*sin(pi/2-2*pi/L)-M		
		quant%xyz(1,0)=-1d0/sqrt(2d0)+Am/sqrt(2d0); quant%xyz(2,0)=1d0/sqrt(2d0)+Am/sqrt(2d0)
		
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge/2
            dx=node*quant%Delta_ll
			Am = M*sin(2*pi*dx/L+pi/2-2*pi/L)-M
			
            quant%xyz(1,node*2)=(dx-1d0)/sqrt(2d0)+Am/sqrt(2d0)
            quant%xyz(2,node*2)=(1d0-dx)/sqrt(2d0)+Am/sqrt(2d0)
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge/2
            dx=node*quant%Delta_ll
			Am = M*sin(2*pi*(dx+1)/L+pi/2-2*pi/L)-M
			
            quant%xyz(1,node*2+Maxedge)=dx/sqrt(2d0)-Am/sqrt(2d0)
            quant%xyz(2,node*2+Maxedge)=dx/sqrt(2d0)+Am/sqrt(2d0)
        enddo
        !$omp end parallel do

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            quant%xyz(1,node*2-1)=(quant%xyz(1,node*2-2)+quant%xyz(1,node*2))/2d0
            quant%xyz(2,node*2-1)=(quant%xyz(2,node*2-2)+quant%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo  		

		quant%Ncorner = 1
		allocate(quant%corner_points(2,quant%Ncorner))
		quant%corner_points(1,1) = 0
		quant%corner_points(2,1) = 0
		
    elseif (quant%model2d==10) then  !*****cup*****

		L1 = 0
		L2 = sqrt(2d0)
		L3 = sqrt(2d0)
		L = 3
        quant%Delta_ll=(2*L1+2*L2+2*L3+L)/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        
		
        !$omp parallel do default(shared) private(node,dx,tt)
        do node=0, Maxedge
            dx=node*quant%Delta_ll
			if(dx<=L1)then
				tt = dx
				quant%xyz(1,node*2)= -(L/2-(L1-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L1+L2+L3)/sqrt(2d0)-tt/sqrt(2d0)			
			else if(dx<=L1+L2)then
				tt = dx - L1
				quant%xyz(1,node*2)= -(L/2-(-L2+L3)/sqrt(2d0))+tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L2+L3)/sqrt(2d0)-tt/sqrt(2d0)
			else if(dx<=L1+L2+L3)then
				tt = dx - (L1+L2)
				quant%xyz(1,node*2)= -(L/2-L3/sqrt(2d0))-tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L3)/sqrt(2d0)-tt/sqrt(2d0)				
			else if(dx<=L1+L2+L3+L)then
				tt = dx - (L1+L2+L3)
				quant%xyz(1,node*2)= -L/2+tt
				quant%xyz(2,node*2)= 0
			else if(dx<=L1+L2+L3+L+L3)then
				tt = dx - (L1+L2+L3+L)
				quant%xyz(1,node*2)= L/2-tt/sqrt(2d0)
				quant%xyz(2,node*2)= tt/sqrt(2d0)		
			else if(dx<=L1+L2+L3+L+L3+L2)then
				tt = dx - (L1+L2+L3+L+L3)
				quant%xyz(1,node*2)= (L/2-L3/sqrt(2d0))+tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L3)/sqrt(2d0)+tt/sqrt(2d0)				
			else
				tt = dx - (L1+L2+L3+L+L3+L2)
				quant%xyz(1,node*2)= (L/2-(-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L3+L2)/sqrt(2d0)+tt/sqrt(2d0)							
			end if 
        enddo
        !$omp end parallel do			

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            quant%xyz(1,node*2-1)=(quant%xyz(1,node*2-2)+quant%xyz(1,node*2))/2d0
            quant%xyz(2,node*2-1)=(quant%xyz(2,node*2-2)+quant%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo  		

		! quant%Ncorner = 4
		! allocate(quant%corner_points(2,quant%Ncorner))
		! quant%corner_points(1,1) = -L/2
		! quant%corner_points(2,1) = 0
		! quant%corner_points(1,2) = L/2
		! quant%corner_points(2,2) = 0
		! quant%corner_points(1,3) = -(L/2-L3/sqrt(2d0))
		! quant%corner_points(2,3) = (L3)/sqrt(2d0)		
		! quant%corner_points(1,4) = (L/2-L3/sqrt(2d0))
		! quant%corner_points(2,4) = (L3)/sqrt(2d0)			

    elseif (quant%model2d==11) then  !*****longer cup*****
		L1 = 1
		L2 = sqrt(2d0)
		L3 = sqrt(2d0)
		L = 3.5
        quant%Delta_ll=(2*L1+2*L2+2*L3+L)/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        
		
        !$omp parallel do default(shared) private(node,dx,tt)
        do node=0, Maxedge
            dx=node*quant%Delta_ll
			if(dx<=L1)then
				tt = dx
				quant%xyz(1,node*2)= -(L/2-(L1-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L1+L2+L3)/sqrt(2d0)-tt/sqrt(2d0)			
			else if(dx<=L1+L2)then
				tt = dx - L1
				quant%xyz(1,node*2)= -(L/2-(-L2+L3)/sqrt(2d0))+tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L2+L3)/sqrt(2d0)-tt/sqrt(2d0)
			else if(dx<=L1+L2+L3)then
				tt = dx - (L1+L2)
				quant%xyz(1,node*2)= -(L/2-L3/sqrt(2d0))-tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L3)/sqrt(2d0)-tt/sqrt(2d0)				
			else if(dx<=L1+L2+L3+L)then
				tt = dx - (L1+L2+L3)
				quant%xyz(1,node*2)= -L/2+tt
				quant%xyz(2,node*2)= 0
			else if(dx<=L1+L2+L3+L+L3)then
				tt = dx - (L1+L2+L3+L)
				quant%xyz(1,node*2)= L/2-tt/sqrt(2d0)
				quant%xyz(2,node*2)= tt/sqrt(2d0)		
			else if(dx<=L1+L2+L3+L+L3+L2)then
				tt = dx - (L1+L2+L3+L+L3)
				quant%xyz(1,node*2)= (L/2-L3/sqrt(2d0))+tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L3)/sqrt(2d0)+tt/sqrt(2d0)				
			else
				tt = dx - (L1+L2+L3+L+L3+L2)
				quant%xyz(1,node*2)= (L/2-(-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L3+L2)/sqrt(2d0)+tt/sqrt(2d0)							
			end if 
        enddo
        !$omp end parallel do			

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            quant%xyz(1,node*2-1)=(quant%xyz(1,node*2-2)+quant%xyz(1,node*2))/2d0
            quant%xyz(2,node*2-1)=(quant%xyz(2,node*2-2)+quant%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo  		

		quant%Ncorner = 6
		allocate(quant%corner_points(2,quant%Ncorner))
		quant%corner_points(1,1) = -(L/2-(-L2+L3)/sqrt(2d0))
		quant%corner_points(2,1) = (L2+L3)/sqrt(2d0)
		quant%corner_points(1,2) = (L/2-(-L2+L3)/sqrt(2d0))
		quant%corner_points(2,2) = (L2+L3)/sqrt(2d0)	
		quant%corner_points(1,3) = -L/2
		quant%corner_points(2,3) = 0
		quant%corner_points(1,4) = L/2
		quant%corner_points(2,4) = 0
		quant%corner_points(1,5) = -(L/2-L3/sqrt(2d0))
		quant%corner_points(2,5) = (L3)/sqrt(2d0)		
		quant%corner_points(1,6) = (L/2-L3/sqrt(2d0))
		quant%corner_points(2,6) = (L3)/sqrt(2d0)		
    endif
    
	quant%maxedgelength = 0
	do edge=1,Maxedge
		quant%maxedgelength = max(quant%maxedgelength,sqrt(sum(abs(quant%xyz(:,quant%info_unk(1,edge))-quant%xyz(:,quant%info_unk(2,edge))))**2))
	end do	
	
	quant%minedgelength = 10000000
	do edge=1,Maxedge
		quant%minedgelength = min(quant%minedgelength,sqrt(sum(abs(quant%xyz(:,quant%info_unk(1,edge))-quant%xyz(:,quant%info_unk(2,edge)))**2)))
	end do		
	
	quant%corner_radius = 2*quant%maxedgelength
	
	! write(*,*)	quant%xyz(1,1:100),sum(quant%xyz(1,:))
	! stop
		
    if(MyID==Main_ID)write (*,*) ''
    if(MyID==Main_ID)write (*,*) 'Maxedge:',Maxedge
	if(MyID==Main_ID)write (*,*) 'dx:',quant%Delta_ll/2
	if(MyID==Main_ID)write (*,*) 'wavelength/maxedgelength:',quant%wavelength/quant%maxedgelength

    
    return
    
end subroutine geo_modeling_CURV
	
	
	
subroutine EM_solve_CURV(bmat,option,msh,quant,ptree,stats)
    
    use z_BPACK_DEFS
	use z_BPACK_Solve_Mul	
    
    implicit none
    
    integer i, j, ii, jj, iii, jjj, ierr
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, n_iter_max, iter, N_unk, N_unk_loc
    real(kind=8) theta, phi, dphi, rcs_V, rcs_H
    real T0
    real(kind=8) n1,n2,rtemp	
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:)
	real(kind=8):: rel_error
	type(z_Hoption)::option
	type(z_mesh)::msh
	type(quant_EMCURV)::quant
	type(z_proctree)::ptree
	class(*)::bmat
	type(z_Hstat)::stats	
	complex(kind=8),allocatable:: current(:),voltage(:)
	
	N_unk_loc = msh%idxe-msh%idxs+1
	
	if(option%ErrSol==1)then
		call z_bpack_test_solve_error(bmat,N_unk_loc,option,ptree,stats)
	endif
	
    if (quant%RCS_static==2) then
    
        phi=180d0
        
        allocate (current(N_unk_loc))
		Current=0
        allocate (voltage(N_unk_loc))
								  
        !$omp parallel do default(shared) private(edge,value_Z)
        do edge=msh%idxs, msh%idxe
            call element_Vinc_VV_CURV(phi,msh%new2old(edge),value_Z,msh,quant)
            voltage(edge-msh%idxs+1)=value_Z
        enddo    
        !$omp end parallel do
        
        n1 = OMP_get_wtime()
        
		call z_bpack_solution(bmat,Current,Voltage,N_unk_loc,1,option,ptree,stats)
		
		n2 = OMP_get_wtime()
		stats%Time_Sol = stats%Time_Sol + n2-n1
		call MPI_ALLREDUCE(stats%Time_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
		if(ptree%MyID==Main_ID)then
			write (*,*) ''
			write (*,*) 'Solving:',rtemp,'Seconds'
			write (*,*) ''
		endif
		call MPI_ALLREDUCE(stats%Flop_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,'(A13Es14.2)') 'Solve flops:',rtemp	
		
		T0=secnds(0.0)
        call RCS_bistatic_CURV(Current,msh,quant,ptree)

		if(ptree%MyID==Main_ID)then
			write (*,*) ''
			write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
			write (*,*) ''
		endif
    
		deallocate(current)
		deallocate(voltage)
	
    elseif (quant%RCS_static==1) then
    
        allocate (current(N_unk_loc))
        num_sample=quant%RCS_Nsample
        dphi=180./num_sample
        
		allocate (b(N_unk_loc,num_sample+1))
		allocate (x(N_unk_loc,num_sample+1))
		x=0
		
        if(ptree%MyID==Main_ID)open (100, file='RCS_monostatic.txt')

        n1=OMP_get_wtime()       
        do j=0, num_sample 
            phi=j*dphi
			!$omp parallel do default(shared) private(edge,value_Z)
			do edge=msh%idxs, msh%idxe
				call element_Vinc_VV_CURV(phi,msh%new2old(edge),value_Z,msh,quant)
				b(edge-msh%idxs+1,j+1)=value_Z
			enddo    
			!$omp end parallel do
		enddo
		
		call z_bpack_solution(bmat,x,b,N_unk_loc,num_sample+1,option,ptree,stats)
			
			
		do j=0, num_sample 	
			phi=j*dphi
			Current=x(:,j+1)
            call RCS_monostatic_VV_CURV(phi,rcs_V,Current,msh,quant,ptree)
!             !$omp parallel do default(shared) private(i)
!             do i=1, N_unk
!                 current(i)=vectors_block(0)%vector(i,2)
!             enddo
!             !$omp end parallel do
!             call RCS_monostatic_HH(theta,phi,rcs_H)
            
            if(ptree%MyID==Main_ID)write (100,*) j,phi,rcs_V !,rcs_H
            
            ! deallocate (vectors_block)
            
        enddo
        
		n2 = OMP_get_wtime()
		stats%Time_Sol = stats%Time_Sol + n2-n1
		call MPI_ALLREDUCE(stats%Time_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)		
		
        if(ptree%MyID==Main_ID .and. option%verbosity>=0)then
			close(100)
			write (*,*) ''
			write (*,*) 'Solving:',rtemp,'Seconds'
			write (*,*) ''     
		endif
		call MPI_ALLREDUCE(stats%Flop_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,'(A13Es14.2)') 'Solve flops:',rtemp	
		
	! call MPI_ALLREDUCE(stats%Time_RedistV,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) '     Time_RedistV:', rtemp			
		
		deallocate(b)
		deallocate(x)
		deallocate(Current)
		
    endif
        
    return
    
end subroutine EM_solve_CURV	
	
	
	
	
	

!**** C interface of initialization of quant_EMCURV
	!Npo: (input) matrix size
	!Locations: (output) coordinates used for clustering  
	!quant_emcurv_Cptr: (output) C-handle of quant_EMCURV
	!ptree: (input) processor tree
	!model2d: (input) choice of 2d geometries, (1=strip; 2=corner reflector; 3=two opposite strips; 4=CR with RRS; 5=cylinder; 6=Rectangle Cavity); 7=half cylinder; 8=corrugated half cylinder; 9=corrugated corner reflector; 10=open polygon; 11=taller open polygon)
	!wavelength: (input) wave-length in free space
	!MPIcomm: user-provided MPI communicator
subroutine C_EMCURV_Init(Npo,Locations,quant_emcurv_Cptr, model2d, wavelength, MPIcomm) bind(c, name="c_emcurv_init")	
	implicit none 
	integer Npo,Dimn
	real(kind=8) Locations(*)
	type(c_ptr) :: quant_emcurv_Cptr
	! type(c_ptr) :: ptree_Cptr
	! type(z_proctree),pointer::ptree	
	type(quant_EMCURV),pointer::quant
	integer model2d
	real(kind=8) wavelength
	integer MPIcomm
	
	real(kind=8),parameter :: cd = 299792458d0
	integer seed_myid(50)
	integer times(8)
	integer edge
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	character(len=1024)  :: strings
	integer MyID,ierr
	
	call MPI_Comm_rank(MPIcomm,MyID,ierr)	
	
	!**** allocate EMCURV structures 
	allocate(quant)
	! call c_f_pointer(ptree_Cptr, ptree)


	Dimn=2
	quant%Nunk = Npo
	quant%model2d = model2d
	quant%wavelength=wavelength
	quant%RCS_static=1
    quant%RCS_Nsample=2000		
    quant%omiga=2*pi/quant%wavelength/sqrt(mu0*eps0)
    quant%wavenum=2*pi/quant%wavelength	
	! quant%rank_approximate_para1=6.0
    ! quant%rank_approximate_para2=6.0
    ! quant%rank_approximate_para3=6.0
	
	
   !***********************************************************************
   if(MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'wavelength:',quant%wavelength
   write (*,*) ''
   endif
   !***********************************************************************
	
	t1 = OMP_get_wtime()
    if(MyID==Main_ID)write(*,*) "geometry modeling......"
    call geo_modeling_CURV(quant,MPIcomm)
	
	! generate the list of points for clustering
	do edge=1, quant%Nunk
		Locations((edge-1)*Dimn+1:edge*Dimn) = quant%xyz(:,edge*2-1)
	enddo	
   
	if(MyID==Main_ID)write(*,*) "modeling finished"
    if(MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1	
	
	!**** return the C address of hodlr structures to C caller
	quant_emcurv_Cptr=c_loc(quant)
	
end subroutine C_EMCURV_Init



!**** C interface of computing mn^th entry of the EMCURV matrix
	!m: (input) matrix row index
	!n: (input) matrix row index  
	!quant_emcurv_Cptr: (output) C-handle of quant_EMCURV
	!value: (output) Z_mn
subroutine C_EMCURV_Sample(m,n,value,quant_emcurv_Cptr) bind(c, name="c_emcurv_sample")	
	implicit none 
	integer m,n
	complex(kind=8) value
	type(c_ptr),value :: quant_emcurv_Cptr
	type(quant_EMCURV),pointer::quant
	class(*),pointer :: QuantApp

	call c_f_pointer(quant_emcurv_Cptr, quant)
	QuantApp=>quant
	!!!!!  m,n changed to m+1,n+1 because Zelem_EMCURV assumes indices staring from 1
	call Zelem_EMCURV(m+1,n+1,value,QuantApp)
end subroutine C_EMCURV_Sample

	

end module EMCURV_MODULE	