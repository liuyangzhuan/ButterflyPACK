module RCS_Bi
use gauss_rule
use analytic_part
use current_mapping
contains 

subroutine RCS_bistatic_SURF(curr,msh,ker)
    !integer flag
	use MODULE_FILE
    implicit none
    
    real*8 rcs
    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1
    real*8 theta,phi,dtheta,dphi
    
    integer i,j,ii,jj,iii,jjj,patch,flag
    real*8 l_edge,l_edgefine
    real*8 a(3)
    real*8,allocatable:: x(:),y(:),z(:),w(:)
    complex(kind=8):: curr(:,:)
	type(mesh)::msh
	type(kernelquant)::ker
	
    integer edge,edge_m,edge_n
    
	allocate(x(msh%integral_points))
	allocate(y(msh%integral_points))
	allocate(z(msh%integral_points))
	allocate(w(msh%integral_points))
	
    open (100, file='VV_bistatic.txt')	
    theta=90.
    dphi=180./ker%RCS_Nsample

    do i=0, ker%RCS_Nsample
        phi=i*dphi 
        ctemp=(0.,0.)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp)   
        do edge=1,msh%Nunk
            call VV_polar_SURF(theta,phi,edge,ctemp_1,curr(:,1),msh,ker)
            ctemp=ctemp+ctemp_1    
        enddo
        !$omp end parallel do
        rcs=(abs(ker%wavenum*ctemp))**2/4/pi
        !rcs=rcs/ker%wavelength
        rcs=10*log10(rcs)
        write(100,*)phi,rcs
    enddo
    close(100)

	call current_node_patch_mapping('V',curr(:,1),msh)    	
	
	
   
    open (1000, file='HH_bistatic.txt')
	   
    do i=0, ker%RCS_Nsample
        phi=i*dphi
        ctemp=(0,0)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp)   
        do edge=1,msh%Nunk
            call HH_polar_SURF(theta,phi,edge,ctemp_1,curr(:,2),msh,ker)
            ctemp=ctemp+ctemp_1    
        enddo
        !$omp end parallel do
        rcs=(abs(ker%wavenum*ctemp))**2/4/pi
        !rcs=rcs/ker%wavelength
        rcs=10*log10(rcs)
        write(1000,*)phi,rcs
    enddo
    close(1000)
    call current_node_patch_mapping('H',curr(:,2),msh)    
	
	deallocate(x,y,z,w)
	
    return
    
end subroutine RCS_bistatic_SURF

subroutine VV_polar_SURF(theta,phi,edge,ctemp_1,curr,msh,ker)
    
    use MODULE_FILE
    implicit none
    
    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1
    real*8 theta,phi
    type(kernelquant)::ker
    integer i,j,ii,jj,iii,jjj,patch,flag
    real*8 l_edge,l_edgefine
    real*8 a(3)
    complex(kind=8)::curr(:)
    integer edge,edge_m,edge_n
    type(mesh)::msh
	real*8,allocatable:: x(:),y(:),z(:),w(:)
		
	allocate(x(msh%integral_points))
	allocate(y(msh%integral_points))
	allocate(z(msh%integral_points))
	allocate(w(msh%integral_points))		
		
            ctemp_rcs(1:3)=(0,0)          
            l_edge=sqrt((msh%xyz(1,msh%info_unk(1,edge))-msh%xyz(1,msh%info_unk(2,edge)))**2+(msh%xyz(2,msh%info_unk(1,edge))-msh%xyz(2,msh%info_unk(2,edge)))**2+(msh%xyz(3,msh%info_unk(1,edge))-msh%xyz(3,msh%info_unk(2,edge)))**2)
			do jj=3,4  ! two triganle
	            !patch=msh%info_unk(jj,edgefine_m)
                call gau_grobal(edge,jj,x,y,z,w,msh)

		        do ii=1,msh%integral_points
		           
				    a(1)=x(ii)-msh%xyz(1,msh%info_unk(jj+2,edge))  
		            a(2)=y(ii)-msh%xyz(2,msh%info_unk(jj+2,edge))
		            a(3)=z(ii)-msh%xyz(3,msh%info_unk(jj+2,edge))
		            phase=junit*ker%wavenum*(x(ii)*sin(theta*pi/180.)*cos(phi*pi/180.)+y(ii)*sin(theta*pi/180.)*sin(phi*pi/180.)+z(ii)*cos(theta*pi/180.))
		            do i=1,3
                        ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*curr(edge)*exp(phase)*w(ii)
                        !ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*vectors_block(0)%vector(edge,1)*exp(phase)*w(ii)  !current
		            enddo
                enddo
            enddo
            ctemp_1=impedence0*(cos(theta*pi/180.)*cos(phi*pi/180.)*ctemp_rcs(1)+cos(theta*pi/180.)*sin(phi*pi/180.)*ctemp_rcs(2)-sin(theta*pi/180.)*ctemp_rcs(3))
    
	deallocate(x,y,z,w)	
	
    return
    
end subroutine VV_polar_SURF

subroutine HH_polar_SURF(theta,phi,edge,ctemp_1,curr,msh,ker)
    
    use MODULE_FILE
    implicit none
    
    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1
    real*8 theta,phi
    type(kernelquant)::ker
    integer i,j,ii,jj,iii,jjj,patch,flag
    real*8 l_edge,l_edgefine
    real*8 a(3)
    complex(kind=8)::curr(:)
	type(mesh)::msh
	real*8,allocatable:: x(:),y(:),z(:),w(:)
	
    integer edge,edge_m,edge_n
    
	allocate(x(msh%integral_points))
	allocate(y(msh%integral_points))
	allocate(z(msh%integral_points))
	allocate(w(msh%integral_points))		
	
        ctemp_rcs(1:3)=(0,0)
        l_edge=sqrt((msh%xyz(1,msh%info_unk(1,edge))-msh%xyz(1,msh%info_unk(2,edge)))**2+(msh%xyz(2,msh%info_unk(1,edge))-msh%xyz(2,msh%info_unk(2,edge)))**2+(msh%xyz(3,msh%info_unk(1,edge))-msh%xyz(3,msh%info_unk(2,edge)))**2)
       
			do jj=3,4  ! two triganle
	            !patch=msh%info_unk(jj+2,edgefine_m)
                call gau_grobal(edge,jj,x,y,z,w,msh)
		        do ii=1,msh%integral_points
		        
				    a(1)=x(ii)-msh%xyz(1,msh%info_unk(jj+2,edge))  !free node coordinate
		            a(2)=y(ii)-msh%xyz(2,msh%info_unk(jj+2,edge))
		            a(3)=z(ii)-msh%xyz(3,msh%info_unk(jj+2,edge))
		            phase=junit*ker%wavenum*(x(ii)*sin(theta*pi/180.)*cos(phi*pi/180.)+y(ii)*sin(theta*pi/180.)*sin(phi*pi/180.)+z(ii)*cos(theta*pi/180.))
		            do i=1,3
                        ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*curr(edge)*exp(phase)*w(ii)
                        !ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*vectors_block(0)%vector(edge,1)*exp(phase)*w(ii)  !current
		            enddo
                enddo
	        enddo
            ctemp_1=impedence0*(-sin(phi*pi/180.)*ctemp_rcs(1)+cos(phi*pi/180.)*ctemp_rcs(2))
	deallocate(x,y,z,w)	
    
    return
    
end subroutine HH_polar_SURF




subroutine RCS_bistatic_CURV(curr,msh,ker)
    !integer flag
	use MODULE_FILE
    implicit none
    complex(kind=8)::curr(:)
    real*8 rcs
    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1
    real*8 ddphi,dphi
    
    integer i,j,ii,jj,iii,jjj,patch,flag
    real*8 l_edge,l_edgefine
	type(mesh)::msh
    integer edge,edge_m,edge_n
    type(kernelquant)::ker 
    open (100, file='VV_bistatic.txt')
  
    ddphi=180./ker%RCS_Nsample
    
    do i=0, ker%RCS_Nsample   !phi=0
        dphi=i*ddphi
        ctemp=0
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp)   
        do edge=1,msh%Nunk
            call VV_polar_CURV(dphi,edge,ctemp_1,curr,msh,ker)
            ctemp=ctemp+ctemp_1    
        enddo
        !$omp end parallel do
        rcs=(abs(impedence0*ctemp))**2/4d0*ker%wavenum
        !rcs=rcs/ker%wavelength
        rcs=10*log10(rcs)
        write(100,*)dphi,rcs
    enddo
    close(100)
    
    return
    
end subroutine RCS_bistatic_CURV

subroutine VV_polar_CURV(dphi,edge,ctemp_1,curr,msh,ker)
    
    use MODULE_FILE
    implicit none
    complex(kind=8)::curr(:)
    complex(kind=8) ctemp,phase,ctemp_1
    real*8 dsita,dphi    
    integer edge
    type(mesh)::msh
	type(kernelquant)::ker
	
    phase=junit*ker%wavenum*(msh%xyz(1,msh%info_unk(0,edge))*cos(dphi*pi/180.)+msh%xyz(2,msh%info_unk(0,edge))*sin(dphi*pi/180.))
    ctemp_1=curr(edge)*msh%Delta_ll*exp(phase)

    return
    
end subroutine VV_polar_CURV



end module RCS_Bi