module RCS_Bi
use gauss_rule
use analytic_part
use current_mapping
contains 

subroutine RCS_bistatic_SURF()
    !integer flag
	use MODULE_FILE
    implicit none
    
    real*8 rcs
    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1
    real*8 theta,phi,dtheta,dphi
    
    integer i,j,ii,jj,iii,jjj,patch,flag
    real*8 l_edge,l_edgefine
    real*8 a(3),x(integral_points),y(integral_points),z(integral_points),w(integral_points)
    
    integer edge,edge_m,edge_n
         
!     if (preconditioner==0) then
!         open(100,file='VV_bistatic_N.txt')
!     elseif (preconditioner==1) then
!         if (Diag_level==1) then
!             open(100,file='VV_bistatic_D1.txt')
!         elseif (Diag_level==2) then
!             open(100,file='VV_bistatic_D2.txt')
!         endif
!     elseif (preconditioner==2) then
!         if (SAI_level==1) then
!             open(100,file='VV_bistatic_S1.txt')
!         elseif (SAI_level==2) then
!             open(100,file='VV_bistatic_S2.txt')
!         elseif (SAI_level==3) then
!             open(100,file='VV_bistatic_S3.txt')
!         endif
!     endif
    open (100, file='VV_bistatic.txt')
	Current = Current2com(:,1)		
    theta=90.
    dphi=180./RCS_sample

    do i=0, RCS_sample
        phi=i*dphi 
        ctemp=(0.,0.)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp)   
        do edge=1,maxedge
            call VV_polar_SURF(theta,phi,edge,ctemp_1)
            ctemp=ctemp+ctemp_1    
        enddo
        !$omp end parallel do
        rcs=(abs(wavenum*ctemp))**2/4/pi
        !rcs=rcs/wavelength
        rcs=10*log10(rcs)
        write(100,*)phi,rcs
    enddo
    close(100)

	call current_node_patch_mapping('V')    	
	
	
   
!     if (preconditioner==0) then
!         open(1000,file='HH_bistatic_N.txt')
!     elseif (preconditioner==1) then
!          if (Diag_level==1) then
!             open(1000,file='HH_bistatic_D1.txt')
!         elseif (Diag_level==2) then
!             open(1000,file='HH_bistatic_D2.txt')
!         endif
!     elseif (preconditioner==2) then
!         if (SAI_level==1) then
!             open(1000,file='HH_bistatic_S1.txt')
!         elseif (SAI_level==2) then
!             open(1000,file='HH_bistatic_S2.txt')
!         elseif (SAI_level==3) then
!             open(1000,file='HH_bistatic_S3.txt')
!         endif
!     endif
    open (1000, file='HH_bistatic.txt')
	Current = Current2com(:,2)	    
    do i=0, RCS_sample
        phi=i*dphi
        ctemp=(0,0)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp)   
        do edge=1,maxedge
            call HH_polar_SURF(theta,phi,edge,ctemp_1)
            ctemp=ctemp+ctemp_1    
        enddo
        !$omp end parallel do
        rcs=(abs(wavenum*ctemp))**2/4/pi
        !rcs=rcs/wavelength
        rcs=10*log10(rcs)
        write(1000,*)phi,rcs
    enddo
    close(1000)
    call current_node_patch_mapping('H')    
	
    return
    
end subroutine RCS_bistatic_SURF

subroutine VV_polar_SURF(theta,phi,edge,ctemp_1)
    
    use MODULE_FILE
    implicit none
    
    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1
    real*8 theta,phi
    
    integer i,j,ii,jj,iii,jjj,patch,flag
    real*8 l_edge,l_edgefine
    real*8 a(3),x(integral_points),y(integral_points),z(integral_points),w(integral_points)
    
    integer edge,edge_m,edge_n
    
            ctemp_rcs(1:3)=(0,0)          
            l_edge=sqrt((xyz(1,node_patch_of_edge(1,edge))-xyz(1,node_patch_of_edge(2,edge)))**2+(xyz(2,node_patch_of_edge(1,edge))-xyz(2,node_patch_of_edge(2,edge)))**2+(xyz(3,node_patch_of_edge(1,edge))-xyz(3,node_patch_of_edge(2,edge)))**2)
			do jj=3,4  ! two triganle
	            !patch=node_patch_of_edge(jj,edgefine_m)
                call gau_grobal(edge,jj,x,y,z,w)

		        do ii=1,integral_points
		           
				    a(1)=x(ii)-xyz(1,node_patch_of_edge(jj+2,edge))  
		            a(2)=y(ii)-xyz(2,node_patch_of_edge(jj+2,edge))
		            a(3)=z(ii)-xyz(3,node_patch_of_edge(jj+2,edge))
		            phase=junit*wavenum*(x(ii)*sin(theta*pi/180.)*cos(phi*pi/180.)+y(ii)*sin(theta*pi/180.)*sin(phi*pi/180.)+z(ii)*cos(theta*pi/180.))
		            do i=1,3
                        ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*current(edge)*exp(phase)*w(ii)
                        !ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*vectors_block(0)%vector(edge,1)*exp(phase)*w(ii)  !current
		            enddo
                enddo
            enddo
            ctemp_1=impedence*(cos(theta*pi/180.)*cos(phi*pi/180.)*ctemp_rcs(1)+cos(theta*pi/180.)*sin(phi*pi/180.)*ctemp_rcs(2)-sin(theta*pi/180.)*ctemp_rcs(3))
    
    return
    
end subroutine VV_polar_SURF

subroutine HH_polar_SURF(theta,phi,edge,ctemp_1)
    
    use MODULE_FILE
    implicit none
    
    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1
    real*8 theta,phi
    
    integer i,j,ii,jj,iii,jjj,patch,flag
    real*8 l_edge,l_edgefine
    real*8 a(3),x(integral_points),y(integral_points),z(integral_points),w(integral_points)
    
    integer edge,edge_m,edge_n
    
        ctemp_rcs(1:3)=(0,0)
        l_edge=sqrt((xyz(1,node_patch_of_edge(1,edge))-xyz(1,node_patch_of_edge(2,edge)))**2+(xyz(2,node_patch_of_edge(1,edge))-xyz(2,node_patch_of_edge(2,edge)))**2+(xyz(3,node_patch_of_edge(1,edge))-xyz(3,node_patch_of_edge(2,edge)))**2)
       
			do jj=3,4  ! two triganle
	            !patch=node_patch_of_edge(jj+2,edgefine_m)
                call gau_grobal(edge,jj,x,y,z,w)
		        do ii=1,integral_points
		        
				    a(1)=x(ii)-xyz(1,node_patch_of_edge(jj+2,edge))  !free node coordinate
		            a(2)=y(ii)-xyz(2,node_patch_of_edge(jj+2,edge))
		            a(3)=z(ii)-xyz(3,node_patch_of_edge(jj+2,edge))
		            phase=junit*wavenum*(x(ii)*sin(theta*pi/180.)*cos(phi*pi/180.)+y(ii)*sin(theta*pi/180.)*sin(phi*pi/180.)+z(ii)*cos(theta*pi/180.))
		            do i=1,3
                        ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*current(edge)*exp(phase)*w(ii)
                        !ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*vectors_block(0)%vector(edge,1)*exp(phase)*w(ii)  !current
		            enddo
                enddo
	        enddo
            ctemp_1=impedence*(-sin(phi*pi/180.)*ctemp_rcs(1)+cos(phi*pi/180.)*ctemp_rcs(2))
    
    return
    
end subroutine HH_polar_SURF




subroutine RCS_bistatic_CURV()
    !integer flag
	use MODULE_FILE
    implicit none
    
    real*8 rcs
    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1
    real*8 ddphi,dphi
    
    integer i,j,ii,jj,iii,jjj,patch,flag
    real*8 l_edge,l_edgefine

    integer edge,edge_m,edge_n
         
    open (100, file='VV_bistatic.txt')
    Current = Current2com(:,1)	  
    ddphi=180./RCS_sample
    
    do i=0, RCS_sample   !phi=0
        dphi=i*ddphi
        ctemp=0
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp)   
        do edge=1,maxedge
            call VV_polar_CURV(dphi,edge,ctemp_1)
            ctemp=ctemp+ctemp_1    
        enddo
        !$omp end parallel do
        rcs=(abs(impedence*ctemp))**2/4d0*wavenum
        !rcs=rcs/wavelength
        rcs=10*log10(rcs)
        write(100,*)dphi,rcs
    enddo
    close(100)
    
    return
    
end subroutine RCS_bistatic_CURV

subroutine VV_polar_CURV(dphi,edge,ctemp_1)
    
    use MODULE_FILE
    implicit none
    
    complex(kind=8) ctemp,phase,ctemp_1
    real*8 dsita,dphi    
    integer edge
    
    phase=junit*wavenum*(xyz(1,node_patch_of_edge(0,edge))*cos(dphi*pi/180.)+xyz(2,node_patch_of_edge(0,edge))*sin(dphi*pi/180.))
    ctemp_1=current(edge)*Delta_ll*exp(phase)

    return
    
end subroutine VV_polar_CURV



end module RCS_Bi