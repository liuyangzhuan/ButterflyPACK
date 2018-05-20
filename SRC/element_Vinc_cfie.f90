module element_Vinc
use gauss_rule
use analytic_part
contains 


subroutine element_Vinc_VV_CURV(phi,edge,value,msh,ker)

    use MODULE_FILE
    implicit none
    
    integer edge
    complex(kind=8) value
    real*8 theta, phi
    complex(kind=8)  phase
	type(mesh)::msh
	type(kernelquant)::ker
	
    phase=junit*ker%wavenum*(msh%xyz(1,msh%info_unk(0,edge))*cos(phi*pi/180.)+msh%xyz(2,msh%info_unk(0,edge))*sin(phi*pi/180.))
    value=exp(phase)
    
    return
    
end subroutine element_Vinc_VV_CURV



subroutine element_Vinc_VV_SURF(theta,phi,edge,value,msh,ker)

    use MODULE_FILE
    implicit none
    
    integer edge
    complex(kind=8) value
    real*8 theta, phi
    integer i,ii,jj,node1,node2,node3
    type(mesh)::msh
	real*8 lm,einc(3),a(3),nr(3),k(3)
	real*8, allocatable::x(:),y(:),z(:),w(:)
    complex(kind=8)  phase,ctemp_e,ctemp_h,value_e,value_h,ee(3),hh(3),hh1(3)
	type(kernelquant)::ker
	
	allocate(x(msh%integral_points))
	allocate(y(msh%integral_points))
	allocate(z(msh%integral_points))
	allocate(w(msh%integral_points))
	
	
    einc(1)=cos(theta*pi/180.)*cos(phi*pi/180.)
    einc(2)=cos(theta*pi/180.)*sin(phi*pi/180.)
    einc(3)=-sin(theta*pi/180.)
    k(1)=sin(theta*pi/180.)*cos(phi*pi/180.)
    k(2)=sin(theta*pi/180.)*sin(phi*pi/180.)
    k(3)=cos(theta*pi/180.)
    
    value_h=(0,0); value_e=(0,0)
        	
    node1=msh%info_unk(1,edge)
    node2=msh%info_unk(2,edge)
    lm=sqrt((msh%xyz(1,node1)-msh%xyz(1,node2))**2+(msh%xyz(2,node1)-msh%xyz(2,node2))**2+(msh%xyz(3,node1)-msh%xyz(3,node2))**2)
    do jj=3,4
	   nr(1:3)=msh%normal_of_patch(1:3,msh%info_unk(jj,edge)) 
	   node3=msh%info_unk(jj+2,edge)
       call gau_grobal(edge,jj,x,y,z,w,msh)
	   do ii=1,msh%integral_points
          phase=junit*ker%wavenum*(sin(theta*pi/180.)*cos(phi*pi/180.)*x(ii)+sin(theta*pi/180.)*sin(phi*pi/180.)*y(ii)+cos(theta*pi/180.)*z(ii))
	      do i=1,3
		     ee(i)=einc(i)*exp(phase)
          end do
          call ccurl(-k,ee,hh1)
          call ccurl(nr,hh1,hh)
          a(1)=x(ii)-msh%xyz(1,node3)
	      a(2)=y(ii)-msh%xyz(2,node3)
	      a(3)=z(ii)-msh%xyz(3,node3)
	      call cscalar(hh,a,ctemp_h) 
	      call cscalar(ee,a,ctemp_e)
	      value_h=value_h+(-1)**(jj+1)*ctemp_h*w(ii)
	      value_e=value_e+(-1)**(jj+1)*ctemp_e*w(ii) 
	   enddo
    end do
    value=lm*(ker%CFIE_alpha*value_e+(1.-ker%CFIE_alpha)*value_h)
    ! value=lm*value_e
    
	deallocate(x,y,z,w)
	
    return
    
end subroutine element_Vinc_VV_SURF

subroutine element_Vinc_HH_SURF(theta,phi,edge,value,msh,ker)

    use MODULE_FILE
    implicit none
	
    type(kernelquant)::ker
    integer edge
    complex(kind=8) value
    real*8 theta, phi
    integer i,ii,jj,node1,node2,node3
	real*8 lm,einc(3),a(3),nr(3),k(3)
	real*8, allocatable::x(:),y(:),z(:),w(:)
    complex(kind=8)  phase,ctemp_e,ctemp_h,value_e,value_h,ee(3),hh(3),hh1(3)
	type(mesh)::msh
	
	allocate(x(msh%integral_points))
	allocate(y(msh%integral_points))
	allocate(z(msh%integral_points))
	allocate(w(msh%integral_points))	
	
    einc(1)=-sin(phi*pi/180.)
    einc(2)=cos(phi*pi/180.)
    einc(3)=0.
    k(1)=sin(theta*pi/180.)*cos(phi*pi/180.)
    k(2)=sin(theta*pi/180.)*sin(phi*pi/180.)
    k(3)=cos(theta*pi/180.)
    
    value_h=(0,0); value_e=(0,0)
        	
    node1=msh%info_unk(1,edge)
    node2=msh%info_unk(2,edge)
    lm=sqrt((msh%xyz(1,node1)-msh%xyz(1,node2))**2+(msh%xyz(2,node1)-msh%xyz(2,node2))**2+(msh%xyz(3,node1)-msh%xyz(3,node2))**2)
    do jj=3,4
	   nr(1:3)=msh%normal_of_patch(1:3,msh%info_unk(jj,edge)) 
	   node3=msh%info_unk(jj+2,edge)
       call gau_grobal(edge,jj,x,y,z,w,msh)
	   do ii=1,msh%integral_points
          phase=junit*ker%wavenum*(sin(theta*pi/180.)*cos(phi*pi/180.)*x(ii)+sin(theta*pi/180.)*sin(phi*pi/180.)*y(ii)+cos(theta*pi/180.)*z(ii))
	      do i=1,3
		     ee(i)=einc(i)*exp(phase)
          end do
          call ccurl(-k,ee,hh1)
          call ccurl(nr,hh1,hh)
          a(1)=x(ii)-msh%xyz(1,node3)
	      a(2)=y(ii)-msh%xyz(2,node3)
	      a(3)=z(ii)-msh%xyz(3,node3)
	      call cscalar(hh,a,ctemp_h) 
	      call cscalar(ee,a,ctemp_e)
	      value_h=value_h+(-1)**(jj+1)*ctemp_h*w(ii)
	      value_e=value_e+(-1)**(jj+1)*ctemp_e*w(ii) 
	   enddo
    end do
    value=lm*(ker%CFIE_alpha*value_e+(1.-ker%CFIE_alpha)*value_h)
    ! value=lm*value_e
	deallocate(x,y,z,w)    
    return
    
end subroutine element_Vinc_HH_SURF
end module element_Vinc