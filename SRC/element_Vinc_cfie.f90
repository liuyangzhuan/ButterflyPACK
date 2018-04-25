module element_Vinc
use gauss_rule
use analytic_part
contains 
subroutine element_Vinc_VV(theta,phi,edge,value)

    use MODULE_FILE
    implicit none
    
    integer edge
    complex(kind=8) value
    real*8 theta, phi
    integer i,ii,jj,node1,node2,node3
    real*8 lm,einc(3),a(3),x(integral_points),y(integral_points),z(integral_points),w(integral_points),nr(3),k(3)
    complex(kind=8)  phase,ctemp_e,ctemp_h,value_e,value_h,ee(3),hh(3),hh1(3)

    einc(1)=cos(theta*pi/180.)*cos(phi*pi/180.)
    einc(2)=cos(theta*pi/180.)*sin(phi*pi/180.)
    einc(3)=-sin(theta*pi/180.)
    k(1)=sin(theta*pi/180.)*cos(phi*pi/180.)
    k(2)=sin(theta*pi/180.)*sin(phi*pi/180.)
    k(3)=cos(theta*pi/180.)
    
    value_h=(0,0); value_e=(0,0)
        	
    node1=node_patch_of_edge(1,edge)
    node2=node_patch_of_edge(2,edge)
    lm=sqrt((xyz(1,node1)-xyz(1,node2))**2+(xyz(2,node1)-xyz(2,node2))**2+(xyz(3,node1)-xyz(3,node2))**2)
    do jj=3,4
	   nr(1:3)=normal_of_patch(1:3,node_patch_of_edge(jj,edge)) 
	   node3=node_patch_of_edge(jj+2,edge)
       call gau_grobal(edge,jj,x,y,z,w)
	   do ii=1,integral_points
          phase=junit*wavenum*(sin(theta*pi/180.)*cos(phi*pi/180.)*x(ii)+sin(theta*pi/180.)*sin(phi*pi/180.)*y(ii)+cos(theta*pi/180.)*z(ii))
	      do i=1,3
		     ee(i)=einc(i)*exp(phase)
          end do
          call ccurl(-k,ee,hh1)
          call ccurl(nr,hh1,hh)
          a(1)=x(ii)-xyz(1,node3)
	      a(2)=y(ii)-xyz(2,node3)
	      a(3)=z(ii)-xyz(3,node3)
	      call cscalar(hh,a,ctemp_h) 
	      call cscalar(ee,a,ctemp_e)
	      value_h=value_h+(-1)**(jj+1)*ctemp_h*w(ii)
	      value_e=value_e+(-1)**(jj+1)*ctemp_e*w(ii) 
	   enddo
    end do
    value=lm*(CFIE_alpha*value_e+(1.-CFIE_alpha)*value_h)
    ! value=lm*value_e
    
    return
    
end subroutine element_Vinc_VV

subroutine element_Vinc_HH(theta,phi,edge,value)

    use MODULE_FILE
    implicit none
    
    integer edge
    complex(kind=8) value
    real*8 theta, phi
    integer i,ii,jj,node1,node2,node3
    real*8 lm,einc(3),a(3),x(integral_points),y(integral_points),z(integral_points),w(integral_points),nr(3),k(3)
    complex(kind=8)  phase,ctemp_e,ctemp_h,value_e,value_h,ee(3),hh(3),hh1(3)

    einc(1)=-sin(phi*pi/180.)
    einc(2)=cos(phi*pi/180.)
    einc(3)=0.
    k(1)=sin(theta*pi/180.)*cos(phi*pi/180.)
    k(2)=sin(theta*pi/180.)*sin(phi*pi/180.)
    k(3)=cos(theta*pi/180.)
    
    value_h=(0,0); value_e=(0,0)
        	
    node1=node_patch_of_edge(1,edge)
    node2=node_patch_of_edge(2,edge)
    lm=sqrt((xyz(1,node1)-xyz(1,node2))**2+(xyz(2,node1)-xyz(2,node2))**2+(xyz(3,node1)-xyz(3,node2))**2)
    do jj=3,4
	   nr(1:3)=normal_of_patch(1:3,node_patch_of_edge(jj,edge)) 
	   node3=node_patch_of_edge(jj+2,edge)
       call gau_grobal(edge,jj,x,y,z,w)
	   do ii=1,integral_points
          phase=junit*wavenum*(sin(theta*pi/180.)*cos(phi*pi/180.)*x(ii)+sin(theta*pi/180.)*sin(phi*pi/180.)*y(ii)+cos(theta*pi/180.)*z(ii))
	      do i=1,3
		     ee(i)=einc(i)*exp(phase)
          end do
          call ccurl(-k,ee,hh1)
          call ccurl(nr,hh1,hh)
          a(1)=x(ii)-xyz(1,node3)
	      a(2)=y(ii)-xyz(2,node3)
	      a(3)=z(ii)-xyz(3,node3)
	      call cscalar(hh,a,ctemp_h) 
	      call cscalar(ee,a,ctemp_e)
	      value_h=value_h+(-1)**(jj+1)*ctemp_h*w(ii)
	      value_e=value_e+(-1)**(jj+1)*ctemp_e*w(ii) 
	   enddo
    end do
    value=lm*(CFIE_alpha*value_e+(1.-CFIE_alpha)*value_h)
    ! value=lm*value_e
    
    return
    
end subroutine element_Vinc_HH
end module element_Vinc