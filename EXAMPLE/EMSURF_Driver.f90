module APPLICATION_MODULE
use MODULE_FILE
use MISC
implicit none

	!**** define your application-related variables here   
	type quant_app
		real*8 wavenum    ! CEM: wave number  
		real*8 wavelength  ! CEM: wave length
		real*8 omiga       ! CEM: angular frequency
		real*8 rank_approximate_para1, rank_approximate_para2, rank_approximate_para3 ! CEM: rank estimation parameter  
		integer RCS_static  ! CEM: 1: monostatic or 2: bistatic RCS
		integer RCS_Nsample ! CEM: number of RCS samples
		real*8:: CFIE_alpha ! CEM: combination parameter in CFIE
	end type quant_app

contains


!**** user-defined subroutine to sample Z_mn
subroutine Z_elem_EMSURF(ker,m,n,value_e,msh,quant)
	
    use MODULE_FILE
    implicit none
	
    integer, INTENT(IN):: m,n
    integer flag,edge_m,edge_n,nodetemp_n,patch
    complex(kind=8) value_e,value_m,value
    integer i,j,ii,jj,iii,jjj
    real*8 ln,lm,am(3),an(3),nr_m(3)
    real*8 nr_n(3)
    complex(kind=8) ctemp,ctemp1,ctemp2,aa(3),bb(1),dg(3),dg1(3),dg2(3)
    complex(kind=8) imp,imp1,imp2,imp3
    real*8 temp
    real*8 distance
    real*8 ianl,ianl1,ianl2
    real*8 area
	type(mesh)::msh
	class(kernelquant)::ker
	class(*),pointer :: quant
	
	
    real*8,allocatable::xm(:),ym(:),zm(:),wm(:),xn(:),yn(:),zn(:),wn(:)
	
	
	select TYPE(quant)
	type is (quant_app)	

		! convert to new indices because msh%info_unk has been reordered
		edge_m = msh%old2new(m)
		edge_n = msh%old2new(n)	
	
		allocate (xm(msh%integral_points), ym(msh%integral_points), zm(msh%integral_points), wm(msh%integral_points))
		allocate (xn(msh%integral_points), yn(msh%integral_points), zn(msh%integral_points), wn(msh%integral_points))
		
		
		lm=sqrt((msh%xyz(1,msh%info_unk(1,edge_m))-msh%xyz(1,msh%info_unk(2,edge_m)))**2+(msh%xyz(2,msh%info_unk(1,edge_m))-msh%xyz(2,msh%info_unk(2,edge_m)))**2+(msh%xyz(3,msh%info_unk(1,edge_m))-msh%xyz(3,msh%info_unk(2,edge_m)))**2)
		ln=sqrt((msh%xyz(1,msh%info_unk(1,edge_n))-msh%xyz(1,msh%info_unk(2,edge_n)))**2+(msh%xyz(2,msh%info_unk(1,edge_n))-msh%xyz(2,msh%info_unk(2,edge_n)))**2+(msh%xyz(3,msh%info_unk(1,edge_n))-msh%xyz(3,msh%info_unk(2,edge_n)))**2)
		
		ctemp1=(0.,0.)
		ctemp2=(0.,0.)
		value_m=(0.,0.)
		
		do ii=3,4
			call gau_grobal(edge_m,ii,xm,ym,zm,wm,msh)                        	       
			nr_m(1:3)=msh%normal_of_patch(1:3,msh%info_unk(ii,edge_m))
			do i=1,msh%integral_points
				am(1)=xm(i)-msh%xyz(1,msh%info_unk(ii+2,edge_m))
				am(2)=ym(i)-msh%xyz(2,msh%info_unk(ii+2,edge_m))
				am(3)=zm(i)-msh%xyz(3,msh%info_unk(ii+2,edge_m))
				bb(1)=(0.,0.)
				aa(1:3)=(0.,0.)
				do jj=3,4                
					call gau_grobal(edge_n,jj,xn,yn,zn,wn,msh)
					nr_n(1:3)=msh%normal_of_patch(1:3,msh%info_unk(jj,edge_n))
					if (msh%info_unk(ii,edge_m)==msh%info_unk(jj,edge_n)) then
						area=triangle_area(msh%info_unk(ii,edge_m),msh)
						imp=(0.,0.);imp1=(0.,0.);imp2=(0.,0.);imp3=(0.,0.)
						an(1)=xm(i)-msh%xyz(1,msh%info_unk(jj+2,edge_n))
						an(2)=ym(i)-msh%xyz(2,msh%info_unk(jj+2,edge_n))
						an(3)=zm(i)-msh%xyz(3,msh%info_unk(jj+2,edge_n))
						call scalar(am,an,temp)    
						value_m=value_m+(-1)**(ii+1)*(-1)**(jj+1)*0.5*temp/(2.*area)*wm(i)
						do j=1,msh%integral_points
							distance=sqrt((xm(i)-xn(j))**2+(ym(i)-yn(j))**2+(zm(i)-zn(j))**2)			                                     
							if(distance==0)then    
								imp=imp+wn(j)*(-junit*quant%wavenum)
								imp1=imp1+msh%ng1(j)*wn(j)*(-junit*quant%wavenum)
								imp2=imp2+msh%ng2(j)*wn(j)*(-junit*quant%wavenum)
								ianl=ianalytic(edge_n,jj,xn(j),yn(j),zn(j),msh)
								ianl1=ianalytic2(edge_n,jj,xn(j),yn(j),zn(j),1,msh)
								ianl2=ianalytic2(edge_n,jj,xn(j),yn(j),zn(j),2,msh)
								imp=imp+ianl  !debug
								imp1=imp1+ianl1
								imp2=imp2+ianl2
							else
								imp=imp+wn(j)*exp(-junit*quant%wavenum*distance)/distance
								imp1=imp1+msh%ng1(j)*wn(j)*exp(-junit*quant%wavenum*distance)/distance
								imp2=imp2+msh%ng2(j)*wn(j)*exp(-junit*quant%wavenum*distance)/distance
							endif                        
						enddo
						imp3=imp-imp1-imp2
						nodetemp_n=msh%info_unk(jj+2,edge_n)
						patch=msh%info_unk(jj,edge_n)
						do jjj=1,3
							aa(jjj)=aa(jjj)+(-1)**(jj+1)*quant%wavenum**2*(msh%xyz(jjj,msh%node_of_patch(1,patch))*imp1+msh%xyz(jjj,msh%node_of_patch(2,patch))*imp2+msh%xyz(jjj,msh%node_of_patch(3,patch))*imp3-msh%xyz(jjj,nodetemp_n)*imp)          
						enddo
						bb(1)=bb(1)+(-1)**(jj+1)*imp
					else
						imp=(0.,0.);imp1=(0.,0.);imp2=(0.,0.);imp3=(0.,0.)
						do j=1,msh%integral_points
							distance=sqrt((xm(i)-xn(j))**2+(ym(i)-yn(j))**2+(zm(i)-zn(j))**2)
							an(1)=xn(j)-msh%xyz(1,msh%info_unk(jj+2,edge_n))
							an(2)=yn(j)-msh%xyz(2,msh%info_unk(jj+2,edge_n))
							an(3)=zn(j)-msh%xyz(3,msh%info_unk(jj+2,edge_n))
							dg(1)=(xm(i)-xn(j))*(1+junit*quant%wavenum*distance)*exp(-junit*quant%wavenum*distance)/(4*pi*distance**3)
							dg(2)=(ym(i)-yn(j))*(1+junit*quant%wavenum*distance)*exp(-junit*quant%wavenum*distance)/(4*pi*distance**3)
							dg(3)=(zm(i)-zn(j))*(1+junit*quant%wavenum*distance)*exp(-junit*quant%wavenum*distance)/(4*pi*distance**3)
							 call ccurl(an,dg,dg1)
							 call ccurl(nr_m,dg1,dg2)
							 call cscalar(dg2,am,ctemp)
							 value_m=value_m-(-1)**(ii+1)*(-1)**(jj+1)*ctemp*wm(i)*wn(j)			             
							imp=imp+wn(j)*exp(-junit*quant%wavenum*distance)/distance
							imp1=imp1+msh%ng1(j)*wn(j)*exp(-junit*quant%wavenum*distance)/distance
							imp2=imp2+msh%ng2(j)*wn(j)*exp(-junit*quant%wavenum*distance)/distance
						enddo                        
						imp3=imp-imp1-imp2
						nodetemp_n=msh%info_unk(jj+2,edge_n)
						patch=msh%info_unk(jj,edge_n)
						do jjj=1,3
						   aa(jjj)=aa(jjj)+(-1)**(jj+1)*quant%wavenum**2*(msh%xyz(jjj,msh%node_of_patch(1,patch))*imp1+msh%xyz(jjj,msh%node_of_patch(2,patch))*imp2+msh%xyz(jjj,msh%node_of_patch(3,patch))*imp3-msh%xyz(jjj,nodetemp_n)*imp)          
						enddo
						bb(1)=bb(1)+(-1)**(jj+1)*imp                                                                                                    
					endif
				enddo	
				call cscalar(aa,am,ctemp)
				ctemp1=ctemp1+(-1)**(ii+1)*ctemp*wm(i) 
				ctemp2=ctemp2+4.*(-1)**(ii+1)*bb(1)*wm(i)    
			enddo
		enddo
		value_e=ln*lm*junit*(ctemp1-ctemp2)/4./pi/quant%omiga/eps0
		value_m=value_m*lm*ln
		
		value=quant%CFIE_alpha*value_e+(1.-quant%CFIE_alpha)*impedence0*value_m
		
		deallocate(xm,ym,zm,wm,xn,yn,zn,wn)
	class default
		write(*,*)"unexpected type"
		stop
	end select	
			
    return
    
end subroutine Z_elem_EMSURF



!***********************************
!	seven points gauss integration
!	provide w(n) and x,y,z(n) (n=7)
!***********************************
  subroutine gau_grobal(nn,j,x,y,z,w,msh)
  
  use MODULE_FILE
  implicit none
  type(mesh)::msh
  integer nn ,flag
  real*8 x(:),y(:),z(:),w(:)
  integer i,j,ii
!	wa=area*9./40
!	wb=area*(155.-sqrt(15.))/1200.
!	wc=area*(155.+sqrt(15.))/1200.

   ! if(flag==1) then

!		ii=msh%info_unkfine(j,nn)
!  do i=1,integral_points
!	x(i)=msh%ng1(i)*msh%xyz(1,node_of_patchfine(1,ii))+msh%ng2(i)*msh%xyz(1,node_of_patchfine(2,ii))+&
!     msh%ng3(i)*msh%xyz(1,node_of_patchfine(3,ii))
!	y(i)=msh%ng1(i)*msh%xyz(2,node_of_patchfine(1,ii))+msh%ng2(i)*msh%xyz(2,node_of_patchfine(2,ii))+&
!     msh%ng3(i)*msh%xyz(2,node_of_patchfine(3,ii))
!	z(i)=msh%ng1(i)*msh%xyz(3,node_of_patchfine(1,ii))+msh%ng2(i)*msh%xyz(3,node_of_patchfine(2,ii))+&
!     msh%ng3(i)*msh%xyz(3,node_of_patchfine(3,ii))
!  enddo
  
!  elseif(flag==2) then

  ii=msh%info_unk(j,nn)
	 do i=1,msh%integral_points
	x(i)=msh%ng1(i)*msh%xyz(1,msh%node_of_patch(1,ii))+msh%ng2(i)*msh%xyz(1,msh%node_of_patch(2,ii))+&
     msh%ng3(i)*msh%xyz(1,msh%node_of_patch(3,ii))
	y(i)=msh%ng1(i)*msh%xyz(2,msh%node_of_patch(1,ii))+msh%ng2(i)*msh%xyz(2,msh%node_of_patch(2,ii))+&
     msh%ng3(i)*msh%xyz(2,msh%node_of_patch(3,ii))
	z(i)=msh%ng1(i)*msh%xyz(3,msh%node_of_patch(1,ii))+msh%ng2(i)*msh%xyz(3,msh%node_of_patch(2,ii))+&
     msh%ng3(i)*msh%xyz(3,msh%node_of_patch(3,ii))
  enddo
  
  w=msh%gauss_w

! 	w(1)=wa
! 	w(2)=wb
! 	w(3)=wb
! 	w(4)=wb
! 	w(5)=wc
! 	w(6)=wc
! 	w(7)=wc
  return
  end subroutine gau_grobal
  

  subroutine gauss_points(msh)
  
      use MODULE_FILE
      implicit none
      
      real*8 v1,v2,v3,v4,v5
      real*8 wa,wb,wc
	  type(mesh)::msh
	  
    !	wa=area*9./40
!	wb=area*(155.-sqrt(15.))/1200.
!	wc=area*(155.+sqrt(15.))/1200.
    if (msh%integral_points==7) then
    
    wa=9./80
	wb=(155.-sqrt(15.))/2400.
	wc=(155.+sqrt(15.))/2400.

	v1=1./3.
	v2=(6.-sqrt(15.))/21.
	v3=(9.+2*sqrt(15.))/21.
	v4=(6.+sqrt(15.))/21.
	v5=(9.-2*sqrt(15.))/21.

	msh%ng1(1)=1-v1-v1
	msh%ng1(2)=1-v2-v2
	msh%ng1(3)=1-v2-v3
	msh%ng1(4)=1-v3-v2
	msh%ng1(5)=1-v4-v4
	msh%ng1(6)=1-v4-v5
	msh%ng1(7)=1-v5-v4

	msh%ng2(1)=v1
	msh%ng2(2)=v2
	msh%ng2(3)=v2
	msh%ng2(4)=v3
	msh%ng2(5)=v4
	msh%ng2(6)=v4
	msh%ng2(7)=v5

	msh%ng3(1)=v1
	msh%ng3(2)=v2
	msh%ng3(3)=v3
	msh%ng3(4)=v2
	msh%ng3(5)=v4
	msh%ng3(6)=v5
	msh%ng3(7)=v4
	
	msh%gauss_w(1)=wa
	msh%gauss_w(2)=wb
	msh%gauss_w(3)=wb
	msh%gauss_w(4)=wb
	msh%gauss_w(5)=wc
	msh%gauss_w(6)=wc
	msh%gauss_w(7)=wc
	
	elseif (msh%integral_points==6) then
	
	v1=0.816847572980459d0
    v2=0.091576213509771d0
    v3=0.108103018168070d0
    v4=0.445948490915965d0
    wa=0.109951743655322d0/2.
    wb=0.223381589678011d0/2.

	msh%ng1(1)=v1
	msh%ng1(2)=v2
	msh%ng1(3)=v2
	msh%ng1(4)=v3
	msh%ng1(5)=v4
	msh%ng1(6)=v4
	!msh%ng1(7)=1-v5-v4

	msh%ng2(1)=v2
	msh%ng2(2)=v1
	msh%ng2(3)=v2
	msh%ng2(4)=v4
	msh%ng2(5)=v3
	msh%ng2(6)=v4
	!msh%ng2(7)=v5

	msh%ng3(1)=v2
	msh%ng3(2)=v2
	msh%ng3(3)=v1
	msh%ng3(4)=v4
	msh%ng3(5)=v4
	msh%ng3(6)=v3
	
	msh%gauss_w(1)=wa
	msh%gauss_w(2)=wa
	msh%gauss_w(3)=wa
	msh%gauss_w(4)=wb
	msh%gauss_w(5)=wb
	msh%gauss_w(6)=wb
	
	elseif (msh%integral_points==4) then
	
	v1=1./3.
	v2=0.2
	v3=0.6
	wa=-27./96.
	wb=25./96.
	
	msh%ng1(1)=v1
	msh%ng1(2)=v2
	msh%ng1(3)=v2
	msh%ng1(4)=v3
	
	msh%ng2(1)=v1
	msh%ng2(2)=v2
	msh%ng2(3)=v3
	msh%ng2(4)=v2
	
	msh%ng3(1)=v1
	msh%ng3(2)=v3
	msh%ng3(3)=v2
	msh%ng3(4)=v2
	
	msh%gauss_w(1)=wa
	msh%gauss_w(2)=wb
	msh%gauss_w(3)=wb
	msh%gauss_w(4)=wb
	
	endif
    
    return
  end subroutine gauss_points


!**********************************
!	mm:commom edge
!	jj:face number(3 or 4)
!**********************************	
function ianalytic(mm,jj,xi,yi,zi,msh)

use     MODULE_FILE
integer mm,jj,j,i
real*8 xi,yi,zi
real*8    temp,ianalytic
integer ii,node1,node2,node3
real*8    u3,v3,u0,v0,w0,l(3)
real*8    u(3),w(3),v(3),a(3),b(3)
real*8    s(2,3),t(-1:1,3)
real*8    f2(3),beta(3)
real*8    r(-1:1,3)
real*8    area
type(mesh)::msh

ii=msh%info_unk(jj,mm)  
node1=msh%node_of_patch(1,ii)
node2=msh%node_of_patch(2,ii)
node3=msh%node_of_patch(3,ii)
do i=1,3
   a(i)=msh%xyz(i,node2)-msh%xyz(i,node1)
   b(i)=msh%xyz(i,node3)-msh%xyz(i,node1)
enddo
 call curl(a,b,w)
area=0.5*sqrt(w(1)**2+w(2)**2+w(3)**2)
do i=1,3
   w(i)=w(i)/2./area
enddo
l(1)=sqrt((msh%xyz(1,node3)-msh%xyz(1,node2))**2+(msh%xyz(2,node3)&
     	-msh%xyz(2,node2))**2+(msh%xyz(3,node3)-msh%xyz(3,node2))**2)
l(2)=sqrt((msh%xyz(1,node3)-msh%xyz(1,node1))**2+(msh%xyz(2,node3)&
     	-msh%xyz(2,node1))**2+(msh%xyz(3,node3)-msh%xyz(3,node1))**2) 
l(3)=sqrt((msh%xyz(1,node1)-msh%xyz(1,node2))**2+(msh%xyz(2,node1)&
     	-msh%xyz(2,node2))**2+(msh%xyz(3,node1)-msh%xyz(3,node2))**2) 
do i=1,3
   u(i)=a(i)/l(3)
enddo
 call curl(w,u,v)
 call scalar(u,b,u3)
v3=2.*area/l(3)
	
b(1)=xi-msh%xyz(1,node1)
b(2)=yi-msh%xyz(2,node1)
b(3)=zi-msh%xyz(3,node1)
 call scalar(u,b,u0)
 call scalar(v,b,v0)
 call scalar(w,b,w0)

s(1,1)=-((l(3)-u0)*(l(3)-u3)+v0*v3)/l(1)
s(2,1)=((u3-u0)*(u3-l(3))+v3*(v3-v0))/l(1)  
s(1,2)=-(u3*(u3-u0)+v3*(v3-v0))/l(2)
s(2,2)=(u0*u3+v0*v3)/l(2)     
s(1,3)=-u0
s(2,3)=l(3)-u0
t(0,1)=(v0*(u3-l(3))+v3*(l(3)-u0))/l(1)
t(0,2)=(u0*v3-v0*u3)/l(2)
t(0,3)=v0
t(-1,1)=sqrt((l(3)-u0)**2+v0**2)
t(1,1)=sqrt((u3-u0)**2+(v3-v0)**2)
t(1,2)=sqrt(u0**2+v0**2)
t(1,3)=t(-1,1)
t(-1,2)=t(1,1)
t(-1,3)=t(1,2)
	
do j=1,3
   do i=-1,1
      r(i,j)=sqrt(t(i,j)**2+w0**2)
   enddo
enddo

if((r(1,i)+s(2,i))/(r(-1,i)+s(1,i)).le.1.e-6)then
   print*,(r(1,i)+s(2,i))/(r(-1,i)+s(1,i))
   print*,"log value error!"
   print*,"ianalytic:",mm
   stop
endif

do i=1,3
   f2(i)=log((r(1,i)+s(2,i))/(r(-1,i)+s(1,i)))
   beta(i)=atan(t(0,i)*s(2,i)/(r(0,i)**2+abs(w0)*r(1,i)))-&
     atan(t(0,i)*s(1,i)/(r(0,i)**2+abs(w0)*r(-1,i)))
enddo

temp=0.
do i=1,3
   temp=temp+t(0,i)*f2(i)-abs(w0)*beta(i)
enddo
ianalytic=temp/area/2
	
return

end function ianalytic



!**********************************
!	mm:commom edge
!	jj:face number(3 or 4)
!	have the ni part(iii)
!**********************************	
function ianalytic2(mm,jj,xi,yi,zi,iii,msh)

use MODULE_FILE
integer mm,jj,j,i
real*8 ianalytic2
integer ii,node1,node2,node3,iii
real*8 xi,yi,zi
real*8 u3,v3,u0,v0,w0,l(3)
real*8 u(3),w(3),v(3),a(3),b(3)
real*8 m1(3),m2(3),m3(3)
real*8 s1(3),s2(3),s3(3)
real*8 s(2,3),t(-1:1,3)
real*8 f2(3),beta(3),f3(3)
real*8 r(-1:1,3)
real*8 temp,temp1,temp2,temp3
real*8 iua,iva
real*8 n0(3)
real*8    area
type(mesh)::msh

ii=msh%info_unk(jj,mm)
node1=msh%node_of_patch(1,ii)
node2=msh%node_of_patch(2,ii)
node3=msh%node_of_patch(3,ii)
do i=1,3
   a(i)=msh%xyz(i,node2)-msh%xyz(i,node1)
   b(i)=msh%xyz(i,node3)-msh%xyz(i,node1)
enddo
call curl(a,b,w)
area=0.5*sqrt(w(1)**2+w(2)**2+w(3)**2)
do i=1,3
   w(i)=w(i)/2./area
enddo
l(1)=sqrt((msh%xyz(1,node3)-msh%xyz(1,node2))**2+(msh%xyz(2,node3)&
     	-msh%xyz(2,node2))**2+(msh%xyz(3,node3)-msh%xyz(3,node2))**2)
l(2)=sqrt((msh%xyz(1,node3)-msh%xyz(1,node1))**2+(msh%xyz(2,node3)&
     	-msh%xyz(2,node1))**2+(msh%xyz(3,node3)-msh%xyz(3,node1))**2) 
l(3)=sqrt((msh%xyz(1,node1)-msh%xyz(1,node2))**2+(msh%xyz(2,node1)&
     	-msh%xyz(2,node2))**2+(msh%xyz(3,node1)-msh%xyz(3,node2))**2) 
do i=1,3
   u(i)=a(i)/l(3)
enddo
 call curl(w,u,v)
 call scalar(u,b,u3)
v3=2.*area/l(3)
	
b(1)=xi-msh%xyz(1,node1)
b(2)=yi-msh%xyz(2,node1)
b(3)=zi-msh%xyz(3,node1)
 call scalar(u,b,u0)
 call scalar(v,b,v0)
 call scalar(w,b,w0)

	
s(1,1)=-((l(3)-u0)*(l(3)-u3)+v0*v3)/l(1)
s(2,1)=((u3-u0)*(u3-l(3))+v3*(v3-v0))/l(1)  
s(1,2)=-(u3*(u3-u0)+v3*(v3-v0))/l(2)
s(2,2)=(u0*u3+v0*v3)/l(2)
s(1,3)=-u0
s(2,3)=l(3)-u0
t(0,1)=(v0*(u3-l(3))+v3*(l(3)-u0))/l(1)
t(0,2)=(u0*v3-v0*u3)/l(2)
t(0,3)=v0
t(-1,1)=sqrt((l(3)-u0)**2+v0**2)
t(1,1)=sqrt((u3-u0)**2+(v3-v0)**2)
t(1,2)=sqrt(u0**2+v0**2)
t(1,3)=t(-1,1)
t(-1,2)=t(1,1)
t(-1,3)=t(1,2)
	
do j=1,3
   do i=-1,1
      r(i,j)=sqrt(t(i,j)**2+w0**2)
   enddo
enddo

if((r(1,i)+s(2,i))/(r(-1,i)+s(1,i)).le.1.e-6)then
	print*,(r(1,i)+s(2,i))/(r(-1,i)+s(1,i))
	print*,"log value error!"
	print*,"ianalytic2:",mm
	stop
endif
	
do i=1,3
   f2(i)=log((r(1,i)+s(2,i))/(r(-1,i)+s(1,i)))
   beta(i)=atan(t(0,i)*s(2,i)/(r(0,i)**2+abs(w0)*r(1,i)))-&
     atan(t(0,i)*s(1,i)/(r(0,i)**2+abs(w0)*r(-1,i)))
enddo
temp=0

do i=1,3
   temp=temp+t(0,i)*f2(i)-abs(w0)*beta(i)  !Ianl
enddo

do i=1,3
   f3(i)=s(2,i)*r(1,i)-s(1,i)*r(-1,i)+r(0,i)**2*f2(i)
enddo
do i=1,3
   s1(i)=(msh%xyz(i,node3)-msh%xyz(i,node2))/l(1)
   s2(i)=(msh%xyz(i,node1)-msh%xyz(i,node3))/l(2)
   s3(i)=(msh%xyz(i,node2)-msh%xyz(i,node1))/l(3)
enddo
call curl(s1,w,m1)
call curl(s2,w,m2)
call curl(s3,w,m3)
call scalar(u,m1,temp1)
temp1=temp1*f3(1)/2
call scalar(u,m2,temp2)
temp2=temp2*f3(2)/2
call scalar(u,m3,temp3)
temp3=temp3*f3(3)/2
iua=temp1+temp2+temp3
	
call scalar(v,m1,temp1)
temp1=temp1*f3(1)/2
call scalar(v,m2,temp2)
temp2=temp2*f3(2)/2
call scalar(v,m3,temp3)
temp3=temp3*f3(3)/2
iva=temp1+temp2+temp3

n0(1)=1.-u0/l(3)+v0*(u3/l(3)-1.)/v3
n0(2)=u0/l(3)-u3*v0/l(3)/v3
n0(3)=v0/v3
	
select case(iii)
case(1)
ianalytic2=(n0(1)*temp-iua/l(3)+iva*(u3/l(3)-1.)/v3)/2/area  !debug  l(2)--->l(3)
case(2)
ianalytic2=(n0(2)*temp+iua/l(3)-iva*u3/l(3)/v3)/2/area
case(3)
ianalytic2=(n0(3)*temp+iva/v3)/2/area
case default
print*,"iii value error!"
end select	

return

end function ianalytic2


subroutine current_node_patch_mapping(chara,curr,msh)
    
    use MODULE_FILE
    implicit none
    
    integer patch, edge, node_patch(3), node_edge, node
    integer i,j,k,ii,jj,kk,flag
    real*8 center(3), current_patch(3,0:3),  current_abs, r, a
    character chara
    character(20) string 
    complex(kind=8)::curr(:)
	type(mesh)::msh
	
    real*8,allocatable :: current_at_patch(:), current_at_node(:)
    integer,allocatable :: edge_of_patch(:,:)
    
    allocate (edge_of_patch(3,msh%maxpatch))
    edge_of_patch = -1
	allocate (current_at_node(msh%maxnode),current_at_patch(msh%maxpatch))
    
    !$omp parallel do default(shared) private(patch,i,edge)
    do patch=1,msh%maxpatch
        i=0
        ! do while (i<3)     !!! Modified by Yang Liu, commented out, this doesn't make sense
            do edge=1, msh%Nunk            
                if (msh%info_unk(3,edge)==patch .or. msh%info_unk(4,edge)==patch) then
                    i=i+1
                    edge_of_patch(i,patch)=edge
                endif
            enddo
        ! enddo                
    enddo
    !$omp end parallel do
    
    !$omp parallel do default(shared) private(patch,i,j,edge,current_abs,current_patch,a,r)
    do patch=1, msh%maxpatch
        do i=1,3
            center(i)=1./3.*(msh%xyz(i,msh%node_of_patch(1,patch))+msh%xyz(i,msh%node_of_patch(2,patch))+msh%xyz(i,msh%node_of_patch(3,patch)))
        enddo
        do edge=1,3
			if(edge_of_patch(edge,patch)==-1)then
				current_patch(:,edge)=0
			else 
				current_abs=dble(curr(edge_of_patch(edge,patch)))
				if (msh%info_unk(3,edge_of_patch(edge,patch))==patch) then
					r=0.                
					do j=1,3
						a=msh%xyz(j,msh%info_unk(5,edge_of_patch(edge,patch)))-center(j)
						r=r+a**2
						current_patch(j,edge)=current_abs*a
					enddo
					r=sqrt(r)
					do j=1,3
						current_patch(j,edge)=current_patch(j,edge)/r
					enddo
				elseif (msh%info_unk(4,edge_of_patch(edge,patch))==patch) then
					r=0.                
					do j=1,3
						a=msh%xyz(j,msh%info_unk(6,edge_of_patch(edge,patch)))-center(j)
						r=r+a**2
						current_patch(j,edge)=-current_abs*a
					enddo
					r=sqrt(r)
					do j=1,3
						current_patch(j,edge)=current_patch(j,edge)/r
					enddo
				endif
			endif	
        enddo
        do i=1,3
            current_patch(i,0)=0.
            do edge=1,3
                current_patch(i,0)=current_patch(i,0)+current_patch(i,edge)
            enddo
        enddo
        current_at_patch(patch)=sqrt(current_patch(1,0)**2+current_patch(2,0)**2+current_patch(3,0)**2)
    enddo
    !$omp end parallel do
    
    !$omp parallel do default(shared) private(patch,i,ii,node,a)
    do node=1, msh%maxnode
        ii=0; a=0.
        do patch=1, msh%maxpatch
            do i=1,3
                if (msh%node_of_patch(i,patch)==node) then
                    a=a+current_at_patch(patch)
                    ii=ii+1
                endif
            enddo
        enddo
        current_at_node(node)=a/dble(ii)
    enddo
    !$omp end parallel do
    
    string='current'//chara//'.out'
    open(30,file=string)
    do node=1, msh%maxnode
        write (30,*) node,current_at_node(node)
    enddo
    close(30)
    
    deallocate (edge_of_patch,current_at_node,current_at_patch)
    return
    
end subroutine current_node_patch_mapping


real*8 function triangle_area(patch,msh)
    
    use MODULE_FILE
    implicit none
    
    integer patch,i
    real*8 a(3),b(3),c(3)
	type(mesh)::msh
    
    do i=1,3
        a(i)=msh%xyz(i,msh%node_of_patch(2,patch))-msh%xyz(i,msh%node_of_patch(1,patch))
        b(i)=msh%xyz(i,msh%node_of_patch(3,patch))-msh%xyz(i,msh%node_of_patch(1,patch))
    enddo
    
    call curl(a,b,c)
    triangle_area=0.5*sqrt(c(1)**2+c(2)**2+c(3)**2)
    
    return
end function triangle_area


subroutine element_Vinc_VV_SURF(theta,phi,edge,value,msh,quant)

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
	type(quant_app)::quant
	
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
          phase=junit*quant%wavenum*(sin(theta*pi/180.)*cos(phi*pi/180.)*x(ii)+sin(theta*pi/180.)*sin(phi*pi/180.)*y(ii)+cos(theta*pi/180.)*z(ii))
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
    value=lm*(quant%CFIE_alpha*value_e+(1.-quant%CFIE_alpha)*value_h)
    ! value=lm*value_e
    
	deallocate(x,y,z,w)
	
    return
    
end subroutine element_Vinc_VV_SURF

subroutine element_Vinc_HH_SURF(theta,phi,edge,value,msh,quant)

    use MODULE_FILE
    implicit none
	
    type(quant_app)::quant
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
          phase=junit*quant%wavenum*(sin(theta*pi/180.)*cos(phi*pi/180.)*x(ii)+sin(theta*pi/180.)*sin(phi*pi/180.)*y(ii)+cos(theta*pi/180.)*z(ii))
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
    value=lm*(quant%CFIE_alpha*value_e+(1.-quant%CFIE_alpha)*value_h)
    ! value=lm*value_e
	deallocate(x,y,z,w)    
    return
    
end subroutine element_Vinc_HH_SURF



subroutine RCS_bistatic_SURF(curr,msh,quant,ptree)
    !integer flag
	use MODULE_FILE
    implicit none
    
    real*8 rcs
    complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1
    real*8 theta,phi,dtheta,dphi
    
    integer i,j,ii,jj,iii,jjj,patch,flag
    real*8 l_edge,l_edgefine
    real*8 a(3)
    real*8,allocatable:: x(:),y(:),z(:),w(:)
    complex(kind=8):: curr(:,:)
	type(mesh)::msh
	type(quant_app)::quant
	type(proctree)::ptree
	
    integer edge,edge_m,edge_n,ierr
    
	allocate(x(msh%integral_points))
	allocate(y(msh%integral_points))
	allocate(z(msh%integral_points))
	allocate(w(msh%integral_points))
	
    if(ptree%MyID==Main_ID)open (100, file='VV_bistatic.txt')	
    theta=90.
    dphi=180./quant%RCS_Nsample

    do i=0, quant%RCS_Nsample
        phi=i*dphi 
        ctemp_loc=(0.,0.)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)   
        do edge=msh%idxs,msh%idxe
            call VV_polar_SURF(theta,phi,edge,ctemp_1,curr(edge-msh%idxs+1,1),msh,quant)
            ctemp_loc=ctemp_loc+ctemp_1    
        enddo
        !$omp end parallel do
		
		call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
		
        rcs=(abs(quant%wavenum*ctemp))**2/4/pi
        !rcs=rcs/quant%wavelength
        rcs=10*log10(rcs)
        if(ptree%MyID==Main_ID)write(100,*)phi,rcs
    enddo
    if(ptree%MyID==Main_ID)close(100)

   
    if(ptree%MyID==Main_ID)open (1000, file='HH_bistatic.txt')
	   
    do i=0, quant%RCS_Nsample
        phi=i*dphi
        ctemp_loc=(0,0)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)   
        do edge=msh%idxs,msh%idxe
            call HH_polar_SURF(theta,phi,edge,ctemp_1,curr(edge-msh%idxs+1,2),msh,quant)
            ctemp_loc=ctemp_loc+ctemp_1    
        enddo
        !$omp end parallel do
		
		call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
        rcs=(abs(quant%wavenum*ctemp))**2/4/pi
        !rcs=rcs/quant%wavelength
        rcs=10*log10(rcs)
        if(ptree%MyID==Main_ID)write(1000,*)phi,rcs
    enddo
    if(ptree%MyID==Main_ID)close(1000)
		
	deallocate(x,y,z,w)
	
    return
    
end subroutine RCS_bistatic_SURF

subroutine VV_polar_SURF(theta,phi,edge,ctemp_1,curr,msh,quant)
    
    use MODULE_FILE
    implicit none
    
    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1
    real*8 theta,phi
    type(quant_app)::quant
    integer i,j,ii,jj,iii,jjj,patch,flag
    real*8 l_edge,l_edgefine
    real*8 a(3)
    complex(kind=8)::curr
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
		            phase=junit*quant%wavenum*(x(ii)*sin(theta*pi/180.)*cos(phi*pi/180.)+y(ii)*sin(theta*pi/180.)*sin(phi*pi/180.)+z(ii)*cos(theta*pi/180.))
		            do i=1,3
                        ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*curr*exp(phase)*w(ii)
                        !ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*vectors_block(0)%vector(edge,1)*exp(phase)*w(ii)  !current
		            enddo
                enddo
            enddo
            ctemp_1=impedence0*(cos(theta*pi/180.)*cos(phi*pi/180.)*ctemp_rcs(1)+cos(theta*pi/180.)*sin(phi*pi/180.)*ctemp_rcs(2)-sin(theta*pi/180.)*ctemp_rcs(3))
    
	deallocate(x,y,z,w)	
	
    return
    
end subroutine VV_polar_SURF

subroutine HH_polar_SURF(theta,phi,edge,ctemp_1,curr,msh,quant)
    
    use MODULE_FILE
    implicit none
    
    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1
    real*8 theta,phi
    type(quant_app)::quant
    integer i,j,ii,jj,iii,jjj,patch,flag
    real*8 l_edge,l_edgefine
    real*8 a(3)
    complex(kind=8)::curr
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
		            phase=junit*quant%wavenum*(x(ii)*sin(theta*pi/180.)*cos(phi*pi/180.)+y(ii)*sin(theta*pi/180.)*sin(phi*pi/180.)+z(ii)*cos(theta*pi/180.))
		            do i=1,3
                        ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*curr*exp(phase)*w(ii)
                        !ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*vectors_block(0)%vector(edge,1)*exp(phase)*w(ii)  !current
		            enddo
                enddo
	        enddo
            ctemp_1=impedence0*(-sin(phi*pi/180.)*ctemp_rcs(1)+cos(phi*pi/180.)*ctemp_rcs(2))
	deallocate(x,y,z,w)	
    
    return
    
end subroutine HH_polar_SURF




subroutine RCS_monostatic_VV_SURF(dsita,dphi,rcs,curr,msh,quant,ptree)

    use MODULE_FILE
    implicit none
    complex(kind=8)::curr(:)
    real*8 rcs
    complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1,ctemp_2
    real*8 dsita,dphi
    integer edge,edge_m,edge_n,ierr
	type(mesh)::msh
    type(quant_app)::quant
	type(proctree)::ptree
	
    ctemp_loc=(0.,0.)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
        do edge=msh%idxs,msh%idxe
            call VV_polar_SURF(dsita,dphi,edge,ctemp_1,curr(edge-msh%idxs+1),msh,quant)
	        ctemp_loc=ctemp_loc+ctemp_1
        enddo
        !$omp end parallel do
		
		call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
		
        rcs=(abs(quant%wavenum*ctemp))**2/4/pi
        !rcs=rcs/quant%wavelength
        rcs=10*log10(rcs)
        
    return
    
end subroutine RCS_monostatic_VV_SURF

subroutine RCS_monostatic_HH_SURF(dsita,dphi,rcs,curr,msh,quant,ptree)

    use MODULE_FILE
    implicit none
    
    real*8 rcs
    complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1,ctemp_2
    real*8 dsita,dphi    
    integer edge,edge_m,edge_n,ierr
    complex(kind=8)::curr(:)
	type(mesh)::msh
	type(quant_app)::quant
	type(proctree)::ptree
	
    ctemp_loc=(0.,0.)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
        do edge=msh%idxs,msh%idxe
            call HH_polar_SURF(dsita,dphi,edge,ctemp_1,curr(edge-msh%idxs+1),msh,quant)
	        ctemp_loc=ctemp_loc+ctemp_1
        enddo
        !$omp end parallel do
		
		call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
		
        rcs=(abs(quant%wavenum*ctemp))**2/4/pi
        !rcs=rcs/quant%wavelength
        rcs=10*log10(rcs)
        
    return
    
end subroutine RCS_monostatic_HH_SURF



end module APPLICATION_MODULE



PROGRAM HODLR_BUTTERFLY_SOLVER_3D
    use MODULE_FILE
	use APPLICATION_MODULE
	! use geometry_model
	use H_structure
	use cascading_factorization
	use matrices_fill
	use omp_lib
	use MISC
    implicit none

	! include "mkl_vml.fi"	 
	
    real*8 para
    real*8 tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seed_myid(12)
	integer times(8)	
	real*8 t1,t2,x,y,z,r,theta,phi
	complex(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:)
	
	character(len=:),allocatable  :: string
	character(len=1024)  :: strings	
	character(len=6)  :: info_env	
	integer :: length
	integer :: ierr
	integer*8 oldmode,newmode
	type(Hoption)::option	
	type(Hstat)::stats	
	type(mesh)::msh	
	type(hobf)::ho_bf,ho_bf_copy
	type(kernelquant)::ker
	type(quant_app),target::quant
	type(proctree)::ptree
	integer,allocatable:: groupmembers(:)	
	integer nmpi
	
	
	! nmpi and groupmembers should be provided by the user 
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo	
	
	call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	
	if(ptree%MyID==Main_ID)write(*,*)'NUMBER_MPI=',nmpi
	
 	threads_num=1
    CALL getenv("OMP_NUM_THREADS", strings)
	strings = TRIM(strings)	
	if(LEN_TRIM(strings)>0)then
		read(strings , *) threads_num
	endif
	if(ptree%MyID==Main_ID)write(*,*)'OMP_NUM_THREADS=',threads_num
	call OMP_set_num_threads(threads_num)
	
		
	call DATE_AND_TIME(values=times)     ! Get the current time 
	seed_myid(1) = times(4) * (360000*times(5) + 6000*times(6) + 100*times(7) + times(8))
	
	! seed_myid(1) = myid*1000
	call RANDOM_SEED(PUT=seed_myid)
	
	! oldmode = vmlsetmode(VML_FTZDAZ_ON)
	! call vmlsetmode(VML_FTZDAZ_ON)
	if(ptree%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "HODLR_BUTTERFLY_SOLVER_3D"
    write(*,*) "   "
	endif
	call InitStat(stats)

	time_tmp = 0
	
 	! register the user-defined function and type in ker 
	ker%FuncZmn=>Z_elem_EMSURF
	ker%QuantZmn=>quant	
	
	msh%Origins=(/0d0,0d0,0d0/)
    msh%integral_points=6
    allocate (msh%ng1(msh%integral_points), msh%ng2(msh%integral_points), msh%ng3(msh%integral_points), msh%gauss_w(msh%integral_points))
    call gauss_points(msh)

	
     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	CALL getarg(1, DATA_DIR)

	option%Nmin_leaf=40
	msh%mesh_normal=1
	! Refined_level=0
	para=0.001
	! levelpara_control=0
	! ACA_tolerance_forward=1d-4
	option%tol_SVD=1d-4
	! SVD_tolerance_factor=1d-4
	option%tol_Rdetect=3d-5
    ! Preset_level_butterfly=0
	msh%scaling=1d0
	quant%wavelength=1.0
	! Discret=0.05
	quant%RCS_static=2
    quant%RCS_Nsample=2000
    ! Optimizing_forward=0
    ! Fast_inverse=0
    ! Add_method_of_base_level=2
    quant%rank_approximate_para1=6.0
    quant%rank_approximate_para2=6.0
    quant%rank_approximate_para3=6.0
	option%tol_LS=1d-12
	! tfqmr_tolerance=1d-6
	option%tol_itersol=3d-3
	option%N_iter=1000
	option%tol_rand=5d-3
	! up_tolerance=1d-4
	! relax_lamda=1d0
	! SolvingMethod=1
	option%level_check=100
	! rank_tmp=7
	! schurinv=1
	! reducelevel_flag=0
	! directschur=1
	option%precon=DIRECT
	! verboselevel=2
	option%xyzsort=3
	option%LnoBP=4000
	option%TwoLayerOnly=0
	quant%CFIE_alpha=1
	option%touch_para=3
    option%schulzorder=3
    option%schulzlevel=3000
	option%LRlevel=100
	option%ErrFillFull=0
	option%RecLR_leaf='A'	
	! call MKL_set_num_threads(NUM_Threads)    ! this overwrites omp_set_num_threads for MKL functions 
	
    ! call OMP_set_dynamic(.true.)
    !call OMP_set_nested(.true.)
	

    !*********************************************************

    quant%omiga=2*pi/quant%wavelength/sqrt(mu0*eps0)
    quant%wavenum=2*pi/quant%wavelength

   !***********************************************************************
	if(ptree%MyID==Main_ID)then							  
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'wavelength:',quant%wavelength
   write (*,*) ''
	endif		
   !***********************************************************************
	
	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID)write(*,*) "geometry modeling......"
    call geo_modeling_SURF(msh,quant,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "modeling finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1

	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "constructing H_matrices formatting......"
    call H_matrix_structuring(ho_bf,para,option,msh,ptree)
	call BPlus_structuring(ho_bf,option,msh,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices formatting finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1
	! stop 
	
    !pause
    
    !call compression_test()
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices filling......"
    call matrices_filling(ho_bf,option,stats,msh,ker,element_Zmn_user,ptree)
	! if(option%precon/=DIRECT)then
		call copy_HOBF(ho_bf,ho_bf_copy)	
	! end if    
	if(ptree%MyID==Main_ID)write(*,*) "H_matrices filling finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	! write(*,*)t2-t1
	
	if(option%precon/=NOPRECON)then								
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing......"
    call cascading_factorizing(ho_bf,option,stats,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	end if

    if(ptree%MyID==Main_ID)write(*,*) "EM_solve......"
    call EM_solve_SURF(ho_bf_copy,ho_bf,option,msh,quant,ptree,stats)
    if(ptree%MyID==Main_ID)write(*,*) "EM_solve finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	

    if(ptree%MyID==Main_ID)write(*,*) "-------------------------------program end-------------------------------------"

    ! ! ! ! pause

end PROGRAM HODLR_BUTTERFLY_SOLVER_3D





subroutine geo_modeling_SURF(msh,quant,ptree)
	use APPLICATION_MODULE
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
	type(quant_app)::quant
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

    if(ptree%MyID==Main_ID)write (*,*) ''
    if(ptree%MyID==Main_ID)write (*,*) 'Maxedge:',Maxedge
	if(ptree%MyID==Main_ID)write (*,*) 'minedgelength:',msh%minedgelength
	if(ptree%MyID==Main_ID)write (*,*) 'wavelength/minedgelength:',quant%wavelength/msh%minedgelength
	if(ptree%MyID==Main_ID)write (*,*) 'maxedgelength:',msh%maxedgelength
	if(ptree%MyID==Main_ID)write (*,*) 'wavelength/maxedgelength:',quant%wavelength/msh%maxedgelength

    if(ptree%MyID==Main_ID)write (*,*) '' 
    
    return
    
end subroutine geo_modeling_SURF




subroutine EM_solve_SURF(ho_bf_for,ho_bf_inv,option,msh,quant,ptree,stats)
    use APPLICATION_MODULE
    use MODULE_FILE
	! use RCS_Bi
	! use RCS_Mono
	! use element_vinc
	use HODLR_Solve
	
    ! use blas95
    implicit none
    
    integer i, j, ii, jj, iii, jjj,ierr
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, N_iter_max, iter ,N_unk, N_unk_loc
    real*8 theta, phi, dphi, rcs_V, rcs_H
    real T0
    real*8 n1,n2,rtemp
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:)
	real*8:: rel_error
	type(Hoption)::option
	type(hobf)::ho_bf_for,ho_bf_inv
	type(mesh)::msh
	type(quant_app)::quant
	type(proctree)::ptree
	type(Hstat)::stats
	complex(kind=8),allocatable:: current(:,:),voltage(:,:)
	
	
	if(option%ErrSol==1)then
		call HODLR_Test_Solve_error(ho_bf_for,ho_bf_inv,option,ptree,stats)
	endif
	
	
	! if(option%PRECON==DIRECT)then
		! msh%idxs = ho_bf_inv%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		! msh%idxe = ho_bf_inv%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)
	! else 
		! msh%idxs = ho_bf_for%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		! msh%idxe = ho_bf_for%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)	
	! endif
	
	! N_unk=msh%Nunk
	N_unk_loc = msh%idxe-msh%idxs+1
	
    if (quant%RCS_static==2) then
    
        theta=90
        phi=0
        
        allocate (current(N_unk_loc,2))
		Current=0
        allocate (voltage(N_unk_loc,2))

		
        !$omp parallel do default(shared) private(edge,value_Z)
        do edge=msh%idxs, msh%idxe
            call element_Vinc_VV_SURF(theta,phi,edge,value_Z,msh,quant)
			voltage(edge-msh%idxs+1,1)=value_Z
			call element_Vinc_HH_SURF(theta,phi,edge,value_Z,msh,quant)
            voltage(edge-msh%idxs+1,2)=value_Z
        enddo    
        !$omp end parallel do
        
        T0=secnds(0.0)
        		
		call HODLR_Solution(ho_bf_for,ho_bf_inv,Current,Voltage,N_unk_loc,2,option,ptree,stats)
					
		
        if(ptree%MyID==Main_ID)write (*,*) ''
        if(ptree%MyID==Main_ID)write (*,*) 'Solving:',secnds(T0),'Seconds'
        if(ptree%MyID==Main_ID)write (*,*) ''

		T0=secnds(0.0)
        call RCS_bistatic_SURF(Current,msh,quant,ptree)
		
		! call current_node_patch_mapping('V',curr(:,1),msh)    		
		! call current_node_patch_mapping('H',curr(:,2),msh)          

        if(ptree%MyID==Main_ID)write (*,*) ''
        if(ptree%MyID==Main_ID)write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
        if(ptree%MyID==Main_ID)write (*,*) ''
		deallocate(Current)
		deallocate(Voltage)
	
    elseif (quant%RCS_static==1) then
    
        allocate (current(N_unk_loc,1))

        
        num_sample=quant%RCS_Nsample
		theta=90.
        dphi=180./num_sample
		allocate (b(N_unk_loc,num_sample+1))
		allocate (x(N_unk_loc,num_sample+1))        
		x=0
		
		
        if(ptree%MyID==Main_ID)open (100, file='bistaticH.out')

        n1=OMP_get_wtime()
		
        do j=0, num_sample 
            phi=j*dphi
			!$omp parallel do default(shared) private(edge,value_Z)
			do edge=msh%idxs, msh%idxe
				call element_Vinc_HH_SURF(theta,phi,edge,value_Z,msh,quant)
				b(edge-msh%idxs+1,j+1)=value_Z
			enddo    
			!$omp end parallel do
        enddo
		
		call HODLR_Solution(ho_bf_for,ho_bf_inv,x,b,N_unk_loc,num_sample+1,option,ptree,stats)
			
		do j=0, num_sample 			
			phi=j*dphi
			
			Current(:,1)=x(:,j+1)
			
            call RCS_monostatic_HH_SURF(theta,phi,rcs_H,Current(:,1),msh,quant,ptree)
!             !$omp parallel do default(shared) private(i)
!             do i=1, N_unk
!                 current(i)=vectors_block(0)%vector(i,2)
!             enddo
!             !$omp end parallel do
!             call RCS_monostatic_HH_SURF(theta,phi,rcs_H)
            
            if(ptree%MyID==Main_ID)write (100,*) phi,rcs_H !,rcs_H
            
            ! deallocate (vectors_block)
            
        enddo
        
		n2 = OMP_get_wtime()
		stats%Time_Sol = stats%Time_Sol + n2-n1
		call MPI_ALLREDUCE(stats%Time_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
		
        if(ptree%MyID==Main_ID)then
			close(100)
			write (*,*) ''
			write (*,*) 'Solving:',rtemp,'Seconds'
			write (*,*) ''   
		endif
		
		call MPI_ALLREDUCE(stats%Flop_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
		if(ptree%MyID==Main_ID)write (*,'(A13Es14.2)') 'Solve flops:',rtemp	
	
		
		deallocate(Current)
		
    endif
        
    return
    
end subroutine EM_solve_SURF



