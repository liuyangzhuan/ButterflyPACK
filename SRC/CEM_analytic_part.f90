module analytic_part
use misc

contains 


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

ii=msh%info_unk(jj,mm)  !ii 表示面
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
s(2,1)=((u3-u0)*(u3-l(3))+v3*(v3-v0))/l(1)  !去掉了个负号
s(1,2)=-(u3*(u3-u0)+v3*(v3-v0))/l(2)
s(2,2)=(u0*u3+v0*v3)/l(2)     !改过了的
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
s(2,1)=((u3-u0)*(u3-l(3))+v3*(v3-v0))/l(1)  !去掉了个负号
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

end module 


