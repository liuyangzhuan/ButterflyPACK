module gauss_rule
use misc
contains

!***********************************
!	seven points gauss integration
!	provide w(n) and x,y,z(n) (n=7)
!***********************************
  subroutine gau_grobal(nn,j,x,y,z,w)
  
  use MODULE_FILE
  implicit none

  integer nn ,flag
  real*8 x(integral_points),y(integral_points),z(integral_points),w(integral_points)
  integer i,j,ii
!	wa=area*9./40
!	wb=area*(155.-sqrt(15.))/1200.
!	wc=area*(155.+sqrt(15.))/1200.

   ! if(flag==1) then

!		ii=node_patch_of_edgefine(j,nn)
!  do i=1,integral_points
!	x(i)=ng1(i)*xyz(1,node_of_patchfine(1,ii))+ng2(i)*xyz(1,node_of_patchfine(2,ii))+&
!     ng3(i)*xyz(1,node_of_patchfine(3,ii))
!	y(i)=ng1(i)*xyz(2,node_of_patchfine(1,ii))+ng2(i)*xyz(2,node_of_patchfine(2,ii))+&
!     ng3(i)*xyz(2,node_of_patchfine(3,ii))
!	z(i)=ng1(i)*xyz(3,node_of_patchfine(1,ii))+ng2(i)*xyz(3,node_of_patchfine(2,ii))+&
!     ng3(i)*xyz(3,node_of_patchfine(3,ii))
!  enddo
  
!  elseif(flag==2) then

  ii=node_patch_of_edge(j,nn)
	 do i=1,integral_points
	x(i)=ng1(i)*xyz(1,node_of_patch(1,ii))+ng2(i)*xyz(1,node_of_patch(2,ii))+&
     ng3(i)*xyz(1,node_of_patch(3,ii))
	y(i)=ng1(i)*xyz(2,node_of_patch(1,ii))+ng2(i)*xyz(2,node_of_patch(2,ii))+&
     ng3(i)*xyz(2,node_of_patch(3,ii))
	z(i)=ng1(i)*xyz(3,node_of_patch(1,ii))+ng2(i)*xyz(3,node_of_patch(2,ii))+&
     ng3(i)*xyz(3,node_of_patch(3,ii))
  enddo
  
  w=gauss_w

! 	w(1)=wa
! 	w(2)=wb
! 	w(3)=wb
! 	w(4)=wb
! 	w(5)=wc
! 	w(6)=wc
! 	w(7)=wc
  return
  end subroutine gau_grobal
  

  subroutine gauss_points()
  
      use MODULE_FILE
      implicit none
      
      real*8 v1,v2,v3,v4,v5
      
    !	wa=area*9./40
!	wb=area*(155.-sqrt(15.))/1200.
!	wc=area*(155.+sqrt(15.))/1200.
    if (integral_points==7) then
    
    wa=9./80
	wb=(155.-sqrt(15.))/2400.
	wc=(155.+sqrt(15.))/2400.

	v1=1./3.
	v2=(6.-sqrt(15.))/21.
	v3=(9.+2*sqrt(15.))/21.
	v4=(6.+sqrt(15.))/21.
	v5=(9.-2*sqrt(15.))/21.

	ng1(1)=1-v1-v1
	ng1(2)=1-v2-v2
	ng1(3)=1-v2-v3
	ng1(4)=1-v3-v2
	ng1(5)=1-v4-v4
	ng1(6)=1-v4-v5
	ng1(7)=1-v5-v4

	ng2(1)=v1
	ng2(2)=v2
	ng2(3)=v2
	ng2(4)=v3
	ng2(5)=v4
	ng2(6)=v4
	ng2(7)=v5

	ng3(1)=v1
	ng3(2)=v2
	ng3(3)=v3
	ng3(4)=v2
	ng3(5)=v4
	ng3(6)=v5
	ng3(7)=v4
	
	gauss_w(1)=wa
	gauss_w(2)=wb
	gauss_w(3)=wb
	gauss_w(4)=wb
	gauss_w(5)=wc
	gauss_w(6)=wc
	gauss_w(7)=wc
	
	elseif (integral_points==6) then
	
	v1=0.816847572980459d0
    v2=0.091576213509771d0
    v3=0.108103018168070d0
    v4=0.445948490915965d0
    wa=0.109951743655322d0/2.
    wb=0.223381589678011d0/2.

	ng1(1)=v1
	ng1(2)=v2
	ng1(3)=v2
	ng1(4)=v3
	ng1(5)=v4
	ng1(6)=v4
	!ng1(7)=1-v5-v4

	ng2(1)=v2
	ng2(2)=v1
	ng2(3)=v2
	ng2(4)=v4
	ng2(5)=v3
	ng2(6)=v4
	!ng2(7)=v5

	ng3(1)=v2
	ng3(2)=v2
	ng3(3)=v1
	ng3(4)=v4
	ng3(5)=v4
	ng3(6)=v3
	
	gauss_w(1)=wa
	gauss_w(2)=wa
	gauss_w(3)=wa
	gauss_w(4)=wb
	gauss_w(5)=wb
	gauss_w(6)=wb
	
	elseif (integral_points==4) then
	
	v1=1./3.
	v2=0.2
	v3=0.6
	wa=-27./96.
	wb=25./96.
	
	ng1(1)=v1
	ng1(2)=v2
	ng1(3)=v2
	ng1(4)=v3
	
	ng2(1)=v1
	ng2(2)=v2
	ng2(3)=v3
	ng2(4)=v2
	
	ng3(1)=v1
	ng3(2)=v3
	ng3(3)=v2
	ng3(4)=v2
	
	gauss_w(1)=wa
	gauss_w(2)=wb
	gauss_w(3)=wb
	gauss_w(4)=wb
	
	endif
    
    return
  end subroutine gauss_points
  
  end module gauss_rule