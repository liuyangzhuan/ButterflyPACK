module element_Z
use gauss_rule
use analytic_part
contains 





subroutine element_Zmn(edge_m,edge_n,value_e)
    use MODULE_FILE
    implicit none
	integer edge_m, edge_n
    complex(kind=8) value_e
	if(Kernel==RBF)then
		call element_Zmn_RBF(edge_m,edge_n,value_e)
	else if(Kernel==EMCURV)then
		call element_Zmn_EMCURV(edge_m,edge_n,value_e)
	else if(Kernel==EMSURF)then
		call element_Zmn_EMSURF(edge_m,edge_n,value_e)
	else if(Kernel==FULL)then
		call element_Zmn_FULL(edge_m,edge_n,value_e)		
	endif
	
end subroutine element_Zmn



subroutine element_Zmn_RBF(edge_m,edge_n,value_e)
    
    use MODULE_FILE
    implicit none
    
    integer edge_m, edge_n, i, j, flag
    real(kind=8) r_mn, rtemp1, rtemp2
    complex(kind=8) value_e
    integer dimn
	dimn = size(xyz,1)
	
	value_e=0
	
	r_mn=sum((xyz(1:dimn,node_patch_of_edge(0,edge_m))-xyz(1:dimn,node_patch_of_edge(0,edge_n)))**2)
	
	value_e = exp(-r_mn/2.0/sigma**2)
	if(r_mn==0)value_e = value_e + lambda

	! if(abs(value_e)<1D-10)write(*,*)r_mn,-r_mn/2.0/sigma**2,exp(-r_mn/2.0/sigma**2),'wocao'
	
	! r_mn=sqrt(r_mn)
	!!value_e=wavenum*impedence/4.0*Delta_ll*Hankel02_Func(wavenum*r_mn)	
	return
    
end subroutine element_Zmn_RBF


subroutine element_Zmn_Full(edge_m,edge_n,value)
    
    use MODULE_FILE
    implicit none
    
    integer edge_m,edge_n
    complex(kind=8) value
    
	value = matZ_glo(new2old(edge_m),new2old(edge_n))
	

end subroutine element_Zmn_Full




subroutine element_Zmn_EMCURV(edge_m,edge_n,value_e)
    
    use MODULE_FILE
    implicit none
    
    integer edge_m, edge_n, i, j, flag
    real*8 r_mn, rtemp1, rtemp2
    complex(kind=8) value_e
       
    if (edge_m/=edge_n) then
    
        flag=0
        do j=1,2
            do i=1,2
                if (node_patch_of_edge(i,edge_m)==node_patch_of_edge(j,edge_n)) then
                    flag=1
                endif
            enddo
        enddo
        
        if (flag==1) then
			! write(*,*)Delta_ll,wavenum,impedence,junit,gamma,'ddd' 
            value_e=wavenum*impedence/4.0*Delta_ll*(1-junit/pi*(3*LOG(3*gamma*wavenum*Delta_ll/4.0)-LOG(gamma*wavenum*Delta_ll/4.0)-2))
        
        else
    
            r_mn=(xyz(1,node_patch_of_edge(0,edge_m))-xyz(1,node_patch_of_edge(0,edge_n)))**2+(xyz(2,node_patch_of_edge(0,edge_m))-xyz(2,node_patch_of_edge(0,edge_n)))**2
            r_mn=sqrt(r_mn)
            value_e=wavenum*impedence/4.0*Delta_ll*Hankel02_Func(wavenum*r_mn)
        endif
    else
        value_e=wavenum*impedence/4.0*Delta_ll*(1.0-junit*2.0/pi*(LOG(gamma*wavenum*Delta_ll/4.0)-1.0))
    endif

	return
    
end subroutine element_Zmn_EMCURV


subroutine element_Zmn_EMSURF(edge_m,edge_n,value_e)
    
    use MODULE_FILE
    implicit none
    
    integer flag,edge_m,edge_n,nodetemp_n,patch
    complex(kind=8) value_e,value_m,value
    integer i,j,ii,jj,iii,jjj
    real*8 ln,lm,am(3),an(3),xm(integral_points),ym(integral_points),zm(integral_points),wm(integral_points),nr_m(3)
    real*8 xn(integral_points),yn(integral_points),zn(integral_points),wn(integral_points),nr_n(3)
    complex(kind=8) ctemp,ctemp1,ctemp2,aa(3),bb(1),dg(3),dg1(3),dg2(3)
    complex(kind=8) imp,imp1,imp2,imp3
    real*8 temp
    real*8 distance
    real*8 ianl,ianl1,ianl2
    real*8 area
    
    lm=sqrt((xyz(1,node_patch_of_edge(1,edge_m))-xyz(1,node_patch_of_edge(2,edge_m)))**2+(xyz(2,node_patch_of_edge(1,edge_m))-xyz(2,node_patch_of_edge(2,edge_m)))**2+(xyz(3,node_patch_of_edge(1,edge_m))-xyz(3,node_patch_of_edge(2,edge_m)))**2)
    ln=sqrt((xyz(1,node_patch_of_edge(1,edge_n))-xyz(1,node_patch_of_edge(2,edge_n)))**2+(xyz(2,node_patch_of_edge(1,edge_n))-xyz(2,node_patch_of_edge(2,edge_n)))**2+(xyz(3,node_patch_of_edge(1,edge_n))-xyz(3,node_patch_of_edge(2,edge_n)))**2)
    
    ctemp1=(0.,0.)
    ctemp2=(0.,0.)
    value_m=(0.,0.)
    
    do ii=3,4
        call gau_grobal(edge_m,ii,xm,ym,zm,wm)                        	       
        nr_m(1:3)=normal_of_patch(1:3,node_patch_of_edge(ii,edge_m))
        do i=1,integral_points
            am(1)=xm(i)-xyz(1,node_patch_of_edge(ii+2,edge_m))
            am(2)=ym(i)-xyz(2,node_patch_of_edge(ii+2,edge_m))
            am(3)=zm(i)-xyz(3,node_patch_of_edge(ii+2,edge_m))
            bb(1)=(0.,0.)
            aa(1:3)=(0.,0.)
            do jj=3,4                
                call gau_grobal(edge_n,jj,xn,yn,zn,wn)
                nr_n(1:3)=normal_of_patch(1:3,node_patch_of_edge(jj,edge_n))
                if (node_patch_of_edge(ii,edge_m)==node_patch_of_edge(jj,edge_n)) then
                    area=triangle_area(node_patch_of_edge(ii,edge_m))
                    imp=(0.,0.);imp1=(0.,0.);imp2=(0.,0.);imp3=(0.,0.)
                    an(1)=xm(i)-xyz(1,node_patch_of_edge(jj+2,edge_n))
                    an(2)=ym(i)-xyz(2,node_patch_of_edge(jj+2,edge_n))
			        an(3)=zm(i)-xyz(3,node_patch_of_edge(jj+2,edge_n))
			        call scalar(am,an,temp)    
                    value_m=value_m+(-1)**(ii+1)*(-1)**(jj+1)*0.5*temp/(2.*area)*wm(i)
                    do j=1,integral_points
                        distance=sqrt((xm(i)-xn(j))**2+(ym(i)-yn(j))**2+(zm(i)-zn(j))**2)			                                     
                        if(distance==0)then    
                            imp=imp+wn(j)*(-junit*wavenum)
                            imp1=imp1+ng1(j)*wn(j)*(-junit*wavenum)
                            imp2=imp2+ng2(j)*wn(j)*(-junit*wavenum)
                            ianl=ianalytic(edge_n,jj,xn(j),yn(j),zn(j))
                            ianl1=ianalytic2(edge_n,jj,xn(j),yn(j),zn(j),1)
                            ianl2=ianalytic2(edge_n,jj,xn(j),yn(j),zn(j),2)
                            imp=imp+ianl  !debug
                            imp1=imp1+ianl1
                            imp2=imp2+ianl2
                        else
                            imp=imp+wn(j)*exp(-junit*wavenum*distance)/distance
                            imp1=imp1+ng1(j)*wn(j)*exp(-junit*wavenum*distance)/distance
                            imp2=imp2+ng2(j)*wn(j)*exp(-junit*wavenum*distance)/distance
                        endif                        
                    enddo
                    imp3=imp-imp1-imp2
                    nodetemp_n=node_patch_of_edge(jj+2,edge_n)
                    patch=node_patch_of_edge(jj,edge_n)
                    do jjj=1,3
                        aa(jjj)=aa(jjj)+(-1)**(jj+1)*wavenum**2*(xyz(jjj,node_of_patch(1,patch))*imp1+xyz(jjj,node_of_patch(2,patch))*imp2+xyz(jjj,node_of_patch(3,patch))*imp3-xyz(jjj,nodetemp_n)*imp)          
                    enddo
                    bb(1)=bb(1)+(-1)**(jj+1)*imp
                else
                    imp=(0.,0.);imp1=(0.,0.);imp2=(0.,0.);imp3=(0.,0.)
                    do j=1,integral_points
                        distance=sqrt((xm(i)-xn(j))**2+(ym(i)-yn(j))**2+(zm(i)-zn(j))**2)
                        an(1)=xn(j)-xyz(1,node_patch_of_edge(jj+2,edge_n))
                        an(2)=yn(j)-xyz(2,node_patch_of_edge(jj+2,edge_n))
			            an(3)=zn(j)-xyz(3,node_patch_of_edge(jj+2,edge_n))
			            dg(1)=(xm(i)-xn(j))*(1+junit*wavenum*distance)*exp(-junit*wavenum*distance)/(4*pi*distance**3)
                        dg(2)=(ym(i)-yn(j))*(1+junit*wavenum*distance)*exp(-junit*wavenum*distance)/(4*pi*distance**3)
                        dg(3)=(zm(i)-zn(j))*(1+junit*wavenum*distance)*exp(-junit*wavenum*distance)/(4*pi*distance**3)
			             call ccurl(an,dg,dg1)
			             call ccurl(nr_m,dg1,dg2)
			             call cscalar(dg2,am,ctemp)
			             value_m=value_m-(-1)**(ii+1)*(-1)**(jj+1)*ctemp*wm(i)*wn(j)			             
                        imp=imp+wn(j)*exp(-junit*wavenum*distance)/distance
                        imp1=imp1+ng1(j)*wn(j)*exp(-junit*wavenum*distance)/distance
                        imp2=imp2+ng2(j)*wn(j)*exp(-junit*wavenum*distance)/distance
                    enddo                        
                    imp3=imp-imp1-imp2
                    nodetemp_n=node_patch_of_edge(jj+2,edge_n)
                    patch=node_patch_of_edge(jj,edge_n)
                    do jjj=1,3
                       aa(jjj)=aa(jjj)+(-1)**(jj+1)*wavenum**2*(xyz(jjj,node_of_patch(1,patch))*imp1+xyz(jjj,node_of_patch(2,patch))*imp2+xyz(jjj,node_of_patch(3,patch))*imp3-xyz(jjj,nodetemp_n)*imp)          
                    enddo
                    bb(1)=bb(1)+(-1)**(jj+1)*imp                                                                                                    
                endif
            enddo	
            call cscalar(aa,am,ctemp)
            ctemp1=ctemp1+(-1)**(ii+1)*ctemp*wm(i) 
            ctemp2=ctemp2+4.*(-1)**(ii+1)*bb(1)*wm(i)    
        enddo
    enddo
    value_e=ln*lm*junit*(ctemp1-ctemp2)/4./pi/omiga/eps0
    value_m=value_m*lm*ln
	
    value=CFIE_alpha*value_e+(1.-CFIE_alpha)*impedence*value_m
    
    return
    
end subroutine element_Zmn_EMSURF

    
real*8 function triangle_area(patch)
    
    use MODULE_FILE
    implicit none
    
    integer patch,i
    real*8 a(3),b(3),c(3)
    
    do i=1,3
        a(i)=xyz(i,node_of_patch(2,patch))-xyz(i,node_of_patch(1,patch))
        b(i)=xyz(i,node_of_patch(3,patch))-xyz(i,node_of_patch(1,patch))
    enddo
    
    call curl(a,b,c)
    triangle_area=0.5*sqrt(c(1)**2+c(2)**2+c(3)**2)
    
    return
end function triangle_area
    

real*8 function arg_thresh_Zmn()
use MODULE_FILE
implicit none 
	if(Kernel==RBF)then
		! write(*,*)-log(SafeUnderflow)*2.0*sigma**2 ,'ri'
		arg_thresh_Zmn=-log(SafeUnderflow)*2.0*sigma**2 
	else 
		arg_thresh_Zmn=1D300
	endif
end function arg_thresh_Zmn


complex(kind=8) function Hankel02_Func(x)

use MODULE_FILE    
implicit none
    
    real*8 x
    complex(kind=8) y
    
    Hankel02_Func=BesselJ0_func(x)-junit*BesselY0_func(x)
    
    return
    
end function Hankel02_Func

real*8 function BesselJ0_func(x)

    implicit none
    
    real*8 x, z, ax
    real*8 y, rtemp1, rtemp2, xx
    
    ax=abs(x)
    if (ax<8d0) then
        y=x*x
        rtemp1=57568490574.0d0+y*(-13362690354.0d0+y*(651619640.7d0+y*(-11214424.18d0+y*(77392.33017d0+y*(-184.9052456d0)))))
        rtemp2=57568490411.0d0+y*(1029532985.0d0+y*(9494680.718d0+y*(59272.64853d0+y*(267.8532712d0+y*1.0d0))))
        BesselJ0_func=rtemp1/rtemp2
    else
        z=8.0d0/ax
        y=z*z

        xx=ax-0.785398164d0
        rtemp1=1.0d0+y*(-0.1098628627d-2+y*(0.2734510407d-4+y*(-0.2073370639d-5+y*0.2093887211d-6)))
        rtemp2=-0.1562499995d-1+y*(0.1430488765d-3+y*(-0.6911147651d-5+y*(0.7621095161d-6-y*0.934935152d-7)))
        BesselJ0_func=sqrt(0.636619772d0/ax)*(cos(xx)*rtemp1-z*sin(xx)*rtemp2)
    endif
    
    return
    
end function BesselJ0_func

real*8 function BesselY0_func(x)

    implicit none
    
    real*8 x, z, ax
    real*8 y, rtemp1, rtemp2, xx
    
    if (x<8.0d0) then
        y=x*x
        rtemp1=-2957821389.0d0+y*(7062834065.0d0+y*(-512359803.6d0+y*(10879881.29d0+y*(-86327.92757d0+y*228.4622733d0))))
        rtemp2=40076544269.0d0+y*(745249964.8d0+y*(7189466.438d0+y*(47447.26470d0+y*(226.1030244d0+y*1.0d0))))
        BesselY0_func=(rtemp1/rtemp2)+0.636619772d0*BesselJ0_func(x)*LOG(x)
    else
        z=8.0d0/x
        y=z*z

        xx=x-0.785398164d0
        rtemp1=1.0d0+y*(-0.1098628627d-2+y*(0.2734510407d-4+y*(-0.2073370639d-5+y*0.2093887211d-6)))
        rtemp2=-0.1562499995d-1+y*(0.1430488765d-3+y*(-0.6911147651d-5+y*(0.7621095161d-6-y*0.934935152d-7)))
        BesselY0_func=sqrt(0.636619772d0/x)*(sin(xx)*rtemp1+z*cos(xx)*rtemp2)
    endif
    
    return
    
end function BesselY0_func
    

end module element_Z
    