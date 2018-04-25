module element_Z
use gauss_rule
use analytic_part
contains 

subroutine element_Zmn(edge_m,edge_n,value_e)
    
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
    !value_m=(0.,0.)
    
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
                !nr_n(1:3)=normal_of_patch(1:3,node_patch_of_edge(jj,edge_n))
                if (node_patch_of_edge(ii,edge_m)==node_patch_of_edge(jj,edge_n)) then
                    !area=triangle_area(node_patch_of_edge(ii,edge_m))
                    imp=(0.,0.);imp1=(0.,0.);imp2=(0.,0.);imp3=(0.,0.)
                    an(1)=xm(i)-xyz(1,node_patch_of_edge(jj+2,edge_n))
                    an(2)=ym(i)-xyz(2,node_patch_of_edge(jj+2,edge_n))
			        an(3)=zm(i)-xyz(3,node_patch_of_edge(jj+2,edge_n))
			        !call scalar(am,an,temp)    
                    !value_m=value_m+(-1)**(ii+1)*(-1)**(jj+1)*0.5*temp/(2.*area)*wm(i)
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
			            !dg(1)=(xm(i)-xn(j))*(1+junit*wavenum*distance)*exp(-junit*wavenum*distance)/(4*pi*distance**3)
                        !dg(2)=(ym(i)-yn(j))*(1+junit*wavenum*distance)*exp(-junit*wavenum*distance)/(4*pi*distance**3)
                        !dg(3)=(zm(i)-zn(j))*(1+junit*wavenum*distance)*exp(-junit*wavenum*distance)/(4*pi*distance**3)
			             !call ccurl(an,dg,dg1)
			             !call ccurl(nr_m,dg1,dg2)
			             !call cscalar(dg2,am,ctemp)
			             !value_m=value_m-(-1)**(ii+1)*(-1)**(jj+1)*ctemp*wm(i)*wn(j)			             
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
    !value_m=value_m*lm*ln
    
    !value=alpha*value_e+(1.-alpha)*impedence*value_m
    
    return
    
end subroutine element_Zmn

    
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
    
    

end module element_Z
    