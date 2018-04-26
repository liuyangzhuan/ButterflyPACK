module RCS_Mono
use RCS_Bi
contains

subroutine RCS_monostatic_VV_SURF(dsita,dphi,rcs)

    use MODULE_FILE
    implicit none
    
    real*8 rcs
    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1,ctemp_2
    real*8 dsita,dphi
    integer edge,edge_m,edge_n
    
    ctemp=(0.,0.)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp)
        do edge=1,maxedge
            call VV_polar_SURF(dsita,dphi,edge,ctemp_1)
	        ctemp=ctemp+ctemp_1
        enddo
        !$omp end parallel do
        rcs=(abs(wavenum*ctemp))**2/4/pi
        !rcs=rcs/wavelength
        rcs=10*log10(rcs)
        
    return
    
end subroutine RCS_monostatic_VV_SURF

subroutine RCS_monostatic_HH_SURF(dsita,dphi,rcs)

    use MODULE_FILE
    implicit none
    
    real*8 rcs
    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1,ctemp_2
    real*8 dsita,dphi    
    integer edge,edge_m,edge_n
    
    ctemp=(0.,0.)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp)
        do edge=1,maxedge
            call HH_polar_SURF(dsita,dphi,edge,ctemp_1)
	        ctemp=ctemp+ctemp_1
        enddo
        !$omp end parallel do
        rcs=(abs(wavenum*ctemp))**2/4/pi
        !rcs=rcs/wavelength
        rcs=10*log10(rcs)
        
    return
    
end subroutine RCS_monostatic_HH_SURF


subroutine RCS_monostatic_VV_CURV(dphi,rcs)

    use MODULE_FILE
    implicit none
    
    real*8 rcs
    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1,ctemp_2
    real*8 dsita,dphi
    integer edge,edge_m,edge_n
    
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
        
    return
    
end subroutine RCS_monostatic_VV_CURV



end module RCS_Mono