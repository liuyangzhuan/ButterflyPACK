module RCS_Mono
use RCS_Bi
contains

subroutine RCS_monostatic_VV_SURF(dsita,dphi,rcs,curr,msh,ker,ptree)

    use MODULE_FILE
    implicit none
    complex(kind=8)::curr(:)
    real*8 rcs
    complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1,ctemp_2
    real*8 dsita,dphi
    integer edge,edge_m,edge_n,ierr
	type(mesh)::msh
    type(kernelquant)::ker
	type(proctree)::ptree
	
    ctemp_loc=(0.,0.)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
        do edge=msh%idxs,msh%idxe
            call VV_polar_SURF(dsita,dphi,edge,ctemp_1,curr(edge-msh%idxs+1),msh,ker)
	        ctemp_loc=ctemp_loc+ctemp_1
        enddo
        !$omp end parallel do
		
		call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
		
        rcs=(abs(ker%wavenum*ctemp))**2/4/pi
        !rcs=rcs/ker%wavelength
        rcs=10*log10(rcs)
        
    return
    
end subroutine RCS_monostatic_VV_SURF

subroutine RCS_monostatic_HH_SURF(dsita,dphi,rcs,curr,msh,ker,ptree)

    use MODULE_FILE
    implicit none
    
    real*8 rcs
    complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1,ctemp_2
    real*8 dsita,dphi    
    integer edge,edge_m,edge_n,ierr
    complex(kind=8)::curr(:)
	type(mesh)::msh
	type(kernelquant)::ker
	type(proctree)::ptree
	
    ctemp_loc=(0.,0.)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
        do edge=msh%idxs,msh%idxe
            call HH_polar_SURF(dsita,dphi,edge,ctemp_1,curr(edge-msh%idxs+1),msh,ker)
	        ctemp_loc=ctemp_loc+ctemp_1
        enddo
        !$omp end parallel do
		
		call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
		
        rcs=(abs(ker%wavenum*ctemp))**2/4/pi
        !rcs=rcs/ker%wavelength
        rcs=10*log10(rcs)
        
    return
    
end subroutine RCS_monostatic_HH_SURF


subroutine RCS_monostatic_VV_CURV(dphi,rcs,curr,msh,ker,ptree)

    use MODULE_FILE
    implicit none
    
    real*8 rcs
    complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1,ctemp_2
    real*8 dsita,dphi
    integer edge,edge_m,edge_n,ierr
    complex(kind=8):: curr(:)
	type(mesh)::msh
	type(kernelquant)::ker
	type(proctree)::ptree
	
    ctemp_loc=0
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
        do edge=msh%idxs,msh%idxe
            call VV_polar_CURV(dphi,edge,ctemp_1,curr(edge-msh%idxs+1),msh,ker)
	        ctemp_loc=ctemp_loc+ctemp_1
        enddo
        !$omp end parallel do
		
		call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
		
        rcs=(abs(impedence0*ctemp))**2/4d0*ker%wavenum
        !rcs=rcs/ker%wavelength
        rcs=10*log10(rcs)
        
    return
    
end subroutine RCS_monostatic_VV_CURV



end module RCS_Mono