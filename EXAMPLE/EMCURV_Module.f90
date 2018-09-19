module EMCURV_MODULE
use z_HODLR_DEFS
use z_misc
implicit none

	!**** define your application-related variables here   
	
	!**** quantities related to geometries, meshes, unknowns and points 

	type quant_EMCURV
		real(kind=8) wavenum    ! CEM: wave number  
		real(kind=8) wavelength  ! CEM: wave length
		real(kind=8) omiga       ! CEM: angular frequency
		real(kind=8) rank_approximate_para1, rank_approximate_para2, rank_approximate_para3 ! CEM: rank estimation parameter  
		integer RCS_static  ! CEM: 1: monostatic or 2: bistatic RCS
		integer RCS_Nsample ! CEM: number of RCS samples
		real(kind=8):: CFIE_alpha ! CEM: combination parameter in CFIE
		
		integer Nunk ! size of the matrix 
		real(kind=8),allocatable:: xyz(:,:)   ! coordinates of the points
		integer,allocatable:: info_unk(:,:)  
		! for 2D mesh: 0 point to coordinates of each edge center (unknown x), 1-2 point to coordinates of each edge vertice  

		! 2D mesh
		integer maxnode ! # of vertices in a mesh
		integer maxedge ! # of edges in a mesh 
		real(kind=8) maxedgelength,minedgelength ! maximum and minimum edge length for 2D and 3D meshes	
		integer model2d ! # shape of 2-D curves: (1=strip; 2=corner reflector; 3=two opposite strips; 4=CR with RRS; 5=cylinder; 6=Rectangle Cavity); 7=half cylinder; 8=corrugated half cylinder; 9=corrugated corner reflector; 10=open polygon; 11=taller open polygon
		real(kind=8) Delta_ll	! edge length of each element in 2D curves	
		integer::Ncorner ! # of corners in 2D curves to facilate partitioning
		real(kind=8),allocatable::corner_points(:,:) ! coordinates of corner points 
		real(kind=8)::corner_radius ! radius of the corner points within which no partitioning is performed 
	end type quant_EMCURV

contains

	!**** user-defined subroutine to sample Z_mn
	subroutine Z_elem_EMCURV(ker,m,n,value_e,msh,quant)
		
		use z_HODLR_DEFS
		implicit none
		
		integer edge_m, edge_n, i, j, flag
		integer, INTENT(IN):: m,n
		real(kind=8) r_mn, rtemp1, rtemp2
		complex(kind=8) value_e
		type(z_mesh)::msh
		class(z_kernelquant)::ker
		class(*),pointer :: quant
		
		select TYPE(quant)
			type is (quant_EMCURV)
				! convert to new indices because msh%info_unk has been reordered
				edge_m = m
				edge_n = n
				
				if (edge_m/=edge_n) then
				
					flag=0
					do j=1,2
						do i=1,2
							if (quant%info_unk(i,edge_m)==quant%info_unk(j,edge_n)) then
								flag=1
							endif
						enddo
					enddo
					
					if (flag==1) then
						
						value_e=quant%wavenum*impedence0/4.0*quant%Delta_ll*(1-junit/pi*(3*LOG(3*gamma*quant%wavenum*quant%Delta_ll/4.0)-LOG(gamma*quant%wavenum*quant%Delta_ll/4.0)-2))
					
					else
				
						r_mn=(quant%xyz(1,quant%info_unk(0,edge_m))-quant%xyz(1,quant%info_unk(0,edge_n)))**2+(quant%xyz(2,quant%info_unk(0,edge_m))-quant%xyz(2,quant%info_unk(0,edge_n)))**2
						r_mn=sqrt(r_mn)
						value_e=quant%wavenum*impedence0/4.0*quant%Delta_ll*z_Hankel02_Func(quant%wavenum*r_mn)
					endif
				else
					value_e=quant%wavenum*impedence0/4.0*quant%Delta_ll*(1.0-junit*2.0/pi*(LOG(gamma*quant%wavenum*quant%Delta_ll/4.0)-1.0))
				endif

			class default
				write(*,*)"unexpected type"
				stop
			end select			
			
		return
		
	end subroutine Z_elem_EMCURV	
	

	subroutine VV_polar_CURV(dphi,edge,ctemp_1,curr,msh,quant)
		
		use z_HODLR_DEFS
		! use EMCURV_MODULE
		implicit none
		complex(kind=8)::curr
		complex(kind=8) ctemp,phase,ctemp_1
		real(kind=8) dsita,dphi    
		integer edge
		type(z_mesh)::msh
		type(quant_EMCURV)::quant
		
		phase=junit*quant%wavenum*(quant%xyz(1,quant%info_unk(0,edge))*cos(dphi*pi/180.)+quant%xyz(2,quant%info_unk(0,edge))*sin(dphi*pi/180.))
		ctemp_1=curr*quant%Delta_ll*exp(phase)

		return
		
	end subroutine VV_polar_CURV

	subroutine RCS_bistatic_CURV(curr,msh,quant,ptree)
		!integer flag
		use z_HODLR_DEFS
		! use EMCURV_MODULE
		implicit none
		complex(kind=8)::curr(:)
		real(kind=8) rcs
		complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1,ctemp_loc
		real(kind=8) ddphi,dphi
		
		integer i,j,ii,jj,iii,jjj,patch,flag
		real(kind=8) l_edge,l_edgefine
		type(z_mesh)::msh
		integer edge,edge_m,edge_n,ierr
		type(quant_EMCURV)::quant 
		type(z_proctree)::ptree
		
		if(ptree%MyID==Main_ID)open (100, file='VV_bistatic.txt')
	  
		ddphi=180./quant%RCS_Nsample
		
		do i=0, quant%RCS_Nsample   !phi=0
			dphi=i*ddphi
			ctemp_loc=0
			rcs=0
			!$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)   
			do edge=msh%idxs,msh%idxe
				call VV_polar_CURV(dphi,msh%new2old(edge),ctemp_1,curr(edge-msh%idxs+1),msh,quant)
				ctemp_loc=ctemp_loc+ctemp_1    
			enddo
			!$omp end parallel do
			
			call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
			
			rcs=(abs(impedence0*ctemp))**2/4d0*quant%wavenum
			!rcs=rcs/quant%wavelength
			rcs=10*log10(rcs)
			if(ptree%MyID==Main_ID)write(100,*)dphi,rcs
		enddo
		if(ptree%MyID==Main_ID) close(100)
		
		return
		
	end subroutine RCS_bistatic_CURV



	subroutine RCS_monostatic_VV_CURV(dphi,rcs,curr,msh,quant,ptree)

		use z_HODLR_DEFS
		! use EMCURV_MODULE
		implicit none
		
		real(kind=8) rcs
		complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1,ctemp_2
		real(kind=8) dsita,dphi
		integer edge,edge_m,edge_n,ierr
		complex(kind=8):: curr(:)
		type(z_mesh)::msh
		type(quant_EMCURV)::quant
		type(z_proctree)::ptree
		
		ctemp_loc=0
			rcs=0
			!$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
			do edge=msh%idxs,msh%idxe
				call VV_polar_CURV(dphi,msh%new2old(edge),ctemp_1,curr(edge-msh%idxs+1),msh,quant)
				ctemp_loc=ctemp_loc+ctemp_1
			enddo
			!$omp end parallel do
			
			call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
			
			rcs=(abs(impedence0*ctemp))**2/4d0*quant%wavenum
			!rcs=rcs/quant%wavelength
			rcs=10*log10(rcs)
			
		return
		
	end subroutine RCS_monostatic_VV_CURV

	subroutine element_Vinc_VV_CURV(phi,edge,value,msh,quant)

		use z_HODLR_DEFS
		! use EMCURV_MODULE
		implicit none
		
		integer edge
		complex(kind=8) value
		real(kind=8) theta, phi
		complex(kind=8)  phase
		type(z_mesh)::msh
		type(quant_EMCURV)::quant
		
		phase=junit*quant%wavenum*(quant%xyz(1,quant%info_unk(0,edge))*cos(phi*pi/180.)+quant%xyz(2,quant%info_unk(0,edge))*sin(phi*pi/180.))
		value=exp(phase)
		
		return
		
	end subroutine element_Vinc_VV_CURV	
		

end module EMCURV_MODULE	