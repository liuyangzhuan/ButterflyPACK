! “ButterflyPACK” Copyright (c) 2018, The Regents of the University of California, through
! Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.

! If you have questions about your rights to use or distribute this software, please contact
! Berkeley Lab's Intellectual Property Office at  IPO@lbl.gov.

! NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the
! U.S. Government consequently retains certain rights. As such, the U.S. Government has been
! granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable
! worldwide license in the Software to reproduce, distribute copies to the public, prepare
! derivative works, and perform publicly and display publicly, and to permit other to do so.

! Developers: Yang Liu
!             (Lawrence Berkeley National Lab, Computational Research Division).

! This exmple works with double-complex precision data
#define DAT 0

#include "ButterflyPACK_config.fi"


module EMSURF_PORT_MODULE
use BPACK_DEFS
use MISC_Utilities
implicit none

	!**** define your application-related variables here


	real(kind=8):: r_TE_nm(3,3)	=	&
	reshape((/3.8317,	1.8412,		3.0542,		&
				7.0156,	5.3314,		6.7061,		&
				10.1735,	8.5363,	 	9.9695	/), (/3,3/))  ! roots of derivative of bessel j
	real(kind=8):: r_TM_nm(3,3)=	&
	reshape((/2.4048,	3.8317,		5.1356,		&
				5.5201,	7.0156,		8.4172,		&
				8.6537,	10.1735,	11.6198	/), (/3,3/))  ! roots of bessel j

	real(kind=8):: A_TE_nm_cir(3,3)	=	&
	reshape((/1.40081,	1.63313,		3.32719,		&
			1.87991,	2.34682,		3.18421,		&
			2.25943,	2.93968,	 	3.59894	/), (/3,3/)) ! normalization factor of TE_nm circular: A_TE_nm_cir/R
	real(kind=8):: A_TM_nm_cir(3,3) =	&
	reshape((/1.08676,	1.98104,		2.96082,		&
			1.65809,	2.65859,		3.38125,		&
			2.07841,	3.19531,	 	3.79779	/), (/3,3/))  ! normalization factor of TM_nm circular A_TM_nm_cir/R

	! normalization factor of TE_nm, rectangular: 1/sqrt(a/b*m^2*A_nm_rec + b/a*n^2*B_nm_rec)
	! normalization factor of TM_nm, rectangular: 1/sqrt(a/b*m^2*B_nm_rec + b/a*n^2*A_nm_rec)
	real(kind=8):: A_nm_rec(3,3)	=	&
	reshape((/0d0,	0d0,			0d0,		&
			0.5d0,	0.125d0,		0.25d0,		&
			0.5d0,	0.25d0,	 		0.125d0	/), (/3,3/))

	real(kind=8):: B_nm_rec(3,3)	=	&
	reshape((/0d0,		0.5d0,		0.5d0,		&
			0d0,		0.125d0,	0.25d0,		&
			0d0,		0.25d0,	 	0.125d0	/), (/3,3/))

	!**** edges of each node
	type edge_node
		integer Nedge
		integer,allocatable::edges(:,:)
	end type edge_node

	!**** define a port
	type port
		integer type   ! 0: circular 1: rectangular
		real(kind=8) origin(3)  ! origin of the local coordinate system
		real(kind=8) x(3),y(3),z(3)   ! local coordinate units
		integer:: Nunk=0            ! number of internal edges on the port
		real(kind=8) R  ! radius of the port if it's circular
		real(kind=8) a, b ! sizes of the port if it's rectangular
		integer:: nmax=1   ! TM_nm and TE_nm
		integer:: mmax=1   ! TM_nm and TE_nm
		real(kind=8):: A_TE_nm(3,3) ! normalization factor of TE_nm
		real(kind=8):: A_TM_nm(3,3) ! normalization factor of TM_nm
		real(kind=8):: impedance_TE_nm(3,3)  ! mode impedance TE_nm
		real(kind=8):: impedance_TM_nm(3,3)  ! mode impedance of TM_nm
		complex(kind=8),allocatable::nxe_dot_rwg(:,:,:,:,:)  ! int_nxe_dot_rwg of shape Nunk * nmax+1 * mmax * 2 * npolar, the third dimension differentiate TM and TE, the last represents polarization degeneracy of circular waveguide
		complex(kind=8),allocatable::e_dot_rwg(:,:,:,:,:)  ! int_e_dot_rwg of shape Nunk * nmax+1 * mmax * 2 * npolar, the third dimension differentiate TM and TE, the last represents polarization degeneracy of circular waveguide
	end type port




	!**** quantities related to geometries, meshes, unknowns and points
	type quant_EMSURF
		real(kind=8) wavenum    ! CEM: wave number
		real(kind=8) wavelength  ! CEM: wave length
		real(kind=8) freq       ! CEM: frequency
		real(kind=8) rank_approximate_para1, rank_approximate_para2, rank_approximate_para3 ! CEM: rank estimation parameter
		integer RCS_static  ! CEM: 1: monostatic or 2: bistatic RCS
		integer RCS_Nsample ! CEM: number of RCS samples
		real(kind=8):: CFIE_alpha ! CEM: combination parameter in CFIE
		real(kind=8):: sigma_s=1/1.7d-8  ! surface conductivity, copper 1/1.7d-8

		integer Nunk ! size of the matrix
		integer Nunk_int ! not port unknowns
		integer Nunk_port ! port unknowns
		real(kind=8),allocatable:: xyz(:,:)   ! coordinates of the points
		integer,allocatable:: info_unk(:,:)
		! for 3D mesh: 0 point to coordinates of each edge center (unknown x), 1-2 point to coordinates of each edge vertice, 3-4 point to two patches that share the same edge, 5-6 point to coordinates of last vertice of the two patches

		integer:: Nport=0 ! number of ports
		integer:: noport=0 ! whether to treat the cavity as closed cavity
		type(port), allocatable:: ports(:)   ! the open ports

		integer:: postprocess = 1 ! whether postprocessing is carried out
		integer:: Nobs = 0 ! number of observation points (currently cannot locate on the walls or ports)
		real(kind=8), allocatable:: obs_points(:,:) ! xyz coordinates, dimensions 3xNobs
		complex(kind=8), allocatable:: obs_Efields(:,:) ! E fields, dimensions 3xNobs

		! 3D mesh
		integer maxnode ! # of vertices in a mesh
		integer maxedge ! # of edges in a mesh
		real(kind=8) maxedgelength,minedgelength ! maximum and minimum edge length for 2D and 3D meshes
		integer integral_points ! #of Gauss quadrature points on a triangular
		integer maxpatch ! # of triangular patches
		integer mesh_normal	 ! flags to flip the unit normal vectors of a triangular patch
		real(kind=8) scaling  ! scaling factor of the coordinates of vertices of a 3D mesh
		real(kind=8), allocatable :: ng1(:),ng2(:),ng3(:),gauss_w(:) ! Gass quadrature and weights
		real(kind=8),allocatable:: normal_of_patch(:,:) ! normal vector of each triangular patch
		integer,allocatable:: node_of_patch(:,:) ! vertices of each triangular patch
		type(edge_node),allocatable:: edge_of_node(:) ! edges of each vertice
		CHARACTER (LEN=1000) DATA_DIR
		integer::CMmode=0 !  1: solve the characteristic mode, 0: solve the eigen mode
		integer::SI=0 ! 0: regular mode 1: shift-invert mode
		complex(kind=8)::shift=0d0 ! the shift value in shift-invert Arnoldi iterations
		integer::nev=1 ! nubmer of requested eigen values
		character(len=2) which ! which portion of eigen spectrum
		real(kind=8) tol_eig ! tolerance in arpack
		real(kind=8):: normalize_factor=1d0 ! normalization factor (default to be the acceleration voltage)
	end type quant_EMSURF

contains

subroutine delete_quant_EMSURF(quant)
implicit none
type(quant_EMSURF):: quant
if(allocated(quant%xyz))deallocate(quant%xyz)
if(allocated(quant%info_unk))deallocate(quant%info_unk)
if(allocated(quant%ng1))deallocate(quant%ng1)
if(allocated(quant%ng2))deallocate(quant%ng2)
if(allocated(quant%ng3))deallocate(quant%ng3)
if(allocated(quant%gauss_w))deallocate(quant%gauss_w)
if(allocated(quant%normal_of_patch))deallocate(quant%normal_of_patch)
if(allocated(quant%node_of_patch))deallocate(quant%node_of_patch)

end subroutine delete_quant_EMSURF


!**** user-defined subroutine to sample Z_mn
subroutine Zelem_EMSURF_T(m,n,value,quant)

    use BPACK_DEFS
    implicit none

    integer, INTENT(IN):: m,n
    integer flag,edge_m,edge_n,nodetemp_n,patch
    complex(kind=8) value_e,value_m,value
    integer i,j,ii,jj,iii,jjj
    real(kind=8) ln,lm,am(3),an(3),nr_m(3)
    real(kind=8) nr_n(3)
    complex(kind=8) ctemp,ctemp1,ctemp2,aa(3),bb(1),dg(3),dg1(3),dg2(3)
    complex(kind=8) imp,imp1,imp2,imp3
    real(kind=8) temp
    real(kind=8) distance
    real(kind=8) ianl,ianl1,ianl2
    real(kind=8) area

	type(quant_EMSURF) :: quant


    real(kind=8),allocatable::xm(:),ym(:),zm(:),wm(:),xn(:),yn(:),zn(:),wn(:)



		! convert to new indices because quant%info_unk has been reordered
		edge_m = m
		edge_n = n

		allocate (xm(quant%integral_points), ym(quant%integral_points), zm(quant%integral_points), wm(quant%integral_points))
		allocate (xn(quant%integral_points), yn(quant%integral_points), zn(quant%integral_points), wn(quant%integral_points))


		lm=sqrt((quant%xyz(1,quant%info_unk(1,edge_m))-quant%xyz(1,quant%info_unk(2,edge_m)))**2+(quant%xyz(2,quant%info_unk(1,edge_m))-quant%xyz(2,quant%info_unk(2,edge_m)))**2+(quant%xyz(3,quant%info_unk(1,edge_m))-quant%xyz(3,quant%info_unk(2,edge_m)))**2)
		ln=sqrt((quant%xyz(1,quant%info_unk(1,edge_n))-quant%xyz(1,quant%info_unk(2,edge_n)))**2+(quant%xyz(2,quant%info_unk(1,edge_n))-quant%xyz(2,quant%info_unk(2,edge_n)))**2+(quant%xyz(3,quant%info_unk(1,edge_n))-quant%xyz(3,quant%info_unk(2,edge_n)))**2)

		ctemp1=(0.,0.)
		ctemp2=(0.,0.)
		value_m=(0.,0.)

			do ii=3,4
				call gau_grobal(edge_m,ii,xm,ym,zm,wm,quant)
				nr_m(1:3)=quant%normal_of_patch(1:3,quant%info_unk(ii,edge_m))
				do i=1,quant%integral_points
					am(1)=xm(i)-quant%xyz(1,quant%info_unk(ii+2,edge_m))
					am(2)=ym(i)-quant%xyz(2,quant%info_unk(ii+2,edge_m))
					am(3)=zm(i)-quant%xyz(3,quant%info_unk(ii+2,edge_m))
					bb(1)=(0.,0.)
					aa(1:3)=(0.,0.)
					do jj=3,4
						call gau_grobal(edge_n,jj,xn,yn,zn,wn,quant)
						nr_n(1:3)=quant%normal_of_patch(1:3,quant%info_unk(jj,edge_n))
						if (quant%info_unk(ii,edge_m)==quant%info_unk(jj,edge_n)) then
							! area=triangle_area(quant%info_unk(ii,edge_m),quant)
							imp=(0.,0.);imp1=(0.,0.);imp2=(0.,0.);imp3=(0.,0.)
							do j=1,quant%integral_points
								distance=sqrt((xm(i)-xn(j))**2+(ym(i)-yn(j))**2+(zm(i)-zn(j))**2)
								if(distance==0)then
									imp=imp+wn(j)*(-junit*quant%wavenum)
									imp1=imp1+quant%ng1(j)*wn(j)*(-junit*quant%wavenum)
									imp2=imp2+quant%ng2(j)*wn(j)*(-junit*quant%wavenum)
									ianl=ianalytic(edge_n,jj,xn(j),yn(j),zn(j),quant)
									ianl1=ianalytic2(edge_n,jj,xn(j),yn(j),zn(j),1,quant)
									ianl2=ianalytic2(edge_n,jj,xn(j),yn(j),zn(j),2,quant)
									imp=imp+ianl  !debug
									imp1=imp1+ianl1
									imp2=imp2+ianl2
								else
									ctemp=exp(-junit*quant%wavenum*distance)/distance
									imp=imp+wn(j)*ctemp
									imp1=imp1+quant%ng1(j)*wn(j)*ctemp
									imp2=imp2+quant%ng2(j)*wn(j)*ctemp
								endif
							enddo
							imp3=imp-imp1-imp2
							nodetemp_n=quant%info_unk(jj+2,edge_n)
							patch=quant%info_unk(jj,edge_n)
							do jjj=1,3
								aa(jjj)=aa(jjj)+(-1)**(jj+1)*quant%wavenum**2*(quant%xyz(jjj,quant%node_of_patch(1,patch))*imp1+quant%xyz(jjj,quant%node_of_patch(2,patch))*imp2+quant%xyz(jjj,quant%node_of_patch(3,patch))*imp3-quant%xyz(jjj,nodetemp_n)*imp)
							enddo
							bb(1)=bb(1)+(-1)**(jj+1)*imp
						else
							imp=(0.,0.);imp1=(0.,0.);imp2=(0.,0.);imp3=(0.,0.)
							do j=1,quant%integral_points
								distance=sqrt((xm(i)-xn(j))**2+(ym(i)-yn(j))**2+(zm(i)-zn(j))**2)
								ctemp=exp(-junit*quant%wavenum*distance)/distance
								imp=imp+wn(j)*ctemp
								imp1=imp1+quant%ng1(j)*wn(j)*ctemp
								imp2=imp2+quant%ng2(j)*wn(j)*ctemp
							enddo
							imp3=imp-imp1-imp2
							nodetemp_n=quant%info_unk(jj+2,edge_n)
							patch=quant%info_unk(jj,edge_n)
							do jjj=1,3
								aa(jjj)=aa(jjj)+(-1)**(jj+1)*quant%wavenum**2*(quant%xyz(jjj,quant%node_of_patch(1,patch))*imp1+quant%xyz(jjj,quant%node_of_patch(2,patch))*imp2+quant%xyz(jjj,quant%node_of_patch(3,patch))*imp3-quant%xyz(jjj,nodetemp_n)*imp)
							enddo
							bb(1)=bb(1)+(-1)**(jj+1)*imp
						endif
					enddo
					call cscalar(aa,am,ctemp)
					ctemp1=ctemp1+(-1)**(ii+1)*ctemp*wm(i)
					ctemp2=ctemp2+4.*(-1)**(ii+1)*bb(1)*wm(i)
				enddo
			enddo
			value_e=ln*lm*junit*(ctemp1-ctemp2)/2./quant%freq/eps0
			value=value_e/impedence0

		deallocate(xm,ym,zm,wm,xn,yn,zn,wn)


    return

end subroutine Zelem_EMSURF_T

subroutine Zelem_EMSURF_K(m,n,value,quant,sign)

    use BPACK_DEFS
    implicit none

    integer, INTENT(IN):: m,n,sign
    integer flag,edge_m,edge_n,nodetemp_n,patch
    complex(kind=8) value_e,value_m,value
    integer i,j,ii,jj,iii,jjj
    real(kind=8) ln,lm,am(3),an(3),nr_m(3),nxan(3)
    real(kind=8) nr_n(3)
    complex(kind=8) ctemp,ctemp1,ctemp2,aa(3),bb(1),dg(3),dg1(3),dg2(3),dg3(3)
    complex(kind=8) imp,imp1,imp2,imp3
    real(kind=8) temp
    real(kind=8) distance
    real(kind=8) ianl,ianl1,ianl2
    real(kind=8) area

	type(quant_EMSURF) :: quant


    real(kind=8),allocatable::xm(:),ym(:),zm(:),wm(:),xn(:),yn(:),zn(:),wn(:)



		! convert to new indices because quant%info_unk has been reordered
		edge_m = m
		edge_n = n

		allocate (xm(quant%integral_points), ym(quant%integral_points), zm(quant%integral_points), wm(quant%integral_points))
		allocate (xn(quant%integral_points), yn(quant%integral_points), zn(quant%integral_points), wn(quant%integral_points))


		lm=sqrt((quant%xyz(1,quant%info_unk(1,edge_m))-quant%xyz(1,quant%info_unk(2,edge_m)))**2+(quant%xyz(2,quant%info_unk(1,edge_m))-quant%xyz(2,quant%info_unk(2,edge_m)))**2+(quant%xyz(3,quant%info_unk(1,edge_m))-quant%xyz(3,quant%info_unk(2,edge_m)))**2)
		ln=sqrt((quant%xyz(1,quant%info_unk(1,edge_n))-quant%xyz(1,quant%info_unk(2,edge_n)))**2+(quant%xyz(2,quant%info_unk(1,edge_n))-quant%xyz(2,quant%info_unk(2,edge_n)))**2+(quant%xyz(3,quant%info_unk(1,edge_n))-quant%xyz(3,quant%info_unk(2,edge_n)))**2)

		ctemp1=(0.,0.)
		ctemp2=(0.,0.)
		value_m=(0.,0.)

		do ii=3,4
			call gau_grobal(edge_m,ii,xm,ym,zm,wm,quant)
			nr_m(1:3)=quant%normal_of_patch(1:3,quant%info_unk(ii,edge_m))
			do i=1,quant%integral_points
				am(1)=xm(i)-quant%xyz(1,quant%info_unk(ii+2,edge_m))
				am(2)=ym(i)-quant%xyz(2,quant%info_unk(ii+2,edge_m))
				am(3)=zm(i)-quant%xyz(3,quant%info_unk(ii+2,edge_m))
				bb(1)=(0.,0.)
				aa(1:3)=(0.,0.)
				do jj=3,4
					call gau_grobal(edge_n,jj,xn,yn,zn,wn,quant)
					nr_n(1:3)=quant%normal_of_patch(1:3,quant%info_unk(jj,edge_n))
					if (quant%info_unk(ii,edge_m)==quant%info_unk(jj,edge_n)) then
						area=triangle_area(quant%info_unk(ii,edge_m),quant)
						an(1)=xm(i)-quant%xyz(1,quant%info_unk(jj+2,edge_n))
						an(2)=ym(i)-quant%xyz(2,quant%info_unk(jj+2,edge_n))
						an(3)=zm(i)-quant%xyz(3,quant%info_unk(jj+2,edge_n))
						call curl(nr_n,an,nxan)
						call scalar(am,nxan,temp)
						value_m=value_m+sign*(-1)**(ii+1)*(-1)**(jj+1)*0.5*temp/(2.*area)*wm(i)  ! note that wm contains the factor of 2
					else
						do j=1,quant%integral_points
							distance=sqrt((xm(i)-xn(j))**2+(ym(i)-yn(j))**2+(zm(i)-zn(j))**2)
							an(1)=xn(j)-quant%xyz(1,quant%info_unk(jj+2,edge_n))
							an(2)=yn(j)-quant%xyz(2,quant%info_unk(jj+2,edge_n))
							an(3)=zn(j)-quant%xyz(3,quant%info_unk(jj+2,edge_n))
							dg(1)=(xm(i)-xn(j))*(1+junit*quant%wavenum*distance)*exp(-junit*quant%wavenum*distance)/(4*pi*distance**3)
							dg(2)=(ym(i)-yn(j))*(1+junit*quant%wavenum*distance)*exp(-junit*quant%wavenum*distance)/(4*pi*distance**3)
							dg(3)=(zm(i)-zn(j))*(1+junit*quant%wavenum*distance)*exp(-junit*quant%wavenum*distance)/(4*pi*distance**3)
							 call ccurl(an,dg,dg2)
							 call cscalar(dg2,am,ctemp)
							 value_m=value_m-(-1)**(ii+1)*(-1)**(jj+1)*ctemp*wm(i)*wn(j)
						enddo
					endif
				enddo
			enddo
		enddo
		value_m=value_m*lm*ln

		value=value_m

		deallocate(xm,ym,zm,wm,xn,yn,zn,wn)


    return

end subroutine Zelem_EMSURF_K



subroutine Zelem_EMSURF_K_Self(m,n,value,quant,sign)

    use BPACK_DEFS
    implicit none

    integer, INTENT(IN):: m,n,sign
    integer flag,edge_m,edge_n,nodetemp_n,patch
    complex(kind=8) value_e,value_m,value
    integer i,j,ii,jj,iii,jjj
    real(kind=8) ln,lm,am(3),an(3),nr_m(3),nxan(3)
    real(kind=8) nr_n(3)
    complex(kind=8) ctemp,ctemp1,ctemp2,aa(3),bb(1),dg(3),dg1(3),dg2(3),dg3(3)
    complex(kind=8) imp,imp1,imp2,imp3
    real(kind=8) temp
    real(kind=8) distance
    real(kind=8) ianl,ianl1,ianl2
    real(kind=8) area

	type(quant_EMSURF) :: quant


    real(kind=8),allocatable::xm(:),ym(:),zm(:),wm(:),xn(:),yn(:),zn(:),wn(:)



		! convert to new indices because quant%info_unk has been reordered
		edge_m = m
		edge_n = n

		allocate (xm(quant%integral_points), ym(quant%integral_points), zm(quant%integral_points), wm(quant%integral_points))
		allocate (xn(quant%integral_points), yn(quant%integral_points), zn(quant%integral_points), wn(quant%integral_points))


		lm=sqrt((quant%xyz(1,quant%info_unk(1,edge_m))-quant%xyz(1,quant%info_unk(2,edge_m)))**2+(quant%xyz(2,quant%info_unk(1,edge_m))-quant%xyz(2,quant%info_unk(2,edge_m)))**2+(quant%xyz(3,quant%info_unk(1,edge_m))-quant%xyz(3,quant%info_unk(2,edge_m)))**2)
		ln=sqrt((quant%xyz(1,quant%info_unk(1,edge_n))-quant%xyz(1,quant%info_unk(2,edge_n)))**2+(quant%xyz(2,quant%info_unk(1,edge_n))-quant%xyz(2,quant%info_unk(2,edge_n)))**2+(quant%xyz(3,quant%info_unk(1,edge_n))-quant%xyz(3,quant%info_unk(2,edge_n)))**2)

		ctemp1=(0.,0.)
		ctemp2=(0.,0.)
		value_m=(0.,0.)

		do ii=3,4
			call gau_grobal(edge_m,ii,xm,ym,zm,wm,quant)
			nr_m(1:3)=quant%normal_of_patch(1:3,quant%info_unk(ii,edge_m))
			do i=1,quant%integral_points
				am(1)=xm(i)-quant%xyz(1,quant%info_unk(ii+2,edge_m))
				am(2)=ym(i)-quant%xyz(2,quant%info_unk(ii+2,edge_m))
				am(3)=zm(i)-quant%xyz(3,quant%info_unk(ii+2,edge_m))
				bb(1)=(0.,0.)
				aa(1:3)=(0.,0.)
				do jj=3,4
					call gau_grobal(edge_n,jj,xn,yn,zn,wn,quant)
					nr_n(1:3)=quant%normal_of_patch(1:3,quant%info_unk(jj,edge_n))
					if (quant%info_unk(ii,edge_m)==quant%info_unk(jj,edge_n)) then
						area=triangle_area(quant%info_unk(ii,edge_m),quant)
						an(1)=xm(i)-quant%xyz(1,quant%info_unk(jj+2,edge_n))
						an(2)=ym(i)-quant%xyz(2,quant%info_unk(jj+2,edge_n))
						an(3)=zm(i)-quant%xyz(3,quant%info_unk(jj+2,edge_n))
						! call curl(nr_n,an,nxan)
						! call scalar(am,nxan,temp)

						call scalar(am,an,temp)
						value_m=value_m+sign*(-1)**(ii+1)*(-1)**(jj+1)*0.5*temp/(2.*area)*wm(i)
					else

					endif
				enddo
			enddo
		enddo
		value_m=value_m*lm*ln

		value=value_m

		deallocate(xm,ym,zm,wm,xn,yn,zn,wn)


    return

end subroutine Zelem_EMSURF_K_Self


subroutine Port_nxe_dot_rwg(m,pp,mm,nn,TETM,rr,value,quant)

    use BPACK_DEFS
    implicit none

    integer, INTENT(IN):: m,pp,mm,nn,TETM,rr
    integer flag,edge_m,patch
    complex(kind=8) value_m,value
    integer i,ii
    real(kind=8) lm,am(3),nr_m(3),nxe(3)
    real(kind=8) temp
	type(quant_EMSURF) :: quant


    real(kind=8),allocatable::xm(:),ym(:),zm(:),wm(:),xn(:),yn(:),zn(:),wn(:)

		edge_m = m
		allocate (xm(quant%integral_points), ym(quant%integral_points), zm(quant%integral_points), wm(quant%integral_points))

		lm=sqrt((quant%xyz(1,quant%info_unk(1,edge_m))-quant%xyz(1,quant%info_unk(2,edge_m)))**2+(quant%xyz(2,quant%info_unk(1,edge_m))-quant%xyz(2,quant%info_unk(2,edge_m)))**2+(quant%xyz(3,quant%info_unk(1,edge_m))-quant%xyz(3,quant%info_unk(2,edge_m)))**2)

		value_m=(0.,0.)

		do ii=3,4
			call gau_grobal(edge_m,ii,xm,ym,zm,wm,quant)
			nr_m(1:3)=quant%normal_of_patch(1:3,quant%info_unk(ii,edge_m))
			do i=1,quant%integral_points
				am(1)=xm(i)-quant%xyz(1,quant%info_unk(ii+2,edge_m))
				am(2)=ym(i)-quant%xyz(2,quant%info_unk(ii+2,edge_m))
				am(3)=zm(i)-quant%xyz(3,quant%info_unk(ii+2,edge_m))

				call Port_nxe(xm(i),ym(i),zm(i),nxe,quant,pp,mm,nn,TETM,rr)
				call scalar(am,nxe,temp)
				value_m=value_m+(-1)**(ii+1)*temp/2.*wm(i)
			enddo
		enddo
		value_m=value_m*lm

		value=value_m

		deallocate(xm,ym,zm,wm)


    return

end subroutine Port_nxe_dot_rwg

subroutine Port_e_dot_rwg(m,pp,mm,nn,TETM,rr,value,quant)

    use BPACK_DEFS
    implicit none

    integer, INTENT(IN):: m,pp,mm,nn,TETM,rr
    integer flag,edge_m,patch
    complex(kind=8) value_m,value
    integer i,ii
    real(kind=8) lm,am(3),nr_m(3),e(3)
    real(kind=8) temp
	type(quant_EMSURF) :: quant


    real(kind=8),allocatable::xm(:),ym(:),zm(:),wm(:),xn(:),yn(:),zn(:),wn(:)

		edge_m = m
		allocate (xm(quant%integral_points), ym(quant%integral_points), zm(quant%integral_points), wm(quant%integral_points))

		lm=sqrt((quant%xyz(1,quant%info_unk(1,edge_m))-quant%xyz(1,quant%info_unk(2,edge_m)))**2+(quant%xyz(2,quant%info_unk(1,edge_m))-quant%xyz(2,quant%info_unk(2,edge_m)))**2+(quant%xyz(3,quant%info_unk(1,edge_m))-quant%xyz(3,quant%info_unk(2,edge_m)))**2)

		value_m=(0.,0.)

		do ii=3,4
			call gau_grobal(edge_m,ii,xm,ym,zm,wm,quant)
			nr_m(1:3)=quant%normal_of_patch(1:3,quant%info_unk(ii,edge_m))
			do i=1,quant%integral_points
				am(1)=xm(i)-quant%xyz(1,quant%info_unk(ii+2,edge_m))
				am(2)=ym(i)-quant%xyz(2,quant%info_unk(ii+2,edge_m))
				am(3)=zm(i)-quant%xyz(3,quant%info_unk(ii+2,edge_m))

				call Port_e(xm(i),ym(i),zm(i),e,quant,pp,mm,nn,TETM,rr)
				call scalar(am,e,temp)
				value_m=value_m+(-1)**(ii+1)*temp/2.*wm(i)
			enddo
		enddo
		value_m=value_m*lm

		value=value_m

		deallocate(xm,ym,zm,wm)


    return

end subroutine Port_e_dot_rwg


subroutine Port_nxe(xm,ym,zm,nxe,quant,pp,mm,nn,TETM,rr)
	real(kind=8) xm,ym,zm,e(3),nxe(3),r,theta,phi,Erho,Ephi,Ex,Ey,kc,x,y,z,nhat(3),rhat(3),phihat(3),rho(3),a,b
	type(quant_EMSURF) :: quant
	integer:: pp,mm,nn,TETM,rr
	real(kind=8):: J(nn+2),Jp(nn+2)

	if(quant%ports(pp)%type==0)then

		Erho=0
		Ephi=0

		call Cart2Sph_Loc(xm, ym, zm, quant%ports(pp)%origin, quant%ports(pp)%x, quant%ports(pp)%y, quant%ports(pp)%z, r, theta, phi)

		nhat = quant%ports(pp)%z
		rhat = 0
		phihat = 0
		if(r>0)then
			rhat = (/xm-quant%ports(pp)%origin(1),ym-quant%ports(pp)%origin(2),zm-quant%ports(pp)%origin(3)/)/r
			call curl(nhat,rhat,phihat)
		endif

		if(TETM==1)then !TE_nm
			kc = r_TE_nm(nn+1,mm)/quant%ports(pp)%R
			x = kc*r
			if(nn==0)then
				call bessjyV(1,x,J,0)
				Jp(1) = -J(2)
			else
				call bessjyV(nn,x,J,0)
				Jp(nn+1) = J(nn)-(nn)/x*J(nn+1)
			endif
			if(rr==1)then
				Erho = nn/x*sin(nn*phi)*J(nn+1)
				Ephi = cos(nn*phi)*Jp(nn+1)
			else
				Erho = -nn/x*cos(nn*phi)*J(nn+1)
				Ephi = sin(nn*phi)*Jp(nn+1)
			endif
			e = (Erho*rhat+Ephi*phihat)*quant%ports(pp)%A_TE_nm(nn+1,mm)/quant%ports(pp)%R
		elseif(TETM==2)then ! TM_nm
			kc = r_TM_nm(nn+1,mm)/quant%ports(pp)%R
			x = kc*r
			if(nn==0)then
				call bessjyV(1,x,J,0)
				Jp(1) = -J(2)
			else
				call bessjyV(nn,x,J,0)
				Jp(nn+1) = J(nn)-(nn)/x*J(nn+1)
			endif
			if(rr==1)then
				Erho = cos(nn*phi)*Jp(nn+1)
				Ephi = -nn/x*sin(nn*phi)*J(nn+1)
			else
				Erho = sin(nn*phi)*Jp(nn+1)
				Ephi = nn/x*cos(nn*phi)*J(nn+1)
			endif
			e = (Erho*rhat+Ephi*phihat)*quant%ports(pp)%A_TM_nm(nn+1,mm)/quant%ports(pp)%R
		endif
		call curl(nhat,e,nxe)
	elseif(quant%ports(pp)%type==1)then
		Ex=0
		Ey=0
		rho(1)=xm
		rho(2)=ym
		rho(3)=zm
		rho = rho - quant%ports(pp)%origin
		x = dot_product(rho,quant%ports(pp)%x)
		y = dot_product(rho,quant%ports(pp)%y)
		z = dot_product(rho,quant%ports(pp)%z)

		a = quant%ports(pp)%a
		b = quant%ports(pp)%b

		if(TETM==1)then !TE_nm
			Ex = mm/b*cos(nn*pi*x/a)*sin(mm*pi*y/b)
			Ey = -nn/a*sin(nn*pi*x/a)*cos(mm*pi*y/b)
			e = (Ex*quant%ports(pp)%x+Ey*quant%ports(pp)%y)*quant%ports(pp)%A_TE_nm(nn+1,mm+1)
		elseif(TETM==2)then ! TM_nm
			Ex = nn/a*cos(nn*pi*x/a)*sin(mm*pi*y/b)
			Ey = mm/b*sin(nn*pi*x/a)*cos(mm*pi*y/b)
			e = (Ex*quant%ports(pp)%x+Ey*quant%ports(pp)%y)*quant%ports(pp)%A_TM_nm(nn+1,mm+1)
		endif
		call curl(quant%ports(pp)%z,e,nxe)
	else
		write(*,*)'unrecognized port type',quant%ports(pp)%type
		stop
	endif

end subroutine Port_nxe


subroutine Port_e(xm,ym,zm,e,quant,pp,mm,nn,TETM,rr)
	real(kind=8) xm,ym,zm,e(3),nxe(3),r,theta,phi,Erho,Ephi,Ex,Ey,kc,x,y,z,nhat(3),rhat(3),phihat(3),rho(3),a,b
	type(quant_EMSURF) :: quant
	integer:: pp,mm,nn,TETM,rr
	real(kind=8):: J(nn+2),Jp(nn+2)

	if(quant%ports(pp)%type==0)then

		Erho=0
		Ephi=0

		call Cart2Sph_Loc(xm, ym, zm, quant%ports(pp)%origin, quant%ports(pp)%x, quant%ports(pp)%y, quant%ports(pp)%z, r, theta, phi)

		nhat = quant%ports(pp)%z
		rhat = 0
		phihat = 0
		if(r>0)then
			rhat = (/xm-quant%ports(pp)%origin(1),ym-quant%ports(pp)%origin(2),zm-quant%ports(pp)%origin(3)/)/r
			call curl(nhat,rhat,phihat)
		endif

		if(TETM==1)then !TE_nm
			kc = r_TE_nm(nn+1,mm)/quant%ports(pp)%R
			x = kc*r
			if(nn==0)then
				call bessjyV(1,x,J,0)
				Jp(1) = -J(2)
			else
				call bessjyV(nn,x,J,0)
				Jp(nn+1) = J(nn)-(nn)/x*J(nn+1)
			endif
			if(rr==1)then
				Erho = nn/x*sin(nn*phi)*J(nn+1)
				Ephi = cos(nn*phi)*Jp(nn+1)
			else
				Erho = -nn/x*cos(nn*phi)*J(nn+1)
				Ephi = sin(nn*phi)*Jp(nn+1)
			endif
			e = (Erho*rhat+Ephi*phihat)*quant%ports(pp)%A_TE_nm(nn+1,mm)/quant%ports(pp)%R
		elseif(TETM==2)then ! TM_nm
			kc = r_TM_nm(nn+1,mm)/quant%ports(pp)%R
			x = kc*r
			if(nn==0)then
				call bessjyV(1,x,J,0)
				Jp(1) = -J(2)
			else
				call bessjyV(nn,x,J,0)
				Jp(nn+1) = J(nn)-(nn)/x*J(nn+1)
			endif
			if(rr==1)then
				Erho = cos(nn*phi)*Jp(nn+1)
				Ephi = -nn/x*sin(nn*phi)*J(nn+1)
			else
				Erho = sin(nn*phi)*Jp(nn+1)
				Ephi = nn/x*cos(nn*phi)*J(nn+1)
			endif
			e = (Erho*rhat+Ephi*phihat)*quant%ports(pp)%A_TM_nm(nn+1,mm)/quant%ports(pp)%R
		endif
		! call curl(nhat,e,nxe)
	elseif(quant%ports(pp)%type==1)then
		Ex=0
		Ey=0
		rho(1)=xm
		rho(2)=ym
		rho(3)=zm
		rho = rho - quant%ports(pp)%origin
		x = dot_product(rho,quant%ports(pp)%x)
		y = dot_product(rho,quant%ports(pp)%y)
		z = dot_product(rho,quant%ports(pp)%z)

		a = quant%ports(pp)%a
		b = quant%ports(pp)%b

		if(TETM==1)then !TE_nm
			Ex = mm/b*cos(nn*pi*x/a)*sin(mm*pi*y/b)
			Ey = -nn/a*sin(nn*pi*x/a)*cos(mm*pi*y/b)
			e = (Ex*quant%ports(pp)%x+Ey*quant%ports(pp)%y)*quant%ports(pp)%A_TE_nm(nn+1,mm+1)
		elseif(TETM==2)then ! TM_nm
			Ex = nn/a*cos(nn*pi*x/a)*sin(mm*pi*y/b)
			Ey = mm/b*sin(nn*pi*x/a)*cos(mm*pi*y/b)
			e = (Ex*quant%ports(pp)%x+Ey*quant%ports(pp)%y)*quant%ports(pp)%A_TM_nm(nn+1,mm+1)
		endif
		! call curl(quant%ports(pp)%z,e,nxe)
	else
		write(*,*)'unrecognized port type',quant%ports(pp)%type
		stop
	endif

end subroutine Port_e



!**** user-defined subroutine to sample Z_mn
subroutine Zelem_EMSURF(m,n,value,quant)

    use BPACK_DEFS
    implicit none

    integer, INTENT(IN):: m,n
    integer flag,edge_m,edge_n,nodetemp_n,patch,cntm,cntn,ppm,ppn
    complex(kind=8) value_e,value_m,value
    integer i,j,ii,jj,iii,jjj,nn,mm,mode,rr
    real(kind=8) ln,lm,am(3),an(3),nr_m(3)
    real(kind=8) nr_n(3)
    complex(kind=8) ctemp,ctemp1,ctemp2,aa(3),bb(1),dg(3),dg1(3),dg2(3)
    complex(kind=8) imp,imp1,imp2,imp3
    real(kind=8) temp
    real(kind=8) distance
    real(kind=8) ianl,ianl1,ianl2
	real(kind=8) area
	integer sign,npolar,off

	class(*),pointer :: quant


    real(kind=8),allocatable::xm(:),ym(:),zm(:),wm(:),xn(:),yn(:),zn(:),wn(:)


	select TYPE(quant)
	type is (quant_EMSURF)

		! convert to new indices because quant%info_unk has been reordered


		if(m<=quant%Nunk-quant%Nunk_port .and. n<=quant%Nunk-quant%Nunk_port)then
			edge_m = m
			edge_n = n
			call Zelem_EMSURF_T(edge_m,edge_n,value,quant)
		elseif(m>quant%Nunk-quant%Nunk_port .and. n>quant%Nunk-quant%Nunk_port)then
			edge_m = m-quant%Nunk_port
			edge_n = n-quant%Nunk_port
			call Zelem_EMSURF_T(edge_m,edge_n,value,quant)

			! find which ports this edge is on
			cntm = quant%Nunk_int
			do ppm=1,quant%Nport
				if(edge_m<=cntm+quant%ports(ppm)%Nunk)exit
				cntm = cntm + quant%ports(ppm)%Nunk
			enddo
			cntn = quant%Nunk_int
			do ppn=1,quant%Nport
				if(edge_n<=cntn+quant%ports(ppn)%Nunk)exit
				cntn = cntn + quant%ports(ppn)%Nunk
			enddo
			if(ppm==ppn)then
				if(quant%ports(ppm)%type==0)then
					npolar=2
					off=0
				elseif(quant%ports(ppm)%type==1)then
					npolar=1
					off=1
				else
					write(*,*)'unrecognized port type',quant%ports(ppm)%type
					stop
				endif
				ctemp=0
				do rr=1,npolar
					do nn=0,quant%ports(ppm)%nmax
						do mm=1-off,quant%ports(ppm)%mmax
							if(mm>0 .or. nn>0)then
							ctemp1=quant%ports(ppm)%nxe_dot_rwg(edge_m-cntm,nn+1,mm+off,1,rr)
							ctemp2=quant%ports(ppm)%nxe_dot_rwg(edge_n-cntn,nn+1,mm+off,1,rr)
							ctemp = ctemp + impedence0/2/quant%ports(ppm)%impedance_TE_nm(nn+1,mm+off)*ctemp1*ctemp2
							endif
						enddo
					enddo
					do nn=0,quant%ports(ppm)%nmax
						do mm=1-off,quant%ports(ppm)%mmax
							if(mm>0 .or. nn>0)then
							ctemp1=quant%ports(ppm)%nxe_dot_rwg(edge_m-cntm,nn+1,mm+off,2,rr)
							ctemp2=quant%ports(ppm)%nxe_dot_rwg(edge_n-cntn,nn+1,mm+off,2,rr)
							ctemp = ctemp + impedence0/2/quant%ports(ppm)%impedance_TM_nm(nn+1,mm+off)*ctemp1*ctemp2
							endif
						enddo
					enddo
				enddo
				value = value +ctemp
			endif


		elseif(m>quant%Nunk-quant%Nunk_port .and. n<=quant%Nunk-quant%Nunk_port)then
			edge_m = m-quant%Nunk_port
			edge_n = n
			sign=0
			call Zelem_EMSURF_K(edge_m,edge_n,value,quant,sign)
			value=-value

		elseif(m<=quant%Nunk-quant%Nunk_port .and. n>quant%Nunk-quant%Nunk_port)then
			edge_m = m
			edge_n = n-quant%Nunk_port
			if(m<=quant%Nunk_int)sign=-1
			if(m>quant%Nunk_int)sign=1
			! sign=-1
			call Zelem_EMSURF_K(edge_m,edge_n,value,quant,sign)
		endif


		value = value*impedence0 ! the solution vector will be J and M/impedence0, this makes it easier to compare with ie3deigen

	class default
		write(*,*)"unexpected type"
		stop
	end select

    return

end subroutine Zelem_EMSURF




!**** user-defined subroutine to sample Z_mn for calculating the normal E fields
subroutine Zelem_EMSURF_Post(m,n,value,quant)

    use BPACK_DEFS
    implicit none

    integer, INTENT(IN):: m,n
    integer flag,edge_m,edge_n,nodetemp_n,patch,cntm,cntn,ppm,ppn
    complex(kind=8) value_e,value_m,value,value0
    integer i,j,ii,jj,iii,jjj,nn,mm,mode,rr
    real(kind=8) ln,lm,am(3),an(3),nr_m(3),point(3),l1,l2,l3
    real(kind=8) nr_n(3)
    complex(kind=8) ctemp,ctemp1,ctemp2,aa(3),bb(1),dg(3),dg1(3),dg2(3),field(3)
    complex(kind=8) imp,imp1,imp2,imp3
    real(kind=8) temp
    real(kind=8) distance
    real(kind=8) ianl,ianl1,ianl2
	real(kind=8) area
	integer sign,npolar,off

	class(*),pointer :: quant


    real(kind=8),allocatable::xm(:),ym(:),zm(:),wm(:),xn(:),yn(:),zn(:),wn(:)


	select TYPE(quant)
	type is (quant_EMSURF)

		if(m>quant%Nunk-quant%Nunk_port)then
			value=0d0
		else
			! Calulate the average normal E field between the center of the two patches (with an offset towards the inside)
			value=0d0
			do ii = 3,4
			edge_m = m
			patch = quant%info_unk(ii,edge_m)
			nr_m(1:3)=quant%normal_of_patch(1:3,patch)
			do i=1,3
				point(i)=1./3.*(quant%xyz(i,quant%node_of_patch(1,patch))+quant%xyz(i,quant%node_of_patch(2,patch))+quant%xyz(i,quant%node_of_patch(3,patch)))
			enddo

			l1=sqrt(sum((quant%xyz(:,quant%node_of_patch(1,patch))-quant%xyz(:,quant%node_of_patch(2,patch)))**2d0))
			l2=sqrt(sum((quant%xyz(:,quant%node_of_patch(1,patch))-quant%xyz(:,quant%node_of_patch(3,patch)))**2d0))
			l3=sqrt(sum((quant%xyz(:,quant%node_of_patch(3,patch))-quant%xyz(:,quant%node_of_patch(2,patch)))**2d0))
			ln= max(l1,max(l2,l3))

			! if(abs(point(3)-0.14454)<1e-4 .or. abs(point(3)+0.14454)<1e-4)write(*,*)point(3),nr_m(3)
			! point = point - nr_m*(ln)  ! use largest edge length as an offset inwards, to avoid singularity
			point = point- nr_m*(ln)/20
			call Field_EMSURF(point,field,n,quant)
			value0 = dot_product(field,nr_m)
			value0 = value0*impedence0
			value = value + value0
			enddo
			! value = value/2
		endif

	class default
		write(*,*)"unexpected type"
		stop
	end select

    return

end subroutine Zelem_EMSURF_Post



! !**** user-defined subroutine to sample Z_mn
! subroutine Zelem_EMSURF(m,n,value,quant)

!     use BPACK_DEFS
!     implicit none

!     integer, INTENT(IN):: m,n
!     integer flag,edge_m,edge_n,nodetemp_n,patch,cntm,cntn,ppm,ppn
!     complex(kind=8) value_e,value_m,value
!     integer i,j,ii,jj,iii,jjj,nn,mm,mode,rr
!     real(kind=8) ln,lm,am(3),an(3),nr_m(3)
!     real(kind=8) nr_n(3)
!     complex(kind=8) ctemp,ctemp1,ctemp2,aa(3),bb(1),dg(3),dg1(3),dg2(3)
!     complex(kind=8) imp,imp1,imp2,imp3
!     real(kind=8) temp
!     real(kind=8) distance
!     real(kind=8) ianl,ianl1,ianl2
! 	real(kind=8) area
! 	integer sign,npolar,off

! 	class(*),pointer :: quant


!     real(kind=8),allocatable::xm(:),ym(:),zm(:),wm(:),xn(:),yn(:),zn(:),wn(:)


! 	select TYPE(quant)
! 	type is (quant_EMSURF)

! 		! convert to new indices because quant%info_unk has been reordered


! 		if(m<=quant%Nunk-quant%Nunk_port .and. n<=quant%Nunk-quant%Nunk_port)then
! 			edge_m = m
! 			edge_n = n
! 			call Zelem_EMSURF_T(edge_m,edge_n,value,quant)
! 		elseif(m>quant%Nunk-quant%Nunk_port .and. n>quant%Nunk-quant%Nunk_port)then
! 			edge_m = m-quant%Nunk_port
! 			edge_n = n-quant%Nunk_port
! 			! call Zelem_EMSURF_T(edge_m,edge_n,value,quant)

! 			! find which ports this edge is on
! 			cntm = quant%Nunk_int
! 			do ppm=1,quant%Nport
! 				if(edge_m<=cntm+quant%ports(ppm)%Nunk)exit
! 				cntm = cntm + quant%ports(ppm)%Nunk
! 			enddo
! 			cntn = quant%Nunk_int
! 			do ppn=1,quant%Nport
! 				if(edge_n<=cntn+quant%ports(ppn)%Nunk)exit
! 				cntn = cntn + quant%ports(ppn)%Nunk
! 			enddo
! 			if(ppm==ppn)then
! 				if(quant%ports(ppm)%type==0)then
! 					npolar=2
! 					off=0
! 				elseif(quant%ports(ppm)%type==1)then
! 					npolar=1
! 					off=1
! 				else
! 					write(*,*)'unrecognized port type',quant%ports(ppm)%type
! 					stop
! 				endif
! 				ctemp=0
! 				do rr=1,npolar
! 					do nn=0,quant%ports(ppm)%nmax
! 						do mm=1-off,quant%ports(ppm)%mmax
! 							if(mm>0 .or. nn>0)then
! 							ctemp1=quant%ports(ppm)%e_dot_rwg(edge_m-cntm,nn+1,mm+off,1,rr)
! 							ctemp2=quant%ports(ppm)%nxe_dot_rwg(edge_n-cntn,nn+1,mm+off,1,rr)
! 							ctemp = ctemp + impedence0/2/quant%ports(ppm)%impedance_TE_nm(nn+1,mm+off)*ctemp1*ctemp2
! 							endif
! 						enddo
! 					enddo
! 					do nn=0,quant%ports(ppm)%nmax
! 						do mm=1-off,quant%ports(ppm)%mmax
! 							if(mm>0 .or. nn>0)then
! 							ctemp1=quant%ports(ppm)%e_dot_rwg(edge_m-cntm,nn+1,mm+off,2,rr)
! 							ctemp2=quant%ports(ppm)%nxe_dot_rwg(edge_n-cntn,nn+1,mm+off,2,rr)
! 							ctemp = ctemp + impedence0/2/quant%ports(ppm)%impedance_TM_nm(nn+1,mm+off)*ctemp1*ctemp2
! 							endif
! 						enddo
! 					enddo
! 				enddo
! 				value = value +ctemp
! 			endif

! 		elseif(m>quant%Nunk-quant%Nunk_port .and. n<=quant%Nunk-quant%Nunk_port)then
! 			edge_m = m-quant%Nunk_port
! 			edge_n = n
! 			sign = 1
! 			call Zelem_EMSURF_K_Self(edge_m,edge_n,value,quant,sign)
! 			value=-value

! 		elseif(m<=quant%Nunk-quant%Nunk_port .and. n>quant%Nunk-quant%Nunk_port)then
! 			edge_m = m
! 			edge_n = n-quant%Nunk_port
! 			if(m<=quant%Nunk_int)sign=-1
! 			if(m>quant%Nunk_int)sign=1
! 			! sign=-1
! 			call Zelem_EMSURF_K(edge_m,edge_n,value,quant,sign)
! 		endif


! 	class default
! 		write(*,*)"unexpected type"
! 		stop
! 	end select

!     return

! end subroutine Zelem_EMSURF



	!**** user-defined subroutine to sample real(Z_mn)
	subroutine Zelem_EMSURF_Real(m,n,value_e,quant)
		use BPACK_DEFS
		implicit none
		integer, INTENT(IN):: m,n
		complex(kind=8) value_e
		class(*),pointer :: quant

		call Zelem_EMSURF(m,n,value_e,quant)
		value_e=dble(value_e)
	end subroutine  Zelem_EMSURF_Real


	!**** user-defined subroutine to sample Z_mn-sigma*Delta_mn or Z_mn-sigma*real(Z_mn)
	subroutine Zelem_EMSURF_Shifted(m,n,value_e,quant)

		use BPACK_DEFS
		implicit none

		integer edge_m, edge_n, i, j, flag
		integer, INTENT(IN):: m,n
		real(kind=8) r_mn, rtemp1, rtemp2
		complex(kind=8) value_e

		class(*),pointer :: quant

		call Zelem_EMSURF(m,n,value_e,quant)
		select TYPE(quant)
			type is (quant_EMSURF)
				if(quant%CMmode==1)then
					value_e = value_e - quant%shift*dble(value_e)
				else
					if(m==n)then
						value_e = value_e - quant%shift
					endif
				endif
			class default
				write(*,*)"unexpected type"
				stop
			end select

		return

	end subroutine Zelem_EMSURF_Shifted



!***********************************
!	seven points gauss integration
!	provide w(n) and x,y,z(n) (n=7)
!***********************************
  subroutine gau_grobal(nn,j,x,y,z,w,quant)

  use BPACK_DEFS
  implicit none
  type(quant_EMSURF)::quant
  integer nn ,flag
  real(kind=8) x(:),y(:),z(:),w(:)
  integer i,j,ii
!	wa=area*9./40
!	wb=area*(155.-sqrt(15.))/1200.
!	wc=area*(155.+sqrt(15.))/1200.

   ! if(flag==1) then

!		ii=quant%info_unkfine(j,nn)
!  do i=1,integral_points
!	x(i)=quant%ng1(i)*quant%xyz(1,node_of_patchfine(1,ii))+quant%ng2(i)*quant%xyz(1,node_of_patchfine(2,ii))+&
!     quant%ng3(i)*quant%xyz(1,node_of_patchfine(3,ii))
!	y(i)=quant%ng1(i)*quant%xyz(2,node_of_patchfine(1,ii))+quant%ng2(i)*quant%xyz(2,node_of_patchfine(2,ii))+&
!     quant%ng3(i)*quant%xyz(2,node_of_patchfine(3,ii))
!	z(i)=quant%ng1(i)*quant%xyz(3,node_of_patchfine(1,ii))+quant%ng2(i)*quant%xyz(3,node_of_patchfine(2,ii))+&
!     quant%ng3(i)*quant%xyz(3,node_of_patchfine(3,ii))
!  enddo

!  elseif(flag==2) then

  ii=quant%info_unk(j,nn)
	 do i=1,quant%integral_points
	x(i)=quant%ng1(i)*quant%xyz(1,quant%node_of_patch(1,ii))+quant%ng2(i)*quant%xyz(1,quant%node_of_patch(2,ii))+&
     quant%ng3(i)*quant%xyz(1,quant%node_of_patch(3,ii))
	y(i)=quant%ng1(i)*quant%xyz(2,quant%node_of_patch(1,ii))+quant%ng2(i)*quant%xyz(2,quant%node_of_patch(2,ii))+&
     quant%ng3(i)*quant%xyz(2,quant%node_of_patch(3,ii))
	z(i)=quant%ng1(i)*quant%xyz(3,quant%node_of_patch(1,ii))+quant%ng2(i)*quant%xyz(3,quant%node_of_patch(2,ii))+&
     quant%ng3(i)*quant%xyz(3,quant%node_of_patch(3,ii))
  enddo

  w=quant%gauss_w

! 	w(1)=wa
! 	w(2)=wb
! 	w(3)=wb
! 	w(4)=wb
! 	w(5)=wc
! 	w(6)=wc
! 	w(7)=wc
  return
  end subroutine gau_grobal


  subroutine gauss_points(quant)

      use BPACK_DEFS
      implicit none

      real(kind=8) v1,v2,v3,v4,v5
      real(kind=8) wa,wb,wc
	  type(quant_EMSURF)::quant

    !	wa=area*9./40
!	wb=area*(155.-sqrt(15.))/1200.
!	wc=area*(155.+sqrt(15.))/1200.
    if (quant%integral_points==7) then

    wa=9./80
	wb=(155.-sqrt(15.))/2400.
	wc=(155.+sqrt(15.))/2400.

	v1=1d0/3d0
	v2=(6d0-sqrt(15d0))/21d0
	v3=(9d0+2*sqrt(15d0))/21d0
	v4=(6d0+sqrt(15d0))/21d0
	v5=(9d0-2*sqrt(15d0))/21d0

	quant%ng1(1)=1-v1-v1
	quant%ng1(2)=1-v2-v2
	quant%ng1(3)=1-v2-v3
	quant%ng1(4)=1-v3-v2
	quant%ng1(5)=1-v4-v4
	quant%ng1(6)=1-v4-v5
	quant%ng1(7)=1-v5-v4

	quant%ng2(1)=v1
	quant%ng2(2)=v2
	quant%ng2(3)=v2
	quant%ng2(4)=v3
	quant%ng2(5)=v4
	quant%ng2(6)=v4
	quant%ng2(7)=v5

	quant%ng3(1)=v1
	quant%ng3(2)=v2
	quant%ng3(3)=v3
	quant%ng3(4)=v2
	quant%ng3(5)=v4
	quant%ng3(6)=v5
	quant%ng3(7)=v4

	quant%gauss_w(1)=wa
	quant%gauss_w(2)=wb
	quant%gauss_w(3)=wb
	quant%gauss_w(4)=wb
	quant%gauss_w(5)=wc
	quant%gauss_w(6)=wc
	quant%gauss_w(7)=wc

	elseif (quant%integral_points==6) then

	v1=0.816847572980459d0
    v2=0.091576213509771d0
    v3=0.108103018168070d0
    v4=0.445948490915965d0
    wa=0.109951743655322d0/2.
    wb=0.223381589678011d0/2.

	quant%ng1(1)=v1
	quant%ng1(2)=v2
	quant%ng1(3)=v2
	quant%ng1(4)=v3
	quant%ng1(5)=v4
	quant%ng1(6)=v4
	!quant%ng1(7)=1-v5-v4

	quant%ng2(1)=v2
	quant%ng2(2)=v1
	quant%ng2(3)=v2
	quant%ng2(4)=v4
	quant%ng2(5)=v3
	quant%ng2(6)=v4
	!quant%ng2(7)=v5

	quant%ng3(1)=v2
	quant%ng3(2)=v2
	quant%ng3(3)=v1
	quant%ng3(4)=v4
	quant%ng3(5)=v4
	quant%ng3(6)=v3

	quant%gauss_w(1)=wa
	quant%gauss_w(2)=wa
	quant%gauss_w(3)=wa
	quant%gauss_w(4)=wb
	quant%gauss_w(5)=wb
	quant%gauss_w(6)=wb

	elseif (quant%integral_points==4) then

	v1=1./3.
	v2=0.2
	v3=0.6
	wa=-27./96.
	wb=25./96.

	quant%ng1(1)=v1
	quant%ng1(2)=v2
	quant%ng1(3)=v2
	quant%ng1(4)=v3

	quant%ng2(1)=v1
	quant%ng2(2)=v2
	quant%ng2(3)=v3
	quant%ng2(4)=v2

	quant%ng3(1)=v1
	quant%ng3(2)=v3
	quant%ng3(3)=v2
	quant%ng3(4)=v2

	quant%gauss_w(1)=wa
	quant%gauss_w(2)=wb
	quant%gauss_w(3)=wb
	quant%gauss_w(4)=wb

	endif

    return
  end subroutine gauss_points


!**********************************
!	mm:commom edge
!	jj:face number(3 or 4)
!**********************************
function ianalytic(mm,jj,xi,yi,zi,quant)

use     BPACK_DEFS
integer mm,jj,j,i
real(kind=8) xi,yi,zi
real(kind=8)    temp,ianalytic
integer ii,node1,node2,node3
real(kind=8)    u3,v3,u0,v0,w0,l(3)
real(kind=8)    u(3),w(3),v(3),a(3),b(3)
real(kind=8)    s(2,3),t(-1:1,3)
real(kind=8)    f2(3),beta(3)
real(kind=8)    r(-1:1,3)
real(kind=8)    area
type(quant_EMSURF)::quant

ii=quant%info_unk(jj,mm)
node1=quant%node_of_patch(1,ii)
node2=quant%node_of_patch(2,ii)
node3=quant%node_of_patch(3,ii)
do i=1,3
   a(i)=quant%xyz(i,node2)-quant%xyz(i,node1)
   b(i)=quant%xyz(i,node3)-quant%xyz(i,node1)
enddo
 call curl(a,b,w)
area=0.5*sqrt(w(1)**2+w(2)**2+w(3)**2)
do i=1,3
   w(i)=w(i)/2./area
enddo
l(1)=sqrt((quant%xyz(1,node3)-quant%xyz(1,node2))**2+(quant%xyz(2,node3)&
     	-quant%xyz(2,node2))**2+(quant%xyz(3,node3)-quant%xyz(3,node2))**2)
l(2)=sqrt((quant%xyz(1,node3)-quant%xyz(1,node1))**2+(quant%xyz(2,node3)&
     	-quant%xyz(2,node1))**2+(quant%xyz(3,node3)-quant%xyz(3,node1))**2)
l(3)=sqrt((quant%xyz(1,node1)-quant%xyz(1,node2))**2+(quant%xyz(2,node1)&
     	-quant%xyz(2,node2))**2+(quant%xyz(3,node1)-quant%xyz(3,node2))**2)
do i=1,3
   u(i)=a(i)/l(3)
enddo
 call curl(w,u,v)
 call scalar(u,b,u3)
v3=2.*area/l(3)

b(1)=xi-quant%xyz(1,node1)
b(2)=yi-quant%xyz(2,node1)
b(3)=zi-quant%xyz(3,node1)
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
function ianalytic2(mm,jj,xi,yi,zi,iii,quant)

use BPACK_DEFS
integer mm,jj,j,i
real(kind=8) ianalytic2
integer ii,node1,node2,node3,iii
real(kind=8) xi,yi,zi
real(kind=8) u3,v3,u0,v0,w0,l(3)
real(kind=8) u(3),w(3),v(3),a(3),b(3)
real(kind=8) m1(3),m2(3),m3(3)
real(kind=8) s1(3),s2(3),s3(3)
real(kind=8) s(2,3),t(-1:1,3)
real(kind=8) f2(3),beta(3),f3(3)
real(kind=8) r(-1:1,3)
real(kind=8) temp,temp1,temp2,temp3
real(kind=8) iua,iva
real(kind=8) n0(3)
real(kind=8)    area
type(quant_EMSURF)::quant

ii=quant%info_unk(jj,mm)
node1=quant%node_of_patch(1,ii)
node2=quant%node_of_patch(2,ii)
node3=quant%node_of_patch(3,ii)
do i=1,3
   a(i)=quant%xyz(i,node2)-quant%xyz(i,node1)
   b(i)=quant%xyz(i,node3)-quant%xyz(i,node1)
enddo
call curl(a,b,w)
area=0.5*sqrt(w(1)**2+w(2)**2+w(3)**2)
do i=1,3
   w(i)=w(i)/2./area
enddo
l(1)=sqrt((quant%xyz(1,node3)-quant%xyz(1,node2))**2+(quant%xyz(2,node3)&
     	-quant%xyz(2,node2))**2+(quant%xyz(3,node3)-quant%xyz(3,node2))**2)
l(2)=sqrt((quant%xyz(1,node3)-quant%xyz(1,node1))**2+(quant%xyz(2,node3)&
     	-quant%xyz(2,node1))**2+(quant%xyz(3,node3)-quant%xyz(3,node1))**2)
l(3)=sqrt((quant%xyz(1,node1)-quant%xyz(1,node2))**2+(quant%xyz(2,node1)&
     	-quant%xyz(2,node2))**2+(quant%xyz(3,node1)-quant%xyz(3,node2))**2)
do i=1,3
   u(i)=a(i)/l(3)
enddo
 call curl(w,u,v)
 call scalar(u,b,u3)
v3=2.*area/l(3)

b(1)=xi-quant%xyz(1,node1)
b(2)=yi-quant%xyz(2,node1)
b(3)=zi-quant%xyz(3,node1)
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
   s1(i)=(quant%xyz(i,node3)-quant%xyz(i,node2))/l(1)
   s2(i)=(quant%xyz(i,node1)-quant%xyz(i,node3))/l(2)
   s3(i)=(quant%xyz(i,node2)-quant%xyz(i,node1))/l(3)
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


subroutine current_node_patch_mapping(JMflag,string,curr,msh,quant,ptree)

    use BPACK_DEFS
    implicit none

    integer JMflag,patch, edge,edge1, node_patch(3), node_edge, node, info_idx
    integer i,j,k,ii,jj,kk,flag,lr
    real(kind=8) center(3), r, a, signs, current_patch(3),  current_abs
	character(*)::string
    complex(kind=8)::curr(:)
	type(quant_EMSURF)::quant
	type(mesh)::msh
	type(proctree)::ptree

    real(kind=8),allocatable :: current_at_patch(:),vec_current_at_patch(:,:), current_at_node(:)
	integer ierr

	allocate (current_at_node(quant%maxnode),current_at_patch(quant%maxpatch))
	allocate (vec_current_at_patch(3,quant%maxpatch))
	current_at_node=0
	current_at_patch=0
	vec_current_at_patch=0

	!$omp parallel do default(shared) private(lr,patch,signs,info_idx,i,j,center,edge,edge1,current_abs,current_patch,a,r)
	do edge=msh%idxs,msh%idxe
		do lr=3,4
			if((msh%new2old(edge)<=quant%Nunk_int+quant%Nunk_port .and. JMflag==0) .or. (msh%new2old(edge)>quant%Nunk_int+quant%Nunk_port .and. JMflag==1))then  ! the electric current (JMflag=0) or magnetic current (JMflag=1)
				if(JMflag==0)edge1=msh%new2old(edge)
				if(JMflag==1)edge1=msh%new2old(edge)-quant%Nunk_port
				patch = quant%info_unk(lr,edge1)
				if(patch/=-1)then
					if(lr==3)then
						signs=1
						info_idx=5
					endif
					if(lr==4)then
						signs=-1
						info_idx=6
					endif

					do i=1,3
						center(i)=1./3.*(quant%xyz(i,quant%node_of_patch(1,patch))+quant%xyz(i,quant%node_of_patch(2,patch))+quant%xyz(i,quant%node_of_patch(3,patch)))
					enddo
					current_abs=dble(curr(edge-msh%idxs+1))
					r=0.
					do j=1,3
						a=quant%xyz(j,quant%info_unk(info_idx,edge1))-center(j)
						r=r+a**2
						current_patch(j)=signs*current_abs*a
					enddo
					r=sqrt(r)
					current_patch = current_patch/r
					do i=1,3
					!$omp atomic
					vec_current_at_patch(i,patch)=vec_current_at_patch(i,patch) + current_patch(i)
					!$omp end atomic
					enddo
				endif
			endif
		enddo
	enddo
	!$omp end parallel do
	call MPI_ALLREDUCE(MPI_IN_PLACE,vec_current_at_patch,quant%maxpatch*3,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)


	do i=1,quant%maxpatch
	current_at_patch(i) = sqrt(sum(vec_current_at_patch(1:3,i)**2))
	enddo

    !$omp parallel do default(shared) private(patch,i,ii,node,a)
    do node=1, quant%maxnode
        ii=0; a=0.
        do patch=1, quant%maxpatch
            do i=1,3
                if (quant%node_of_patch(i,patch)==node) then
                    a=a+current_at_patch(patch)
                    ii=ii+1
                endif
            enddo
        enddo
        current_at_node(node)=a/dble(ii)
    enddo
    !$omp end parallel do

	if(ptree%MyID == Main_ID)then
    ! string='current'//chara//'.out'
    open(30,file=string)
    do node=1, quant%maxnode
        write (30,*) current_at_node(node)
    enddo
    close(30)
	endif

    deallocate (current_at_node,current_at_patch,vec_current_at_patch)
    return

end subroutine current_node_patch_mapping


real(kind=8) function triangle_area(patch,quant)

    use BPACK_DEFS
    implicit none

    integer patch,i
    real(kind=8) a(3),b(3),c(3)
	type(quant_EMSURF)::quant

    do i=1,3
        a(i)=quant%xyz(i,quant%node_of_patch(2,patch))-quant%xyz(i,quant%node_of_patch(1,patch))
        b(i)=quant%xyz(i,quant%node_of_patch(3,patch))-quant%xyz(i,quant%node_of_patch(1,patch))
    enddo

    call curl(a,b,c)
    triangle_area=0.5*sqrt(c(1)**2+c(2)**2+c(3)**2)

    return
end function triangle_area


logical function in_triangle(point,patch,quant)

    use BPACK_DEFS
    implicit none

    integer patch,i
    real(kind=8) a(3),b(3),c(3),area,point(3),alpha1,beta1,gamma1
	type(quant_EMSURF)::quant

    do i=1,3
        a(i)=quant%xyz(i,quant%node_of_patch(2,patch))-quant%xyz(i,quant%node_of_patch(1,patch))
        b(i)=quant%xyz(i,quant%node_of_patch(3,patch))-quant%xyz(i,quant%node_of_patch(1,patch))
    enddo
    call curl(a,b,c)
	area=0.5*sqrt(sum(c**2d0))

    do i=1,3
        a(i)=quant%xyz(i,quant%node_of_patch(2,patch))-point(i)
        b(i)=quant%xyz(i,quant%node_of_patch(3,patch))-point(i)
    enddo
    call curl(a,b,c)
    alpha1=0.5*sqrt(sum(c**2d0))/area

    do i=1,3
        a(i)=quant%xyz(i,quant%node_of_patch(2,patch))-point(i)
        b(i)=quant%xyz(i,quant%node_of_patch(1,patch))-point(i)
    enddo
    call curl(a,b,c)
	beta1=0.5*sqrt(sum(c**2d0))/area


    do i=1,3
        a(i)=quant%xyz(i,quant%node_of_patch(1,patch))-point(i)
        b(i)=quant%xyz(i,quant%node_of_patch(3,patch))-point(i)
    enddo
    call curl(a,b,c)
	gamma1=0.5*sqrt(sum(c**2d0))/area

	in_triangle = alpha1>=0 .and. alpha1<=1 .and. beta1>=0 .and. beta1<=1 .and. gamma1>=0 .and. gamma1<=1 .and. abs(alpha1+beta1+gamma1-1)<=1d-13

    return
end function in_triangle


subroutine element_Vinc_VV_SURF(theta,phi,edge,value,quant)

    use BPACK_DEFS
    implicit none

    integer edge
    complex(kind=8) value
    real(kind=8) theta, phi
    integer i,ii,jj,node1,node2,node3
	real(kind=8) lm,einc(3),a(3),nr(3),k(3)
	real(kind=8), allocatable::x(:),y(:),z(:),w(:)
    complex(kind=8)  phase,ctemp_e,ctemp_h,value_e,value_h,ee(3),hh(3),hh1(3)
	type(quant_EMSURF)::quant

	allocate(x(quant%integral_points))
	allocate(y(quant%integral_points))
	allocate(z(quant%integral_points))
	allocate(w(quant%integral_points))


    einc(1)=cos(theta*pi/180.)*cos(phi*pi/180.)
    einc(2)=cos(theta*pi/180.)*sin(phi*pi/180.)
    einc(3)=-sin(theta*pi/180.)
    k(1)=sin(theta*pi/180.)*cos(phi*pi/180.)
    k(2)=sin(theta*pi/180.)*sin(phi*pi/180.)
    k(3)=cos(theta*pi/180.)

    value_h=(0,0); value_e=(0,0)

    node1=quant%info_unk(1,edge)
    node2=quant%info_unk(2,edge)
    lm=sqrt((quant%xyz(1,node1)-quant%xyz(1,node2))**2+(quant%xyz(2,node1)-quant%xyz(2,node2))**2+(quant%xyz(3,node1)-quant%xyz(3,node2))**2)
    do jj=3,4
	   nr(1:3)=quant%normal_of_patch(1:3,quant%info_unk(jj,edge))
	   node3=quant%info_unk(jj+2,edge)
       call gau_grobal(edge,jj,x,y,z,w,quant)
	   do ii=1,quant%integral_points
          phase=junit*quant%wavenum*(sin(theta*pi/180.)*cos(phi*pi/180.)*x(ii)+sin(theta*pi/180.)*sin(phi*pi/180.)*y(ii)+cos(theta*pi/180.)*z(ii))
	      do i=1,3
		     ee(i)=einc(i)*exp(phase)
          end do
          call ccurl(-k,ee,hh1)
          call ccurl(nr,hh1,hh)
          a(1)=x(ii)-quant%xyz(1,node3)
	      a(2)=y(ii)-quant%xyz(2,node3)
	      a(3)=z(ii)-quant%xyz(3,node3)
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

subroutine element_Vinc_HH_SURF(theta,phi,edge,value,quant)

    use BPACK_DEFS
    implicit none

    type(quant_EMSURF)::quant
    integer edge
    complex(kind=8) value
    real(kind=8) theta, phi
    integer i,ii,jj,node1,node2,node3
	real(kind=8) lm,einc(3),a(3),nr(3),k(3)
	real(kind=8), allocatable::x(:),y(:),z(:),w(:)
    complex(kind=8)  phase,ctemp_e,ctemp_h,value_e,value_h,ee(3),hh(3),hh1(3)

	allocate(x(quant%integral_points))
	allocate(y(quant%integral_points))
	allocate(z(quant%integral_points))
	allocate(w(quant%integral_points))

    einc(1)=-sin(phi*pi/180.)
    einc(2)=cos(phi*pi/180.)
    einc(3)=0.
    k(1)=sin(theta*pi/180.)*cos(phi*pi/180.)
    k(2)=sin(theta*pi/180.)*sin(phi*pi/180.)
    k(3)=cos(theta*pi/180.)

    value_h=(0,0); value_e=(0,0)

    node1=quant%info_unk(1,edge)
    node2=quant%info_unk(2,edge)
    lm=sqrt((quant%xyz(1,node1)-quant%xyz(1,node2))**2+(quant%xyz(2,node1)-quant%xyz(2,node2))**2+(quant%xyz(3,node1)-quant%xyz(3,node2))**2)
    do jj=3,4
	   nr(1:3)=quant%normal_of_patch(1:3,quant%info_unk(jj,edge))
	   node3=quant%info_unk(jj+2,edge)
       call gau_grobal(edge,jj,x,y,z,w,quant)
	   do ii=1,quant%integral_points
          phase=junit*quant%wavenum*(sin(theta*pi/180.)*cos(phi*pi/180.)*x(ii)+sin(theta*pi/180.)*sin(phi*pi/180.)*y(ii)+cos(theta*pi/180.)*z(ii))
	      do i=1,3
		     ee(i)=einc(i)*exp(phase)
          end do
          call ccurl(-k,ee,hh1)
          call ccurl(nr,hh1,hh)
          a(1)=x(ii)-quant%xyz(1,node3)
	      a(2)=y(ii)-quant%xyz(2,node3)
	      a(3)=z(ii)-quant%xyz(3,node3)
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
	use BPACK_DEFS
    implicit none

    real(kind=8) rcs
    complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1
    real(kind=8) theta,phi,dtheta,dphi

    integer i,j,ii,jj,iii,jjj,patch,flag
    real(kind=8) l_edge,l_edgefine
    real(kind=8) a(3)
    real(kind=8),allocatable:: x(:),y(:),z(:),w(:)
    complex(kind=8):: curr(:,:)
	type(mesh)::msh
	type(quant_EMSURF)::quant
	type(proctree)::ptree

    integer edge,edge_m,edge_n,ierr

    if(ptree%MyID==Main_ID)open (100, file='VV_bistatic.txt')
    theta=90.
    dphi=180./quant%RCS_Nsample

    do i=0, quant%RCS_Nsample
        phi=i*dphi
        ctemp_loc=(0.,0.)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
        do edge=msh%idxs,msh%idxe
            call VV_polar_SURF(theta,phi,msh%new2old(edge),ctemp_1,curr(edge-msh%idxs+1,1),quant)
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
            call HH_polar_SURF(theta,phi,msh%new2old(edge),ctemp_1,curr(edge-msh%idxs+1,2),quant)
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



    return

end subroutine RCS_bistatic_SURF

subroutine VV_polar_SURF(theta,phi,edge,ctemp_1,curr,quant)

    use BPACK_DEFS
    implicit none

    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1
    real(kind=8) theta,phi
    type(quant_EMSURF)::quant
    integer i,j,ii,jj,iii,jjj,patch,flag
    real(kind=8) l_edge,l_edgefine
    real(kind=8) a(3)
    complex(kind=8)::curr
    integer edge,edge_m,edge_n
    ! type(mesh)::quant
	real(kind=8),allocatable:: x(:),y(:),z(:),w(:)

	allocate(x(quant%integral_points))
	allocate(y(quant%integral_points))
	allocate(z(quant%integral_points))
	allocate(w(quant%integral_points))

            ctemp_rcs(1:3)=(0,0)
            l_edge=sqrt((quant%xyz(1,quant%info_unk(1,edge))-quant%xyz(1,quant%info_unk(2,edge)))**2+(quant%xyz(2,quant%info_unk(1,edge))-quant%xyz(2,quant%info_unk(2,edge)))**2+(quant%xyz(3,quant%info_unk(1,edge))-quant%xyz(3,quant%info_unk(2,edge)))**2)
			do jj=3,4  ! two triganle
	            !patch=quant%info_unk(jj,edgefine_m)
                call gau_grobal(edge,jj,x,y,z,w,quant)

		        do ii=1,quant%integral_points

				    a(1)=x(ii)-quant%xyz(1,quant%info_unk(jj+2,edge))
		            a(2)=y(ii)-quant%xyz(2,quant%info_unk(jj+2,edge))
		            a(3)=z(ii)-quant%xyz(3,quant%info_unk(jj+2,edge))
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

subroutine HH_polar_SURF(theta,phi,edge,ctemp_1,curr,quant)

    use BPACK_DEFS
    implicit none

    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1
    real(kind=8) theta,phi
    type(quant_EMSURF)::quant
    integer i,j,ii,jj,iii,jjj,patch,flag
    real(kind=8) l_edge,l_edgefine
    real(kind=8) a(3)
    complex(kind=8)::curr
	real(kind=8),allocatable:: x(:),y(:),z(:),w(:)

    integer edge,edge_m,edge_n

	allocate(x(quant%integral_points))
	allocate(y(quant%integral_points))
	allocate(z(quant%integral_points))
	allocate(w(quant%integral_points))

        ctemp_rcs(1:3)=(0,0)
        l_edge=sqrt((quant%xyz(1,quant%info_unk(1,edge))-quant%xyz(1,quant%info_unk(2,edge)))**2+(quant%xyz(2,quant%info_unk(1,edge))-quant%xyz(2,quant%info_unk(2,edge)))**2+(quant%xyz(3,quant%info_unk(1,edge))-quant%xyz(3,quant%info_unk(2,edge)))**2)

			do jj=3,4  ! two triganle
	            !patch=quant%info_unk(jj+2,edgefine_m)
                call gau_grobal(edge,jj,x,y,z,w,quant)
		        do ii=1,quant%integral_points

				    a(1)=x(ii)-quant%xyz(1,quant%info_unk(jj+2,edge))  !free node coordinate
		            a(2)=y(ii)-quant%xyz(2,quant%info_unk(jj+2,edge))
		            a(3)=z(ii)-quant%xyz(3,quant%info_unk(jj+2,edge))
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

    use BPACK_DEFS
    implicit none
    complex(kind=8)::curr(:)
    real(kind=8) rcs
    complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1,ctemp_2
    real(kind=8) dsita,dphi
    integer edge,edge_m,edge_n,ierr
	type(mesh)::msh
    type(quant_EMSURF)::quant
	type(proctree)::ptree

    ctemp_loc=(0.,0.)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
        do edge=msh%idxs,msh%idxe
            call VV_polar_SURF(dsita,dphi,msh%new2old(edge),ctemp_1,curr(edge-msh%idxs+1),quant)
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

    use BPACK_DEFS
    implicit none

    real(kind=8) rcs
    complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1,ctemp_2
    real(kind=8) dsita,dphi
    integer edge,edge_m,edge_n,ierr
    complex(kind=8)::curr(:)
	type(mesh)::msh
	type(quant_EMSURF)::quant
	type(proctree)::ptree

    ctemp_loc=(0.,0.)
        rcs=0
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
        do edge=msh%idxs,msh%idxe
            call HH_polar_SURF(dsita,dphi,msh%new2old(edge),ctemp_1,curr(edge-msh%idxs+1),quant)
	        ctemp_loc=ctemp_loc+ctemp_1
        enddo
        !$omp end parallel do

		call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)

        rcs=(abs(quant%wavenum*ctemp))**2/4/pi
        !rcs=rcs/quant%wavelength
        rcs=10*log10(rcs)

    return

end subroutine RCS_monostatic_HH_SURF




subroutine geo_modeling_SURF(quant,MPIcomm,DATA_DIR)
    use BPACK_DEFS
	use MISC_Utilities
    implicit none

	real T0
    integer i,j,ii,jj,nn,iii,jjj,pp,cnt,mode,mm,rr
    integer intemp
    integer node, patch, edge,edge1, flag
    integer node1, node2,found,npolar,off
    integer node_temp(2)
    real(kind=8) a(3),b(3),c(3),r0,x1,x2,y1,y2,z1,z2
	! type(mesh)::msh
	type(quant_EMSURF)::quant
	! type(proctree)::ptree
	integer MPIcomm,MyID,ierr
	CHARACTER (*) DATA_DIR
	integer Maxedge
	integer,parameter:: NperNode=10
	integer,allocatable::tmpint(:,:)
	integer,allocatable:: info_unk(:,:)
	real(kind=8),allocatable::distance(:)
	integer,allocatable::order(:)
	integer v1,v2

	call MPI_Comm_rank(MPIcomm,MyID,ierr)

    open(11,file=trim(DATA_DIR)//'_node.inp')
    open(111,file=trim(DATA_DIR)//'_elem.inp')

    read(11,*)quant%maxnode
    read(111,*)quant%maxpatch
    Maxedge=quant%maxpatch*3



    allocate(quant%xyz(3,quant%maxnode+Maxedge))
    allocate(quant%node_of_patch(0:3,quant%maxpatch),quant%info_unk(0:6,maxedge))
	quant%info_unk=-1
	quant%node_of_patch=-1
    allocate(quant%normal_of_patch(3,quant%maxpatch))
	allocate(quant%edge_of_node(quant%maxnode))
	do i=1,quant%maxnode
		quant%edge_of_node(i)%Nedge=0
	enddo

    !************quant%xyz****************
    i=1
    do while(i<=quant%maxnode)
        read(11,*)intemp,quant%xyz(1:3,i)
        quant%xyz(1:3,i)=quant%xyz(1:3,i)/quant%scaling
        i=i+1
    enddo
    close(11)

    i=1
    if (quant%mesh_normal==1) then
        do while(i<=quant%maxpatch)
            read(111,*)intemp,quant%node_of_patch(1:3,i)
            i=i+1
        enddo
    elseif (quant%mesh_normal==-1) then
        do while(i<=quant%maxpatch)
            read(111,*)intemp,quant%node_of_patch(3,i),quant%node_of_patch(2,i),quant%node_of_patch(1,i)
            i=i+1
        enddo
    endif
    close(111)

    !************quant%normal_of_patch****************

    !$omp parallel do default(shared) private(patch,i,a,b,c,r0)
    do patch=1,quant%maxpatch
        do i=1,3
            a(i)=(quant%xyz(i,quant%node_of_patch(2,patch))-quant%xyz(i,quant%node_of_patch(1,patch)))
            b(i)=(quant%xyz(i,quant%node_of_patch(3,patch))-quant%xyz(i,quant%node_of_patch(1,patch)))
        enddo
        call curl(a,b,c)
        r0=sqrt(c(1)**2+c(2)**2+c(3)**2)
        c(1)=c(1)/r0
        c(2)=c(2)/r0
        c(3)=c(3)/r0
        quant%normal_of_patch(1:3,patch)=c(1:3)
    enddo
    !$omp end parallel do

    !************quant%info_unk****************

    edge=0
    do i=1,quant%maxpatch
		do ii=1,3
			jj=ii+1
			if(jj==4)jj=1
			node_temp(1)=quant%node_of_patch(ii,i)
			node_temp(2)=quant%node_of_patch(jj,i)

			found=0
			do nn=1,quant%edge_of_node(node_temp(1))%Nedge
			if(quant%edge_of_node(node_temp(1))%edges(nn,2)==node_temp(2))then
				found=1
				edge1=quant%edge_of_node(node_temp(1))%edges(nn,1)
				exit
			endif
			enddo
			if(found==0)then
				edge=edge+1
                if(node_temp(1)<node_temp(2)) then
                    quant%info_unk(1,edge)=node_temp(1)
                    quant%info_unk(2,edge)=node_temp(2)
                else
                    quant%info_unk(1,edge)=node_temp(2)
                    quant%info_unk(2,edge)=node_temp(1)
                endif
                quant%info_unk(3,edge)=i
                quant%info_unk(0,edge)=0

				if(mod(quant%edge_of_node(node_temp(1))%Nedge,NperNode)==0)then
					allocate(tmpint(quant%edge_of_node(node_temp(1))%Nedge+NperNode,2))
					if(quant%edge_of_node(node_temp(1))%Nedge>0)then
						tmpint(1:quant%edge_of_node(node_temp(1))%Nedge,:)=quant%edge_of_node(node_temp(1))%edges
						deallocate(quant%edge_of_node(node_temp(1))%edges)
					endif
					allocate(quant%edge_of_node(node_temp(1))%edges(quant%edge_of_node(node_temp(1))%Nedge+NperNode,2))
					quant%edge_of_node(node_temp(1))%edges=tmpint
					deallocate(tmpint)
				endif
				quant%edge_of_node(node_temp(1))%Nedge=quant%edge_of_node(node_temp(1))%Nedge+1
				quant%edge_of_node(node_temp(1))%edges(quant%edge_of_node(node_temp(1))%Nedge,1)=edge
				quant%edge_of_node(node_temp(1))%edges(quant%edge_of_node(node_temp(1))%Nedge,2)=node_temp(2)
			else
				quant%info_unk(4,edge1)=i
			endif


			found=0
			do nn=1,quant%edge_of_node(node_temp(2))%Nedge
			if(quant%edge_of_node(node_temp(2))%edges(nn,2)==node_temp(1))then
				found=1
				edge1=quant%edge_of_node(node_temp(2))%edges(nn,1)
				exit
			endif
			enddo
			if(found==0)then
				if(mod(quant%edge_of_node(node_temp(2))%Nedge,NperNode)==0)then
					allocate(tmpint(quant%edge_of_node(node_temp(2))%Nedge+NperNode,2))
					if(quant%edge_of_node(node_temp(2))%Nedge>0)then
						tmpint(1:quant%edge_of_node(node_temp(2))%Nedge,:)=quant%edge_of_node(node_temp(2))%edges
						deallocate(quant%edge_of_node(node_temp(2))%edges)
					endif
					allocate(quant%edge_of_node(node_temp(2))%edges(quant%edge_of_node(node_temp(2))%Nedge+NperNode,2))
					quant%edge_of_node(node_temp(2))%edges=tmpint
					deallocate(tmpint)
				endif
				quant%edge_of_node(node_temp(2))%Nedge=quant%edge_of_node(node_temp(2))%Nedge+1
				quant%edge_of_node(node_temp(2))%edges(quant%edge_of_node(node_temp(2))%Nedge,1)=edge
				quant%edge_of_node(node_temp(2))%edges(quant%edge_of_node(node_temp(2))%Nedge,2)=node_temp(1)
			endif

		enddo
	enddo
	do ii=1,quant%maxnode
		if(quant%edge_of_node(ii)%Nedge>0)deallocate(quant%edge_of_node(ii)%edges)
	enddo
	deallocate(quant%edge_of_node)

	allocate(info_unk(0:6,edge))
	edge=0
	do ii=1,size(info_unk,2)
		if(quant%info_unk(4,ii)/=-1)then
			edge=edge+1
			info_unk(0:6,edge)=quant%info_unk(0:6,ii)
			if(info_unk(3,edge)>info_unk(4,edge))then
				intemp= info_unk(3,edge)
				info_unk(3,edge)=info_unk(4,edge)
				info_unk(4,edge)=intemp
			endif
		endif
	enddo
	Maxedge=edge
	deallocate(quant%info_unk)
	allocate(quant%info_unk(0:6,Maxedge))
	quant%info_unk(0:6,1:Maxedge)=info_unk(0:6,1:Maxedge)
	deallocate(info_unk)


    ! edge=0
    ! do i=1,quant%maxpatch-1
        ! do j=i+1,quant%maxpatch
            ! flag=0;node1=0;node2=0;iii=1
            ! do ii=1,3
                ! do jj=1,3
	     	         ! if(quant%node_of_patch(ii,i)==quant%node_of_patch(jj,j))then
                        ! flag=flag+1
                        ! node_temp(iii)=quant%node_of_patch(ii,i)
                        ! iii=iii+1
                    ! endif
                ! enddo
            ! enddo
            ! if(flag==2)then
                ! edge=edge+1
                ! if(node_temp(1)<node_temp(2)) then
                    ! quant%info_unk(1,edge)=node_temp(1)
                    ! quant%info_unk(2,edge)=node_temp(2)
                ! else
                    ! quant%info_unk(1,edge)=node_temp(2)
                    ! quant%info_unk(2,edge)=node_temp(1)
                ! endif
                ! quant%info_unk(3,edge)=i
                ! quant%info_unk(4,edge)=j       ! notice that : i<j
                ! quant%info_unk(0,edge)=0
            ! endif
        ! enddo
    ! enddo

    ! Maxedge=edge

    !$omp parallel do default(shared) private(edge,node_temp,jj,iii,jjj)
    do edge=1,maxedge
	    node_temp(1)=0
	    node_temp(2)=0
	    do jj=3,4
             do iii=1,3
                 do jjj=1,2
        	            if(quant%node_of_patch(iii,quant%info_unk(jj,edge))==quant%info_unk(jjj,edge)) then
                           node_temp(jj-2)=node_temp(jj-2)+iii
        	            endif
                 enddo
             enddo
         enddo
         node_temp(1)=6-node_temp(1)
         node_temp(2)=6-node_temp(2)
         quant%info_unk(5,edge)=quant%node_of_patch(node_temp(1),quant%info_unk(3,edge))
         quant%info_unk(6,edge)=quant%node_of_patch(node_temp(2),quant%info_unk(4,edge))
    enddo
    !$omp end parallel do


	! handle the ports
	allocate(distance(Maxedge))
	distance=0
	allocate(order(Maxedge))
	do edge=1,Maxedge
		v1 = quant%info_unk(5,edge)
		v2 = quant%info_unk(6,edge)
		do ii =1,quant%Nport
			if(quant%ports(ii)%type==0)then
				if(abs(dot_product(quant%xyz(:,v1)-quant%ports(ii)%origin, quant%ports(ii)%z))<SafeEps .and. sqrt(sum((quant%xyz(:,v1)-quant%ports(ii)%origin)**2d0))<quant%ports(ii)%R*(1+SafeEps)  .and. abs(dot_product(quant%xyz(:,v2)-quant%ports(ii)%origin, quant%ports(ii)%z))<SafeEps .and. sqrt(sum((quant%xyz(:,v2)-quant%ports(ii)%origin)**2d0))<quant%ports(ii)%R*(1+SafeEps))then
					distance(edge)=ii
					quant%ports(ii)%Nunk=quant%ports(ii)%Nunk+1
					quant%Nunk_port = quant%Nunk_port+1
				endif
			else if(quant%ports(ii)%type==1)then
				x1 = dot_product(quant%xyz(:,v1)-quant%ports(ii)%origin, quant%ports(ii)%x)
				y1 = dot_product(quant%xyz(:,v1)-quant%ports(ii)%origin, quant%ports(ii)%y)
				z1 = dot_product(quant%xyz(:,v1)-quant%ports(ii)%origin, quant%ports(ii)%z)
				x2 = dot_product(quant%xyz(:,v2)-quant%ports(ii)%origin, quant%ports(ii)%x)
				y2 = dot_product(quant%xyz(:,v2)-quant%ports(ii)%origin, quant%ports(ii)%y)
				z2 = dot_product(quant%xyz(:,v2)-quant%ports(ii)%origin, quant%ports(ii)%z)

				if(abs(z1)<SafeEps .and. x1<quant%ports(ii)%a*(1+SafeEps) .and. y1<quant%ports(ii)%b*(1+SafeEps) .and. abs(z2)<SafeEps .and. x2<quant%ports(ii)%a*(1+SafeEps) .and. y2<quant%ports(ii)%b*(1+SafeEps))then
					distance(edge)=ii
					quant%ports(ii)%Nunk=quant%ports(ii)%Nunk+1
					quant%Nunk_port = quant%Nunk_port+1
				endif
			else
				write(*,*)'unrecognized port type',quant%ports(ii)%type
				stop
			endif
		enddo
	enddo
	call quick_sort(distance, order, Maxedge)
	allocate(info_unk(0:6,edge))
	info_unk(0:6,1:Maxedge)=quant%info_unk(0:6,1:Maxedge)
	quant%info_unk(0:6,:)=info_unk(0:6,order)
	deallocate(info_unk)
	deallocate(distance)
	deallocate(order)
	quant%Nunk_int = Maxedge - quant%Nunk_port
	quant%Nunk = quant%Nunk_int + quant%Nunk_port*2



    node=quant%maxnode
    do edge=1, Maxedge
        node=node+1
        quant%info_unk(0,edge)=node
        do i=1,3
            quant%xyz(i,node)=1./2.*(quant%xyz(i,quant%info_unk(1,edge))+quant%xyz(i,quant%info_unk(2,edge)))
        enddo
    enddo

	quant%maxedgelength = 0
	do edge=1,Maxedge
		quant%maxedgelength = max(quant%maxedgelength,sqrt(sum(abs(quant%xyz(:,quant%info_unk(1,edge))-quant%xyz(:,quant%info_unk(2,edge)))**2)))
		if(sqrt(sum(abs(quant%xyz(:,quant%info_unk(1,edge))-quant%xyz(:,quant%info_unk(2,edge)))**2))>0.9d0)write(*,*)edge,sqrt(sum(abs(quant%xyz(:,quant%info_unk(1,edge))-quant%xyz(:,quant%info_unk(2,edge)))**2)),quant%info_unk(1,edge),quant%info_unk(2,edge),quant%xyz(:,quant%info_unk(1,edge)),quant%xyz(:,quant%info_unk(2,edge))
	end do


	quant%minedgelength = BigValue
	do edge=1,Maxedge
		quant%minedgelength = min(quant%minedgelength,sqrt(sum(abs(quant%xyz(:,quant%info_unk(1,edge))-quant%xyz(:,quant%info_unk(2,edge)))**2)))
	end do

	! write(*,*)	quant%xyz(1,1:100),sum(quant%xyz(1,:))
	! stop

    if(MyID==Main_ID)write (*,*) ''
    if(MyID==Main_ID)write (*,*) 'Maxedge:',Maxedge
    if(MyID==Main_ID)write (*,*) 'Nunk:',quant%Nunk
	if(MyID==Main_ID)write (*,*) 'minedgelength:',quant%minedgelength
	if(MyID==Main_ID)write (*,*) 'wavelength/minedgelength:',quant%wavelength/quant%minedgelength
	if(MyID==Main_ID)write (*,*) 'maxedgelength:',quant%maxedgelength
	if(MyID==Main_ID)write (*,*) 'wavelength/maxedgelength:',quant%wavelength/quant%maxedgelength


	if(MyID==Main_ID)write (*,*) 'Pre-computing the nxe_dot_rwg vectors at the ports ...'
	T0=secnds(0.0)
	edge=quant%Nunk_int
	do pp=1,quant%Nport
		if(quant%ports(pp)%type==0)then
			npolar=2
			off=0
		elseif(quant%ports(pp)%type==1)then
			npolar=1
			off=1
		else
			write(*,*)'unrecognized port type',quant%ports(pp)%type
			stop
		endif

		allocate(quant%ports(pp)%nxe_dot_rwg(quant%ports(pp)%Nunk,quant%ports(pp)%nmax+1,quant%ports(pp)%mmax+off,2,npolar))
		allocate(quant%ports(pp)%e_dot_rwg(quant%ports(pp)%Nunk,quant%ports(pp)%nmax+1,quant%ports(pp)%mmax+off,2,npolar))
		do cnt=1,quant%ports(pp)%Nunk
			do rr=1,npolar
				do nn=0,quant%ports(pp)%nmax
					do mm=1-off,quant%ports(pp)%mmax
						if(mm>0 .or. nn>0)then
						call Port_nxe_dot_rwg(edge+cnt,pp,mm,nn,1,rr,quant%ports(pp)%nxe_dot_rwg(cnt,nn+1,mm+off,1,rr),quant)
						call Port_e_dot_rwg(edge+cnt,pp,mm,nn,1,rr,quant%ports(pp)%e_dot_rwg(cnt,nn+1,mm+off,1,rr),quant)
						endif
					enddo
				enddo
				do nn=0,quant%ports(pp)%nmax
					do mm=1-off,quant%ports(pp)%mmax
						if(mm>0 .or. nn>0)then
						call Port_nxe_dot_rwg(edge+cnt,pp,mm,nn,2,rr,quant%ports(pp)%nxe_dot_rwg(cnt,nn+1,mm+off,2,rr),quant)
						call Port_e_dot_rwg(edge+cnt,pp,mm,nn,2,rr,quant%ports(pp)%e_dot_rwg(cnt,nn+1,mm+off,2,rr),quant)
						endif
					enddo
				enddo
			enddo
		enddo
		edge =edge+quant%ports(pp)%Nunk
	enddo
	if(MyID==Main_ID)write (*,*) 'Pre-computing the nxe_dot_rwg vectors at the ports:',secnds(T0),'Seconds'
	if(MyID==Main_ID)write (*,*) ''

    return

end subroutine geo_modeling_SURF


subroutine EM_solve_SURF(bmat,option,msh,quant,ptree,stats)
    use BPACK_DEFS
	use BPACK_Solve_Mul


    implicit none

    integer i, j, ii, jj, iii, jjj,ierr
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, n_iter_max, iter ,N_unk, N_unk_loc
    real(kind=8) theta, phi, dphi, rcs_V, rcs_H
    real T0
    real(kind=8) n1,n2,rtemp
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:)
	real(kind=8):: rel_error
	type(Hoption)::option
	type(Bmatrix)::bmat
	type(mesh)::msh
	type(quant_EMSURF)::quant
	type(proctree)::ptree
	type(Hstat)::stats
	complex(kind=8),allocatable:: current(:,:),voltage(:,:)

	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve......"

	N_unk_loc = msh%idxe-msh%idxs+1

    if (quant%RCS_static==2) then

        theta=90
        phi=0

        allocate (current(N_unk_loc,2))
		Current=0
        allocate (voltage(N_unk_loc,2))


        !$omp parallel do default(shared) private(edge,value_Z)
        do edge=msh%idxs, msh%idxe
            call element_Vinc_VV_SURF(theta,phi,msh%new2old(edge),value_Z,quant)
			voltage(edge-msh%idxs+1,1)=value_Z
			call element_Vinc_HH_SURF(theta,phi,msh%new2old(edge),value_Z,quant)
            voltage(edge-msh%idxs+1,2)=value_Z
        enddo
        !$omp end parallel do

        T0=secnds(0.0)

		call BPACK_Solution(bmat,Current,Voltage,N_unk_loc,2,option,ptree,stats)


        if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''
        if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'Solving:',secnds(T0),'Seconds'
        if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''

		T0=secnds(0.0)
        call RCS_bistatic_SURF(Current,msh,quant,ptree)

		! call current_node_patch_mapping('V',curr(:,1),msh)
		! call current_node_patch_mapping('H',curr(:,2),msh)

        if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''
        if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
        if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''
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
				call element_Vinc_HH_SURF(theta,phi,msh%new2old(edge),value_Z,quant)
				b(edge-msh%idxs+1,j+1)=value_Z
			enddo
			!$omp end parallel do
        enddo

		call BPACK_Solution(bmat,x,b,N_unk_loc,num_sample+1,option,ptree,stats)

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
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,'(A13Es14.2)') 'Solve flops:',rtemp


		deallocate(Current)

    endif

	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve finished"
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "

    return

end subroutine EM_solve_SURF





subroutine EM_cavity_postprocess(option,msh,quant,ptree,stats,eigvec,nth,norm,eigval,Enormal_GF,ith,model)
    use BPACK_DEFS
	use BPACK_Solve_Mul


    implicit none

    integer i, j, ii, jj, iii, jjj,ierr,nth, ith
    integer level, blocks, edge, edge_n, patch, node, group, pp, ppn, port, cntn
    integer rank, index_near, m, n, length, flag, num_sample, n_iter_max, iter ,N_unk, N_unk_loc
    real(kind=8) theta, phi, dphi, rcs_V, rcs_H, freq, norm
    real T0
    real(kind=8) n1,n2,rtemp, center(3),an(3)
    complex(kind=8) value_Z,ctemp,a,eigval
    complex(kind=8) eigvec(:),Enormal_GF(:)
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:)
	real(kind=8):: rel_error,Rs, max_Enormal_GF
	type(Hoption)::option
	type(mesh)::msh
	type(quant_EMSURF)::quant
	type(proctree)::ptree
	type(Hstat)::stats
	integer,allocatable:: port_of_patch(:),cnt_patch(:)
	complex(kind=8),allocatable:: current(:,:),voltage(:,:),Enormal_at_patch(:),Enormal_at_node(:),Ht_at_patch(:,:),Et_at_patch(:,:)
	complex(kind=8):: field(3),cpoint(3), volt_acc
	real(kind=8) ln, area,point(3), Ht2int, dx(3)
	character(4096)::string
	real(kind=8),allocatable::xn(:),yn(:),zn(:),wn(:), weight_at_patch(:),ExH_at_ports(:,:)
	character(len=1024)  :: substring,substring1,model

	write(substring , *) nth
	write(substring1 , *) quant%freq


	n1 = OMP_get_wtime()
	quant%obs_Efields=0

	do n=msh%idxs, msh%idxe
		do ii =1,quant%Nobs
			field=0
			edge = msh%new2old(n)
			call Field_EMSURF(quant%obs_points(:,ii),field,edge,quant)
			quant%obs_Efields(:,ii) = quant%obs_Efields(:,ii) + field*eigvec(n-msh%idxs+1)*impedence0 ! the solution vector is J and M/impedence0, this makes it easier to compare with ie3deigen
		enddo
	enddo

	call MPI_ALLREDUCE(MPI_IN_PLACE,quant%obs_Efields,quant%Nobs*3,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)

	if(ptree%MyID==Main_ID)then
		volt_acc = 0
		do ii =2,quant%Nobs
			dx =  quant%obs_points(:,ii) - quant%obs_points(:,ii-1)
			volt_acc = volt_acc + dot_product(quant%obs_Efields(:,ii),dx)
		enddo
		quant%normalize_factor=abs(volt_acc)


		quant%obs_Efields=quant%obs_Efields/quant%normalize_factor
		volt_acc=volt_acc/quant%normalize_factor

	! write(*,*)abs(quant%obs_Efields(1,quant%Nobs/2))
	! write(*,*)'x'

	! write(*,*)abs(quant%obs_Efields(2,quant%Nobs/2))
	! write(*,*)'y'

	! write(*,*)abs(quant%obs_Efields(3,quant%Nobs/2))
	! write(*,*)'z'
	open(28,file=trim(adjustl(model))//'_Eobs_'//trim(adjustl(substring))//'_freq_'//trim(adjustl(substring1))//'.out')
	do ii=1,quant%Nobs
		write(28,*) quant%obs_points(:,ii),abs(quant%obs_Efields(1,ii)), abs(quant%obs_Efields(2,ii)), abs(quant%obs_Efields(3,ii))
	enddo
	close(28)
	write(*,*)ith,': Computing_nth_'//trim(adjustl(substring))//'_freq_'//trim(adjustl(substring1))
	write(*,*)'   eigval: ', eigval
	write(*,*)'   norm_1/norm_inf: ', norm
	write(*,*)'   acceleration voltage: ', abs(volt_acc)

	endif
	n2 = OMP_get_wtime()
	! if(ptree%MyID==Main_ID)write(*,*)n2-n1,' seconds'




	!!!!! Compute the maximum normal electric fields, this is not very accurate due to definition of RWG
	n1 = OMP_get_wtime()
	allocate(Enormal_at_patch(quant%maxpatch))
	Enormal_at_patch=0

	do edge=msh%idxs,msh%idxe
		edge_n = msh%new2old(edge)
		if(edge_n<=quant%Nunk-quant%Nunk_port)then
			ln=sqrt((quant%xyz(1,quant%info_unk(1,edge_n))-quant%xyz(1,quant%info_unk(2,edge_n)))**2+(quant%xyz(2,quant%info_unk(1,edge_n))-quant%xyz(2,quant%info_unk(2,edge_n)))**2+(quant%xyz(3,quant%info_unk(1,edge_n))-quant%xyz(3,quant%info_unk(2,edge_n)))**2)

			do jj=3,4
			patch = quant%info_unk(jj,edge_n)
			if(patch/=-1)then
				area=triangle_area(patch,quant)
				Enormal_at_patch(patch) = Enormal_at_patch(patch)-ln/area*(-1)**(jj+1)*eigvec(edge-msh%idxs+1)*cd**2d0*mu0/junit/quant%freq/(2*pi)
			endif
			enddo
		endif
	enddo

	call MPI_ALLREDUCE(MPI_IN_PLACE,Enormal_at_patch,quant%maxpatch,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
	Enormal_at_patch = Enormal_at_patch/quant%normalize_factor

	if(ptree%MyID==Main_ID)then  ! check the normal E field on the wall
	do ii=1,quant%Nobs
	do patch =1,quant%maxpatch
	point = quant%obs_points(:,ii)
	point(3) = 0.1d0
	if(in_triangle(point,patch,quant))then
		write(666,*) point,abs(Enormal_at_patch(patch))
	endif
	enddo
	enddo
	endif


	allocate(Enormal_at_node(quant%maxnode))
	Enormal_at_node=0

    !$omp parallel do default(shared) private(patch,i,ii,node,a)
    do node=1, quant%maxnode
        ii=0; a=0.
        do patch=1, quant%maxpatch
            do i=1,3
                if (quant%node_of_patch(i,patch)==node) then
                    a=a+Enormal_at_patch(patch)
                    ii=ii+1
                endif
            enddo
        enddo
        Enormal_at_node(node)=a/dble(ii)
    enddo
    !$omp end parallel do

	if(ptree%MyID == Main_ID)then
    ! string='current'//chara//'.out'
    open(30,file=trim(adjustl(model))//'_EigEn_'//trim(adjustl(substring))//'_freq_'//trim(adjustl(substring1))//'.out')
    do node=1, quant%maxnode
        write (30,*) abs(Enormal_at_node(node))
    enddo
    close(30)
	endif
	if(ptree%MyID==Main_ID)write(*,*)'   max normal E:',maxval(abs(Enormal_at_patch))


	!!!!! Compute the maximum normal electric fields using GF, but there is a delta offset at each patch to avoid singularity
	n1 = OMP_get_wtime()

	Enormal_GF = Enormal_GF/quant%normalize_factor
	Enormal_at_patch=0
	allocate(cnt_patch(quant%maxpatch))
	cnt_patch=0
	do edge=msh%idxs,msh%idxe
		edge_n = msh%new2old(edge)
		if(edge_n<=quant%Nunk-quant%Nunk_port)then
			do jj=3,4
			patch = quant%info_unk(jj,edge_n)
			if(patch/=-1)then
				Enormal_at_patch(patch) = Enormal_at_patch(patch) + Enormal_GF(edge-msh%idxs+1)
				cnt_patch(patch) = cnt_patch(patch)+1
			endif
			enddo
		endif
	enddo

	call MPI_ALLREDUCE(MPI_IN_PLACE,Enormal_at_patch,quant%maxpatch,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
	call MPI_ALLREDUCE(MPI_IN_PLACE,cnt_patch,quant%maxpatch,MPI_INTEGER,MPI_SUM,ptree%Comm,ierr)
	do patch =1,quant%maxpatch
		if(cnt_patch(patch)>0)then
			Enormal_at_patch(patch)=Enormal_at_patch(patch)/cnt_patch(patch)
		endif
	enddo



	if(ptree%MyID==Main_ID)then  ! check the normal E field on the wall
	do ii=1,quant%Nobs
	do patch =1,quant%maxpatch
	point = quant%obs_points(:,ii)
	point(3) = 0.1d0
	if(in_triangle(point,patch,quant))then
		write(667,*) point,abs(Enormal_at_patch(patch))
	endif
	enddo
	enddo
	endif
	deallocate(cnt_patch)

	Enormal_at_node=0

    !$omp parallel do default(shared) private(patch,i,ii,node,a)
    do node=1, quant%maxnode
        ii=0; a=0.
        do patch=1, quant%maxpatch
            do i=1,3
                if (quant%node_of_patch(i,patch)==node) then
                    a=a+Enormal_at_patch(patch)
                    ii=ii+1
                endif
            enddo
        enddo
        Enormal_at_node(node)=a/dble(ii)
    enddo
    !$omp end parallel do

	if(ptree%MyID == Main_ID)then
    ! string='current'//chara//'.out'
    open(30,file=trim(adjustl(model))//'_EigEn1_'//trim(adjustl(substring))//'_freq_'//trim(adjustl(substring1))//'.out')
    do node=1, quant%maxnode
        write (30,*) abs(Enormal_at_node(node))
    enddo
    close(30)
	endif





	max_Enormal_GF = maxval(abs(Enormal_GF))
	call MPI_ALLREDUCE(MPI_IN_PLACE,max_Enormal_GF,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)


	if(ptree%MyID==Main_ID)write(*,*)'   max normal E (GF):',max_Enormal_GF
	n2 = OMP_get_wtime()
	! if(ptree%MyID==Main_ID)write(*,*)n2-n1,' seconds'



	!!!!! Compute the surface integral of |Ht|^2
	n1 = OMP_get_wtime()
	allocate(Ht_at_patch(3,quant%maxpatch*quant%integral_points))
	allocate(weight_at_patch(quant%maxpatch*quant%integral_points))
	Ht_at_patch=0
	weight_at_patch=0
	allocate (xn(quant%integral_points), yn(quant%integral_points), zn(quant%integral_points), wn(quant%integral_points))


	do edge_n=1,quant%Nunk-quant%Nunk_port
		do jj=3,4
			patch = quant%info_unk(jj,edge_n)
			area=triangle_area(patch,quant)
			call gau_grobal(edge_n,jj,xn,yn,zn,wn,quant)
			do j=1,quant%integral_points
				weight_at_patch((patch-1)*quant%integral_points+j) = wn(j)/area*2d0
			enddo
		enddo
	enddo


	do edge=msh%idxs,msh%idxe
		edge_n = msh%new2old(edge)
		if(edge_n<=quant%Nunk_int)then
			ln=sqrt((quant%xyz(1,quant%info_unk(1,edge_n))-quant%xyz(1,quant%info_unk(2,edge_n)))**2+(quant%xyz(2,quant%info_unk(1,edge_n))-quant%xyz(2,quant%info_unk(2,edge_n)))**2+(quant%xyz(3,quant%info_unk(1,edge_n))-quant%xyz(3,quant%info_unk(2,edge_n)))**2)
			do jj=3,4
				patch = quant%info_unk(jj,edge_n)
				! area=triangle_area(patch,quant)
				call gau_grobal(edge_n,jj,xn,yn,zn,wn,quant)
				do j=1,quant%integral_points
					an(1)=xn(j)-quant%xyz(1,quant%info_unk(jj+2,edge_n))
					an(2)=yn(j)-quant%xyz(2,quant%info_unk(jj+2,edge_n))
					an(3)=zn(j)-quant%xyz(3,quant%info_unk(jj+2,edge_n))
					Ht_at_patch(:,(patch-1)*quant%integral_points+j) = Ht_at_patch(:,(patch-1)*quant%integral_points+j)+ (-1)**(jj+1)*an*ln/(2)*eigvec(edge-msh%idxs+1)
				enddo
			enddo
		endif
	enddo

	call MPI_ALLREDUCE(MPI_IN_PLACE,Ht_at_patch,quant%maxpatch*quant%integral_points*3,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
    Ht_at_patch = Ht_at_patch/quant%normalize_factor

	if(ptree%MyID==Main_ID)then
		Ht2int=0
		do ii=1,quant%maxpatch*quant%integral_points
			patch = int((ii-1)/quant%integral_points)+1
			rtemp = abs(Ht_at_patch(1,ii))**2d0+abs(Ht_at_patch(2,ii))**2d0+abs(Ht_at_patch(3,ii))**2d0
			Ht2int = Ht2int + weight_at_patch(ii)* rtemp
		enddo
		Rs = 1/quant%sigma_s/sqrt(2d0/(quant%sigma_s*2*pi*quant%freq*mu0))
		Ht2int = Ht2int/2*Rs
		write(*,*)'   int Rs|H_t|^2/2 and Rs:',Ht2int,Rs
	endif

	n2 = OMP_get_wtime()
	! if(ptree%MyID==Main_ID)write(*,*)n2-n1,' seconds'



	!!!!! Compute power at ports
	if(quant%Nport>0)then
		n1 = OMP_get_wtime()
		allocate(Et_at_patch(3,quant%maxpatch*quant%integral_points))
		Et_at_patch=0
		allocate(ExH_at_ports(3,quant%Nport))
		ExH_at_ports=0
		allocate(port_of_patch(quant%maxpatch))
		port_of_patch=0
		Ht_at_patch=0

		do edge=msh%idxs,msh%idxe
			edge_n = msh%new2old(edge)
			if(edge_n>quant%Nunk-quant%Nunk_port)then
				edge_n = edge_n - quant%Nunk_port
				ln=sqrt((quant%xyz(1,quant%info_unk(1,edge_n))-quant%xyz(1,quant%info_unk(2,edge_n)))**2+(quant%xyz(2,quant%info_unk(1,edge_n))-quant%xyz(2,quant%info_unk(2,edge_n)))**2+(quant%xyz(3,quant%info_unk(1,edge_n))-quant%xyz(3,quant%info_unk(2,edge_n)))**2)
				do jj=3,4
					patch = quant%info_unk(jj,edge_n)
					call gau_grobal(edge_n,jj,xn,yn,zn,wn,quant)
					do j=1,quant%integral_points
						an(1)=xn(j)-quant%xyz(1,quant%info_unk(jj+2,edge_n))
						an(2)=yn(j)-quant%xyz(2,quant%info_unk(jj+2,edge_n))
						an(3)=zn(j)-quant%xyz(3,quant%info_unk(jj+2,edge_n))
						Ht_at_patch(:,(patch-1)*quant%integral_points+j) = Ht_at_patch(:,(patch-1)*quant%integral_points+j)+ (-1)**(jj+1)*an*ln/(2)*eigvec(edge-msh%idxs+1)
					enddo
				enddo
			endif
		enddo
		Ht_at_patch = Ht_at_patch/quant%normalize_factor

		do edge=msh%idxs,msh%idxe
			edge_n = msh%new2old(edge)
			if(edge_n>quant%Nunk-quant%Nunk_port)then
				edge_n = edge_n - quant%Nunk_port
				ln=sqrt((quant%xyz(1,quant%info_unk(1,edge_n))-quant%xyz(1,quant%info_unk(2,edge_n)))**2+(quant%xyz(2,quant%info_unk(1,edge_n))-quant%xyz(2,quant%info_unk(2,edge_n)))**2+(quant%xyz(3,quant%info_unk(1,edge_n))-quant%xyz(3,quant%info_unk(2,edge_n)))**2)
				do jj=3,4
					patch = quant%info_unk(jj,edge_n)
					area=triangle_area(patch,quant)

					cntn = quant%Nunk_int
					do ppn=1,quant%Nport
						if(edge_n<=cntn+quant%ports(ppn)%Nunk)then
							port_of_patch(patch)=ppn
							exit
						endif
						cntn = cntn + quant%ports(ppn)%Nunk
					enddo

					call gau_grobal(edge_n,jj,xn,yn,zn,wn,quant)
					do j=1,quant%integral_points
						an(1)=xn(j)-quant%xyz(1,quant%info_unk(jj+2,edge_n))
						an(2)=yn(j)-quant%xyz(2,quant%info_unk(jj+2,edge_n))
						an(3)=zn(j)-quant%xyz(3,quant%info_unk(jj+2,edge_n))
						Et_at_patch(:,(patch-1)*quant%integral_points+j) = Et_at_patch(:,(patch-1)*quant%integral_points+j)+ (-1)**(jj+1)*an*ln/(2)*eigvec(edge-msh%idxs+1)*impedence0

						! if(port_of_patch(patch)==2)then
						! 	write(*,*)port_of_patch(patch),patch, abs(Et_at_patch(:,(patch-1)*quant%integral_points+j)),abs(Ht_at_patch(:,(patch-1)*quant%integral_points+j))
						! endif

					enddo
				enddo
			endif
		enddo
		Et_at_patch = Et_at_patch/quant%normalize_factor

		call MPI_ALLREDUCE(MPI_IN_PLACE,Ht_at_patch,quant%maxpatch*quant%integral_points*3,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
		call MPI_ALLREDUCE(MPI_IN_PLACE,Et_at_patch,quant%maxpatch*quant%integral_points*3,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
		call MPI_ALLREDUCE(MPI_IN_PLACE,port_of_patch,quant%maxpatch,MPI_INTEGER,MPI_MAX,ptree%Comm,ierr)



		if(ptree%MyID==Main_ID)then
			do ii=1,quant%maxpatch*quant%integral_points
				patch = int((ii-1)/quant%integral_points)+1
				port = port_of_patch(patch)
				if(port>0)then
					call cccurl(Et_at_patch(:,ii),conjg(Ht_at_patch(:,ii)),cpoint)
					cpoint = cpoint*weight_at_patch(ii)
					point = dble(cpoint)/2d0
					ExH_at_ports(:,port) = ExH_at_ports(:,port) + point
				endif
			enddo
			do pp=1,quant%Nport
				write(*,*)'   power at port ',pp, dot_product(ExH_at_ports(:,pp),quant%ports(pp)%z)
			enddo
		endif
		deallocate(Et_at_patch)
		deallocate(port_of_patch)
		deallocate(ExH_at_ports)
		n2 = OMP_get_wtime()
		! if(ptree%MyID==Main_ID)write(*,*)n2-n1,' seconds'

	endif





	deallocate(Enormal_at_patch)
	deallocate(Enormal_at_node)

	deallocate(weight_at_patch)
	deallocate(Ht_at_patch)
	deallocate(xn,yn,zn,wn)




end subroutine EM_cavity_postprocess




!**** Compute the fields at a given point due to the K operator on a rwg basis
subroutine Field_EMSURF_K(point,field,n,quant)

    use BPACK_DEFS
    implicit none

    integer, INTENT(IN):: n
    integer flag,edge_m,edge_n,nodetemp_n,patch
    complex(kind=8) value_e,value_m,value
    integer i,j,ii,jj,iii,jjj
    real(kind=8) ln,lm,am(3),an(3),nr_m(3),nxan(3)
    real(kind=8) nr_n(3)
    complex(kind=8) ctemp,ctemp1,ctemp2,aa(3),bb(1),dg(3),dg1(3),dg2(3),dg3(3)
    complex(kind=8) imp,imp1,imp2,imp3
    real(kind=8) temp
    real(kind=8) distance
    real(kind=8) ianl,ianl1,ianl2
    real(kind=8) point(3)
    complex(kind=8) field(3)

	type(quant_EMSURF) :: quant


    real(kind=8),allocatable::xm(:),ym(:),zm(:),wm(:),xn(:),yn(:),zn(:),wn(:)

		! convert to new indices because quant%info_unk has been reordered
		edge_n = n

		allocate (xn(quant%integral_points), yn(quant%integral_points), zn(quant%integral_points), wn(quant%integral_points))

		ln=sqrt((quant%xyz(1,quant%info_unk(1,edge_n))-quant%xyz(1,quant%info_unk(2,edge_n)))**2+(quant%xyz(2,quant%info_unk(1,edge_n))-quant%xyz(2,quant%info_unk(2,edge_n)))**2+(quant%xyz(3,quant%info_unk(1,edge_n))-quant%xyz(3,quant%info_unk(2,edge_n)))**2)

		field = 0d0

		do jj=3,4
			call gau_grobal(edge_n,jj,xn,yn,zn,wn,quant)
			do j=1,quant%integral_points
				distance=sqrt((point(1)-xn(j))**2+(point(2)-yn(j))**2+(point(3)-zn(j))**2)
				an(1)=xn(j)-quant%xyz(1,quant%info_unk(jj+2,edge_n))
				an(2)=yn(j)-quant%xyz(2,quant%info_unk(jj+2,edge_n))
				an(3)=zn(j)-quant%xyz(3,quant%info_unk(jj+2,edge_n))
				dg(1)=(point(1)-xn(j))*(1+junit*quant%wavenum*distance)*exp(-junit*quant%wavenum*distance)/(4*pi*distance**3)
				dg(2)=(point(2)-yn(j))*(1+junit*quant%wavenum*distance)*exp(-junit*quant%wavenum*distance)/(4*pi*distance**3)
				dg(3)=(point(3)-zn(j))*(1+junit*quant%wavenum*distance)*exp(-junit*quant%wavenum*distance)/(4*pi*distance**3)
				call ccurl(an,dg,dg2)
				field = field - (-1)**(jj+1)*dg2*wn(j)*ln
			enddo
		enddo
		deallocate(xn,yn,zn,wn)

    return

end subroutine Field_EMSURF_K



!**** Compute the fields at a given point due to the T operator on a rwg basis
subroutine Field_EMSURF_T(point,field,n,quant)

    use BPACK_DEFS
    implicit none

    integer, INTENT(IN):: n
    integer flag,edge_m,edge_n,nodetemp_n,patch
    integer i,j,ii,jj,iii,jjj
    real(kind=8) ln,an(3)
    complex(kind=8) ctemp,ctemp1,ctemp2,ctemp3
    real(kind=8) temp
	real(kind=8) distance,distance1,distance2,distance3
	complex(kind=8) phi,phi1,phi2,phi3,phase
    real(kind=8) delta
    real(kind=8) point(3)
    complex(kind=8) field(3),field_A(3),dg(3),field_phi(3)

	type(quant_EMSURF) :: quant


    real(kind=8),allocatable::xn(:),yn(:),zn(:),wn(:)



		! convert to new indices because quant%info_unk has been reordered
		edge_n = n
		allocate (xn(quant%integral_points), yn(quant%integral_points), zn(quant%integral_points), wn(quant%integral_points))

		ln=sqrt((quant%xyz(1,quant%info_unk(1,edge_n))-quant%xyz(1,quant%info_unk(2,edge_n)))**2+(quant%xyz(2,quant%info_unk(1,edge_n))-quant%xyz(2,quant%info_unk(2,edge_n)))**2+(quant%xyz(3,quant%info_unk(1,edge_n))-quant%xyz(3,quant%info_unk(2,edge_n)))**2)

		field_A = 0d0
		field_phi = 0d0
		do jj=3,4
			call gau_grobal(edge_n,jj,xn,yn,zn,wn,quant)
			do j=1,quant%integral_points
				an(1)=xn(j)-quant%xyz(1,quant%info_unk(jj+2,edge_n))
				an(2)=yn(j)-quant%xyz(2,quant%info_unk(jj+2,edge_n))
				an(3)=zn(j)-quant%xyz(3,quant%info_unk(jj+2,edge_n))
				an = an*(-1)**(jj+1)
				distance=sqrt((point(1)-xn(j))**2+(point(2)-yn(j))**2+(point(3)-zn(j))**2)

				phase = -quant%wavenum*distance
				ctemp=exp(-aimag(phase))*(cos(dble(phase))+junit*sin(dble(phase)))/(4*pi)
				field_A = field_A + an*wn(j)*ctemp/distance


				dg(1)=(point(1)-xn(j))*(1+junit*quant%wavenum*distance)*ctemp/(distance**3)
				dg(2)=(point(2)-yn(j))*(1+junit*quant%wavenum*distance)*ctemp/(distance**3)
				dg(3)=(point(3)-zn(j))*(1+junit*quant%wavenum*distance)*ctemp/(distance**3)

				field_phi = field_phi + 2d0*wn(j)*dg*(-1)**(jj+1)

			enddo
		enddo

		field_A = field_A * junit*quant%wavenum
		field_phi = field_phi * junit/quant%wavenum

		! write(777,*)'x', field_A(1), field_phi(1)
		! write(777,*)'y', field_A(2), field_phi(2)
		! write(777,*)'z', field_A(3), field_phi(3)

		field = field_A -  field_phi
		field = field  * ln

		deallocate(xn,yn,zn,wn)



    return

end subroutine Field_EMSURF_T




!**** Compute the fields at a given point due to the T or K operator on a rwg basis
subroutine Field_EMSURF(point,field,n,quant)

    use BPACK_DEFS
    implicit none

	integer, INTENT(IN):: n
	integer edge
	type(quant_EMSURF) :: quant
	complex(kind=8) field(3)
	real(kind=8) point(3)
	edge=n
	if(edge<=quant%Nunk-quant%Nunk_port)then
		call Field_EMSURF_T(point,field,edge,quant)
	else
		edge = edge-quant%Nunk_port
		call Field_EMSURF_K(point,field,edge,quant)
	endif
	! field=field*2d0
end subroutine Field_EMSURF

end module EMSURF_PORT_MODULE