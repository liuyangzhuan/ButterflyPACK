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
!> @file
!> @brief This files contains functions and data types for the 3D EFIE/MFIE/CIFE examples
!> @details Note that the use of the following \n
!> #define DAT 0 \n
!> #include "zButterflyPACK_config.fi" \n
!> which will macro replace subroutine, function, type names with those defined in SRC_DOUBLECOMLEX with double-complex precision

! This exmple works with double-complex precision data
module EMSURF_MODULE
use z_BPACK_DEFS
use z_MISC_Utilities
implicit none

	!**** define your application-related variables here


	!**** edges of each node
	type edge_node
		integer Nedge
		integer,allocatable::edges(:,:)
	end type edge_node

	!**** quantities related to geometries, meshes, unknowns and points
	type quant_EMSURF
		real(kind=8) wavenum    ! CEM: wave number
		real(kind=8) wavelength  ! CEM: wave length
		real(kind=8) freq       ! CEM: frequency
		real(kind=8) rank_approximate_para1, rank_approximate_para2, rank_approximate_para3 ! CEM: rank estimation parameter
		integer RCS_static  ! CEM: 1: monostatic or 2: bistatic RCS
		integer RCS_Nsample ! CEM: number of RCS samples
		real(kind=8):: CFIE_alpha ! CEM: combination parameter in CFIE

		integer Nunk ! size of the matrix
		real(kind=8),allocatable:: xyz(:,:)   ! coordinates of the points
		integer,allocatable:: info_unk(:,:)
		! for 3D z_mesh: 0 point to coordinates of each edge center (unknown x), 1-2 point to coordinates of each edge vertice, 3-4 point to two patches that share the same edge, 5-6 point to coordinates of last vertice of the two patches

		! 3D z_mesh
		integer maxnode ! # of vertices in a z_mesh
		integer maxedge ! # of edges in a z_mesh
		real(kind=8) maxedgelength,minedgelength ! maximum and minimum edge length for 2D and 3D meshes
		integer integral_points ! #of Gauss quadrature points on a triangular
		integer maxpatch ! # of triangular patches
		integer mesh_normal	 ! flags to flip the unit normal vectors of a triangular patch
		real(kind=8) scaling  ! scaling factor of the coordinates of vertices of a 3D z_mesh
		real(kind=8), allocatable :: ng1(:),ng2(:),ng3(:),gauss_w(:) ! Gass quadrature and weights
		real(kind=8),allocatable:: normal_of_patch(:,:) ! normal vector of each triangular patch
		integer,allocatable:: node_of_patch(:,:) ! vertices of each triangular patch
		integer,allocatable:: edge_of_patch(:,:,:) ! unknown ids of each triangular patch
		type(edge_node),allocatable:: edge_of_node(:) ! edges of each vertice
		CHARACTER (LEN=1000) DATA_DIR
		integer::CMmode=0 !  1: solve the characteristic mode, 0: solve the eigen mode
		integer::SI=0 ! 0: regular mode 1: shift-invert mode
		complex(kind=8)::shift=0d0 ! the shift value in shift-invert Arnoldi iterations
		integer::nev=1 ! nubmer of requested eigen values
		character(len=2) which ! which portion of eigen spectrum
		real(kind=8) tol_eig ! tolerance in arpack

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
subroutine Zelem_EMSURF(m,n,value,quant)

    use z_BPACK_DEFS
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

	class(*),pointer :: quant


    real(kind=8),allocatable::xm(:),ym(:),zm(:),wm(:),xn(:),yn(:),zn(:),wn(:)


	select TYPE(quant)
	type is (quant_EMSURF)

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
						imp=(0.,0.);imp1=(0.,0.);imp2=(0.,0.);imp3=(0.,0.)
						an(1)=xm(i)-quant%xyz(1,quant%info_unk(jj+2,edge_n))
						an(2)=ym(i)-quant%xyz(2,quant%info_unk(jj+2,edge_n))
						an(3)=zm(i)-quant%xyz(3,quant%info_unk(jj+2,edge_n))
						call z_scalar(am,an,temp)
						value_m=value_m+(-1)**(ii+1)*(-1)**(jj+1)*0.5*temp/(2.*area)*wm(i)
						do j=1,quant%integral_points
							distance=sqrt((xm(i)-xn(j))**2+(ym(i)-yn(j))**2+(zm(i)-zn(j))**2)
							if(distance==0)then
								imp=imp+wn(j)*(-BPACK_junit*quant%wavenum)
								imp1=imp1+quant%ng1(j)*wn(j)*(-BPACK_junit*quant%wavenum)
								imp2=imp2+quant%ng2(j)*wn(j)*(-BPACK_junit*quant%wavenum)
								ianl=ianalytic(edge_n,jj,xn(j),yn(j),zn(j),quant)
								ianl1=ianalytic2(edge_n,jj,xn(j),yn(j),zn(j),1,quant)
								ianl2=ianalytic2(edge_n,jj,xn(j),yn(j),zn(j),2,quant)
								imp=imp+ianl  !debug
								imp1=imp1+ianl1
								imp2=imp2+ianl2
							else
								imp=imp+wn(j)*exp(-BPACK_junit*quant%wavenum*distance)/distance
								imp1=imp1+quant%ng1(j)*wn(j)*exp(-BPACK_junit*quant%wavenum*distance)/distance
								imp2=imp2+quant%ng2(j)*wn(j)*exp(-BPACK_junit*quant%wavenum*distance)/distance
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
							an(1)=xn(j)-quant%xyz(1,quant%info_unk(jj+2,edge_n))
							an(2)=yn(j)-quant%xyz(2,quant%info_unk(jj+2,edge_n))
							an(3)=zn(j)-quant%xyz(3,quant%info_unk(jj+2,edge_n))
							dg(1)=(xm(i)-xn(j))*(1+BPACK_junit*quant%wavenum*distance)*exp(-BPACK_junit*quant%wavenum*distance)/(4*BPACK_pi*distance**3)
							dg(2)=(ym(i)-yn(j))*(1+BPACK_junit*quant%wavenum*distance)*exp(-BPACK_junit*quant%wavenum*distance)/(4*BPACK_pi*distance**3)
							dg(3)=(zm(i)-zn(j))*(1+BPACK_junit*quant%wavenum*distance)*exp(-BPACK_junit*quant%wavenum*distance)/(4*BPACK_pi*distance**3)
								call z_rccurl(an,dg,dg1)
								call z_rccurl(nr_m,dg1,dg2)
								call z_cscalar(dg2,am,ctemp)
								value_m=value_m-(-1)**(ii+1)*(-1)**(jj+1)*ctemp*wm(i)*wn(j)
							imp=imp+wn(j)*exp(-BPACK_junit*quant%wavenum*distance)/distance
							imp1=imp1+quant%ng1(j)*wn(j)*exp(-BPACK_junit*quant%wavenum*distance)/distance
							imp2=imp2+quant%ng2(j)*wn(j)*exp(-BPACK_junit*quant%wavenum*distance)/distance
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
				call z_cscalar(aa,am,ctemp)
				ctemp1=ctemp1+(-1)**(ii+1)*ctemp*wm(i)
				ctemp2=ctemp2+4.*(-1)**(ii+1)*bb(1)*wm(i)
			enddo
		enddo
		value_e=ln*lm*BPACK_junit*(ctemp1-ctemp2)/8./BPACK_pi**2d0/quant%freq/BPACK_eps0
		value_m=value_m*lm*ln

		value=quant%CFIE_alpha*value_e+(1.-quant%CFIE_alpha)*BPACK_impedence0*value_m

		deallocate(xm,ym,zm,wm,xn,yn,zn,wn)

	class default
		write(*,*)"unexpected type"
		stop
	end select

    return

end subroutine Zelem_EMSURF




subroutine Zelem_EMSURF_block(Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps, quant) ! interface of user-defined element extraction routine in Fortran. allrows,allcols represents indices in natural order

    use z_BPACK_DEFS
    implicit none

    integer i,j,ii,jj

	class(*), pointer :: quant
	integer:: Ninter
	integer:: allrows(:), allcols(:)
	complex(kind=8),target::alldat_loc(:)
	integer::colidx(Ninter), rowidx(Ninter), pgidx(Ninter)
	integer::Npmap, pmaps(Npmap, 3)
	integer::myid,nmpi,ierr,nn,pp,nprow,npcol,pid,nr,nc,nn1
	integer::idxcnt_m,idxcnt_n,idx_m_s,idx_m_e,idx_n_s,idx_n_e,iii,jjj
	integer(kind=8)::idx_row=0,idx_col=0,idx_val=0,idx_temp
	integer::NinterNew=0,nrmax=0,ncmax=0,nvalmax=0,idxr,idxc,m_1,n_1,m,n
	integer(kind=8)::idx_row_map(Ninter)
	integer(kind=8)::idx_col_map(Ninter)
	integer(kind=8)::idx_val_map(Ninter)
	integer::inter_map(Ninter)
	complex(kind=8) value
	integer,allocatable::M_elemlist(:,:),N_elemlist(:,:),M_element_offset(:),N_element_offset(:)
	idx_row_map=0
	idx_col_map=0
	idx_val_map=0
	inter_map=0


	call MPI_comm_rank(MPI_COMM_WORLD,myid,ierr)
	call MPI_comm_size(MPI_COMM_WORLD,nmpi,ierr)

	select TYPE(quant)
	type is (quant_EMSURF)
		NinterNew=0
		idx_row=0
		idx_col=0
		idx_val=0
		nrmax=0
		ncmax=0
		nvalmax=0
		do nn=1,Ninter
			pp=pgidx(nn)
			nprow = pmaps(pp,1)
			npcol = pmaps(pp,2)
			pid = pmaps(pp,3)
			nr = rowidx(nn)
			nc = colidx(nn)
			if(nprow*npcol==1)then
				if(myid==pid)then
					! write(*,*)Ninter,NinterNew
					idx_row_map(NinterNew+1)=idx_row
					idx_col_map(NinterNew+1)=idx_col
					idx_val_map(NinterNew+1)=idx_val
					idx_val=idx_val+nr*nc
					inter_map(NinterNew+1)=nn
					NinterNew=NinterNew+1
				endif
				idx_row = idx_row+nr
				idx_col = idx_col+nc
				nrmax = max(nr,nrmax)
				ncmax = max(nc,ncmax)
				nvalmax = max(nr*nc,nvalmax)
			else
				write(*,*)"nprow*npcol>1 in Zelem_EMSURF_block"
				stop
			endif
		enddo
		if(idx_val>0)then
			alldat_loc(1:idx_val)=0
		endif

		!$omp parallel do default(shared) private(nn1)
		do nn1=1,NinterNew
			call Zelem_EMSURF_oneblock(nn1,inter_map,allrows, allcols, alldat_loc, rowidx, colidx, idx_row_map,idx_col_map,idx_val_map, quant)

		enddo
		!$omp end parallel do
	class default
		write(*,*)"unexpected type"
		stop
	end select

    return

end subroutine Zelem_EMSURF_block



subroutine Zelem_EMSURF_oneblock(nn1,inter_map,allrows, allcols, alldat_loc, rowidx, colidx, idx_row_map,idx_col_map,idx_val_map,quant)

    use z_BPACK_DEFS
    implicit none

    integer flag,edge_m,edge_n,nodetemp_n,patch
    complex(kind=8) value_e,value_m,value
    integer i,j,ii,jj,iii,jjj
    real(kind=8) ln,lm,am(3),an(3),nr_m(3)
    real(kind=8) nr_n(3)
    complex(kind=8) ctemp,ctemp1,ctemp2,aa(3),bb(1),dg(3),dg1(3),dg2(3)
    complex(kind=8) imp,imp1,imp2,imp3, consts
    real(kind=8) temp
    real(kind=8) distance
    real(kind=8) ianl,ianl1,ianl2
    real(kind=8) area

	type(quant_EMSURF):: quant
    real(kind=8),allocatable::xm(:),ym(:),zm(:),wm(:),xn(:),yn(:),zn(:),wn(:)
    integer kk,ss,sss,ttt
	integer:: Tord,Sord,N_Tnodes,N_Snodes
	integer M_elem,N_elem,Selem_id,Telem_id,test_unknown_id,source_unknown_id, M_elem_tot, N_elem_tot, lastelem
	real(kind=8):: row_sign,col_sign
	integer,parameter::Nnodes_max=3
	complex(kind=8),dimension(1:Nnodes_max,1:Nnodes_max,1:3):: temp_matrix
	complex(kind=8),dimension(1:Nnodes_max,1:Nnodes_max):: final_matrix


	integer:: Ninter
	integer:: allrows(:), allcols(:)
	complex(kind=8),target::alldat_loc(:)
	integer::colidx(:), rowidx(:)
	integer::myid,nmpi,ierr,nn,pp,nprow,npcol,pid,nr,nc,nn1
	integer::idxcnt_m,idxcnt_n,idx_m_s,idx_m_e,idx_n_s,idx_n_e
	integer(kind=8)::idx_row=0,idx_col=0,idx_val=0,idx_temp
	integer::NinterNew=0,nrmax=0,ncmax=0,nvalmax=0,idxr,idxc,m_1,n_1,m,n
	integer(kind=8)::idx_row_map(:)
	integer(kind=8)::idx_col_map(:)
	integer(kind=8)::idx_val_map(:)
	integer::inter_map(:)
	integer,allocatable::M_elemlist(:,:),N_elemlist(:,:),M_element_offset(:),N_element_offset(:)
	integer::trange,srange


	nn=inter_map(nn1)
	nr = rowidx(nn)
	nc = colidx(nn)


	!!!! build the mapping from all active elements to the unknowns
	allocate(M_elemlist(nr*2,3))
	M_elemlist=-1
	idxcnt_m=0
	do idxr=1,nr
		m = allrows(idx_row_map(nn1)+idxr)
		edge_m = m
		do ii=3,4
			Telem_id = quant%info_unk(ii,edge_m)
			if(Telem_id/=-1)then
				idxcnt_m=idxcnt_m+1
				M_elemlist(idxcnt_m,1)=Telem_id
				M_elemlist(idxcnt_m,3)=idxr
				do iii=1,3
					if(quant%edge_of_patch(Telem_id,iii,1)==edge_m)then
						M_elemlist(idxcnt_m,2)=iii
					endif
				enddo
			endif
		enddo
	enddo
	call z_PIKSRT_INT_Multi(idxcnt_m, 3, M_elemlist(1:idxcnt_m,:))
	M_elem_tot=0
	allocate(M_element_offset(nr*2))
	M_element_offset=0
	lastelem=-1
	do ii=1,idxcnt_m
		if(M_elemlist(ii,1)/=lastelem)then
			M_elem_tot=M_elem_tot+1
			M_element_offset(M_elem_tot)=ii
		endif
		lastelem = M_elemlist(ii,1)
	enddo


	allocate(N_elemlist(nc*2,3))
	N_elemlist=-1
	idxcnt_n=0
	do idxc=1,nc
		n = allcols(idx_col_map(nn1)+idxc)
		edge_n = n
		do ii=3,4
			Selem_id = quant%info_unk(ii,edge_n)
			if(Selem_id/=-1)then
				idxcnt_n=idxcnt_n+1
				N_elemlist(idxcnt_n,1)=Selem_id
				N_elemlist(idxcnt_n,3)=idxc
				do iii=1,3
					if(quant%edge_of_patch(Selem_id,iii,1)==edge_n)then
						N_elemlist(idxcnt_n,2)=iii
					endif
				enddo
			endif
		enddo
	enddo
	call z_PIKSRT_INT_Multi(idxcnt_n, 3, N_elemlist(1:idxcnt_n,:))
	N_elem_tot=0
	allocate(N_element_offset(nc*2))
	N_element_offset=0
	lastelem=-1
	do ii=1,idxcnt_n
		if(N_elemlist(ii,1)/=lastelem)then
			N_elem_tot=N_elem_tot+1
			N_element_offset(N_elem_tot)=ii
		endif
		lastelem = N_elemlist(ii,1)
	enddo



	allocate (xm(quant%integral_points), ym(quant%integral_points), zm(quant%integral_points), wm(quant%integral_points))
	allocate (xn(quant%integral_points), yn(quant%integral_points), zn(quant%integral_points), wn(quant%integral_points))



	!!!! loop all element pairs
	do kk=1,M_elem_tot
		idx_m_s = M_element_offset(kk)
		if(kk==M_elem_tot)then
			idx_m_e = idxcnt_m
		else
			idx_m_e = M_element_offset(kk+1)-1
		endif
		Telem_id = M_elemlist(idx_m_s,1)
		call gau_grobal_patchid(Telem_id,xm,ym,zm,wm,quant)
		nr_m(1:3)=quant%normal_of_patch(1:3,Telem_id)

		do ss=1,N_elem_tot
			idx_n_s = N_element_offset(ss)
			if(ss==N_elem_tot)then
				idx_n_e = idxcnt_n
			else
				idx_n_e = N_element_offset(ss+1)-1
			endif
			Selem_id = N_elemlist(idx_n_s,1)

			call gau_grobal_patchid(Selem_id,xn,yn,zn,wn,quant)
			nr_n(1:3)=quant%normal_of_patch(1:3,Selem_id)

			trange=0
			do ttt=1,3
				if(quant%edge_of_patch(Telem_id,ttt,1)/=-1)trange=trange+1
			enddo
			srange=0
			do sss=1,3
				if(quant%edge_of_patch(Selem_id,sss,1)/=-1)srange=srange+1
			enddo

			temp_matrix = 0
			final_matrix = 0
			if (Telem_id==Selem_id) then
				area=triangle_area(Telem_id,quant)

				do ttt=1,trange
					edge_m = quant%edge_of_patch(Telem_id,ttt,1)
					ii = quant%edge_of_patch(Telem_id,ttt,2)
					do sss=1,srange
						edge_n = quant%edge_of_patch(Selem_id,sss,1)
						jj = quant%edge_of_patch(Selem_id,sss,2)
						ctemp1=(0.,0.)
						ctemp2=(0.,0.)
						value_m=(0.,0.)
						do i=1,quant%integral_points
							am(1)=xm(i)-quant%xyz(1,quant%info_unk(ii+2,edge_m))
							am(2)=ym(i)-quant%xyz(2,quant%info_unk(ii+2,edge_m))
							am(3)=zm(i)-quant%xyz(3,quant%info_unk(ii+2,edge_m))
							bb(1)=(0.,0.)
							aa(1:3)=(0.,0.)
							imp=(0.,0.);imp1=(0.,0.);imp2=(0.,0.);imp3=(0.,0.)
							an(1)=xm(i)-quant%xyz(1,quant%info_unk(jj+2,edge_n))
							an(2)=ym(i)-quant%xyz(2,quant%info_unk(jj+2,edge_n))
							an(3)=zm(i)-quant%xyz(3,quant%info_unk(jj+2,edge_n))
							call z_scalar(am,an,temp)
							temp_matrix(ttt,sss,3)=temp_matrix(ttt,sss,3)+(-1)**(ii+1)*(-1)**(jj+1)*0.5*temp/(2.*area)*wm(i)
							do j=1,quant%integral_points
								distance=sqrt((xm(i)-xn(j))**2+(ym(i)-yn(j))**2+(zm(i)-zn(j))**2)
								if(distance==0)then
									imp=imp+wn(j)*(-BPACK_junit*quant%wavenum)
									imp1=imp1+quant%ng1(j)*wn(j)*(-BPACK_junit*quant%wavenum)
									imp2=imp2+quant%ng2(j)*wn(j)*(-BPACK_junit*quant%wavenum)
									ianl=ianalytic(edge_n,jj,xn(j),yn(j),zn(j),quant)
									ianl1=ianalytic2(edge_n,jj,xn(j),yn(j),zn(j),1,quant)
									ianl2=ianalytic2(edge_n,jj,xn(j),yn(j),zn(j),2,quant)
									imp=imp+ianl  !debug
									imp1=imp1+ianl1
									imp2=imp2+ianl2
								else
									imp=imp+wn(j)*exp(-BPACK_junit*quant%wavenum*distance)/distance
									imp1=imp1+quant%ng1(j)*wn(j)*exp(-BPACK_junit*quant%wavenum*distance)/distance
									imp2=imp2+quant%ng2(j)*wn(j)*exp(-BPACK_junit*quant%wavenum*distance)/distance
								endif
							enddo
							imp3=imp-imp1-imp2
							nodetemp_n=quant%info_unk(jj+2,edge_n)
							patch=Selem_id
							do jjj=1,3
								aa(jjj)=aa(jjj)+(-1)**(jj+1)*quant%wavenum**2*(quant%xyz(jjj,quant%node_of_patch(1,patch))*imp1+quant%xyz(jjj,quant%node_of_patch(2,patch))*imp2+quant%xyz(jjj,quant%node_of_patch(3,patch))*imp3-quant%xyz(jjj,nodetemp_n)*imp)
							enddo
							bb(1)=bb(1)+(-1)**(jj+1)*imp
							call z_cscalar(aa,am,ctemp)
							temp_matrix(ttt,sss,1)=temp_matrix(ttt,sss,1)+(-1)**(ii+1)*ctemp*wm(i)
							temp_matrix(ttt,sss,2)=temp_matrix(ttt,sss,2)+4.*(-1)**(ii+1)*bb(1)*wm(i)
						enddo

						lm=sqrt((quant%xyz(1,quant%info_unk(1,edge_m))-quant%xyz(1,quant%info_unk(2,edge_m)))**2+(quant%xyz(2,quant%info_unk(1,edge_m))-quant%xyz(2,quant%info_unk(2,edge_m)))**2+(quant%xyz(3,quant%info_unk(1,edge_m))-quant%xyz(3,quant%info_unk(2,edge_m)))**2)
						ln=sqrt((quant%xyz(1,quant%info_unk(1,edge_n))-quant%xyz(1,quant%info_unk(2,edge_n)))**2+(quant%xyz(2,quant%info_unk(1,edge_n))-quant%xyz(2,quant%info_unk(2,edge_n)))**2+(quant%xyz(3,quant%info_unk(1,edge_n))-quant%xyz(3,quant%info_unk(2,edge_n)))**2)
						value_e=ln*lm*BPACK_junit*(temp_matrix(ttt,sss,1)-temp_matrix(ttt,sss,2))/8./BPACK_pi**2d0/quant%freq/BPACK_eps0
						value_m=temp_matrix(ttt,sss,3)*lm*ln
						final_matrix(ttt,sss)=quant%CFIE_alpha*value_e+(1.-quant%CFIE_alpha)*BPACK_impedence0*value_m
					enddo
				enddo
			else

				do i=1,quant%integral_points
					imp=(0.,0.);imp1=(0.,0.);imp2=(0.,0.);imp3=(0.,0.)
					do j=1,quant%integral_points
						distance=sqrt((xm(i)-xn(j))**2+(ym(i)-yn(j))**2+(zm(i)-zn(j))**2)
						consts = exp(-BPACK_junit*quant%wavenum*distance)

						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						! The following block is only used for MFIE
						if(abs((1.-quant%CFIE_alpha))>1d-14)then
							dg(1)=(xm(i)-xn(j))*(1+BPACK_junit*quant%wavenum*distance)*consts/(4*BPACK_pi*distance**3)
							dg(2)=(ym(i)-yn(j))*(1+BPACK_junit*quant%wavenum*distance)*consts/(4*BPACK_pi*distance**3)
							dg(3)=(zm(i)-zn(j))*(1+BPACK_junit*quant%wavenum*distance)*consts/(4*BPACK_pi*distance**3)

							do ttt=1,trange
								edge_m = quant%edge_of_patch(Telem_id,ttt,1)
								ii = quant%edge_of_patch(Telem_id,ttt,2)
								am(1)=xm(i)-quant%xyz(1,quant%info_unk(ii+2,edge_m))
								am(2)=ym(i)-quant%xyz(2,quant%info_unk(ii+2,edge_m))
								am(3)=zm(i)-quant%xyz(3,quant%info_unk(ii+2,edge_m))
								do sss=1,srange
									edge_n = quant%edge_of_patch(Selem_id,sss,1)
									jj = quant%edge_of_patch(Selem_id,sss,2)

									an(1)=xn(j)-quant%xyz(1,quant%info_unk(jj+2,edge_n))
									an(2)=yn(j)-quant%xyz(2,quant%info_unk(jj+2,edge_n))
									an(3)=zn(j)-quant%xyz(3,quant%info_unk(jj+2,edge_n))

									call z_rccurl(an,dg,dg1)
									call z_rccurl(nr_m,dg1,dg2)
									call z_cscalar(dg2,am,ctemp)
									temp_matrix(ttt,sss,3)=temp_matrix(ttt,sss,3)-(-1)**(ii+1)*(-1)**(jj+1)*ctemp*wm(i)*wn(j)

								enddo
							enddo
						endif
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						! The following block is only used for EFIE
						if(abs(quant%CFIE_alpha)>1d-14)then
							imp=imp+wn(j)*consts/distance
							imp1=imp1+quant%ng1(j)*wn(j)*consts/distance
							imp2=imp2+quant%ng2(j)*wn(j)*consts/distance
						endif
					enddo

					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					! The following block is only used for EFIE
					if(abs(quant%CFIE_alpha)>1d-14)then
						imp3=imp-imp1-imp2
						do ttt=1,trange
							edge_m = quant%edge_of_patch(Telem_id,ttt,1)
							ii = quant%edge_of_patch(Telem_id,ttt,2)
							am(1)=xm(i)-quant%xyz(1,quant%info_unk(ii+2,edge_m))
							am(2)=ym(i)-quant%xyz(2,quant%info_unk(ii+2,edge_m))
							am(3)=zm(i)-quant%xyz(3,quant%info_unk(ii+2,edge_m))
							do sss=1,srange
								edge_n = quant%edge_of_patch(Selem_id,sss,1)
								jj = quant%edge_of_patch(Selem_id,sss,2)
								nodetemp_n=quant%info_unk(jj+2,edge_n)
								patch=Selem_id
								do jjj=1,3
								aa(jjj)=(-1)**(jj+1)*quant%wavenum**2*(quant%xyz(jjj,quant%node_of_patch(1,patch))*imp1+quant%xyz(jjj,quant%node_of_patch(2,patch))*imp2+quant%xyz(jjj,quant%node_of_patch(3,patch))*imp3-quant%xyz(jjj,nodetemp_n)*imp)
								enddo
								bb(1)=(-1)**(jj+1)*imp
								call z_cscalar(aa,am,ctemp)
								temp_matrix(ttt,sss,1)=temp_matrix(ttt,sss,1)+(-1)**(ii+1)*ctemp*wm(i)
								temp_matrix(ttt,sss,2)=temp_matrix(ttt,sss,2)+4.*(-1)**(ii+1)*bb(1)*wm(i)
							enddo
						enddo
					endif
				enddo

				do ttt=1,trange
					edge_m = quant%edge_of_patch(Telem_id,ttt,1)
					ii = quant%edge_of_patch(Telem_id,ttt,2)
					do sss=1,srange
						edge_n = quant%edge_of_patch(Selem_id,sss,1)
						jj = quant%edge_of_patch(Selem_id,sss,2)

						lm=sqrt((quant%xyz(1,quant%info_unk(1,edge_m))-quant%xyz(1,quant%info_unk(2,edge_m)))**2+(quant%xyz(2,quant%info_unk(1,edge_m))-quant%xyz(2,quant%info_unk(2,edge_m)))**2+(quant%xyz(3,quant%info_unk(1,edge_m))-quant%xyz(3,quant%info_unk(2,edge_m)))**2)
						ln=sqrt((quant%xyz(1,quant%info_unk(1,edge_n))-quant%xyz(1,quant%info_unk(2,edge_n)))**2+(quant%xyz(2,quant%info_unk(1,edge_n))-quant%xyz(2,quant%info_unk(2,edge_n)))**2+(quant%xyz(3,quant%info_unk(1,edge_n))-quant%xyz(3,quant%info_unk(2,edge_n)))**2)

						value_e=ln*lm*BPACK_junit*(temp_matrix(ttt,sss,1)-temp_matrix(ttt,sss,2))/8./BPACK_pi**2d0/quant%freq/BPACK_eps0
						value_m=temp_matrix(ttt,sss,3)*lm*ln
						final_matrix(ttt,sss)=quant%CFIE_alpha*value_e+(1.-quant%CFIE_alpha)*BPACK_impedence0*value_m
					enddo
				enddo
			endif

			do iii=idx_m_s,idx_m_e
				ttt=M_elemlist(iii,2)
				idxr= M_elemlist(iii,3)
				do jjj=idx_n_s,idx_n_e
					sss=N_elemlist(jjj,2)
					idxc= N_elemlist(jjj,3)
					alldat_loc(idx_val_map(nn1)+idxr+(idxc-1)*nr) = alldat_loc(idx_val_map(nn1)+idxr+(idxc-1)*nr) &
					&				+ final_matrix(ttt,sss)
				enddo
			enddo
		enddo
	enddo
	deallocate(M_elemlist)
	deallocate(N_elemlist)
	deallocate(M_element_offset)
	deallocate(N_element_offset)

	deallocate(xm,ym,zm,wm,xn,yn,zn,wn)


end subroutine Zelem_EMSURF_oneblock




	!**** user-defined subroutine to sample real(Z_mn)
	subroutine Zelem_EMSURF_Real(m,n,value_e,quant)
		use z_BPACK_DEFS
		implicit none
		integer, INTENT(IN):: m,n
		complex(kind=8) value_e
		class(*),pointer :: quant

		call Zelem_EMSURF(m,n,value_e,quant)
		value_e=dble(value_e)
	end subroutine  Zelem_EMSURF_Real


	!**** user-defined subroutine to sample Z_mn-sigma*Delta_mn or Z_mn-sigma*real(Z_mn)
	subroutine Zelem_EMSURF_Shifted(m,n,value_e,quant)

		use z_BPACK_DEFS
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

  use z_BPACK_DEFS
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



  subroutine gau_grobal_patchid(ii,x,y,z,w,quant)

  use z_BPACK_DEFS
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
  end subroutine gau_grobal_patchid




  subroutine gauss_points(quant)

      use z_BPACK_DEFS
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

	v1=1./3.
	v2=(6.-sqrt(15.))/21.
	v3=(9.+2*sqrt(15.))/21.
	v4=(6.+sqrt(15.))/21.
	v5=(9.-2*sqrt(15.))/21.

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

use     z_BPACK_DEFS
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
 call z_rrcurl(a,b,w)
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
 call z_rrcurl(w,u,v)
 call z_scalar(u,b,u3)
v3=2.*area/l(3)

b(1)=xi-quant%xyz(1,node1)
b(2)=yi-quant%xyz(2,node1)
b(3)=zi-quant%xyz(3,node1)
 call z_scalar(u,b,u0)
 call z_scalar(v,b,v0)
 call z_scalar(w,b,w0)

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

use z_BPACK_DEFS
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
call z_rrcurl(a,b,w)
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
 call z_rrcurl(w,u,v)
 call z_scalar(u,b,u3)
v3=2.*area/l(3)

b(1)=xi-quant%xyz(1,node1)
b(2)=yi-quant%xyz(2,node1)
b(3)=zi-quant%xyz(3,node1)
 call z_scalar(u,b,u0)
 call z_scalar(v,b,v0)
 call z_scalar(w,b,w0)


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
call z_rrcurl(s1,w,m1)
call z_rrcurl(s2,w,m2)
call z_rrcurl(s3,w,m3)
call z_scalar(u,m1,temp1)
temp1=temp1*f3(1)/2
call z_scalar(u,m2,temp2)
temp2=temp2*f3(2)/2
call z_scalar(u,m3,temp3)
temp3=temp3*f3(3)/2
iua=temp1+temp2+temp3

call z_scalar(v,m1,temp1)
temp1=temp1*f3(1)/2
call z_scalar(v,m2,temp2)
temp2=temp2*f3(2)/2
call z_scalar(v,m3,temp3)
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


subroutine current_node_patch_mapping(string,curr,msh,quant,ptree)

    use z_BPACK_DEFS
    implicit none

    integer patch, edge, node_patch(3), node_edge, node, info_idx
    integer i,j,k,ii,jj,kk,flag,lr
    real(kind=8) center(3), r, a, signs, current_patch(3),  current_abs
	character(*)::string
    complex(kind=8)::curr(:)
	type(quant_EMSURF)::quant
	type(z_mesh)::msh
	type(z_proctree)::ptree

    real(kind=8),allocatable :: current_at_patch(:),vec_current_at_patch(:,:), current_at_node(:)
	complex(kind=8),allocatable:: current1(:,:)
	integer ierr

	allocate (current_at_node(quant%maxnode),current_at_patch(quant%maxpatch))
	allocate (current1(quant%maxnode,1))
	allocate (vec_current_at_patch(3,quant%maxpatch))
	current_at_node=0
	current_at_patch=0
	vec_current_at_patch=0
#ifdef HAVE_OPENMP
	!$omp parallel do default(shared) private(lr,patch,signs,info_idx,i,j,center,edge,current_abs,current_patch,a,r)
#endif
	do edge=msh%idxs,msh%idxe
		do lr=3,4
		patch = quant%info_unk(lr,msh%new2old(edge))
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
				a=quant%xyz(j,quant%info_unk(info_idx,msh%new2old(edge)))-center(j)
				r=r+a**2
				current_patch(j)=signs*current_abs*a
			enddo
			r=sqrt(r)
			current_patch = current_patch/r
			do i=1,3
#ifdef HAVE_TASKLOOP
			!$omp atomic
#endif
			vec_current_at_patch(i,patch)=vec_current_at_patch(i,patch) + current_patch(i)
#ifdef HAVE_TASKLOOP
			!$omp end atomic
#endif
			enddo
		endif
		enddo
	enddo
#ifdef HAVE_OPENMP
	!$omp end parallel do
#endif
	call MPI_ALLREDUCE(MPI_IN_PLACE,vec_current_at_patch,quant%maxpatch*3,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)


	do i=1,quant%maxpatch
	current_at_patch(i) = sqrt(sum(vec_current_at_patch(1:3,i)**2))
	enddo
#ifdef HAVE_OPENMP
    !$omp parallel do default(shared) private(patch,i,ii,node,a)
#endif
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
#ifdef HAVE_OPENMP
    !$omp end parallel do
#endif
	if(ptree%MyID == Main_ID)then
    ! string='current'//chara//'.out'
    open(30,file=string)
    do node=1, quant%maxnode
        write (30,*) current_at_node(node)
    enddo
    close(30)
	endif
	current1(:,1) = current_at_node

!	if(ptree%MyID == Main_ID)then
!	write(*,*)'J norm: ', z_fnorm(current1,quant%maxnode,1,'1')/z_fnorm(current1,quant%maxnode,1,'I')
!	endif

    deallocate (current1,current_at_node,current_at_patch,vec_current_at_patch)
    return

end subroutine current_node_patch_mapping


real(kind=8) function triangle_area(patch,quant)

    use z_BPACK_DEFS
    implicit none

    integer patch,i
    real(kind=8) a(3),b(3),c(3)
	type(quant_EMSURF)::quant

    do i=1,3
        a(i)=quant%xyz(i,quant%node_of_patch(2,patch))-quant%xyz(i,quant%node_of_patch(1,patch))
        b(i)=quant%xyz(i,quant%node_of_patch(3,patch))-quant%xyz(i,quant%node_of_patch(1,patch))
    enddo

    call z_rrcurl(a,b,c)
    triangle_area=0.5*sqrt(c(1)**2+c(2)**2+c(3)**2)

    return
end function triangle_area


subroutine element_Vinc_VV_SURF(theta,phi,edge,value,quant)

    use z_BPACK_DEFS
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


    einc(1)=cos(theta*BPACK_pi/180.)*cos(phi*BPACK_pi/180.)
    einc(2)=cos(theta*BPACK_pi/180.)*sin(phi*BPACK_pi/180.)
    einc(3)=-sin(theta*BPACK_pi/180.)
    k(1)=sin(theta*BPACK_pi/180.)*cos(phi*BPACK_pi/180.)
    k(2)=sin(theta*BPACK_pi/180.)*sin(phi*BPACK_pi/180.)
    k(3)=cos(theta*BPACK_pi/180.)

    value_h=(0,0); value_e=(0,0)

    node1=quant%info_unk(1,edge)
    node2=quant%info_unk(2,edge)
    lm=sqrt((quant%xyz(1,node1)-quant%xyz(1,node2))**2+(quant%xyz(2,node1)-quant%xyz(2,node2))**2+(quant%xyz(3,node1)-quant%xyz(3,node2))**2)
    do jj=3,4
	   nr(1:3)=quant%normal_of_patch(1:3,quant%info_unk(jj,edge))
	   node3=quant%info_unk(jj+2,edge)
       call gau_grobal(edge,jj,x,y,z,w,quant)
	   do ii=1,quant%integral_points
          phase=BPACK_junit*quant%wavenum*(sin(theta*BPACK_pi/180.)*cos(phi*BPACK_pi/180.)*x(ii)+sin(theta*BPACK_pi/180.)*sin(phi*BPACK_pi/180.)*y(ii)+cos(theta*BPACK_pi/180.)*z(ii))
	      do i=1,3
		     ee(i)=einc(i)*exp(phase)
          end do
          call z_rccurl(-k,ee,hh1)
          call z_rccurl(nr,hh1,hh)
          a(1)=x(ii)-quant%xyz(1,node3)
	      a(2)=y(ii)-quant%xyz(2,node3)
	      a(3)=z(ii)-quant%xyz(3,node3)
	      call z_cscalar(hh,a,ctemp_h)
	      call z_cscalar(ee,a,ctemp_e)
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

    use z_BPACK_DEFS
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

    einc(1)=-sin(phi*BPACK_pi/180.)
    einc(2)=cos(phi*BPACK_pi/180.)
    einc(3)=0.
    k(1)=sin(theta*BPACK_pi/180.)*cos(phi*BPACK_pi/180.)
    k(2)=sin(theta*BPACK_pi/180.)*sin(phi*BPACK_pi/180.)
    k(3)=cos(theta*BPACK_pi/180.)

    value_h=(0,0); value_e=(0,0)

    node1=quant%info_unk(1,edge)
    node2=quant%info_unk(2,edge)
    lm=sqrt((quant%xyz(1,node1)-quant%xyz(1,node2))**2+(quant%xyz(2,node1)-quant%xyz(2,node2))**2+(quant%xyz(3,node1)-quant%xyz(3,node2))**2)
    do jj=3,4
	   nr(1:3)=quant%normal_of_patch(1:3,quant%info_unk(jj,edge))
	   node3=quant%info_unk(jj+2,edge)
       call gau_grobal(edge,jj,x,y,z,w,quant)
	   do ii=1,quant%integral_points
          phase=BPACK_junit*quant%wavenum*(sin(theta*BPACK_pi/180.)*cos(phi*BPACK_pi/180.)*x(ii)+sin(theta*BPACK_pi/180.)*sin(phi*BPACK_pi/180.)*y(ii)+cos(theta*BPACK_pi/180.)*z(ii))
	      do i=1,3
		     ee(i)=einc(i)*exp(phase)
          end do
          call z_rccurl(-k,ee,hh1)
          call z_rccurl(nr,hh1,hh)
          a(1)=x(ii)-quant%xyz(1,node3)
	      a(2)=y(ii)-quant%xyz(2,node3)
	      a(3)=z(ii)-quant%xyz(3,node3)
	      call z_cscalar(hh,a,ctemp_h)
	      call z_cscalar(ee,a,ctemp_e)
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
	use z_BPACK_DEFS
    implicit none

    real(kind=8) rcs
    complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1
    real(kind=8) theta,phi,dtheta,dphi

    integer i,j,ii,jj,iii,jjj,patch,flag
    real(kind=8) l_edge,l_edgefine
    real(kind=8) a(3)
    real(kind=8),allocatable:: x(:),y(:),z(:),w(:)
    complex(kind=8):: curr(:,:)
	type(z_mesh)::msh
	type(quant_EMSURF)::quant
	type(z_proctree)::ptree

    integer edge,edge_m,edge_n,ierr

    if(ptree%MyID==Main_ID)open (100, file='VV_bistatic.txt')
    theta=90.
    dphi=180./quant%RCS_Nsample

    do i=0, quant%RCS_Nsample
        phi=i*dphi
        ctemp_loc=(0.,0.)
        rcs=0
#ifdef HAVE_OPENMP
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
#endif
		do edge=msh%idxs,msh%idxe
            call VV_polar_SURF(theta,phi,msh%new2old(edge),ctemp_1,curr(edge-msh%idxs+1,1),quant)
            ctemp_loc=ctemp_loc+ctemp_1
        enddo
#ifdef HAVE_OPENMP
        !$omp end parallel do
#endif
		call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)

        rcs=(abs(quant%wavenum*ctemp))**2/4/BPACK_pi
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
#ifdef HAVE_OPENMP
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
#endif
		do edge=msh%idxs,msh%idxe
            call HH_polar_SURF(theta,phi,msh%new2old(edge),ctemp_1,curr(edge-msh%idxs+1,2),quant)
            ctemp_loc=ctemp_loc+ctemp_1
        enddo
#ifdef HAVE_OPENMP
        !$omp end parallel do
#endif
		call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)
        rcs=(abs(quant%wavenum*ctemp))**2/4/BPACK_pi
        !rcs=rcs/quant%wavelength
        rcs=10*log10(rcs)
        if(ptree%MyID==Main_ID)write(1000,*)phi,rcs
    enddo
    if(ptree%MyID==Main_ID)close(1000)



    return

end subroutine RCS_bistatic_SURF

subroutine VV_polar_SURF(theta,phi,edge,ctemp_1,curr,quant)

    use z_BPACK_DEFS
    implicit none

    complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1
    real(kind=8) theta,phi
    type(quant_EMSURF)::quant
    integer i,j,ii,jj,iii,jjj,patch,flag
    real(kind=8) l_edge,l_edgefine
    real(kind=8) a(3)
    complex(kind=8)::curr
    integer edge,edge_m,edge_n
    ! type(z_mesh)::quant
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
		            phase=BPACK_junit*quant%wavenum*(x(ii)*sin(theta*BPACK_pi/180.)*cos(phi*BPACK_pi/180.)+y(ii)*sin(theta*BPACK_pi/180.)*sin(phi*BPACK_pi/180.)+z(ii)*cos(theta*BPACK_pi/180.))
		            do i=1,3
                        ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*curr*exp(phase)*w(ii)
                        !ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*vectors_block(0)%vector(edge,1)*exp(phase)*w(ii)  !current
		            enddo
                enddo
            enddo
            ctemp_1=BPACK_impedence0*(cos(theta*BPACK_pi/180.)*cos(phi*BPACK_pi/180.)*ctemp_rcs(1)+cos(theta*BPACK_pi/180.)*sin(phi*BPACK_pi/180.)*ctemp_rcs(2)-sin(theta*BPACK_pi/180.)*ctemp_rcs(3))

	deallocate(x,y,z,w)

    return

end subroutine VV_polar_SURF

subroutine HH_polar_SURF(theta,phi,edge,ctemp_1,curr,quant)

    use z_BPACK_DEFS
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
		            phase=BPACK_junit*quant%wavenum*(x(ii)*sin(theta*BPACK_pi/180.)*cos(phi*BPACK_pi/180.)+y(ii)*sin(theta*BPACK_pi/180.)*sin(phi*BPACK_pi/180.)+z(ii)*cos(theta*BPACK_pi/180.))
		            do i=1,3
                        ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*curr*exp(phase)*w(ii)
                        !ctemp_rcs(i)=ctemp_rcs(i)+(-1)**(jj+1)*l_edge*a(i)*vectors_block(0)%vector(edge,1)*exp(phase)*w(ii)  !current
		            enddo
                enddo
	        enddo
            ctemp_1=BPACK_impedence0*(-sin(phi*BPACK_pi/180.)*ctemp_rcs(1)+cos(phi*BPACK_pi/180.)*ctemp_rcs(2))
	deallocate(x,y,z,w)

    return

end subroutine HH_polar_SURF




subroutine RCS_monostatic_VV_SURF(dsita,dphi,rcs,curr,msh,quant,ptree)

    use z_BPACK_DEFS
    implicit none
    complex(kind=8)::curr(:)
    real(kind=8) rcs
    complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1,ctemp_2
    real(kind=8) dsita,dphi
    integer edge,edge_m,edge_n,ierr
	type(z_mesh)::msh
    type(quant_EMSURF)::quant
	type(z_proctree)::ptree

    ctemp_loc=(0.,0.)
        rcs=0
#ifdef HAVE_OPENMP
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
#endif
		do edge=msh%idxs,msh%idxe
            call VV_polar_SURF(dsita,dphi,msh%new2old(edge),ctemp_1,curr(edge-msh%idxs+1),quant)
	        ctemp_loc=ctemp_loc+ctemp_1
        enddo
#ifdef HAVE_OPENMP
        !$omp end parallel do
#endif
		call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)

        rcs=(abs(quant%wavenum*ctemp))**2/4/BPACK_pi
        !rcs=rcs/quant%wavelength
        rcs=10*log10(rcs)

    return

end subroutine RCS_monostatic_VV_SURF

subroutine RCS_monostatic_HH_SURF(dsita,dphi,rcs,curr,msh,quant,ptree)

    use z_BPACK_DEFS
    implicit none

    real(kind=8) rcs
    complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1,ctemp_2
    real(kind=8) dsita,dphi
    integer edge,edge_m,edge_n,ierr
    complex(kind=8)::curr(:)
	type(z_mesh)::msh
	type(quant_EMSURF)::quant
	type(z_proctree)::ptree

    ctemp_loc=(0.,0.)
        rcs=0
#ifdef HAVE_OPENMP
        !$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
#endif
		do edge=msh%idxs,msh%idxe
            call HH_polar_SURF(dsita,dphi,msh%new2old(edge),ctemp_1,curr(edge-msh%idxs+1),quant)
	        ctemp_loc=ctemp_loc+ctemp_1
        enddo
#ifdef HAVE_OPENMP
        !$omp end parallel do
#endif
		call MPI_ALLREDUCE(ctemp_loc,ctemp,1,MPI_DOUBLE_COMPLEX,MPI_SUM,ptree%Comm,ierr)

        rcs=(abs(quant%wavenum*ctemp))**2/4/BPACK_pi
        !rcs=rcs/quant%wavelength
        rcs=10*log10(rcs)

    return

end subroutine RCS_monostatic_HH_SURF




subroutine geo_modeling_SURF(quant,MPIcomm,DATA_DIR)
    use z_BPACK_DEFS
	use z_MISC_Utilities
    implicit none

    integer i,j,ii,jj,nn,iii,jjj
    integer intemp
    integer node, patch, edge,edge1, flag
    integer node1, node2,found
    integer node_temp(2)
    real(kind=8) a(3),b(3),c(3),r0
	! type(z_mesh)::msh
	type(quant_EMSURF)::quant
	! type(z_proctree)::ptree
	integer MPIcomm,MyID,ierr
	CHARACTER (*) DATA_DIR
	integer Maxedge
	integer,parameter:: NperNode=10
	integer,allocatable::tmpint(:,:)
	integer,allocatable:: info_unk(:,:)
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
#ifdef HAVE_OPENMP
    !$omp parallel do default(shared) private(patch,i,a,b,c,r0)
#endif
	do patch=1,quant%maxpatch
        do i=1,3
            a(i)=(quant%xyz(i,quant%node_of_patch(2,patch))-quant%xyz(i,quant%node_of_patch(1,patch)))
            b(i)=(quant%xyz(i,quant%node_of_patch(3,patch))-quant%xyz(i,quant%node_of_patch(1,patch)))
        enddo
        call z_rrcurl(a,b,c)
        r0=sqrt(c(1)**2+c(2)**2+c(3)**2)
        c(1)=c(1)/r0
        c(2)=c(2)/r0
        c(3)=c(3)/r0
        quant%normal_of_patch(1:3,patch)=c(1:3)
    enddo
#ifdef HAVE_OPENMP
    !$omp end parallel do
#endif
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
#ifdef HAVE_OPENMP
    !$omp parallel do default(shared) private(edge,node_temp,jj,iii,jjj)
#endif
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
#ifdef HAVE_OPENMP
    !$omp end parallel do
#endif



	allocate(quant%edge_of_patch(quant%maxpatch,3,2))
	quant%edge_of_patch=-1
	do edge=1,maxedge
	    do jj=3,4
			patch=quant%info_unk(jj,edge)
			do iii=1,3
				if(quant%edge_of_patch(patch,iii,1)==edge .or. quant%edge_of_patch(patch,iii,1)==-1)then
					quant%edge_of_patch(patch,iii,1) = edge
					quant%edge_of_patch(patch,iii,2) = jj
					exit
				endif
			enddo
		enddo
    enddo


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


	quant%minedgelength = BPACK_Bigvalue
	do edge=1,Maxedge
		quant%minedgelength = min(quant%minedgelength,sqrt(sum(abs(quant%xyz(:,quant%info_unk(1,edge))-quant%xyz(:,quant%info_unk(2,edge)))**2)))
	end do

	! write(*,*)	quant%xyz(1,1:100),sum(quant%xyz(1,:))
	! stop
	quant%Nunk = Maxedge

    if(MyID==Main_ID)write (*,*) ''
    if(MyID==Main_ID)write (*,*) 'Maxedge:',Maxedge
	if(MyID==Main_ID)write (*,*) 'minedgelength:',quant%minedgelength
	if(MyID==Main_ID)write (*,*) 'wavelength/minedgelength:',quant%wavelength/quant%minedgelength
	if(MyID==Main_ID)write (*,*) 'maxedgelength:',quant%maxedgelength
	if(MyID==Main_ID)write (*,*) 'wavelength/maxedgelength:',quant%wavelength/quant%maxedgelength

    if(MyID==Main_ID)write (*,*) ''

    return

end subroutine geo_modeling_SURF


subroutine EM_solve_SURF(bmat,option,msh,quant,ptree,stats)
    use z_BPACK_DEFS
	use z_BPACK_Solve_Mul


    implicit none

    integer i, j, ii, jj, iii, jjj,ierr
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, n_iter_max, iter ,N_unk, N_unk_loc
    real(kind=8) theta, phi, dphi, rcs_V, rcs_H
    real T0
    real(kind=8) n1,n2,rtemp,n3,n4
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:)
	real(kind=8):: rel_error
	type(z_Hoption)::option
	type(z_Bmatrix)::bmat
	type(z_mesh)::msh
	type(quant_EMSURF)::quant
	type(z_proctree)::ptree
	type(z_Hstat)::stats
	complex(kind=8),allocatable:: current(:,:),voltage(:,:)

	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve......"

	N_unk_loc = msh%idxe-msh%idxs+1

    if (quant%RCS_static==2) then

        theta=90
        phi=0

        allocate (current(N_unk_loc,2))
		Current=0
        allocate (voltage(N_unk_loc,2))

#ifdef HAVE_OPENMP
        !$omp parallel do default(shared) private(edge,value_Z)
#endif
		do edge=msh%idxs, msh%idxe
            call element_Vinc_VV_SURF(theta,phi,msh%new2old(edge),value_Z,quant)
			voltage(edge-msh%idxs+1,1)=value_Z
			call element_Vinc_HH_SURF(theta,phi,msh%new2old(edge),value_Z,quant)
            voltage(edge-msh%idxs+1,2)=value_Z
        enddo
#ifdef HAVE_OPENMP
        !$omp end parallel do
#endif
        n1 = MPI_Wtime()
		call z_BPACK_Solution(bmat,Current,Voltage,N_unk_loc,2,option,ptree,stats)
		n2 = MPI_Wtime()

        if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''
        if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'Solving:',n2-n1,'Seconds'
        if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''

		n1 = MPI_Wtime()
        call RCS_bistatic_SURF(Current,msh,quant,ptree)
		n2 = MPI_Wtime()
        if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''
        if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'Bistatic RCS',n2-n1,'Seconds'
        if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''

		! n3 = MPI_Wtime()
		! call current_node_patch_mapping('V_current',Current(:,1),msh,quant,ptree)
		! call current_node_patch_mapping('H_current',Current(:,2),msh,quant,ptree)
		! n4 = MPI_Wtime()
        ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''
        ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) 'Current dumping',n4-n3,'Seconds'
        ! if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,*) ''

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


        if(ptree%MyID==Main_ID)open (100, file='monostaticH.out')

        n1=MPI_Wtime()

        do j=0, num_sample
            phi=j*dphi
#ifdef HAVE_OPENMP
			!$omp parallel do default(shared) private(edge,value_Z)
#endif
			do edge=msh%idxs, msh%idxe
				call element_Vinc_HH_SURF(theta,phi,msh%new2old(edge),value_Z,quant)
				b(edge-msh%idxs+1,j+1)=value_Z
			enddo
#ifdef HAVE_OPENMP
			!$omp end parallel do
#endif
        enddo

		call z_BPACK_Solution(bmat,x,b,N_unk_loc,num_sample+1,option,ptree,stats)

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

		n2 = MPI_Wtime()
		stats%Time_Sol = stats%Time_Sol + n2-n1
		call MPI_ALLREDUCE(stats%Time_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)

        if(ptree%MyID==Main_ID)then
			close(100)
			write (*,*) ''
			write (*,*) 'Solving:',rtemp,'Seconds'
			write (*,*) ''
		endif

		call MPI_ALLREDUCE(stats%Flop_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
		if(ptree%MyID==Main_ID .and. option%verbosity>=0)write (*,'(A13,Es14.2)') 'Solve flops:',rtemp


		deallocate(Current)

    endif

	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "EM_solve finished"
	if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "    "

    return

end subroutine EM_solve_SURF





end module EMSURF_MODULE
