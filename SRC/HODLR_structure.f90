#include "HODLR_config.fi"
module HODLR_structure 
use HODLR_Utilities
contains 


real(kind=8) function euclidean_distance(node1,node2,msh)
    
    use HODLR_DEFS
    implicit none

    integer node1, node2
    real(kind=8) dis
    integer i, j
    integer Dimn
	type(mesh)::msh
	
	Dimn = 0
	if(allocated(msh%xyz))Dimn = size(msh%xyz,1)
	
    dis=0d0
    do i=1,Dimn
        dis=dis+(msh%xyz(i,node1)-msh%xyz(i,node2))**2
    enddo
    
    euclidean_distance=dis
    
    return
    
end function euclidean_distance


!**** l2 gram distance^2 between element edgem and edgen is 
!     defined as: Re{Z_ii+Z_jj-Z_ij-Z_ji} for SPD, HPD, general symmetric real, and hermitian matrices
!     undefined otherwise
!**** angular gram distance^2 is  
!     defined as 1-Z_ij^2/(Z_iiZ_jj)
!     undefined otherwise
!     Use with caution !!!
real(kind=8) function gram_distance(edgem,edgen,ker,msh,element_Zmn)
    
    use HODLR_DEFS
    implicit none

    integer edgem, edgen
	type(mesh)::msh
	type(kernelquant)::ker
	procedure(Zelem)::element_Zmn
	DT r1,r2,r3,r4
! l2 distance	
#if 0	
	call element_Zmn(edgem,edgem,r1,msh,ker)
	call element_Zmn(edgen,edgen,r2,msh,ker)
	call element_Zmn(edgem,edgen,r3,msh,ker)
	call element_Zmn(edgen,edgem,r4,msh,ker)
    gram_distance=dble(r1+r2-r3-r4)
! angular distance
#else   
	call element_Zmn(edgem,edgem,r1,msh,ker)
	call element_Zmn(edgen,edgen,r2,msh,ker)
	call element_Zmn(edgem,edgen,r3,msh,ker)
    gram_distance=dble(1d0-r3**2d0/(r1*r2))	
#endif    
    return
    
end function gram_distance


subroutine HODLR_structuring(ho_bf1,option,msh,ker,element_Zmn,ptree)
    
    use HODLR_DEFS
	use misc
    implicit none
    
	integer Cflag
    integer i, j, ii, jj, iii, jjj
    integer level, edge, node, patch, group, group_m, group_n,col_group,row_group,fidx
    integer blocks
    integer center_edge
    
    integer index_temp
    DT r1,r2,r3,r4
    real(kind=8) a, b, c, d, para, xmax,xmin,ymax,ymin,zmax,zmin,seperator,r,theta,phi,phi_tmp
    real(kind=8) radius, radiusmax, radius2, radiusmax2
    real(kind=8),allocatable:: xyzrange(:),xyzmin(:),xyzmax(:),auxpoint(:),groupcenter(:)
    real(kind=8), allocatable :: distance(:),array(:,:)
	integer level_c,sortdirec,mm,phi_end,Ninfo_edge,ind_i,ind_j
	real(kind=8) t1,t2
	integer Maxgroup,nlevel_pre
	character(len=1024)  :: strings	
    integer, allocatable :: order(:), edge_temp(:,:),map_temp(:)
	integer dimn,groupsize,idxstart,Nsmp
	type(Hoption)::option
	type(hobf)::ho_bf1
	integer Maxlevel,Maxgrp
	type(mesh)::msh
	type(kernelquant)::ker
	type(proctree)::ptree
    procedure(Zelem)::element_Zmn
	integer,allocatable:: perms(:)
	
	

	!*************Initialize permutation vector ********	   
	allocate(msh%new2old(msh%Nunk))    
	do ii=1,msh%Nunk
		msh%new2old(ii) = ii
	end do
	
	
	!************Compute Maxlevel of hodlr tree*******************	
	nlevel_pre=0
	if(allocated(msh%pretree))then
		nlevel_pre = ceiling_safe(log(dble(size(msh%pretree,1))) / log(2d0))
	endif
	level=0; i=1
	do while (int(msh%Nunk/i)>option%Nmin_leaf)
		level=level+1
		i=2**level
	enddo
	Maxlevel=level
	if(Maxlevel<nlevel_pre)Maxlevel=nlevel_pre
	ho_bf1%Maxlevel=Maxlevel
	
	
	!************** check whether the sorting option is valid
	if(Maxlevel>nlevel_pre)then
		if(.not. allocated(msh%xyz))then
			if(option%xyzsort==CKD .or. option%xyzsort==TM)then
				write(*,*)'Geometrical information is not provided. Try use NATRUAL or TM_GRAM ordering'
				stop
			endif
		endif
		
		if(option%xyzsort==TM_GRAM)then
			call random_number(a)
			ind_i =floor_safe(a*(msh%Nunk-1))+1
			ind_j = ind_i
			do while(ind_i==ind_j)
				call random_number(a)
				ind_j =floor_safe(a*(msh%Nunk-1))+1
			enddo
			call element_Zmn(ind_i,ind_i,r1,msh,ker)
			call element_Zmn(ind_j,ind_j,r2,msh,ker)
			call element_Zmn(ind_i,ind_j,r3,msh,ker)
			call element_Zmn(ind_j,ind_i,r4,msh,ker)
			if(aimag(cmplx(r1,kind=8))/=0 .or. aimag(cmplx(r2,kind=8))/=0 .or. abs(r3-conjg(cmplx(r3,kind=8)))>abs(r3)*SafeEps)then
				write(*,*)'Matrix not hermitian. The gram distance is undefined'
			endif
		endif		
	endif
	
	!***************************************************

	Maxgroup=2**(Maxlevel+1)-1
	msh%Maxgroup = Maxgroup
	allocate (msh%basis_group(Maxgroup))
	if (ptree%MyID==Main_ID .and. option%verbosity>0)then
		write (*,*) ''
		write (*,*) 'Maxlevel_for_blocks:',ho_bf1%Maxlevel
		write (*,*) 'N_leaf:',int(msh%Nunk/(2**Maxlevel))
		write (*,*) ''
		write (*,*) 'Constructing basis groups...'
	endif
	
	!**** construct the top few levels whose ordering is provided by the user
	msh%basis_group(1)%head=1 ; msh%basis_group(1)%tail=msh%Nunk; msh%basis_group(1)%pgno=1
	do level=nlevel_pre,0,-1
		idxstart=1
		do group=2**level, 2**(level+1)-1
			! msh%basis_group(group)%level=level

			if(level==nlevel_pre)then
				if(nlevel_pre==0)then
					groupsize = msh%Nunk
				else
					groupsize = msh%pretree(group-2**nlevel_pre+1)				
				endif
				call assert(groupsize>0,'zero leafsize may not be handled')
				msh%basis_group(group)%head = idxstart
				msh%basis_group(group)%tail = idxstart + groupsize -1
				idxstart = idxstart + groupsize
			else
				msh%basis_group(group)%head = msh%basis_group(2*group)%head
				msh%basis_group(group)%tail = msh%basis_group(2*group+1)%tail
			endif					
		enddo
	enddo	
	
	
	!**** if necessary, continue ordering the sub-trees using clustering method specified by option%xyzsort
	dimn=0
	if(allocated(msh%xyz))Dimn = size(msh%xyz,1)
	if(dimn>0)then
	allocate(xyzrange(dimn))
	allocate(xyzmin(dimn))
	allocate(xyzmax(dimn))
	allocate(auxpoint(dimn))
	allocate(groupcenter(dimn))
	endif	

	do level=nlevel_pre, Maxlevel
		do group=2**level, 2**(level+1)-1
			! msh%basis_group(group)%level=level
								
			if(option%xyzsort==NATURAL)then !natural ordering		 
				mm = msh%basis_group(group)%tail - msh%basis_group(group)%head + 1
				allocate (distance(mm))	
				do i=msh%basis_group(group)%head, msh%basis_group(group)%tail
					distance(i-msh%basis_group(group)%head+1)=dble(i)
				enddo
		
			else if(option%xyzsort==CKD)then !msh%xyz sort		 
				xyzmin= 1d300
				xyzmax= -1d300
				do edge=msh%basis_group(group)%head, msh%basis_group(group)%tail
					do ii=1,Dimn
						xyzmax(ii) = max(xyzmax(ii),msh%xyz(ii,msh%new2old(edge)))
						xyzmin(ii) = min(xyzmin(ii),msh%xyz(ii,msh%new2old(edge)))
					enddo
				enddo
				xyzrange(1:Dimn) = xyzmax(1:Dimn)-xyzmin(1:Dimn)
				
				mm = msh%basis_group(group)%tail - msh%basis_group(group)%head + 1
				allocate (distance(mm))	
				sortdirec = maxloc(xyzrange(1:Dimn),1)
				! write(*,*)'gaw',sortdirec,xyzrange(1:Dimn)
				
				! if(ker%Kernel==EMSURF)then
				! if(mod(level,2)==1)then           !!!!!!!!!!!!!!!!!!!!!!!!! note: applys only to plates
					! sortdirec=2
				! else 
					! sortdirec=3
				! end if
				! endif
				
				
				!$omp parallel do default(shared) private(i)
				do i=msh%basis_group(group)%head, msh%basis_group(group)%tail
					distance(i-msh%basis_group(group)%head+1)=msh%xyz(sortdirec,msh%new2old(i))
				enddo
				!$omp end parallel do

			else if(option%xyzsort==TM)then !cobblestone sort

				groupcenter(1:dimn)=0.0d0
				! !$omp parallel do default(shared) private(edge,ii) reduction(+:groupcenter)
				do edge=msh%basis_group(group)%head, msh%basis_group(group)%tail
					do ii=1,dimn
						groupcenter(ii)=groupcenter(ii)+msh%xyz(ii,msh%new2old(edge))
					enddo
				enddo
				! !$omp end parallel do
				do ii=1,dimn
					groupcenter(ii)=groupcenter(ii)/(msh%basis_group(group)%tail-msh%basis_group(group)%head+1)
				enddo

				radiusmax=0.
				do edge=msh%basis_group(group)%head, msh%basis_group(group)%tail
					radius=0
					do ii=1,dimn
						radius=radius+(msh%xyz(ii,msh%new2old(edge))-groupcenter(ii))**2
					enddo
					radius=sqrt(radius)
					if (radius>radiusmax) then
						radiusmax=radius
						center_edge=edge
					endif			
				enddo

				mm = msh%basis_group(group)%tail - msh%basis_group(group)%head + 1
				allocate (distance(mm))	
				
				distance(1:mm)=Bigvalue
				!$omp parallel do default(shared) private(i)
				do i=msh%basis_group(group)%head, msh%basis_group(group)%tail
					distance(i-msh%basis_group(group)%head+1)=euclidean_distance(msh%new2old(i),msh%new2old(center_edge),msh)
				enddo
				!$omp end parallel do					

			else if(option%xyzsort==TM_GRAM)then !GRAM-distance-based cobblestone sort

				Nsmp = min(msh%basis_group(group)%tail-msh%basis_group(group)%head+1,50)
				allocate(perms(msh%basis_group(group)%tail-msh%basis_group(group)%head+1))
				call rperm(msh%basis_group(group)%tail-msh%basis_group(group)%head+1, perms)
				
				radiusmax2=0.
				do edge=msh%basis_group(group)%head, msh%basis_group(group)%tail
					radius2=0
					do ii=1,Nsmp  ! take average of distance^2 to Nsmp samples as the distance^2 to the group center 
						radius2 = radius2 + gram_distance(edge,perms(ii)+msh%basis_group(group)%head-1,ker,msh,element_Zmn)
					enddo
					! call assert(radius2>0,'radius2<0 cannot take square root')
					! radius2 = sqrt(radius2)
					radius2 = radius2/Nsmp
					if (radius2>radiusmax2) then
						radiusmax2=radius2
						center_edge=edge
					endif			
				enddo

				mm = msh%basis_group(group)%tail - msh%basis_group(group)%head + 1
				allocate (distance(mm))	
				
				distance(1:mm)=Bigvalue
				
				
				!$omp parallel do default(shared) private(i)
				do i=msh%basis_group(group)%head, msh%basis_group(group)%tail
					distance(i-msh%basis_group(group)%head+1)=gram_distance(i,center_edge,ker,msh,element_Zmn)
				enddo
				!$omp end parallel do

				deallocate(perms)
				
			end if
			
			allocate (order(mm))
			allocate(map_temp(mm))

			call quick_sort(distance,order,mm)       
			!$omp parallel do default(shared) private(ii)     
			do ii=1, mm
				map_temp(ii) = msh%new2old(order(ii)+msh%basis_group(group)%head-1)
			enddo
			!$omp end parallel do

			!$omp parallel do default(shared) private(ii)     
			do ii=1, mm
				msh%new2old(ii+msh%basis_group(group)%head-1) = map_temp(ii)
			enddo
			!$omp end parallel do			
			

			! deallocate(edge_temp)
			deallocate(map_temp)
			deallocate(order)

			deallocate(distance)
			
			
			if (level<Maxlevel) then
				msh%basis_group(2*group)%head=msh%basis_group(group)%head
				msh%basis_group(2*group)%tail=int((msh%basis_group(group)%head+msh%basis_group(group)%tail)/2)
				msh%basis_group(2*group+1)%head=msh%basis_group(2*group)%tail+1
				msh%basis_group(2*group+1)%tail=msh%basis_group(group)%tail
								
				if(option%xyzsort==CKD)then 
					seperator = msh%xyz(sortdirec,msh%new2old(msh%basis_group(2*group)%tail))
				end if
				
				msh%basis_group(group)%boundary(1) = sortdirec
				msh%basis_group(group)%boundary(2) = seperator
													
			endif	
		enddo
	enddo	
	
   
	if(dimn>0)then
	deallocate(xyzrange)
	deallocate(xyzmin)
	deallocate(xyzmax)
	deallocate(auxpoint)
	deallocate(groupcenter)	
	endif
	
	allocate(msh%old2new(msh%Nunk))
	do ii=1,msh%Nunk
		msh%old2new(msh%new2old(ii)) = ii
	end do		
		
	
	!**********Dump the ordering into a file********************************
	
#if	0	
	write(strings , *) Dimn
	do level=0, Maxlevel
		do group=2**level, 2**(level+1)-1
			do edge=msh%basis_group(group)%head, msh%basis_group(group)%tail
				write(113,'(I5,I8,'//TRIM(strings)//'Es16.8)')level,group,msh%xyz(1:Dimn,msh%new2old(edge))
			enddo
		enddo
	enddo	   
#endif


    return
    
end subroutine HODLR_structuring


subroutine BPlus_structuring(ho_bf1,option,msh,ptree)
	use HODLR_DEFS
	use misc
	implicit none 
	
    integer i, j, ii, jj, kk, iii, jjj,ll,bb,sortdirec,ii_sch
    integer level, edge, patch, node, group, group_touch
    integer rank, index_near, m, n, length, flag, itemp,cnt,detection
    real T0
	real(kind=8):: tolerance, rtemp,rel_error,seperator,dist
    real(kind=8) Memory_direct_forward,Memory_butterfly_forward
	integer mm,nn,header_m,header_n,edge_m,edge_n,group_m,group_n,group_m1,group_n1,group_m2,group_n2,levelm,groupm_start,index_i_m,index_j_m
	integer level_c,iter,level_cc,level_BP,Nboundall,level_butterfly	
	type(matrixblock),pointer::blocks,block_f,block_sch,block_inv
	real(kind=8)::minbound,theta,phi,r,rmax,phi_tmp,measure
	real(kind=8),allocatable::Centroid_M(:,:),Centroid_N(:,:)
	integer,allocatable::Isboundary_M(:),Isboundary_N(:)
	integer Dimn,col_group,row_group,Maxgrp
	type(Hoption)::option
	type(mesh)::msh
	type(hobf)::ho_bf1
	character(len=1024)  :: strings
	type(proctree)::ptree

	Maxgrp=2**(ptree%nlevel)-1
	
	msh%basis_group(1)%pgno=1
	do level=0, ho_bf1%Maxlevel
		do group=2**level, 2**(level+1)-1
			if(level<ho_bf1%Maxlevel)then
			if(msh%basis_group(group)%pgno*2<=Maxgrp)then
				msh%basis_group(2*group)%pgno=msh%basis_group(group)%pgno*2
			else
				msh%basis_group(2*group)%pgno=msh%basis_group(group)%pgno
			endif
			if(msh%basis_group(group)%pgno*2+1<=Maxgrp)then
				msh%basis_group(2*group+1)%pgno=msh%basis_group(group)%pgno*2+1
			else
				msh%basis_group(2*group+1)%pgno=msh%basis_group(group)%pgno
			endif
			endif
		enddo
	enddo

	
	ho_bf1%N=msh%Nunk
	allocate(ho_bf1%levels(ho_bf1%Maxlevel+1))

	
	if(2**ho_bf1%Maxlevel<2**(ptree%nlevel-1))then
		write(*,*)'too many processes for paralleling leaf boxes!'
		stop
	endif	
	
	do level_c = 1,ho_bf1%Maxlevel+1
		ho_bf1%levels(level_c)%level = level_c
		if(level_c == ho_bf1%Maxlevel+1)then
			ho_bf1%levels(level_c)%N_block_forward = 2**(level_c-1)
		else
			ho_bf1%levels(level_c)%N_block_forward = 2**level_c			
		endif
		ho_bf1%levels(level_c)%N_block_inverse = 2**(level_c-1)
		ho_bf1%levels(level_c)%Bidxs = 2**(ho_bf1%Maxlevel+1)
		ho_bf1%levels(level_c)%Bidxe = -2**(ho_bf1%Maxlevel+1)

		allocate(ho_bf1%levels(level_c)%BP(ho_bf1%levels(level_c)%N_block_forward))			
		allocate(ho_bf1%levels(level_c)%BP_inverse(ho_bf1%levels(level_c)%N_block_inverse))
		allocate(ho_bf1%levels(level_c)%BP_inverse_update(ho_bf1%levels(level_c)%N_block_forward))	
		allocate(ho_bf1%levels(level_c)%BP_inverse_schur(ho_bf1%levels(level_c)%N_block_inverse))
	end do

	
	ho_bf1%levels(1)%BP_inverse(1)%level = 0
	ho_bf1%levels(1)%BP_inverse(1)%col_group = 1
	ho_bf1%levels(1)%BP_inverse(1)%row_group = 1
	ho_bf1%levels(1)%BP_inverse(1)%pgno = 1
	! ho_bf1%levels(1)%BP_inverse(1)%style = 2
	
	! treat hodlr as a full matrix if Maxlevel=0
	if(ho_bf1%Maxlevel==0)then 
		ho_bf1%levels(1)%BP(1)%level = 0
		ho_bf1%levels(1)%BP(1)%col_group = 1
		ho_bf1%levels(1)%BP(1)%row_group = 1
		ho_bf1%levels(1)%BP(1)%pgno = 1		
	endif
	
	do level_c = 1,ho_bf1%Maxlevel
		do ii = 1, ho_bf1%levels(level_c)%N_block_inverse
			col_group = ho_bf1%levels(level_c)%BP_inverse(ii)%col_group
			row_group = ho_bf1%levels(level_c)%BP_inverse(ii)%row_group
	
			allocate(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(LplusMax))
			do ll=1,LplusMax
			ho_bf1%levels(level_c)%BP_inverse(ii)%LL(ll)%Nbound = 0
			end do			
			ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%Nbound=1
			
			allocate(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1))
			block_inv =>ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)
			block_inv%col_group=col_group
			block_inv%row_group=row_group
			block_inv%level=ho_bf1%levels(level_c)%BP_inverse(ii)%level
			block_inv%pgno=ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
			block_inv%headm=msh%basis_group(row_group)%head
			block_inv%headn=msh%basis_group(col_group)%head
			block_inv%M=msh%basis_group(row_group)%tail-msh%basis_group(row_group)%head+1
			block_inv%N=msh%basis_group(col_group)%tail-msh%basis_group(col_group)%head+1
			
			block_inv%level_butterfly = ho_bf1%Maxlevel - block_inv%level
			
			call ComputeParallelIndices(ho_bf1%Maxlevel,block_inv,block_inv%pgno,ptree,msh,0)

			if(IOwnPgrp(ptree,block_inv%pgno))then
				ho_bf1%levels(level_c)%Bidxs = min(ho_bf1%levels(level_c)%Bidxs,ii)
				ho_bf1%levels(level_c)%Bidxe = max(ho_bf1%levels(level_c)%Bidxe,ii)
			endif
			
			if(GetTreelevel(ho_bf1%levels(level_c)%BP_inverse(ii)%pgno)==ptree%nlevel)then
				ho_bf1%levels(level_c)%BP(ii*2-1)%pgno=ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
				ho_bf1%levels(level_c)%BP(ii*2)%pgno=ho_bf1%levels(level_c)%BP_inverse(ii)%pgno				
				ho_bf1%levels(level_c+1)%BP_inverse(ii*2-1)%pgno=ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
				ho_bf1%levels(level_c+1)%BP_inverse(ii*2)%pgno=ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
				ho_bf1%levels(level_c)%BP_inverse_schur(ii)%pgno=ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
			else
				ho_bf1%levels(level_c)%BP(ii*2-1)%pgno=ho_bf1%levels(level_c)%BP_inverse(ii)%pgno*2
				ho_bf1%levels(level_c)%BP(ii*2)%pgno=ho_bf1%levels(level_c)%BP_inverse(ii)%pgno*2+1
				ho_bf1%levels(level_c+1)%BP_inverse(ii*2-1)%pgno=ho_bf1%levels(level_c)%BP_inverse(ii)%pgno*2
				ho_bf1%levels(level_c+1)%BP_inverse(ii*2)%pgno=ho_bf1%levels(level_c)%BP_inverse(ii)%pgno*2+1				
				ho_bf1%levels(level_c)%BP_inverse_schur(ii)%pgno=ho_bf1%levels(level_c)%BP_inverse(ii)%pgno					
			endif
	
			! off-diagonal blocks and their updates
			ho_bf1%levels(level_c)%BP(ii*2-1)%level = level_c
			ho_bf1%levels(level_c)%BP(ii*2-1)%col_group = col_group*2+1
			ho_bf1%levels(level_c)%BP(ii*2-1)%row_group = row_group*2
			ho_bf1%levels(level_c)%BP(ii*2)%level = level_c
			ho_bf1%levels(level_c)%BP(ii*2)%col_group = col_group*2
			ho_bf1%levels(level_c)%BP(ii*2)%row_group = row_group*2+1
			
			! schur complement of every two off-diagonal blocks
			ho_bf1%levels(level_c)%BP_inverse_schur(ii)%level = level_c+1
			ho_bf1%levels(level_c)%BP_inverse_schur(ii)%col_group = col_group * 2
			ho_bf1%levels(level_c)%BP_inverse_schur(ii)%row_group = row_group * 2			
			ho_bf1%levels(level_c)%BP_inverse_schur(ii)%Lplus = 1
	
			! diagonal blocks and their inverses at bottom level
			if(level_c==ho_bf1%Maxlevel)then
				ho_bf1%levels(level_c+1)%BP(ii*2-1)%level = level_c+1
				ho_bf1%levels(level_c+1)%BP(ii*2-1)%col_group = col_group*2
				ho_bf1%levels(level_c+1)%BP(ii*2-1)%row_group = row_group*2
				ho_bf1%levels(level_c+1)%BP(ii*2)%level = level_c+1
				ho_bf1%levels(level_c+1)%BP(ii*2)%col_group = col_group*2+1
				ho_bf1%levels(level_c+1)%BP(ii*2)%row_group = row_group*2+1				
				if(GetTreelevel(ho_bf1%levels(level_c)%BP_inverse(ii)%pgno)==ptree%nlevel)then
					ho_bf1%levels(level_c+1)%BP(ii*2-1)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
					ho_bf1%levels(level_c+1)%BP(ii*2)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
				else
					ho_bf1%levels(level_c+1)%BP(ii*2-1)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno*2
					ho_bf1%levels(level_c+1)%BP(ii*2)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno*2+1					
				endif
			end if
	
			
			ho_bf1%levels(level_c+1)%BP_inverse(ii*2-1)%level = level_c
			ho_bf1%levels(level_c+1)%BP_inverse(ii*2-1)%col_group = col_group*2
			ho_bf1%levels(level_c+1)%BP_inverse(ii*2-1)%row_group = row_group*2
			ho_bf1%levels(level_c+1)%BP_inverse(ii*2)%level = level_c
			ho_bf1%levels(level_c+1)%BP_inverse(ii*2)%col_group = col_group*2+1
			ho_bf1%levels(level_c+1)%BP_inverse(ii*2)%row_group = row_group*2+1				
		end do
	end do

	! do level_c = 1,ho_bf1%Maxlevel+1
	! deallocate(ho_bf1%levels(level_c)%BP_inverse)
	! enddo	
	
	
	Dimn = 0
	if(allocated(msh%xyz))Dimn = size(msh%xyz,1)
	
	do level_c = 1,ho_bf1%Maxlevel+1
		do ii =1,ho_bf1%levels(level_c)%N_block_forward
			! if(IOwnPgrp(ptree,ho_bf1%levels(level_c)%BP(ii)%pgno))then	
			if(level_c==ho_bf1%Maxlevel+1)then
			
				! bottom level dense blocks 
				ho_bf1%levels(level_c)%BP(ii)%Lplus=1				
				allocate(ho_bf1%levels(level_c)%BP(ii)%LL(LplusMax))
				do ll=1,LplusMax
					ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound=0
				end do
				ho_bf1%levels(level_c)%BP(ii)%LL(1)%Nbound = 1
				allocate(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1))
				block_f => ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)
				block_f%level = ho_bf1%levels(level_c)%BP(ii)%level
				block_f%col_group = ho_bf1%levels(level_c)%BP(ii)%col_group
				block_f%row_group = ho_bf1%levels(level_c)%BP(ii)%row_group
				block_f%style = 1  !!!!! be careful here
				block_f%pgno = msh%basis_group(block_f%row_group)%pgno
				
				block_f%M = msh%basis_group(block_f%row_group)%tail - msh%basis_group(block_f%row_group)%head + 1
				block_f%N = msh%basis_group(block_f%col_group)%tail - msh%basis_group(block_f%col_group)%head + 1
				block_f%headm = msh%basis_group(block_f%row_group)%head
				block_f%headn = msh%basis_group(block_f%col_group)%head
				block_f%level_butterfly=0
				call ComputeParallelIndices(ho_bf1%Maxlevel+1,block_f,block_f%pgno,ptree,msh,0)
				
				! bottom level dense blocks' inverse 
				ho_bf1%levels(level_c)%BP_inverse(ii)%Lplus=1				
				allocate(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(LplusMax))
				do ll=1,LplusMax
					ho_bf1%levels(level_c)%BP_inverse(ii)%LL(ll)%Nbound=0
				end do
				ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%Nbound = 1
				allocate(ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1))
				block_inv => ho_bf1%levels(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)
				block_inv%level = ho_bf1%levels(level_c)%BP_inverse(ii)%level
				block_inv%col_group = ho_bf1%levels(level_c)%BP_inverse(ii)%col_group
				block_inv%row_group = ho_bf1%levels(level_c)%BP_inverse(ii)%row_group
				block_inv%style = 1  !!!!! be careful here
				block_inv%pgno = msh%basis_group(block_inv%row_group)%pgno
				
				block_inv%M = msh%basis_group(block_inv%row_group)%tail - msh%basis_group(block_inv%row_group)%head + 1
				block_inv%N = msh%basis_group(block_inv%col_group)%tail - msh%basis_group(block_inv%col_group)%head + 1
				block_inv%headm = msh%basis_group(block_inv%row_group)%head
				block_inv%headn = msh%basis_group(block_inv%col_group)%head
				block_inv%level_butterfly=0
				call ComputeParallelIndices(ho_bf1%Maxlevel+1,block_inv,block_inv%pgno,ptree,msh,0)			
				if(IOwnPgrp(ptree,block_inv%pgno))then
					ho_bf1%levels(level_c)%Bidxs = min(ho_bf1%levels(level_c)%Bidxs,ii)
					ho_bf1%levels(level_c)%Bidxe = max(ho_bf1%levels(level_c)%Bidxe,ii)
				endif					
			else 
				allocate(ho_bf1%levels(level_c)%BP(ii)%LL(LplusMax))
				do ll=1,LplusMax
					ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound=0
				end do
				
				ho_bf1%levels(level_c)%BP(ii)%LL(1)%Nbound = 1
				allocate(ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1))
				block_f => ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)
				block_f%level = ho_bf1%levels(level_c)%BP(ii)%level
				
				if(level_c>option%LRlevel)then
					block_f%level_butterfly = 0 ! low rank below LRlevel
				else 
					if(ho_bf1%Maxlevel - block_f%level<option%lnoBP)then				
						block_f%level_butterfly = ho_bf1%Maxlevel - ho_bf1%levels(level_c)%BP(ii)%level   ! butterfly 
					else
						block_f%level_butterfly = int((ho_bf1%Maxlevel - ho_bf1%levels(level_c)%BP(ii)%level)/2)*2 ! butterfly plus needs even number of levels				
					endif
				endif
				
				block_f%col_group = ho_bf1%levels(level_c)%BP(ii)%col_group
				block_f%row_group = ho_bf1%levels(level_c)%BP(ii)%row_group
				block_f%pgno=msh%basis_group(block_f%row_group)%pgno

				
				
				! compute the partial indices when BP is shared by double number of processes
				ii_sch = ceiling_safe(ii/2d0)
				block_inv => ho_bf1%levels(level_c)%BP_inverse(ii_sch)%LL(1)%matrices_block(1)
				block_f%pgno_db = block_inv%pgno	 	
	
						
		
				
				block_f%M = msh%basis_group(block_f%row_group)%tail - msh%basis_group(block_f%row_group)%head + 1
				block_f%N = msh%basis_group(block_f%col_group)%tail - msh%basis_group(block_f%col_group)%head + 1
				block_f%headm = msh%basis_group(block_f%row_group)%head
				block_f%headn = msh%basis_group(block_f%col_group)%head				
				
				call ComputeParallelIndices(ho_bf1%Maxlevel,block_f,block_f%pgno,ptree,msh,0)
				call ComputeParallelIndices(ho_bf1%Maxlevel,block_f,block_f%pgno_db,ptree,msh,1)
				! if(block_f%M==2500)write(*,*)ptree%myID,block_f%pgno,block_f%pgno_db,block_f%N_loc,block_f%N_loc_db,'eref'
				
				block_f%style = 2
				allocate(ho_bf1%levels(level_c)%BP(ii)%LL(1)%boundary_map(1))
				ho_bf1%levels(level_c)%BP(ii)%LL(1)%boundary_map(1) = block_f%col_group
				ho_bf1%levels(level_c)%BP(ii)%Lplus=0
				
				
				group = floor((ii-1+2**level_c)/2d0)
				sortdirec = NINT(msh%basis_group(group)%boundary(1))
				seperator = msh%basis_group(group)%boundary(2)
				
				do ll=1,LplusMax-1
					if(ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound>0)then
						ho_bf1%levels(level_c)%BP(ii)%Lplus = ho_bf1%levels(level_c)%BP(ii)%Lplus + 1
						call assert(ho_bf1%levels(level_c)%BP(ii)%Lplus<=LplusMax,'increase LplusMax')
						! write(*,*)'nini',level_c,ho_bf1%Maxlevel - ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(1)%level,option%lnoBP,ll

						block_f => ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(1)
						
						if(ho_bf1%Maxlevel - block_f%level<option%lnoBP .or. ll==LplusMax-1 .or. block_f%level_butterfly==0)then
							ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%Nbound=0
						else
							! write(*,*)'gggggg'
							! level_butterfly = int((ho_bf1%Maxlevel - ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(1)%level)/2)*2 
							level_butterfly = block_f%level_butterfly
							level_BP = ho_bf1%levels(level_c)%BP(ii)%level
							levelm = ceiling_safe(dble(level_butterfly)/2d0)						
							groupm_start=block_f%row_group*2**levelm						
							Nboundall = 2**(block_f%level+levelm-level_BP)
							allocate(ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%boundary_map(Nboundall))
							ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%boundary_map = -1
							ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%Nbound=0
							
							do bb = 1,ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound
								blocks => ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(bb)
								
								allocate(Centroid_M(2**levelm,Dimn))
								allocate(Isboundary_M(2**levelm))
								Isboundary_M = 0
								Centroid_M = 0
								
								do index_i_m=1, 2**levelm
									group_m=blocks%row_group   ! Note: row_group and col_group interchanged here 
									group_m=group_m*2**levelm-1+index_i_m  										

									CNT = 0
									if(option%xyzsort==CKD)then	
										do nn = msh%basis_group(group_m)%head,msh%basis_group(group_m)%tail
											measure = abs(msh%xyz(sortdirec,msh%new2old(nn))-seperator)
											if(measure<option%touch_para)then
												Isboundary_M(index_i_m) = 1 
												CNT = CNT + 1
												Centroid_M(index_i_m,1:Dimn) = Centroid_M(index_i_m,1:Dimn)+ msh%xyz(1:Dimn,msh%new2old(nn))
											end if
										end do
										if(Isboundary_M(index_i_m)==1)Centroid_M(index_i_m,:) = Centroid_M(index_i_m,:)/CNT
										
										! if(blocks%col_group==8 .or. blocks%col_group==9)then
											! write(*,*)'wocaoo',group_m,Isboundary_M(index_i_m),CNT,sortdirec,seperator
										! endif

										write(*,*)'Bplus for other sorting not considered yet'
										stop
									end if									
								end do
														
								
								allocate(Centroid_N(2**(level_butterfly-levelm),Dimn))
								allocate(Isboundary_N(2**(level_butterfly-levelm)))
								Isboundary_N = 0
								Centroid_N = 0
								
								
								do index_j_m=1, 2**(level_butterfly-levelm)	
									group_n=blocks%col_group  
									group_n=group_n*2**(level_butterfly-levelm)-1+index_j_m
								
									CNT = 0
									if(option%xyzsort==CKD)then	
										do nn = msh%basis_group(group_n)%head,msh%basis_group(group_n)%tail
											measure = abs(msh%xyz(sortdirec,msh%new2old(nn))-seperator)
											if(measure<option%touch_para)then
												Isboundary_N(index_j_m) = 1 
												CNT = CNT + 1
												
												Centroid_N(index_j_m,1:Dimn) = Centroid_N(index_j_m,1:Dimn) + msh%xyz(1:Dimn,msh%new2old(nn))
											end if
										end do
										if(Isboundary_N(index_j_m)==1)Centroid_N(index_j_m,:) = Centroid_N(index_j_m,:)/CNT
										write(*,*)'Bplus for other sorting not considered yet'
										stop
									end if	
								end do								
								
								
								! if(level_c==1)then
								! ! write(*,*)Isboundary_N,Isboundary_M
								! do kk=1,2**levelm
									! if(Isboundary_M(kk)==1)then
									! ! write(*,*)Centroid_M(kk,1),Centroid_M(kk,2),Centroid_M(kk,3)
									! write(777,*)Centroid_M(kk,1),Centroid_M(kk,2),Centroid_M(kk,3)
									! end if
								! end do
								
								! do kk=1,2**(level_butterfly-levelm)
									! if(Isboundary_N(kk)==1)then
									! write(777,*)Centroid_N(kk,1),Centroid_N(kk,2),Centroid_N(kk,3)
									! end if
								! end do								
								! end if
								
								do index_i_m=1, 2**levelm
									group_m=blocks%row_group   ! Note: row_group and col_group interchanged here 
									group_m=group_m*2**levelm-1+index_i_m  								

									if(Isboundary_M(index_i_m)==1)then
										ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%Nbound = ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%Nbound + 1									
										dist = 100000000d0
										do index_j_m=1, 2**(level_butterfly-levelm)	
											group_n=blocks%col_group  
											group_n=group_n*2**(level_butterfly-levelm)-1+index_j_m
											
											! if(blocks%col_group==8 .or. blocks%col_group==9)then
												! write(*,*)group_m,group_n,sqrt(sum((Centroid_N(index_j_m,:)-Centroid_M(index_i_m,:))**2d0)),'nima'
											! end	if
											
											if(Isboundary_N(index_j_m)==1)then
												if(dist > sqrt(sum((Centroid_N(index_j_m,:)-Centroid_M(index_i_m,:))**2d0)))then
													! if(level_c==1)write(*,*)index_i_m,index_j_m
													dist = sqrt(sum((Centroid_N(index_j_m,:)-Centroid_M(index_i_m,:))**2d0))
													ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%boundary_map(group_m - groupm_start + 1) = group_n
												end if
											end if
										end do
									end if	
								enddo
								deallocate(Isboundary_M)
								deallocate(Isboundary_N)
								deallocate(Centroid_M)
								deallocate(Centroid_N)
								
							end do	
							
							if(ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%Nbound>1)then
								! write(*,*)level_c,ii,ll,ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%col_group,ho_bf1%levels(level_c)%BP(ii)%LL(1)%matrices_block(1)%row_group,ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%Nbound,'niamaa'
							endif
							
							
							call assert(ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%Nbound>0,'why is no boundary group detected')	
								
							allocate(ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%matrices_block(ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%Nbound))
							
							cnt = 0
							do bb = 1,	Nboundall
								if(ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%boundary_map(bb)/=-1)then
									cnt = cnt + 1
									group_m = bb+groupm_start-1
									group_n = ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%boundary_map(bb)
									blocks => ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%matrices_block(cnt)
									blocks%row_group = group_m
									blocks%col_group = group_n
									blocks%level = GetTreelevel(group_m)-1

									blocks%pgno = msh%basis_group(group_m)%pgno										
									blocks%M = msh%basis_group(group_m)%tail - msh%basis_group(group_m)%head + 1
									blocks%N = msh%basis_group(group_n)%tail - msh%basis_group(group_n)%head + 1
									blocks%headm = msh%basis_group(group_m)%head
									blocks%headn = msh%basis_group(group_n)%head						
									
									
									! blocks%level_butterfly = int((ho_bf1%Maxlevel - blocks%level)/2)*2
									blocks%level_butterfly = 0 ! only two layer butterfly plus here
									
									blocks%style = 2
									call ComputeParallelIndices(ho_bf1%Maxlevel,blocks,blocks%pgno,ptree,msh,0)
								end if
							end do
						end if		
					else 
						exit
					end if
				end do
				
				! write(*,*)level_c,ii,ho_bf1%levels(level_c)%BP(ii)%Lplus,'gaogao '
				
				if(mod(ii,2)==1)then  ! in the beginning only even block hold information about the schur complement
					ii_sch = ceiling_safe(ii/2d0)
					
					allocate(ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(LplusMax))
					ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%Lplus=ho_bf1%levels(level_c)%BP(ii)%Lplus
					do ll=1,LplusMax
						ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%Nbound=0
					end do	
									
					ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(1)%Nbound = 1
			
					allocate(ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(1)%boundary_map(1))
					ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(1)%boundary_map(1) = ho_bf1%levels(level_c)%BP(ii)%row_group		
					
								
					do ll=1,LplusMax-1
						if(ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound>0)then

							ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%rankmax = 0
							ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%Nbound = ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound
							
							allocate(ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%Nbound))
							
							do bb =1,ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%Nbound
								block_f => ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(bb)
								block_sch => ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(bb)
								
								row_group = block_f%row_group
								
								block_sch%row_group = row_group
								block_sch%col_group = row_group
								
								if(msh%basis_group(row_group)%pgno/=msh%basis_group(INT(row_group/2d0))%pgno)then
									block_sch%pgno = msh%basis_group(INT(row_group/2d0))%pgno								
								else 
									block_sch%pgno = msh%basis_group(row_group)%pgno
								end if
								
								
								block_sch%style = block_f%style
								block_sch%level = block_f%level
								block_sch%level_butterfly = block_f%level_butterfly
								
								block_sch%M = msh%basis_group(row_group)%tail - msh%basis_group(row_group)%head + 1
								block_sch%N = msh%basis_group(row_group)%tail - msh%basis_group(row_group)%head + 1
								block_sch%headm = msh%basis_group(row_group)%head 
								block_sch%headn = msh%basis_group(row_group)%head
								call ComputeParallelIndices(ho_bf1%Maxlevel,block_sch,block_sch%pgno,ptree,msh,0)
							end do
							
							
							if(ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%Nbound==0)then		
								ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll+1)%Nbound=0
							else 
								level_butterfly = ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(1)%level_butterfly
								level_BP = ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%level
								levelm = ceiling_safe(dble(level_butterfly)/2d0)						
								groupm_start=ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(1)%row_group*2**levelm		
								Nboundall = 2**(ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(1)%level+levelm-level_BP)				
								
								allocate(ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll+1)%boundary_map(Nboundall))
								
								ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll+1)%boundary_map = ho_bf1%levels(level_c)%BP(ii)%LL(ll+1)%boundary_map
								do bb=1,Nboundall
									if(ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll+1)%boundary_map(bb)/=-1)then
										ho_bf1%levels(level_c)%BP_inverse_schur(ii_sch)%LL(ll+1)%boundary_map(bb) = bb + groupm_start - 1
									end if
								end do
							end if
						else 
							exit
						end if
					end do
					
				end if			
								
				! ! if(level_c==1 .and. ii==1)then
				 
				! write(strings , *) 2*dimn
				! ! write(177,*)'Bplus:', level_c,ii
				! do ll=1,ho_bf1%levels(level_c)%BP(ii)%Lplus
				! ! write(*,*)ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound,'ddd'
					! do bb = 1,ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound
						! write(177,'(I3,I7,I3,I3,'//TRIM(strings)//'Es16.7)')level_c,ii,ll,ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%level,msh%basis_group(ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%row_group)%center(1:dimn),msh%basis_group(ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%col_group)%center(1:dimn)
					! end do
				! end do
				! ! end if
			end if
			! end if
			
			call copy_Bplus(ho_bf1%levels(level_c)%BP(ii),ho_bf1%levels(level_c)%BP_inverse_update(ii))		
			
		end do
	end do	

	

		
	msh%idxs = ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
	msh%idxe = ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)	
	
	
	
	
end subroutine BPlus_structuring


end module HODLR_structure 
