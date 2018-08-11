module H_structure 
use Utilities
contains 

real*8 function group_dist(group_m,group_n)

    use MODULE_FILE
    implicit none
    
    integer group_m, group_n,farblock, level
    integer i, j, ii, jj
    real*8 dis, rad1, rad2, para
	real*8,allocatable::a(:), b(:)
    integer Dimn
	Dimn = size(basis_group(group_m)%center,1)
    allocate(a(Dimn))
    allocate(b(Dimn))
	
    do i=1,Dimn
        a(i)=basis_group(group_m)%center(i)
        b(i)=basis_group(group_n)%center(i)
    enddo
        
    dis=0d0
    do i=1,Dimn
        dis=dis+(a(i)-b(i))**2
    enddo
    group_dist=sqrt(dis)
    
	deallocate(a,b)

end function group_dist


logical function near_or_far(group_m,group_n,para)

    use MODULE_FILE
    implicit none
    
    integer group_m, group_n,farblock, level
    integer i, j, ii, jj
    real*8 dis, rad1, rad2, para
    real*8,allocatable:: a(:), b(:)
    integer Dimn
	
	Dimn = size(basis_group(group_m)%center,1)
    allocate(a(Dimn))
    allocate(b(Dimn))
    do i=1,Dimn
        a(i)=basis_group(group_m)%center(i)
        b(i)=basis_group(group_n)%center(i)
    enddo
    
    rad1=basis_group(group_m)%radius
    rad2=basis_group(group_n)%radius
    
    dis=0d0
    do i=1,Dimn
        dis=dis+(a(i)-b(i))**2
    enddo
    dis=sqrt(dis)
      
		! write(*,*)dis/((rad1+rad2)/2)
	 
    ! if (dis>para*max(rad1,rad2)) then
    if (dis>para*(rad1+rad2)/2) then
        near_or_far=.true.
    else
        near_or_far=.false.
    endif
    deallocate(a,b)
	
end function

real*8 function func_distance(node1,node2,msh)
    
    use MODULE_FILE
    implicit none

    integer node1, node2
    real*8 dis
    integer i, j
    integer Dimn
	type(mesh)::msh
	
	Dimn = size(msh%xyz,1)
	
    dis=0d0
    do i=1,Dimn
        dis=dis+(msh%xyz(i,node1)-msh%xyz(i,node2))**2
    enddo
    
    func_distance=dis
    
    return
    
end function


subroutine H_matrix_structuring(ho_bf1,para,option,msh,ptree)
    
    use MODULE_FILE
	use misc
    implicit none
    
	integer Cflag
    integer i, j, ii, jj, iii, jjj
    integer level, edge, node, patch, group, group_m, group_n,col_group,row_group,fidx
    integer blocks
    integer center_edge
    
    integer index_temp
    
    real*8 a, b, c, d, para, xmax,xmin,ymax,ymin,zmax,zmin,seperator,r,theta,phi,phi_tmp
    real*8 radius, radiusmax
    real*8,allocatable:: xyzrange(:),xyzmin(:),xyzmax(:),auxpoint(:),groupcenter(:)
    real*8, allocatable :: distance(:),array(:,:)
	integer level_c,sortdirec,mm,phi_end,Ninfo_edge
	real*8 t1,t2
	integer Maxgroup
	character(len=1024)  :: strings	
    integer, allocatable :: order(:), edge_temp(:,:),map_temp(:)
	integer dimn
	type(Hoption)::option
	type(hobf)::ho_bf1
	integer Maxlevel,Maxgrp
	type(mesh)::msh
	type(proctree)::ptree
    !**********************Maxlevel*******************
    level=0; i=1
    do while (int(msh%Nunk/i)>option%Nmin_leaf)
        level=level+1
        i=2**level
    enddo
    
    Maxlevel=level
    ho_bf1%Maxlevel=Maxlevel
    !***************************************************
    
    Maxgroup=2**(Maxlevel+1)-1
    allocate (basis_group(Maxgroup))
    
	dimn=size(msh%xyz,1) 
	
	allocate(xyzrange(dimn))
	allocate(xyzmin(dimn))
	allocate(xyzmax(dimn))
	allocate(auxpoint(dimn))
	allocate(groupcenter(dimn))
		
    if(ptree%MyID==Main_ID)write (*,*) ''
    if(ptree%MyID==Main_ID)write (*,*) 'Maxlevel_for_blocks:',ho_bf1%Maxlevel
    if(ptree%MyID==Main_ID)write (*,*) 'N_leaf:',int(msh%Nunk/i)
 	if(ptree%MyID==Main_ID)write (*,*) ''
    if(ptree%MyID==Main_ID)write (*,*) 'Constructing basis groups...'

       
 
    !***************************************************************************************	   
	allocate(msh%new2old(msh%Nunk))    

	do ii=1,msh%Nunk
		msh%new2old(ii) = ii
	end do
	
	! write(*,*)'gan', msh%info_unk(0,100)   
	
	   
    ! allocate (distance(msh%Nunk))     
    
    !********************************index_of_group**************************************
    
	Maxgrp=2**(ptree%nlevel)-1
    basis_group(1)%head=1 ; basis_group(1)%tail=msh%Nunk; basis_group(1)%pgno=1
    do level=0, Maxlevel
        do group=2**level, 2**(level+1)-1
            basis_group(group)%level=level

            groupcenter(1:dimn)=0.0d0
            ! !$omp parallel do default(shared) private(edge,ii) reduction(+:groupcenter)
            do edge=basis_group(group)%head, basis_group(group)%tail
                do ii=1,dimn
					! write(*,*)edge,msh%info_unk(0,edge)
                    groupcenter(ii)=groupcenter(ii)+msh%xyz(ii,msh%info_unk(0,edge))
                enddo
				! if(group==24)write(*,*)edge,groupcenter(:),msh%xyz(1,msh%info_unk(0,edge))
			enddo
            ! !$omp end parallel do
            do ii=1,dimn
                groupcenter(ii)=groupcenter(ii)/(basis_group(group)%tail-basis_group(group)%head+1)
            enddo
            allocate(basis_group(group)%center(dimn))
			basis_group(group)%center(1:dimn)=groupcenter(1:dimn)

			radiusmax=0.
			do edge=basis_group(group)%head, basis_group(group)%tail
				radius=0
				do ii=1,dimn
					radius=radius+(msh%xyz(ii,msh%info_unk(0,edge))-groupcenter(ii))**2
				enddo
				! write(*,*)'really',edge, msh%info_unk(0,edge)
				radius=sqrt(radius)
				if (radius>radiusmax) then
					radiusmax=radius
					center_edge=edge
				endif			
			enddo
			basis_group(group)%radius=radiusmax		
				
			if(size(basis_group_pre,1)<2*group+1)then ! if groups not predefined, need to order the points
				if(option%xyzsort==1)then !msh%xyz sort		 
					xyzmin= 1d300
					xyzmax= -1d300
					do edge=basis_group(group)%head, basis_group(group)%tail
						do ii=1,Dimn
							xyzmax(ii) = max(xyzmax(ii),msh%xyz(ii,msh%info_unk(0,edge)))
							xyzmin(ii) = min(xyzmin(ii),msh%xyz(ii,msh%info_unk(0,edge)))
						enddo
					enddo
					xyzrange(1:Dimn) = xyzmax(1:Dimn)-xyzmin(1:Dimn)
					
					mm = basis_group(group)%tail - basis_group(group)%head + 1
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
					do i=basis_group(group)%head, basis_group(group)%tail
						distance(i-basis_group(group)%head+1)=msh%xyz(sortdirec,msh%info_unk(0,i))
					enddo
					!$omp end parallel do

				else if(option%xyzsort==2)then ! spherical sort
					phi_end = 0
					do edge=basis_group(group)%head, basis_group(group)%tail
						call Cart2Sph(msh%xyz(1,msh%info_unk(0,edge)),msh%xyz(2,msh%info_unk(0,edge)),msh%xyz(3,msh%info_unk(0,edge)),msh%Origins,r,theta,phi)					
						if(abs(phi-2*pi)< option%touch_para*msh%minedgelength/r)phi_end=1      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! this group contains points near phi=2pi
					enddo
					mm = basis_group(group)%tail - basis_group(group)%head + 1
					allocate (distance(mm))	
				
					if(mod(level,2)==1)then           !!!!!!!!!!!!!!!!!!!!!!!!! note: applys only to spheres sort Phi first
						sortdirec=1
					else 
						sortdirec=2
					end if
					
					!$omp parallel do default(shared) private(i,theta,phi,r,phi_tmp)
					do i=basis_group(group)%head, basis_group(group)%tail
						call Cart2Sph(msh%xyz(1,msh%info_unk(0,i)),msh%xyz(2,msh%info_unk(0,i)),msh%xyz(3,msh%info_unk(0,i)),msh%Origins,r,theta,phi)
						if(sortdirec==1)distance(i-basis_group(group)%head+1)=theta
						if(sortdirec==2)then
							distance(i-basis_group(group)%head+1)=phi
							if(phi_end==1)then   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! if the group contains points near phi=2pi, substract phi by pi
								phi_tmp = phi-pi
								if(phi_tmp<=0) phi_tmp = phi_tmp + 2*pi
								distance(i-basis_group(group)%head+1)=phi_tmp
							end if
						end if
					enddo
					!$omp end parallel do
					
				else if(option%xyzsort==3)then !cobblestone sort
					
					mm = basis_group(group)%tail - basis_group(group)%head + 1
					allocate (distance(mm))	
					
					distance(1:mm)=Bigvalue
					!$omp parallel do default(shared) private(i)
					do i=basis_group(group)%head, basis_group(group)%tail
						distance(i-basis_group(group)%head+1)=func_distance(msh%info_unk(0,i),msh%info_unk(0,center_edge),msh)
					enddo
					!$omp end parallel do					

				end if
				
				Ninfo_edge=size(msh%info_unk,1)-1
				allocate (order(mm),edge_temp(0:Ninfo_edge,mm)) 
				allocate(map_temp(mm))

				call quick_sort(distance,order,mm)       
				!$omp parallel do default(shared) private(ii)     
				do ii=1, mm
					edge_temp(:,ii)=msh%info_unk(:,order(ii)+basis_group(group)%head-1)
					map_temp(ii) = msh%new2old(order(ii)+basis_group(group)%head-1)
				enddo
				!$omp end parallel do

				!$omp parallel do default(shared) private(ii)     
				do ii=1, mm
					msh%info_unk(:,ii+basis_group(group)%head-1)=edge_temp(:,ii)
					msh%new2old(ii+basis_group(group)%head-1) = map_temp(ii)
				enddo
				!$omp end parallel do			
				

				deallocate(edge_temp)
				deallocate(map_temp)
				deallocate(order)

				deallocate(distance)
				
				
				if (level<Maxlevel) then
					basis_group(2*group)%head=basis_group(group)%head
					basis_group(2*group)%tail=int((basis_group(group)%head+basis_group(group)%tail)/2)
					basis_group(2*group+1)%head=basis_group(2*group)%tail+1
					basis_group(2*group+1)%tail=basis_group(group)%tail
					
					if(basis_group(group)%pgno*2<=Maxgrp)then
						basis_group(2*group)%pgno=basis_group(group)%pgno*2
					else
						basis_group(2*group)%pgno=basis_group(group)%pgno
					endif
					if(basis_group(group)%pgno*2+1<=Maxgrp)then
						basis_group(2*group+1)%pgno=basis_group(group)%pgno*2+1
					else
						basis_group(2*group+1)%pgno=basis_group(group)%pgno
					endif					
					
					if(option%xyzsort==1)then 
						seperator = msh%xyz(sortdirec,msh%info_unk(0,basis_group(2*group)%tail))
						
					else if(option%xyzsort==2)then  
						call Cart2Sph(msh%xyz(1,msh%info_unk(0,basis_group(2*group)%tail)),msh%xyz(2,msh%info_unk(0,basis_group(2*group)%tail)),msh%xyz(3,msh%info_unk(0,basis_group(2*group)%tail)),msh%Origins,r,theta,phi)
						if(sortdirec==1)seperator = theta
						if(sortdirec==2)then
							seperator = phi
							! write(*,*)level,phi*180/pi,basis_group(2*group)%tail,'ganni',phi_end
						end if
					end if
					
					
					basis_group(group)%boundary(1) = sortdirec
					basis_group(group)%boundary(2) = seperator
					
					! fidx = (2*group)- 2**(level+1)+1
					! ho_bf1%levels(level+1)%BP(fidx)%boundary(1) = sortdirec
					! ho_bf1%levels(level+1)%BP(fidx)%boundary(2) = seperator
					
					! fidx = (2*group+1)- 2**(level+1)+1
					! ho_bf1%levels(level+1)%BP(fidx)%boundary(1) = sortdirec
					! ho_bf1%levels(level+1)%BP(fidx)%boundary(2) = seperator				 					
				endif	

			else 
				if (level<Maxlevel) then
					basis_group(2*group)%head=basis_group_pre(2*group,1)
					basis_group(2*group)%tail=basis_group_pre(2*group,2)
					basis_group(2*group+1)%head=basis_group_pre(2*group+1,1)
					basis_group(2*group+1)%tail=basis_group_pre(2*group+1,2)
					
					basis_group(group)%boundary(1) = 0
					basis_group(group)%boundary(2) = 0
					
					! fidx = (2*group)- 2**(level+1)+1
					! ho_bf1%levels(level+1)%BP(fidx)%boundary(1) = 0   ! dummy parameters
					! ho_bf1%levels(level+1)%BP(fidx)%boundary(2) = 0
					
					! fidx = (2*group+1)- 2**(level+1)+1
					! ho_bf1%levels(level+1)%BP(fidx)%boundary(1) = 0
					! ho_bf1%levels(level+1)%BP(fidx)%boundary(2) = 0				 
				endif					
			end if
        enddo
    enddo
	
	
	
	! write(strings , *) Dimn
	! do level=0, Maxlevel
        ! do group=2**level, 2**(level+1)-1
            ! do edge=basis_group(group)%head, basis_group(group)%tail
				! ! write(*,*)edge,msh%info_unk(0,edge)
				! ! write(113,'(I5,I8,Es16.8,Es16.8,Es16.8)')level,group,msh%xyz(1:Dimn,msh%info_unk(0,edge))
				! write(113,'(I5,I8,'//TRIM(strings)//'Es16.8)')level,group,msh%xyz(1:Dimn,msh%info_unk(0,edge))
			! enddo
        ! enddo
    ! enddo	   
    !*************************************************************************************
   
	
	deallocate(xyzrange)
	deallocate(xyzmin)
	deallocate(xyzmax)
	deallocate(auxpoint)
	deallocate(groupcenter)	
    
	
	allocate(msh%old2new(msh%Nunk))
	do ii=1,msh%Nunk
		msh%old2new(msh%new2old(ii)) = ii
	end do	
	
    return
    
end subroutine H_matrix_structuring


subroutine BPlus_structuring(ho_bf1,option,msh,ptree)
	use MODULE_FILE
	use misc
	implicit none 
	
    integer i, j, ii, jj, kk, iii, jjj,ll,bb,sortdirec,ii_sch
    integer level, edge, patch, node, group, group_touch
    integer rank, index_near, m, n, length, flag, itemp,cnt,detection
    real T0
	real*8:: tolerance, rtemp,rel_error,seperator,dist
    real*8 Memory_direct_forward,Memory_butterfly_forward
	integer mm,nn,header_m,header_n,edge_m,edge_n,group_m,group_n,group_m1,group_n1,group_m2,group_n2,levelm,groupm_start,index_i_m,index_j_m
	complex(kind=8)::ctemp,ctemp1,ctemp2
	integer level_c,iter,level_cc,level_BP,Nboundall,level_butterfly	
	type(matrixblock),pointer::blocks,block_f,block_sch,block_inv
	real*8::minbound,theta,phi,r,rmax,phi_tmp,measure
	real*8,allocatable::Centroid_M(:,:),Centroid_N(:,:)
	integer,allocatable::Isboundary_M(:),Isboundary_N(:)
	integer Dimn,col_group,row_group,Maxgrp
	type(Hoption)::option
	type(mesh)::msh
	type(hobf)::ho_bf1
	character(len=1024)  :: strings
	type(proctree)::ptree

	Maxgrp=2**(ptree%nlevel)-1
	
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
		allocate(ho_bf1%levels(level_c)%BP_inverse_schur(ho_bf1%levels(level_c)%N_block_inverse))
	end do

	
	ho_bf1%levels(1)%BP_inverse(1)%level = 0
	ho_bf1%levels(1)%BP_inverse(1)%col_group = 1
	ho_bf1%levels(1)%BP_inverse(1)%row_group = 1
	ho_bf1%levels(1)%BP_inverse(1)%pgno = 1
	! ho_bf1%levels(1)%BP_inverse(1)%style = 2
	
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
			block_inv%headm=basis_group(row_group)%head
			block_inv%headn=basis_group(col_group)%head
			block_inv%M=basis_group(row_group)%tail-basis_group(row_group)%head+1
			block_inv%N=basis_group(col_group)%tail-basis_group(col_group)%head+1
			
			block_inv%level_butterfly = ho_bf1%Maxlevel - block_inv%level
			
			call ComputeParallelIndices(ho_bf1%Maxlevel,block_inv,block_inv%pgno,ptree,0)

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
				if(GetTreelevel(ho_bf1%levels(level_c)%BP_inverse(ii)%pgno)==ptree%nlevel)then
					ho_bf1%levels(level_c+1)%BP(ii*2-1)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
					ho_bf1%levels(level_c+1)%BP(ii*2)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno
				else
					ho_bf1%levels(level_c+1)%BP(ii*2-1)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno*2
					ho_bf1%levels(level_c+1)%BP(ii*2)%pgno = ho_bf1%levels(level_c)%BP_inverse(ii)%pgno*2+1					
				endif
				ho_bf1%levels(level_c+1)%BP(ii*2)%level = level_c+1
				ho_bf1%levels(level_c+1)%BP(ii*2)%col_group = col_group*2+1
				ho_bf1%levels(level_c+1)%BP(ii*2)%row_group = row_group*2+1
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
	
	
	Dimn = size(msh%xyz,1)
	
	do level_c = 1,ho_bf1%Maxlevel+1
		do ii =1,ho_bf1%levels(level_c)%N_block_forward
            ! if(ptree%MyID >=ptree%pgrp(ho_bf1%levels(level_c)%BP(ii)%pgno)%head .and. ptree%MyID <=ptree%pgrp(ho_bf1%levels(level_c)%BP(ii)%pgno)%tail)then
			
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
				block_f%pgno = basis_group(block_f%row_group)%pgno
				
				block_f%M = basis_group(block_f%row_group)%tail - basis_group(block_f%row_group)%head + 1
				block_f%N = basis_group(block_f%col_group)%tail - basis_group(block_f%col_group)%head + 1
				block_f%headm = basis_group(block_f%row_group)%head
				block_f%headn = basis_group(block_f%col_group)%head
				block_f%level_butterfly=0
				call ComputeParallelIndices(ho_bf1%Maxlevel+1,block_f,block_f%pgno,ptree,0)
				
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
				block_inv%pgno = basis_group(block_inv%row_group)%pgno
				
				block_inv%M = basis_group(block_inv%row_group)%tail - basis_group(block_inv%row_group)%head + 1
				block_inv%N = basis_group(block_inv%col_group)%tail - basis_group(block_inv%col_group)%head + 1
				block_inv%headm = basis_group(block_inv%row_group)%head
				block_inv%headn = basis_group(block_inv%col_group)%head
				block_inv%level_butterfly=0
				call ComputeParallelIndices(ho_bf1%Maxlevel+1,block_inv,block_inv%pgno,ptree,0)			
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
					if(ho_bf1%Maxlevel - block_f%level<option%LnoBP)then				
						block_f%level_butterfly = ho_bf1%Maxlevel - ho_bf1%levels(level_c)%BP(ii)%level   ! butterfly 
					else
						block_f%level_butterfly = int((ho_bf1%Maxlevel - ho_bf1%levels(level_c)%BP(ii)%level)/2)*2 ! butterfly plus needs even number of levels				
					endif
				endif
				
				block_f%col_group = ho_bf1%levels(level_c)%BP(ii)%col_group
				block_f%row_group = ho_bf1%levels(level_c)%BP(ii)%row_group
				block_f%pgno=basis_group(block_f%row_group)%pgno

				
				
				! compute the partial indices when BP is shared by double number of processes
				ii_sch = ceiling_safe(ii/2d0)
				block_inv => ho_bf1%levels(level_c)%BP_inverse(ii_sch)%LL(1)%matrices_block(1)
				block_f%pgno_db = block_inv%pgno	 	
	
						
		
				
				block_f%M = basis_group(block_f%row_group)%tail - basis_group(block_f%row_group)%head + 1
				block_f%N = basis_group(block_f%col_group)%tail - basis_group(block_f%col_group)%head + 1
				block_f%headm = basis_group(block_f%row_group)%head
				block_f%headn = basis_group(block_f%col_group)%head				
				
				call ComputeParallelIndices(ho_bf1%Maxlevel,block_f,block_f%pgno,ptree,0)
				call ComputeParallelIndices(ho_bf1%Maxlevel,block_f,block_f%pgno_db,ptree,1)
				! if(block_f%M==2500)write(*,*)ptree%myID,block_f%pgno,block_f%pgno_db,block_f%N_loc,block_f%N_loc_db,'eref'
				
				block_f%style = 2
				allocate(ho_bf1%levels(level_c)%BP(ii)%LL(1)%boundary_map(1))
				ho_bf1%levels(level_c)%BP(ii)%LL(1)%boundary_map(1) = block_f%col_group
				ho_bf1%levels(level_c)%BP(ii)%Lplus=0
				
				
				group = floor((ii-1+2**level_c)/2d0)
				sortdirec = NINT(basis_group(group)%boundary(1))
				seperator = basis_group(group)%boundary(2)
				
				do ll=1,LplusMax-1
					if(ho_bf1%levels(level_c)%BP(ii)%LL(ll)%Nbound>0)then
						ho_bf1%levels(level_c)%BP(ii)%Lplus = ho_bf1%levels(level_c)%BP(ii)%Lplus + 1
						call assert(ho_bf1%levels(level_c)%BP(ii)%Lplus<=LplusMax,'increase LplusMax')
						! write(*,*)'nini',level_c,ho_bf1%Maxlevel - ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(1)%level,option%LnoBP,ll

						block_f => ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(1)
						
						if(ho_bf1%Maxlevel - block_f%level<option%LnoBP .or. ll==LplusMax-1 .or. block_f%level_butterfly==0)then
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
									if(option%xyzsort==1)then	
										do nn = basis_group(group_m)%head,basis_group(group_m)%tail
											measure = abs(msh%xyz(sortdirec,msh%info_unk(0,nn))-seperator)
											if(measure<option%touch_para*msh%minedgelength)then
												Isboundary_M(index_i_m) = 1 
												CNT = CNT + 1
												Centroid_M(index_i_m,1:Dimn) = Centroid_M(index_i_m,1:Dimn)+ msh%xyz(1:Dimn,msh%info_unk(0,nn))
											end if
										end do
										if(Isboundary_M(index_i_m)==1)Centroid_M(index_i_m,:) = Centroid_M(index_i_m,:)/CNT
										
										! if(blocks%col_group==8 .or. blocks%col_group==9)then
											! write(*,*)'wocaoo',group_m,Isboundary_M(index_i_m),CNT,sortdirec,seperator
										! endif
										
										
									else if(option%xyzsort==2)then
										do nn = basis_group(group_m)%head,basis_group(group_m)%tail
											call Cart2Sph(msh%xyz(1,msh%info_unk(0,nn)),msh%xyz(2,msh%info_unk(0,nn)),msh%xyz(3,msh%info_unk(0,nn)),msh%Origins,r,theta,phi)
											if(sortdirec==1)then
												measure = abs(theta-seperator)
												if(measure< option%touch_para*msh%minedgelength/r)then
													Isboundary_M(index_i_m) = 1 
													CNT = CNT + 1
													Centroid_M(index_i_m,1) = Centroid_M(index_i_m,1) + msh%xyz(1,msh%info_unk(0,nn))
													Centroid_M(index_i_m,2) = Centroid_M(index_i_m,2) + msh%xyz(2,msh%info_unk(0,nn))
													Centroid_M(index_i_m,3) = Centroid_M(index_i_m,3) + msh%xyz(3,msh%info_unk(0,nn))													
												end if
											end if
											if(sortdirec==2)then
												measure = abs(phi-seperator)
												measure = min(measure, abs(phi-(seperator+2*pi)))  !!! handle phi end point here
												measure = min(measure, abs(phi-(seperator-2*pi)))  !!! handle phi end point here
												if(level_c==1)then	
													measure = min(measure, abs(phi-pi))     !!! treat phi=pi as extra seperators for the first level groups
													! minbound = min(minbound, abs(phi-2*pi))  
												end if
												if(measure< option%touch_para*msh%minedgelength/r)then
													Isboundary_M(index_i_m) = 1 
													CNT = CNT + 1
													Centroid_M(index_i_m,1) = Centroid_M(index_i_m,1) + msh%xyz(1,msh%info_unk(0,nn))
													Centroid_M(index_i_m,2) = Centroid_M(index_i_m,2) + msh%xyz(2,msh%info_unk(0,nn))
													Centroid_M(index_i_m,3) = Centroid_M(index_i_m,3) + msh%xyz(3,msh%info_unk(0,nn))	
													! if(level_c==1 .and. index_i_m==1)write(*,*)CNT, msh%xyz(:,msh%info_unk(0,nn))												
												end if												
											end if
										end do										
										! if(level_c==1 .and. index_i_m==1)write(*,*)Centroid_M(index_i_m,:),CNT,'dddddd'
										
										if(Isboundary_M(index_i_m)==1)Centroid_M(index_i_m,:) = Centroid_M(index_i_m,:)/CNT
										! if(level_c==1)write(*,*)Centroid_M(1,:),CNT,'dddddd'
									else 
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
									if(option%xyzsort==1)then	
										do nn = basis_group(group_n)%head,basis_group(group_n)%tail
											measure = abs(msh%xyz(sortdirec,msh%info_unk(0,nn))-seperator)
											if(measure<option%touch_para*msh%minedgelength)then
												Isboundary_N(index_j_m) = 1 
												CNT = CNT + 1
												! write(*,*)nn,index_j_m,'ok'
												! write(*,*)Centroid_N(index_j_m,1),msh%xyz(1,msh%info_unk(0,nn))
												Centroid_N(index_j_m,1:Dimn) = Centroid_N(index_j_m,1:Dimn) + msh%xyz(1:Dimn,msh%info_unk(0,nn))
											end if
										end do
										if(Isboundary_N(index_j_m)==1)Centroid_N(index_j_m,:) = Centroid_N(index_j_m,:)/CNT
										
									else if(option%xyzsort==2)then
										do nn = basis_group(group_n)%head,basis_group(group_n)%tail
											call Cart2Sph(msh%xyz(1,msh%info_unk(0,nn)),msh%xyz(2,msh%info_unk(0,nn)),msh%xyz(3,msh%info_unk(0,nn)),msh%Origins,r,theta,phi)
											if(sortdirec==1)then
												measure = abs(theta-seperator)
												if(measure< option%touch_para*msh%minedgelength/r)then
													Isboundary_N(index_j_m) = 1 
													CNT = CNT + 1
													Centroid_N(index_j_m,1) = Centroid_N(index_j_m,1) + msh%xyz(1,msh%info_unk(0,nn))
													Centroid_N(index_j_m,2) = Centroid_N(index_j_m,2) + msh%xyz(2,msh%info_unk(0,nn))
													Centroid_N(index_j_m,3) = Centroid_N(index_j_m,3) + msh%xyz(3,msh%info_unk(0,nn))													
												end if
											end if
											if(sortdirec==2)then
												measure = abs(phi-seperator)
												measure = min(measure, abs(phi-(seperator+2*pi)))  !!! handle phi end point here
												measure = min(measure, abs(phi-(seperator-2*pi)))  !!! handle phi end point here
												if(level_c==1)then	
													measure = min(measure, abs(phi-pi))     !!! treat phi=pi as extra seperators for the first level groups
													! minbound = min(minbound, abs(phi-2*pi))  
												end if
												if(measure< option%touch_para*msh%minedgelength/r)then
													Isboundary_N(index_j_m) = 1 
													CNT = CNT + 1
													Centroid_N(index_j_m,1) = Centroid_N(index_j_m,1) + msh%xyz(1,msh%info_unk(0,nn))
													Centroid_N(index_j_m,2) = Centroid_N(index_j_m,2) + msh%xyz(2,msh%info_unk(0,nn))
													Centroid_N(index_j_m,3) = Centroid_N(index_j_m,3) + msh%xyz(3,msh%info_unk(0,nn))													
												end if												
											end if
										end do										
										if(Isboundary_N(index_j_m)==1)Centroid_N(index_j_m,:) = Centroid_N(index_j_m,:)/CNT
									else 
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
									blocks%level = basis_group(group_m)%level

									blocks%pgno = basis_group(group_m)%pgno										
									blocks%M = basis_group(group_m)%tail - basis_group(group_m)%head + 1
									blocks%N = basis_group(group_n)%tail - basis_group(group_n)%head + 1
									blocks%headm = basis_group(group_m)%head
									blocks%headn = basis_group(group_n)%head						
									
									
									! blocks%level_butterfly = int((ho_bf1%Maxlevel - blocks%level)/2)*2
									blocks%level_butterfly = 0 ! only two layer butterfly plus here
									
									blocks%style = 2
									call ComputeParallelIndices(ho_bf1%Maxlevel,blocks,blocks%pgno,ptree,0)
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
								
								if(basis_group(row_group)%pgno/=basis_group(INT(row_group/2d0))%pgno)then
									block_sch%pgno = basis_group(INT(row_group/2d0))%pgno								
								else 
									block_sch%pgno = basis_group(row_group)%pgno
								end if
								
								
								block_sch%style = block_f%style
								block_sch%level = block_f%level
								block_sch%level_butterfly = block_f%level_butterfly
								
								block_sch%M = basis_group(row_group)%tail - basis_group(row_group)%head + 1
								block_sch%N = basis_group(row_group)%tail - basis_group(row_group)%head + 1
								block_sch%headm = basis_group(row_group)%head 
								block_sch%headn = basis_group(row_group)%head
								call ComputeParallelIndices(ho_bf1%Maxlevel,block_sch,block_sch%pgno,ptree,0)
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
						! write(177,'(I3,I7,I3,I3,'//TRIM(strings)//'Es16.7)')level_c,ii,ll,ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%level,basis_group(ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%row_group)%center(1:dimn),basis_group(ho_bf1%levels(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%col_group)%center(1:dimn)
					! end do
				! end do
				! ! end if
			end if
			! end if
		end do
	end do	

	

		
	msh%idxs = ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
	msh%idxe = ho_bf1%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)	
	
	
	
	
end subroutine BPlus_structuring


end module H_structure 