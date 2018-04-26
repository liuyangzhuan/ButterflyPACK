module H_structure 

contains 

real*8 function group_dist(group_m,group_n)

    use MODULE_FILE
    implicit none
    
    integer group_m, group_n,farblock, level
    integer i, j, ii, jj
    real*8 a(MaxDim), b(MaxDim), dis, rad1, rad2, para
    integer Dimn
	Dimn = size(xyz,1)
    
    do i=1,Dimn
        a(i)=basis_group(group_m)%center(i)
        b(i)=basis_group(group_n)%center(i)
    enddo
        
    dis=0d0
    do i=1,Dimn
        dis=dis+(a(i)-b(i))**2
    enddo
    group_dist=sqrt(dis)
      

end function group_dist


logical function near_or_far(group_m,group_n,para)

    use MODULE_FILE
    implicit none
    
    integer group_m, group_n,farblock, level
    integer i, j, ii, jj
    real*8 a(MaxDim), b(MaxDim), dis, rad1, rad2, para
    integer Dimn
	Dimn = size(xyz,1)   
    
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
    
end function

real*8 function func_distance(node1,node2)
    
    use MODULE_FILE
    implicit none

    integer node1, node2
    real*8 dis
    integer i, j
    integer Dimn
	Dimn = size(xyz,1)
	
    dis=0d0
    do i=1,Dimn
        dis=dis+(xyz(i,node1)-xyz(i,node2))**2
    enddo
    
    func_distance=dis
    
    return
    
end function

integer function rank_approximate_func(group_m, group_n, flag)

    use MODULE_FILE
    implicit none

    integer i, j, k, mm, nn, edge_head, edge_tail, rank, group_m, group_n, flag
    real*8 a, b, c, aa(2), bb(2), cc(2), angle, distance

	if(group_m/=group_n)then
		
		distance=group_dist(group_m,group_n)
		distance=distance**2d0
		! distance=(basis_group(group_m)%center(1)-basis_group(group_n)%center(1))**2+(basis_group(group_m)%center(2)-basis_group(group_n)%center(2))**2+(basis_group(group_m)%center(3)-basis_group(group_n)%center(3))**2
		! ! distance=sqrt(distance)
		angle=4*pi*(basis_group(group_m)%radius)**2/distance
		rank=int(4*pi*(basis_group(group_n)%radius)**2*angle/wavelength**2)+1
		! if(group_m==4 .and. group_n==24)write(*,*)int(rank*rank_approximate_para1),rank,basis_group(group_n)%radius,basis_group(group_n)%radius,angle,distance
		if (flag==1) then
			rank_approximate_func=int(rank*rank_approximate_para1**2)
		elseif (flag==2) then
			rank_approximate_func=int(rank*rank_approximate_para2**2)
		elseif (flag==3) then
			rank_approximate_func=int(rank*rank_approximate_para3**2)
		endif
	else 
		rank_approximate_func = 100000
	end if
    !if (rank==0) then
    !    pause
    !    continue
    !endif

    return

end function rank_approximate_func





subroutine H_matrix_structuring(para)
    
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
    real*8 groupcenter(MaxDim), radius, radiusmax,auxpoint(MaxDim)
    real*8:: xyzrange(MaxDim),xyzmin(MaxDim),xyzmax(MaxDim)
    real*8, allocatable :: distance(:),array(:,:)
	integer level_c,sortdirec,mm,phi_end,Ninfo_edge
	real*8 t1,t2

    integer, allocatable :: order(:), edge_temp(:,:),map_temp(:)
	integer dimn
  
    !**********************Maxlevel*******************
    level=0; i=1
    do while (int(Maxedge/i)>Nmin_leaf)
        level=level+1
        i=2**level
    enddo
    
    Maxlevel=level
    Maxlevel_for_blocks=Maxlevel-Refined_level
    !***************************************************
    
    Maxgroup=2**(Maxlevel+1)-1
    allocate (basis_group(Maxgroup))
    
	dimn=size(xyz,1) 
	if(dimn>MaxDim)then
	! if(MyID==Main_ID)then
		write(*,*)'increase MaxDim to: ', dimn
		stop
	! endif
	endif    
	
	
    write (*,*) ''
    write (*,*) 'Maxlevel_for_blocks:',Maxlevel_for_blocks
    write (*,*) 'N_leaf:',int(Maxedge/i)
 	write (*,*) ''
    write (*,*) 'Constructing basis groups...'
    open (256,file='Info.txt',position='append')
    write (256,*) 'Maxlevel_for_blocks:',Maxlevel_for_blocks
    close (256) 
       
 
    !***************************************************************************************	   
    
    Maxblock=(4**Maxlevel_for_blocks-1)*4/3+1
	allocate(cascading_factors(Maxlevel_for_blocks+1))

	do level_c = 1,Maxlevel_for_blocks
		cascading_factors(level_c)%level = level_c
		cascading_factors(level_c)%N_block_forward = 2**level_c
		cascading_factors(level_c)%N_block_inverse = 2**(level_c-1)
		allocate(cascading_factors(level_c)%matrices_block(cascading_factors(level_c)%N_block_forward))
		allocate(cascading_factors(level_c)%matrices_block_inverse(cascading_factors(level_c)%N_block_inverse))
		allocate(cascading_factors(level_c)%matrices_block_inverse_schur(cascading_factors(level_c)%N_block_inverse))
		
		allocate(cascading_factors(level_c)%BP(cascading_factors(level_c)%N_block_forward))		
		allocate(cascading_factors(level_c)%BP_inverse(cascading_factors(level_c)%N_block_inverse))
		allocate(cascading_factors(level_c)%BP_inverse_schur(cascading_factors(level_c)%N_block_inverse))
		
		
	end do
	level_c = Maxlevel_for_blocks+1
	cascading_factors(level_c)%level = level_c
	cascading_factors(level_c)%N_block_forward = 2**(level_c-1)
	cascading_factors(level_c)%N_block_inverse = 2**(level_c-1)	
	allocate(cascading_factors(level_c)%matrices_block(cascading_factors(level_c)%N_block_forward))
	allocate(cascading_factors(level_c)%matrices_block_inverse(cascading_factors(level_c)%N_block_inverse))
	
	allocate(cascading_factors(level_c)%BP(cascading_factors(level_c)%N_block_forward))	
	allocate(cascading_factors(level_c)%BP_inverse(cascading_factors(level_c)%N_block_inverse))
	
	cascading_factors(1)%BP_inverse(1)%level = 0
	cascading_factors(1)%BP_inverse(1)%col_group = 1
	cascading_factors(1)%BP_inverse(1)%row_group = 1
	! cascading_factors(1)%BP_inverse(1)%style = 2
	
	do level_c = 1,Maxlevel_for_blocks
		do ii = 1, cascading_factors(level_c)%N_block_inverse
			col_group = cascading_factors(level_c)%BP_inverse(ii)%col_group
			row_group = cascading_factors(level_c)%BP_inverse(ii)%row_group
			
			cascading_factors(level_c)%matrices_block(ii*2-1)%level = level_c
			cascading_factors(level_c)%matrices_block(ii*2-1)%col_group = col_group*2+1
			cascading_factors(level_c)%matrices_block(ii*2-1)%row_group = row_group*2
			cascading_factors(level_c)%matrices_block(ii*2-1)%style = 2			
			cascading_factors(level_c)%matrices_block(ii*2)%level = level_c
			cascading_factors(level_c)%matrices_block(ii*2)%col_group = col_group*2
			cascading_factors(level_c)%matrices_block(ii*2)%row_group = row_group*2+1
			cascading_factors(level_c)%matrices_block(ii*2)%style = 2		

			cascading_factors(level_c)%BP(ii*2-1)%level = level_c
			cascading_factors(level_c)%BP(ii*2-1)%col_group = col_group*2+1
			cascading_factors(level_c)%BP(ii*2-1)%row_group = row_group*2
		
			cascading_factors(level_c)%BP(ii*2)%level = level_c
			cascading_factors(level_c)%BP(ii*2)%col_group = col_group*2
			cascading_factors(level_c)%BP(ii*2)%row_group = row_group*2+1
			
			

			
			
			if(level_c/=Maxlevel_for_blocks)then
			
				cascading_factors(level_c+1)%BP_inverse(ii*2-1)%level = level_c
				cascading_factors(level_c+1)%BP_inverse(ii*2-1)%col_group = col_group*2
				cascading_factors(level_c+1)%BP_inverse(ii*2-1)%row_group = row_group*2
				! cascading_factors(level_c+1)%BP_inverse(ii*2-1)%style = 2
				cascading_factors(level_c+1)%BP_inverse(ii*2)%level = level_c
				cascading_factors(level_c+1)%BP_inverse(ii*2)%col_group = col_group*2+1
				cascading_factors(level_c+1)%BP_inverse(ii*2)%row_group = row_group*2+1
				! cascading_factors(level_c+1)%BP_inverse(ii*2)%style = 2	

			else 

				
				cascading_factors(level_c+1)%matrices_block(ii*2-1)%level = level_c+1
				cascading_factors(level_c+1)%matrices_block(ii*2-1)%col_group = col_group*2
				cascading_factors(level_c+1)%matrices_block(ii*2-1)%row_group = row_group*2
				cascading_factors(level_c+1)%matrices_block(ii*2-1)%style = 1
				cascading_factors(level_c+1)%matrices_block(ii*2)%level = level_c+1
				cascading_factors(level_c+1)%matrices_block(ii*2)%col_group = col_group*2+1
				cascading_factors(level_c+1)%matrices_block(ii*2)%row_group = row_group*2+1
				cascading_factors(level_c+1)%matrices_block(ii*2)%style = 1						

				cascading_factors(level_c+1)%BP(ii*2-1)%level = level_c+1
				cascading_factors(level_c+1)%BP(ii*2-1)%col_group = col_group*2
				cascading_factors(level_c+1)%BP(ii*2-1)%row_group = row_group*2
				cascading_factors(level_c+1)%BP(ii*2)%level = level_c+1
				cascading_factors(level_c+1)%BP(ii*2)%col_group = col_group*2+1
				cascading_factors(level_c+1)%BP(ii*2)%row_group = row_group*2+1

				cascading_factors(level_c+1)%BP_inverse(ii*2-1)%level = level_c+1
				cascading_factors(level_c+1)%BP_inverse(ii*2-1)%col_group = col_group*2
				cascading_factors(level_c+1)%BP_inverse(ii*2-1)%row_group = row_group*2
				! cascading_factors(level_c+1)%BP_inverse(ii*2-1)%style = 1
				cascading_factors(level_c+1)%BP_inverse(ii*2)%level = level_c+1
				cascading_factors(level_c+1)%BP_inverse(ii*2)%col_group = col_group*2+1
				cascading_factors(level_c+1)%BP_inverse(ii*2)%row_group = row_group*2+1
				! cascading_factors(level_c+1)%BP_inverse(ii*2)%style = 1	
				
			end if

		end do
	end do
	

	if(schurinv==1)then	
		do level_c = 1,Maxlevel_for_blocks
			do ii = 1, cascading_factors(level_c)%N_block_inverse
				cascading_factors(level_c)%BP_inverse_schur(ii)%level = cascading_factors(level_c)%BP_inverse(ii)%level+1
				! cascading_factors(level_c)%BP_inverse_schur(ii)%style = 2
				cascading_factors(level_c)%BP_inverse_schur(ii)%col_group = cascading_factors(level_c)%BP_inverse(ii)%col_group * 2
				cascading_factors(level_c)%BP_inverse_schur(ii)%row_group = cascading_factors(level_c)%BP_inverse(ii)%row_group * 2
				
			end do
		end do	
	end if
  
	allocate(new2old(Maxedge))    

	do ii=1,Maxedge
		new2old(ii) = ii
	end do
	
	! write(*,*)'gan', node_patch_of_edge(0,100)   
	
	   
    ! allocate (distance(Maxedge))     
    
    !********************************index_of_group**************************************
    
    basis_group(1)%head=1 ; basis_group(1)%tail=Maxedge
    do level=0, Maxlevel
        do group=2**level, 2**(level+1)-1
            basis_group(group)%level=level

            groupcenter(1:dimn)=0.0d0
            ! !$omp parallel do default(shared) private(edge,ii) reduction(+:groupcenter)
            do edge=basis_group(group)%head, basis_group(group)%tail
                do ii=1,dimn
                    groupcenter(ii)=groupcenter(ii)+xyz(ii,node_patch_of_edge(0,edge))
                enddo
				! if(group==24)write(*,*)edge,groupcenter(:),xyz(1,node_patch_of_edge(0,edge))
			enddo
            ! !$omp end parallel do
            do ii=1,dimn
                groupcenter(ii)=groupcenter(ii)/(basis_group(group)%tail-basis_group(group)%head+1)
            enddo
            basis_group(group)%center(1:dimn)=groupcenter(1:dimn)

			radiusmax=0.
			do edge=basis_group(group)%head, basis_group(group)%tail
				radius=0
				do ii=1,dimn
					radius=radius+(xyz(ii,node_patch_of_edge(0,edge))-groupcenter(ii))**2
				enddo
				! write(*,*)'really',edge, node_patch_of_edge(0,edge)
				radius=sqrt(radius)
				if (radius>radiusmax) then
					radiusmax=radius
				endif			
			enddo
			basis_group(group)%radius=radiusmax		
				
			if(size(basis_group_pre,1)<2*group+1)then ! if groups not predefined, need to order the points
				if(xyzsort==1)then !xyz sort		 
					xyzmin= 1d300
					xyzmax= -1d300
					do edge=basis_group(group)%head, basis_group(group)%tail
						do ii=1,Dimn
							xyzmax(ii) = max(xyzmax(ii),xyz(ii,node_patch_of_edge(0,edge)))
							xyzmin(ii) = min(xyzmin(ii),xyz(ii,node_patch_of_edge(0,edge)))
						enddo
					enddo
					xyzrange(1:Dimn) = xyzmax(1:Dimn)-xyzmin(1:Dimn)
					
					mm = basis_group(group)%tail - basis_group(group)%head + 1
					allocate (distance(mm))	
					sortdirec = maxloc(xyzrange(1:Dimn),1)
					! write(*,*)'gaw',sortdirec,xyzrange(1:Dimn)
					
					! if(Kernel==EMSURF)then
					! if(mod(level,2)==1)then           !!!!!!!!!!!!!!!!!!!!!!!!! note: applys only to plates
						! sortdirec=2
					! else 
						! sortdirec=3
					! end if
					! endif
					
					
					!$omp parallel do default(shared) private(i)
					do i=basis_group(group)%head, basis_group(group)%tail
						distance(i-basis_group(group)%head+1)=xyz(sortdirec,node_patch_of_edge(0,i))
					enddo
					!$omp end parallel do

				else if(xyzsort==2)then ! spherical sort
					phi_end = 0
					do edge=basis_group(group)%head, basis_group(group)%tail
						call Cart2Sph(xyz(1,node_patch_of_edge(0,edge)),xyz(2,node_patch_of_edge(0,edge)),xyz(3,node_patch_of_edge(0,edge)),Origins,r,theta,phi)					
						if(abs(phi-2*pi)< 3*minedgelength/r)phi_end=1      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! this group contains points near phi=2pi
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
						call Cart2Sph(xyz(1,node_patch_of_edge(0,i)),xyz(2,node_patch_of_edge(0,i)),xyz(3,node_patch_of_edge(0,i)),Origins,r,theta,phi)
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
				end if
				
				Ninfo_edge=size(node_patch_of_edge,1)-1
				allocate (order(mm),edge_temp(0:Ninfo_edge,mm)) 
				allocate(map_temp(mm))

				call quick_sort(distance,order,mm)       
				!$omp parallel do default(shared) private(ii)     
				do ii=1, mm
					edge_temp(:,ii)=node_patch_of_edge(:,order(ii)+basis_group(group)%head-1)
					map_temp(ii) = new2old(order(ii)+basis_group(group)%head-1)
				enddo
				!$omp end parallel do

				!$omp parallel do default(shared) private(ii)     
				do ii=1, mm
					node_patch_of_edge(:,ii+basis_group(group)%head-1)=edge_temp(:,ii)
					new2old(ii+basis_group(group)%head-1) = map_temp(ii)
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
					
					
					if(xyzsort==1)then 
						seperator = xyz(sortdirec,node_patch_of_edge(0,basis_group(2*group)%tail))
						
					else 
						call Cart2Sph(xyz(1,node_patch_of_edge(0,basis_group(2*group)%tail)),xyz(2,node_patch_of_edge(0,basis_group(2*group)%tail)),xyz(3,node_patch_of_edge(0,basis_group(2*group)%tail)),Origins,r,theta,phi)
						if(sortdirec==1)seperator = theta
						if(sortdirec==2)then
							seperator = phi
							! write(*,*)level,phi*180/pi,basis_group(2*group)%tail,'ganni',phi_end
						end if
					end if
					
					
					fidx = (2*group)- 2**(level+1)+1
					cascading_factors(level+1)%BP(fidx)%boundary(1) = sortdirec
					cascading_factors(level+1)%BP(fidx)%boundary(2) = seperator
					
					fidx = (2*group+1)- 2**(level+1)+1
					cascading_factors(level+1)%BP(fidx)%boundary(1) = sortdirec
					cascading_factors(level+1)%BP(fidx)%boundary(2) = seperator				 					
				endif	

			else 
				if (level<Maxlevel) then
					basis_group(2*group)%head=basis_group_pre(2*group,1)
					basis_group(2*group)%tail=basis_group_pre(2*group,2)
					basis_group(2*group+1)%head=basis_group_pre(2*group+1,1)
					basis_group(2*group+1)%tail=basis_group_pre(2*group+1,2)
					
					fidx = (2*group)- 2**(level+1)+1
					cascading_factors(level+1)%BP(fidx)%boundary(1) = 0   ! dummy parameters
					cascading_factors(level+1)%BP(fidx)%boundary(2) = 0
					
					fidx = (2*group+1)- 2**(level+1)+1
					cascading_factors(level+1)%BP(fidx)%boundary(1) = 0
					cascading_factors(level+1)%BP(fidx)%boundary(2) = 0				 
				endif					
			end if
        enddo
    enddo
	
	
	do level=0, Maxlevel
        do group=2**level, 2**(level+1)-1
            do edge=basis_group(group)%head, basis_group(group)%tail
				! write(*,*)edge,node_patch_of_edge(0,edge)
				! write(113,'(I5,I8,Es16.8,Es16.8,Es16.8)')level,group,xyz(1:Dimn,node_patch_of_edge(0,edge))
				write(113,'(I5,I8,Es16.8,Es16.8)')level,group,xyz(1:Dimn,node_patch_of_edge(0,edge))
			enddo
        enddo
    enddo	   
    !*************************************************************************************
   
	
	
    
!     write (*,*) basis_group(matrices_block(9,0)%col_group)%center(1:2),basis_group(matrices_block(9,0)%row_group)%center(1:2)
!     write (*,*) basis_group(matrices_block(10,0)%col_group)%center(1:2),basis_group(matrices_block(10,0)%row_group)%center(1:2)
!     write (*,*) basis_group(matrices_block(11,0)%col_group)%center(1:2),basis_group(matrices_block(11,0)%row_group)%center(1:2)
!     pause
    
    return
    
end subroutine H_matrix_structuring


subroutine BPlus_structuring()
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
	type(matrixblock),pointer::blocks
	real*8::minbound,theta,phi,r,rmax,phi_tmp,measure
	real*8,allocatable::Centroid_M(:,:),Centroid_N(:,:)
	integer,allocatable::Isboundary_M(:),Isboundary_N(:)
	integer Dimn

	Dimn = size(xyz,1)
	
	do level_c = 1,Maxlevel_for_blocks+1
		do ii =1,cascading_factors(level_c)%N_block_forward
           
			if(level_c==Maxlevel_for_blocks+1)then


				cascading_factors(level_c)%BP(ii)%Lplus=1				
				allocate(cascading_factors(level_c)%BP(ii)%LL(LplusMax))
				do ll=1,LplusMax
					cascading_factors(level_c)%BP(ii)%LL(ll)%Nbound=0
				end do
				cascading_factors(level_c)%BP(ii)%LL(1)%Nbound = 1
				allocate(cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1))
				cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%level = cascading_factors(level_c)%BP(ii)%level
				cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%col_group = cascading_factors(level_c)%BP(ii)%col_group
				cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%row_group = cascading_factors(level_c)%BP(ii)%row_group
				cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%style = 1  !!!!! be careful here
				
				cascading_factors(level_c)%BP_inverse(ii)%Lplus=1
				allocate(cascading_factors(level_c)%BP_inverse(ii)%LL(LplusMax))
				do ll=1,LplusMax
					cascading_factors(level_c)%BP_inverse(ii)%LL(ll)%Nbound=0
				end do	
				cascading_factors(level_c)%BP_inverse(ii)%LL(1)%Nbound = 1
				allocate(cascading_factors(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1))
				cascading_factors(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%level = cascading_factors(level_c)%BP_inverse(ii)%level
				cascading_factors(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%col_group = cascading_factors(level_c)%BP_inverse(ii)%col_group
				cascading_factors(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%row_group = cascading_factors(level_c)%BP_inverse(ii)%row_group
				cascading_factors(level_c)%BP_inverse(ii)%LL(1)%matrices_block(1)%style = 1  !!!!! be careful here				
				
			
			else 
				allocate(cascading_factors(level_c)%BP(ii)%LL(LplusMax))
				do ll=1,LplusMax
					cascading_factors(level_c)%BP(ii)%LL(ll)%Nbound=0
				end do
				
				cascading_factors(level_c)%BP(ii)%LL(1)%Nbound = 1
				allocate(cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1))
				cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%level = cascading_factors(level_c)%BP(ii)%level
				cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%col_group = cascading_factors(level_c)%BP(ii)%col_group
				cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%row_group = cascading_factors(level_c)%BP(ii)%row_group			
				cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%style = 2
				allocate(cascading_factors(level_c)%BP(ii)%LL(1)%boundary_map(1))
				cascading_factors(level_c)%BP(ii)%LL(1)%boundary_map(1) = cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%col_group
				cascading_factors(level_c)%BP(ii)%Lplus=0
				
				sortdirec = NINT(cascading_factors(level_c)%BP(ii)%boundary(1))
				seperator = cascading_factors(level_c)%BP(ii)%boundary(2)
				
				do ll=1,LplusMax-1
					if(cascading_factors(level_c)%BP(ii)%LL(ll)%Nbound>0)then
						cascading_factors(level_c)%BP(ii)%Lplus = cascading_factors(level_c)%BP(ii)%Lplus + 1
						call assert(cascading_factors(level_c)%BP(ii)%Lplus<=LplusMax,'increase LplusMax')
						
						level_butterfly = int((Maxlevel_for_blocks - cascading_factors(level_c)%BP(ii)%LL(ll)%matrices_block(1)%level)/2)*2 
						! write(*,*)'nini',Maxlevel_for_blocks - cascading_factors(level_c)%BP(ii)%LL(ll)%matrices_block(1)%level,LnoBP
						if(Maxlevel_for_blocks - cascading_factors(level_c)%BP(ii)%LL(ll)%matrices_block(1)%level<LnoBP)then
							cascading_factors(level_c)%BP(ii)%LL(ll+1)%Nbound=0
						else 
							level_BP = cascading_factors(level_c)%BP(ii)%level
							levelm = ceiling_safe(dble(level_butterfly)/2d0)						
							groupm_start=cascading_factors(level_c)%BP(ii)%LL(ll)%matrices_block(1)%row_group*2**levelm						
							Nboundall = 2**(cascading_factors(level_c)%BP(ii)%LL(ll)%matrices_block(1)%level+levelm-level_BP)
							allocate(cascading_factors(level_c)%BP(ii)%LL(ll+1)%boundary_map(Nboundall))
							cascading_factors(level_c)%BP(ii)%LL(ll+1)%boundary_map = -1
							cascading_factors(level_c)%BP(ii)%LL(ll+1)%Nbound=0
							
							do bb = 1,cascading_factors(level_c)%BP(ii)%LL(ll)%Nbound
								blocks => cascading_factors(level_c)%BP(ii)%LL(ll)%matrices_block(bb)
								
								allocate(Centroid_M(2**levelm,Dimn))
								allocate(Isboundary_M(2**levelm))
								Isboundary_M = 0
								Centroid_M = 0
								
								do index_i_m=1, 2**levelm
									group_m=blocks%row_group   ! Note: row_group and col_group interchanged here 
									group_m=group_m*2**levelm-1+index_i_m  										

									CNT = 0
									if(xyzsort==1)then	
										do nn = basis_group(group_m)%head,basis_group(group_m)%tail
											measure = abs(xyz(sortdirec,node_patch_of_edge(0,nn))-seperator)
											if(measure<3*minedgelength)then
												Isboundary_M(index_i_m) = 1 
												CNT = CNT + 1
												Centroid_M(index_i_m,1:Dimn) = Centroid_M(index_i_m,1:Dimn)+ xyz(1:Dimn,node_patch_of_edge(0,nn))
											end if
										end do
										if(Isboundary_M(index_i_m)==1)Centroid_M(index_i_m,:) = Centroid_M(index_i_m,:)/CNT
										
										if(blocks%col_group==8 .or. blocks%col_group==9)then
											write(*,*)'wocaoo',group_m,Isboundary_M(index_i_m),CNT
										endif
										
										
									else if(xyzsort==2)then
										do nn = basis_group(group_m)%head,basis_group(group_m)%tail
											call Cart2Sph(xyz(1,node_patch_of_edge(0,nn)),xyz(2,node_patch_of_edge(0,nn)),xyz(3,node_patch_of_edge(0,nn)),Origins,r,theta,phi)
											if(sortdirec==1)then
												measure = abs(theta-seperator)
												if(measure< 3*minedgelength/r)then
													Isboundary_M(index_i_m) = 1 
													CNT = CNT + 1
													Centroid_M(index_i_m,1) = Centroid_M(index_i_m,1) + xyz(1,node_patch_of_edge(0,nn))
													Centroid_M(index_i_m,2) = Centroid_M(index_i_m,2) + xyz(2,node_patch_of_edge(0,nn))
													Centroid_M(index_i_m,3) = Centroid_M(index_i_m,3) + xyz(3,node_patch_of_edge(0,nn))													
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
												if(measure< 3*minedgelength/r)then
													Isboundary_M(index_i_m) = 1 
													CNT = CNT + 1
													Centroid_M(index_i_m,1) = Centroid_M(index_i_m,1) + xyz(1,node_patch_of_edge(0,nn))
													Centroid_M(index_i_m,2) = Centroid_M(index_i_m,2) + xyz(2,node_patch_of_edge(0,nn))
													Centroid_M(index_i_m,3) = Centroid_M(index_i_m,3) + xyz(3,node_patch_of_edge(0,nn))	
													! if(level_c==1 .and. index_i_m==1)write(*,*)CNT, xyz(:,node_patch_of_edge(0,nn))												
												end if												
											end if
										end do										
										! if(level_c==1 .and. index_i_m==1)write(*,*)Centroid_M(index_i_m,:),CNT,'dddddd'
										
										if(Isboundary_M(index_i_m)==1)Centroid_M(index_i_m,:) = Centroid_M(index_i_m,:)/CNT
										! if(level_c==1)write(*,*)Centroid_M(1,:),CNT,'dddddd'
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
									if(xyzsort==1)then	
										do nn = basis_group(group_n)%head,basis_group(group_n)%tail
											measure = abs(xyz(sortdirec,node_patch_of_edge(0,nn))-seperator)
											if(measure<3*minedgelength)then
												Isboundary_N(index_j_m) = 1 
												CNT = CNT + 1
												! write(*,*)nn,index_j_m,'ok'
												! write(*,*)Centroid_N(index_j_m,1),xyz(1,node_patch_of_edge(0,nn))
												Centroid_N(index_j_m,1:Dimn) = Centroid_N(index_j_m,1:Dimn) + xyz(1:Dimn,node_patch_of_edge(0,nn))
											end if
										end do
										if(Isboundary_N(index_j_m)==1)Centroid_N(index_j_m,:) = Centroid_N(index_j_m,:)/CNT
										
									else if(xyzsort==2)then
										do nn = basis_group(group_n)%head,basis_group(group_n)%tail
											call Cart2Sph(xyz(1,node_patch_of_edge(0,nn)),xyz(2,node_patch_of_edge(0,nn)),xyz(3,node_patch_of_edge(0,nn)),Origins,r,theta,phi)
											if(sortdirec==1)then
												measure = abs(theta-seperator)
												if(measure< 3*minedgelength/r)then
													Isboundary_N(index_j_m) = 1 
													CNT = CNT + 1
													Centroid_N(index_j_m,1) = Centroid_N(index_j_m,1) + xyz(1,node_patch_of_edge(0,nn))
													Centroid_N(index_j_m,2) = Centroid_N(index_j_m,2) + xyz(2,node_patch_of_edge(0,nn))
													Centroid_N(index_j_m,3) = Centroid_N(index_j_m,3) + xyz(3,node_patch_of_edge(0,nn))													
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
												if(measure< 3*minedgelength/r)then
													Isboundary_N(index_j_m) = 1 
													CNT = CNT + 1
													Centroid_N(index_j_m,1) = Centroid_N(index_j_m,1) + xyz(1,node_patch_of_edge(0,nn))
													Centroid_N(index_j_m,2) = Centroid_N(index_j_m,2) + xyz(2,node_patch_of_edge(0,nn))
													Centroid_N(index_j_m,3) = Centroid_N(index_j_m,3) + xyz(3,node_patch_of_edge(0,nn))													
												end if												
											end if
										end do										
										if(Isboundary_N(index_j_m)==1)Centroid_N(index_j_m,:) = Centroid_N(index_j_m,:)/CNT
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
										cascading_factors(level_c)%BP(ii)%LL(ll+1)%Nbound = cascading_factors(level_c)%BP(ii)%LL(ll+1)%Nbound + 1									
										dist = 100000000d0
										do index_j_m=1, 2**(level_butterfly-levelm)	
											group_n=blocks%col_group  
											group_n=group_n*2**(level_butterfly-levelm)-1+index_j_m
											
											if(blocks%col_group==8 .or. blocks%col_group==9)then
												write(*,*)group_m,group_n,sqrt(sum((Centroid_N(index_j_m,:)-Centroid_M(index_i_m,:))**2d0)),'nima'
											end	if
											
											if(Isboundary_N(index_j_m)==1)then
												if(dist > sqrt(sum((Centroid_N(index_j_m,:)-Centroid_M(index_i_m,:))**2d0)))then
													! if(level_c==1)write(*,*)index_i_m,index_j_m
													dist = sqrt(sum((Centroid_N(index_j_m,:)-Centroid_M(index_i_m,:))**2d0))
													cascading_factors(level_c)%BP(ii)%LL(ll+1)%boundary_map(group_m - groupm_start + 1) = group_n
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
							
							if(cascading_factors(level_c)%BP(ii)%LL(ll+1)%Nbound>1)then
								write(*,*)level_c,ii,ll,cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%col_group,cascading_factors(level_c)%BP(ii)%LL(1)%matrices_block(1)%row_group,cascading_factors(level_c)%BP(ii)%LL(ll+1)%Nbound,'niamaa'
							endif
							
							call assert(cascading_factors(level_c)%BP(ii)%LL(ll+1)%Nbound>0,'why is no boundary group detected')	
								
							allocate(cascading_factors(level_c)%BP(ii)%LL(ll+1)%matrices_block(cascading_factors(level_c)%BP(ii)%LL(ll+1)%Nbound))
							
							cnt = 0
							do bb = 1,	Nboundall
								if(cascading_factors(level_c)%BP(ii)%LL(ll+1)%boundary_map(bb)/=-1)then
									cnt = cnt + 1
									group_m = bb+groupm_start-1
									group_n = cascading_factors(level_c)%BP(ii)%LL(ll+1)%boundary_map(bb)
									
									cascading_factors(level_c)%BP(ii)%LL(ll+1)%matrices_block(cnt)%row_group = group_m
									cascading_factors(level_c)%BP(ii)%LL(ll+1)%matrices_block(cnt)%col_group = group_n
									cascading_factors(level_c)%BP(ii)%LL(ll+1)%matrices_block(cnt)%level = basis_group(group_m)%level
									cascading_factors(level_c)%BP(ii)%LL(ll+1)%matrices_block(cnt)%style = 2
									
								end if
							end do
						end if		
					else 
						exit
					end if
				end do

				
				! write(*,*)level_c,ii,cascading_factors(level_c)%BP(ii)%Lplus,'gaogao '
				
				if(mod(ii,2)==1)then
					ii_sch = ceiling_safe(ii/2d0)
					
					allocate(cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(LplusMax))
					cascading_factors(level_c)%BP_inverse_schur(ii_sch)%Lplus=cascading_factors(level_c)%BP(ii)%Lplus
					do ll=1,LplusMax
						cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%Nbound=0
					end do	
									
					cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(1)%Nbound = 1
					cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(1)%rankmax = 0
					allocate(cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(1)%matrices_block(1))
					cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(1)%matrices_block(1)%level = cascading_factors(level_c)%BP(ii)%level
					cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(1)%matrices_block(1)%col_group = cascading_factors(level_c)%BP(ii)%row_group
					cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(1)%matrices_block(1)%row_group = cascading_factors(level_c)%BP(ii)%row_group			
					cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(1)%matrices_block(1)%style = 2
					allocate(cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(1)%boundary_map(1))
					cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(1)%boundary_map(1) = cascading_factors(level_c)%BP(ii)%row_group		
					
								
					do ll=1,LplusMax-1
						if(cascading_factors(level_c)%BP(ii)%LL(ll)%Nbound>0)then

							cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%rankmax = 0
							cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%Nbound = cascading_factors(level_c)%BP(ii)%LL(ll)%Nbound
							
							allocate(cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%Nbound))
							
							do bb =1,cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%Nbound
								cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(bb)%row_group = cascading_factors(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%row_group
								cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(bb)%col_group = cascading_factors(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%row_group
								cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(bb)%style = cascading_factors(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%style
								cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(bb)%level = cascading_factors(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%level
							end do
							
							
							if(cascading_factors(level_c)%BP(ii)%LL(ll+1)%Nbound==0)then		
								cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll+1)%Nbound=0
							else 
								level_butterfly = int((Maxlevel_for_blocks - cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(1)%level)/2)*2 
								level_BP = cascading_factors(level_c)%BP_inverse_schur(ii_sch)%level
								levelm = ceiling_safe(dble(level_butterfly)/2d0)						
								groupm_start=cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(1)%row_group*2**levelm		
								Nboundall = 2**(cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll)%matrices_block(1)%level+levelm-level_BP)				
								
								allocate(cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll+1)%boundary_map(Nboundall))
								! write(*,*)shape(Bplus%LL(ll+1)%boundary_map),shape(Bplus_randomized(1)%LL(ll+1)%boundary_map),'didi',ll
								
								cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll+1)%boundary_map = cascading_factors(level_c)%BP(ii)%LL(ll+1)%boundary_map
								do bb=1,Nboundall
									if(cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll+1)%boundary_map(bb)/=-1)then
										cascading_factors(level_c)%BP_inverse_schur(ii_sch)%LL(ll+1)%boundary_map(bb) = bb + groupm_start - 1
									end if
								end do
							end if
						else 
							exit
						end if
					end do				
				end if			
								
				! if(level_c==1 .and. ii==1)then
				 
				! write(177,*)'Bplus:', level_c,ii
				do ll=1,cascading_factors(level_c)%BP(ii)%Lplus
				! write(*,*)cascading_factors(level_c)%BP(ii)%LL(ll)%Nbound,'ddd'
					do bb = 1,cascading_factors(level_c)%BP(ii)%LL(ll)%Nbound
						write(177,'(I3,I7,I3,I3,Es16.7,Es16.7,Es16.7,Es16.7,Es16.7,Es16.7)')level_c,ii,ll,cascading_factors(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%level,basis_group(cascading_factors(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%row_group)%center(1:dimn),basis_group(cascading_factors(level_c)%BP(ii)%LL(ll)%matrices_block(bb)%col_group)%center(1:dimn)
					end do
				end do
				! end if
			end if
			
		end do
	end do	

	
	
	
end subroutine BPlus_structuring


end module H_structure 