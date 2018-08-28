module APPLICATION_MODULE
use z_HODLR_DEFS
use z_misc
implicit none

	!**** define your application-related variables here   
	type quant_app
		real(kind=8) wavenum    ! CEM: wave number  
		real(kind=8) wavelength  ! CEM: wave length
		real(kind=8) omiga       ! CEM: angular frequency
		real(kind=8) rank_approximate_para1, rank_approximate_para2, rank_approximate_para3 ! CEM: rank estimation parameter  
		integer RCS_static  ! CEM: 1: monostatic or 2: bistatic RCS
		integer RCS_Nsample ! CEM: number of RCS samples
		real(kind=8):: CFIE_alpha ! CEM: combination parameter in CFIE
	end type quant_app

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
			type is (quant_app)
				! convert to new indices because msh%info_unk has been reordered
				edge_m = msh%old2new(m)
				edge_n = msh%old2new(n)
				
				if (edge_m/=edge_n) then
				
					flag=0
					do j=1,2
						do i=1,2
							if (msh%info_unk(i,edge_m)==msh%info_unk(j,edge_n)) then
								flag=1
							endif
						enddo
					enddo
					
					if (flag==1) then
						! write(*,*)msh%Delta_ll,quant%wavenum,impedence0,junit,gamma,'ddd' 
						value_e=quant%wavenum*impedence0/4.0*msh%Delta_ll*(1-junit/pi*(3*LOG(3*gamma*quant%wavenum*msh%Delta_ll/4.0)-LOG(gamma*quant%wavenum*msh%Delta_ll/4.0)-2))
					
					else
				
						r_mn=(msh%xyz(1,msh%info_unk(0,edge_m))-msh%xyz(1,msh%info_unk(0,edge_n)))**2+(msh%xyz(2,msh%info_unk(0,edge_m))-msh%xyz(2,msh%info_unk(0,edge_n)))**2
						r_mn=sqrt(r_mn)
						value_e=quant%wavenum*impedence0/4.0*msh%Delta_ll*z_Hankel02_Func(quant%wavenum*r_mn)
					endif
				else
					value_e=quant%wavenum*impedence0/4.0*msh%Delta_ll*(1.0-junit*2.0/pi*(LOG(gamma*quant%wavenum*msh%Delta_ll/4.0)-1.0))
				endif

			class default
				write(*,*)"unexpected type"
				stop
			end select			
			
		return
		
	end subroutine Z_elem_EMCURV	
	
	
	

	subroutine VV_polar_CURV(dphi,edge,ctemp_1,curr,msh,quant)
		
		use z_HODLR_DEFS
		! use APPLICATION_MODULE
		implicit none
		complex(kind=8)::curr
		complex(kind=8) ctemp,phase,ctemp_1
		real(kind=8) dsita,dphi    
		integer edge
		type(z_mesh)::msh
		type(quant_app)::quant
		
		phase=junit*quant%wavenum*(msh%xyz(1,msh%info_unk(0,edge))*cos(dphi*pi/180.)+msh%xyz(2,msh%info_unk(0,edge))*sin(dphi*pi/180.))
		ctemp_1=curr*msh%Delta_ll*exp(phase)

		return
		
	end subroutine VV_polar_CURV

	subroutine RCS_bistatic_CURV(curr,msh,quant,ptree)
		!integer flag
		use z_HODLR_DEFS
		! use APPLICATION_MODULE
		implicit none
		complex(kind=8)::curr(:)
		real(kind=8) rcs
		complex(kind=8) ctemp_rcs(3),ctemp,phase,ctemp_1,ctemp_loc
		real(kind=8) ddphi,dphi
		
		integer i,j,ii,jj,iii,jjj,patch,flag
		real(kind=8) l_edge,l_edgefine
		type(z_mesh)::msh
		integer edge,edge_m,edge_n,ierr
		type(quant_app)::quant 
		type(z_proctree)::ptree
		
		if(ptree%MyID==Main_ID)open (100, file='VV_bistatic.txt')
	  
		ddphi=180./quant%RCS_Nsample
		
		do i=0, quant%RCS_Nsample   !phi=0
			dphi=i*ddphi
			ctemp_loc=0
			rcs=0
			!$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)   
			do edge=msh%idxs,msh%idxe
				call VV_polar_CURV(dphi,edge,ctemp_1,curr(edge-msh%idxs+1),msh,quant)
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
		! use APPLICATION_MODULE
		implicit none
		
		real(kind=8) rcs
		complex(kind=8) ctemp_rcs(3),ctemp,ctemp_loc,phase,ctemp_1,ctemp_2
		real(kind=8) dsita,dphi
		integer edge,edge_m,edge_n,ierr
		complex(kind=8):: curr(:)
		type(z_mesh)::msh
		type(quant_app)::quant
		type(z_proctree)::ptree
		
		ctemp_loc=0
			rcs=0
			!$omp parallel do default(shared) private(edge,ctemp_1) reduction(+:ctemp_loc)
			do edge=msh%idxs,msh%idxe
				call VV_polar_CURV(dphi,edge,ctemp_1,curr(edge-msh%idxs+1),msh,quant)
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
		! use APPLICATION_MODULE
		implicit none
		
		integer edge
		complex(kind=8) value
		real(kind=8) theta, phi
		complex(kind=8)  phase
		type(z_mesh)::msh
		type(quant_app)::quant
		
		phase=junit*quant%wavenum*(msh%xyz(1,msh%info_unk(0,edge))*cos(phi*pi/180.)+msh%xyz(2,msh%info_unk(0,edge))*sin(phi*pi/180.))
		value=exp(phase)
		
		return
		
	end subroutine element_Vinc_VV_CURV	
		

end module APPLICATION_MODULE	



PROGRAM HODLR_BUTTERFLY_SOLVER_2D
    use z_HODLR_DEFS
    use APPLICATION_MODULE
	! use geometry_model
	use z_H_structure
	use z_cascading_factorization
	use z_HODLR_construction
	use omp_lib
	use z_misc
    implicit none

	! include "mkl_vml.fi"	 
	
    real(kind=8) para
    real(kind=8) tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seed_myid(12)
	integer times(8)	
	real(kind=8) t1,t2,x,y,z,r,theta,phi
	complex(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:)
	
	character(len=:),allocatable  :: string
	character(len=1024)  :: strings	
	character(len=6)  :: info_env	
	integer :: length
	integer :: ierr
	integer*8 oldmode,newmode
	type(z_Hoption)::option	
	type(z_Hstat)::stats
	type(z_mesh)::msh
	type(z_kernelquant)::ker
	type(z_hobf)::ho_bf,ho_bf_copy
	integer,allocatable:: groupmembers(:)
	integer nmpi
	type(z_proctree)::ptree
	type(quant_app),target::quant
	CHARACTER (LEN=1000) DATA_DIR	
	
	
	! nmpi and groupmembers should be provided by the user 
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo	
	
	call z_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	
	if(ptree%MyID==Main_ID)write(*,*)'NUMBER_MPI=',nmpi
	
 	threads_num=1
    CALL getenv("OMP_NUM_THREADS", strings)
	strings = TRIM(strings)	
	if(LEN_TRIM(strings)>0)then
		read(strings , *) threads_num
	endif
	if(ptree%MyID==Main_ID)write(*,*)'OMP_NUM_THREADS=',threads_num
	
	call OMP_set_num_threads(threads_num)		
		
		
	call DATE_AND_TIME(values=times)     ! Get the current time 
	seed_myid(1) = times(4) * (360000*times(5) + 6000*times(6) + 100*times(7) + times(8))
	
	! seed_myid(1) = myid*1000
	call RANDOM_SEED(PUT=seed_myid)
	
	! oldmode = vmlsetmode(VML_FTZDAZ_ON)
	! call vmlsetmode(VML_FTZDAZ_ON)
	if(ptree%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "HODLR_BUTTERFLY_SOLVER_2D"
    write(*,*) "   "
	endif
	
	call z_InitStat(stats)
	call z_SetDefaultOptions(option)
	
	time_tmp = 0
	
	msh%Origins=(/0d0,0d0,0d0/)
	
 	! register the user-defined function and type in ker 
	ker%FuncZmn=>Z_elem_EMCURV
	ker%QuantZmn=>quant

     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	CALL getarg(1, DATA_DIR)
	
	
	msh%model2d=10 ! Model type (1=strip; 2=corner reflector; 3=two opposite strips; 4=CR with RRS; 5=cylinder; 6=Rectangle Cavity); 7=half cylinder; 8=corrugated half cylinder; 9=corrugated corner reflector; 10=open polygon; 11=taller open polygon 
	! msh%Nunk=1280000
	! msh%Nunk=320000
	! msh%Nunk=80000
	! msh%Nunk=20000
	 msh%Nunk=5000
    ! msh%Nunk=160000	
	! Refined_level=0
	
	
	
	
	! ker%Kernel = FULL	
	! allocate(ker%matZ_glo(msh%Nunk,msh%Nunk))
	! call RandomMat(msh%Nunk,msh%Nunk,msh%Nunk,ker%matZ_glo,0)
	! call MPI_Bcast(ker%matZ_glo,msh%Nunk*msh%Nunk,MPI_DOUBLE_COMPLEX,0,MPI_Comm_World,ierr)
	
	
	
	
	option%preorder=0 
	msh%scaling=1d0
	! quant%wavelength=0.0006
	!quant%wavelength=0.0003
	quant%wavelength=0.08

! quant%wavelength=0.08
! Discret=0.05
	quant%RCS_static=1
    quant%RCS_Nsample=2000

	
	
	quant%rank_approximate_para1=6.0
    quant%rank_approximate_para2=6.0
    quant%rank_approximate_para3=6.0
	
	
	option%Nmin_leaf=50
	option%tol_comp=1d-4
	option%tol_Rdetect=3d-5	
	option%tol_LS=1d-12
	option%tol_itersol=1d-6
	option%n_iter=1000
	option%tol_rand=1d-3
	option%level_check=100
	option%precon=DIRECT
	option%xyzsort=TM
	option%lnoBP=40000
	option%TwoLayerOnly=1
	option%touch_para=3
    option%schulzorder=3
    option%schulzlevel=3000
	option%LRlevel=0
	option%ErrFillFull=1
	option%RecLR_leaf=BACA
	option%ErrSol=1

	! call MKL_set_num_threads(NUM_Threads)    ! this overwrites omp_set_num_threads for MKL functions 
	
    ! call OMP_set_dynamic(.true.)
    !call OMP_set_nested(.true.)
	

    !*********************************************************

    quant%omiga=2*pi/quant%wavelength/sqrt(mu0*eps0)
    quant%wavenum=2*pi/quant%wavelength

   !***********************************************************************
   if(ptree%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'wavelength:',quant%wavelength
   write (*,*) ''
   endif
   !***********************************************************************
	
	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID)write(*,*) "geometry modeling......"
    call geo_modeling_CURV(msh,quant,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "modeling finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1

	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "constructing H_matrices formatting......"
    call z_H_matrix_structuring(ho_bf,option,msh,ptree)
	call z_BPlus_structuring(ho_bf,option,msh,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices formatting finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1
	! stop 
	
    !pause
    
    !call compression_test()
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices filling......"
    ! call matrices_construction(ho_bf,option,stats,msh,ker,element_Zmn_FULL,ptree)
    call z_matrices_construction(ho_bf,option,stats,msh,ker,z_element_Zmn_user,ptree)
	! if(option%precon/=DIRECT)then
		! call copy_HOBF(ho_bf,ho_bf_copy)	
	! end if
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices filling finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	! write(*,*)t2-t1
	
	if(option%precon/=NOPRECON)then
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing......"
    call z_cascading_factorizing(ho_bf,option,stats,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	end if

    if(ptree%MyID==Main_ID)write(*,*) "EM_solve......"
    call EM_solve_CURV(ho_bf,option,msh,quant,ptree,stats)
    if(ptree%MyID==Main_ID)write(*,*) "EM_solve finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	
	
    if(ptree%MyID==Main_ID)write(*,*) "-------------------------------program end-------------------------------------"

    ! ! ! ! pause

end PROGRAM HODLR_BUTTERFLY_SOLVER_2D





subroutine geo_modeling_CURV(msh,quant,ptree)

    use z_HODLR_DEFS
	use APPLICATION_MODULE
    implicit none
    type(z_mesh)::msh
	type(quant_app)::quant
    integer i,j,ii,jj,iii,jjj
    integer intemp
    integer node, patch, edge, flag
    integer node1, node2, num_node,Maxedge
    integer node_temp(2)
    real(kind=8) dx, xx, yy, rr, theta,L,M,Am,tt,L1,L2,L3
    
    real(kind=8),allocatable :: node_xy_original(:,:)
    integer,allocatable :: num_edge_of_node(:)
    
    real(kind=8) a(3),b(3),c(3),r0, phi_start
	type(z_proctree)::ptree
	
	Maxedge	= msh%Nunk

	msh%Ncorner = 0
	
    if (msh%model2d==1) then  !*****single strip*****
        
        msh%Delta_ll=2d0/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        msh%xyz(1,0)=-1d0 ; msh%xyz(2,0)=0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, msh%maxnode-1
            dx=node*msh%Delta_ll/2
            msh%xyz(1,node)=dx-1d0
            msh%xyz(2,node)=0
        enddo
        !$omp end parallel do
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1   ! 1 start 0 center 2 end
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo           
        
    elseif (msh%model2d==2) then  !*****corner reflector*****
    
        msh%Delta_ll=2d0/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        
        msh%xyz(1,0)=-1d0/sqrt(2d0) ; msh%xyz(2,0)=1d0/sqrt(2d0)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*msh%Delta_ll/2
            msh%xyz(1,node)=(dx-1d0)/sqrt(2d0)
            msh%xyz(2,node)=(1d0-dx)/sqrt(2d0)
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*msh%Delta_ll/2
            msh%xyz(1,node+Maxedge)=dx/sqrt(2d0)
            msh%xyz(2,node+Maxedge)=dx/sqrt(2d0)
        enddo
        !$omp end parallel do
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo
		
		msh%Ncorner = 1
		allocate(msh%corner_points(2,msh%Ncorner))
		msh%corner_points(1,1) = 0
		msh%corner_points(2,1) = 0        
    
    elseif (msh%model2d==3) then  !*****two opposite strips*****
    
        msh%Delta_ll=2d0/Maxedge
        msh%maxnode=2*Maxedge+2
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        
        msh%xyz(1,0)=-0.5d0 ; msh%xyz(2,0)=1.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*msh%Delta_ll/2
            msh%xyz(1,node)=dx-0.5d0
            msh%xyz(2,node)=1.0d0
        enddo
        !$omp end parallel do
        msh%xyz(1,Maxedge+1)=-0.5d0 ; msh%xyz(2,Maxedge+1)=-1.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*msh%Delta_ll/2
            msh%xyz(1,node+Maxedge+1)=dx-0.5d0
            msh%xyz(2,node+Maxedge+1)=-1.0d0
        enddo
        !$omp end parallel do
        
        intemp=int(Maxedge/2)
        
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, intemp
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo
        
        msh%info_unk(1,intemp+1)=msh%maxnode/2 ; msh%info_unk(2,intemp+1)=msh%maxnode/2+2 ; msh%info_unk(0,intemp+1)=msh%maxnode/2+1
        do edge=intemp+2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo
    
    
    elseif (msh%model2d==5) then !************cylinder*****************
        
        msh%Delta_ll=  2.0d0*pi/Maxedge
        msh%maxnode=2*Maxedge
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        msh%xyz(1,0)=1.0d0 ; msh%xyz(2,0)=0.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, msh%maxnode-1
            dx=node*msh%Delta_ll/2
            msh%xyz(1,node)=cos(dx)
            msh%xyz(2,node)=sin(dx)
        enddo
        !$omp end parallel do
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge-1
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo     
        msh%info_unk(1,Maxedge)=msh%info_unk(2,Maxedge-1)
        msh%info_unk(2,Maxedge)=0
        msh%info_unk(0,Maxedge)=msh%info_unk(0,Maxedge-1)+2      
        
    elseif (msh%model2d==6) then !***********cavity******************
    
        msh%Delta_ll=11d0/3d0/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        dx=msh%Delta_ll/2d0
        msh%xyz(1,0)=-1d0/6. ; msh%xyz(2,0)=1.0d0
        node=0
        do while (flag==0) 
            node=node+1
            msh%xyz(1,node)=msh%xyz(1,node-1)-dx
            msh%xyz(2,node)=msh%xyz(2,node-1)
            if (msh%xyz(1,node)<=-0.5d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            msh%xyz(1,node)=msh%xyz(1,node-1)
            msh%xyz(2,node)=msh%xyz(2,node-1)-dx
            if (msh%xyz(2,node)<=0.) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            msh%xyz(1,node)=msh%xyz(1,node-1)+dx
            msh%xyz(2,node)=msh%xyz(2,node-1)
            if (msh%xyz(1,node)>=0.5d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            msh%xyz(1,node)=msh%xyz(1,node-1)
            msh%xyz(2,node)=msh%xyz(2,node-1)+dx
            if (msh%xyz(2,node)>=1.0d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            msh%xyz(1,node)=msh%xyz(1,node-1)-dx
            msh%xyz(2,node)=msh%xyz(2,node-1)
            if (node>=msh%maxnode-1) then
                flag=1
            endif
        enddo
        
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo

		msh%Ncorner = 4
		allocate(msh%corner_points(2,msh%Ncorner))
		msh%corner_points(1,1) = -0.5
		msh%corner_points(2,1) = 0     
		msh%corner_points(1,2) = 0.5
		msh%corner_points(2,2) = 0
		msh%corner_points(1,3) = 0.5
		msh%corner_points(2,3) = 1
		msh%corner_points(1,4) = -0.5
		msh%corner_points(2,4) = 1		
		
    elseif (msh%model2d==7) then   !************open cylinder*****************
    
	
	
        msh%Delta_ll=1d0*pi/Maxedge !2.0d0*pi*5d0/6d0/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        
        msh%xyz(1,0)=cos(0*pi) ; msh%xyz(2,0)=sin(0*pi)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*msh%Delta_ll
            msh%xyz(1,node*2)=cos(0*pi+dx)
            msh%xyz(2,node*2)=sin(0*pi+dx)
        enddo
        !$omp end parallel do
		
		
        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            msh%xyz(1,node*2-1)=(msh%xyz(1,node*2-2)+msh%xyz(1,node*2))/2d0
            msh%xyz(2,node*2-1)=(msh%xyz(2,node*2-2)+msh%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do		
		
		
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo           

    elseif (msh%model2d==8) then   !************corrugated open cylinder*****************
		M = 0.2d0*quant%wavelength
		L = 1.5d0*quant%wavelength
        phi_start = 3d0/2d0*pi
		
		msh%Delta_ll=1d0*pi/Maxedge !2.0d0*pi*5d0/6d0/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        
        msh%xyz(1,0)=cos(phi_start) ; msh%xyz(2,0)=sin(phi_start)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*msh%Delta_ll
            msh%xyz(1,node*2)=(1+M*sin(2*pi*dx/L))*cos(phi_start+dx)
            msh%xyz(2,node*2)=(1+M*sin(2*pi*dx/L))*sin(phi_start+dx)
        enddo
        !$omp end parallel do
		
		
        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            msh%xyz(1,node*2-1)=(msh%xyz(1,node*2-2)+msh%xyz(1,node*2))/2d0
            msh%xyz(2,node*2-1)=(msh%xyz(2,node*2-2)+msh%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do		
		
		
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo    

		
    elseif (msh%model2d==9) then  !*****corrugated corner reflector*****
		M = 0.2d0*quant%wavelength
		L = 1.5d0*quant%wavelength
		
        msh%Delta_ll=2d0/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        

		Am = M*sin(pi/2-2*pi/L)-M		
		msh%xyz(1,0)=-1d0/sqrt(2d0)+Am/sqrt(2d0); msh%xyz(2,0)=1d0/sqrt(2d0)+Am/sqrt(2d0)
		
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge/2
            dx=node*msh%Delta_ll
			Am = M*sin(2*pi*dx/L+pi/2-2*pi/L)-M
			
            msh%xyz(1,node*2)=(dx-1d0)/sqrt(2d0)+Am/sqrt(2d0)
            msh%xyz(2,node*2)=(1d0-dx)/sqrt(2d0)+Am/sqrt(2d0)
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge/2
            dx=node*msh%Delta_ll
			Am = M*sin(2*pi*(dx+1)/L+pi/2-2*pi/L)-M
			
            msh%xyz(1,node*2+Maxedge)=dx/sqrt(2d0)-Am/sqrt(2d0)
            msh%xyz(2,node*2+Maxedge)=dx/sqrt(2d0)+Am/sqrt(2d0)
        enddo
        !$omp end parallel do

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            msh%xyz(1,node*2-1)=(msh%xyz(1,node*2-2)+msh%xyz(1,node*2))/2d0
            msh%xyz(2,node*2-1)=(msh%xyz(2,node*2-2)+msh%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo  		

		msh%Ncorner = 1
		allocate(msh%corner_points(2,msh%Ncorner))
		msh%corner_points(1,1) = 0
		msh%corner_points(2,1) = 0
		
    elseif (msh%model2d==10) then  !*****cup*****

		L1 = 0
		L2 = sqrt(2d0)
		L3 = sqrt(2d0)
		L = 3
        msh%Delta_ll=(2*L1+2*L2+2*L3+L)/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        
		
        !$omp parallel do default(shared) private(node,dx,tt)
        do node=0, Maxedge
            dx=node*msh%Delta_ll
			if(dx<=L1)then
				tt = dx
				msh%xyz(1,node*2)= -(L/2-(L1-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L1+L2+L3)/sqrt(2d0)-tt/sqrt(2d0)			
			else if(dx<=L1+L2)then
				tt = dx - L1
				msh%xyz(1,node*2)= -(L/2-(-L2+L3)/sqrt(2d0))+tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L2+L3)/sqrt(2d0)-tt/sqrt(2d0)
			else if(dx<=L1+L2+L3)then
				tt = dx - (L1+L2)
				msh%xyz(1,node*2)= -(L/2-L3/sqrt(2d0))-tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L3)/sqrt(2d0)-tt/sqrt(2d0)				
			else if(dx<=L1+L2+L3+L)then
				tt = dx - (L1+L2+L3)
				msh%xyz(1,node*2)= -L/2+tt
				msh%xyz(2,node*2)= 0
			else if(dx<=L1+L2+L3+L+L3)then
				tt = dx - (L1+L2+L3+L)
				msh%xyz(1,node*2)= L/2-tt/sqrt(2d0)
				msh%xyz(2,node*2)= tt/sqrt(2d0)		
			else if(dx<=L1+L2+L3+L+L3+L2)then
				tt = dx - (L1+L2+L3+L+L3)
				msh%xyz(1,node*2)= (L/2-L3/sqrt(2d0))+tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L3)/sqrt(2d0)+tt/sqrt(2d0)				
			else
				tt = dx - (L1+L2+L3+L+L3+L2)
				msh%xyz(1,node*2)= (L/2-(-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L3+L2)/sqrt(2d0)+tt/sqrt(2d0)							
			end if 
        enddo
        !$omp end parallel do			

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            msh%xyz(1,node*2-1)=(msh%xyz(1,node*2-2)+msh%xyz(1,node*2))/2d0
            msh%xyz(2,node*2-1)=(msh%xyz(2,node*2-2)+msh%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo  		

		! msh%Ncorner = 4
		! allocate(msh%corner_points(2,msh%Ncorner))
		! msh%corner_points(1,1) = -L/2
		! msh%corner_points(2,1) = 0
		! msh%corner_points(1,2) = L/2
		! msh%corner_points(2,2) = 0
		! msh%corner_points(1,3) = -(L/2-L3/sqrt(2d0))
		! msh%corner_points(2,3) = (L3)/sqrt(2d0)		
		! msh%corner_points(1,4) = (L/2-L3/sqrt(2d0))
		! msh%corner_points(2,4) = (L3)/sqrt(2d0)			

    elseif (msh%model2d==11) then  !*****longer cup*****
		L1 = 1
		L2 = sqrt(2d0)
		L3 = sqrt(2d0)
		L = 3.5
        msh%Delta_ll=(2*L1+2*L2+2*L3+L)/Maxedge
        msh%maxnode=2*Maxedge+1
        allocate (msh%xyz(2,0:msh%maxnode-1), msh%info_unk(0:2,Maxedge))
        
		
        !$omp parallel do default(shared) private(node,dx,tt)
        do node=0, Maxedge
            dx=node*msh%Delta_ll
			if(dx<=L1)then
				tt = dx
				msh%xyz(1,node*2)= -(L/2-(L1-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L1+L2+L3)/sqrt(2d0)-tt/sqrt(2d0)			
			else if(dx<=L1+L2)then
				tt = dx - L1
				msh%xyz(1,node*2)= -(L/2-(-L2+L3)/sqrt(2d0))+tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L2+L3)/sqrt(2d0)-tt/sqrt(2d0)
			else if(dx<=L1+L2+L3)then
				tt = dx - (L1+L2)
				msh%xyz(1,node*2)= -(L/2-L3/sqrt(2d0))-tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L3)/sqrt(2d0)-tt/sqrt(2d0)				
			else if(dx<=L1+L2+L3+L)then
				tt = dx - (L1+L2+L3)
				msh%xyz(1,node*2)= -L/2+tt
				msh%xyz(2,node*2)= 0
			else if(dx<=L1+L2+L3+L+L3)then
				tt = dx - (L1+L2+L3+L)
				msh%xyz(1,node*2)= L/2-tt/sqrt(2d0)
				msh%xyz(2,node*2)= tt/sqrt(2d0)		
			else if(dx<=L1+L2+L3+L+L3+L2)then
				tt = dx - (L1+L2+L3+L+L3)
				msh%xyz(1,node*2)= (L/2-L3/sqrt(2d0))+tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L3)/sqrt(2d0)+tt/sqrt(2d0)				
			else
				tt = dx - (L1+L2+L3+L+L3+L2)
				msh%xyz(1,node*2)= (L/2-(-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				msh%xyz(2,node*2)= (L3+L2)/sqrt(2d0)+tt/sqrt(2d0)							
			end if 
        enddo
        !$omp end parallel do			

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            msh%xyz(1,node*2-1)=(msh%xyz(1,node*2-2)+msh%xyz(1,node*2))/2d0
            msh%xyz(2,node*2-1)=(msh%xyz(2,node*2-2)+msh%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        msh%info_unk(1,1)=0 ; msh%info_unk(2,1)=2 ; msh%info_unk(0,1)=1
        do edge=2, Maxedge
            msh%info_unk(1,edge)=msh%info_unk(2,edge-1)
            msh%info_unk(2,edge)=msh%info_unk(2,edge-1)+2
            msh%info_unk(0,edge)=msh%info_unk(0,edge-1)+2
        enddo  		

		msh%Ncorner = 6
		allocate(msh%corner_points(2,msh%Ncorner))
		msh%corner_points(1,1) = -(L/2-(-L2+L3)/sqrt(2d0))
		msh%corner_points(2,1) = (L2+L3)/sqrt(2d0)
		msh%corner_points(1,2) = (L/2-(-L2+L3)/sqrt(2d0))
		msh%corner_points(2,2) = (L2+L3)/sqrt(2d0)	
		msh%corner_points(1,3) = -L/2
		msh%corner_points(2,3) = 0
		msh%corner_points(1,4) = L/2
		msh%corner_points(2,4) = 0
		msh%corner_points(1,5) = -(L/2-L3/sqrt(2d0))
		msh%corner_points(2,5) = (L3)/sqrt(2d0)		
		msh%corner_points(1,6) = (L/2-L3/sqrt(2d0))
		msh%corner_points(2,6) = (L3)/sqrt(2d0)		
    endif
    
	msh%maxedgelength = 0
	do edge=1,Maxedge
		msh%maxedgelength = max(msh%maxedgelength,sqrt(sum(abs(msh%xyz(:,msh%info_unk(1,edge))-msh%xyz(:,msh%info_unk(2,edge))))**2))
	end do	
	
	msh%minedgelength = 10000000
	do edge=1,Maxedge
		msh%minedgelength = min(msh%minedgelength,sqrt(sum(abs(msh%xyz(:,msh%info_unk(1,edge))-msh%xyz(:,msh%info_unk(2,edge)))**2)))
	end do		
	
	msh%corner_radius = 2*msh%maxedgelength
	
	! write(*,*)	msh%xyz(1,1:100),sum(msh%xyz(1,:))
	! stop

	
		
    if(ptree%MyID==Main_ID)write (*,*) ''
    if(ptree%MyID==Main_ID)write (*,*) 'Maxedge:',Maxedge
	if(ptree%MyID==Main_ID)write (*,*) 'dx:',msh%Delta_ll/2
	if(ptree%MyID==Main_ID)write (*,*) 'wavelength/maxedgelength:',quant%wavelength/msh%maxedgelength

    
    return
    
end subroutine geo_modeling_CURV




subroutine EM_solve_CURV(ho_bf_inv,option,msh,quant,ptree,stats)
    
    use z_HODLR_DEFS
	use APPLICATION_MODULE
	! use RCS_Bi
	! use RCS_Mono
	! use element_vinc
	use z_HODLR_Solve	
    ! use blas95
    implicit none
    
    integer i, j, ii, jj, iii, jjj, ierr
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, n_iter_max, iter, N_unk, N_unk_loc
    real(kind=8) theta, phi, dphi, rcs_V, rcs_H
    real T0
    real(kind=8) n1,n2,rtemp	
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:)
	real(kind=8):: rel_error
	type(z_Hoption)::option
	type(z_mesh)::msh
	type(quant_app)::quant
	type(z_proctree)::ptree
	type(z_hobf)::ho_bf_inv
	type(z_Hstat)::stats	
	complex(kind=8),allocatable:: current(:),voltage(:)
	
	if(option%ErrSol==1)then
		call z_hodlr_test_solve_error(ho_bf_inv,option,ptree,stats)
	endif
	
	
	! if(option%PRECON==DIRECT)then
		! msh%idxs = ho_bf_inv%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		! msh%idxe = ho_bf_inv%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)
	! else 
		! ! write(*,*)associated(ho_bf_for%levels(1)%BP_inverse),'dd' !%matrices_block(1)%N_p),'nima'
		! msh%idxs = ho_bf_for%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,1)
		! msh%idxe = ho_bf_for%levels(1)%BP_inverse(1)%LL(1)%matrices_block(1)%N_p(ptree%MyID - ptree%pgrp(1)%head + 1,2)	
	! endif
	
	! N_unk=msh%Nunk
	N_unk_loc = msh%idxe-msh%idxs+1
	
    if (quant%RCS_static==2) then
    
        phi=180d0
        
        allocate (current(N_unk_loc))
		Current=0
        allocate (voltage(N_unk_loc))
								  
        !$omp parallel do default(shared) private(edge,value_Z)
        do edge=msh%idxs, msh%idxe
            call element_Vinc_VV_CURV(phi,edge,value_Z,msh,quant)
            voltage(edge-msh%idxs+1)=value_Z
        enddo    
        !$omp end parallel do
        
        n1 = OMP_get_wtime()
        
		call z_hodlr_solution(ho_bf_inv,Current,Voltage,N_unk_loc,1,option,ptree,stats)
		
		n2 = OMP_get_wtime()
		stats%Time_Sol = stats%Time_Sol + n2-n1
		call MPI_ALLREDUCE(stats%Time_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
		if(ptree%MyID==Main_ID)then
			write (*,*) ''
			write (*,*) 'Solving:',rtemp,'Seconds'
			write (*,*) ''
		endif
		call MPI_ALLREDUCE(stats%Flop_Sol,rtemp,1,MPI_DOUBLE_PRECISION,MPI_SUM,ptree%Comm,ierr)
		if(ptree%MyID==Main_ID)write (*,'(A13Es14.2)') 'Solve flops:',rtemp	
		
		T0=secnds(0.0)
        call RCS_bistatic_CURV(Current,msh,quant,ptree)

		if(ptree%MyID==Main_ID)then
			write (*,*) ''
			write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
			write (*,*) ''
		endif
    
		deallocate(current)
		deallocate(voltage)
	
    elseif (quant%RCS_static==1) then
    
        allocate (current(N_unk_loc))
        num_sample=quant%RCS_Nsample
        dphi=180./num_sample
        
		allocate (b(N_unk_loc,num_sample+1))
		allocate (x(N_unk_loc,num_sample+1))
		x=0
		
        if(ptree%MyID==Main_ID)open (100, file='RCS_monostatic.txt')

        n1=OMP_get_wtime()       
        do j=0, num_sample 
            phi=j*dphi
			!$omp parallel do default(shared) private(edge,value_Z)
			do edge=msh%idxs, msh%idxe
				call element_Vinc_VV_CURV(phi,edge,value_Z,msh,quant)
				b(edge-msh%idxs+1,j+1)=value_Z
			enddo    
			!$omp end parallel do
		enddo
		
		call z_hodlr_solution(ho_bf_inv,x,b,N_unk_loc,num_sample+1,option,ptree,stats)
			
			
		do j=0, num_sample 	
			phi=j*dphi
			Current=x(:,j+1)
            call RCS_monostatic_VV_CURV(phi,rcs_V,Current,msh,quant,ptree)
!             !$omp parallel do default(shared) private(i)
!             do i=1, N_unk
!                 current(i)=vectors_block(0)%vector(i,2)
!             enddo
!             !$omp end parallel do
!             call RCS_monostatic_HH(theta,phi,rcs_H)
            
            if(ptree%MyID==Main_ID)write (100,*) j,phi,rcs_V !,rcs_H
            
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
		if(ptree%MyID==Main_ID)write (*,'(A13Es14.2)') 'Solve flops:',rtemp	
		
	! call MPI_ALLREDUCE(stats%Time_RedistV,rtemp,1,MPI_DOUBLE_PRECISION,MPI_MAX,ptree%Comm,ierr)
	! if(ptree%MyID==Main_ID)write (*,*) '     Time_RedistV:', rtemp			
		
		deallocate(b)
		deallocate(x)
		deallocate(Current)
		
    endif
        
    return
    
end subroutine EM_solve_CURV






