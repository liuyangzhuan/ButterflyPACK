PROGRAM HODLR_BUTTERFLY_SOLVER_2D
    use MODULE_FILE
	! use geometry_model
	use H_structure
	use cascading_factorization
	use matrices_fill
	use omp_lib
	use MISC
    implicit none

	! include "mkl_vml.fi"	 
	
    real*8 para
    real*8 tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seed_myid(12)
	integer times(8)	
	real*8 t1,t2,x,y,z,r,theta,phi
	complex(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:)
	
	character(len=:),allocatable  :: string
	character(len=1024)  :: strings	
	character(len=6)  :: info_env	
	integer :: length
	integer :: ierr
	integer*8 oldmode,newmode
	type(Hoption)::option	
	type(Hstat)::stats
	type(mesh)::msh
	type(kernelquant)::ker
	type(hobf)::ho_bf,ho_bf_copy
	integer,allocatable:: groupmembers(:)
	integer nmpi
	type(proctree)::ptree
	
	
	
	
	! nmpi and groupmembers should be provided by the user 
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo	
	
	call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	
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
	
	call InitStat(stats)
	
	! time_indexarray = 0
	! time_leastsquare = 0
	! time_buttermul = 0
	! time_buttermulinv = 0
	! time_kernelupdate = 0
	! time_memcopy = 0
	! time_gemm = 0
	! time_gemm1 = 0		   
    ! time_getvec = 0
	! time_resolve = 0
	! time_halfbuttermul = 0
	time_tmp = 0
	
	msh%Origins=(/0d0,0d0,0d0/)
 

     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	CALL getarg(1, DATA_DIR)
	
	para=0.001
	
	! ker%Kernel = EMCURV	
	
	
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
	
	
	
	
	 
	msh%scaling=1d0
	! ker%wavelength=0.0006
	!ker%wavelength=0.0003
	ker%wavelength=0.08

! ker%wavelength=0.08
! Discret=0.05
	ker%RCS_static=1
    ker%RCS_Nsample=2000

	
	
	ker%rank_approximate_para1=6.0
    ker%rank_approximate_para2=6.0
    ker%rank_approximate_para3=6.0
	
	
	option%Nmin_leaf=50
	option%tol_SVD=1d-4
	option%tol_Rdetect=3d-5	
	option%tol_LS=1d-12
	option%tol_itersol=1d-6
	option%N_iter=1000
	option%tol_rand=1d-3
	option%level_check=100
	option%precon=DIRECT
	option%xyzsort=3
	option%LnoBP=40000
	option%TwoLayerOnly=1
	option%touch_para=3
    option%schulzorder=3
    option%schulzlevel=3000
	option%LRlevel=0
	option%ErrFillFull=0
	option%RecLR_leaf='A'
	option%ErrSol=1

	! call MKL_set_num_threads(NUM_Threads)    ! this overwrites omp_set_num_threads for MKL functions 
	
    ! call OMP_set_dynamic(.true.)
    !call OMP_set_nested(.true.)
	

    !*********************************************************

    ker%omiga=2*pi/ker%wavelength/sqrt(mu0*eps0)
    ker%wavenum=2*pi/ker%wavelength

   !***********************************************************************
   if(ptree%MyID==Main_ID)then
   write (*,*) ''
   write (*,*) 'EFIE computing'
   write (*,*) 'wavelength:',ker%wavelength
   write (*,*) ''
   endif
   !***********************************************************************
	
	t1 = OMP_get_wtime()
    if(ptree%MyID==Main_ID)write(*,*) "geometry modeling......"
    call geo_modeling_CURV(msh,ker,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "modeling finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1

	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "constructing H_matrices formatting......"
    call H_matrix_structuring(ho_bf,para,option,msh,ptree)
	call BPlus_structuring(ho_bf,option,msh,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices formatting finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1
	! stop 
	
    !pause
    
    !call compression_test()
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices filling......"
    ! call matrices_filling(ho_bf,option,stats,msh,ker,element_Zmn_FULL,ptree)
    call matrices_filling(ho_bf,option,stats,msh,ker,element_Zmn_EMCURV,ptree)
	! if(option%precon/=DIRECT)then
		call copy_HOBF(ho_bf,ho_bf_copy)	
	! end if
    if(ptree%MyID==Main_ID)write(*,*) "H_matrices filling finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	! write(*,*)t2-t1
	
	if(option%precon/=NOPRECON)then
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing......"
    call cascading_factorizing(ho_bf,option,stats,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	end if

    if(ptree%MyID==Main_ID)write(*,*) "EM_solve......"
    call EM_solve_CURV(ho_bf_copy,ho_bf,option,msh,ker,ptree,stats)
    if(ptree%MyID==Main_ID)write(*,*) "EM_solve finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "	
	
	
    if(ptree%MyID==Main_ID)write(*,*) "-------------------------------program end-------------------------------------"

    ! ! ! ! pause

end PROGRAM HODLR_BUTTERFLY_SOLVER_2D





subroutine geo_modeling_CURV(msh,ker,ptree)

    use MODULE_FILE
    implicit none
    type(mesh)::msh
	type(kernelquant)::ker
    integer i,j,ii,jj,iii,jjj
    integer intemp
    integer node, patch, edge, flag
    integer node1, node2, num_node,Maxedge
    integer node_temp(2)
    real(kind=8) dx, xx, yy, rr, theta,L,M,Am,tt,L1,L2,L3
    
    real(kind=8),allocatable :: node_xy_original(:,:)
    integer,allocatable :: num_edge_of_node(:)
    
    real(kind=8) a(3),b(3),c(3),r0, phi_start
	type(proctree)::ptree
	
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
		M = 0.2d0*ker%wavelength
		L = 1.5d0*ker%wavelength
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
		M = 0.2d0*ker%wavelength
		L = 1.5d0*ker%wavelength
		
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
	if(ptree%MyID==Main_ID)write (*,*) 'wavelength/maxedgelength:',ker%wavelength/msh%maxedgelength

    
    return
    
end subroutine geo_modeling_CURV




subroutine EM_solve_CURV(ho_bf_for,ho_bf_inv,option,msh,ker,ptree,stats)
    
    use MODULE_FILE
	use RCS_Bi
	use RCS_Mono
	use element_vinc
	use HODLR_Solve	
    ! use blas95
    implicit none
    
    integer i, j, ii, jj, iii, jjj, ierr
    integer level, blocks, edge, patch, node, group
    integer rank, index_near, m, n, length, flag, num_sample, N_iter_max, iter, N_unk, N_unk_loc
    real*8 theta, phi, dphi, rcs_V, rcs_H
    real T0
    real*8 n1,n2,rtemp	
    complex(kind=8) value_Z
    complex(kind=8),allocatable:: Voltage_pre(:),x(:,:),b(:,:)
	real*8:: rel_error
	type(Hoption)::option
	type(mesh)::msh
	type(kernelquant)::ker
	type(proctree)::ptree
	type(hobf)::ho_bf_for,ho_bf_inv
	type(Hstat)::stats	
	complex(kind=8),allocatable:: current(:),voltage(:)
	
	if(option%ErrSol==1)then
		call HODLR_Test_Solve_error(ho_bf_for,ho_bf_inv,option,ptree,stats)
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
	
    if (ker%RCS_static==2) then
    
        phi=180d0
        
        allocate (current(N_unk_loc))
		Current=0
        allocate (voltage(N_unk_loc))
								  
        !$omp parallel do default(shared) private(edge,value_Z)
        do edge=msh%idxs, msh%idxe
            call element_Vinc_VV_CURV(phi,edge,value_Z,msh,ker)
            voltage(edge-msh%idxs+1)=value_Z
        enddo    
        !$omp end parallel do
        
        n1 = OMP_get_wtime()
        
		call HODLR_Solution(ho_bf_for,ho_bf_inv,Current,Voltage,N_unk_loc,1,option,ptree,stats)
		
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
        call RCS_bistatic_CURV(Current,msh,ker,ptree)

		if(ptree%MyID==Main_ID)then
			write (*,*) ''
			write (*,*) 'Bistatic RCS',secnds(T0),'Seconds'
			write (*,*) ''
		endif
    
		deallocate(current)
		deallocate(voltage)
	
    elseif (ker%RCS_static==1) then
    
        allocate (current(N_unk_loc))
        num_sample=ker%RCS_Nsample
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
				call element_Vinc_VV_CURV(phi,edge,value_Z,msh,ker)
				b(edge-msh%idxs+1,j+1)=value_Z
			enddo    
			!$omp end parallel do
		enddo
		
		call HODLR_Solution(ho_bf_for,ho_bf_inv,x,b,N_unk_loc,num_sample+1,option,ptree,stats)
			
			
		do j=0, num_sample 	
			phi=j*dphi
			Current=x(:,j+1)
            call RCS_monostatic_VV_CURV(phi,rcs_V,Current,msh,ker,ptree)
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

