PROGRAM HODLR_BUTTERFLY_SOLVER_2D
    use z_HODLR_DEFS
    use EMCURV_MODULE
	
	use z_HODLR_structure
	use z_HODLR_factor
	use z_HODLR_constr
	use omp_lib
	use z_misc
    implicit none

	! include "mkl_vml.fi"	 
	
    real(kind=8) para
    real(kind=8) tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seed_myid(50)
	integer times(8)	
	integer edge
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
	type(quant_EMCURV),target::quant
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

	
 	! register the user-defined function and type in ker 
	ker%FuncZmn=>Z_elem_EMCURV
	ker%QuantZmn=>quant

     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	CALL getarg(1, DATA_DIR)
	
	
	quant%model2d=10 ! Model type (1=strip; 2=corner reflector; 3=two opposite strips; 4=CR with RRS; 5=cylinder; 6=Rectangle Cavity); 7=half cylinder; 8=corrugated half cylinder; 9=corrugated corner reflector; 10=open polygon; 11=taller open polygon 
	! msh%Nunk=1280000
	! msh%Nunk=320000
	! msh%Nunk=80000
	! msh%Nunk=20000
	 msh%Nunk=5000
    ! msh%Nunk=160000	
	! Refined_level=0
	quant%Nunk = msh%Nunk
	
	
	
	! ker%Kernel = FULL	
	! allocate(ker%matZ_glo(msh%Nunk,msh%Nunk))
	! call RandomMat(msh%Nunk,msh%Nunk,msh%Nunk,ker%matZ_glo,0)
	! call MPI_Bcast(ker%matZ_glo,msh%Nunk*msh%Nunk,MPI_DOUBLE_COMPLEX,0,MPI_Comm_World,ierr)
	
	
	
	
	option%nogeo=0 
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
    call geo_modeling_CURV(quant,ptree)
	
	! generate the list of points for clustering
	allocate(msh%xyz(2,quant%Nunk))
	do edge=1, quant%Nunk
		msh%xyz(:,edge) = quant%xyz(:,edge*2-1)
	enddo	
    option%touch_para = 3* quant%minedgelength
	
	if(ptree%MyID==Main_ID)write(*,*) "modeling finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1

	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "constructing HODLR formatting......"
    call z_HODLR_structuring(ho_bf,option,msh,ker,z_element_Zmn_user,ptree)
	call z_BPlus_structuring(ho_bf,option,msh,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "HODLR formatting finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	! write(*,*)t2-t1
	! stop 
	
    !pause
    
    !call compression_test()
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "HODLR construction......"
    ! call HODLR_construction(ho_bf,option,stats,msh,ker,element_Zmn_FULL,ptree)
    call z_HODLR_construction(ho_bf,option,stats,msh,ker,z_element_Zmn_user,ptree)
	! if(option%precon/=DIRECT)then
		! call copy_HOBF(ho_bf,ho_bf_copy)	
	! end if
    if(ptree%MyID==Main_ID)write(*,*) "HODLR construction finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
 	t2 = OMP_get_wtime()   
	! write(*,*)t2-t1
	
	if(option%precon/=NOPRECON)then
    if(ptree%MyID==Main_ID)write(*,*) "Cascading factorizing......"
    call z_HODLR_Factorization(ho_bf,option,stats,ptree)
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





subroutine geo_modeling_CURV(quant,ptree)

    use z_HODLR_DEFS
	use EMCURV_MODULE
    implicit none
    ! type(z_mesh)::msh
	type(quant_EMCURV)::quant
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
	
	Maxedge	= quant%Nunk

	quant%Ncorner = 0
	
    if (quant%model2d==1) then  !*****single strip*****
        
        quant%Delta_ll=2d0/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        quant%xyz(1,0)=-1d0 ; quant%xyz(2,0)=0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, quant%maxnode-1
            dx=node*quant%Delta_ll/2
            quant%xyz(1,node)=dx-1d0
            quant%xyz(2,node)=0
        enddo
        !$omp end parallel do
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1   ! 1 start 0 center 2 end
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo           
        
    elseif (quant%model2d==2) then  !*****corner reflector*****
    
        quant%Delta_ll=2d0/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        
        quant%xyz(1,0)=-1d0/sqrt(2d0) ; quant%xyz(2,0)=1d0/sqrt(2d0)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*quant%Delta_ll/2
            quant%xyz(1,node)=(dx-1d0)/sqrt(2d0)
            quant%xyz(2,node)=(1d0-dx)/sqrt(2d0)
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*quant%Delta_ll/2
            quant%xyz(1,node+Maxedge)=dx/sqrt(2d0)
            quant%xyz(2,node+Maxedge)=dx/sqrt(2d0)
        enddo
        !$omp end parallel do
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo
		
		quant%Ncorner = 1
		allocate(quant%corner_points(2,quant%Ncorner))
		quant%corner_points(1,1) = 0
		quant%corner_points(2,1) = 0        
    
    elseif (quant%model2d==3) then  !*****two opposite strips*****
    
        quant%Delta_ll=2d0/Maxedge
        quant%maxnode=2*Maxedge+2
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        
        quant%xyz(1,0)=-0.5d0 ; quant%xyz(2,0)=1.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*quant%Delta_ll/2
            quant%xyz(1,node)=dx-0.5d0
            quant%xyz(2,node)=1.0d0
        enddo
        !$omp end parallel do
        quant%xyz(1,Maxedge+1)=-0.5d0 ; quant%xyz(2,Maxedge+1)=-1.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*quant%Delta_ll/2
            quant%xyz(1,node+Maxedge+1)=dx-0.5d0
            quant%xyz(2,node+Maxedge+1)=-1.0d0
        enddo
        !$omp end parallel do
        
        intemp=int(Maxedge/2)
        
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, intemp
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo
        
        quant%info_unk(1,intemp+1)=quant%maxnode/2 ; quant%info_unk(2,intemp+1)=quant%maxnode/2+2 ; quant%info_unk(0,intemp+1)=quant%maxnode/2+1
        do edge=intemp+2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo
    
    
    elseif (quant%model2d==5) then !************cylinder*****************
        
        quant%Delta_ll=  2.0d0*pi/Maxedge
        quant%maxnode=2*Maxedge
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        quant%xyz(1,0)=1.0d0 ; quant%xyz(2,0)=0.0d0
        !$omp parallel do default(shared) private(node,dx)
        do node=1, quant%maxnode-1
            dx=node*quant%Delta_ll/2
            quant%xyz(1,node)=cos(dx)
            quant%xyz(2,node)=sin(dx)
        enddo
        !$omp end parallel do
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge-1
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo     
        quant%info_unk(1,Maxedge)=quant%info_unk(2,Maxedge-1)
        quant%info_unk(2,Maxedge)=0
        quant%info_unk(0,Maxedge)=quant%info_unk(0,Maxedge-1)+2      
        
    elseif (quant%model2d==6) then !***********cavity******************
    
        quant%Delta_ll=11d0/3d0/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        dx=quant%Delta_ll/2d0
        quant%xyz(1,0)=-1d0/6. ; quant%xyz(2,0)=1.0d0
        node=0
        do while (flag==0) 
            node=node+1
            quant%xyz(1,node)=quant%xyz(1,node-1)-dx
            quant%xyz(2,node)=quant%xyz(2,node-1)
            if (quant%xyz(1,node)<=-0.5d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            quant%xyz(1,node)=quant%xyz(1,node-1)
            quant%xyz(2,node)=quant%xyz(2,node-1)-dx
            if (quant%xyz(2,node)<=0.) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            quant%xyz(1,node)=quant%xyz(1,node-1)+dx
            quant%xyz(2,node)=quant%xyz(2,node-1)
            if (quant%xyz(1,node)>=0.5d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            quant%xyz(1,node)=quant%xyz(1,node-1)
            quant%xyz(2,node)=quant%xyz(2,node-1)+dx
            if (quant%xyz(2,node)>=1.0d0) then
                flag=1
            endif
        enddo
        flag=0
        do while (flag==0) 
            node=node+1
            quant%xyz(1,node)=quant%xyz(1,node-1)-dx
            quant%xyz(2,node)=quant%xyz(2,node-1)
            if (node>=quant%maxnode-1) then
                flag=1
            endif
        enddo
        
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo

		quant%Ncorner = 4
		allocate(quant%corner_points(2,quant%Ncorner))
		quant%corner_points(1,1) = -0.5
		quant%corner_points(2,1) = 0     
		quant%corner_points(1,2) = 0.5
		quant%corner_points(2,2) = 0
		quant%corner_points(1,3) = 0.5
		quant%corner_points(2,3) = 1
		quant%corner_points(1,4) = -0.5
		quant%corner_points(2,4) = 1		
		
    elseif (quant%model2d==7) then   !************open cylinder*****************
    
	
	
        quant%Delta_ll=1d0*pi/Maxedge !2.0d0*pi*5d0/6d0/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        
        quant%xyz(1,0)=cos(0*pi) ; quant%xyz(2,0)=sin(0*pi)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*quant%Delta_ll
            quant%xyz(1,node*2)=cos(0*pi+dx)
            quant%xyz(2,node*2)=sin(0*pi+dx)
        enddo
        !$omp end parallel do
		
		
        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            quant%xyz(1,node*2-1)=(quant%xyz(1,node*2-2)+quant%xyz(1,node*2))/2d0
            quant%xyz(2,node*2-1)=(quant%xyz(2,node*2-2)+quant%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do		
		
		
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo           

    elseif (quant%model2d==8) then   !************corrugated open cylinder*****************
		M = 0.2d0*quant%wavelength
		L = 1.5d0*quant%wavelength
        phi_start = 3d0/2d0*pi
		
		quant%Delta_ll=1d0*pi/Maxedge !2.0d0*pi*5d0/6d0/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        
        quant%xyz(1,0)=cos(phi_start) ; quant%xyz(2,0)=sin(phi_start)
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge
            dx=node*quant%Delta_ll
            quant%xyz(1,node*2)=(1+M*sin(2*pi*dx/L))*cos(phi_start+dx)
            quant%xyz(2,node*2)=(1+M*sin(2*pi*dx/L))*sin(phi_start+dx)
        enddo
        !$omp end parallel do
		
		
        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            quant%xyz(1,node*2-1)=(quant%xyz(1,node*2-2)+quant%xyz(1,node*2))/2d0
            quant%xyz(2,node*2-1)=(quant%xyz(2,node*2-2)+quant%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do		
		
		
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo    

		
    elseif (quant%model2d==9) then  !*****corrugated corner reflector*****
		M = 0.2d0*quant%wavelength
		L = 1.5d0*quant%wavelength
		
        quant%Delta_ll=2d0/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        

		Am = M*sin(pi/2-2*pi/L)-M		
		quant%xyz(1,0)=-1d0/sqrt(2d0)+Am/sqrt(2d0); quant%xyz(2,0)=1d0/sqrt(2d0)+Am/sqrt(2d0)
		
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge/2
            dx=node*quant%Delta_ll
			Am = M*sin(2*pi*dx/L+pi/2-2*pi/L)-M
			
            quant%xyz(1,node*2)=(dx-1d0)/sqrt(2d0)+Am/sqrt(2d0)
            quant%xyz(2,node*2)=(1d0-dx)/sqrt(2d0)+Am/sqrt(2d0)
        enddo
        !$omp end parallel do
        !$omp parallel do default(shared) private(node,dx)
        do node=1, Maxedge/2
            dx=node*quant%Delta_ll
			Am = M*sin(2*pi*(dx+1)/L+pi/2-2*pi/L)-M
			
            quant%xyz(1,node*2+Maxedge)=dx/sqrt(2d0)-Am/sqrt(2d0)
            quant%xyz(2,node*2+Maxedge)=dx/sqrt(2d0)+Am/sqrt(2d0)
        enddo
        !$omp end parallel do

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            quant%xyz(1,node*2-1)=(quant%xyz(1,node*2-2)+quant%xyz(1,node*2))/2d0
            quant%xyz(2,node*2-1)=(quant%xyz(2,node*2-2)+quant%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo  		

		quant%Ncorner = 1
		allocate(quant%corner_points(2,quant%Ncorner))
		quant%corner_points(1,1) = 0
		quant%corner_points(2,1) = 0
		
    elseif (quant%model2d==10) then  !*****cup*****

		L1 = 0
		L2 = sqrt(2d0)
		L3 = sqrt(2d0)
		L = 3
        quant%Delta_ll=(2*L1+2*L2+2*L3+L)/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        
		
        !$omp parallel do default(shared) private(node,dx,tt)
        do node=0, Maxedge
            dx=node*quant%Delta_ll
			if(dx<=L1)then
				tt = dx
				quant%xyz(1,node*2)= -(L/2-(L1-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L1+L2+L3)/sqrt(2d0)-tt/sqrt(2d0)			
			else if(dx<=L1+L2)then
				tt = dx - L1
				quant%xyz(1,node*2)= -(L/2-(-L2+L3)/sqrt(2d0))+tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L2+L3)/sqrt(2d0)-tt/sqrt(2d0)
			else if(dx<=L1+L2+L3)then
				tt = dx - (L1+L2)
				quant%xyz(1,node*2)= -(L/2-L3/sqrt(2d0))-tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L3)/sqrt(2d0)-tt/sqrt(2d0)				
			else if(dx<=L1+L2+L3+L)then
				tt = dx - (L1+L2+L3)
				quant%xyz(1,node*2)= -L/2+tt
				quant%xyz(2,node*2)= 0
			else if(dx<=L1+L2+L3+L+L3)then
				tt = dx - (L1+L2+L3+L)
				quant%xyz(1,node*2)= L/2-tt/sqrt(2d0)
				quant%xyz(2,node*2)= tt/sqrt(2d0)		
			else if(dx<=L1+L2+L3+L+L3+L2)then
				tt = dx - (L1+L2+L3+L+L3)
				quant%xyz(1,node*2)= (L/2-L3/sqrt(2d0))+tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L3)/sqrt(2d0)+tt/sqrt(2d0)				
			else
				tt = dx - (L1+L2+L3+L+L3+L2)
				quant%xyz(1,node*2)= (L/2-(-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L3+L2)/sqrt(2d0)+tt/sqrt(2d0)							
			end if 
        enddo
        !$omp end parallel do			

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            quant%xyz(1,node*2-1)=(quant%xyz(1,node*2-2)+quant%xyz(1,node*2))/2d0
            quant%xyz(2,node*2-1)=(quant%xyz(2,node*2-2)+quant%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo  		

		! quant%Ncorner = 4
		! allocate(quant%corner_points(2,quant%Ncorner))
		! quant%corner_points(1,1) = -L/2
		! quant%corner_points(2,1) = 0
		! quant%corner_points(1,2) = L/2
		! quant%corner_points(2,2) = 0
		! quant%corner_points(1,3) = -(L/2-L3/sqrt(2d0))
		! quant%corner_points(2,3) = (L3)/sqrt(2d0)		
		! quant%corner_points(1,4) = (L/2-L3/sqrt(2d0))
		! quant%corner_points(2,4) = (L3)/sqrt(2d0)			

    elseif (quant%model2d==11) then  !*****longer cup*****
		L1 = 1
		L2 = sqrt(2d0)
		L3 = sqrt(2d0)
		L = 3.5
        quant%Delta_ll=(2*L1+2*L2+2*L3+L)/Maxedge
        quant%maxnode=2*Maxedge+1
        allocate (quant%xyz(2,0:quant%maxnode-1), quant%info_unk(0:2,Maxedge))
        
		
        !$omp parallel do default(shared) private(node,dx,tt)
        do node=0, Maxedge
            dx=node*quant%Delta_ll
			if(dx<=L1)then
				tt = dx
				quant%xyz(1,node*2)= -(L/2-(L1-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L1+L2+L3)/sqrt(2d0)-tt/sqrt(2d0)			
			else if(dx<=L1+L2)then
				tt = dx - L1
				quant%xyz(1,node*2)= -(L/2-(-L2+L3)/sqrt(2d0))+tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L2+L3)/sqrt(2d0)-tt/sqrt(2d0)
			else if(dx<=L1+L2+L3)then
				tt = dx - (L1+L2)
				quant%xyz(1,node*2)= -(L/2-L3/sqrt(2d0))-tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L3)/sqrt(2d0)-tt/sqrt(2d0)				
			else if(dx<=L1+L2+L3+L)then
				tt = dx - (L1+L2+L3)
				quant%xyz(1,node*2)= -L/2+tt
				quant%xyz(2,node*2)= 0
			else if(dx<=L1+L2+L3+L+L3)then
				tt = dx - (L1+L2+L3+L)
				quant%xyz(1,node*2)= L/2-tt/sqrt(2d0)
				quant%xyz(2,node*2)= tt/sqrt(2d0)		
			else if(dx<=L1+L2+L3+L+L3+L2)then
				tt = dx - (L1+L2+L3+L+L3)
				quant%xyz(1,node*2)= (L/2-L3/sqrt(2d0))+tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L3)/sqrt(2d0)+tt/sqrt(2d0)				
			else
				tt = dx - (L1+L2+L3+L+L3+L2)
				quant%xyz(1,node*2)= (L/2-(-L2+L3)/sqrt(2d0))-tt/sqrt(2d0)
				quant%xyz(2,node*2)= (L3+L2)/sqrt(2d0)+tt/sqrt(2d0)							
			end if 
        enddo
        !$omp end parallel do			

        !$omp parallel do default(shared) private(node)
        do node=1, Maxedge
            quant%xyz(1,node*2-1)=(quant%xyz(1,node*2-2)+quant%xyz(1,node*2))/2d0
            quant%xyz(2,node*2-1)=(quant%xyz(2,node*2-2)+quant%xyz(2,node*2))/2d0
        enddo
        !$omp end parallel do			
		
        quant%info_unk(1,1)=0 ; quant%info_unk(2,1)=2 ; quant%info_unk(0,1)=1
        do edge=2, Maxedge
            quant%info_unk(1,edge)=quant%info_unk(2,edge-1)
            quant%info_unk(2,edge)=quant%info_unk(2,edge-1)+2
            quant%info_unk(0,edge)=quant%info_unk(0,edge-1)+2
        enddo  		

		quant%Ncorner = 6
		allocate(quant%corner_points(2,quant%Ncorner))
		quant%corner_points(1,1) = -(L/2-(-L2+L3)/sqrt(2d0))
		quant%corner_points(2,1) = (L2+L3)/sqrt(2d0)
		quant%corner_points(1,2) = (L/2-(-L2+L3)/sqrt(2d0))
		quant%corner_points(2,2) = (L2+L3)/sqrt(2d0)	
		quant%corner_points(1,3) = -L/2
		quant%corner_points(2,3) = 0
		quant%corner_points(1,4) = L/2
		quant%corner_points(2,4) = 0
		quant%corner_points(1,5) = -(L/2-L3/sqrt(2d0))
		quant%corner_points(2,5) = (L3)/sqrt(2d0)		
		quant%corner_points(1,6) = (L/2-L3/sqrt(2d0))
		quant%corner_points(2,6) = (L3)/sqrt(2d0)		
    endif
    
	quant%maxedgelength = 0
	do edge=1,Maxedge
		quant%maxedgelength = max(quant%maxedgelength,sqrt(sum(abs(quant%xyz(:,quant%info_unk(1,edge))-quant%xyz(:,quant%info_unk(2,edge))))**2))
	end do	
	
	quant%minedgelength = 10000000
	do edge=1,Maxedge
		quant%minedgelength = min(quant%minedgelength,sqrt(sum(abs(quant%xyz(:,quant%info_unk(1,edge))-quant%xyz(:,quant%info_unk(2,edge)))**2)))
	end do		
	
	quant%corner_radius = 2*quant%maxedgelength
	
	! write(*,*)	quant%xyz(1,1:100),sum(quant%xyz(1,:))
	! stop
		
    if(ptree%MyID==Main_ID)write (*,*) ''
    if(ptree%MyID==Main_ID)write (*,*) 'Maxedge:',Maxedge
	if(ptree%MyID==Main_ID)write (*,*) 'dx:',quant%Delta_ll/2
	if(ptree%MyID==Main_ID)write (*,*) 'wavelength/maxedgelength:',quant%wavelength/quant%maxedgelength

    
    return
    
end subroutine geo_modeling_CURV




subroutine EM_solve_CURV(ho_bf_inv,option,msh,quant,ptree,stats)
    
    use z_HODLR_DEFS
	use EMCURV_MODULE
	! use RCS_Bi
	! use RCS_Mono
	! use element_vinc
	use z_HODLR_Solve_Mul	
    
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
	type(quant_EMCURV)::quant
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
            call element_Vinc_VV_CURV(phi,msh%new2old(edge),value_Z,msh,quant)
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
				call element_Vinc_VV_CURV(phi,msh%new2old(edge),value_Z,msh,quant)
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






