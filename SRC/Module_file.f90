module MODULE_FILE
	implicit none
    INCLUDE 'mpif.h'   
	
	integer,parameter:: LplusMax=3
	real*8,parameter :: pi = 4d0*atan(1d0)
	real*8,parameter :: cd = 299792458d0
	complex(kind=8),parameter :: junit=(0d0,1d0)
	real*8,parameter :: Bigvalue=1d300
	real(kind=8),parameter:: SafeUnderflow=1D-150
	real(kind=8),parameter:: SafeAbsoulte=1D-14
	complex(kind=8),parameter :: cone=CMPLX(1d0,0d0,kind=8)
	complex(kind=8),parameter :: czero=CMPLX(0d0,0d0,kind=8)
	
	real*8,parameter :: eps0=1d7/(4d0*pi*cd**2)
    real*8,parameter :: mu0=pi*4d-7
    real*8,parameter :: gamma=1.781072418d0
    real*8,parameter :: impedence0=sqrt(mu0/eps0)
	
	 integer,parameter:: EMCURV=1
	 integer,parameter:: EMSURF=2
	 integer,parameter:: RBF=3
	 integer,parameter:: FULL=4	
	
	 integer,parameter:: NOPRECON=2
	 integer,parameter:: SAIPRECON=3	
	
	 integer,parameter:: DIRECT=1
	 integer,parameter :: Main_ID = 0
	 
	 integer,parameter :: nbslpk=16
	
	
	 
	type grid
		integer :: nsprow,nspcol ! number of process rows/columns as square as possible, it's designed that nspcol>=nsprow
		integer :: hprow,hpcol ! head process in the row/column process tree		
		integer :: ctxt ! blacs context
		integer :: Comm ! MPI communicator for this grid
		integer :: gprow,gpcol ! the group number in the row and column dimension, no longer needed once constructed 
		type(grid),pointer::gdc(:)=>null() ! pointer to its two children
	end type grid
	
	
	! type parACAblock
		 ! integer headm,headn
		 ! integer tailm,tailn
		 ! complex(kind=8),allocatable :: Umat(:,:),Vmat(:,:)
		 ! type(grid),pointer::gd=>null()
		 ! type(parACAblock),pointer::bk1=>null(), bk2=>null()
     ! end type parACAblock
	
	
	type procgroup
		!****Define MPI variables
		integer :: nprow,npcol,nproc ! number of processors
		integer :: head,tail
		integer :: ctxt ! blacs context 
		integer :: ctxt1D ! blacs context 1D Row noncyclic
		integer :: ctxt_head ! blacs context only involving the head process 
		integer :: Comm ! MPI communicator
		type(grid),pointer::gd=>null() ! the hierarchical process grid structure associated with each process group
	end type procgroup
	
	type proctree
		integer nlevel
		integer Comm
		integer nproc
		integer :: MyID
		type(procgroup),allocatable::pgrp(:)
	end type proctree
	
	type commquant1D
		integer offset ! offset in my local array
		integer size ! size of the message along first dimension 
		complex(kind=8),allocatable::dat(:,:) ! communication buffer
	end type commquant1D
	

     type basisgroup
		 integer pgno
         integer head
         integer tail
         integer level
         real*8 radius
		 real*8:: boundary(2)
         real*8,allocatable:: center(:)
     end type basisgroup

     type vectorsblock
         integer style
         integer head
         integer tail
         complex(kind=8),allocatable :: vector(:,:)
     end type vectorsblock
     
     type vectorset
         complex(kind=8),allocatable :: vector(:)
     end type vectorset

     type butterfly_vector
         complex(kind=8),allocatable :: vector(:)
     end type butterfly_vector

     type butterflymatrix
         complex(kind=8),allocatable :: matrix(:,:)
         complex(kind=8),allocatable :: null(:,:)
		 integer,allocatable :: list(:) 		 
		 integer::nulldim,mdim,ndim		 
     end type butterflymatrix
     
     type butterfly_col_select
         integer,allocatable :: select_columns(:)
     end type butterfly_col_select
     
     type butterfly_Kerl
         integer num_col
         integer num_row
         type(butterflymatrix),allocatable :: blocks(:,:)
     end type butterfly_Kerl

     type butterfly_UV
         integer num_blk
         type(butterflymatrix),allocatable :: blocks(:)
     end type butterfly_UV
     
     type matrixblock
		 integer pgno,pgno_db
         integer level
         integer col_group
         integer row_group
         integer nested_num
         integer style
         integer data_type
         integer level_butterfly
		 integer rankmax,rankmin
		 integer dimension_rank
		 integer M,N
		 integer M_loc,N_loc,M_loc_db,N_loc_db
		 integer headm,headn
         integer,pointer:: M_p(:,:)=>null()
         integer,pointer:: N_p(:,:)=>null()
         integer,pointer:: M_p_db(:,:)=>null()
         integer,pointer:: N_p_db(:,:)=>null()
		 complex(kind=8),allocatable :: fullmat(:,:)
         complex(kind=8),allocatable:: KerInv(:,:)	
		 type(butterfly_UV) :: ButterflyU
         type(butterfly_UV) :: ButterflyV
         type(butterflymatrix),allocatable :: ButterflyMiddle(:,:)		 
         type(butterfly_Kerl),allocatable :: ButterflyKerl(:)
         type(butterfly_col_select),allocatable :: ButterflyColSelect(:,:)
     end type matrixblock
	 
	 
	 type onelplus
		integer Nbound
		integer rankmax
		type(matrixblock),pointer:: matrices_block(:)=>null()
		integer,allocatable::boundary_map(:)
	 end type onelplus
	 
	 
	 type blockplus
         integer level
         integer col_group
         integer row_group
		 integer pgno
         integer Lplus
		 real*8:: boundary(2)
		 type(onelplus),pointer:: LL(:)=>null() 
	 end type blockplus
	 

	 type schulz_operand
		type(matrixblock):: matrices_block
		real*8::A2norm,scale
		real*8,allocatable::diags(:)
		integer order
	 end type schulz_operand
	
	 type cascadingfactors
		integer level
		integer N_block_forward
		integer N_block_inverse
		integer Bidxs,Bidxe   ! indice range of my local groups
		! type(matrixblock),pointer:: matrices_block(:)	
		! type(matrixblock),pointer:: matrices_block_inverse(:)		
		! type(matrixblock),pointer:: matrices_block_inverse_schur(:)		
		type(blockplus),pointer:: BP(:)=>null()
		type(blockplus),pointer:: BP_inverse(:)=>null()
		type(blockplus),pointer:: BP_inverse_schur(:)=>null()
	 end type cascadingfactors

	type hobf
		integer Maxlevel,N
		integer ind_lv,ind_bk
		type(cascadingfactors),allocatable::levels(:)	 
	end type hobf 
	 
	 type partitionedblocks
		integer level
		type(matrixblock),pointer:: blocks_A=>null()
		type(matrixblock),pointer:: blocks_B=>null()
		type(matrixblock),pointer:: blocks_C=>null()
		type(matrixblock),pointer:: blocks_D=>null()
	 end type partitionedblocks

     type RandomBlock
         integer level_butterfly
         type(butterfly_Kerl), allocatable :: RandomVectorRR(:)
         type(butterfly_Kerl), allocatable :: RandomVectorLL(:)
     end type RandomBlock
	 
	 type Hoption
		real*8 tol_Rdetect
		real*8 tol_SVD
		real*8 tol_LS
		real*8 tol_rand
		real*8 tol_itersol
		real*8 touch_para
		
		integer N_iter
		integer precon
		integer:: TwoLayerOnly 
		integer:: LnoBP 
		integer xyzsort
		integer::level_check
		integer::Nmin_leaf
		integer::schulzorder
		integer::schulzlevel
		integer::LRlevel
	 end type Hoption
	 
	type Hstat
		real*8 Time_random(5)  ! Intialization, MVP, Reconstruction, Reconstruction of one subblock 
		real*8 Time_Sblock,Time_Inv,Time_SMW,Time_Fill,Time_RedistB,Time_RedistV,Time_Sol
		real*8 Mem_peak,Mem_Sblock,Mem_SMW,Mem_Direct,Mem_int_vec,Mem_For
		real*8 Flop_Fill, Flop_Factor, Flop_Sol, Flop_Tmp
		integer, allocatable:: rankmax_of_level(:),rankmin_of_level(:)
	 end type Hstat	

	type Mesh
		integer integral_points
		integer maxpatch
		integer maxnode
		integer maxedge
		integer model2d
		integer Nunk
		integer idxs,idxe  ! range of local indices after reordering
		integer mesh_normal	 !
		real*8 Delta_ll		
		real*8 scaling
		real*8, allocatable :: ng1(:),ng2(:),ng3(:),gauss_w(:)
		real*8,allocatable:: normal_of_patch(:,:) ! normal vector of each patch
		integer,allocatable:: node_of_patch(:,:) ! vertices of each patch
		real*8,allocatable:: xyz(:,:)   ! coordinates of the points
		real*8:: Origins(3) ! only used for spherical coordinate transform
		real*8 maxedgelength,minedgelength
		integer::Ncorner
		real*8,allocatable::corner_points(:,:)
		real*8::corner_radius		
					! for 2D mesh: 0 point to coordinates of each edge center (unknown x), 1-2 point to coordinates of each edge vertice  
					! for 3D mesh: 0 point to coordinates of each edge center (unknown x), 1-2 point to coordinates of each edge vertice, 3-4 point to two patches that share the same edge, 5-6 point to coordinates of last vertice of the two patches
					! for general: 0 point to coordinates of each unknown x
		integer,allocatable:: info_unk(:,:)  
	end type Mesh
	 
	 
	 
	type kernelquant
		integer Kernel
		real*8 wavenum
		real*8 wavelength
		real*8 omiga 
		real*8 rank_approximate_para1, rank_approximate_para2, rank_approximate_para3 ! rank estimation parameter
		integer RCS_static  ! monostatic or bistatic RCS
		integer RCS_Nsample ! number of RCS samples
		real*8:: CFIE_alpha ! combination parameter in CFIE
		real*8 sigma, lambda ! parameters in ker%Kernel matrices		
		
		complex(kind=8), allocatable :: matZ_glo(:,:)
	end type kernelquant
	
	
	abstract interface
		subroutine HOBF_MVP_blk(trans,N,num_vect_sub,Vin,Vout,operand)
			implicit none
			character trans
			integer M, N, num_vect_sub
			complex(kind=8) :: Vin(:,:), Vout(:,:)
			class(*)::operand
		end subroutine HOBF_MVP_blk

		subroutine BF_MVP_blk(operand,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,ptree,stats,operand1)
			import :: matrixblock,proctree,Hstat
			implicit none
			class(*)::operand	
			class(*),optional::operand1	
			type(matrixblock)::block_o
			character trans
			integer M, N, num_vect_sub
			type(proctree)::ptree
			type(Hstat)::stats
			complex(kind=8) :: Vin(:,:), Vout(:,:),a,b
		end subroutine BF_MVP_blk

		subroutine Z_elem(edge_m,edge_n,value,msh,ker)
			import::mesh,kernelquant
			implicit none
			integer edge_m, edge_n
			complex(kind=8) value
			type(mesh)::msh
			type(kernelquant)::ker
		end subroutine Z_elem		
	end interface
	

	 real*8::time_tmp
     type(basisgroup),allocatable:: basis_group(:)
	 integer vecCNT
	 integer,allocatable:: new2old(:)
	 integer,allocatable:: basis_group_pre(:,:)
	 CHARACTER (LEN=1000) DATA_DIR	
	 



	
	 
end module

