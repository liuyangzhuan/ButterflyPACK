module MODULE_FILE
	use iso_c_binding    
	implicit none
    INCLUDE 'mpif.h'   
	
	!**** common parameters
	integer,parameter::dp=kind(0.0d0),sp=kind(0.0)
	real*8,parameter :: pi = 4d0*atan(1d0)
	complex(kind=8),parameter :: junit=(0d0,1d0)
	real*8,parameter :: Bigvalue=1d300
	real(kind=8),parameter:: SafeUnderflow=1D-30
	real(kind=8),parameter:: SafeAbsoulte=1D-14
	complex(kind=8),parameter :: cone=CMPLX(1d0,0d0,kind=8)
	complex(kind=8),parameter :: czero=CMPLX(0d0,0d0,kind=8)
    integer,parameter :: Main_ID = 0 ! Head MPI rank
	integer,parameter :: nbslpk=16 ! blacs/scalapack block size
	
	
	!**** parameters for CEM
	real*8,parameter :: cd = 299792458d0 ! free-space speed of light
	real*8,parameter :: eps0=1d7/(4d0*pi*cd**2) ! free-space permittivity 
    real*8,parameter :: mu0=pi*4d-7 ! free-space permeability
    real*8,parameter :: gamma=1.781072418d0 ! gamma constant 
    real*8,parameter :: impedence0=sqrt(mu0/eps0) ! free-space wave impedance 
	
	
	!**** solver parameters
	integer,parameter:: DIRECT=1	 ! use factored HODLR as direct solver
	integer,parameter:: NOPRECON=2  ! use compressed HODLR as fast matvec 
	integer,parameter:: HODLRPRECON=3	! use factored HODLR as preconditioner 
	integer,parameter:: LplusMax=3


	
	!**** hierarchical process grid associated with one process node in the process tree (used for parallel recursive LR compression) 
	type grid
		integer :: nsprow,nspcol ! number of process rows/columns as square as possible, it's designed that nspcol>=nsprow
		integer :: hprow,hpcol ! head process in the row/column process tree		
		integer :: ctxt ! blacs context
		integer :: Comm ! MPI communicator for this grid
		integer :: gprow,gpcol ! the group number in the row and column dimension, no longer needed once constructed 
		type(grid),pointer::gdc(:)=>null() ! pointer to its two children
	end type grid
	
	
	
	!**** process groups/nodes in the process tree
	type procgroup
		integer :: nprow,npcol,nproc ! number of processors and 2D grids
		integer :: head,tail ! start and end process in the Comm of proctree
		integer :: ctxt ! blacs context 
		integer :: ctxt1D ! blacs context 1D Row noncyclic (used to distribute from 2D grids to customized noncyclic 1D grid)
		integer :: ctxt_head ! blacs context only involving the head process (used to gather and broadcast)
		integer :: Comm ! MPI communicator for all processes in this node
		type(grid),pointer::gd=>null() ! the hierarchical process grid structure associated with each process group
	end type procgroup
	
	
	
	!**** binary process tree
	type proctree
		integer nlevel ! number of tree levels 
		integer Comm ! MPI communicator for all processes in this tree
		integer nproc ! # of processes in this tree
		integer :: MyID ! MPI Rank in Comm
		type(procgroup),allocatable::pgrp(:) ! tree nodes 
	end type proctree
	
	
	
	!**** communication buffer for all to all communication
	type commquant1D
		integer offset ! offset in my local array
		integer size ! size of the message along first dimension 
		complex(kind=8),allocatable::dat(:,:) ! communication buffer
	end type commquant1D
	
	
	
	!**** cluster of points/indices  	
     type basisgroup
		 integer pgno ! process group number of this cluster
         integer head ! head index
         integer tail ! tail index
         integer level ! level of this cluster group 
         real*8 radius ! geomerical radius of this group 
		 real*8:: boundary(2) ! seperators used to split this group into children group
         real*8,allocatable:: center(:) ! geometrical center of this group
     end type basisgroup

	 
	 
	 !**** input and output vectors for applying a Bplus
     type vectorsblock
         ! integer style
         ! integer head
         ! integer tail
         complex(kind=8),allocatable :: vector(:,:)
     end type vectorsblock
     
	 
	 !**** a vector used to extract one element of a butterfly
     type vectorset
         complex(kind=8),allocatable :: vector(:)
     end type vectorset

 
	 !**** one rank*rank butterfly block
     type butterflymatrix
         complex(kind=8),allocatable :: matrix(:,:) ! entries of the block
		 integer::mdim,ndim	 ! dimensions of the block	 
     end type butterflymatrix
	 
	 
     !**** keep track of skeleton columns
     type butterfly_col_select
         integer,allocatable :: select_columns(:)
     end type butterfly_col_select
     
	 
	 !**** one interior factor
     type butterfly_Kerl
         integer num_col ! # block columns 
         integer num_row ! # block rows 
         type(butterflymatrix),allocatable :: blocks(:,:)
     end type butterfly_Kerl

	 
	 !**** one outter most factor
     type butterfly_UV
         integer num_blk
         type(butterflymatrix),allocatable :: blocks(:)
     end type butterfly_UV
     
	 
	 !**** butterfly or LR structure 
     type matrixblock
		 integer pgno ! process group 
 		 integer pgno_db ! process group when #MPI is doubled 
         integer level ! level in HODLR
         integer col_group ! column group number 
         integer row_group ! row group number 
         integer style ! 1: full block 2: compressed block 
         integer level_butterfly ! butterfly levels 
		 integer rankmax,rankmin ! maximum and minimum butterfly ranks
		 integer dimension_rank ! estimated maximum rank
		 integer M,N ! size of the block
		 integer M_loc,N_loc,M_loc_db,N_loc_db ! local size of the block 
		 integer headm,headn ! header indices in row and column dimension
         integer,pointer:: M_p(:,:)=>null() ! row sizes of all processes sharing this block 
         integer,pointer:: N_p(:,:)=>null() ! column sizes of all processes sharing this block 
         integer,pointer:: M_p_db(:,:)=>null()
         integer,pointer:: N_p_db(:,:)=>null()
		 complex(kind=8),allocatable :: fullmat(:,:) ! full matrix entries
         complex(kind=8),allocatable:: KerInv(:,:)	! a small random rank*rank block used in randomized construction 
		 type(butterfly_UV) :: ButterflyU ! leftmost factor 
         type(butterfly_UV) :: ButterflyV ! rightmost factor 
         type(butterflymatrix),allocatable :: ButterflyMiddle(:,:) ! middle factor 		 
         type(butterfly_Kerl),allocatable :: ButterflyKerl(:) ! interior factors 
         type(butterfly_col_select),allocatable :: ButterflyColSelect(:,:) ! keep track of skeleton columns
     end type matrixblock
	 
	 
	 !**** one layer in a Bplus
	 type onelplus
		integer Nbound ! # of corrected blocks that are further decomposed into deeper layers
		integer rankmax ! maximum butterfly rank on this layer
		type(matrixblock),pointer:: matrices_block(:)=>null() 
		integer,allocatable::boundary_map(:) ! closest subgroup for each subgroup
	 end type onelplus
	 
	 
	 
	 !**** Bplus structure 
	 type blockplus
         integer level ! block level in HODLR
         integer col_group ! column group number 
         integer row_group ! row group number 
		 integer pgno   ! process group number 
         integer Lplus  ! Number of Bplus layers, 2 if option%TwoLayerOnly=1
		 real*8:: boundary(2) ! A analytical seperator defined by one out of three coordinates, boundary(1): direction, boundary(2): value  
		 type(onelplus),pointer:: LL(:)=>null() !  
	 end type blockplus
	 
	 
	 
	 !**** Structure holding operand in Schulz iteration
	 type schulz_operand
		type(matrixblock):: matrices_block ! the original butterfly B in I+B
		real*8::A2norm ! largest singular value of B in I+B
		real*8::scale ! scalar carried on Identities
		real*8,allocatable::diags(:)
		integer order ! order of schulz iteration 
	 end type schulz_operand
	
	
	!**** One level in HODLR
	 type cascadingfactors
		integer level  ! level number 
		integer N_block_forward ! # of forward blocks 
		integer N_block_inverse ! # of inverse blocks 
		integer Bidxs,Bidxe   ! indice range of my local groups		
		type(blockplus),pointer:: BP(:)=>null()  ! forward blocks 
		type(blockplus),pointer:: BP_inverse(:)=>null() ! inverse blocks 
		type(blockplus),pointer:: BP_inverse_schur(:)=>null() ! schur complement blocks 
	 end type cascadingfactors

	 
	!**** HODLR structure  
	type hobf
		integer Maxlevel,N ! HODLR levels and sizes 
		integer ind_lv,ind_bk ! iterator of level and block number in a HODLR
		type(cascadingfactors),allocatable::levels(:) ! 	 
	end type hobf 
	 
	 
	!**** partitioned blocks for reursive computing (I+B)^-1
	 type partitionedblocks
		integer level
		type(matrixblock),pointer:: blocks_A=>null()
		type(matrixblock),pointer:: blocks_B=>null()
		type(matrixblock),pointer:: blocks_C=>null()
		type(matrixblock),pointer:: blocks_D=>null()
	 end type partitionedblocks

	 
	 !**** intermidate vectors for applying a butterfly
     type RandomBlock
         integer level_butterfly
         type(butterfly_Kerl), allocatable :: RandomVectorRR(:)
         type(butterfly_Kerl), allocatable :: RandomVectorLL(:)
     end type RandomBlock
	 
	 
	 !**** HODLR solver options
	 type Hoption
		
		! options for Bplus, Butterfly or LR
		integer::LRlevel  ! The top LRlevel level blocks are butterfly or Bplus		
		integer:: LnoBP ! the bottom LnoBP levels are either Butterfly or LR, but not Bplus  		
		integer:: TwoLayerOnly  ! restrict Bplus as Butterfly + LR
		real*8 touch_para   ! parameters used to determine whether one patch is closer to seperator		
		
		! options for matrix fill
		real*8 tol_SVD      ! matrix fill tolerance
		integer::Nmin_leaf ! leaf sizes of HODLR tree		
		integer xyzsort ! clustering methods given geometrical points: 1: cartesian kd tree 2: spherical kd tree (only for 3D points) 3: 2 means (no recursive)
		character::RecLR_leaf ! bottom level operations in a recursive merge-based LR compression 'Q': RRQR 'S':SVD 'A': ACA
		
		! options for inversion
		real*8 tol_LS       ! tolerance in pseudo inverse
		real*8 tol_Rdetect  ! tolerance to detect numerical ranks
		real*8 tol_rand     ! tolerance for randomized contruction 
		integer::schulzorder ! order (2 or 3) of schultz iteration 
		integer::schulzlevel ! (I+B)^-1 is computed by schultz iteration for butterfly with more than schulzlevel levels 
		
		! options for solve 
		real*8 tol_itersol  ! tolerance for iterative solvers 
		integer N_iter  ! maximum number of iterations for iterative solver
		integer precon  ! DIRECT: use factored HODLR as direct solver, HODLRPRECON: use factored HODLR as preconditioner, NOPRECON: use forward HODLR as fast matvec,  
		
		
		! options for error checking 
		integer::level_check ! check compression quality by picking random entries at level_check (only work for nmpi=1 now)		
		integer::ErrFillFull ! check compression quality by computing all block elements
		integer::ErrSol ! check solution quality by using artificially generated true solution vector

	 end type Hoption
	

	!**** statistics 	
	type Hstat
		real*8 Time_random(5)  ! Intialization, MVP, Reconstruction, Reconstruction of one subblock 
		real*8 Time_Sblock,Time_Inv,Time_SMW,Time_Fill,Time_RedistB,Time_RedistV,Time_Sol
		real*8 Mem_peak,Mem_Sblock,Mem_SMW,Mem_Direct,Mem_int_vec,Mem_For
		real*8 Flop_Fill, Flop_Factor, Flop_Sol, Flop_Tmp
		integer, allocatable:: rankmax_of_level(:),rankmin_of_level(:) ! maximum and minimum ranks observed at each level of HODLR during matrix fill and factorization
	 end type Hstat	

	 
	!**** quantities related to geometries, meshes, unknowns and points 
	type Mesh	
		integer Nunk ! size of the matrix 
		integer idxs,idxe  ! range of local row/column indices after reordering
		real*8,allocatable:: xyz(:,:)   ! coordinates of the points
		integer,allocatable:: info_unk(:,:)  
		! for 2D mesh: 0 point to coordinates of each edge center (unknown x), 1-2 point to coordinates of each edge vertice  
		! for 3D mesh: 0 point to coordinates of each edge center (unknown x), 1-2 point to coordinates of each edge vertice, 3-4 point to two patches that share the same edge, 5-6 point to coordinates of last vertice of the two patches
		! for general: 0 point to coordinates of each unknown x
		
		
		! 2D/3D mesh
		integer maxnode ! # of vertices in a mesh
		integer maxedge ! # of edges in a mesh 
		real*8 maxedgelength,minedgelength ! maximum and minimum edge length for 2D and 3D meshes	
		
		! 3D mesh
		integer integral_points ! #of Gauss quadrature points on a triangular
		integer maxpatch ! # of triangular patches
		integer mesh_normal	 ! flags to flip the unit normal vectors of a triangular patch
		real*8 scaling  ! scaling factor of the coordinates of vertices of a 3D mesh 
		real*8, allocatable :: ng1(:),ng2(:),ng3(:),gauss_w(:) ! Gass quadrature and weights
		real*8,allocatable:: normal_of_patch(:,:) ! normal vector of each triangular patch
		integer,allocatable:: node_of_patch(:,:) ! vertices of each triangular patch
		real*8:: Origins(3) ! orgin of the 3D points, only used for spherical coordinate transform	
		
		! 2D mesh
		integer model2d ! # shape of 2-D curves: (1=strip; 2=corner reflector; 3=two opposite strips; 4=CR with RRS; 5=cylinder; 6=Rectangle Cavity); 7=half cylinder; 8=corrugated half cylinder; 9=corrugated corner reflector; 10=open polygon; 11=taller open polygon
		real*8 Delta_ll	! edge length of each element in 2D curves	
		integer::Ncorner ! # of corners in 2D curves to facilate partitioning
		real*8,allocatable::corner_points(:,:) ! coordinates of corner points 
		real*8::corner_radius ! radius of the corner points within which no partitioning is performed 
		
	end type Mesh
	 
	 
	!**** quantities related to specific matrix kernels  
	type kernelquant
		real*8 wavenum    ! CEM: wave number  
		real*8 wavelength  ! CEM: wave length
		real*8 omiga       ! CEM: angular frequency
		real*8 rank_approximate_para1, rank_approximate_para2, rank_approximate_para3 ! CEM: rank estimation parameter  
		integer RCS_static  ! CEM: 1: monostatic or 2: bistatic RCS
		integer RCS_Nsample ! CEM: number of RCS samples
		real*8:: CFIE_alpha ! CEM: combination parameter in CFIE
		
		real*8 sigma, lambda ! Kernel Regression: RBF parameters 		
		character(LEN=500)::trainfile_p,trainfile_tree,trainfile_l,testfile_p,testfile_l !Kernel Regression: file pointers to train data, preordered tree, train labels, test data, and test labels
		
		complex(kind=8), allocatable :: matZ_glo(:,:) ! Full Matrix: the full matrix to sample its entries
		
		type(c_ptr),pointer :: QuantZmn ! Kernels Defined in C: c_pointer to the user-supplied object for computing one element of Z 
		type(c_funptr),pointer :: FuncZmn ! Kernels Defined in C: c_function_pointer to the user-supplied function for computing one element of Z
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
	

	
	!*** need to further work on the following: 
	 real*8::time_tmp
     type(basisgroup),allocatable:: basis_group(:)
	 integer vecCNT
	 integer,allocatable:: new2old(:)
	 integer,allocatable:: basis_group_pre(:,:)
	 CHARACTER (LEN=1000) DATA_DIR	

end module

