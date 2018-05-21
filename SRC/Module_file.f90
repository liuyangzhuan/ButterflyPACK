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
	
	
	type procgroup
		!****Define MPI variables
		integer :: nprow,npcol,nproc ! number of processors
		integer :: head,tail
		integer :: ctxt
	end type procgroup
	
	type proctree
		integer nlevel
		integer Comm
		integer nproc
		integer :: MyID
		type(procgroup),allocatable::pgrp(:)
	end type proctree
	

     type basisgroup
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
		 integer::nulldim		 
     end type butterflymatrix
     
     type butterfly_col_select
         integer,allocatable :: select_columns(:)
     end type butterfly_col_select
     
     type butterfly_Kerl
         integer num_col
         integer num_row
         type(butterflymatrix),allocatable :: blocks(:,:)
     end type butterfly_Kerl
     
     type butterfly_Kerl_Global
         integer maxlevel
         integer num_col
         integer num_row
         type(butterflymatrix),allocatable :: blocks(:,:)
     end type butterfly_Kerl_Global
     
     
     type butterfly_multivectors
         integer num
         type(vectorset), allocatable :: blocks(:)
     end type butterfly_multivectors

     type matrixblock
         integer level
         integer col_group
         integer row_group
         integer nested_num
         integer style
         integer data_type
         integer level_butterfly
		 integer rankmax,rankmin
		 integer dimension_rank
         complex(kind=8),allocatable :: fullmat(:,:)
         complex(kind=8),allocatable:: KerInv(:,:)	
		 type(butterflymatrix),allocatable :: ButterflyU(:)
         type(butterflymatrix),allocatable :: ButterflyV(:)
         type(butterflymatrix),allocatable :: ButterflyMiddle(:,:)		 
         type(butterfly_Kerl),allocatable :: ButterflyKerl(:)
         type(butterfly_col_select),allocatable :: ButterflyColSelect(:,:)
     end type matrixblock
	 
	 
	 type onelplus
		integer Nbound
		integer rankmax
		type(matrixblock),pointer:: matrices_block(:)
		integer,allocatable::boundary_map(:)
	 end type onelplus
	 
	 
	 type blockplus
         integer level
         integer col_group
         integer row_group
         integer Lplus
		 real*8:: boundary(2)
		 type(onelplus),pointer:: LL(:) 
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
		! type(matrixblock),pointer:: matrices_block(:)	
		! type(matrixblock),pointer:: matrices_block_inverse(:)		
		! type(matrixblock),pointer:: matrices_block_inverse_schur(:)		
		type(blockplus),pointer:: BP(:)
		type(blockplus),pointer:: BP_inverse(:)
		type(blockplus),pointer:: BP_inverse_schur(:)
	 end type cascadingfactors

	 
	type hobf
		integer Maxlevel,N
		integer ind_lv,ind_bk
		type(cascadingfactors),allocatable::levels(:)	 
	end type hobf 
	 

	 type partitionedblocks
		integer level
		type(matrixblock),pointer:: blocks_A,blocks_B,blocks_C,blocks_D
	 end type partitionedblocks

	 
         
     type blockindex
         integer num_row
         integer num_col
         integer demension_row
         integer demension_col
         integer, allocatable :: IndicesForward(:,:,:)
         integer, allocatable :: IndicesBackward(:,:,:)
     end type blockindex
     
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
		real*8 Time_Sblock,Time_SMW
		real*8 Mem_peak,Mem_Sblock,Mem_SMW,Mem_Direct,Mem_int_vec
		integer, allocatable:: rankmax_of_level(:),rankmin_of_level(:)
	 end type Hstat	

	type Mesh
		integer integral_points
		integer maxpatch
		integer maxnode
		integer maxedge
		integer model2d
		integer Nunk
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

		subroutine BF_MVP_blk(operand,block_o,trans,M,N,num_vect_sub,Vin,Vout,a,b,operand1)
			import :: matrixblock
			implicit none
			class(*)::operand	
			class(*),optional::operand1	
			type(matrixblock)::block_o
			character trans
			integer M, N, num_vect_sub
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

