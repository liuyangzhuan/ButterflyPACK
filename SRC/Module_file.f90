module MODULE_FILE
	
	integer,parameter :: MaxDim = 100
	
     type basisgroup
         integer head
         integer tail
         integer level
         real*8 radius
         real*8 center(MaxDim)
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
	 
	 
	
	 type cascadingfactors
		integer level
		integer N_block_forward
		integer N_block_inverse
		type(matrixblock),pointer:: matrices_block(:)	
		type(matrixblock),pointer:: matrices_block_inverse(:)		
		type(matrixblock),pointer:: matrices_block_inverse_schur(:)		
		type(blockplus),pointer:: BP(:)
		type(blockplus),pointer:: BP_inverse(:)
		type(blockplus),pointer:: BP_inverse_schur(:)
	 end type cascadingfactors

	 
	type hobf
		integer Maxlevel_for_blocks
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
     
     ! type matrixblock
         ! integer dimension_m
         ! integer dimension_n
         ! integer dimension_rank
		 ! integer rankmax,rankmin
         ! integer level_butterfly
		 ! complex(kind=8),allocatable:: KerInv(:,:)									
         ! type(butterflymatrix),allocatable :: ButterflyU(:)
         ! type(butterflymatrix),allocatable :: ButterflyV(:)
         ! type(butterflymatrix),allocatable :: ButterflyMiddle(:,:)			 
         ! type(butterfly_Kerl),allocatable :: ButterflyKerl(:)
     ! end type matrixblock

     real*8, allocatable :: ng1(:),ng2(:),ng3(:),gauss_w(:)
     integer,allocatable:: node_of_patch(:,:),node_patch_of_edge(:,:)
     real*8,allocatable:: normal_of_patch(:,:)
     real*8,allocatable:: xyz(:,:)
     complex(kind=8),allocatable:: current2com(:,:),current(:),current_pre(:),vector_Vinc(:,:),voltage(:)
	 
	 integer error_cnt
	 real*8 error_avr_glo
     real*8 wa,wb,wc
     real*8 pi,impedence
     real*8 wavenum,wavelength
     real*8 mu0,eps0,omiga
     real*8 alpha, Delta_ll, gamma
     real*8 Scale
     complex(kind=8) junit
     integer maxpatch,maxnode,maxedge, geo_model, RCS_sample
	 real*8 maxedgelength,minedgelength
     integer Maxlevel, Maxlevel_for_blocks, Refined_level, mesh_normal
     integer RCS_static, Recompress_method, Forward_recompress
     integer integral_points, SAI_level, Diag_level
     integer maxvector_near, maxvector_SAI, Maxvector_Diagblock
     integer preconditioner, Forward_compress_elem,adaptive_flag,reducelevel_flag,verboselevel,directschur  
     complex(kind=8) Test_value
     integer Test_blocks, Test_i, Test_j, Preset_level_butterfly, Optimizing_forward
	 real*8 Bigvalue, ACA_tolerance_forward, SVD_tolerance_forward,SVD_tolerance_factor, ACA_tolerance_split, SVD_tolerance_split,relax_lamda,Rank_detection_factor
	 real*8 SVD_tolerance_add, ACA_tolerance_add, Discret, levelpara_control, ACA_tolerance_aggregate, SVD_tolerance_aggregate
	 real*8 touch_para, rank_approximate_para1, rank_approximate_para2, rank_approximate_para3
     integer Maxgroup,Maxblock, Split_method, Forward_method, Add_method, Aggregate_method, Add_method_of_base_level
     integer Nmin_leaf, Static, Fast_inverse, rank_control_forward, rank_control_add, rank_control_split, rank_control_aggregate
     real*8 Time_Direct_LLD, Time_Add_Multiply_Baselevel, Time_Multiply_Butterfly, Time_Solve_XLM
     real*8 Time_Split_Blocks, Time_Aggregate_Blocks, Time_Aggregate_blocks_OPT, Time_Add_Butterfly,Time_Rightmultiply_inverse_randomized,Time_Inversion_diagnal_randomized
	 real*8 Time_Init_forward,Time_Vector_forward,Time_Reconstruct_forward,Time_Oneblock_forward,Time_Init_inverse,Time_Vector_inverse,Time_Reconstruct_inverse,Time_Oneblock_inverse,Time_InvertBlock,Time_Init_HODLR_MVP
     real*8::mem_target,mem_recons,time_indexarray,time_leastsquare,time_kernelupdate,time_buttermul,time_buttermulinv
	 real*8::time_getvec,time_halfbuttermul,time_resolve,time_memcopy,time_gemm,time_tmp,time_gemm1
	 
     type(RandomBlock),pointer :: Random_Block(:)
     type(basisgroup),allocatable:: basis_group(:)
     type(matrixblock),allocatable:: matrices_block(:,:)
	 type(matrixblock),allocatable:: agent_block(:)
	 type(blockplus),allocatable:: agent_bplus(:)
     type(matrixblock),pointer :: butterfly_block_randomized(:)
     type(blockplus),pointer :: Bplus_randomized(:)
	 ! type(cascadingfactors),allocatable::cascading_factors(:),cascading_factors_copy(:)
	 type(hobf)::ho_bf,ho_bf_copy
	 type(partitionedblocks),allocatable::partitioned_blocks(:)
     type(vectorsblock),allocatable:: vectors_block(:)
	 type(vectorsblock),pointer:: RandomVectors_InOutput(:)
	 ! type(vectorsblock),pointer:: RandomVectors_InOutput_tmp(:)
     type(butterfly_Kerl),allocatable:: ButterflyKerl(:), Kerlmat_Global(:), Kerlmat_Global_Shrink(:)
     type(butterfly_Kerl),allocatable:: Kerlmat_Global_Shrink1(:), Kerlmat_Global_Shrink2(:), Kerlmat_Global_Shrink3(:), Kerlmat_Global_Shrink4(:)
     type(butterflymatrix),allocatable:: ButterflyU_Global(:), ButterflyV_Global(:)
     type(matrixblock),allocatable:: matricesblocktemp(:)
     integer, allocatable:: leafs_of_level(:), leafindex_of_level(:),rankmax_of_level(:),rankmin_of_level(:)
     integer, allocatable:: index_of_Diagblock(:,:), index_of_SAIblock(:,:)
     complex(kind=8),allocatable:: matrixtemp_U(:,:),matrixtemp_V(:,:),fullmatrix(:,:)
     complex(kind=8),allocatable:: matrixtemp_UU(:,:),matrixtemp_VV(:,:), MatrixSubselection(:,:)
     complex(kind=8),allocatable:: UU_1(:,:),UU_2(:,:),UU_3(:,:),UU_4(:,:),VV_1(:,:),VV_2(:,:),VV_3(:,:),VV_4(:,:)
     real*8,allocatable:: Single_1(:),Single_2(:),Single_3(:),Single_4(:)
     integer(kind=8),allocatable:: CPU_clocks(:)
     integer, allocatable :: rankmax_for_butterfly(:),rankmin_for_butterfly(:)
	 real*8:: LS_tolerance, iter_tolerance,tfqmr_tolerance,tfqmr_tolerance_solving,up_tolerance
	 integer:: schurinv,SolvingMethod,level_tmp,rank_tmp,ForwardSymmetricFlag,SblockSymmetricFlag
	 real*8:: Memory_int_vec,Memory_tfqmr_vec
	 complex(kind=8),allocatable:: RandVectInR(:,:),RandVectOutR(:,:),RandVectInL(:,:),RandVectOutL(:,:)
	 complex(kind=8), allocatable :: Basis_split(:,:)
	 integer::Ncorner, ranktmp_glo
	 real*8,allocatable::corner_points(:,:),edge_cen(:,:)
	 real*8::corner_radius
	 integer xyzsort
	 integer,parameter:: LplusMax=5
	 integer:: LnoBP 
	 integer:: TwoLayerOnly 
	 real*8:: CFIE_alpha
		 
	 real*8:: Origins(3)
	 complex(kind=8), allocatable :: matU_glo(:,:), matV_glo(:,:),matSub_glo(:,:),matZ_glo(:,:)
	 integer explicitflag, fullmatflag, Nin1,Nout1,Nin2,Nout2, vecCNT
	 integer,allocatable:: new2old(:)
	 integer,allocatable:: basis_group_pre(:,:)
	 
	 real(kind=8),parameter:: SafeUnderflow=1D-150
	 real(kind=8),parameter:: SafeAbsoulte=1D-14
	 
	 real*8 sigma, lambda ! parameters in Kernel matrices
	 
	 integer,parameter:: EMCURV=1
	 integer,parameter:: EMSURF=2
	 integer,parameter:: RBF=3
	 integer,parameter:: FULL=4
	 integer Kernel
	 
	 
	 integer,parameter:: DIRECT=1
	 integer,parameter:: NOPRECON=2
	 integer,parameter:: SAIPRECON=3
	 
	 
	 CHARACTER (LEN=1000) DATA_DIR	 
	 
end module

module Super_Block
    integer,allocatable :: index_of_superblock(:,:)
    logical,allocatable :: nearorfar_of_superblock(:,:)
    integer Maxsuperblock
end module
