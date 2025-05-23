!> @file
!> @brief This example generates a random LR product, or reads a full matrix from disk, and compress it using entry-valuation-based APIs
!> @details Note that instead of the use of precision dependent subroutine/module/type names "z_", one can also use the following \n
!> #define DAT 0 \n
!> #include "zButterflyPACK_config.fi" \n
!> which will macro replace precision-independent subroutine/module/type names "X" with "z_X" defined in SRC_DOUBLECOMLEX with double-complex precision

module APPLICATION_MODULE_FULL_SIMPLE
use z_BPACK_DEFS
implicit none

	!**** define your application-related variables here
	type quant_app
		complex(kind=8), allocatable :: matU_glo(:,:),matV_glo(:,:) ! Full Matrix: the random LR matrix to sample its entries
		complex(kind=8), allocatable :: matZ_glo(:,:) ! Full Matrix: Full matrix read from files
		real(kind=8), allocatable :: locations(:,:),locations_m(:,:),locations_n(:,:) ! geometrical points
		integer,allocatable:: permutation_m(:),permutation_n(:)
		integer:: rank
		integer:: Nunk,Nunk_m,Nunk_n
		integer:: Ndim
		real(kind=8):: lambda
		CHARACTER (LEN=1000) DATA_DIR
		CHARACTER (LEN=1000) GEO_DIR,GEO_DIR_m,GEO_DIR_n
		integer:: tst=1
	end type quant_app

contains

	!**** user-defined subroutine to sample Z_mn as two LR products
	subroutine Zelem_LR(m,n,value_e,quant)
		use z_BPACK_DEFS
		implicit none

		class(*),pointer :: quant
		integer, INTENT(IN):: m,n
		complex(kind=8)::value_e
		integer ii

		integer dimn

		select TYPE(quant)
		type is (quant_app)
			value_e = 0
			do ii=1,quant%rank
				value_e = value_e + quant%matU_glo(m,ii)*quant%matV_glo(ii,n)
			enddo
			if(m==n)then
				value_e = value_e + quant%lambda
			endif

			! value_e = quant%matZ_glo(m,n)
		class default
			write(*,*)"unexpected type"
			stop
		end select
	end subroutine Zelem_LR


	!**** user-defined subroutine to sample Z_mn as two LR products (note that this is for the BF interface not for the BPACK interface)
	subroutine ZBelem_LR(m,n,value_e,quant)
		use z_BPACK_DEFS
		implicit none

		class(*),pointer :: quant
		integer, INTENT(IN):: m,n
		complex(kind=8)::value_e
		integer ii,m1,n1

		if(m>0)then
			m1=m
			n1=-n
		else
			m1=n
			n1=-m
		endif

		!!! m,n still need to convert to the original order, using new2old of mshr and mshc
		select TYPE(quant)
		type is (quant_app)
			m1=quant%permutation_m(m1)
			n1=quant%permutation_n(n1)
		class default
			write(*,*)"unexpected type"
			stop
		end select


		call Zelem_LR(m1,n1,value_e,quant)
	end subroutine ZBelem_LR


	!**** user-defined subroutine to sample Z_mn as full matrix
	subroutine Zelem_FULL(m,n,value,quant)
		use z_BPACK_DEFS
		implicit none

		class(*),pointer :: quant
		integer, INTENT(IN):: m,n
		complex(kind=8)::value
		integer ii

		integer dimn

		select TYPE(quant)
		type is (quant_app)
			value = quant%matZ_glo(m,n)
		class default
			write(*,*)"unexpected type"
			stop
		end select
	end subroutine Zelem_FULL


	!**** user-defined subroutine to sample Z_mn as full matrix (note that this is for the BF interface not for the BPACK interface)
	subroutine ZBelem_FULL(m,n,value_e,quant)
		use z_BPACK_DEFS
		implicit none

		class(*),pointer :: quant
		integer, INTENT(IN):: m,n
		complex(kind=8)::value_e
		integer ii,m1,n1

		if(m>0)then
			m1=m
			n1=-n
		else
			m1=n
			n1=-m
		endif

		!!! m,n still need to convert to the original order, using new2old of mshr and mshc
		select TYPE(quant)
		type is (quant_app)
			m1=quant%permutation_m(m1)
			n1=quant%permutation_n(n1)
		class default
			write(*,*)"unexpected type"
			stop
		end select

		call Zelem_FULL(m1,n1,value_e,quant)
	end subroutine ZBelem_FULL


end module APPLICATION_MODULE_FULL_SIMPLE


PROGRAM ButterflyPACK_TEMPLATE
    use z_BPACK_DEFS
    use APPLICATION_MODULE_FULL_SIMPLE
	use z_BPACK_Solve_Mul

	use z_BPACK_structure
	use z_BPACK_factor
	use z_BPACK_constr
#ifdef HAVE_OPENMP
	use omp_lib
#endif
	use z_MISC_Utilities
	use z_BPACK_constr
	use z_BPACK_utilities
    implicit none

    integer rank,ii
	real(kind=8),allocatable:: datain(:)

	integer :: ierr
	type(z_Hoption),target::option
	type(z_Hstat),target::stats
	type(z_mesh),target::msh,mshr,mshc
	type(z_kernelquant),target::ker
	type(quant_app),target::quant
	type(z_Bmatrix),target::bmat
	integer,allocatable:: groupmembers(:)
	integer nmpi
	integer level,Maxlevel
	type(z_proctree),target::ptree
	integer,allocatable::Permutation(:)
	integer Nunk_loc,Nunk_m_loc, Nunk_n_loc
	integer,allocatable::tree(:),tree_m(:),tree_n(:)
	complex(kind=8),allocatable::rhs_glo(:,:),rhs_loc(:,:),x_glo(:,:),x_loc(:,:)
	integer nrhs
	type(z_matrixblock) ::blocks
	type(z_blockplus) ::BP
	character(len=1024)  :: strings,strings1
	integer flag,nargs

	!**** nmpi and groupmembers should be provided by the user
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo

	!**** create the process tree
	call z_CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	deallocate(groupmembers)
	!**** initialize stats and option
	call z_InitStat(stats)
	call z_SetDefaultOptions(option)


	!**** set solver parameters
	option%ErrSol=1  ! whether or not checking the factorization accuracy
	option%format=  HODLR! HMAT!   ! the hierarhical format
	option%near_para=2.01d0        ! admissibiltiy condition, not referenced if option%format=  HODLR
	option%verbosity=1             ! verbosity level
	option%LRlevel=0             ! 0: low-rank compression 100: butterfly compression


	quant%tst = 1 ! 1: use a LR product as the kernel, 2: read the full matrix from a file as the kernel

	! quant%DATA_DIR = '../EXAMPLE/FULLMAT_DATA/K05N4096.csv' ! file storing the full matrix
	! quant%Nunk = 4096  ! matrix size
	! quant%DATA_DIR = '../EXAMPLE/FULLMAT_DATA/A_alpha_N64.csv' ! file storing the full matrix
	! quant%Nunk = 63  ! matrix size
	! quant%DATA_DIR = '../EXAMPLE/FULLMAT_DATA/A_alpha_N128.csv' ! file storing the full matrix
	! quant%Nunk = 127  ! matrix size
	! quant%DATA_DIR = '../EXAMPLE/FULLMAT_DATA/A_alpha_N256.csv' ! file storing the full matrix
	! quant%Nunk = 255  ! matrix size
	quant%DATA_DIR = '../EXAMPLE/FULLMAT_DATA/A_alpha_N512.csv' ! file storing the full matrix
	quant%Nunk = 511  ! matrix size

	quant%rank = 6   ! rank of the LR product kernel
	option%nogeo = 1  ! 1. no geometry info available. 2. geometry info available


	if(option%nogeo==1)then
		! no geometry points available
		option%xyzsort=NATURAL ! no reordering will be perfomed
	else
		! geometry points available
		option%xyzsort=TM ! no reordering will be perfomed
		option%knn=20   ! neareat neighbour points per geometry point, which helps improving the compression accuracy
		quant%GEO_DIR = '../EXAMPLE/FULLMAT_DATA/Geometry_3D.csv' ! file storing the geometry
		quant%Ndim = 3 ! dimension of the geometry information, not referenced if option%nogeo=1
	endif



	nargs = iargc()
	ii=1
	do while(ii<=nargs)
		call getarg(ii,strings)
		if(trim(strings)=='-quant')then ! user-defined quantity parameters
			flag=1
			do while(flag==1)
				ii=ii+1
				if(ii<=nargs)then
					call getarg(ii,strings)
					if(strings(1:2)=='--')then
						ii=ii+1
						call getarg(ii,strings1)
						if(trim(strings)=='--tst')then
							read(strings1,*)quant%tst
						else
							if(ptree%MyID==Main_ID)write(*,*)'ignoring unknown quant: ', trim(strings)
						endif
					else
						flag=0
					endif
				else
					flag=0
				endif
			enddo
		else if(trim(strings)=='-option')then ! options of ButterflyPACK
			call z_ReadOption(option,ptree,ii)
		else
			if(ptree%MyID==Main_ID)write(*,*)'ignoring unknown argument: ',trim(strings)
			ii=ii+1
		endif
	enddo



	call z_PrintOptions(option,ptree)

!******************************************************************************!
! generate a LR matrix as two matrix product and do a H matrix compression and factorization
	if(quant%tst==1)then
		allocate(tree(1))
		tree=quant%Nunk
		!**** Get matrix size and rank and create the matrix
		quant%lambda = 1d5
		allocate(quant%matU_glo(quant%Nunk,quant%rank))
		call z_RandomMat(quant%Nunk,quant%rank,quant%rank,quant%matU_glo,0)
		call MPI_Bcast(quant%matU_glo,quant%Nunk*quant%rank,MPI_DOUBLE_COMPLEX,Main_ID,ptree%Comm,ierr)

		allocate(quant%matV_glo(quant%rank,quant%Nunk))
		call z_RandomMat(quant%rank,quant%Nunk,quant%rank,quant%matV_glo,0)
		call MPI_Bcast(quant%matV_glo,quant%Nunk*quant%rank,MPI_DOUBLE_COMPLEX,Main_ID,ptree%Comm,ierr)
	   !***********************************************************************
	   if(ptree%MyID==Main_ID)then
	   write (*,*) ''
	   write (*,*) 'Random LR Kernel computing'
	   write (*,*) 'Matrix size:', quant%Nunk
       write (*,*) ''
	   endif
	   !***********************************************************************

		!**** register the user-defined function and type in ker
		ker%QuantApp => quant
		ker%FuncZmn => Zelem_LR

		!**** initialization of the construction phase
		allocate(Permutation(quant%Nunk))
		call z_BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree)
		call MPI_Bcast(Permutation,quant%Nunk,MPI_integer,0,ptree%comm,ierr)

		!**** computation of the construction phase
		call z_BPACK_construction_Element(bmat,option,stats,msh,ker,ptree)

		!**** factorization phase
		call z_BPACK_Factorization(bmat,option,stats,ptree,msh)

		!**** generate testing RHSs stored globally on each rank
		nrhs=2
		allocate(rhs_glo(quant%Nunk,nrhs))
		allocate(x_glo(quant%Nunk,nrhs))
		rhs_glo=1d0
		call MPI_Bcast(rhs_glo,quant%Nunk*nrhs,MPI_DOUBLE_COMPLEX,0,ptree%comm,ierr)

		!**** convert global RHSs to local RHSs
		allocate(rhs_loc(msh%idxe-msh%idxs+1,nrhs))
		allocate(x_loc(msh%idxe-msh%idxs+1,nrhs))
		x_loc=0
		do ii=1,msh%idxe-msh%idxs+1
			rhs_loc(ii,:) = rhs_glo(Permutation(ii+msh%idxs-1),:)
		end do

		!**** call the solve routine
		call z_BPACK_Solution(bmat,x_loc,rhs_loc,msh%idxe-msh%idxs+1,1,option,ptree,stats)

		!**** convert local solutions to global solutions
		x_glo=0
		do ii=1,msh%idxe-msh%idxs+1
			x_glo(Permutation(ii+msh%idxs-1),:) = x_loc(ii,:)
		end do
		call MPI_ALLREDUCE(MPI_IN_PLACE,x_glo,quant%Nunk*nrhs,MPI_DOUBLE_COMPLEX,MPI_SUM, ptree%Comm, ierr)

	endif

!******************************************************************************!
! generate a full matrix stored in a files and do a H matrix compression and factorization
	if(quant%tst==2)then
		allocate(tree(1))
		tree=quant%Nunk

		!**** Get matrix size and rank and create the matrix
		allocate(quant%matZ_glo(quant%Nunk,quant%Nunk))

		!***** assuming reading one row every time, reading a real matrix stored in file as a complex matrix
		allocate(datain(quant%Nunk))
		if(ptree%MyID==Main_ID)then
		open(10, file=quant%DATA_DIR)
		do ii=1,quant%Nunk
			read(10,*) datain
			quant%matZ_glo(ii,:)=datain(:)
		enddo
		close(10)
		endif
		deallocate(datain)

		call MPI_Bcast(quant%matZ_glo,quant%Nunk*quant%Nunk,MPI_DOUBLE_COMPLEX,Main_ID,ptree%Comm,ierr)

		if(option%nogeo==0)then
			!***** assuming reading one dimension each time, the geometry is used to partition the matrix
			allocate(quant%locations(quant%Ndim,quant%Nunk))
			quant%locations=0
			allocate(datain(quant%Nunk))
			open(10, file=quant%GEO_DIR)
			do ii=1,quant%Ndim
				read(10,*) datain(:)
				quant%locations(ii,:)=datain(:)
			enddo
			close(10)
			deallocate(datain)
		endif

	   !***********************************************************************
		if(ptree%MyID==Main_ID)then
			write (*,*) ''
			write (*,*) 'FullMat computing'
			write (*,*) 'Matrix size:', quant%Nunk
			write (*,*) ''
		endif
		!***********************************************************************

		!**** register the user-defined function and type in ker
		ker%QuantApp => quant
		ker%FuncZmn => Zelem_FULL

		!**** initialization of the construction phase
		allocate(Permutation(quant%Nunk))
		call z_BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree,Coordinates=quant%locations,tree=tree)
		call MPI_Bcast(Permutation,quant%Nunk,MPI_integer,0,ptree%comm,ierr)
		deallocate(tree)

		!**** computation of the construction phase
		call z_BPACK_construction_Element(bmat,option,stats,msh,ker,ptree)

		deallocate(quant%matZ_glo)
		!**** factorization phase
		call z_BPACK_Factorization(bmat,option,stats,ptree,msh)

		!**** generate testing RHSs stored globally on each rank
		nrhs=2
		allocate(rhs_glo(quant%Nunk,nrhs))
		allocate(x_glo(quant%Nunk,nrhs))
		rhs_glo=1d0
		call MPI_Bcast(rhs_glo,quant%Nunk*nrhs,MPI_DOUBLE_COMPLEX,0,ptree%comm,ierr)

		!**** convert global RHSs to local RHSs
		allocate(rhs_loc(msh%idxe-msh%idxs+1,nrhs))
		allocate(x_loc(msh%idxe-msh%idxs+1,nrhs))
		x_loc=0
		do ii=1,msh%idxe-msh%idxs+1
			rhs_loc(ii,:) = rhs_glo(Permutation(ii+msh%idxs-1),:)
		end do

		!**** call the solve routine
		call z_BPACK_Solution(bmat,x_loc,rhs_loc,msh%idxe-msh%idxs+1,1,option,ptree,stats)

		!**** convert local solutions to global solutions
		x_glo=0
		do ii=1,msh%idxe-msh%idxs+1
			x_glo(Permutation(ii+msh%idxs-1),:) = x_loc(ii,:)
		end do
		call MPI_ALLREDUCE(MPI_IN_PLACE,x_glo,quant%Nunk*nrhs,MPI_DOUBLE_COMPLEX,MPI_SUM, ptree%Comm, ierr)
	endif


!******************************************************************************!
! generate a LR matrix as two matrix product and do a BF compression
	if(quant%tst==3)then
		quant%Nunk_m=1000
		quant%Nunk_n=800

		!**** Get matrix size and rank and create the matrix
		quant%lambda = 0
		allocate(quant%matU_glo(quant%Nunk_m,quant%rank))
		call z_RandomMat(quant%Nunk_m,quant%rank,quant%rank,quant%matU_glo,0)
		call MPI_Bcast(quant%matU_glo,quant%Nunk_m*quant%rank,MPI_DOUBLE_COMPLEX,Main_ID,ptree%Comm,ierr)

		allocate(quant%matV_glo(quant%rank,quant%Nunk_n))
		call z_RandomMat(quant%rank,quant%Nunk_n,quant%rank,quant%matV_glo,0)
		call MPI_Bcast(quant%matV_glo,quant%Nunk_n*quant%rank,MPI_DOUBLE_COMPLEX,Main_ID,ptree%Comm,ierr)
	   !***********************************************************************
	   if(ptree%MyID==Main_ID)then
	   write (*,*) ''
	   write (*,*) 'Random LR Kernel computing'
	   write (*,*) 'Matrix size:', quant%Nunk_m, quant%Nunk_n
       write (*,*) ''
	   endif
	   !***********************************************************************

		!**** register the user-defined function and type in ker
		ker%QuantApp => quant
		ker%FuncZmn => ZBelem_LR

		allocate(quant%Permutation_m(quant%Nunk_m))
		allocate(quant%Permutation_n(quant%Nunk_n))
		call z_BP_Construct_Init(quant%Nunk_m, quant%Nunk_n, Nunk_m_loc, Nunk_n_loc, quant%Permutation_m, quant%Permutation_n, BP, option, stats, msh, ker, ptree)
		call MPI_Bcast(quant%Permutation_m,quant%Nunk_m,MPI_integer,0,ptree%comm,ierr)
		call MPI_Bcast(quant%Permutation_n,quant%Nunk_n,MPI_integer,0,ptree%comm,ierr)

		call z_BP_Construct_Element_Compute(BP, option, stats, msh, ker, ptree)
	endif

!******************************************************************************!
! Read a full non-square matrix and do a BF compression
	if(quant%tst==4)then
		quant%Nunk_m=26040
		quant%Nunk_n=15930
		quant%Ndim=2
		quant%DATA_DIR = '../EXAMPLE/FULLMAT_DATA/Mat_redatuming_Freq200.csv' ! file storing the full matrix
		quant%GEO_DIR_m = '../EXAMPLE/FULLMAT_DATA/Geometry_redatuming_src.csv' ! file storing the geometry
		quant%GEO_DIR_n = '../EXAMPLE/FULLMAT_DATA/Geometry_redatuming_rec.csv' ! file storing the geometry


		!**** Get matrix size and rank and create the matrix
		allocate(quant%matZ_glo(quant%Nunk_m,quant%Nunk_n))
		quant%matZ_glo=1d0

		!***** assuming reading one row every time, the first half of a row is the real part, the second half is the imaginary part
		allocate(datain(quant%Nunk_n*2))
		if(ptree%MyID==Main_ID)then
		open(10, file=quant%DATA_DIR)
		do ii=1,quant%Nunk_m
			read(10,*) datain
			quant%matZ_glo(ii,:)=datain(1:quant%Nunk_n) + BPACK_junit*datain(1+quant%Nunk_n:quant%Nunk_n*2)
		enddo
		close(10)
		endif
		deallocate(datain)

		call MPI_Bcast(quant%matZ_glo,quant%Nunk_m*quant%Nunk_n,MPI_DOUBLE_COMPLEX,Main_ID,ptree%Comm,ierr)

		if(option%nogeo==0)then
			!***** assuming reading one dimension each time, the geometry is used to partition the matrix
			allocate(quant%locations_m(quant%Ndim,quant%Nunk_m))
			quant%locations_m=0
			allocate(datain(quant%Nunk_m))
			open(10, file=quant%GEO_DIR_m)
			do ii=1,quant%Ndim
				read(10,*) datain(:)
				quant%locations_m(ii,:)=datain(:)
			enddo
			close(10)
			deallocate(datain)

			allocate(quant%locations_n(quant%Ndim,quant%Nunk_n))
			quant%locations_n=0
			allocate(datain(quant%Nunk_n))
			open(10, file=quant%GEO_DIR_n)
			do ii=1,quant%Ndim
				read(10,*) datain(:)
				quant%locations_n(ii,:)=datain(:)
			enddo
			close(10)
			deallocate(datain)
		endif


	   if(ptree%MyID==Main_ID)then
	   write (*,*) ''
	   write (*,*) 'FullMat Kernel computing'
	   write (*,*) 'Matrix size:', quant%Nunk_m, quant%Nunk_n
       write (*,*) ''
	   endif
	   !***********************************************************************

		!**** register the user-defined function and type in ker
		ker%QuantApp => quant
		ker%FuncZmn => ZBelem_FULL

		allocate(quant%Permutation_m(quant%Nunk_m))
		allocate(quant%Permutation_n(quant%Nunk_n))
		call z_BP_Construct_Init(quant%Nunk_m, quant%Nunk_n, Nunk_m_loc, Nunk_n_loc, quant%Permutation_m, quant%Permutation_n, BP, option, stats, msh, ker, ptree, Coordinates_m=quant%locations_m,Coordinates_n=quant%locations_n)
		call MPI_Bcast(quant%Permutation_m,quant%Nunk_m,MPI_integer,0,ptree%comm,ierr)
		call MPI_Bcast(quant%Permutation_n,quant%Nunk_n,MPI_integer,0,ptree%comm,ierr)

		call z_BP_Construct_Element_Compute(BP, option, stats, msh, ker, ptree)

		nrhs=1
		allocate(rhs_loc(Nunk_m_loc,nrhs))
		rhs_loc=0
		allocate(x_loc(Nunk_n_loc,nrhs))
		x_loc=1 ! a simple test vector
		call z_BP_Mult(BP, 'N', x_loc, rhs_loc, Nunk_n_loc, Nunk_m_loc, nrhs, ptree, stats)

	endif

!******************************************************************************!


	!**** print statistics
	call z_PrintStat(stats,ptree)

	if(allocated(quant%matU_glo))deallocate(quant%matU_glo)
	if(allocated(quant%matV_glo))deallocate(quant%matV_glo)
	if(allocated(quant%matZ_glo))deallocate(quant%matZ_glo)

	call z_delete_proctree(ptree)
	call z_delete_Hstat(stats)
	call z_delete_mesh(msh)
	call z_delete_kernelquant(ker)
	call z_BPACK_delete(bmat)

	if(allocated(Permutation))deallocate(Permutation)

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	call z_blacs_exit_wrp(1)
	call MPI_Finalize(ierr)

end PROGRAM ButterflyPACK_TEMPLATE



