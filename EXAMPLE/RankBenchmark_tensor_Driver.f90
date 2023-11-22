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
		real(kind=8), allocatable :: locations_m(:,:),locations_n(:,:) ! geometrical points
		integer,allocatable:: permutation_m(:,:),permutation_n(:,:)
		integer:: rank
		integer,allocatable:: Nunk_m(:),Nunk_n(:)
		integer:: Ndim
		integer:: tst=1
		real(kind=8)::wavelen,zdist
	end type quant_app

contains

	!**** user-defined subroutine to sample Z_mn as full matrix
	subroutine Zelem_MD_User(Ndim, m,n,value,quant)
		use z_BPACK_DEFS
		implicit none

		class(*),pointer :: quant
		integer, INTENT(IN):: Ndim
		integer, INTENT(IN):: m(Ndim),n(Ndim)
		complex(kind=8)::value
		integer ii, dim_i

		real(kind=8)::pos_o(Ndim),pos_s(Ndim), dist, waven

		select TYPE(quant)
		type is (quant_app)
			do dim_i=1,Ndim
				pos_o(dim_i) = quant%locations_m(dim_i,m(dim_i))
				pos_s(dim_i) = quant%locations_n(dim_i,n(dim_i))
			enddo
			if(quant%tst==1)then
				dist = sqrt(sum((pos_o-pos_s)**2d0))
			elseif(quant%tst==2)then
				dist = sqrt(sum((pos_o-pos_s)**2d0) + quant%zdist**2d0)
			elseif(quant%tst==3)then
				dist = sqrt(sum((pos_o-pos_s)**2d0))
			else
				write(*,*)'tst unknown'
			endif

			waven=2*BPACK_pi/quant%wavelen
			value = EXP(-BPACK_junit*waven*dist)/dist

		class default
			write(*,*)"unexpected type"
			stop
		end select
	end subroutine Zelem_MD_User


	!**** user-defined subroutine to sample Z_mn as full matrix (note that this is for the BF interface not for the BPACK interface)
	subroutine ZBelem_MD_User(Ndim, m, n,value_e,quant)
		use z_BPACK_DEFS
		implicit none

		class(*),pointer :: quant
		integer, INTENT(IN):: Ndim
		integer, INTENT(IN):: m(Ndim),n(Ndim)
		complex(kind=8)::value_e
		integer ii,dim_i, m1(Ndim),n1(Ndim)

		if(m(1)>0)then
			m1=m
			n1=-n
		else
			m1=n
			n1=-m
		endif

		!!! m,n still need to convert to the original order, using new2old of mshr and mshc
		select TYPE(quant)
		type is (quant_app)
			do dim_i=1,Ndim
				m1(dim_i)=quant%permutation_m(m1(dim_i),dim_i)
				n1(dim_i)=quant%permutation_n(n1(dim_i),dim_i)
			enddo
		class default
			write(*,*)"unexpected type"
			stop
		end select

		call Zelem_MD_User(Ndim,m1,n1,value_e,quant)
	end subroutine ZBelem_MD_User


end module APPLICATION_MODULE_FULL_SIMPLE


PROGRAM ButterflyPACK_RankBenchmark
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

    integer rank,ii,jj,kk
	real(kind=8),allocatable:: datain(:)
	real(kind=8) :: wavelen, ds, ppw
	integer :: ierr
	type(z_Hoption),target::option
	type(z_Hstat),target::stats
	type(z_mesh),target,allocatable::msh(:)
	type(z_kernelquant),target::ker
	type(quant_app),target::quant
	type(z_Bmatrix),target::bmat
	integer,allocatable:: groupmembers(:)
	integer nmpi, Nperdim, dims(3), inds(3)
	integer level,Maxlevel,m,n
	type(z_proctree),target::ptree
	integer,allocatable::Permutation(:)
	integer,allocatable:: Nunk_m_loc(:), Nunk_n_loc(:)
	integer,allocatable::tree(:),tree_m(:),tree_n(:)
	complex(kind=8),allocatable::rhs_glo(:,:),rhs_loc(:,:),x_glo(:,:),x_loc(:,:)
	integer nrhs
	type(z_matrixblock_MD) ::blocks
	character(len=1024)  :: strings,strings1
	integer flag,nargs,dim_i

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
	! option%format=  HODLR! HMAT!   ! the hierarhical format
	option%near_para=2.01d0        ! admissibiltiy condition, not referenced if option%format=  HODLR
	option%verbosity=1             ! verbosity level
	option%LRlevel=0             ! 0: low-rank compression 100: butterfly compression
	option%format=HSS_MD           ! currently this is the only format supported in MD

	! geometry points available
	option%xyzsort=TM ! no reordering will be perfomed
	option%knn=0   ! neareat neighbour points per geometry point, which helps improving the compression accuracy

	quant%tst = 2
	quant%wavelen = 0.25d0/8d0


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
						elseif(trim(strings)=='--wavelen')then
							read(strings1,*)quant%wavelen
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
! Read a full non-square matrix and do a BF compression

	ppw=5

    ds = quant%wavelen/ppw
    if(quant%tst==1)then ! two colinear plate
	  quant%Ndim = 2
	  quant%zdist = 0
      Nperdim = z_ceiling_safe(1d0/ds)
	  allocate(quant%Nunk_m(quant%Ndim))
	  allocate(quant%Nunk_n(quant%Ndim))
      quant%Nunk_m = Nperdim
      quant%Nunk_n = Nperdim
	  allocate(quant%locations_m(quant%Ndim,Nperdim))
	  allocate(quant%locations_n(quant%Ndim,Nperdim))
	  do m=1,Nperdim
		quant%locations_m(1,m)=m*ds+2
		quant%locations_m(2,m)=m*ds
	  enddo
	  do n=1,Nperdim
		quant%locations_n(1,n)=n*ds
		quant%locations_n(2,n)=n*ds
	  enddo

    elseif(quant%tst==2)then ! two parallel plate
	  quant%Ndim = 2
	  quant%zdist = 1d0
      Nperdim = z_ceiling_safe(1d0/ds)
	  allocate(quant%Nunk_m(quant%Ndim))
	  allocate(quant%Nunk_n(quant%Ndim))
      quant%Nunk_m = Nperdim
      quant%Nunk_n = Nperdim
	  allocate(quant%locations_m(quant%Ndim,Nperdim))
	  allocate(quant%locations_n(quant%Ndim,Nperdim))
	  dims = Nperdim
	  do m=1,Nperdim
		quant%locations_m(1,m)=m*ds
		quant%locations_m(2,m)=m*ds
	  enddo
	  do n=1,Nperdim
		quant%locations_n(1,n)=n*ds
		quant%locations_n(2,n)=n*ds
	  enddo

	elseif(quant%tst==3)then ! two 3D cubes

	  quant%Ndim = 3
      Nperdim = z_ceiling_safe(1d0/ds)
	  allocate(quant%Nunk_m(quant%Ndim))
	  allocate(quant%Nunk_n(quant%Ndim))
      quant%Nunk_m = Nperdim
      quant%Nunk_n = Nperdim
	  allocate(quant%locations_m(quant%Ndim,Nperdim))
	  allocate(quant%locations_n(quant%Ndim,Nperdim))
	  do m=1,Nperdim
		quant%locations_m(1,m)=m*ds
		quant%locations_m(2,m)=m*ds
		quant%locations_m(3,m)=m*ds
	  enddo
	  do n=1,Nperdim
		quant%locations_n(1,n)=n*ds+2
		quant%locations_n(2,n)=n*ds
		quant%locations_n(3,n)=n*ds
	  enddo
	endif

	allocate(Nunk_m_loc(quant%Ndim))
	allocate(Nunk_n_loc(quant%Ndim))
	allocate(msh(quant%Ndim))

	if(ptree%MyID==Main_ID)then
	write (*,*) ''
	write (*,*) 'RankBenchmark(Tensor) computing'
	write (*,*) 'Tensor size:', quant%Nunk_m, quant%Nunk_n
	write (*,*) ''
	endif
	!***********************************************************************

	!**** register the user-defined function and type in ker
	ker%QuantApp => quant
	ker%FuncZmn_MD => ZBelem_MD_User

	allocate(quant%Permutation_m(maxval(quant%Nunk_m),quant%Ndim))
	allocate(quant%Permutation_n(maxval(quant%Nunk_n),quant%Ndim))

	call z_BF_MD_Construct_Init(quant%Ndim, quant%Nunk_m, quant%Nunk_n, Nunk_m_loc, Nunk_n_loc, quant%Permutation_m, quant%Permutation_n, blocks, option, stats, msh, ker, ptree, Coordinates_m=quant%locations_m,Coordinates_n=quant%locations_n)
	call MPI_Bcast(quant%Permutation_m,maxval(quant%Nunk_m)*quant%Ndim,MPI_integer,0,ptree%comm,ierr)
	call MPI_Bcast(quant%Permutation_n,maxval(quant%Nunk_n)*quant%Ndim,MPI_integer,0,ptree%comm,ierr)

	call z_BF_MD_Construct_Element_Compute(quant%Ndim, blocks, option, stats, msh, ker, ptree)

!******************************************************************************!

	!**** print statistics
	call z_PrintStat(stats,ptree)

	call z_delete_proctree(ptree)
	call z_delete_Hstat(stats)
	do dim_i=1,quant%Ndim
		call z_delete_mesh(msh(dim_i))
	enddo
	deallocate(msh)
	call z_delete_kernelquant(ker)


    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	call z_blacs_exit_wrp(1)
	call MPI_Finalize(ierr)

end PROGRAM ButterflyPACK_RankBenchmark



