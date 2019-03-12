! “ButterflyPACK” Copyright (c) 2018, The Regents of the University of California, through
! Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.

! If you have questions about your rights to use or distribute this software, please contact
! Berkeley Lab's Intellectual Property Office at  IPO@lbl.gov.

! NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the
! U.S. Government consequently retains certain rights. As such, the U.S. Government has been
! granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable
! worldwide license in the Software to reproduce, distribute copies to the public, prepare
! derivative works, and perform publicly and display publicly, and to permit other to do so.

! Developers: Yang Liu
!             (Lawrence Berkeley National Lab, Computational Research Division).


! This exmple works with double-complex precision data
#define DAT 0

#include "ButterflyPACK_config.fi"



PROGRAM ButterflyPACK_IE_2D
    use BPACK_DEFS
    use EMCURV_MODULE

	use BPACK_structure
	use BPACK_Solve_Mul
	use BPACK_factor
	use BPACK_constr
	use BPACK_Utilities
	use omp_lib
	use misc
    implicit none

	! include "mkl_vml.fi"

    real(kind=8) para
    real(kind=8) tolerance
    integer Primary_block, nn, mm,kk,mn,rank,ii,jj
    integer i,j,k, threads_num
	integer seed_myid(50)
	integer times(8)
	integer edge
	real(kind=8) t1,t2,t3, x,y,z,r,theta,phi
	complex(kind=8),allocatable:: matU(:,:),matV(:,:),matZ(:,:),LL(:,:),RR(:,:),matZ1(:,:)

	character(len=:),allocatable  :: string
	character(len=1024)  :: strings
	character(len=6)  :: info_env
	integer :: length
	integer :: ierr
	integer*8 oldmode,newmode
	type(Hoption)::option,option1,option2
	type(Hstat)::stats,stats1,stats2
	type(mesh)::msh,msh1,msh2
	type(kernelquant)::ker,ker1,ker2
	type(Bmatrix)::bmat,bmat1,bmat2
	integer,allocatable:: groupmembers(:)
	integer nmpi
	type(proctree)::ptree,ptree1,ptree2
	type(quant_EMCURV),target::quant
	CHARACTER (LEN=1000) DATA_DIR
	integer:: randsize=50
	real(kind=8),allocatable::xyz(:,:)
	integer,allocatable::Permutation(:),tree(:)
	integer Nunk_loc,Maxlevel

	integer maxn, maxnev, maxncv, ldv
	integer iparam(11), ipntr(14)
    logical,allocatable:: select(:)
    complex(kind=8),allocatable:: ax(:), mx(:), d(:), v(:,:), workd(:), workev(:), resid(:), workl(:)
    real(kind=8),allocatable:: rwork(:), rd(:,:)

    character bmattype
	character(len=10) which
    integer ido, n, nx, nev, ncv, lworkl, info, nconv, maxitr, ishfts, mode
    complex(kind=8) sigma
    real(kind=8) tol
    logical rvec
	real(kind=8),external :: pdznorm2, dlapy2


	! nmpi and groupmembers should be provided by the user
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo

	! generate the process tree
	call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree)
	deallocate(groupmembers)


	if(ptree%MyID==Main_ID)then
    write(*,*) "-------------------------------Program Start----------------------------------"
    write(*,*) "ButterflyPACK_IE_2D"
    write(*,*) "   "
	endif

	!**** initialize stats and option
	call InitStat(stats)
	call SetDefaultOptions(option)

	!**** intialize the user-defined derived type quant
	quant%RCS_static=1
    quant%RCS_Nsample=2000
	quant%model2d=10
	quant%wavelength=0.08
	quant%Nunk=5000

	!**** default parameters for the eigen solvers
	quant%CMmode=0
	quant%SI=0
	quant%shift=(0d0, 0d0)	
	
	option%ErrSol=1
	option%format=  HODLR !HMAT!  HODLR !
	option%near_para=0.01d0
	option%verbosity=2
	option%ILU=0
	option%forwardN15flag=0
        ! option%schulzlevel=0
        ! option%LRlevel=100
       ! option%level_check=1
    option%tol_itersol=1d-5
    ! option%sample_para=4d0

	
	
	
	if(iargc()>=1)then
		call getarg(1,strings)
		read(strings,*)option%LR_BLK_NUM
	endif
	if(iargc()>=2)then
		call getarg(2,strings)
		read(strings,*)quant%model2d
	endif
	if(iargc()>=3)then
		call getarg(3,strings)
		read(strings,*)quant%Nunk
	endif
	if(iargc()>=4)then
		call getarg(4,strings)
		read(strings,*)quant%wavelength
	endif
	if(iargc()>=5)then
		call getarg(5,strings)
		read(strings,*)option%tol_comp
		option%tol_rand=option%tol_comp
		option%tol_Rdetect=option%tol_comp*1d-1
	endif
	if(iargc()>=6)then
		call getarg(6,strings)
		read(strings,*)option%ErrFillFull
	endif
	if(iargc()>=7)then
		call getarg(7,strings)
		read(strings,*)option%RecLR_leaf
	endif
	if(iargc()>=8)then
		call getarg(8,strings)
		read(strings,*)option%BACA_Batch
	endif
	if(iargc()>=9)then
		call getarg(9,strings)
		read(strings,*)option%LRlevel
	endif
	if(iargc()>=10)then
		call getarg(10,strings)
		read(strings,*)option%precon
	endif
	if(iargc()>=11)then
		call getarg(11,strings)
		read(strings,*)option%xyzsort
	endif
	if(iargc()>=12)then
		call getarg(12,strings)
		read(strings,*)option%Nmin_leaf
	endif
	if(iargc()>=13)then
		call getarg(13,strings)
		read(strings,*)option%near_para
	endif

    if(iargc()>=14)then
        call getarg(14,strings)
        read(strings,*)option%pat_comp
    endif

    if(iargc()>=15)then
        call getarg(15,strings)
        read(strings,*)option%schulzlevel
    endif

    if(iargc()>=16)then
        call getarg(16,strings)
        read(strings,*)option%Nbundle
    endif

    if(iargc()>=17)then
        call getarg(17,strings)
        read(strings,*)option%format
    endif

    if(iargc()>=18)then
        call getarg(18,strings)
        read(strings,*)quant%CMmode
    endif	
	
    if(iargc()>=19)then
        call getarg(19,strings)
        read(strings,*)quant%SI
    endif	
	
    ! if(iargc()>=20)then
        ! call getarg(20,strings)
        ! read(strings,*)quant%shift
    ! endif	
		
	


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


    !***********************************************************************
	!**** construct compressed A - sigma I or  A - sigma real(A)

	!**** geometry generalization and discretization
    call geo_modeling_CURV(quant,ptree%Comm)
	!**** register the user-defined function and type in ker
	ker%QuantApp => quant
	ker%FuncZmn => Zelem_EMCURV_Shifted
	!**** initialization of the construction phase
	allocate(xyz(2,quant%Nunk))
	do edge=1, quant%Nunk
		xyz(:,edge) = quant%xyz(:,edge*2-1)
	enddo
    allocate(Permutation(quant%Nunk))
	call BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat,option,stats,msh,ker,ptree,Coordinates=xyz)
	deallocate(Permutation) ! caller can use this permutation vector if needed
	deallocate(xyz)
	!**** computation of the construction phase
    call BPACK_construction_Element(bmat,option,stats,msh,ker,element_Zmn_user,ptree)
	!**** factorization phase
    call BPACK_Factorization(bmat,option,stats,ptree,msh)
	!**** print statistics
	call PrintStat(stats,ptree)

    !***********************************************************************
	!**** construct compressed A

	call CopyOptions(option,option1)
	ker1%QuantApp => quant
	ker1%FuncZmn => Zelem_EMCURV
	msh1%Nunk = msh%Nunk
	!**** initialize stats and option
	call InitStat(stats1)
	!**** create the process tree, can use larger number of mpis if needed
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo
	call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree1)
	deallocate(groupmembers)
	!**** use the clustering tree from the first HODLR
	select case(option%format)
	case(HODLR)
		Maxlevel=bmat%ho_bf%Maxlevel
	case(HMAT)
		Maxlevel=bmat%h_mat%Maxlevel
	end select
	allocate (tree(2**Maxlevel))
	do ii=1,2**Maxlevel
		tree(ii)=msh%basis_group(2**Maxlevel+ii-1)%tail-msh%basis_group(2**Maxlevel+ii-1)%head+1
	enddo

	allocate(xyz(2,quant%Nunk))
	do edge=1, quant%Nunk
		xyz(:,edge) = quant%xyz(:,edge*2-1)
	enddo
	allocate(Permutation(quant%Nunk))
	call BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat1,option1,stats1,msh1,ker1,ptree1,Coordinates=xyz)
	deallocate(Permutation)
	deallocate(tree)
	deallocate(xyz)
	call BPACK_construction_Element(bmat1,option1,stats1,msh1,ker1,element_Zmn_user,ptree1)
	!**** print statistics
	call PrintStat(stats1,ptree1)

	! call FULLMAT_Element(option1,stats1,msh1,ker1,element_Zmn_user,ptree)


	!***********************************************************************



    !***********************************************************************
	!**** construct compressed real(A)
	if(quant%CMmode==1)then ! solve the characteristic mode
		call CopyOptions(option,option2)
		ker2%QuantApp => quant
		ker2%FuncZmn => Zelem_EMCURV_Real
		msh2%Nunk = msh%Nunk
		!**** initialize stats and option
		call InitStat(stats2)
		!**** create the process tree, can use larger number of mpis if needed
		allocate(groupmembers(nmpi))
		do ii=1,nmpi
			groupmembers(ii)=(ii-1)
		enddo
		call CreatePtree(nmpi,groupmembers,MPI_Comm_World,ptree2)
		deallocate(groupmembers)
		!**** use the clustering tree from the first HODLR
		select case(option%format)
		case(HODLR)
			Maxlevel=bmat%ho_bf%Maxlevel
		case(HMAT)
			Maxlevel=bmat%h_mat%Maxlevel
		end select
		allocate (tree(2**Maxlevel))
		do ii=1,2**Maxlevel
			tree(ii)=msh%basis_group(2**Maxlevel+ii-1)%tail-msh%basis_group(2**Maxlevel+ii-1)%head+1
		enddo

		allocate(xyz(2,quant%Nunk))
		do edge=1, quant%Nunk
			xyz(:,edge) = quant%xyz(:,edge*2-1)
		enddo
		allocate(Permutation(quant%Nunk))
		call BPACK_construction_Init(quant%Nunk,Permutation,Nunk_loc,bmat2,option2,stats2,msh2,ker2,ptree2,Coordinates=xyz)
		deallocate(Permutation)
		deallocate(tree)
		deallocate(xyz)
		call BPACK_construction_Element(bmat2,option2,stats2,msh2,ker2,element_Zmn_user,ptree2)
		!**** factorization phase
		call BPACK_Factorization(bmat2,option2,stats2,ptree2,msh2)
		!**** print statistics
		call PrintStat(stats2,ptree2)
	endif
	!***********************************************************************




	! !**** solve phase
   ! call EM_solve_CURV(bmat,option,msh,quant,ptree,stats)






#ifdef HAVE_ARPACK

	  maxn = quant%Nunk
	  ldv = Nunk_loc
	  maxnev = 4
	  maxncv = 20

	  n = maxn
      nev   = maxnev
      ncv   = maxncv
      lworkl  = 3*ncv**2+5*ncv
      tol    = 1D-13 !option%tol_comp
      ido    = 0
      info   = 0

      ishfts = 1
      maxitr = 3000
      ! mode   = 3
      iparam(1) = ishfts
      iparam(3) = maxitr
      ! iparam(7) = mode


	  allocate(select(maxncv))
	  allocate(ax(Nunk_loc))
	  allocate(mx(Nunk_loc))
	  allocate(d(maxn))
	  allocate(v(ldv,maxncv))
	  allocate(workd(3*Nunk_loc))
	  allocate(workev(3*maxncv))
	  allocate(resid(maxn))
	  allocate(workl(3*maxncv*maxncv+5*maxncv))
	  allocate(rwork(maxncv))
	  allocate(rd(maxncv,3))

	if(quant%CMmode==0)then ! solve the eigen mode

		do while(ido/=99)

			if(quant%SI==0)then  ! regular mode
				bmattype  = 'I'
				which='LM'  ! largest eigenvalues
				iparam(7) = 1
				sigma=0d0
				call pznaupd(ptree%Comm, ido, bmattype, Nunk_loc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,rwork, info )
				if(ido == -1 .or. ido == 1)then !Perform  y <--- OP*x = A*x
					call BPACK_Mult('N',Nunk_loc,1,workd(ipntr(1)),workd(ipntr(2)),bmat1,ptree,option1,stats1)
				endif
			else                ! shift-invert mode
				bmattype  = 'I'
				which='LM' ! eigenvalues closest to the shifts
				iparam(7) = 3
				sigma = quant%shift
				call pznaupd(ptree%Comm, ido, bmattype, Nunk_loc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,rwork, info )
				if(ido == -1 .or. ido == 1)then !Perform  y <--- OP*x = inv[A-SIGMA*I]*x
					call BPACK_Solution(bmat,workd(ipntr(2)),workd(ipntr(1)),Nunk_loc,1,option,ptree,stats)
				endif
			endif

		enddo

      if ( info < 0 ) then
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation of _naupd.'
         print *, ' '
      else
         rvec = .true.
         call pzneupd  (ptree%Comm,rvec, 'A', select, d, v, ldv, sigma, workev, bmattype, Nunk_loc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, ierr)
         if ( ierr /= 0) then
             print *, ' '
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
         else
             nconv = iparam(5)
             do j=1, nconv

				call BPACK_Mult('N',Nunk_loc,1,v(1,j),ax,bmat1,ptree1,option1,stats1)
				mx(1:Nunk_loc) = v(1:Nunk_loc,j)
                call zaxpy (Nunk_loc, -d(j), mx, 1, ax, 1)
                rd(j,1) = dble (d(j))
                rd(j,2) = dimag (d(j))
                rd(j,3) = pdznorm2(ptree%Comm, Nunk_loc, ax, 1)
                rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
			enddo
             call pdmout(ptree%Comm, 6, nconv, 3, rd, maxncv, -6,'Ritz values (Real, Imag) and direct residuals')
          end if

		 if (ptree%MyID==Main_ID)then
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' '
             print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
             print *, ' '
         end if
         print *, ' '
         print *, '_NDRV1'
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of processors is ', ptree%nproc
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', nconv
         print *, ' The number of Implicit Arnoldi update',' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
         endif
      end if

	else if(quant%CMmode==1)then ! solve the characteristic mode
		do while(ido/=99)

			if(quant%SI==0)then ! regular mode
				bmattype  = 'G'
				which='LM'	 ! largest eigenvalues
				iparam(7) = 2
				sigma = 0
				call pznaupd(ptree%Comm, ido, bmattype, Nunk_loc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,rwork, info )
				if(ido == -1 .or. ido == -1)then !Perform  y <--- OP*x = inv[M]*A*x
					call BPACK_Mult('N',Nunk_loc,1,workd(ipntr(1)),workd(ipntr(3)),bmat1,ptree1,option1,stats1)
					call BPACK_Solution(bmat2,workd(ipntr(2)),workd(ipntr(3)),Nunk_loc,1,option2,ptree2,stats2)
				else if(ido==2)then !Perform  y <--- M*x
					call BPACK_Mult('N',Nunk_loc,1,workd(ipntr(1)),workd(ipntr(2)),bmat2,ptree2,option2,stats2)
				endif
			else ! shift-invert mode
				bmattype  = 'G'
				which='LM'	! eigenvalues closest to the shifts
				iparam(7) = 3
				sigma = quant%shift
				call pznaupd(ptree%Comm, ido, bmattype, Nunk_loc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,rwork, info )
				if(ido == -1)then !Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x
					call BPACK_Mult('N',Nunk_loc,1,workd(ipntr(1)),workd(ipntr(3)),bmat2,ptree2,option2,stats2)
					call BPACK_Solution(bmat,workd(ipntr(2)),workd(ipntr(3)),Nunk_loc,1,option,ptree,stats)
				else if(ido==1)then !Perform y <-- OP*x = inv[A-sigma*M]*M*x, M*x has been saved in workd(ipntr(3))
					call BPACK_Solution(bmat,workd(ipntr(2)),workd(ipntr(3)),Nunk_loc,1,option,ptree,stats)
				else if(ido==2)then !Perform  y <--- M*x
					call BPACK_Mult('N',Nunk_loc,1,workd(ipntr(1)),workd(ipntr(2)),bmat2,ptree2,option2,stats2)
				endif
			endif
		enddo

      if ( info < 0 ) then
         print *, ' '
         print *, ' Error with _naupd, info = ', info
         print *, ' Check the documentation of _naupd.'
         print *, ' '
      else
         rvec = .true.
         call pzneupd  (ptree%Comm,rvec, 'A', select, d, v, ldv, sigma, workev, bmattype, Nunk_loc, which, nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, rwork, ierr)
         if ( ierr /= 0) then
             print *, ' '
             print *, ' Error with _neupd, info = ', ierr
             print *, ' Check the documentation of _neupd. '
             print *, ' '
         else
             nconv = iparam(5)
             do j=1, nconv

				call BPACK_Mult('N',Nunk_loc,1,v(1,j),ax,bmat1,ptree1,option1,stats1)
				call BPACK_Mult('N',Nunk_loc,1,v(1,j),mx,bmat2,ptree2,option2,stats2)
                call zaxpy (Nunk_loc, -d(j), mx, 1, ax, 1)
                rd(j,1) = dble (d(j))
                rd(j,2) = dimag (d(j))
                rd(j,3) = pdznorm2(ptree%Comm, Nunk_loc, ax, 1)
                rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
			enddo
             call pdmout(ptree%Comm, 6, nconv, 3, rd, maxncv, -6,'Ritz values (Real, Imag) and direct residuals')
          end if

		 if (ptree%MyID==Main_ID)then
         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' '
             print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
             print *, ' '
         end if
         print *, ' '
         print *, '_NDRV1'
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of processors is ', ptree%nproc
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', nconv
         print *, ' The number of Implicit Arnoldi update',' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
         endif
      end if


	endif


	  deallocate(select)
	  deallocate(ax)
	  deallocate(mx)
	  deallocate(d)
	  deallocate(v)
	  deallocate(workd)
	  deallocate(workev)
	  deallocate(resid)
	  deallocate(workl)
	  deallocate(rwork)
	  deallocate(rd)

#endif


	!**** deletion of quantities
	call delete_quant_EMCURV(quant)
	call delete_proctree(ptree)
	call delete_Hstat(stats)
	call delete_mesh(msh)
	call delete_kernelquant(ker)
	call BPACK_delete(bmat)

	call delete_proctree(ptree1)
	call delete_Hstat(stats1)
	call delete_mesh(msh1)
	call delete_kernelquant(ker1)
	call BPACK_delete(bmat1)

	if(quant%CMmode==1)then ! solve the characteristic mode
		call delete_proctree(ptree2)
		call delete_Hstat(stats2)
		call delete_mesh(msh2)
		call delete_kernelquant(ker2)
		call BPACK_delete(bmat2)
	endif

    if(ptree%MyID==Main_ID .and. option%verbosity>=0)write(*,*) "-------------------------------program end-------------------------------------"

	call blacs_exit(1)
	call MPI_Finalize(ierr)

end PROGRAM ButterflyPACK_IE_2D




