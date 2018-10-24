module APPLICATION_MODULE
use z_HODLR_DEFS
implicit none

	!**** define your application-related variables here   
	type quant_app
		complex(kind=8), allocatable :: matZ_glo(:,:) ! Full Matrix: Full matrix read from files	
		complex(kind=8), allocatable :: matZ_loc(:,:) ! Local Matrix: Loccal matrix in a npx1 blasc grid	
		integer,pointer :: N_p(:,:) ! column sizes of all processes sharing this hodlr
		type(z_hobf),pointer::ho_bf ! Use this precomputed hodbf as matvec
	end type quant_app

contains
	
	!**** user-defined subroutine to sample Z_mn as full matrix
	subroutine Z_elem_FULL(m,n,value_e,quant)
		use z_HODLR_DEFS
		implicit none 
		
		class(*),pointer :: quant
		integer, INTENT(IN):: m,n
		complex(kind=8)::value_e 
		integer ii
		
		select TYPE(quant)
		type is (quant_app)
			value_e = quant%matZ_glo(m,n)
		class default
			write(*,*)"unexpected type"
			stop
		end select	
	end subroutine Z_elem_FULL
	


	

	subroutine HODLR_MVP_OneHODLR(trans,Nloc,num_vect,Vin,Vout,msh,ptree,stats,quant)
		use z_HODLR_DEFS
		use z_DenseLA
		use z_misc
		use z_HODLR_Solve_Mul
		implicit none 
		character trans
		complex(kind=8) Vin(:,:),Vout(:,:)
		complex(kind=8),allocatable:: Vin_tmp(:,:),Vout_tmp(:,:),Vin_tmp_2D(:,:),Vout_tmp_2D(:,:)
		complex(kind=8) ctemp,a,b
		integer ii,jj,nn,fl_transpose,kk,black_step
		integer, INTENT(in)::Nloc,num_vect
		real(kind=8) n1,n2,tmp(2)
		type(z_mesh)::msh
		type(z_proctree)::ptree
		integer idxs_o,idxe_o,N
		integer nproc,ctxt,info,nb1Dc, nb1Dr, level_p,pgno,num_blocks,ii_new,gg,proc,myi,myj,myAcols,myArows,nprow,npcol,myrow,mycol,Nrow,Ncol
		integer::descsVin(9),descsVout(9),descsMat2D(9),descsVin2D(9),descsVout2D(9)
		class(*),pointer :: quant
		type(z_hobf),pointer::ho_bf
		type(z_Hstat)::stats
		
		pgno=1
		nproc = ptree%pgrp(pgno)%nproc
		
		select TYPE(quant)   
		type is (quant_app)
			ho_bf=>quant%ho_bf
			call z_MVM_Z_forward(trans,Nloc,num_vect,1,ho_bf%Maxlevel+1,Vin,Vout,ho_bf,ptree,stats)	
		end select
		
	end subroutine HODLR_MVP_OneHODLR




	
	subroutine HODLR_MVP_Fullmat(trans,Nloc,num_vect,Vin,Vout,msh,ptree,stats,quant)
		use z_HODLR_DEFS
		use z_DenseLA
		use z_misc
		implicit none 
		character trans
		complex(kind=8) Vin(:,:),Vout(:,:)
		complex(kind=8),allocatable:: Vin_tmp(:,:),Vout_tmp(:,:),Vin_tmp_2D(:,:),Vout_tmp_2D(:,:)
		complex(kind=8) ctemp,a,b
		integer ii,jj,nn,fl_transpose,kk,black_step
		integer, INTENT(in)::Nloc,num_vect
		real(kind=8) n1,n2,tmp(2)
		type(z_mesh)::msh
		type(z_proctree)::ptree
		type(z_Hstat)::stats
		integer idxs_o,idxe_o,N
		integer nproc,ctxt,info,nb1Dc, nb1Dr, level_p,pgno,num_blocks,ii_new,gg,proc,myi,myj,myAcols,myArows,nprow,npcol,myrow,mycol,Nrow,Ncol
		integer::descsVin(9),descsVout(9),descsMat2D(9),descsVin2D(9),descsVout2D(9)
		class(*),pointer :: quant

		pgno=1
		nproc = ptree%pgrp(pgno)%nproc
		
		
		select TYPE(quant)   
		type is (quant_app)
		
		N = quant%N_p(nproc,2)
		
					
		!!!!**** generate 2D grid blacs quantities
		ctxt = ptree%pgrp(pgno)%ctxt		
		call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)	
		if(myrow/=-1 .and. mycol/=-1)then
			myArows = z_numroc_wp(N, nbslpk, myrow, 0, nprow)
			myAcols = z_numroc_wp(num_vect, nbslpk, mycol, 0, npcol)
			call descinit( descsVin2D, N, num_vect, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
			call descinit( descsVout2D, N, num_vect, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )
			myArows = z_numroc_wp(N, nbslpk, myrow, 0, nprow)
			myAcols = z_numroc_wp(N, nbslpk, mycol, 0, npcol)			
			call descinit( descsMat2D, N, N, nbslpk, nbslpk, 0, 0, ctxt, max(myArows,1), info )					
			allocate(Vin_tmp_2D(myArows,myAcols))
			allocate(Vout_tmp_2D(myArows,myAcols))
			Vout_tmp_2D=0			
		else 
			descsVin2D(2)=-1
			descsVout2D(2)=-1
			descsMat2D(2)=-1
		endif

		
		!!!!**** redistribution of input vectors
		call z_Redistribute1Dto2D(Vin,quant%N_p,0,pgno,Vin_tmp_2D,N,0,pgno,num_vect,ptree)	
			

		!!!!**** perform gemm on 2d grid
		if(myrow/=-1 .and. mycol/=-1)then
			call z_pgemmf90(trans,'N',N,num_vect,N,cone, quant%matZ_loc,1,1,descsMat2D,Vin_tmp_2D,1,1,descsVin2D,czero,Vout_tmp_2D,1,1,descsVout2D)
		endif
		
	
		!!!!**** redistribution of output vectors
		call z_Redistribute2Dto1D(Vout_tmp_2D,N,0,pgno,Vout,quant%N_p,0,pgno,num_vect,ptree)	
		
		
		!!!!**** deallocation buffers
		if(myrow/=-1 .and. mycol/=-1)then
			deallocate(Vin_tmp_2D)
			deallocate(Vout_tmp_2D)
		endif
		
		
		end select
		
	end subroutine HODLR_MVP_Fullmat


	subroutine CreateDistDenseMat(N,msh,ptree,quant)
		use z_HODLR_DEFS
		use z_DenseLA
		use z_misc
		implicit none 
		complex(kind=8),allocatable:: Vin_tmp(:,:),Vout_tmp(:,:)
		complex(kind=8) ctemp,a,b
		integer ii,jj,nn,fl_transpose,kk,black_step
		integer, INTENT(in)::N
		real(kind=8) n1,n2,tmp(2)
		type(z_mesh)::msh
		type(z_proctree)::ptree
		type(quant_app) :: quant
		integer nproc,ctxt,nb1Dc, nb1Dr, level_p,pgno,num_blocks,ii_new,gg,proc,myi,myj,myAcols,myArows,nprow,npcol,myrow,mycol
		
		
		pgno=1
		nproc = ptree%pgrp(pgno)%nproc
	
		!!!!****** allocate index array for 1D HODLR layout
		level_p = ptree%nlevel-z_GetTreelevel(pgno)
		num_blocks = 2**level_p
		allocate(quant%N_p(nproc,2))
		quant%N_p(:,1) = N+1	
		quant%N_p(:,2) = -N-1	
		do ii=1,num_blocks
			ii_new=ii
			gg = 2**level_p+ii_new-1
			proc = ptree%pgrp(pgno*2**level_p+ii-1)%head - ptree%pgrp(pgno)%head
			quant%N_p(proc+1,1) = min(quant%N_p(proc+1,1),msh%basis_group(gg)%head)	
			quant%N_p(proc+1,2) = max(quant%N_p(proc+1,2),msh%basis_group(gg)%tail)				
		enddo

		
		!!!!****** assemble the full matrix in 1D blacs layout
		ctxt = ptree%pgrp(pgno)%ctxt
		call blacs_gridinfo(ctxt, nprow, npcol, myrow, mycol)			
		if(myrow/=-1 .and. mycol/=-1)then
		myArows = z_numroc_wp(N, nbslpk, myrow, 0, nprow)
		myAcols = z_numroc_wp(N, nbslpk, mycol, 0, npcol)	
		allocate(quant%matZ_loc(myArows,myAcols))
		quant%matZ_loc=0			
		do myi=1,myArows
			call z_l2g(myi,myrow,N,nprow,nbslpk,ii)
			do myj=1,myAcols
				call z_l2g(myj,mycol,N,npcol,nbslpk,jj)
				quant%matZ_loc(myi,myj) = quant%matZ_glo(msh%new2old(ii),msh%new2old(jj))			
			enddo
		enddo	
		endif
	
	end subroutine CreateDistDenseMat


	
	
end module APPLICATION_MODULE	


PROGRAM MLMDA_DIRECT_SOLVER_3D_CFIE

    use z_HODLR_DEFS
	
	use z_HODLR_structure
	use z_HODLR_factor
	use z_HODLR_constr
	use omp_lib
	use z_Bplus_compress
	use z_HODLR_randomMVP
	use APPLICATION_MODULE
	
    implicit none

    real(kind=8) para,error,tmp1,tmp2,norm1,norm2
    real(kind=8) tolerance,rankrate
    integer Primary_block, nn, mm
    integer i,j,k, threads_num,ii,jj
	integer seed_myid(50)
	integer times(8)	
	real(kind=8) t1,t2,t3,t4,x,y,z,r,theta,phi,tmp(3),Memory
	complex(kind=8),allocatable:: InputVec(:)
	complex(kind=8):: ctemp
	integer Ntunnel,kk,black_step,rank0
	complex(kind=8),allocatable::Vout1(:,:),Vout2(:,:),Vin(:,:)
	character(len=1024)  :: strings
	type(z_Hoption):: option
	type(z_Hstat)::stats,stats1	
	type(z_mesh)::msh,msh1	
	type(z_kernelquant)::ker
	integer:: explicitflag
	type(z_hobf),target::ho_bf,ho_bf1
	integer Nin1,Nout1,Nin2,Nout2	
	type(z_proctree)::ptree,ptree1
	integer,allocatable:: groupmembers(:)	
	integer :: ierr
	integer :: nmpi
	CHARACTER (LEN=1000) DATA_DIR	
	type(quant_app),target::quant
	integer N_unk_loc
	
	! nmpi and groupmembers should be provided by the user 
	call MPI_Init(ierr)
	call MPI_Comm_size(MPI_Comm_World,nmpi,ierr)
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo	
	
	call z_createptree(nmpi,groupmembers,MPI_Comm_World,ptree)
	deallocate(groupmembers)
	
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
	
    !real scale


    if(ptree%MyID==Main_ID)write(*,*) "-------------------------------Program Start----------------------------------"
    if(ptree%MyID==Main_ID)write(*,*) "MLMDA_DIRECT_SOLVER_3D_CFIE_NEW"
    if(ptree%MyID==Main_ID)write(*,*) "FOR X64 COMPILER"
    if(ptree%MyID==Main_ID)write(*,*) "   "

	call z_initstat(stats)
	call z_setdefaultoptions(option)
	
     !*************************input******************************

    !tol=0.000001
	!itmax=10000
	
	CALL getarg(1, DATA_DIR)
	

	option%Nmin_leaf=100
	option%tol_comp=1d-4
	option%tol_Rdetect=1d-4 !3d-5
	option%nogeo=0
	option%tol_LS=1d-10
	option%tol_itersol=3d-3
	option%n_iter=1000
	option%tol_rand=option%tol_comp !1d-4
	option%level_check=100
	option%precon=DIRECT
	option%xyzsort=CKD
	option%lnoBP=600
	option%TwoLayerOnly=1
	option%LRlevel=0
	option%RecLR_leaf=BACA
	option%rank0 = 64
	option%rankrate = 1.5d0	
	
	explicitflag=0



	
	
	!Nmin_leaf=250
	!para=0.01
	!tolerance=0.001
	!alpha=0.5
    !ker%wavelength=2.
    ! tolerance=ACA_tolerance_forward
    ! call OMP_set_dynamic(.true.)
    !call OMP_set_nested(.true.)

    !*********************************************************

    ! ker%omiga=2*pi/ker%wavelength/sqrt(mu0*eps0)
    ! ker%wavenum=2*pi/ker%wavelength

	!**** register the user-defined function and type in ker 
	ker%FuncZmn=>Z_elem_FULL
	ker%FuncMatVec=>HODLR_MVP_Fullmat
	ker%QuantZmn=>quant	
	
	Ntunnel = 15600
	msh%Nunk = 3720
	Nin1 = 320*2 
	Nout1 = 610*2
	Nin2 = 320*2 
	Nout2 = 610*2	
	call z_assert(Nin1+Nout1+Nin2+Nout2==msh%Nunk,'The two surfaces have mismatched number of unknowns')
	
	! predefine the first three levels of tree due to the physical meanings
	allocate(msh%pretree(4))
	msh%pretree(1) = Nin1
	msh%pretree(2) = Nout1
	msh%pretree(3) = Nin2
	msh%pretree(4) = Nout2
	
	
	if(ptree%MyID==Main_ID)write(*,*)'Blackbox HODLR for scattering matrix compression'
	
	if(ptree%MyID==Main_ID)write(*,'(A10,I9,A11,I9)')' Ntunnel: ',Ntunnel,' Nsurface: ',msh%Nunk
	
	allocate(msh%xyz(3,msh%Nunk))

	!!!!!! for simplicity, I'm duplicating locations of each point

	open(unit=521,file=trim(DATA_DIR)//'/edge_cen.out',status='old')
	do kk=1,Ntunnel/2
		read(521,*) tmp(1),tmp(2),tmp(3)
	end do          
	do kk=1,msh%Nunk/2
		read(521,*) msh%xyz(1,2*kk-1),msh%xyz(2,2*kk-1),msh%xyz(3,2*kk-1)
		msh%xyz(:,2*kk) = msh%xyz(:,2*kk-1)
	end do
	close(521)
	
	
	
	t1 = OMP_get_wtime()	
    if(ptree%MyID==Main_ID)write(*,*) "constructing HODLR formatting......"
    call z_HODLR_structuring(ho_bf,option,msh,ker,z_element_Zmn_user,ptree)
	call z_BPlus_structuring(ho_bf,option,msh,ptree)
    if(ptree%MyID==Main_ID)write(*,*) "HODLR formatting finished"
    if(ptree%MyID==Main_ID)write(*,*) "    "
	t2 = OMP_get_wtime()
	if(ptree%MyID==Main_ID)write(*,*)t2-t1,'secnds'	
	! stop 

	t1 = OMP_get_wtime()
	if(ptree%MyID==Main_ID)write(*,*) "Generating fullmat ......"
	allocate(quant%matZ_glo(msh%Nunk,msh%Nunk))
	quant%matZ_glo = 0
	
	open(unit=888,file=trim(DATA_DIR)//'/fullmat.out',status='old')
	do ii=1,msh%Nunk
	do kk=1,msh%Nunk
		read(888,*)tmp(1),tmp(2)
		quant%matZ_glo(kk,ii) = cmplx(tmp(1),tmp(2),dp) 
	end do
	end do
	close(unit=888)
	
	call CreateDistDenseMat(msh%Nunk,msh,ptree,quant)
	
	
	if(ptree%MyID==Main_ID)write(*,*) "Generating fullmat finished"
	t2 = OMP_get_wtime()   
	if(ptree%MyID==Main_ID)write(*,*)t2-t1, 'secnds'	

	
	if(explicitflag ==1)then

		t1 = OMP_get_wtime()	
		if(ptree%MyID==Main_ID)write(*,*) "HODLR construction......"
		call z_HODLR_construction(ho_bf,option,stats,msh,ker,z_element_Zmn_user,ptree)		
		if(ptree%MyID==Main_ID)write(*,*) "HODLR construction finished"
		if(ptree%MyID==Main_ID)write(*,*) "    "
		t2 = OMP_get_wtime()   
		if(ptree%MyID==Main_ID)write(*,*)t2-t1, 'secnds'			
	
	
		N_unk_loc = msh%idxe-msh%idxs+1
		allocate(Vin(N_unk_loc,1))
		allocate(Vout1(N_unk_loc,1))
		allocate(Vout2(N_unk_loc,1))
		Vout2=0
		do ii=1,N_unk_loc
			call z_random_dp_number(Vin(ii,1))
		end do
		
		
		call z_MVM_Z_forward('N',N_unk_loc,1,1,ho_bf%Maxlevel+1,Vin,Vout1,ho_bf,ptree,stats)
	
		call z_matvec_user('N',N_unk_loc,1,Vin,Vout2,msh,ptree,stats,ker)
		
		
		tmp1 = z_fnorm(Vout2-Vout1,N_unk_loc,1)**2d0
		call MPI_ALLREDUCE(tmp1, norm1, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
		tmp2 = z_fnorm(Vout2,N_unk_loc,1)**2d0
		call MPI_ALLREDUCE(tmp2, norm2, 1,MPI_double_precision, MPI_SUM, ptree%Comm,ierr)
		error = sqrt(norm1)/sqrt(norm2)
		
		
		deallocate(Vin,Vout1,Vout2)
		
		if(ptree%MyID==Main_ID)write(*,*)error,'accuracy of construction'

		
	else if(explicitflag ==0)then
		N_unk_loc = msh%idxe-msh%idxs+1
		t1 = OMP_get_wtime()	
		if(ptree%MyID==Main_ID)write(*,*) "FullMATVEC-based HODLR construction......"		
		call z_HODLR_randomized(ho_bf,z_matvec_user,N_unk_loc,Memory,error,option,stats,ker,ptree,msh)
		t2 = OMP_get_wtime()  
		if(ptree%MyID==Main_ID)write(*,*) "FullMATVEC-based HODLR construction finished",t2-t1, 'secnds. Error: ', error	

	end if	
	
	
	


	option%nogeo=1
	option%xyzsort=NATURAL
	msh1%Nunk = msh%Nunk
	ker%FuncMatVec=>HODLR_MVP_OneHODLR
	quant%ho_bf=>ho_bf
	
	call z_initstat(stats1)
	
	allocate(groupmembers(nmpi))
	do ii=1,nmpi
		groupmembers(ii)=(ii-1)
	enddo	
	call z_createptree(nmpi,groupmembers,MPI_Comm_World,ptree1)
	deallocate(groupmembers)
	
	allocate (msh1%pretree(2**ho_bf%Maxlevel))	
	do ii=1,2**ho_bf%Maxlevel
		msh1%pretree(ii)=msh%basis_group(2**ho_bf%Maxlevel+ii-1)%tail-msh%basis_group(2**ho_bf%Maxlevel+ii-1)%head+1
	enddo
	
    call z_HODLR_structuring(ho_bf1,option,msh1,ker,z_element_Zmn_user,ptree1)
	call z_BPlus_structuring(ho_bf1,option,msh1,ptree1)	
	
	N_unk_loc = msh1%idxe-msh1%idxs+1
	t1 = OMP_get_wtime()	
	if(ptree%MyID==Main_ID)write(*,*) "FastMATVEC-based HODLR construction......"		
	call z_HODLR_randomized(ho_bf1,z_matvec_user,N_unk_loc,Memory,error,option,stats1,ker,ptree1,msh1)
	t2 = OMP_get_wtime()  
	if(ptree%MyID==Main_ID)write(*,*) "FastMATVEC-based HODLR construction finished",t2-t1, 'secnds. Error: ', error		
	
	
	
	if(allocated(quant%matZ_glo))deallocate(quant%matZ_glo)
	
	call z_delete_proctree(ptree)
	call z_delete_Hstat(stats)
	call z_delete_mesh(msh)
	call z_delete_kernelquant(ker)
	call z_delete_HOBF(ho_bf)
	
	
	call z_delete_proctree(ptree1)
	call z_delete_Hstat(stats1)
	call z_delete_mesh(msh1)	
	call z_delete_HOBF(ho_bf1)
	
	
    if(ptree%MyID==Main_ID)write(*,*) "-------------------------------program end-------------------------------------"
	
	call blacs_exit(1)
	call MPI_Finalize(ierr)
    ! ! ! ! pause

end PROGRAM MLMDA_DIRECT_SOLVER_3D_CFIE



