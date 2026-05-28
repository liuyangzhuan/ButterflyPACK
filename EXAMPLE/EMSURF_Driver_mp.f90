! “ButterflyPACK” Copyright (c) 2018, The Regents of the University of California, through
! Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the
! U.S. Dept. of Energy). All rights reserved.

!> @file
!> @brief Mixed-precision 3D EFIE/MFIE/CFIE driver.
!> @details This driver constructs the double-complex forward matrix, copies the compressed
!> structure to single-complex precision for factorization, and solves with double-complex
!> outer iterations preconditioned by the single-complex inverse.

PROGRAM ButterflyPACK_IE_3D
   use c_BPACK_DEFS
   use EMSURF_MODULE, only: quant_EMSURF, delete_quant_EMSURF, gauss_points, &
      geo_modeling_SURF, Zelem_EMSURF, Zelem_EMSURF_block, element_Vinc_VV_SURF, &
      element_Vinc_HH_SURF, RCS_bistatic_SURF, RCS_monostatic_HH_SURF
   use c_BPACK_factor, only: c_BPACK_Factorization
   use c_BPACK_structure, only: c_BPACK_delete
   use c_BPACK_utilities
   use c_MISC_Utilities, only: c_CreatePtree, c_blacs_exit_wrp
   use m_BPACK_utilities, only: z2c_BPACK_copy, z2c_BPACK_Solution, z2c_mesh_copy, z2c_CopyOptions
   use z_BPACK_DEFS, only: z_Hoption, z_Hstat, z_mesh, z_Bmatrix, z_kernelquant, z_proctree
   use z_BPACK_structure, only: z_BPACK_delete
   use z_BPACK_constr, only: z_BPACK_construction_Init, z_BPACK_construction_Element
   use z_MISC_Utilities, only: z_CreatePtree
   use z_BPACK_utilities, only: z_InitStat, z_SetDefaultOptions, z_ReadOption, z_PrintOptions, &
      z_PrintStat, z_delete_proctree, z_delete_Hstat, z_delete_kernelquant, z_delete_mesh
#ifdef HAVE_OPENMP
   use omp_lib
#endif
   implicit none

   integer provided, ierr, nmpi, ii, nargs, flag
   integer v_major, v_minor, v_bugfix
   integer Nunk_loc_a, Nunk_loc_m
   integer, allocatable::groupmembers(:), Permutation_a(:)
   real(kind=8), allocatable::xyz_a(:, :)
   real(kind=8) t1, t2
   character(len=1024)::strings, strings1

   type(c_Hoption)::option_m
   type(c_Hstat)::stats_m
   type(c_mesh)::msh_m
   type(c_Bmatrix)::bmat_m
   type(c_proctree)::ptree_m

   type(z_Hoption)::option_a
   type(z_Hstat)::stats_a
   type(z_mesh)::msh_a
   type(z_Bmatrix)::bmat_a
   type(z_kernelquant)::ker_a
   type(quant_EMSURF), target::quant_a
   type(z_proctree)::ptree_a

   call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided, ierr)

   call MPI_Comm_size(MPI_Comm_World, nmpi, ierr)
   allocate(groupmembers(nmpi))
   do ii = 1, nmpi
      groupmembers(ii) = ii - 1
   enddo

   call c_CreatePtree(nmpi, groupmembers, MPI_Comm_World, ptree_m)
   call z_CreatePtree(nmpi, groupmembers, MPI_Comm_World, ptree_a)
   deallocate(groupmembers)

   if (ptree_m%MyID == Main_ID) then
      write(*,*) "-------------------------------Program Start----------------------------------"
      write(*,*) "ButterflyPACK_IE_3D_MP"
      call c_BPACK_GetVersionNumber(v_major, v_minor, v_bugfix)
      write(*,'(A23,I1,A1,I1,A1,I1,A1)') " ButterflyPACK Version:", v_major, ".", v_minor, ".", v_bugfix
      write(*,*) "   "
   endif

   call c_InitStat(stats_m)
   call z_InitStat(stats_a)
   call z_SetDefaultOptions(option_a)

   quant_a%integral_points = 6
   allocate(quant_a%ng1(quant_a%integral_points), quant_a%ng2(quant_a%integral_points), &
      quant_a%ng3(quant_a%integral_points), quant_a%gauss_w(quant_a%integral_points))
   call gauss_points(quant_a)

   quant_a%DATA_DIR = '../EXAMPLE/EM3D_DATA/sphere_2300'
   quant_a%mesh_normal = 1
   quant_a%scaling = 1d0
   quant_a%wavelength = 2.0d0
   quant_a%freq = 1/quant_a%wavelength/sqrt(BPACK_mu0*BPACK_eps0)
   quant_a%RCS_static = 2
   quant_a%RCS_Nsample = 1000
   quant_a%CFIE_alpha = 1.0d0

   option_a%ErrSol = 1
   option_a%format = HODLR
   option_a%near_para = 2.01d0
   option_a%verbosity = 1
   option_a%ILU = 0
   option_a%forwardN15flag = 0
   option_a%LRlevel = 100
   option_a%tol_itersol = 1d-5
   option_a%sample_para = 4d0
   option_a%knn = 50

   nargs = iargc()
   ii = 1
   do while (ii <= nargs)
      call getarg(ii, strings)
      if (trim(strings) == '-quant') then
         flag = 1
         do while (flag == 1)
            ii = ii + 1
            if (ii <= nargs) then
               call getarg(ii, strings)
               if (strings(1:2) == '--') then
                  ii = ii + 1
                  call getarg(ii, strings1)
                  if (trim(strings) == '--data_dir') then
                     quant_a%data_dir = trim(strings1)
                  else if (trim(strings) == '--wavelength') then
                     read(strings1,*) quant_a%wavelength
                     quant_a%freq = 1/quant_a%wavelength/sqrt(BPACK_mu0*BPACK_eps0)
                  else if (trim(strings) == '--freq') then
                     read(strings1,*) quant_a%freq
                     quant_a%wavelength = 1/quant_a%freq/sqrt(BPACK_mu0*BPACK_eps0)
                  else
                     if (ptree_m%MyID == Main_ID) write(*,*) 'ignoring unknown quant: ', trim(strings)
                  endif
               else
                  flag = 0
               endif
            else
               flag = 0
            endif
         enddo
      else if (trim(strings) == '-option') then
         call z_ReadOption(option_a, ptree_a, ii)
      else
         if (ptree_m%MyID == Main_ID) write(*,*) 'ignoring unknown argument: ', trim(strings)
         ii = ii + 1
      endif
   enddo

   quant_a%wavenum = 2*BPACK_pi/quant_a%wavelength

   if (ptree_m%MyID == Main_ID) then
      write(*,*) ''
      write(*,*) 'EFIE computing'
      write(*,*) 'frequency:', quant_a%freq
      write(*,*) 'wavelength:', quant_a%wavelength
      write(*,*) ''
   endif

   call geo_modeling_SURF(quant_a, ptree_a%Comm, quant_a%DATA_DIR)
   option_a%touch_para = 3*quant_a%minedgelength
   call z2c_CopyOptions(option_a, option_m)

   option_m%precon = DIRECT ! direct factorization
   option_a%precon = BPACKPRECON ! preconditioned iterative solver


   ker_a%QuantApp => quant_a
   ker_a%FuncZmn => Zelem_EMSURF
   ker_a%FuncZmnBlock => Zelem_EMSURF_block

   allocate(xyz_a(3, quant_a%Nunk))
   do ii = 1, quant_a%Nunk
      xyz_a(:, ii) = quant_a%xyz(:, quant_a%maxnode + ii)
   enddo

   t1 = MPI_Wtime()
   allocate(Permutation_a(quant_a%Nunk))
   if (ptree_m%MyID == Main_ID .and. option_m%verbosity >= 0) write(*,*) 'Constructing double-complex A'
   call z_PrintOptions(option_a, ptree_a)
   call z_BPACK_construction_Init(quant_a%Nunk, Permutation_a, Nunk_loc_a, bmat_a, option_a, &
      stats_a, msh_a, ker_a, ptree_a, Coordinates=xyz_a)
   deallocate(Permutation_a, xyz_a)

   call z2c_mesh_copy(msh_a, msh_m)
   Nunk_loc_m = msh_m%idxe - msh_m%idxs + 1
   if (Nunk_loc_a /= Nunk_loc_m) then
      if (ptree_m%MyID == Main_ID) write(*,*) 'double A and single M have different local ownership'
      stop
   endif
   t2 = MPI_Wtime()

   if (ptree_m%MyID == Main_ID .and. option_m%verbosity >= 0) write(*,*) 'Double-complex initialization:', t2 - t1, 'Seconds'

   if (ptree_m%MyID == Main_ID .and. option_m%verbosity >= 0) write(*,*) 'Single-complex M options'
   call c_PrintOptions(option_m, ptree_m)

   call z_BPACK_construction_Element(bmat_a, option_a, stats_a, msh_a, ker_a, ptree_a)

   if (ptree_m%MyID == Main_ID .and. option_m%verbosity >= 0) write(*,*) 'Copying double-complex A to single-complex M'
   call z2c_BPACK_copy(bmat_a, bmat_m, ptree_a)

   if (ptree_m%MyID == Main_ID .and. option_m%verbosity >= 0) write(*,*) 'Factorizing single-complex M'
   call c_BPACK_Factorization(bmat_m, option_m, stats_m, ptree_m, msh_m)

   call EM_solve_SURF(bmat_a, bmat_m, option_a, option_m, msh_a, quant_a, &
      ptree_a, ptree_m, stats_a, stats_m)

   if (ptree_m%MyID == Main_ID .and. option_m%verbosity >= 0) write(*,*) 'Double-complex A statistics'
   call z_PrintStat(stats_a, ptree_a)
   if (ptree_m%MyID == Main_ID .and. option_m%verbosity >= 0) write(*,*) 'Single-complex M statistics'
   call c_PrintStat(stats_m, ptree_m)

   if (ptree_m%MyID == Main_ID .and. option_m%verbosity >= 0) write(*,*) &
      "-------------------------------program end-------------------------------------"

   call delete_quant_EMSURF(quant_a)
   call z_delete_proctree(ptree_a)
   call c_delete_proctree(ptree_m)
   call z_delete_Hstat(stats_a)
   call c_delete_Hstat(stats_m)
   call z_delete_mesh(msh_a)
   call c_delete_mesh(msh_m)
   call z_delete_kernelquant(ker_a)
   call z_BPACK_delete(bmat_a)
   call c_BPACK_delete(bmat_m)

   call c_blacs_exit_wrp(1)
   call MPI_Finalize(ierr)

contains

   subroutine EM_solve_SURF(amat, mmat, option_a, option_m, msh_a, quant, ptree_a, ptree_m, stats_a, stats_m)
      implicit none
      type(z_Bmatrix), target::amat
      type(c_Bmatrix), target::mmat
      type(z_Hoption), target::option_a
      type(c_Hoption), target::option_m
      type(z_mesh)::msh_a
      type(quant_EMSURF)::quant
      type(z_proctree), target::ptree_a
      type(c_proctree), target::ptree_m
      type(z_Hstat), target::stats_a
      type(c_Hstat), target::stats_m

      integer edge, j, ierr
      integer num_sample, N_unk_loc
      real(kind=8) theta, phi, dphi, rcs_H
      real(kind=8) n1, n2, rtemp
      complex(kind=8) value_Z
      complex(kind=8), allocatable::current(:, :), voltage(:, :), x(:, :), b(:, :)

      if (ptree_m%MyID == Main_ID .and. option_m%verbosity >= 0) write(*,*) 'EM_solve......'

      N_unk_loc = msh_a%idxe - msh_a%idxs + 1

      if (quant%RCS_static == 2) then
         theta = 90d0
         phi = 0d0

         allocate(current(N_unk_loc, 2))
         current = 0d0
         allocate(voltage(N_unk_loc, 2))

#ifdef HAVE_OPENMP
!$omp parallel do default(shared) private(edge,value_Z)
#endif
         do edge = msh_a%idxs, msh_a%idxe
            call element_Vinc_VV_SURF(theta, phi, msh_a%new2old(edge), value_Z, quant)
            voltage(edge - msh_a%idxs + 1, 1) = value_Z
            call element_Vinc_HH_SURF(theta, phi, msh_a%new2old(edge), value_Z, quant)
            voltage(edge - msh_a%idxs + 1, 2) = value_Z
         enddo
#ifdef HAVE_OPENMP
!$omp end parallel do
#endif

         n1 = MPI_Wtime()
         call z2c_BPACK_Solution(amat, mmat, current, voltage, N_unk_loc, 2, option_a, ptree_a, &
            stats_a, option_m, ptree_m, stats_m)
         n2 = MPI_Wtime()
         stats_m%Time_Sol = stats_m%Time_Sol + n2 - n1

         if (ptree_m%MyID == Main_ID .and. option_m%verbosity >= 0) then
            write(*,*) ''
            write(*,*) 'Mixed-precision solving:', n2 - n1, 'Seconds'
            write(*,*) ''
         endif

         n1 = MPI_Wtime()
         call RCS_bistatic_SURF(current, msh_a, quant, ptree_a)
         n2 = MPI_Wtime()

         if (ptree_m%MyID == Main_ID .and. option_m%verbosity >= 0) then
            write(*,*) ''
            write(*,*) 'Bistatic RCS', n2 - n1, 'Seconds'
            write(*,*) ''
         endif
         deallocate(current, voltage)

      elseif (quant%RCS_static == 1) then
         allocate(current(N_unk_loc, 1))

         num_sample = quant%RCS_Nsample
         theta = 90d0
         dphi = 180d0/num_sample
         allocate(b(N_unk_loc, num_sample + 1))
         allocate(x(N_unk_loc, num_sample + 1))
         x = 0d0

         if (ptree_m%MyID == Main_ID) open(100, file='monostaticH.out')

         n1 = MPI_Wtime()
         do j = 0, num_sample
            phi = j*dphi
#ifdef HAVE_OPENMP
!$omp parallel do default(shared) private(edge,value_Z)
#endif
            do edge = msh_a%idxs, msh_a%idxe
               call element_Vinc_HH_SURF(theta, phi, msh_a%new2old(edge), value_Z, quant)
               b(edge - msh_a%idxs + 1, j + 1) = value_Z
            enddo
#ifdef HAVE_OPENMP
!$omp end parallel do
#endif
         enddo

         call z2c_BPACK_Solution(amat, mmat, x, b, N_unk_loc, num_sample + 1, option_a, ptree_a, &
            stats_a, option_m, ptree_m, stats_m)

         do j = 0, num_sample
            phi = j*dphi
            current(:, 1) = x(:, j + 1)
            call RCS_monostatic_HH_SURF(theta, phi, rcs_H, current(:, 1), msh_a, quant, ptree_a)
            if (ptree_m%MyID == Main_ID) write(100,*) phi, rcs_H
         enddo

         n2 = MPI_Wtime()
         stats_m%Time_Sol = stats_m%Time_Sol + n2 - n1
         call MPI_ALLREDUCE(stats_m%Time_Sol, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree_m%Comm, ierr)

         if (ptree_m%MyID == Main_ID) then
            close(100)
            write(*,*) ''
            write(*,*) 'Mixed-precision solving:', rtemp, 'Seconds'
            write(*,*) ''
         endif

         call MPI_ALLREDUCE(stats_m%Flop_Sol, rtemp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree_m%Comm, ierr)
         if (ptree_m%MyID == Main_ID .and. option_m%verbosity >= 0) write(*,'(A13,Es14.2)') 'Solve flops:', rtemp

         deallocate(current, b, x)
      endif

      if (ptree_m%MyID == Main_ID .and. option_m%verbosity >= 0) then
         write(*,*) 'EM_solve finished'
         write(*,*) '    '
      endif
   end subroutine EM_solve_SURF

end PROGRAM ButterflyPACK_IE_3D
