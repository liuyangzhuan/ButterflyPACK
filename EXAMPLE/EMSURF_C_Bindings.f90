module EMSURF_C_BINDINGS
   use, intrinsic :: iso_c_binding
   use EMSURF_MODULE
   use z_BPACK_DEFS, only: BPACK_eps0, BPACK_mu0, BPACK_pi
   implicit none

   type(quant_EMSURF), target, save :: quant_c
   logical, save :: quant_c_initialized = .false.

contains

   subroutine emsurf_initialize_c(data_dir_c, wavelength, frequency, use_frequency, &
      cfie_alpha, rcs_static, rcs_nsample, mesh_normal, scaling, mpi_fcomm) bind(C)
      character(kind=c_char), intent(in) :: data_dir_c(*)
      real(c_double), value :: wavelength, frequency, cfie_alpha, scaling
      integer(c_int), value :: use_frequency, rcs_static, rcs_nsample, mesh_normal
      integer(c_int), value :: mpi_fcomm
      character(len=:), allocatable :: data_dir
      integer :: ii, path_length
      if (quant_c_initialized) call emsurf_finalize_c()

      path_length = 0
      do while (data_dir_c(path_length + 1) /= c_null_char)
         path_length = path_length + 1
      end do
      allocate(character(len=path_length) :: data_dir)
      do ii = 1, path_length
         data_dir(ii:ii) = data_dir_c(ii)
      end do

      quant_c%integral_points = 6
      allocate (quant_c%ng1(quant_c%integral_points))
      allocate (quant_c%ng2(quant_c%integral_points))
      allocate (quant_c%ng3(quant_c%integral_points))
      allocate (quant_c%gauss_w(quant_c%integral_points))
      call gauss_points(quant_c)

      quant_c%mesh_normal = mesh_normal
      quant_c%scaling = scaling
      quant_c%RCS_static = rcs_static
      quant_c%RCS_Nsample = rcs_nsample
      quant_c%CFIE_alpha = cfie_alpha

      if (use_frequency /= 0) then
         quant_c%freq = frequency
         quant_c%wavelength = 1d0/frequency/dsqrt(BPACK_mu0*BPACK_eps0)
      else
         quant_c%wavelength = wavelength
         quant_c%freq = 1d0/wavelength/dsqrt(BPACK_mu0*BPACK_eps0)
      end if
      quant_c%wavenum = 2d0*BPACK_pi/quant_c%wavelength

      call geo_modeling_SURF(quant_c, mpi_fcomm, data_dir)
      quant_c_initialized = .true.
   end subroutine emsurf_initialize_c

   subroutine emsurf_get_problem_size_c(nunk) bind(C)
      integer(c_int), intent(out) :: nunk
      nunk = quant_c%Nunk
   end subroutine emsurf_get_problem_size_c

   subroutine emsurf_get_coordinates_c(coordinates) bind(C)
      real(c_double), intent(out) :: coordinates(*)
      integer :: edge, dim

      do edge = 1, quant_c%Nunk
         do dim = 1, 3
            coordinates(3*(edge - 1) + dim) = quant_c%xyz(dim, quant_c%maxnode + edge)
         end do
      end do
   end subroutine emsurf_get_coordinates_c

   subroutine emsurf_get_minedgelength_c(minedgelength) bind(C)
      real(c_double), intent(out) :: minedgelength
      minedgelength = quant_c%minedgelength
   end subroutine emsurf_get_minedgelength_c

   subroutine emsurf_get_wavenumber_c(wavenumber) bind(C)
      real(c_double), intent(out) :: wavenumber
      wavenumber = quant_c%wavenum
   end subroutine emsurf_get_wavenumber_c

   subroutine emsurf_entry_c(m, n, value, context) bind(C)
      integer(c_int), intent(in) :: m, n
      complex(c_double_complex), intent(out) :: value
      type(c_ptr), value :: context
      class(*), pointer :: quant_ptr

      quant_ptr => quant_c
      call Zelem_EMSURF(m, n, value, quant_ptr)
   end subroutine emsurf_entry_c

   subroutine emsurf_block_c(Ninter, Nallrows, Nallcols, Nalldat_loc, allrows, allcols, &
      alldat, rowidx, colidx, pgidx, Npmap, pmaps, context) bind(C)
      integer(c_int), intent(in) :: Ninter, Nallrows, Nallcols, Npmap
      integer(c_int64_t), intent(in) :: Nalldat_loc
      integer(c_int), intent(in), target :: allrows(*), allcols(*), rowidx(*), colidx(*), pgidx(*)
      integer(c_int), intent(in), target :: pmaps(*)
      complex(c_double_complex), intent(out), target :: alldat(*)
      type(c_ptr), value :: context
      integer(c_int), pointer :: pmaps_2d(:, :)
      class(*), pointer :: quant_ptr
      integer :: local_data_count

      if (Nalldat_loc > int(huge(local_data_count), c_int64_t)) then
         error stop 'EMSURF block callback exceeds the 32-bit global-array interface'
      end if
      local_data_count = int(Nalldat_loc)
      call c_f_pointer(c_loc(pmaps(1)), pmaps_2d, [Npmap, 3])
      quant_ptr => quant_c
      call Zelem_EMSURF_block(Ninter, allrows(1:Nallrows), allcols(1:Nallcols), &
         alldat(1:local_data_count), rowidx(1:Ninter), colidx(1:Ninter), &
         pgidx(1:Ninter), Npmap, pmaps_2d, quant_ptr)
   end subroutine emsurf_block_c

   subroutine emsurf_incident_c(edge, polarization, theta, phi, value) bind(C)
      integer(c_int), value :: edge, polarization
      real(c_double), value :: theta, phi
      complex(c_double_complex), intent(out) :: value

      if (polarization == 0) then
         call element_Vinc_VV_SURF(theta, phi, edge, value, quant_c)
      else
         call element_Vinc_HH_SURF(theta, phi, edge, value, quant_c)
      end if
   end subroutine emsurf_incident_c

   subroutine emsurf_rcs_contribution_c(edge, polarization, theta, phi, current, value) bind(C)
      integer(c_int), value :: edge, polarization
      real(c_double), value :: theta, phi
      complex(c_double_complex), value :: current
      complex(c_double_complex), intent(out) :: value

      if (polarization == 0) then
         call VV_polar_SURF(theta, phi, edge, value, current, quant_c)
      else
         call HH_polar_SURF(theta, phi, edge, value, current, quant_c)
      end if
   end subroutine emsurf_rcs_contribution_c

   subroutine emsurf_finalize_c() bind(C)
      integer :: node

      if (.not. quant_c_initialized) return
      call delete_quant_EMSURF(quant_c)
      if (allocated(quant_c%edge_of_patch)) deallocate (quant_c%edge_of_patch)
      if (allocated(quant_c%edge_of_node)) then
         do node = 1, size(quant_c%edge_of_node)
            if (allocated(quant_c%edge_of_node(node)%edges)) &
               deallocate (quant_c%edge_of_node(node)%edges)
         end do
         deallocate (quant_c%edge_of_node)
      end if
      quant_c_initialized = .false.
   end subroutine emsurf_finalize_c

end module EMSURF_C_BINDINGS
