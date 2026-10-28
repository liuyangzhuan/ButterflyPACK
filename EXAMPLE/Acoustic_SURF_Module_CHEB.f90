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
!> @file
!> @brief This file contains functions and data types for the 3D sound-soft acoustic scattering example discretized by the Chebyshev-based rectangular-polar Nystrom method
!> @details The discretization follows O.P. Bruno and E. Garza, "A Chebyshev-based rectangular-polar integral solver for scattering by geometries described by non-overlapping patches", J. Comput. Phys. 421 (2020) 109740. \n
!> The surface is a union of non-overlapping logically-quadrilateral patches r(s,t), (s,t) in [-1,1]^2. Each patch carries N x N Fejer (first-kind Chebyshev) nodes in (u,v), mapped to (s,t) by the edge change of variables s=eta_s(u), t=eta_t(v). \n
!> Closed surfaces solve 1/2 phi + D[phi] - ik S[phi] = -u_inc, open surfaces solve S[phi] = -u_inc, with G(x,y)=exp(ik|x-y|)/(4 pi |x-y|). \n
!> Far interactions use Fejer's first rule. Self-patch and near-singular interactions use rectangular-polar product-integration weights against the Chebyshev interpolant of the density, precomputed once and replicated on all MPI ranks.

! This exmple works with double-complex precision data
module Acoustic_SURF_MODULE_CHEB
use z_BPACK_DEFS
use z_MISC_Utilities
#ifdef HAVE_OPENMP
use omp_lib
#endif
implicit none

	!**** analytic geometries
	integer, parameter:: CHEB_SPHERE = 1, CHEB_CUBE = 2, CHEB_DISK = 3, CHEB_PLATE = 4
	!**** patch parametrizations
	integer, parameter:: PATCH_AFFINE = 1, PATCH_CUBEDSPHERE = 2, PATCH_DISKSIDE = 3
	!**** unknowns: physical density phi, or edge-resolved density psi = eta_s'(u)*eta_t'(v)*phi
	integer, parameter:: UNK_PHI = 0, UNK_PSI = 1

	!**** per-thread cache of one weight row for single-entry requests in the on-the-fly mode
	integer, save:: cache_m = -1, cache_kk = -1
	complex(kind=8), allocatable, save:: cache_W(:)
#ifdef HAVE_OPENMP
	!$omp threadprivate(cache_m, cache_kk, cache_W)
#endif

	!**** one logically-quadrilateral patch r(s,t), (s,t) in [-1,1]^2
	type patch_CHEB
		integer:: kind = PATCH_AFFINE
		integer:: edge(4) = 0 ! 1 if the side s=-1, s=+1, t=-1, t=+1 (in this order) is a geometric edge
		logical:: planar = .true.
		real(kind=8):: c(3) = 0d0 ! affine: r=c+s*a1+t*a2; cubed sphere: face axis; disk side: origin
		real(kind=8):: a1(3) = 0d0, a2(3) = 0d0 ! affine: tangent vectors; cubed sphere and disk side: local axes
		real(kind=8):: s0 = 0d0, hs = 1d0, t0 = 0d0, ht = 1d0 ! sub-patch of the parent face: (s0+s*hs, t0+t*ht)
		real(kind=8):: rad = 1d0 ! sphere radius or disk radius
		real(kind=8):: rin = 0.5d0 ! disk: half width of the central square
		real(kind=8):: bcenter(3) = 0d0, bradius = 0d0 ! bounding sphere
	end type patch_CHEB

	!**** quantities related to geometries, discretization, near-field weights and tests
	type quant_ACOUSTIC_CHEB
		real(kind=8):: wavenum = 2d0*BPACK_pi ! wave number
		real(kind=8):: wavelength = 1d0 ! wave length
		integer:: geo = CHEB_SPHERE ! one of CHEB_SPHERE, CHEB_CUBE, CHEB_DISK, CHEB_PLATE
		real(kind=8):: radius = 1d0 ! sphere radius, cube/plate half side, or disk radius
		integer:: split = 1 ! each face is split into split x split patches
		logical:: closed = .true. ! .true.: combined-field equation, .false.: single-layer equation
		integer:: unknown = -1 ! UNK_PHI, UNK_PSI, or -1: phi for closed and psi for open surfaces
		real(kind=8):: inc_theta = 0d0, inc_phi = 0d0 ! incident direction in degrees
		real(kind=8):: dinc(3) = (/0d0, 0d0, 1d0/) ! incident direction

		integer:: N = 12 ! Fejer nodes per patch per dimension
		integer:: Nbeta = 100 ! rectangular-polar quadrature nodes per dimension; geometries with edges need ~200 for errors below 1e-5
		real(kind=8):: pxi = 6d0 ! exponent p of the rectangular-polar change of variables xi
		real(kind=8):: pedge = 2d0 ! exponent p of the edge change of variables eta
		real(kind=8):: delta = -1d0 ! proximity distance; <0: half the typical patch size
		integer:: naive_chord = 0 ! 1: form self-patch chords by subtracting coordinates (for testing only)
		integer:: forest = 0 ! 1: order the patches by the patch hierarchy (each base patch split 2x2 recursively, Morton order) and pass it to ButterflyPACK as a forest of cluster trees (option%ntree=nbase)
		integer:: myid = 0 ! MPI rank of this process
		integer:: nbase = 1 ! number of base patches (sphere/cube 6, disk 5, plate 1), i.e. the roots of the forest
		integer:: onthefly = 0 ! 1: compute each self/near-singular weight row when its matrix entries are requested (no replicated weights); 0: precompute and replicate all of them

		integer:: npatch = 0
		type(patch_CHEB), allocatable:: patches(:)
		real(kind=8), allocatable:: u(:), wu(:) ! Fejer nodes and weights (N)
		real(kind=8), allocatable:: tb(:), wb(:) ! Fejer nodes and weights (Nbeta)
		real(kind=8), allocatable:: Ccheb(:, :) ! Ccheb(n,i) = c_n T_n(u_i), maps node values to Chebyshev coefficients

		integer:: Nunk = 0 ! size of the matrix
		real(kind=8), allocatable:: xyz(:, :) ! coordinates of the nodes
		real(kind=8), allocatable:: nrm(:, :) ! unit normal at the nodes
		real(kind=8), allocatable:: st(:, :) ! patch coordinates (s,t) of the nodes
		real(kind=8), allocatable:: Jw(:) ! surface Jacobian times Fejer weights, J(s,t)*w_i*w_j
		real(kind=8), allocatable:: e(:) ! edge factors eta_s'(u_i)*eta_t'(v_j)

		integer, allocatable:: near_ptr(:) ! near pairs of target m are near_ptr(m):near_ptr(m+1)-1
		integer, allocatable:: near_q(:) ! patch of each near pair
		complex(kind=8), allocatable:: near_W(:) ! N*N weights of each near pair (precomputed mode)
		real(kind=8), allocatable:: near_ub(:), near_vb(:) ! (u,v) of the target's closest point on the patch of each near pair (on-the-fly mode)
		integer:: npair = 0
		integer(kind=8):: nweights_computed = 0 ! number of weight rows computed on the fly (diagnostic)
		integer(kind=8):: nentries = 0, nsubmat = 0 ! number of matrix entries and submatrices requested through Zsubmatrix_CHEB (diagnostic)
		integer(kind=8):: nweights_cat(3) = 0 ! ... requested by: 1 whole leaf-patch blocks, 2 other submatrices, 3 single entries
		integer, allocatable:: paircount(:) ! diagnostic: number of on-the-fly computations of each near pair within whole leaf-patch blocks

		integer:: eigtest = 5 ! degree of the spherical-harmonic forward-map test (sphere only), 0: skip
		integer:: nfar = 181 ! number of far-field directions
		integer:: gmres_restart = 0 ! >0: GMRES restart length when precon/=DIRECT (0: ButterflyPACK's default solve path, GMRES(30))
		integer:: faceorder = 0 ! order of the six faces of the sphere/cube (the forest's base patches): 0: +x,-x,+y,-y,+z,-z; 1: sweep +z,+x,+y,-x,-y,-z
		integer:: printstruct = 0 ! 1: write every block's rank before (fort.100000+rank) and after (fort.200000+rank) the factorization, and the cluster geometry (clusters.txt)
		integer:: tensor = 0 ! 1: iterative solve with a forest of H-TENSOR blocks, one per pair of base patches on their (split*N) x (split*N) parameter grids (split must be a power of two)
		integer:: tlev = 0 ! log2(split): levels of the per-dimension trees of a base patch down to single patches
		real(kind=8):: tnear_para = 2.2d0 ! admissibility parameter of the tensor forest (option%near_para)
		integer, allocatable:: sp_face(:), sp_is(:), sp_js(:) ! base patch and split position (is,js) of each patch
		integer, allocatable:: nsub_ptr(:), nsub_q(:) ! patches holding self/near pairs of the targets of each patch (CSR)
		real(kind=8), allocatable:: tnode_c(:, :, :), tnode_r(:, :) ! bounding sphere (center, radius) of every quadtree node of every base patch (tnode_id)
	end type quant_ACOUSTIC_CHEB

	!**** one block of the tensor forest: Z(base patch p, base patch q) as a 4-mode tensor on the two parameter grids
	type quant_pair_CHEB
		integer:: p = 0, q = 0
		type(quant_ACOUSTIC_CHEB), pointer:: g => null()
	end type quant_pair_CHEB

	!**** the tensor forest: nbase x nbase H-TENSOR blocks sharing the parameter grids, the per-dimension cluster trees and
	! the vector layout. A local vector holds the local box of every base patch, base patch after base patch.
	type tforest_CHEB
		integer:: nface = 0
		integer:: Ns(2) = 0, Ns_loc(2) = 0, nloc = 0 ! grid points per dimension, local box, local entries per base patch
		type(z_Bmatrix), allocatable:: bmat(:, :)
		type(z_mesh), allocatable:: msh(:, :, :) ! msh(dim, p, q)
		type(z_kernelquant), allocatable:: ker(:, :)
		type(z_Hstat), allocatable:: stats(:, :)
		type(quant_pair_CHEB), allocatable:: pair(:, :)
		integer, allocatable:: gidx(:, :) ! gidx(i, p): node (natural order) of local entry i of base patch p
		type(z_Hoption), pointer:: option => null()
		type(z_proctree), pointer:: ptree => null()
		real(kind=8):: tmvp = 0d0 ! time of the products
		integer:: nmvp = 0 ! number of products (vectors)
	end type tforest_CHEB

contains

	subroutine delete_quant_Acoustic_CHEB(quant)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		if (allocated(quant%patches)) deallocate (quant%patches)
		if (allocated(quant%u)) deallocate (quant%u)
		if (allocated(quant%wu)) deallocate (quant%wu)
		if (allocated(quant%tb)) deallocate (quant%tb)
		if (allocated(quant%wb)) deallocate (quant%wb)
		if (allocated(quant%Ccheb)) deallocate (quant%Ccheb)
		if (allocated(quant%xyz)) deallocate (quant%xyz)
		if (allocated(quant%nrm)) deallocate (quant%nrm)
		if (allocated(quant%st)) deallocate (quant%st)
		if (allocated(quant%Jw)) deallocate (quant%Jw)
		if (allocated(quant%e)) deallocate (quant%e)
		if (allocated(quant%near_ptr)) deallocate (quant%near_ptr)
		if (allocated(quant%near_q)) deallocate (quant%near_q)
		if (allocated(quant%near_W)) deallocate (quant%near_W)
		if (allocated(quant%near_ub)) deallocate (quant%near_ub)
		if (allocated(quant%near_vb)) deallocate (quant%near_vb)
		if (allocated(quant%paircount)) deallocate (quant%paircount)
		if (allocated(quant%sp_face)) deallocate (quant%sp_face)
		if (allocated(quant%sp_is)) deallocate (quant%sp_is)
		if (allocated(quant%sp_js)) deallocate (quant%sp_js)
		if (allocated(quant%nsub_ptr)) deallocate (quant%nsub_ptr)
		if (allocated(quant%nsub_q)) deallocate (quant%nsub_q)
		if (allocated(quant%tnode_c)) deallocate (quant%tnode_c)
		if (allocated(quant%tnode_r)) deallocate (quant%tnode_r)
	end subroutine delete_quant_Acoustic_CHEB

	!**** user-defined subroutine to sample Z_mn
	subroutine Zelem_Acoustic_CHEB(m, n, value, quant)
		implicit none
		integer, INTENT(IN):: m, n
		complex(kind=8) value
		class(*), pointer :: quant

		select TYPE (quant)
		type is (quant_ACOUSTIC_CHEB)
			value = Zentry_CHEB(quant, m, n)
		class default
			write (*, *) "unexpected type"
			stop
		end select
	end subroutine Zelem_Acoustic_CHEB

	!**** matrix entry Z_mn: row m is the target node, column n the source node (both in natural order)
	complex(kind=8) function Zentry_CHEB(quant, m, n)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		integer m, n, q, ij, N2, kk
		complex(kind=8) base

		N2 = quant%N*quant%N
		q = (n - 1)/N2 + 1
		ij = n - (q - 1)*N2
		kk = near_pair_CHEB(quant, m, q)
		if (kk > 0) then ! self or near-singular: rectangular-polar weights
			if (quant%onthefly == 1) then ! a one-row cache per thread serves consecutive columns of the same patch
				if (cache_m /= m .or. cache_kk /= kk) then
					if (.not. allocated(cache_W)) allocate (cache_W(N2))
					if (size(cache_W) /= N2) then
						deallocate (cache_W)
						allocate (cache_W(N2))
					endif
					call near_weights_row_CHEB(quant, m, kk, cache_W, 3)
					cache_m = m
					cache_kk = kk
				endif
				base = cache_W(ij)
			else
				base = quant%near_W(int(kk - 1, 8)*N2 + ij)
			endif
		else ! far: Fejer's rule
			base = Zfar_CHEB(quant, m, n)
		endif
		Zentry_CHEB = Zscale_CHEB(quant, m, n, base)
	end function Zentry_CHEB

	!**** index of the near pair (m, patch q) in near_q, or 0 if patch q is far from target m
	integer function near_pair_CHEB(quant, m, q)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		integer m, q, kk
		near_pair_CHEB = 0
		do kk = quant%near_ptr(m), quant%near_ptr(m + 1) - 1
			if (quant%near_q(kk) == q) then
				near_pair_CHEB = kk
				exit
			endif
		enddo
	end function near_pair_CHEB

	!**** far interaction: kernel times the Fejer weight of source node n
	complex(kind=8) function Zfar_CHEB(quant, m, n)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		integer m, n
		real(kind=8) d(3), R, dn
		d = quant%xyz(:, m) - quant%xyz(:, n)
		R = sqrt(sum(d**2))
		dn = dot_product(d, quant%nrm(:, n))
		Zfar_CHEB = kernel_H(quant%wavenum, R, dn, quant%closed)*quant%Jw(n)
	end function Zfar_CHEB

	!**** apply the edge factor of the chosen unknown and the 1/2 identity term (closed surfaces)
	complex(kind=8) function Zscale_CHEB(quant, m, n, base)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		integer m, n
		complex(kind=8) base, value
		value = base
		if (quant%unknown == UNK_PHI) value = value*quant%e(n)
		if (quant%closed .and. m == n) then
			if (quant%unknown == UNK_PHI) then
				value = value + 0.5d0
			else
				value = value + 0.5d0/quant%e(n)
			endif
		endif
		Zscale_CHEB = value
	end function Zscale_CHEB

	!**** on-the-fly mode: the N*N rectangular-polar weights of near pair kk (target m)
	subroutine near_weights_row_CHEB(quant, m, kk, W, cat)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		integer m, kk, q, own, cat
		complex(kind=8) W(quant%N*quant%N)
		q = quant%near_q(kk)
		own = merge(1, 0, q == (m - 1)/(quant%N*quant%N) + 1)
		call pair_weights_CHEB(quant, m, q, quant%near_ub(kk), quant%near_vb(kk), own, W)
#ifdef HAVE_OPENMP
		!$omp atomic
#endif
		quant%nweights_computed = quant%nweights_computed + 1
#ifdef HAVE_OPENMP
		!$omp atomic
#endif
		quant%nweights_cat(cat) = quant%nweights_cat(cat) + 1
		if (cat == 1 .and. allocated(quant%paircount)) then
#ifdef HAVE_OPENMP
			!$omp atomic
#endif
			quant%paircount(kk) = quant%paircount(kk) + 1
		endif
	end subroutine near_weights_row_CHEB

	!**** user-defined subroutine to sample a list of submatrices (ButterflyPACK's FuncZmnBlock interface, used with option%elem_extract>=1).
	! allrows, allcols hold the (natural-order) row and column indices of all Ninter submatrices; alldat_loc receives the owned submatrices, each column major.
	subroutine Zelem_Acoustic_CHEB_block(Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps, quant)
		implicit none
		class(*), pointer :: quant
		integer:: Ninter
		integer:: allrows(:), allcols(:)
		complex(kind=8), target:: alldat_loc(:)
		integer:: colidx(Ninter), rowidx(Ninter), pgidx(Ninter)
		integer:: Npmap, pmaps(Npmap, 3)
		integer nn, nn1, nloc, pp, nr, nc, nthreads
		integer(kind=8) idx_row, idx_col, idx_val
		integer(kind=8), allocatable:: roff(:), coff(:), voff(:)
		integer, allocatable:: inter_map(:)

		select TYPE (quant)
		type is (quant_ACOUSTIC_CHEB)
			allocate (roff(Ninter), coff(Ninter), voff(Ninter), inter_map(Ninter))
			nloc = 0
			idx_row = 0
			idx_col = 0
			idx_val = 0
			do nn = 1, Ninter
				pp = pgidx(nn)
				if (pmaps(pp, 1)*pmaps(pp, 2) /= 1) then
					write (*, *) 'Zelem_Acoustic_CHEB_block: submatrices shared by several processes are not supported'
					stop
				endif
				if (pmaps(pp, 3) == quant%myid) then
					nloc = nloc + 1
					inter_map(nloc) = nn
					roff(nloc) = idx_row
					coff(nloc) = idx_col
					voff(nloc) = idx_val
					idx_val = idx_val + int(rowidx(nn), 8)*colidx(nn)
				endif
				idx_row = idx_row + rowidx(nn)
				idx_col = idx_col + colidx(nn)
			enddo
			nthreads = 1
#ifdef HAVE_OPENMP
			nthreads = omp_get_max_threads()
#endif
			if (nloc >= nthreads) then ! many submatrices: one thread per submatrix
#ifdef HAVE_OPENMP
			!$omp parallel do default(shared) private(nn1, nr, nc) schedule(dynamic,1)
#endif
			do nn1 = 1, nloc
				nr = rowidx(inter_map(nn1))
				nc = colidx(inter_map(nn1))
				call Zsubmatrix_CHEB(quant, nr, nc, allrows(roff(nn1) + 1:roff(nn1) + nr), allcols(coff(nn1) + 1:coff(nn1) + nc), &
					&	alldat_loc(voff(nn1) + 1:voff(nn1) + int(nr, 8)*nc), 0)
			enddo
#ifdef HAVE_OPENMP
			!$omp end parallel do
#endif
			else ! few submatrices (ButterflyPACK's H construction passes one block at a time): threads share the rows of each
			do nn1 = 1, nloc
				nr = rowidx(inter_map(nn1))
				nc = colidx(inter_map(nn1))
				call Zsubmatrix_CHEB(quant, nr, nc, allrows(roff(nn1) + 1:roff(nn1) + nr), allcols(coff(nn1) + 1:coff(nn1) + nc), &
					&	alldat_loc(voff(nn1) + 1:voff(nn1) + int(nr, 8)*nc), 1)
			enddo
			endif
			deallocate (roff, coff, voff, inter_map)
		class default
			write (*, *) "unexpected type"
			stop
		end select
	end subroutine Zelem_Acoustic_CHEB_block

	!**** one submatrix Z(rows, cols): each needed weight row (target, near patch) is computed once
	subroutine Zsubmatrix_CHEB(quant, nr, nc, rows, cols, dat, par)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		integer nr, nc, rows(nr), cols(nc), par
		complex(kind=8) dat(nr, nc)
		integer i, j, m, n, q, k, kk, k0, nk, N2, cat
		complex(kind=8) base
		complex(kind=8), allocatable:: Wrows(:, :)
		integer, allocatable:: have(:)

		N2 = quant%N*quant%N
#ifdef HAVE_OPENMP
		!$omp atomic
#endif
		quant%nentries = quant%nentries + int(nr, 8)*nc
#ifdef HAVE_OPENMP
		!$omp atomic
#endif
		quant%nsubmat = quant%nsubmat + 1
		cat = 2 ! a whole leaf-patch block: all N*N rows of one patch against all N*N columns of one patch
		if (nr == N2 .and. nc == N2) then
			if (minval(rows) == maxval(rows) - N2 + 1 .and. mod(minval(rows) - 1, N2) == 0 .and. &
				&	minval(cols) == maxval(cols) - N2 + 1 .and. mod(minval(cols) - 1, N2) == 0) cat = 1
		endif
#ifdef HAVE_OPENMP
		!$omp parallel if(par == 1) default(shared) private(i, j, m, n, q, k, kk, k0, nk, base, Wrows, have) firstprivate(cat)
#endif
		allocate (Wrows(merge(N2, 1, quant%onthefly == 1), 8), have(8))
#ifdef HAVE_OPENMP
		!$omp do schedule(dynamic,1)
#endif
		do i = 1, nr
			m = rows(i)
			k0 = quant%near_ptr(m)
			nk = quant%near_ptr(m + 1) - k0
			if (nk > size(have)) then
				deallocate (Wrows, have)
				allocate (Wrows(merge(N2, 1, quant%onthefly == 1), nk), have(nk))
			endif
			have(1:nk) = 0
			do j = 1, nc
				n = cols(j)
				q = (n - 1)/N2 + 1
				kk = 0
				do k = 1, nk
					if (quant%near_q(k0 + k - 1) == q) then
						kk = k
						exit
					endif
				enddo
				if (kk > 0) then
					if (quant%onthefly == 1) then
						if (have(kk) == 0) then
							call near_weights_row_CHEB(quant, m, k0 + kk - 1, Wrows(:, kk), cat)
							have(kk) = 1
						endif
						base = Wrows(n - (q - 1)*N2, kk)
					else
						base = quant%near_W(int(k0 + kk - 2, 8)*N2 + n - (q - 1)*N2)
					endif
				else
					base = Zfar_CHEB(quant, m, n)
				endif
				dat(i, j) = Zscale_CHEB(quant, m, n, base)
			enddo
		enddo
#ifdef HAVE_OPENMP
		!$omp end do
#endif
		deallocate (Wrows, have)
#ifdef HAVE_OPENMP
		!$omp end parallel
#endif
	end subroutine Zsubmatrix_CHEB

	!**** kernel H(x,y): dG/dn_y - ik G (closed) or G (open); dn = (x-y).n_y
	complex(kind=8) function kernel_H(k, R, dn, closed)
		implicit none
		real(kind=8) k, R, dn
		logical closed
		complex(kind=8) G
		if (R < 1d-14) then
			kernel_H = 0d0
			return
		endif
		G = exp(BPACK_junit*k*R)/(4d0*BPACK_pi*R)
		if (closed) then
			kernel_H = G*(1d0 - BPACK_junit*k*R)*dn/R**2 - BPACK_junit*k*G
		else
			kernel_H = G
		endif
	end function kernel_H

	!**** Fejer's first quadrature rule on [-1,1], eqs. (28)-(29)
	subroutine fejer_rule(N, x, w)
		implicit none
		integer N, j, l
		real(kind=8) x(N), w(N)
		do j = 1, N
			x(j) = cos(BPACK_pi*(2*j - 1)/(2d0*N))
			w(j) = 1d0
			do l = 1, N/2
				w(j) = w(j) - 2d0*cos(l*BPACK_pi*(2*j - 1)/dble(N))/(4d0*l*l - 1d0)
			enddo
			w(j) = 2d0*w(j)/N
		enddo
	end subroutine fejer_rule

	!**** the change of variables w(tau), eqs. (19)-(20), and its derivative
	real(kind=8) function vfun(t, p)
		implicit none
		real(kind=8) t, p
		vfun = max(0d0, (1d0/p - 0.5d0)*((BPACK_pi - t)/BPACK_pi)**3 + (1d0/p)*((t - BPACK_pi)/BPACK_pi) + 0.5d0)
	end function vfun

	real(kind=8) function dvfun(t, p)
		implicit none
		real(kind=8) t, p
		dvfun = -3d0*(1d0/p - 0.5d0)*(BPACK_pi - t)**2/BPACK_pi**3 + 1d0/(p*BPACK_pi)
	end function dvfun

	real(kind=8) function wmap(t, p)
		implicit none
		real(kind=8) t, p, tt, a, b
		tt = min(2d0*BPACK_pi, max(0d0, t))
		a = vfun(tt, p)**p
		b = vfun(2d0*BPACK_pi - tt, p)**p
		wmap = 2d0*BPACK_pi*a/(a + b)
	end function wmap

	real(kind=8) function dwmap(t, p)
		implicit none
		real(kind=8) t, p, tt, va, vb, a, b, da, db
		tt = min(2d0*BPACK_pi, max(0d0, t))
		va = vfun(tt, p)
		vb = vfun(2d0*BPACK_pi - tt, p)
		a = va**p
		b = vb**p
		da = p*va**(p - 1d0)*dvfun(tt, p)
		db = -p*vb**(p - 1d0)*dvfun(2d0*BPACK_pi - tt, p)
		dwmap = 2d0*BPACK_pi*(da*b - a*db)/(a + b)**2
	end function dwmap

	!**** edge change of variables s=eta(u), eqs. (21)-(22); el/er flag edges at s=-1/s=+1
	subroutine eta_map(u, el, er, p, s, ds)
		implicit none
		real(kind=8) u, p, s, ds
		integer el, er
		if (el == 1 .and. er == 1) then
			s = -1d0 + wmap(BPACK_pi*(u + 1d0), p)/BPACK_pi
			ds = dwmap(BPACK_pi*(u + 1d0), p)
		elseif (el == 1) then
			s = -1d0 + 2d0/BPACK_pi*wmap(BPACK_pi/2d0*(u + 1d0), p)
			ds = dwmap(BPACK_pi/2d0*(u + 1d0), p)
		elseif (er == 1) then
			s = -3d0 + 2d0/BPACK_pi*wmap(BPACK_pi + BPACK_pi/2d0*(u + 1d0), p)
			ds = dwmap(BPACK_pi + BPACK_pi/2d0*(u + 1d0), p)
		else
			s = u
			ds = 1d0
		endif
		s = min(1d0, max(-1d0, s))
	end subroutine eta_map

	!**** inverse of eta_map by bisection (eta is monotone)
	real(kind=8) function eta_inv(s, el, er, p)
		implicit none
		real(kind=8) s, p, lo, hi, mid, sm, dummy
		integer el, er, it
		if (el == 0 .and. er == 0) then
			eta_inv = s
			return
		endif
		if (s <= -1d0) then
			eta_inv = -1d0
			return
		endif
		if (s >= 1d0) then
			eta_inv = 1d0
			return
		endif
		lo = -1d0
		hi = 1d0
		do it = 1, 200
			mid = 0.5d0*(lo + hi)
			if (mid <= lo .or. mid >= hi) exit
			call eta_map(mid, el, er, p, sm, dummy)
			if (sm < s) then
				lo = mid
			else
				hi = mid
			endif
		enddo
		eta_inv = 0.5d0*(lo + hi)
	end function eta_inv

	!**** rectangular-polar change of variables xi_alpha(tau), eq. (40); returns xi-alpha (formed without subtraction) and xi'
	subroutine xi_map(tau, alpha, p, off, dxi)
		implicit none
		real(kind=8) tau, alpha, p, off, dxi, sg
		if (alpha >= 1d0) then
			off = -2d0/BPACK_pi*wmap(BPACK_pi*(1d0 - tau)/2d0, p)
			dxi = dwmap(BPACK_pi*(1d0 - tau)/2d0, p)
		elseif (alpha <= -1d0) then
			off = 2d0/BPACK_pi*wmap(BPACK_pi*(tau + 1d0)/2d0, p)
			dxi = dwmap(BPACK_pi*(tau + 1d0)/2d0, p)
		else
			sg = sign(1d0, tau)
			off = (sg - alpha)/BPACK_pi*wmap(BPACK_pi*abs(tau), p)
			dxi = (1d0 - alpha*sg)*dwmap(BPACK_pi*abs(tau), p)
		endif
	end subroutine xi_map

	!**** L(a,i) = l_i(x_a), the Lagrange basis on the first-kind Chebyshev nodes (via discrete orthogonality)
	subroutine cheb_lagrange(N, Ccheb, M, x, L)
		implicit none
		integer N, M, a, nn
		real(kind=8) Ccheb(0:N - 1, N), x(M), L(M, N)
		real(kind=8) T(M, 0:N - 1)
		do a = 1, M
			T(a, 0) = 1d0
			if (N > 1) T(a, 1) = x(a)
			do nn = 1, N - 2
				T(a, nn + 1) = 2d0*x(a)*T(a, nn) - T(a, nn - 1)
			enddo
		enddo
		L = matmul(T, Ccheb)
	end subroutine cheb_lagrange

	!**** position r(s,t) and tangent vectors r_s, r_t of a patch
	subroutine patch_eval(P, s, t, r, rs, rt)
		implicit none
		type(patch_CHEB):: P
		real(kind=8) s, t, r(3), rs(3), rt(3)
		real(kind=8) f1, f2, PP(3), Ps(3), Pt(3), np, sig, tau, lam, cth, sth, sq(2), cr(2), q(2), qs(2), qt(2)

		select case (P%kind)
		case (PATCH_AFFINE)
			r = P%c + s*P%a1 + t*P%a2
			rs = P%a1
			rt = P%a2
		case (PATCH_CUBEDSPHERE) ! equiangular cubed sphere: r = rad*P/|P|, P = c + tan(f1)*a1 + tan(f2)*a2
			f1 = BPACK_pi/4d0*(P%s0 + s*P%hs)
			f2 = BPACK_pi/4d0*(P%t0 + t*P%ht)
			PP = P%c + tan(f1)*P%a1 + tan(f2)*P%a2
			Ps = BPACK_pi/4d0*P%hs/cos(f1)**2*P%a1
			Pt = BPACK_pi/4d0*P%ht/cos(f2)**2*P%a2
			np = sqrt(sum(PP**2))
			r = P%rad*PP/np
			rs = P%rad*(Ps/np - PP*dot_product(PP, Ps)/np**3)
			rt = P%rad*(Pt/np - PP*dot_product(PP, Pt)/np**3)
		case (PATCH_DISKSIDE) ! blend between a side of the central square and a quarter of the circle
			sig = P%s0 + s*P%hs
			tau = P%t0 + t*P%ht
			lam = (1d0 + sig)/2d0
			cth = cos(BPACK_pi/4d0*tau)
			sth = sin(BPACK_pi/4d0*tau)
			sq = (/P%rin, P%rin*tau/)
			cr = (/P%rad*cth, P%rad*sth/)
			q = (1d0 - lam)*sq + lam*cr
			qs = P%hs/2d0*(cr - sq)
			qt = P%ht*((1d0 - lam)*(/0d0, P%rin/) + lam*P%rad*BPACK_pi/4d0*(/-sth, cth/))
			r = P%c + q(1)*P%a1 + q(2)*P%a2
			rs = qs(1)*P%a1 + qs(2)*P%a2
			rt = qt(1)*P%a1 + qt(2)*P%a2
		end select
	end subroutine patch_eval

	!**** position, unit normal and surface Jacobian of a patch
	subroutine patch_geom(P, s, t, r, nr, J)
		implicit none
		type(patch_CHEB):: P
		real(kind=8) s, t, r(3), nr(3), J, rs(3), rt(3), cr(3)
		call patch_eval(P, s, t, r, rs, rt)
		call z_rrcurl(rs, rt, cr)
		J = sqrt(sum(cr**2))
		nr = cr/J
	end subroutine patch_geom

	!**** chord D = r(s+ds,t+dt) - r(s,t) computed from the parameter offsets, avoiding the cancellation of subtracting O(1) coordinates
	subroutine patch_chord(P, s, t, ds, dt, D)
		implicit none
		type(patch_CHEB):: P
		real(kind=8) s, t, ds, dt, D(3)
		real(kind=8) f1, f2, df1, df2, dtan1, dtan2, PP(3), dP(3), PP1(3), np, np1, dnp, r0(3), r1(3), rs(3), rt(3)

		select case (P%kind)
		case (PATCH_AFFINE)
			D = ds*P%a1 + dt*P%a2
		case (PATCH_CUBEDSPHERE)
			f1 = BPACK_pi/4d0*(P%s0 + s*P%hs)
			f2 = BPACK_pi/4d0*(P%t0 + t*P%ht)
			df1 = BPACK_pi/4d0*P%hs*ds
			df2 = BPACK_pi/4d0*P%ht*dt
			dtan1 = sin(df1)/(cos(f1)*cos(f1 + df1)) ! tan(f1+df1)-tan(f1)
			dtan2 = sin(df2)/(cos(f2)*cos(f2 + df2))
			PP = P%c + tan(f1)*P%a1 + tan(f2)*P%a2
			dP = dtan1*P%a1 + dtan2*P%a2
			PP1 = PP + dP
			np = sqrt(sum(PP**2))
			np1 = sqrt(sum(PP1**2))
			dnp = -(2d0*dot_product(PP, dP) + dot_product(dP, dP))/(np + np1) ! |P|-|P1|
			D = P%rad*(dP/np1 + PP*dnp/(np*np1))
		case default ! planar patches: the normal component vanishes identically, plain subtraction suffices
			call patch_eval(P, s + ds, t + dt, r1, rs, rt)
			call patch_eval(P, s, t, r0, rs, rt)
			D = r1 - r0
		end select
	end subroutine patch_chord

	!**** closest point (sb,tb) of a patch to x, and the distance; Gauss-Newton with box constraints from the best point of a coarse grid
	subroutine patch_project(P, x, sb, tb, dist)
		implicit none
		type(patch_CHEB):: P
		real(kind=8) x(3), sb, tb, dist
		integer, parameter:: ng = 21
		integer ii, jj, it, ib
		real(kind=8) s, t, r(3), rs(3), rt(3), f(3), g1, g2, a11, a12, a22, det, ds, dt, snew, tnew, dnew, dcur
		logical fs, ft

		dist = BPACK_Bigvalue
		do jj = 1, ng
		do ii = 1, ng
			s = -1d0 + 2d0*(ii - 1)/(ng - 1d0)
			t = -1d0 + 2d0*(jj - 1)/(ng - 1d0)
			call patch_eval(P, s, t, r, rs, rt)
			if (sqrt(sum((r - x)**2)) < dist) then
				dist = sqrt(sum((r - x)**2))
				sb = s
				tb = t
			endif
		enddo
		enddo

		do it = 1, 100
			call patch_eval(P, sb, tb, r, rs, rt)
			f = r - x
			dcur = sqrt(sum(f**2))
			g1 = dot_product(rs, f)
			g2 = dot_product(rt, f)
			a11 = dot_product(rs, rs)
			a12 = dot_product(rs, rt)
			a22 = dot_product(rt, rt)
			det = a11*a22 - a12**2
			ds = -(a22*g1 - a12*g2)/det
			dt = -(-a12*g1 + a11*g2)/det
			fs = (sb <= -1d0 .and. ds < 0d0) .or. (sb >= 1d0 .and. ds > 0d0)
			ft = (tb <= -1d0 .and. dt < 0d0) .or. (tb >= 1d0 .and. dt > 0d0)
			if (fs .and. ft) then
				ds = -g1/a11
				dt = -g2/a22
			elseif (fs) then
				ds = 0d0
				dt = -g2/a22
			elseif (ft) then
				dt = 0d0
				ds = -g1/a11
			endif
			do ib = 1, 40 ! backtracking
				snew = min(1d0, max(-1d0, sb + ds))
				tnew = min(1d0, max(-1d0, tb + dt))
				call patch_eval(P, snew, tnew, r, rs, rt)
				dnew = sqrt(sum((r - x)**2))
				if (dnew <= dcur) exit
				ds = ds/2d0
				dt = dt/2d0
			enddo
			if (dnew > dcur .or. abs(snew - sb) + abs(tnew - tb) < 1d-15) exit
			sb = snew
			tb = tnew
		enddo
		call patch_eval(P, sb, tb, r, rs, rt)
		dist = sqrt(sum((r - x)**2))
	end subroutine patch_project

	!**** position (1..S*S) of sub-patch (is,js) within its base patch: row-major, or for the forest the Morton order of the
	! 2x2 patch hierarchy written as a binary tree that splits s first and then t at each level (so every tree node is a contiguous range)
	integer function subpatch_pos(is, js, S, forest)
		implicit none
		integer is, js, S, forest, b
		if (forest == 0) then
			subpatch_pos = (js - 1)*S + is
		else
			subpatch_pos = 1
			b = 0
			do while (2**b < S)
				subpatch_pos = subpatch_pos + ibits(is - 1, b, 1)*2**(2*b + 1) + ibits(js - 1, b, 1)*2**(2*b)
				b = b + 1
			enddo
		endif
	end function subpatch_pos

	!**** build the patches of an analytic geometry
	subroutine geo_modeling_CHEB(quant, ptree)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		type(z_proctree):: ptree
		real(kind=8) fc(3, 6), fa1(3, 6), fa2(3, 6), h, s, t, r(3), rs(3), rt(3), rot(2, 4)
		integer f, fpos, forder(6), is, js, q, ii, jj, k
		type(patch_CHEB):: P

		! face axes c, a1, a2 with a1 x a2 = c (outward normals)
		fc = reshape((/1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1/), (/3, 6/))
		fa1 = reshape((/0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0/), (/3, 6/))
		fa2 = reshape((/0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0/), (/3, 6/))
		rot = reshape((/1, 0, 0, 1, -1, 0, 0, -1/), (/2, 4/))
		forder = (/1, 2, 3, 4, 5, 6/)
		if (quant%faceorder == 1) forder = (/5, 1, 3, 2, 4, 6/)
		h = 1d0/quant%split
		if (quant%forest == 1 .and. iand(quant%split, quant%split - 1) /= 0) then
			if (ptree%MyID == Main_ID) write (*, *) 'the forest ordering requires split to be a power of two'
			stop
		endif

		select case (quant%geo)
		case (CHEB_SPHERE)
			quant%closed = .true.
			quant%npatch = 6*quant%split**2
			quant%nbase = 6
			allocate (quant%patches(quant%npatch))
			q = 0
			do fpos = 1, 6
			f = forder(fpos)
			do js = 1, quant%split
			do is = 1, quant%split
				q = (fpos - 1)*quant%split**2 + subpatch_pos(is, js, quant%split, quant%forest)
				P = patch_CHEB()
				P%kind = PATCH_CUBEDSPHERE
				P%planar = .false.
				P%c = fc(:, f)
				P%a1 = fa1(:, f)
				P%a2 = fa2(:, f)
				P%s0 = -1d0 + (2*is - 1)*h
				P%hs = h
				P%t0 = -1d0 + (2*js - 1)*h
				P%ht = h
				P%rad = quant%radius
				quant%patches(q) = P
			enddo
			enddo
			enddo
		case (CHEB_CUBE)
			quant%closed = .true.
			quant%npatch = 6*quant%split**2
			quant%nbase = 6
			allocate (quant%patches(quant%npatch))
			q = 0
			do fpos = 1, 6
			f = forder(fpos)
			do js = 1, quant%split
			do is = 1, quant%split
				q = (fpos - 1)*quant%split**2 + subpatch_pos(is, js, quant%split, quant%forest)
				P = patch_CHEB()
				P%kind = PATCH_AFFINE
				P%c = quant%radius*fc(:, f) + quant%radius*(-1d0 + (2*is - 1)*h)*fa1(:, f) + quant%radius*(-1d0 + (2*js - 1)*h)*fa2(:, f)
				P%a1 = quant%radius*h*fa1(:, f)
				P%a2 = quant%radius*h*fa2(:, f)
				P%edge = (/merge(1, 0, is == 1), merge(1, 0, is == quant%split), merge(1, 0, js == 1), merge(1, 0, js == quant%split)/)
				quant%patches(q) = P
			enddo
			enddo
			enddo
		case (CHEB_PLATE)
			quant%closed = .false.
			quant%npatch = quant%split**2
			quant%nbase = 1
			allocate (quant%patches(quant%npatch))
			q = 0
			do js = 1, quant%split
			do is = 1, quant%split
				q = subpatch_pos(is, js, quant%split, quant%forest)
				P = patch_CHEB()
				P%kind = PATCH_AFFINE
				P%c = (/quant%radius*(-1d0 + (2*is - 1)*h), quant%radius*(-1d0 + (2*js - 1)*h), 0d0/)
				P%a1 = (/quant%radius*h, 0d0, 0d0/)
				P%a2 = (/0d0, quant%radius*h, 0d0/)
				P%edge = (/merge(1, 0, is == 1), merge(1, 0, is == quant%split), merge(1, 0, js == 1), merge(1, 0, js == quant%split)/)
				quant%patches(q) = P
			enddo
			enddo
		case (CHEB_DISK) ! central square of half width radius/2 plus four blended side patches; the circle is the only edge
			quant%closed = .false.
			quant%npatch = 5*quant%split**2
			quant%nbase = 5
			allocate (quant%patches(quant%npatch))
			q = 0
			do js = 1, quant%split
			do is = 1, quant%split
				q = subpatch_pos(is, js, quant%split, quant%forest)
				P = patch_CHEB()
				P%kind = PATCH_AFFINE
				P%c = (/0.5d0*quant%radius*(-1d0 + (2*is - 1)*h), 0.5d0*quant%radius*(-1d0 + (2*js - 1)*h), 0d0/)
				P%a1 = (/0.5d0*quant%radius*h, 0d0, 0d0/)
				P%a2 = (/0d0, 0.5d0*quant%radius*h, 0d0/)
				quant%patches(q) = P
			enddo
			enddo
			do k = 1, 4
			do js = 1, quant%split
			do is = 1, quant%split
				q = k*quant%split**2 + subpatch_pos(is, js, quant%split, quant%forest)
				P = patch_CHEB()
				P%kind = PATCH_DISKSIDE
				P%a1 = (/rot(1, k), rot(2, k), 0d0/)
				P%a2 = (/-rot(2, k), rot(1, k), 0d0/)
				P%s0 = -1d0 + (2*is - 1)*h
				P%hs = h
				P%t0 = -1d0 + (2*js - 1)*h
				P%ht = h
				P%rad = quant%radius
				P%rin = 0.5d0*quant%radius
				P%edge = (/0, merge(1, 0, is == quant%split), 0, 0/)
				quant%patches(q) = P
			enddo
			enddo
			enddo
		end select

		! bounding spheres
		do q = 1, quant%npatch
			quant%patches(q)%bcenter = 0d0
			do jj = 1, 11
			do ii = 1, 11
				s = -1d0 + 0.2d0*(ii - 1)
				t = -1d0 + 0.2d0*(jj - 1)
				call patch_eval(quant%patches(q), s, t, r, rs, rt)
				quant%patches(q)%bcenter = quant%patches(q)%bcenter + r/121d0
			enddo
			enddo
			quant%patches(q)%bradius = 0d0
			do jj = 1, 11
			do ii = 1, 11
				s = -1d0 + 0.2d0*(ii - 1)
				t = -1d0 + 0.2d0*(jj - 1)
				call patch_eval(quant%patches(q), s, t, r, rs, rt)
				quant%patches(q)%bradius = max(quant%patches(q)%bradius, sqrt(sum((r - quant%patches(q)%bcenter)**2)))
			enddo
			enddo
			quant%patches(q)%bradius = 1.1d0*quant%patches(q)%bradius
		enddo

		if (quant%unknown < 0) quant%unknown = merge(UNK_PHI, UNK_PSI, quant%closed)
		if (ptree%MyID == Main_ID) then
			write (*, *) 'number of patches:', quant%npatch, ' closed surface:', quant%closed
		endif
	end subroutine geo_modeling_CHEB

	!**** Fejer nodes on all patches
	subroutine build_nodes_CHEB(quant, ptree)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		type(z_proctree):: ptree
		integer N, N2, q, i, j, idx, nn
		real(kind=8) s, t, es, et, r(3), nr(3), Jac, area, hpatch

		N = quant%N
		N2 = N*N
		quant%myid = ptree%MyID
		allocate (quant%u(N), quant%wu(N))
		call fejer_rule(N, quant%u, quant%wu)
		allocate (quant%tb(quant%Nbeta), quant%wb(quant%Nbeta))
		call fejer_rule(quant%Nbeta, quant%tb, quant%wb)
		allocate (quant%Ccheb(0:N - 1, N))
		do i = 1, N
		do nn = 0, N - 1
			quant%Ccheb(nn, i) = merge(1d0, 2d0, nn == 0)/N*cos(nn*BPACK_pi*(2*i - 1)/(2d0*N))
		enddo
		enddo

		quant%Nunk = quant%npatch*N2
		allocate (quant%xyz(3, quant%Nunk), quant%nrm(3, quant%Nunk), quant%st(2, quant%Nunk), quant%Jw(quant%Nunk), quant%e(quant%Nunk))
		do q = 1, quant%npatch
		do j = 1, N
		do i = 1, N
			idx = (q - 1)*N2 + (j - 1)*N + i
			call eta_map(quant%u(i), quant%patches(q)%edge(1), quant%patches(q)%edge(2), quant%pedge, s, es)
			call eta_map(quant%u(j), quant%patches(q)%edge(3), quant%patches(q)%edge(4), quant%pedge, t, et)
			call patch_geom(quant%patches(q), s, t, r, nr, Jac)
			quant%xyz(:, idx) = r
			quant%nrm(:, idx) = nr
			quant%st(:, idx) = (/s, t/)
			quant%Jw(idx) = Jac*quant%wu(i)*quant%wu(j)
			quant%e(idx) = es*et
		enddo
		enddo
		enddo

		area = sum(quant%Jw*quant%e)
		hpatch = sqrt(area/quant%npatch)
		if (quant%delta < 0) quant%delta = 0.5d0*hpatch
		if (ptree%MyID == Main_ID) then
			write (*, *) 'Nunk:', quant%Nunk, ' N:', N, ' Nbeta:', quant%Nbeta
			write (*, *) 'surface area:', area, ' typical patch size:', hpatch, ' delta:', quant%delta
			write (*, *) 'points per wavelength:', N/(hpatch/quant%wavelength), ' min edge factor:', minval(quant%e)
			write (*, *) 'unknown: ', merge('phi', 'psi', quant%unknown == UNK_PHI), '  p_xi:', quant%pxi, ' p_edge:', quant%pedge
		endif
	end subroutine build_nodes_CHEB

	!**** rectangular-polar weights W(i,j) = int int H(x_l, r(u,v)) J(u,v) l_i(u) l_j(v) du dv for target l and patch q, eqs. (37),(41)
	subroutine pair_weights_CHEB(quant, l, q, ubar, vbar, own, W)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		integer l, q, own, a, b, N, Nb
		real(kind=8) ubar, vbar, x(3), y(3), nr(3), J, D(3), R, dn, s0, t0, dummy
		complex(kind=8) W(quant%N*quant%N)
		real(kind=8), allocatable:: offu(:), muu(:), offv(:), muv(:), ut(:), vt(:), sa(:), tv(:), Lu(:, :), Lv(:, :)
		complex(kind=8), allocatable:: F(:, :), tmp(:, :)
		type(patch_CHEB):: P

		N = quant%N
		Nb = quant%Nbeta
		P = quant%patches(q)
		allocate (offu(Nb), muu(Nb), offv(Nb), muv(Nb), ut(Nb), vt(Nb), sa(Nb), tv(Nb), Lu(Nb, N), Lv(Nb, N), F(Nb, Nb), tmp(Nb, N))
		do a = 1, Nb
			call xi_map(quant%tb(a), ubar, quant%pxi, offu(a), muu(a))
			call xi_map(quant%tb(a), vbar, quant%pxi, offv(a), muv(a))
			ut(a) = min(1d0, max(-1d0, ubar + offu(a)))
			vt(a) = min(1d0, max(-1d0, vbar + offv(a)))
			call eta_map(ut(a), P%edge(1), P%edge(2), quant%pedge, sa(a), dummy)
			call eta_map(vt(a), P%edge(3), P%edge(4), quant%pedge, tv(a), dummy)
		enddo
		call cheb_lagrange(N, quant%Ccheb, Nb, ut, Lu)
		call cheb_lagrange(N, quant%Ccheb, Nb, vt, Lv)

		x = quant%xyz(:, l)
		s0 = quant%st(1, l)
		t0 = quant%st(2, l)
		do b = 1, Nb
		do a = 1, Nb
			call patch_geom(P, sa(a), tv(b), y, nr, J)
			if (own == 1 .and. quant%naive_chord == 0) then
				call patch_chord(P, s0, t0, sa(a) - s0, tv(b) - t0, D) ! D = y - x
				R = sqrt(sum(D**2))
				dn = -dot_product(D, nr)
				if (P%planar) dn = 0d0
			else
				D = x - y
				R = sqrt(sum(D**2))
				dn = dot_product(D, nr)
			endif
			F(a, b) = kernel_H(quant%wavenum, R, dn, quant%closed)*J*muu(a)*muv(b)*quant%wb(a)*quant%wb(b)
		enddo
		enddo
		tmp = matmul(F, Lv)
		W = reshape(matmul(transpose(Lu), tmp), (/N*N/))
		deallocate (offu, muu, offv, muv, ut, vt, sa, tv, Lu, Lv, F, tmp)
	end subroutine pair_weights_CHEB

	!**** find all (target, patch) self and near pairs, compute their weights in parallel, and replicate them on all ranks
	subroutine precompute_nearfield_CHEB(quant, ptree)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		type(z_proctree):: ptree
		integer N2, nproc, myid, ls, le, l, q, ql, npl, cap, k, ierr, r
		integer, allocatable:: cnt_loc(:), pl(:), pq(:), pown(:), tcounts(:), tdispls(:), pcounts(:), pdispls(:), wcounts(:), wdispls(:), itmp(:), nnear(:)
		real(kind=8), allocatable:: pub(:), pvb(:), rtmp(:)
		complex(kind=8), allocatable:: Wloc(:, :)
		real(kind=8) sb, tb, dist, t1, t2, bytes
		type(patch_CHEB):: P

		t1 = MPI_Wtime()
		N2 = quant%N*quant%N
		nproc = ptree%nproc
		myid = ptree%MyID
		ls = int((int(quant%Nunk, 8)*myid)/nproc) + 1
		le = int((int(quant%Nunk, 8)*(myid + 1))/nproc)

		!**** pass 1: detect the pairs of the local targets (self pair first)
		cap = max(16, 10*(le - ls + 1))
		allocate (pl(cap), pq(cap), pown(cap), pub(cap), pvb(cap))
		allocate (cnt_loc(max(1, le - ls + 1)))
		cnt_loc = 0
		npl = 0
		do l = ls, le
			ql = (l - 1)/N2 + 1
			do q = 1, quant%npatch
				P = quant%patches(q)
				if (q == ql) then
					k = l - (ql - 1)*N2
					sb = quant%u(mod(k - 1, quant%N) + 1)
					tb = quant%u((k - 1)/quant%N + 1)
				else
					if (sqrt(sum((quant%xyz(:, l) - P%bcenter)**2)) - P%bradius >= quant%delta) cycle
					call patch_project(P, quant%xyz(:, l), sb, tb, dist)
					if (dist >= quant%delta) cycle
					sb = eta_inv(sb, P%edge(1), P%edge(2), quant%pedge)
					tb = eta_inv(tb, P%edge(3), P%edge(4), quant%pedge)
				endif
				if (npl == cap) then
					cap = 2*cap
					allocate (itmp(cap)); itmp(1:npl) = pl(1:npl); call move_alloc(itmp, pl)
					allocate (itmp(cap)); itmp(1:npl) = pq(1:npl); call move_alloc(itmp, pq)
					allocate (itmp(cap)); itmp(1:npl) = pown(1:npl); call move_alloc(itmp, pown)
					allocate (rtmp(cap)); rtmp(1:npl) = pub(1:npl); call move_alloc(rtmp, pub)
					allocate (rtmp(cap)); rtmp(1:npl) = pvb(1:npl); call move_alloc(rtmp, pvb)
				endif
				npl = npl + 1
				pl(npl) = l
				pq(npl) = q
				pown(npl) = merge(1, 0, q == ql)
				pub(npl) = sb
				pvb(npl) = tb
				cnt_loc(l - ls + 1) = cnt_loc(l - ls + 1) + 1
			enddo
		enddo

		!**** pass 2: rectangular-polar weights of the local pairs (precomputed mode only)
		if (quant%onthefly == 0) then
		allocate (Wloc(N2, max(1, npl)))
#ifdef HAVE_OPENMP
		!$omp parallel do default(shared) private(k) schedule(dynamic,4)
#endif
		do k = 1, npl
			call pair_weights_CHEB(quant, pl(k), pq(k), pub(k), pvb(k), pown(k), Wloc(:, k))
		enddo
#ifdef HAVE_OPENMP
		!$omp end parallel do
#endif
		endif

		!**** replicate on all ranks: pair counts per target, pair patches, and weights
		allocate (tcounts(nproc), tdispls(nproc), pcounts(nproc), pdispls(nproc), wcounts(nproc), wdispls(nproc))
		do r = 0, nproc - 1
			tdispls(r + 1) = int((int(quant%Nunk, 8)*r)/nproc)
			tcounts(r + 1) = int((int(quant%Nunk, 8)*(r + 1))/nproc) - tdispls(r + 1)
		enddo
		allocate (nnear(quant%Nunk))
		call MPI_ALLGATHERV(cnt_loc, le - ls + 1, MPI_INTEGER, nnear, tcounts, tdispls, MPI_INTEGER, ptree%Comm, ierr)
		call MPI_ALLGATHER(npl, 1, MPI_INTEGER, pcounts, 1, MPI_INTEGER, ptree%Comm, ierr)
		pdispls(1) = 0
		do r = 2, nproc
			pdispls(r) = pdispls(r - 1) + pcounts(r - 1)
		enddo
		quant%npair = pdispls(nproc) + pcounts(nproc)
		if (quant%onthefly == 0 .and. dble(quant%npair)*N2 > dble(huge(0))) then
			if (myid == Main_ID) write (*, *) 'near-field weights exceed the 32-bit MPI count limit; use fewer nodes per patch'
			stop
		endif
		allocate (quant%near_ptr(quant%Nunk + 1))
		quant%near_ptr(1) = 1
		do l = 1, quant%Nunk
			quant%near_ptr(l + 1) = quant%near_ptr(l) + nnear(l)
		enddo
		allocate (quant%near_q(quant%npair))
		call MPI_ALLGATHERV(pq, npl, MPI_INTEGER, quant%near_q, pcounts, pdispls, MPI_INTEGER, ptree%Comm, ierr)
		if (quant%onthefly == 0) then
			wcounts = pcounts*N2
			wdispls = pdispls*N2
			allocate (quant%near_W(int(quant%npair, 8)*N2))
			call MPI_ALLGATHERV(Wloc, npl*N2, MPI_DOUBLE_COMPLEX, quant%near_W, wcounts, wdispls, MPI_DOUBLE_COMPLEX, ptree%Comm, ierr)
			deallocate (Wloc)
		else ! keep only the closest points; the weights are computed when their entries are requested
			allocate (quant%near_ub(quant%npair), quant%near_vb(quant%npair))
			allocate (quant%paircount(quant%npair))
			quant%paircount = 0
			call MPI_ALLGATHERV(pub, npl, MPI_DOUBLE_PRECISION, quant%near_ub, pcounts, pdispls, MPI_DOUBLE_PRECISION, ptree%Comm, ierr)
			call MPI_ALLGATHERV(pvb, npl, MPI_DOUBLE_PRECISION, quant%near_vb, pcounts, pdispls, MPI_DOUBLE_PRECISION, ptree%Comm, ierr)
		endif

		deallocate (pl, pq, pown, pub, pvb, cnt_loc, tcounts, tdispls, pcounts, pdispls, wcounts, wdispls, nnear)
		t2 = MPI_Wtime()
		if (myid == Main_ID) then
			write (*, *) 'near-field pairs:', quant%npair, ' (', dble(quant%npair)/quant%Nunk, ' per target)'
			if (quant%onthefly == 0) then
				bytes = dble(quant%npair)*N2*16d0
				write (*, '(A,F10.3,A,F10.3,A)') ' near-field weights: ', bytes/1024d0**2, ' MB per rank, precomputed in ', t2 - t1, ' seconds'
			else
				bytes = dble(quant%npair)*20d0
				write (*, '(A,F10.3,A,F10.3,A)') ' near-field pairs (weights computed on the fly): ', bytes/1024d0**2, ' MB per rank, found in ', t2 - t1, ' seconds'
			endif
		endif
	end subroutine precompute_nearfield_CHEB

	!**** spherical Bessel functions j_l(x), y_l(x), l=0..lmax, for real x>0
	subroutine sph_bessel_jy(lmax, x, jl, yl)
		implicit none
		integer lmax, l, Ls
		real(kind=8) x, jl(0:lmax), yl(0:lmax), j0, j1, sc
		real(kind=8), allocatable:: f(:)

		yl(0) = -cos(x)/x
		if (lmax >= 1) yl(1) = -cos(x)/x**2 - sin(x)/x
		do l = 1, lmax - 1
			yl(l + 1) = (2*l + 1)/x*yl(l) - yl(l - 1)
		enddo

		Ls = lmax + 30 + int(x) + int(sqrt(40d0*(lmax + x)))
		allocate (f(0:Ls + 1))
		f(Ls + 1) = 0d0
		f(Ls) = 1d0
		do l = Ls, 1, -1 ! Miller's downward recurrence
			f(l - 1) = (2*l + 1)/x*f(l) - f(l + 1)
			if (abs(f(l - 1)) > 1d200) f(l - 1:Ls + 1) = f(l - 1:Ls + 1)*1d-200
		enddo
		j0 = sin(x)/x
		j1 = sin(x)/x**2 - cos(x)/x
		if (abs(j0) > abs(j1)) then
			sc = j0/f(0)
		else
			sc = j1/f(1)
		endif
		jl(0:lmax) = f(0:lmax)*sc
		deallocate (f)
	end subroutine sph_bessel_jy

	real(kind=8) function legendreP(l, x)
		implicit none
		integer l, nn
		real(kind=8) x, p0, p1, p2
		p0 = 1d0
		p1 = x
		if (l == 0) then
			legendreP = p0
			return
		endif
		do nn = 1, l - 1
			p2 = ((2*nn + 1)*x*p1 - nn*p0)/(nn + 1)
			p0 = p1
			p1 = p2
		enddo
		legendreP = p1
	end function legendreP

	!**** exact eigenvalue of 1/2 I + D - ik S on the sphere of radius a for spherical harmonics of degree l
	complex(kind=8) function sphere_eigenvalue(k, a, l)
		implicit none
		real(kind=8) k, a
		integer l
		real(kind=8) jl(0:l + 1), yl(0:l + 1), dj, dy
		complex(kind=8) h, dh, lamS, lamD
		call sph_bessel_jy(l + 1, k*a, jl, yl)
		dj = -jl(l + 1) + l/(k*a)*jl(l) ! j_l'(x) = l/x j_l - j_{l+1}
		dy = -yl(l + 1) + l/(k*a)*yl(l)
		h = jl(l) + BPACK_junit*yl(l)
		dh = dj + BPACK_junit*dy
		lamS = BPACK_junit*k*a**2*jl(l)*h
		lamD = BPACK_junit*k**2*a**2/2d0*(jl(l)*dh + h*dj)
		sphere_eigenvalue = 0.5d0 + lamD - BPACK_junit*k*lamS
	end function sphere_eigenvalue

	!**** far-field pattern of a sound-soft sphere of radius a for the incident wave exp(ik d.x); cosg = xhat.d
	complex(kind=8) function mie_farfield_soft(k, a, cosg)
		implicit none
		real(kind=8) k, a, cosg
		integer l, lmax
		real(kind=8), allocatable:: jl(:), yl(:)
		complex(kind=8) val
		lmax = int(k*a + 4.05d0*(k*a)**(1d0/3d0)) + 15
		allocate (jl(0:lmax), yl(0:lmax))
		call sph_bessel_jy(lmax, k*a, jl, yl)
		val = 0d0
		do l = 0, lmax
			val = val + (2*l + 1)*jl(l)/(jl(l) + BPACK_junit*yl(l))*legendreP(l, cosg)
		enddo
		mie_farfield_soft = BPACK_junit/k*val
		deallocate (jl, yl)
	end function mie_farfield_soft

	!**** forward-map test on the sphere: Y = P_l(ehat.xhat) at all nodes and the exact eigenvalue of the operator for degree l
	subroutine eigtest_Y_CHEB(quant, Y, lam)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		real(kind=8), allocatable:: Y(:)
		complex(kind=8) lam
		real(kind=8) ehat(3)
		integer n
		ehat = (/0.3d0, -0.5d0, 0.8d0/)
		ehat = ehat/sqrt(sum(ehat**2))
		lam = sphere_eigenvalue(quant%wavenum, quant%radius, quant%eigtest)
		allocate (Y(quant%Nunk))
		do n = 1, quant%Nunk
			Y(n) = legendreP(quant%eigtest, dot_product(ehat, quant%xyz(:, n))/sqrt(sum(quant%xyz(:, n)**2)))
		enddo
	end subroutine eigtest_Y_CHEB

	!**** forward-map test report: vout(i) is the compressed operator applied to Y at the local node gi(i); also applies the
	! uncompressed entries on a sample of the local rows
	subroutine forward_map_report_CHEB(quant, ptree, Nloc, gi, vout, Y, lam)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		type(z_proctree):: ptree
		integer Nloc, gi(Nloc), ii, n, nsamp, kk, ierr
		complex(kind=8) vout(Nloc), lam, val
		real(kind=8) Y(:), errs(4), errg(4)

		errs = 0d0
		do ii = 1, Nloc
			errs(1) = max(errs(1), abs(vout(ii) - lam*Y(gi(ii))))
			errs(2) = max(errs(2), abs(lam*Y(gi(ii))))
		enddo
		nsamp = min(Nloc, 16)
		do kk = 1, nsamp
			ii = 1 + ((kk - 1)*Nloc)/nsamp
			val = 0d0
#ifdef HAVE_OPENMP
			!$omp parallel do default(shared) private(n) reduction(+:val)
#endif
			do n = 1, quant%Nunk
				if (quant%unknown == UNK_PSI) then
					val = val + Zentry_CHEB(quant, gi(ii), n)*Y(n)*quant%e(n)
				else
					val = val + Zentry_CHEB(quant, gi(ii), n)*Y(n)
				endif
			enddo
#ifdef HAVE_OPENMP
			!$omp end parallel do
#endif
			errs(3) = max(errs(3), abs(val - lam*Y(gi(ii))))
			errs(4) = max(errs(4), abs(lam*Y(gi(ii))))
		enddo
		call MPI_ALLREDUCE(errs, errg, 4, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
		if (ptree%MyID == Main_ID) then
			write (*, *) ''
			write (*, '(A,I3,A,2Es14.6)') ' forward-map test with Y of degree', quant%eigtest, ', exact eigenvalue:', dble(lam), aimag(lam)
			write (*, '(A,Es12.4)') '   max error, uncompressed entries (sampled rows): ', errg(3)/errg(4)
			write (*, '(A,Es12.4)') '   max error, compressed operator (all rows):      ', errg(1)/errg(2)
			write (*, *) ''
		endif
	end subroutine forward_map_report_CHEB

	!**** forward-map test on the sphere: apply the operator to Y = P_l(ehat.xhat) and compare with the exact eigenvalue
	subroutine forward_map_test_CHEB(bmat, option, msh, quant, ptree, stats)
		use z_BPACK_Solve_Mul
		implicit none
		type(z_Bmatrix):: bmat
		type(z_Hoption):: option
		type(z_mesh):: msh
		type(quant_ACOUSTIC_CHEB):: quant
		type(z_proctree):: ptree
		type(z_Hstat):: stats
		integer Nloc, ii
		integer, allocatable:: gi(:)
		real(kind=8), allocatable:: Y(:)
		complex(kind=8), allocatable:: vin(:, :), vout(:, :)
		complex(kind=8) lam

		if (quant%geo /= CHEB_SPHERE .or. quant%eigtest <= 0) return
		call eigtest_Y_CHEB(quant, Y, lam)
		Nloc = msh%idxe - msh%idxs + 1
		allocate (gi(Nloc), vin(Nloc, 1), vout(Nloc, 1))
		do ii = 1, Nloc
			gi(ii) = msh%new2old(msh%idxs + ii - 1)
			vin(ii, 1) = Y(gi(ii))
			if (quant%unknown == UNK_PSI) vin(ii, 1) = vin(ii, 1)*quant%e(gi(ii))
		enddo
		call z_BPACK_Mult('N', Nloc, 1, vin, vout, bmat, ptree, option, stats, use_blockcopy=0) ! forward blocks, not their post-factorization copy
		call forward_map_report_CHEB(quant, ptree, Nloc, gi, vout(:, 1), Y, lam)
		deallocate (Y, gi, vin, vout)
	end subroutine forward_map_test_CHEB

	!**** GMRES(restart) on the compressed operator (precon=NOPRECON) or preconditioned by its factorization (precon=BPACKPRECON); same as BPACK_Z_iter but with a user-chosen restart length
	subroutine gmres_solve_CHEB(bmat, option, ptree, stats, Nloc, b, x, restart, iter, relres)
		use z_BPACK_Solve_Mul
		implicit none
		type(z_Bmatrix), target:: bmat
		type(z_Hoption), target:: option
		type(z_proctree), target:: ptree
		type(z_Hstat), target:: stats
		integer Nloc, restart, iter
		real(kind=8) relres
		complex(kind=8) b(Nloc, 1), x(Nloc, 1)
		type(z_kernelquant), target:: kerb
		type(z_quant_bmat), target:: qb

		nullify (qb%msh, qb%msh_md)
		qb%bmat => bmat
		qb%ptree => ptree
		qb%option => option
		qb%stats => stats
		qb%ker => kerb
		kerb%QuantApp => qb
		iter = 0
		relres = option%tol_itersol
		call z_BPACK_Zgmres_usermatvec_precon(option%n_iter, Nloc, b, x, relres, iter, z_blackbox_BPACK_MVP, z_blackbox_BPACK_precon_MVP, &
			&	ptree, option, stats, kerb, restart_in=restart)
	end subroutine gmres_solve_CHEB

	!**** solve for a plane wave and evaluate the far field on the great circle containing the incident direction
	subroutine solve_scattering_CHEB(bmat, option, msh, quant, ptree, stats)
		use z_BPACK_Solve_Mul
		implicit none
		type(z_Bmatrix):: bmat
		type(z_Hoption):: option
		type(z_mesh):: msh
		type(quant_ACOUSTIC_CHEB):: quant
		type(z_proctree):: ptree
		type(z_Hstat):: stats
		integer Nloc, ii, iter
		integer, allocatable:: gi(:)
		real(kind=8) t1, t2, relres
		complex(kind=8), allocatable:: x(:, :), b(:, :)

		Nloc = msh%idxe - msh%idxs + 1
		allocate (gi(Nloc), x(Nloc, 1), b(Nloc, 1))
		do ii = 1, Nloc
			gi(ii) = msh%new2old(msh%idxs + ii - 1)
			b(ii, 1) = -exp(BPACK_junit*quant%wavenum*dot_product(quant%dinc, quant%xyz(:, gi(ii))))
		enddo
		x = 0d0
		t1 = MPI_Wtime()
		if (option%precon /= DIRECT .and. quant%gmres_restart > 0) then
			call gmres_solve_CHEB(bmat, option, ptree, stats, Nloc, b, x, quant%gmres_restart, iter, relres)
			t2 = MPI_Wtime()
			stats%Time_Sol = stats%Time_Sol + t2 - t1
			if (ptree%MyID == Main_ID) write (*, '(A,I5,A,I6,A,Es12.4)') ' GMRES(', quant%gmres_restart, ') iterations:', iter, '   relative residual:', relres
		else
			call z_BPACK_Solution(bmat, x, b, Nloc, 1, option, ptree, stats)
			t2 = MPI_Wtime()
		endif
		if (ptree%MyID == Main_ID .and. option%verbosity >= 0) write (*, *) 'Solving:', t2 - t1, 'Seconds'
		call farfield_CHEB(quant, ptree, Nloc, gi, x(:, 1))
		deallocate (gi, x, b)
	end subroutine solve_scattering_CHEB

	!**** far field of the solution x (x(i) at the local node gi(i)) on the great circle containing the incident direction,
	! written to farfield_CHEB.out, and its error against the Mie series on the sphere
	subroutine farfield_CHEB(quant, ptree, Nloc, gi, x)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		type(z_proctree):: ptree
		integer Nloc, gi(Nloc), ii, m, ierr
		complex(kind=8) x(Nloc)
		real(kind=8) e1(3), xh(3), gam, errmax, fmax, k
		complex(kind=8), allocatable:: ffloc(:), ff(:)
		complex(kind=8) coef, umie

		k = quant%wavenum
		! a unit vector orthogonal to the incident direction
		if (abs(quant%dinc(1)) < 0.9d0) then
			call z_rrcurl(quant%dinc, (/1d0, 0d0, 0d0/), e1)
		else
			call z_rrcurl(quant%dinc, (/0d0, 1d0, 0d0/), e1)
		endif
		e1 = e1/sqrt(sum(e1**2))

		allocate (ffloc(quant%nfar), ff(quant%nfar))
		ffloc = 0d0
		do m = 1, quant%nfar
			gam = BPACK_pi*(m - 1)/max(1, quant%nfar - 1)
			xh = cos(gam)*quant%dinc + sin(gam)*e1
			do ii = 1, Nloc
				if (quant%closed) then
					coef = -BPACK_junit*k*(dot_product(xh, quant%nrm(:, gi(ii))) + 1d0)
				else
					coef = 1d0
				endif
				coef = coef*exp(-BPACK_junit*k*dot_product(xh, quant%xyz(:, gi(ii))))*quant%Jw(gi(ii))/(4d0*BPACK_pi)
				if (quant%unknown == UNK_PHI) coef = coef*quant%e(gi(ii)) ! psi = e*phi
				ffloc(m) = ffloc(m) + coef*x(ii)
			enddo
		enddo
		call MPI_ALLREDUCE(ffloc, ff, quant%nfar, MPI_DOUBLE_COMPLEX, MPI_SUM, ptree%Comm, ierr)

		if (ptree%MyID == Main_ID) then
			open (100, file='farfield_CHEB.out')
			errmax = 0d0
			fmax = 0d0
			do m = 1, quant%nfar
				gam = BPACK_pi*(m - 1)/max(1, quant%nfar - 1)
				if (quant%geo == CHEB_SPHERE) then
					umie = mie_farfield_soft(k, quant%radius, cos(gam))
					errmax = max(errmax, abs(ff(m) - umie))
					fmax = max(fmax, abs(umie))
					write (100, '(F10.4,5Es18.9)') gam*180d0/BPACK_pi, dble(ff(m)), aimag(ff(m)), abs(ff(m)), dble(umie), aimag(umie)
				else
					write (100, '(F10.4,3Es18.9)') gam*180d0/BPACK_pi, dble(ff(m)), aimag(ff(m)), abs(ff(m))
				endif
			enddo
			close (100)
			write (*, *) 'far field written to farfield_CHEB.out'
			if (any(ff /= ff)) write (*, *) 'WARNING: the far field contains NaN, the solve failed'
			if (quant%geo == CHEB_SPHERE) write (*, '(A,Es12.4)') ' far-field error vs. Mie series (relative to max): ', errmax/fmax
		endif
		deallocate (ffloc, ff)
	end subroutine farfield_CHEB

	!**** on-the-fly near field: how often the weight rows were computed during the construction, and whether every
	! self/near pair was computed exactly once in a dense leaf block (forest and tensor forest)
	subroutine nearfield_report_CHEB(quant, ptree)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		type(z_proctree):: ptree
		integer(kind=8) nweights, ncat(3), nent(2)
		integer ii, kk, nshow, ierr

		call MPI_ALLREDUCE((/quant%nentries, quant%nsubmat/), nent, 2, MPI_INTEGER8, MPI_SUM, ptree%Comm, ierr)
		if (ptree%MyID == Main_ID) write (*, '(A,I14,A,F10.4,A,I12,A)') ' matrix entries requested during construction:', nent(1), &
			&	' (', dble(nent(1))/(dble(quant%Nunk)**2), ' of the dense matrix) in', nent(2), ' submatrices'
		if (quant%onthefly /= 1) return
		call MPI_ALLREDUCE(quant%nweights_computed, nweights, 1, MPI_INTEGER8, MPI_SUM, ptree%Comm, ierr)
		call MPI_ALLREDUCE(quant%nweights_cat, ncat, 3, MPI_INTEGER8, MPI_SUM, ptree%Comm, ierr)
		if (ptree%MyID == Main_ID) write (*, '(A,I12,A,F10.3,A)') ' near-field weight rows computed during construction:', nweights, &
			&	' (', dble(nweights)/quant%npair, ' per near pair)'
		if (ptree%MyID == Main_ID) write (*, '(A,3I12)') '   requested by leaf-patch blocks / other submatrices / single entries:', ncat
		call MPI_ALLREDUCE(MPI_IN_PLACE, quant%paircount, quant%npair, MPI_INTEGER, MPI_SUM, ptree%Comm, ierr)
		if (ptree%MyID == Main_ID) then
			write (*, '(A,4I10)') '   near pairs computed 0/1/2/>2 times in leaf-patch blocks:', count(quant%paircount == 0), count(quant%paircount == 1), &
				&	count(quant%paircount == 2), count(quant%paircount > 2)
			if ((quant%forest == 1 .or. quant%tensor == 1) .and. count(quant%paircount /= 1) > 0) write (*, '(A,I10,A)') ' WARNING:', count(quant%paircount /= 1), &
				&	' near pairs are not in exactly one dense leaf block (some lie in admissible blocks); increase -option --near_para'
			nshow = 0
			do ii = 1, quant%Nunk
				do kk = quant%near_ptr(ii), quant%near_ptr(ii + 1) - 1
					if (quant%paircount(kk) > 1 .and. nshow < 5) then
						nshow = nshow + 1
						write (*, '(A,I8,A,I6,A,I6,A,I3)') '     e.g. target', ii, ' (patch', (ii - 1)/quant%N**2 + 1, ') and patch', quant%near_q(kk), ': computed', quant%paircount(kk)
					endif
				enddo
			enddo
		endif
	end subroutine nearfield_report_CHEB

	!======================================================================================== tensor forest (quant%tensor=1)
	! Each base patch is split into split x split patches with N x N Fejer nodes each, so its nodes form a product grid of
	! split*N points along s times split*N points along t of the base patch. The interaction of base patches p and q is a
	! 4-mode tensor on these grids, compressed as one H-TENSOR (option%format=HTENSOR) per pair (p,q): the per-dimension
	! cluster trees are the binary trees of the split (a leaf is one patch), so every quadtree node of the patch hierarchy is
	! a product of one node per dimension. The product with the whole matrix sums the nbase x nbase H-TENSOR products.

	!**** (is,js) of the patch at position pos (1..S*S) of its base patch: the inverse of subpatch_pos
	subroutine subpatch_ij(pos, S, forest, is, js)
		implicit none
		integer pos, S, forest, is, js, b
		if (forest == 0) then
			is = mod(pos - 1, S) + 1
			js = (pos - 1)/S + 1
		else
			is = 1
			js = 1
			b = 0
			do while (2**b < S)
				is = is + ibits(pos - 1, 2*b + 1, 1)*2**b
				js = js + ibits(pos - 1, 2*b, 1)*2**b
				b = b + 1
			enddo
		endif
	end subroutine subpatch_ij

	!**** node (natural order) of the grid point (isig,itau) of base patch p, isig = (is-1)*N+i along s and itau = (js-1)*N+j along t
	integer function tnode_CHEB(quant, p, isig, itau)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		integer p, isig, itau, is, js, N
		N = quant%N
		is = (isig - 1)/N + 1
		js = (itau - 1)/N + 1
		tnode_CHEB = ((p - 1)*quant%split**2 + subpatch_pos(is, js, quant%split, quant%forest) - 1)*N*N &
			&	+ (itau - (js - 1)*N - 1)*N + isig - (is - 1)*N
	end function tnode_CHEB

	!**** index of the quadtree node at level l and position (a,b) (along s and t) of a base patch
	integer function tnode_id(l, a, b)
		implicit none
		integer l, a, b
		tnode_id = (4**l - 1)/3 + (b - 1)*2**l + a
	end function tnode_id

	!**** tensor forest: the split position of every patch, the patches holding self/near pairs of the targets of every
	! patch, and the bounding spheres of the quadtree nodes of every base patch
	subroutine tensor_setup_CHEB(quant, ptree)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		type(z_proctree):: ptree
		integer S, Lt, N2, q, pos, is, js, t, m, kk, k, l, a, b, c, id, p
		integer, allocatable:: mark(:)
		real(kind=8) cen(3), r

		S = quant%split
		if (iand(S, S - 1) /= 0) then
			if (ptree%MyID == Main_ID) write (*, *) 'the tensor forest requires split to be a power of two'
			stop
		endif
		Lt = 0
		do while (2**Lt < S)
			Lt = Lt + 1
		enddo
		quant%tlev = Lt
		N2 = quant%N**2

		allocate (quant%sp_face(quant%npatch), quant%sp_is(quant%npatch), quant%sp_js(quant%npatch))
		do q = 1, quant%npatch
			quant%sp_face(q) = (q - 1)/S**2 + 1
			pos = q - (quant%sp_face(q) - 1)*S**2
			call subpatch_ij(pos, S, quant%forest, is, js)
			quant%sp_is(q) = is
			quant%sp_js(q) = js
		enddo

		allocate (mark(quant%npatch), quant%nsub_ptr(quant%npatch + 1))
		mark = 0
		quant%nsub_ptr(1) = 1
		do t = 1, quant%npatch
			k = 0
			do m = (t - 1)*N2 + 1, t*N2
				do kk = quant%near_ptr(m), quant%near_ptr(m + 1) - 1
					if (mark(quant%near_q(kk)) /= t) then
						mark(quant%near_q(kk)) = t
						k = k + 1
					endif
				enddo
			enddo
			quant%nsub_ptr(t + 1) = quant%nsub_ptr(t) + k
		enddo
		allocate (quant%nsub_q(quant%nsub_ptr(quant%npatch + 1) - 1))
		mark = 0
		do t = 1, quant%npatch
			k = quant%nsub_ptr(t) - 1
			do m = (t - 1)*N2 + 1, t*N2
				do kk = quant%near_ptr(m), quant%near_ptr(m + 1) - 1
					if (mark(quant%near_q(kk)) /= t) then
						mark(quant%near_q(kk)) = t
						k = k + 1
						quant%nsub_q(k) = quant%near_q(kk)
					endif
				enddo
			enddo
		enddo
		deallocate (mark)

		allocate (quant%tnode_c(3, (4**(Lt + 1) - 1)/3, quant%nbase), quant%tnode_r((4**(Lt + 1) - 1)/3, quant%nbase))
		do p = 1, quant%nbase
			do js = 1, S
				do is = 1, S
					q = (p - 1)*S**2 + subpatch_pos(is, js, S, quant%forest)
					quant%tnode_c(:, tnode_id(Lt, is, js), p) = quant%patches(q)%bcenter
					quant%tnode_r(tnode_id(Lt, is, js), p) = quant%patches(q)%bradius
				enddo
			enddo
			do l = Lt - 1, 0, -1
				do b = 1, 2**l
					do a = 1, 2**l
						cen = 0d0
						do c = 0, 3
							cen = cen + quant%tnode_c(:, tnode_id(l + 1, 2*a - 1 + mod(c, 2), 2*b - 1 + c/2), p)/4d0
						enddo
						r = 0d0
						do c = 0, 3
							id = tnode_id(l + 1, 2*a - 1 + mod(c, 2), 2*b - 1 + c/2)
							r = max(r, sqrt(sum((quant%tnode_c(:, id, p) - cen)**2)) + quant%tnode_r(id, p))
						enddo
						quant%tnode_c(:, tnode_id(l, a, b), p) = cen
						quant%tnode_r(tnode_id(l, a, b), p) = r
					enddo
				enddo
			enddo
		enddo
	end subroutine tensor_setup_CHEB

	!**** range lo..hi of patch positions (along one dimension of a base patch) and quadtree level l of the group grp of a
	! per-dimension cluster tree (a complete binary tree whose level tlev holds single patches; deeper groups lie in one patch)
	subroutine tgroup_range(quant, grp, lo, hi, l)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		integer grp, lo, hi, l, g
		l = 0
		do while (2**(l + 1) <= grp)
			l = l + 1
		enddo
		g = grp
		do while (l > quant%tlev)
			g = g/2
			l = l - 1
		enddo
		lo = (g - 2**l)*2**(quant%tlev - l) + 1
		hi = lo + 2**(quant%tlev - l) - 1
	end subroutine tgroup_range

	!**** bounding sphere of the patches lo(1)..hi(1) x lo(2)..hi(2) of base patch p: their quadtree node if both dimensions
	! are at the same level, otherwise merged from the patches
	subroutine tcluster_sphere(quant, p, lo, hi, lv, cen, r)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		integer p, lo(2), hi(2), lv(2), is, js, id
		real(kind=8) cen(3), r
		if (lv(1) == lv(2)) then
			id = tnode_id(lv(1), (lo(1) - 1)/(hi(1) - lo(1) + 1) + 1, (lo(2) - 1)/(hi(2) - lo(2) + 1) + 1)
			cen = quant%tnode_c(:, id, p)
			r = quant%tnode_r(id, p)
		else
			cen = 0d0
			do js = lo(2), hi(2)
				do is = lo(1), hi(1)
					cen = cen + quant%tnode_c(:, tnode_id(quant%tlev, is, js), p)
				enddo
			enddo
			cen = cen/((hi(1) - lo(1) + 1)*(hi(2) - lo(2) + 1))
			r = 0d0
			do js = lo(2), hi(2)
				do is = lo(1), hi(1)
					id = tnode_id(quant%tlev, is, js)
					r = max(r, sqrt(sum((quant%tnode_c(:, id, p) - cen)**2)) + quant%tnode_r(id, p))
				enddo
			enddo
		endif
	end subroutine tcluster_sphere

	!**** near (0) or far (1) for the target groups gm (on base patch p) and the source groups gn (on base patch q): far if
	! the bounding spheres are farther apart than tnear_para*(r1+r2)/2 and no target in gm has a self/near pair with a patch in gn
	integer function tensor_nearfar_CHEB(quant, p, q, gm, gn)
		implicit none
		type(quant_ACOUSTIC_CHEB):: quant
		integer p, q, gm(2), gn(2), lo(2, 2), hi(2, 2), lv(2, 2), d, is, js, t, kk, sp
		real(kind=8) cen(3, 2), r(2)
		do d = 1, 2
			call tgroup_range(quant, gm(d), lo(d, 1), hi(d, 1), lv(d, 1))
			call tgroup_range(quant, gn(d), lo(d, 2), hi(d, 2), lv(d, 2))
		enddo
		call tcluster_sphere(quant, p, lo(:, 1), hi(:, 1), lv(:, 1), cen(:, 1), r(1))
		call tcluster_sphere(quant, q, lo(:, 2), hi(:, 2), lv(:, 2), cen(:, 2), r(2))
		tensor_nearfar_CHEB = 0
		if (sqrt(sum((cen(:, 1) - cen(:, 2))**2)) <= quant%tnear_para*(r(1) + r(2))/2d0) return
		do js = lo(2, 1), hi(2, 1)
			do is = lo(1, 1), hi(1, 1)
				t = (p - 1)*quant%split**2 + subpatch_pos(is, js, quant%split, quant%forest)
				do kk = quant%nsub_ptr(t), quant%nsub_ptr(t + 1) - 1
					sp = quant%nsub_q(kk)
					if (quant%sp_face(sp) == q .and. quant%sp_is(sp) >= lo(1, 2) .and. quant%sp_is(sp) <= hi(1, 2) &
						&	.and. quant%sp_js(sp) >= lo(2, 2) .and. quant%sp_js(sp) <= hi(2, 2)) return
				enddo
			enddo
		enddo
		tensor_nearfar_CHEB = 1
	end function tensor_nearfar_CHEB

	!**** user-defined subroutine to sample one entry of the H-TENSOR block of base patches (p,q) (ButterflyPACK's
	! FuncZmn_MD interface); m and n are grid indices (along s, along t) in natural order
	subroutine Zelem_tensor_CHEB(Ndim, m, n, value, quant)
		implicit none
		integer, INTENT(IN):: Ndim
		integer, INTENT(IN):: m(Ndim), n(Ndim)
		complex(kind=8) value
		class(*), pointer :: quant
		select TYPE (quant)
		type is (quant_pair_CHEB)
			value = Zentry_CHEB(quant%g, tnode_CHEB(quant%g, quant%p, m(1), m(2)), tnode_CHEB(quant%g, quant%q, n(1), n(2)))
		class default
			write (*, *) "unexpected type"
			stop
		end select
	end subroutine Zelem_tensor_CHEB

	!**** user-defined subroutine to sample a list of subtensors of the H-TENSOR block of base patches (p,q) (ButterflyPACK's
	! FuncZmnBlock_MD interface, used with option%elem_extract=2). allrows(:,d), allcols(:,d) hold the concatenated grid
	! indices along dimension d of all subtensors; each subtensor is returned as a column-major matrix whose rows and columns
	! are the target and source multi-indices, dimension 1 fastest
	subroutine Zelem_tensor_CHEB_block(Ndim, Ninter, allrows, allcols, alldat_loc, rowidx, colidx, pgidx, Npmap, pmaps, quant)
		implicit none
		integer, INTENT(IN):: Ndim
		class(*), pointer :: quant
		integer:: Ninter
		integer:: allrows(:, :), allcols(:, :)
		complex(kind=8), target:: alldat_loc(:)
		integer:: colidx(Ndim, Ninter), rowidx(Ndim, Ninter), pgidx(Ninter)
		integer:: Npmap, pmaps(Npmap, 3)
		integer nn, nthreads
		integer, allocatable:: roff(:, :), coff(:, :)
		integer(kind=8), allocatable:: voff(:)

		select TYPE (quant)
		type is (quant_pair_CHEB)
			allocate (roff(Ndim, Ninter), coff(Ndim, Ninter), voff(Ninter))
			roff(:, 1) = 0
			coff(:, 1) = 0
			voff(1) = 0
			do nn = 2, Ninter
				roff(:, nn) = roff(:, nn - 1) + rowidx(:, nn - 1)
				coff(:, nn) = coff(:, nn - 1) + colidx(:, nn - 1)
				voff(nn) = voff(nn - 1) + int(product(rowidx(:, nn - 1)), 8)*product(colidx(:, nn - 1))
			enddo
			nthreads = 1
#ifdef HAVE_OPENMP
			nthreads = omp_get_max_threads()
#endif
			if (Ninter >= nthreads) then ! many subtensors: one thread per subtensor
#ifdef HAVE_OPENMP
				!$omp parallel do default(shared) private(nn) schedule(dynamic,1)
#endif
				do nn = 1, Ninter
					call tensor_subtensor_CHEB(quant, rowidx(:, nn), colidx(:, nn), allrows(roff(1, nn) + 1:, 1), allrows(roff(2, nn) + 1:, 2), &
						&	allcols(coff(1, nn) + 1:, 1), allcols(coff(2, nn) + 1:, 2), alldat_loc(voff(nn) + 1:), 0)
				enddo
#ifdef HAVE_OPENMP
				!$omp end parallel do
#endif
			else ! few subtensors: the threads share the rows of each
				do nn = 1, Ninter
					call tensor_subtensor_CHEB(quant, rowidx(:, nn), colidx(:, nn), allrows(roff(1, nn) + 1:, 1), allrows(roff(2, nn) + 1:, 2), &
						&	allcols(coff(1, nn) + 1:, 1), allcols(coff(2, nn) + 1:, 2), alldat_loc(voff(nn) + 1:), 1)
				enddo
			endif
			deallocate (roff, coff, voff)
		class default
			write (*, *) "unexpected type"
			stop
		end select
	end subroutine Zelem_tensor_CHEB_block

	!**** one subtensor: rows r1 x r2 of base patch p against columns c1 x c2 of base patch q, as the submatrix of the nodes
	subroutine tensor_subtensor_CHEB(pq, nr, nc, r1, r2, c1, c2, dat, par)
		implicit none
		type(quant_pair_CHEB):: pq
		integer nr(2), nc(2), r1(*), r2(*), c1(*), c2(*), par, a, b
		complex(kind=8) dat(*)
		integer, allocatable:: rows(:), cols(:)
		allocate (rows(nr(1)*nr(2)), cols(nc(1)*nc(2)))
		do b = 1, nr(2)
			do a = 1, nr(1)
				rows(a + (b - 1)*nr(1)) = tnode_CHEB(pq%g, pq%p, r1(a), r2(b))
			enddo
		enddo
		do b = 1, nc(2)
			do a = 1, nc(1)
				cols(a + (b - 1)*nc(1)) = tnode_CHEB(pq%g, pq%q, c1(a), c2(b))
			enddo
		enddo
		call Zsubmatrix_CHEB(pq%g, nr(1)*nr(2), nc(1)*nc(2), rows, cols, dat, par)
		deallocate (rows, cols)
	end subroutine tensor_subtensor_CHEB

	!**** user-defined compressibility of the H-TENSOR block of base patches (p,q) (ButterflyPACK's FuncNearFar_MD interface,
	! option%nogeo=2): groupm, groupn are groups of the per-dimension cluster trees
	subroutine NearFar_tensor_CHEB(Ndim, groupm, groupn, val, quant)
		implicit none
		integer, INTENT(IN):: Ndim
		integer, INTENT(IN):: groupm(Ndim), groupn(Ndim)
		integer:: val
		class(*), pointer :: quant
		select TYPE (quant)
		type is (quant_pair_CHEB)
			val = tensor_nearfar_CHEB(quant%g, quant%p, quant%q, groupm, groupn)
		class default
			write (*, *) "unexpected type"
			stop
		end select
	end subroutine NearFar_tensor_CHEB

	!**** construct the nbase x nbase H-TENSOR blocks, and the node of every local entry
	subroutine tensor_build_CHEB(tf, quant, option, ptree)
		use z_BPACK_constr
		use z_BPACK_utilities
		implicit none
		type(tforest_CHEB), target:: tf
		type(quant_ACOUSTIC_CHEB), target:: quant
		type(z_Hoption), target:: option
		type(z_proctree), target:: ptree
		integer P, p1, q1, verb, a1, a2, i, ierr, l, lmax
		integer Ns_loc(2), idxs(2)
		integer, allocatable:: tree(:, :), perm(:, :)
		real(kind=8) t1, t2, tsum, mem(2), memg(2), memtot(2), tp, te, tesum
		integer rkmax, rk(0:63)

		P = quant%nbase
		tf%nface = P
		tf%Ns = quant%split*quant%N
		tf%option => option
		tf%ptree => ptree
		allocate (tf%bmat(P, P), tf%msh(2, P, P), tf%ker(P, P), tf%stats(P, P), tf%pair(P, P))
		allocate (tree(quant%split, 2), perm(maxval(tf%Ns), 2))
		verb = option%verbosity
		tsum = 0d0
		tesum = 0d0
		memtot = 0d0
		rk = 0
		if (ptree%MyID == Main_ID) then
			write (*, *) ''
			write (*, '(A,I3,A,I3,A,I6,A,I6,A)') ' tensor forest: ', P, ' x', P, ' H-TENSOR blocks on', tf%Ns(1), ' x', tf%Ns(2), ' parameter grids'
			write (*, '(A)') '   p   q   construction (s)   entry evaluation (s)   compressed (MB)   dense (MB)   max rank'
		endif
		do q1 = 1, P
			do p1 = 1, P
				tf%pair(p1, q1)%p = p1
				tf%pair(p1, q1)%q = q1
				tf%pair(p1, q1)%g => quant
				call z_InitStat(tf%stats(p1, q1))
				tf%ker(p1, q1)%QuantApp => tf%pair(p1, q1)
				tf%ker(p1, q1)%FuncZmn_MD => Zelem_tensor_CHEB
				tf%ker(p1, q1)%FuncZmnBlock_MD => Zelem_tensor_CHEB_block
				tf%ker(p1, q1)%FuncNearFar_MD => NearFar_tensor_CHEB
				if (p1 /= 1 .or. q1 /= 1) option%verbosity = -1 ! ButterflyPACK's output for the first block only
				tree = quant%N
				t1 = MPI_Wtime()
				call z_BPACK_MD_construction_Init(tf%Ns, 2, perm, Ns_loc, tf%bmat(p1, q1), option, tf%stats(p1, q1), tf%msh(:, p1, q1), &
					&	tf%ker(p1, q1), ptree, tree=tree)
				call z_BPACK_MD_construction_Element(2, tf%bmat(p1, q1), option, tf%stats(p1, q1), tf%msh(:, p1, q1), tf%ker(p1, q1), ptree)
				t2 = MPI_Wtime()
				option%verbosity = verb
				if (p1 == 1 .and. q1 == 1) then
					tf%Ns_loc = Ns_loc
				else
					call z_assert(all(Ns_loc == tf%Ns_loc), 'tensor forest: the H-TENSOR blocks have different vector layouts')
				endif
				call MPI_ALLREDUCE(t2 - t1, tp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
				call MPI_ALLREDUCE(tf%stats(p1, q1)%Time_Entry, te, 1, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
				tsum = tsum + tp
				tesum = tesum + te
				mem = (/tf%stats(p1, q1)%Mem_Comp_for, tf%stats(p1, q1)%Mem_Direct_for/)
				call MPI_ALLREDUCE(mem, memg, 2, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
				memtot = memtot + memg
				lmax = ubound(tf%stats(p1, q1)%rankmax_of_level_global, 1)
				rkmax = maxval(tf%stats(p1, q1)%rankmax_of_level_global)
				do l = 0, min(lmax, 63)
					rk(l) = max(rk(l), tf%stats(p1, q1)%rankmax_of_level_global(l))
				enddo
				if (ptree%MyID == Main_ID) write (*, '(2I4,F16.3,F20.3,2F18.3,I11)') p1, q1, tp, te, memg(1), memg(2), rkmax
			enddo
		enddo
		if (ptree%MyID == Main_ID) then
			write (*, '(A,F12.3,A,F12.3,A)') ' tensor forest construction time:', tsum, ' s, of which entry evaluation', tesum, ' s'
			write (*, '(A,F12.3,A,F12.3,A,F12.3,A)') ' tensor forest memory: ', memtot(1), ' MB compressed + ', memtot(2), ' MB dense = ', sum(memtot), ' MB'
			write (*, '(A,64I5)') ' tensor forest max rank per level:', rk(0:min(lmax, 63))
		endif
		tf%Ns_loc = Ns_loc
		tf%nloc = product(tf%Ns_loc)
		do l = 1, 2
			idxs(l) = tf%msh(l, 1, 1)%idxs
		enddo
		allocate (tf%gidx(tf%nloc, P))
		do p1 = 1, P
			do a2 = 1, tf%Ns_loc(2)
				do a1 = 1, tf%Ns_loc(1)
					i = a1 + (a2 - 1)*tf%Ns_loc(1)
					tf%gidx(i, p1) = tnode_CHEB(quant, p1, tf%msh(1, 1, 1)%new2old(idxs(1) + a1 - 1), tf%msh(2, 1, 1)%new2old(idxs(2) + a2 - 1))
				enddo
			enddo
		enddo
		deallocate (tree, perm)
	end subroutine tensor_build_CHEB

	!**** Vout = Z*Vin with the tensor forest; local vectors hold the local boxes of the base patches one after another
	subroutine tensor_mult_CHEB(tf, nvec, Vin, Vout)
		use z_BPACK_Solve_Mul
		implicit none
		type(tforest_CHEB):: tf
		integer nvec, p1, q1, nl, verb
		complex(kind=8) Vin(:, :), Vout(:, :)
		complex(kind=8), allocatable:: vq(:, :), vp(:, :)
		real(kind=8) t1
		t1 = MPI_Wtime()
		verb = tf%option%verbosity
		if (verb < 2) tf%option%verbosity = min(verb, 0) ! per-level timing of every block product only at verbosity >= 2
		nl = tf%nloc
		allocate (vq(max(nl, 1), nvec), vp(max(nl, 1), nvec))
		Vout(1:nl*tf%nface, 1:nvec) = 0d0
		do q1 = 1, tf%nface
			vq(1:nl, :) = Vin((q1 - 1)*nl + 1:q1*nl, 1:nvec)
			do p1 = 1, tf%nface
				call z_BPACK_MD_Mult(2, 'N', tf%Ns_loc, nvec, vq, vp, tf%bmat(p1, q1), tf%ptree, tf%option, tf%stats(p1, q1), tf%msh(:, p1, q1))
				Vout((p1 - 1)*nl + 1:p1*nl, 1:nvec) = Vout((p1 - 1)*nl + 1:p1*nl, 1:nvec) + vp(1:nl, :)
			enddo
		enddo
		deallocate (vq, vp)
		tf%option%verbosity = verb
		tf%tmvp = tf%tmvp + MPI_Wtime() - t1
		tf%nmvp = tf%nmvp + nvec
	end subroutine tensor_mult_CHEB

	!**** product timers of ButterflyPACK on this rank, summed over the face pairs: redistribution in, levels, redistribution out, final update,
	!**** within the levels the dense blocks and the stages of the tensor-butterfly blocks (init, right, all2all, middle, left, cleanup),
	!**** the stages of the cached-plan redistributions (setup/alloc/free, self copy, pack, MPI, unpack, total),
	!**** the wall time of the concurrent products with single-rank blocks, and the flops
	subroutine tensor_mvp_timers_CHEB(tf, t)
		implicit none
		type(tforest_CHEB):: tf
		real(kind=8) t(19)
		integer p1, q1
		t = 0d0
		do q1 = 1, tf%nface
			do p1 = 1, tf%nface
				associate (s => tf%stats(p1, q1))
					t = t + [s%Time_C_Mult_RedistIn, s%Time_C_Mult_Level, s%Time_C_Mult_RedistOut, s%Time_C_Mult_Final, s%Time_C_Mult_Full, &
						s%Time_C_Mult_Init, s%Time_C_Mult_Right, s%Time_C_Mult_All2All, s%Time_C_Mult_Middle, s%Time_C_Mult_Left, s%Time_C_Mult_Cleanup, &
						s%Time_C_Mult_RedistSetup + s%Time_C_Mult_RedistAlloc + s%Time_C_Mult_RedistFree, s%Time_C_Mult_RedistSelf, &
						s%Time_C_Mult_RedistPack, s%Time_C_Mult_RedistMPI, s%Time_C_Mult_RedistUnpack, s%Time_C_Mult_RedistTotal, &
						s%Time_C_Mult_Local, s%Flop_C_Mult]
				end associate
			enddo
		enddo
	end subroutine tensor_mvp_timers_CHEB

	!**** per-product breakdown of the tensor-forest product between two snapshots of tensor_mvp_timers_CHEB (max over ranks)
	subroutine tensor_mvp_profile_CHEB(tf, ptree, t0, t1, twall, nprod)
		implicit none
		type(tforest_CHEB):: tf
		type(z_proctree):: ptree
		real(kind=8) t0(19), t1(19), twall, tl(20), tg(20), fl
		integer nprod, p1, q1, ll, bb, nlev, nbf, nfull, nshared, cnt(4), cntg(4), ierr

		tl(1:19) = (t1 - t0)/max(1, nprod)
		tl(20) = twall/max(1, nprod)
		call MPI_ALLREDUCE(tl, tg, 20, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
		call MPI_ALLREDUCE(tl(19), fl, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
		nlev = 0
		nbf = 0
		nfull = 0
		nshared = 0
		do q1 = 1, tf%nface
			do p1 = 1, tf%nface
				associate (bp => tf%bmat(p1, q1)%h_mat_md%BP)
					nlev = nlev + (bp%Lplus - 1 + max(1, tf%option%htensor_mvp_level_batch) - 1)/max(1, tf%option%htensor_mvp_level_batch)
					do ll = 2, bp%Lplus
						do bb = 1, bp%LL(ll)%Nbound_loc
							if (bp%LL(ll)%matrices_block(bb)%style == 1) then
								nfull = nfull + 1
							else
								nbf = nbf + 1
							endif
							if (tf%ptree%pgrp(bp%LL(ll)%matrices_block(bb)%pgno)%nproc > 1) nshared = nshared + 1
						enddo
					enddo
				end associate
			enddo
		enddo
		cnt = [nlev, nbf, nfull, nshared]
		call MPI_ALLREDUCE(cnt, cntg, 4, MPI_INTEGER, MPI_MAX, ptree%Comm, ierr)
		if (ptree%MyID == Main_ID) then
			write (*, '(A,I6,A,I6,A,I6,A)') ' tensor forest product profile (s per product, max over ranks;', cntg(1), ' level calls,', &
				cntg(2), ' butterfly and', cntg(3), ' dense blocks on the busiest rank):'
			write (*, '(A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4)') '   wall', tg(20), '  redist in', tg(1), '  levels', tg(2), '  redist out', tg(3), '  final', tg(4)
			write (*, '(A,F10.4,A,I8,A,F10.3,A)') '   levels: single-rank blocks (concurrent, wall)', tg(18), '   blocks shared by several ranks:', cntg(4), &
				'   flops per product:', fl/1d9, ' GFlop (all ranks)'
			write (*, '(A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4)') '   levels (summed over threads): dense', tg(5), '  BF init', tg(6), '  right', tg(7), &
				'  all2all', tg(8), '  middle', tg(9), '  left', tg(10), '  cleanup', tg(11)
			write (*, '(A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4,A,F10.4)') '   redistribution (in + out): setup', tg(12), '  self copy', tg(13), &
				'  pack', tg(14), '  MPI', tg(15), '  unpack', tg(16), '  total', tg(17)
		endif
	end subroutine tensor_mvp_profile_CHEB

	!**** ButterflyPACK's HMatVec interface for GMRES: the product with the tensor forest
	subroutine tensor_blackbox_MVP_CHEB(trans, M, N, num_vect, Vin, Vout, ker)
		implicit none
		character trans
		integer, INTENT(IN):: M, N, num_vect
		complex(kind=8):: Vin(:, :), Vout(:, :)
		type(z_kernelquant):: ker
		select TYPE (tf => ker%QuantApp)
		type is (tforest_CHEB)
			call tensor_mult_CHEB(tf, num_vect, Vin, Vout)
		class default
			write (*, *) "unexpected type"
			stop
		end select
	end subroutine tensor_blackbox_MVP_CHEB

	!**** ButterflyPACK's HMatVec interface: no preconditioner
	subroutine identity_MVP_CHEB(trans, M, N, num_vect, Vin, Vout, ker)
		implicit none
		character trans
		integer, INTENT(IN):: M, N, num_vect
		complex(kind=8):: Vin(:, :), Vout(:, :)
		type(z_kernelquant):: ker
		Vout(1:M, 1:num_vect) = Vin(1:M, 1:num_vect)
	end subroutine identity_MVP_CHEB

	!**** accuracy of the tensor-forest product: Z*x for a random x against the uncompressed entries on a sample of rows
	subroutine tensor_mvp_check_CHEB(tf, quant, ptree)
		implicit none
		type(tforest_CHEB):: tf
		type(quant_ACOUSTIC_CHEB):: quant
		type(z_proctree):: ptree
		integer nl, nsamp, kk, ii, n, ierr
		real(kind=8), allocatable:: xr(:, :)
		complex(kind=8), allocatable:: x(:, :), y(:, :), xg(:)
		complex(kind=8) val
		real(kind=8) errs(2), errg(2)

		nl = tf%nloc*tf%nface
		allocate (xr(max(nl, 1), 2), x(max(nl, 1), 1), y(max(nl, 1), 1), xg(quant%Nunk))
		call random_number(xr)
		x(:, 1) = cmplx(xr(:, 1) - 0.5d0, xr(:, 2) - 0.5d0, kind=8)
		xg = 0d0
		do ii = 1, nl
			xg(tf%gidx(mod(ii - 1, tf%nloc) + 1, (ii - 1)/tf%nloc + 1)) = x(ii, 1)
		enddo
		call MPI_ALLREDUCE(MPI_IN_PLACE, xg, quant%Nunk, MPI_DOUBLE_COMPLEX, MPI_SUM, ptree%Comm, ierr)
		call tensor_mult_CHEB(tf, 1, x, y)
		errs = 0d0
		nsamp = min(nl, 8)
		do kk = 1, nsamp
			ii = 1 + ((kk - 1)*nl)/nsamp
			val = 0d0
#ifdef HAVE_OPENMP
			!$omp parallel do default(shared) private(n) reduction(+:val)
#endif
			do n = 1, quant%Nunk
				val = val + Zentry_CHEB(quant, tf%gidx(mod(ii - 1, tf%nloc) + 1, (ii - 1)/tf%nloc + 1), n)*xg(n)
			enddo
#ifdef HAVE_OPENMP
			!$omp end parallel do
#endif
			errs(1) = errs(1) + abs(y(ii, 1) - val)**2
			errs(2) = errs(2) + abs(val)**2
		enddo
		call MPI_ALLREDUCE(errs, errg, 2, MPI_DOUBLE_PRECISION, MPI_SUM, ptree%Comm, ierr)
		if (ptree%MyID == Main_ID) write (*, '(A,Es12.4)') ' tensor forest product vs. uncompressed entries (random vector, sampled rows):', sqrt(errg(1)/errg(2))
		deallocate (xr, x, y, xg)
	end subroutine tensor_mvp_check_CHEB

	!**** forward-map test of the tensor forest on the sphere
	subroutine tensor_forward_map_test_CHEB(tf, quant, ptree)
		implicit none
		type(tforest_CHEB):: tf
		type(quant_ACOUSTIC_CHEB):: quant
		type(z_proctree):: ptree
		integer nl, ii
		integer, allocatable:: gi(:)
		real(kind=8), allocatable:: Y(:)
		complex(kind=8), allocatable:: vin(:, :), vout(:, :)
		complex(kind=8) lam

		if (quant%geo /= CHEB_SPHERE .or. quant%eigtest <= 0) return
		call eigtest_Y_CHEB(quant, Y, lam)
		nl = tf%nloc*tf%nface
		allocate (gi(max(nl, 1)), vin(max(nl, 1), 1), vout(max(nl, 1), 1))
		do ii = 1, nl
			gi(ii) = tf%gidx(mod(ii - 1, tf%nloc) + 1, (ii - 1)/tf%nloc + 1)
			vin(ii, 1) = Y(gi(ii))
			if (quant%unknown == UNK_PSI) vin(ii, 1) = vin(ii, 1)*quant%e(gi(ii))
		enddo
		call tensor_mult_CHEB(tf, 1, vin, vout)
		call forward_map_report_CHEB(quant, ptree, nl, gi, vout(:, 1), Y, lam)
		deallocate (Y, gi, vin, vout)
	end subroutine tensor_forward_map_test_CHEB

	!**** solve for a plane wave by GMRES with the tensor-forest product (no preconditioner) and evaluate the far field
	subroutine tensor_solve_CHEB(tf, quant, option, ptree)
		use z_BPACK_Solve_Mul
		use z_BPACK_utilities
		implicit none
		type(tforest_CHEB), target:: tf
		type(quant_ACOUSTIC_CHEB):: quant
		type(z_Hoption):: option
		type(z_proctree):: ptree
		type(z_kernelquant):: kert
		type(z_Hstat):: stats_sol
		integer nl, ii, iter, restart, nmvp0
		integer, allocatable:: gi(:)
		real(kind=8) t1, t2, relres, tmvp0, tp0(19), tp1(19)
		complex(kind=8), allocatable:: x(:, :), b(:, :)

		nl = tf%nloc*tf%nface
		allocate (gi(max(nl, 1)), x(max(nl, 1), 1), b(max(nl, 1), 1))
		do ii = 1, nl
			gi(ii) = tf%gidx(mod(ii - 1, tf%nloc) + 1, (ii - 1)/tf%nloc + 1)
			b(ii, 1) = -exp(BPACK_junit*quant%wavenum*dot_product(quant%dinc, quant%xyz(:, gi(ii))))
		enddo
		x = 0d0
		call z_InitStat(stats_sol)
		kert%QuantApp => tf
		restart = 30
		if (quant%gmres_restart > 0) restart = quant%gmres_restart
		iter = 0
		relres = option%tol_itersol
		tmvp0 = tf%tmvp
		nmvp0 = tf%nmvp
		call tensor_mvp_timers_CHEB(tf, tp0)
		t1 = MPI_Wtime()
		call z_BPACK_Zgmres_usermatvec_precon(option%n_iter, nl, b, x, relres, iter, tensor_blackbox_MVP_CHEB, identity_MVP_CHEB, &
			&	ptree, option, stats_sol, kert, restart_in=restart)
		t2 = MPI_Wtime()
		call tensor_mvp_timers_CHEB(tf, tp1)
		if (ptree%MyID == Main_ID) then
			write (*, '(A,I5,A,I6,A,Es12.4)') ' GMRES(', restart, ') iterations:', iter, '   relative residual:', relres
			write (*, '(A,F12.3,A,F10.4,A)') ' tensor forest solve time:', t2 - t1, ' s,', (tf%tmvp - tmvp0)/max(1, tf%nmvp - nmvp0), ' s per product'
		endif
		call tensor_mvp_profile_CHEB(tf, ptree, tp0, tp1, tf%tmvp - tmvp0, tf%nmvp - nmvp0)
		call farfield_CHEB(quant, ptree, nl, gi, x(:, 1))
		call z_delete_Hstat(stats_sol)
		deallocate (gi, x, b)
	end subroutine tensor_solve_CHEB

	!**** delete the tensor forest
	subroutine tensor_delete_CHEB(tf)
		use z_BPACK_utilities
		implicit none
		type(tforest_CHEB):: tf
		integer p1, q1, d
		do q1 = 1, tf%nface
			do p1 = 1, tf%nface
				call z_BPACK_delete(tf%bmat(p1, q1))
				do d = 1, 2
					call z_delete_mesh(tf%msh(d, p1, q1))
				enddo
				call z_delete_kernelquant(tf%ker(p1, q1))
				call z_delete_Hstat(tf%stats(p1, q1))
			enddo
		enddo
		deallocate (tf%bmat, tf%msh, tf%ker, tf%stats, tf%pair)
		if (allocated(tf%gidx)) deallocate (tf%gidx)
	end subroutine tensor_delete_CHEB

	!**** the tensor-forest pipeline: construction, accuracy checks, iterative solve and far field (quant%tensor=1)
	subroutine tensor_run_CHEB(quant, option, ptree)
		use z_BPACK_utilities
		implicit none
		type(quant_ACOUSTIC_CHEB), target:: quant
		type(z_Hoption), target:: option
		type(z_proctree), target:: ptree
		type(tforest_CHEB), target:: tf

		option%format = HTENSOR
		option%nogeo = 2 ! admissibility by NearFar_tensor_CHEB: the geometric test of the HTENSOR format assumes a Cartesian grid
		option%xyzsort = NATURAL
		option%Nmin_leaf = quant%N ! a leaf of a per-dimension tree is one patch
		option%knn = 0
		option%trans_invariant = 0
		if (option%elem_extract == 0) option%elem_extract = 2 ! subtensor extraction: each near-field weight row once per subtensor
		quant%tnear_para = option%near_para
		call tensor_setup_CHEB(quant, ptree)
		call z_PrintOptions(option, ptree)
		call tensor_build_CHEB(tf, quant, option, ptree)
		call nearfield_report_CHEB(quant, ptree)
		call tensor_forward_map_test_CHEB(tf, quant, ptree)
		call tensor_mvp_check_CHEB(tf, quant, ptree)
		call tensor_solve_CHEB(tf, quant, option, ptree)
		call tensor_delete_CHEB(tf)
	end subroutine tensor_run_CHEB

end module Acoustic_SURF_MODULE_CHEB
