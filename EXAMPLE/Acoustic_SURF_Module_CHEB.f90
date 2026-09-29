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
		integer(kind=8):: nweights_cat(3) = 0 ! ... requested by: 1 whole leaf-patch blocks, 2 other submatrices, 3 single entries
		integer, allocatable:: paircount(:) ! diagnostic: number of on-the-fly computations of each near pair within whole leaf-patch blocks

		integer:: eigtest = 5 ! degree of the spherical-harmonic forward-map test (sphere only), 0: skip
		integer:: nfar = 181 ! number of far-field directions
		integer:: gmres_restart = 0 ! >0: GMRES restart length when precon/=DIRECT (0: ButterflyPACK's default solve path, GMRES(30))
		integer:: faceorder = 0 ! order of the six faces of the sphere/cube (the forest's base patches): 0: +x,-x,+y,-y,+z,-z; 1: sweep +z,+x,+y,-x,-y,-z
		integer:: printstruct = 0 ! 1: write every block's rank before (fort.100000+rank) and after (fort.200000+rank) the factorization, and the cluster geometry (clusters.txt)
	end type quant_ACOUSTIC_CHEB

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
		integer Nloc, ii, gi, n, nsamp, kk, ierr
		real(kind=8) ehat(3), errs(4), errg(4)
		real(kind=8), allocatable:: Y(:)
		complex(kind=8), allocatable:: vin(:, :), vout(:, :)
		complex(kind=8) lam, val

		if (quant%geo /= CHEB_SPHERE .or. quant%eigtest <= 0) return
		ehat = (/0.3d0, -0.5d0, 0.8d0/)
		ehat = ehat/sqrt(sum(ehat**2))
		lam = sphere_eigenvalue(quant%wavenum, quant%radius, quant%eigtest)
		allocate (Y(quant%Nunk))
		do n = 1, quant%Nunk
			Y(n) = legendreP(quant%eigtest, dot_product(ehat, quant%xyz(:, n))/sqrt(sum(quant%xyz(:, n)**2)))
		enddo

		! compressed operator
		Nloc = msh%idxe - msh%idxs + 1
		allocate (vin(Nloc, 1), vout(Nloc, 1))
		do ii = 1, Nloc
			gi = msh%new2old(msh%idxs + ii - 1)
			vin(ii, 1) = Y(gi)
			if (quant%unknown == UNK_PSI) vin(ii, 1) = vin(ii, 1)*quant%e(gi)
		enddo
		call z_BPACK_Mult('N', Nloc, 1, vin, vout, bmat, ptree, option, stats, use_blockcopy=0) ! forward blocks, not their post-factorization copy
		errs = 0d0
		do ii = 1, Nloc
			gi = msh%new2old(msh%idxs + ii - 1)
			errs(1) = max(errs(1), abs(vout(ii, 1) - lam*Y(gi)))
			errs(2) = max(errs(2), abs(lam*Y(gi)))
		enddo

		! uncompressed entries on a sample of rows
		nsamp = min(Nloc, 16)
		do kk = 1, nsamp
			ii = 1 + ((kk - 1)*Nloc)/nsamp
			gi = msh%new2old(msh%idxs + ii - 1)
			val = 0d0
#ifdef HAVE_OPENMP
			!$omp parallel do default(shared) private(n) reduction(+:val)
#endif
			do n = 1, quant%Nunk
				if (quant%unknown == UNK_PSI) then
					val = val + Zentry_CHEB(quant, gi, n)*Y(n)*quant%e(n)
				else
					val = val + Zentry_CHEB(quant, gi, n)*Y(n)
				endif
			enddo
#ifdef HAVE_OPENMP
			!$omp end parallel do
#endif
			errs(3) = max(errs(3), abs(val - lam*Y(gi)))
			errs(4) = max(errs(4), abs(lam*Y(gi)))
		enddo
		call MPI_ALLREDUCE(errs, errg, 4, MPI_DOUBLE_PRECISION, MPI_MAX, ptree%Comm, ierr)
		if (ptree%MyID == Main_ID) then
			write (*, *) ''
			write (*, '(A,I3,A,2Es14.6)') ' forward-map test with Y of degree', quant%eigtest, ', exact eigenvalue:', dble(lam), aimag(lam)
			write (*, '(A,Es12.4)') '   max error, uncompressed entries (sampled rows): ', errg(3)/errg(4)
			write (*, '(A,Es12.4)') '   max error, compressed operator (all rows):      ', errg(1)/errg(2)
			write (*, *) ''
		endif
		deallocate (Y, vin, vout)
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
		integer Nloc, ii, gi, m, ierr, iter
		real(kind=8) e1(3), xh(3), gam, t1, t2, errmax, fmax, relres
		complex(kind=8), allocatable:: x(:, :), b(:, :), ffloc(:), ff(:)
		complex(kind=8) coef, umie
		real(kind=8) k

		k = quant%wavenum
		Nloc = msh%idxe - msh%idxs + 1
		allocate (x(Nloc, 1), b(Nloc, 1))
		do ii = 1, Nloc
			gi = msh%new2old(msh%idxs + ii - 1)
			b(ii, 1) = -exp(BPACK_junit*k*dot_product(quant%dinc, quant%xyz(:, gi)))
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
				gi = msh%new2old(msh%idxs + ii - 1)
				if (quant%closed) then
					coef = -BPACK_junit*k*(dot_product(xh, quant%nrm(:, gi)) + 1d0)
				else
					coef = 1d0
				endif
				coef = coef*exp(-BPACK_junit*k*dot_product(xh, quant%xyz(:, gi)))*quant%Jw(gi)/(4d0*BPACK_pi)
				if (quant%unknown == UNK_PHI) coef = coef*quant%e(gi) ! psi = e*phi
				ffloc(m) = ffloc(m) + coef*x(ii, 1)
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
		deallocate (x, b, ffloc, ff)
	end subroutine solve_scattering_CHEB

end module Acoustic_SURF_MODULE_CHEB
