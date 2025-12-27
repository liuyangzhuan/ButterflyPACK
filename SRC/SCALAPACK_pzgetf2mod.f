      SUBROUTINE pzgetf2mod( M, N, A, IA, JA, DESCA, THRESH, IPIV, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            DESCA( * ), IPIV( * )
      COMPLEX*16         A( * )
      DOUBLE PRECISION   THRESH
*     ..
*
*  Purpose
*  =======
*
*  pzgetf2mod computes an LU factorization of a general M-by-N
*  distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) using
*  partial pivoting with row interchanges.
*
*  The factorization has the form sub( A ) = P * L * U, where P is a
*  permutation matrix, L is lower triangular with unit diagonal
*  elements (lower trapezoidal if m > n), and U is upper triangular
*  (upper trapezoidal if m < n).
*
*  This is the right-looking Parallel Level 2 BLAS version of the
*  algorithm.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  This routine requires N <= NB_A-MOD(JA-1, NB_A) and square block
*  decomposition ( MB_A = NB_A ).
*
*  Arguments
*  =========
*
*  M       (global input) INTEGER
*          The number of rows to be operated on, i.e. the number of rows
*          of the distributed submatrix sub( A ). M >= 0.
*
*  N       (global input) INTEGER
*          The number of columns to be operated on, i.e. the number of
*          columns of the distributed submatrix sub( A ).
*          NB_A-MOD(JA-1, NB_A) >= N >= 0.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the M-by-N
*          distributed matrix sub( A ). On exit, this array contains
*          the local pieces of the factors L and U from the factoriza-
*          tion sub( A ) = P*L*U; the unit diagonal elements of L are
*          not stored.
*
*  IA      (global input) INTEGER
*          The row index in the global array A indicating the first
*          row of sub( A ).
*
*  JA      (global input) INTEGER
*          The column index in the global array A indicating the
*          first column of sub( A ).
*
*  DESCA   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix A.
*
*  IPIV    (local output) INTEGER array, dimension ( LOCr(M_A)+MB_A )
*          This array contains the pivoting information.
*          IPIV(i) -> The global row local row i was swapped with.
*          This array is tied to the distributed matrix A.
*
*  INFO    (local output) INTEGER
*          = 0:  successful exit
*          < 0:  If the i-th argument is an array and the j-entry had
*                an illegal value, then INFO = -(i*100+j), if the i-th
*                argument is a scalar and had an illegal value, then
*                INFO = -i.
*          > 0:  If INFO = K, U(IA+K-1,JA+K-1) is exactly zero.
*                The factorization has been completed, but the factor U
*                is exactly singular, and division by zero will occur if
*                it is used to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      parameter( block_cyclic_2d = 1, dlen_ = 9, dtype_ = 1,
     $                     ctxt_ = 2, m_ = 3, n_ = 4, mb_ = 5, nb_ = 6,
     $                     rsrc_ = 7, csrc_ = 8, lld_ = 9 )
      COMPLEX*16         ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          ROWBTOP
      INTEGER            I, IACOL, IAROW, ICOFF, ICTXT, IIA, IROFF, J,
     $                   JJA, MN, MYCOL, MYROW, NPCOL, NPROW
      COMPLEX*16         GMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           blacs_abort, blacs_gridinfo, chk1mat, igebr2d,
     $                   igebs2d, infog2l, pb_topget, pxerbla, pzamax,
     $                   pzgeru, pzscal, pzswap, pzelset
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          min, mod
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters.
*
      ictxt = desca( ctxt_ )
      CALL blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )
*
*     Test the input parameters.
*
      info = 0
      IF( nprow.EQ.-1 ) THEN
         info = -(600+ctxt_)
      ELSE
         CALL chk1mat( m, 1, n, 2, ia, ja, desca, 6, info )
         IF( info.EQ.0 ) THEN
            iroff = mod( ia-1, desca( mb_ ) )
            icoff = mod( ja-1, desca( nb_ ) )
            IF( n+icoff.GT.desca( nb_ ) ) THEN
               info = -2
            ELSE IF( iroff.NE.0 ) THEN
               info = -4
            ELSE IF( icoff.NE.0 ) THEN
               info = -5
            ELSE IF( desca( mb_ ).NE.desca( nb_ ) ) THEN
               info = -(600+nb_)
            END IF
         END IF
      END IF
*
      IF( info.NE.0 ) THEN
         CALL pxerbla( ictxt, 'pzgetf2mod', -info )
         CALL blacs_abort( ictxt, 1 )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 )
     $   RETURN
*
      mn = min( m, n )
      CALL infog2l( ia, ja, desca, nprow, npcol, myrow, mycol, iia, jja,
     $              iarow, iacol )
      CALL pb_topget( ictxt, 'Broadcast', 'Rowwise', rowbtop )
*
      IF( mycol.EQ.iacol ) THEN
         DO 10 j = ja, ja+mn-1
            i = ia + j - ja
*
*           Find pivot and test for singularity.
*
            CALL pzamax( m-j+ja, gmax, ipiv( iia+j-ja ), a, i, j,
     $                   desca, 1 )

*           Tiny pivot replacement (preserve sign)
            IF ( ABS( gmax ) .LT. THRESH ) THEN
               IF ( gmax .EQ. ZERO ) THEN
                  gmax = THRESH
               ELSE
                  gmax = gmax/abs(gmax) * THRESH
               END IF

*              Make sure the distributed matrix pivot entry matches gmax
*              Pivot entry is at global (I,J) after the row swap step logic.
               CALL pzelset( a, i, j, desca, gmax )

*              (optional) count tiny pivots
*              ntiny = ntiny + 1
            END IF

            IF( gmax.NE.zero ) THEN
*
*              Apply the row interchanges to columns JA:JA+N-1
*
               CALL pzswap( n, a, i, ja, desca, desca( m_ ), a,
     $                      ipiv( iia+j-ja ), ja, desca, desca( m_ ) )
*
*              Compute elements I+1:IA+M-1 of J-th column.
*
               IF( j-ja+1.LT.m )
     $            CALL pzscal( m-j+ja-1, one / gmax, a, i+1, j,
     $                         desca, 1 )
            ELSE IF( info.EQ.0 ) THEN
               info = j - ja + 1
            END IF
*
*           Update trailing submatrix
*
            IF( j-ja+1.LT.mn ) THEN
               CALL pzgeru( m-j+ja-1, n-j+ja-1, -one, a, i+1, j, desca,
     $                      1, a, i, j+1, desca, desca( m_ ), a, i+1,
     $                      j+1, desca )
            END IF
   10    CONTINUE
*
         CALL igebs2d( ictxt, 'Rowwise', rowbtop, mn, 1, ipiv( iia ),
     $                 mn )
*
      ELSE
*
         CALL igebr2d( ictxt, 'Rowwise', rowbtop, mn, 1, ipiv( iia ),
     $                 mn, myrow, iacol )
*
      END IF
*
      RETURN
*
*     End of pzgetf2mod
*
      END