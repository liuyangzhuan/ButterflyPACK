      SUBROUTINE pzgetrfmod( M, N, A, IA, JA, DESCA, THRESH, IPIV, INFO )
*
*  -- ScaLAPACK routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 25, 2001
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
*  pzgetrfmod computes an LU factorization of a general M-by-N distributed
*  matrix sub( A ) = (IA:IA+M-1,JA:JA+N-1) using partial pivoting with
*  row interchanges.
*
*  The factorization has the form sub( A ) = P * L * U, where P is a
*  permutation matrix, L is lower triangular with unit diagonal ele-
*  ments (lower trapezoidal if m > n), and U is upper triangular
*  (upper trapezoidal if m < n). L and U are stored in sub( A ).
*
*  This is the right-looking Parallel Level 3 BLAS version of the
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
*  This routine requires square block decomposition ( MB_A = NB_A ).
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
*          columns of the distributed submatrix sub( A ). N >= 0.
*
*  A       (local input/local output) COMPLEX*16 pointer into the
*          local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
*          On entry, this array contains the local pieces of the M-by-N
*          distributed matrix sub( A ) to be factored. On exit, this
*          array contains the local pieces of the factors L and U from
*          the factorization sub( A ) = P*L*U; the unit diagonal ele-
*          ments of L are not stored.
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
*  INFO    (global output) INTEGER
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
      COMPLEX*16         ONE
      parameter( one = 1.0d+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          COLBTOP, COLCTOP, ROWBTOP
      INTEGER            I, ICOFF, ICTXT, IINFO, IN, IROFF, J, JB, JN,
     $                   MN, MYCOL, MYROW, NPCOL, NPROW
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM1( 1 ), IDUM2( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           blacs_gridinfo, chk1mat, igamn2d, pchk1mat,
     $                   pb_topget, pb_topset, pxerbla, pzgemm, pzgetf2mod,
     $                   pzlaswp, pztrsm
*     ..
*     .. External Functions ..
      INTEGER            ICEIL
      EXTERNAL           iceil
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          min, mod
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters
*
      ictxt = desca( ctxt_ )
      CALL blacs_gridinfo( ictxt, nprow, npcol, myrow, mycol )
*
*     Test the input parameters
*
      info = 0
      IF( nprow.EQ.-1 ) THEN
         info = -(600+ctxt_)
      ELSE
         CALL chk1mat( m, 1, n, 2, ia, ja, desca, 6, info )
         IF( info.EQ.0 ) THEN
            iroff = mod( ia-1, desca( mb_ ) )
            icoff = mod( ja-1, desca( nb_ ) )
            IF( iroff.NE.0 ) THEN
               info = -4
            ELSE IF( icoff.NE.0 ) THEN
               info = -5
            ELSE IF( desca( mb_ ).NE.desca( nb_ ) ) THEN
               info = -(600+nb_)
            END IF
         END IF
         CALL pchk1mat( m, 1, n, 2, ia, ja, desca, 6, 0, idum1,
     $                  idum2, info )
      END IF
*
      IF( info.NE.0 ) THEN
         CALL pxerbla( ictxt, 'pzgetrfmod', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( desca( m_ ).EQ.1 ) THEN
         ipiv( 1 ) = 1
         RETURN
      ELSE IF( m.EQ.0 .OR. n.EQ.0 ) THEN
         RETURN
      END IF
*
*     Split-ring topology for the communication along process rows
*
      CALL pb_topget( ictxt, 'Broadcast', 'Rowwise', rowbtop )
      CALL pb_topget( ictxt, 'Broadcast', 'Columnwise', colbtop )
      CALL pb_topget( ictxt, 'Combine', 'Columnwise', colctop )
      CALL pb_topset( ictxt, 'Broadcast', 'Rowwise', 'S-ring' )
      CALL pb_topset( ictxt, 'Broadcast', 'Columnwise', ' ' )
      CALL pb_topset( ictxt, 'Combine', 'Columnwise', ' ' )
*
*     Handle the first block of columns separately
*
      mn = min( m, n )
      in = min( iceil( ia, desca( mb_ ) )*desca( mb_ ), ia+m-1 )
      jn = min( iceil( ja, desca( nb_ ) )*desca( nb_ ), ja+mn-1 )
      jb = jn - ja + 1
*
*     Factor diagonal and subdiagonal blocks and test for exact
*     singularity.
*
      CALL pzgetf2mod( m, jb, a, ia, ja, desca, THRESH, ipiv, info )
*
      IF( jb+1.LE.n ) THEN
*
*        Apply interchanges to columns JN+1:JA+N-1.
*
         CALL pzlaswp( 'Forward', 'Rows', n-jb, a, ia, jn+1, desca,
     $                 ia, in, ipiv )
*
*        Compute block row of U.
*
         CALL pztrsm( 'Left', 'Lower', 'No transpose', 'Unit', jb,
     $                n-jb, one, a, ia, ja, desca, a, ia, jn+1, desca )
*
         IF( jb+1.LE.m ) THEN
*
*           Update trailing submatrix.
*
            CALL pzgemm( 'No transpose', 'No transpose', m-jb, n-jb, jb,
     $                   -one, a, in+1, ja, desca, a, ia, jn+1, desca,
     $                   one, a, in+1, jn+1, desca )
*
         END IF
      END IF
*
*     Loop over the remaining blocks of columns.
*
      DO 10 j = jn+1, ja+mn-1, desca( nb_ )
         jb = min( mn-j+ja, desca( nb_ ) )
         i = ia + j - ja
*
*        Factor diagonal and subdiagonal blocks and test for exact
*        singularity.
*
         CALL pzgetf2mod( m-j+ja, jb, a, i, j, desca, THRESH, ipiv, iinfo )
*
         IF( info.EQ.0 .AND. iinfo.GT.0 )
     $      info = iinfo + j - ja
*
*        Apply interchanges to columns JA:J-JA.
*
         CALL pzlaswp( 'Forward', 'Rowwise', j-ja, a, ia, ja, desca,
     $                 i, i+jb-1, ipiv )
*
         IF( j-ja+jb+1.LE.n ) THEN
*
*           Apply interchanges to columns J+JB:JA+N-1.
*
            CALL pzlaswp( 'Forward', 'Rowwise', n-j-jb+ja, a, ia, j+jb,
     $                    desca, i, i+jb-1, ipiv )
*
*           Compute block row of U.
*
            CALL pztrsm( 'Left', 'Lower', 'No transpose', 'Unit', jb,
     $                   n-j-jb+ja, one, a, i, j, desca, a, i, j+jb,
     $                   desca )
*
            IF( j-ja+jb+1.LE.m ) THEN
*
*              Update trailing submatrix.
*
               CALL pzgemm( 'No transpose', 'No transpose', m-j-jb+ja,
     $                      n-j-jb+ja, jb, -one, a, i+jb, j, desca, a,
     $                      i, j+jb, desca, one, a, i+jb, j+jb, desca )
*
            END IF
         END IF
*
   10 CONTINUE
*
      IF( info.EQ.0 )
     $   info = mn + 1
      CALL igamn2d( ictxt, 'Rowwise', ' ', 1, 1, info, 1, idum1, idum2,
     $              -1, -1, mycol )
      IF( info.EQ.mn+1 )
     $   info = 0
*
      CALL pb_topset( ictxt, 'Broadcast', 'Rowwise', rowbtop )
      CALL pb_topset( ictxt, 'Broadcast', 'Columnwise', colbtop )
      CALL pb_topset( ictxt, 'Combine', 'Columnwise', colctop )
*
      RETURN
*
*     End of pzgetrfmod
*
      END