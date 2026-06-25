*> \brief \b CGETRFmod
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CGETRFmod + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgetrf.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgetrf.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgetrf.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CGETRFmod( M, N, A, LDA, IPIV, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX            A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CGETRFmod computes an LU factorization of a general M-by-N matrix A
*> using partial pivoting with row interchanges.
*>
*> The factorization has the form
*>    A = P * L * U
*> where P is a permutation matrix, L is lower triangular with unit
*> diagonal elements (lower trapezoidal if m > n), and U is upper
*> triangular (upper trapezoidal if m < n).
*>
*> This is the right-looking Level 3 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          On entry, the M-by-N matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (min(M,N))
*>          The pivot indices; for 1 <= i <= min(M,N), row i of the
*>          matrix was interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*>                has been completed, but the factor U is exactly
*>                singular, and division by zero will occur if it is used
*>                to solve a system of equations.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup getrf
*
*  =====================================================================

      SUBROUTINE cgetrfmod( M, N, A, LDA, THRESH, IPIV, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * )
      REAL               THRESH
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX            ONE
      parameter( one = ( 1.0e+0, 0.0e+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           cgemm, cgetrf2mod, claswp, ctrsm,
     $                   xerbla
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ilaenv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      IF( m.LT.0 ) THEN
         info = -1
      ELSE IF( n.LT.0 ) THEN
         info = -2
      ELSE IF( lda.LT.max( 1, m ) ) THEN
         info = -4
      END IF
      IF( info.NE.0 ) THEN
         CALL xerbla( 'CGETRFmod', -info )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      nb = ilaenv( 1, 'CGETRFmod', ' ', m, n, -1, -1 )
      IF( nb.LE.1 .OR. nb.GE.min( m, n ) ) THEN
*
*        Use unblocked code.
*
         CALL cgetrf2mod( m, n, a, lda, THRESH, ipiv, info )
      ELSE
*
*        Use blocked code.
*
         DO 20 j = 1, min( m, n ), nb
            jb = min( min( m, n )-j+1, nb )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL cgetrf2mod( m-j+1, jb, a( j, j ), lda, THRESH, ipiv( j ),
     $                    iinfo )
*
*           Adjust INFO and the pivot indices.
*
            IF( info.EQ.0 .AND. iinfo.GT.0 )
     $         info = iinfo + j - 1
            DO 10 i = j, min( m, j+jb-1 )
               ipiv( i ) = j - 1 + ipiv( i )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL claswp( j-1, a, lda, j, j+jb-1, ipiv, 1 )
*
            IF( j+jb.LE.n ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL claswp( n-j-jb+1, a( 1, j+jb ), lda, j, j+jb-1,
     $                      ipiv, 1 )
*
*              Compute block row of U.
*
               CALL ctrsm( 'Left', 'Lower', 'No transpose', 'Unit',
     $                     jb,
     $                     n-j-jb+1, one, a( j, j ), lda, a( j, j+jb ),
     $                     lda )
               IF( j+jb.LE.m ) THEN
*
*                 Update trailing submatrix.
*
                  CALL cgemm( 'No transpose', 'No transpose',
     $                        m-j-jb+1,
     $                        n-j-jb+1, jb, -one, a( j+jb, j ), lda,
     $                        a( j, j+jb ), lda, one, a( j+jb, j+jb ),
     $                        lda )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of CGETRFmod
*

      END