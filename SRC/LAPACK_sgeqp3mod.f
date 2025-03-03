*  =====================================================================
      subroutine SGEQP3mod( M, N, A, LDA, JPVT, TAU, WORK, LWORK, INFO,
     $                      RANK, RTOL, ATOL )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*  -- Modified by F-H Rouet, Lawrence Berkeley National Lab, August 2014
*
      IMPLICIT NONE
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N, RANK
*     ..
*     .. Array Arguments ..
      INTEGER            JPVT( * )
      REAL               A( LDA, * ), TAU( * ), WORK( * ), RTOL, ATOL
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            INB, INBMIN, IXOVER
      PARAMETER          ( INB = 1, INBMIN = 2, IXOVER = 3 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            FJB, IWS, J, JB, LWKOPT, MINMN, MINWS, NA, NB,
     $                   NBMIN, NFXD, NX, SM, SMINMN, SN, TOPBMN, C
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEQRF, SLAQP2, SLAQPS, SORMQR, SSWAP, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      REAL               SNRM2
      EXTERNAL           ILAENV, SNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
*
      IF( INFO.EQ.0 ) THEN
         MINMN = MIN( M, N )
         IF( MINMN.EQ.0 ) THEN
            IWS = 1
            LWKOPT = 1
         ELSE
            IWS = 3*N + 1
            NB = ILAENV( INB, 'SGEQRF', ' ', M, N, -1, -1 )
            LWKOPT = 2*N + ( N + 1 )*NB
         END IF
         WORK( 1 ) = LWKOPT
*
         IF( ( LWORK.LT.IWS ) .AND. .NOT.LQUERY ) THEN
            INFO = -8
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGEQP3', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( MINMN.EQ.0 ) THEN
         RETURN
      END IF
*
*     Move initial columns up front.
*
      NFXD = 1
      DO 10 J = 1, N
         IF( JPVT( J ).NE.0 ) THEN
            IF( J.NE.NFXD ) THEN
               CALL SSWAP( M, A( 1, J ), 1, A( 1, NFXD ), 1 )
               JPVT( J ) = JPVT( NFXD )
               JPVT( NFXD ) = J
            ELSE
               JPVT( J ) = J
            END IF
            NFXD = NFXD + 1
         ELSE
            JPVT( J ) = J
         END IF
   10 CONTINUE
      NFXD = NFXD - 1
*
*     Factorize fixed columns
*  =======================
*
*     Compute the QR factorization of fixed columns and update
*     remaining columns.
*
      IF( NFXD.GT.0 ) THEN
         NA = MIN( M, NFXD )
*CC      CALL SGEQR2( M, NA, A, LDA, TAU, WORK, INFO )
         CALL SGEQRF( M, NA, A, LDA, TAU, WORK, LWORK, INFO )
         IWS = MAX( IWS, INT( WORK( 1 ) ) )
         IF( NA.LT.N ) THEN
*CC         CALL SORM2R( 'Left', 'Transpose', M, N-NA, NA, A, LDA,
*CC  $                   TAU, A( 1, NA+1 ), LDA, WORK, INFO )
            CALL SORMQR( 'Left', 'Transpose', M, N-NA, NA, A, LDA, TAU,
     $                   A( 1, NA+1 ), LDA, WORK, LWORK, INFO )
            IWS = MAX( IWS, INT( WORK( 1 ) ) )
         END IF
      END IF
*
*     Factorize free columns
*  ======================
*
      RANK=0

      IF( NFXD.LT.MINMN ) THEN
*
         SM = M - NFXD
         SN = N - NFXD
         SMINMN = MINMN - NFXD
*
*        Determine the block size.
*
         NB = ILAENV( INB, 'SGEQRF', ' ', SM, SN, -1, -1 )
         NBMIN = 2
         NX = 0
*
         IF( ( NB.GT.1 ) .AND. ( NB.LT.SMINMN ) ) THEN
*
*           Determine when to cross over from blocked to unblocked code.
*
            NX = MAX( 0, ILAENV( IXOVER, 'SGEQRF', ' ', SM, SN, -1,
     $           -1 ) )
*
*
            IF( NX.LT.SMINMN ) THEN
*
*              Determine if workspace is large enough for blocked code.
*
               MINWS = 2*SN + ( SN+1 )*NB
               IWS = MAX( IWS, MINWS )
               IF( LWORK.LT.MINWS ) THEN
*
*                 Not enough workspace to use optimal NB: Reduce NB and
*                 determine the minimum value of NB.
*
                  NB = ( LWORK-2*SN ) / ( SN+1 )
                  NBMIN = MAX( 2, ILAENV( INBMIN, 'SGEQRF', ' ', SM, SN,
     $                    -1, -1 ) )
*
*
               END IF
            END IF
         END IF
*
*        Initialize partial column norms. The first N elements of work
*        store the exact column norms.
*
         DO 20 J = NFXD + 1, N
            WORK( J ) = SNRM2( SM, A( NFXD+1, J ), 1 )
            WORK( N+J ) = WORK( J )
   20    CONTINUE
*
         IF( ( NB.GE.NBMIN ) .AND. ( NB.LT.SMINMN ) .AND.
     $       ( NX.LT.SMINMN ) ) THEN
*
*           Use blocked code initially.
*
            J = NFXD + 1
*
*           Compute factorization: while loop.
*
*
            TOPBMN = MINMN - NX
   30       CONTINUE
            IF( J.LE.TOPBMN ) THEN
               JB = MIN( NB, TOPBMN-J+1 )
*
*              Factorize JB columns among columns J:N.
*
               CALL SLAQPS( M, N-J+1, J-1, JB, FJB, A( 1, J ), LDA,
     $                      JPVT( J ), TAU( J ), WORK( J ), WORK( N+J ),
     $                      WORK( 2*N+1 ), WORK( 2*N+JB+1 ), N-J+1 )
               DO C=J,J+FJB-1
                 IF(ABS(A(1,1))<ATOL) THEN
                   GOTO 99
                 ELSEIF(ABS(A(C,C))/ABS(A(1,1))<RTOL)THEN
                   GOTO 99
                 ELSE
                   RANK=RANK+1
                 ENDIF
               END DO
               J = J + FJB
               GO TO 30
            END IF
         ELSE
            J = NFXD + 1
         END IF
*
*        Use unblocked code to factor the last or only block.
*
*
         IF( J.LE.MINMN )
     $      CALL SLAQP2( M, N-J+1, J-1, A( 1, J ), LDA, JPVT( J ),
     $                   TAU( J ), WORK( J ), WORK( N+J ),
     $                   WORK( 2*N+1 ) )
            DO C=J,MINMN
               IF(ABS(A(1,1))<ATOL) THEN
                  GOTO 99
                ELSEIF(ABS(A(C,C))/ABS(A(1,1))<RTOL)THEN
                  GOTO 99
                ELSE
                RANK=RANK+1
              ENDIF
            END DO
      END IF
*
      WORK( 1 ) = IWS
   99 CONTINUE
      RETURN
*
*     End of SGEQP3
*
      end subroutine SGEQP3mod
