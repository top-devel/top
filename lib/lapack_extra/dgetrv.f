      SUBROUTINE DGETRV( TRANS, N, A, LDA, IPIV, B, INFO )
*
*     Adapted from DGETRS, but using DTRSV instead of DTRSM.
*     This gains time compared to DGETRS with NHRS=1.
*     November 4, 2011
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRV solves a system of linear equations
*     A * x = b  or  A' * x = b
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF. x and b are vectors.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * x = b  (No transpose)
*          = 'T':  A'* x = b  (Transpose)
*          = 'C':  A'* x = b  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array of dimension at least  N.
*          On entry, the right hand side vector b.
*          On exit, the solution vector x.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRV', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * x = b.
*
*        Apply row interchanges to the right hand sides.
*
         CALL DLASWP( 1, B, N, 1, N, IPIV, 1 )
*
*        Solve L*x = b, overwriting b with x.
*
         CALL DTRSV( 'L', 'N', 'U', N, A, LDA, B, 1 )
*
*        Solve U*x = b, overwriting b with x.
*
         CALL DTRSV( 'U', 'N', 'N', N, A, LDA, B, 1 )
      ELSE
*
*        Solve A' * x = b.
*
*        Solve U'*x = b, overwriting b with x.
*
         CALL DTRSV( 'U', 'T', 'N', N, A, LDA, B, 1 )
*
*        Solve L'*x = b, overwriting b with x.
*
         CALL DTRSV( 'L', 'T', 'U', N, A, LDA, B, 1 )
*
*        Apply row interchanges to the solution vectors.
*
         CALL DLASWP( 1, B, N, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of DGETRV
*
      END
