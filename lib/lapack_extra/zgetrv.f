      SUBROUTINE ZGETRV( TRANS, N, A, LDA, IPIV, B, INFO )
*
*     Adapted from ZGETRS, but using ZTRSV instead of ZTRSM.
*     This gains time compared to ZGETRS with NHRS=1.
*     November 4, 2011
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETRV solves a system of linear equations
*     A * x = b,  A**T * x = b,  or  A**H * x = b
*  with a general N-by-N matrix A using the LU factorization computed
*  by ZGETRF. x and b are vectors.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * x = b     (No transpose)
*          = 'T':  A**T * x = b  (Transpose)
*          = 'C':  A**H * x = b  (Conjugate transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by ZGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) COMPLEX*16 array, dimension >= N 
*          On entry, the right hand side vector b.
*          On exit, the solution vector b.
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
      EXTERNAL           XERBLA, ZLASWP, ZTRSV
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
         CALL XERBLA( 'ZGETRV', -INFO )
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
         CALL ZLASWP( 1, B, N, 1, N, IPIV, 1 )
*
*        Solve L*x = b, overwriting b with x.
*
         CALL ZTRSV( 'L', 'N', 'U', N, A, LDA, B, 1 )
*
*        Solve U*x = b, overwriting b with x.
*
         CALL ZTRSV( 'U', 'N', 'N', N, A, LDA, B, 1 )
      ELSE
*
*        Solve A**T * x = b  or A**H * x = b.
*
*        Solve U'*x = b, overwriting b with x.
*
         CALL ZTRSV( 'U', TRANS, 'N', N, A, LDA, B, 1 )
*
*        Solve L'*x = b, overwriting b with x.
*
         CALL ZTRSV( 'L', TRANS, 'U', N, A, LDA, B, 1 )
*
*        Apply row interchanges to the solution vectors.
*
         CALL ZLASWP( 1, B, N, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of ZGETRV
*
      END
