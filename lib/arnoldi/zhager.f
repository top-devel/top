      SUBROUTINE ZHAGER (N, X, ISGN, EST, KASE, REVCOM)
      INTEGER N, ISGN(N), KASE, REVCOM
      DOUBLE PRECISION EST
      COMPLEX *16 X(N,4)
C
C     ZHAGER ESTIMATES THE 1-NORM OF A SQUARE, COMPLEX MATRIX  A.
C     REVERSE COMMUNICATION IS USED FOR EVALUATING
C     MATRIX-VECTOR PRODUCTS. 
C
C     ON ENTRY
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX.  N .GE. 1.
C
C        ISGN    INTEGER(N)
C                USED AS WORKSPACE.
C
C        KASE    INTEGER
C                = 0.
C
C     ON INTERMEDIATE RETURNS 
C
C        KASE    = 1 OR 2.
C
C        X       COMPLEX(N,4)
C                MUST BE OVERWRITTEN BY 
C
C                    X(:,3)=A*X(:,1),             IF KASE=1, 
C                    X(:,4)=TRANSPOSE(A)*X(:,3),  IF KASE=2, 
C
C                AND SONEST MUST BE RE-CALLED, WITH ALL THE OTHER
C                PARAMETERS UNCHANGED.
C
C     ON FINAL RETURN
C
C        KASE    = 3.
C
C        EST     REAL
C                CONTAINS AN ESTIMATE (A LOWER BOUND) FOR NORM(A).
C
C     THIS VERSION DATED MARCH 16, 1988.
C     NICK HIGHAM, UNIVERSITY OF MANCHESTER.
C
C     REFERENCE
C     N.J. HIGHAM (1987) FORTRAN CODES FOR ESTIMATING
C     THE ONE-NORM OF A REAL OR COMPLEX MATRIX, WITH APPLICATIONS
C     TO CONDITION  ESTIMATION, NUMERICAL ANALYSIS REPORT NO. 135,
C     UNIVERSITY OF MANCHESTER, MANCHESTER M13 9PL, ENGLAND.
C
C     SUBROUTINES AND FUNCTIONS
C     BLAS     IDAMAX, DASUM, DCOPY
C     GENERIC  DABS, DNINT, DREAL, DSIGN
C
      INTEGER ITMAX
      PARAMETER (ITMAX = 5)
      DOUBLE PRECISION ZERO, ONE, TWO
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      COMPLEX *16 CZERO, CONE, CTWO
      PARAMETER (CZERO=(0.0D0,0.D0),CONE =(1.0D0,0.d0),
     &                CTWO = (2.0D0,0.D0))

C
C     INTERNAL VARIABLES
      INTEGER I, ITER, J, JLAST, JUMP, ICMAX1
      DOUBLE PRECISION ALTSGN, ESTOLD, TEMP, DCSUM1
C
      SAVE
C
      IF (KASE .EQ. 0) THEN
         DO 10,I = 1,N
            X(I,1) =( CONE/DBLE(N) )	
   10    CONTINUE
         KASE = 1
         REVCOM = 1
         JUMP = 1
         RETURN
      ENDIF
C
      GOTO (100, 200, 300, 400, 500) JUMP
C
C     ................ ENTRY   (JUMP = 1)
C     FIRST ITERATION.  X(:,3) HAS BEEN OVERWRITTEN BY A*X(:,1).
C
  100 CONTINUE
      IF (N.EQ.1) THEN
	EST = ABS(X(1,3))
c ..    QUIT
 	GOTO 510
      ENDIF
      EST = DCSUM1(N,X(1,3),1)
C
      DO 110,I = 1,N
      IF ( ABS(X(I,3)).NE.ZERO) THEN 
	X(I,3)=X(I,3)/DCMPLX(ABS(X(I,3)),0.D0)
      ELSE
	X(I,3)=CONE
      ENDIF
  110 CONTINUE
      KASE = 2
      REVCOM = 2
      JUMP = 2
      RETURN
C
C     ................ ENTRY   (JUMP = 2)
C     FIRST ITERATION.  X(:,3) HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X(:,2).
C
  200 CONTINUE
      J = ICMAX1(N,X(1,4),1)
      ITER = 2
C
C     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
C
  220 CONTINUE
      DO 230,I = 1,N
         X(I,1) = CZERO 
  230 CONTINUE
      X(J,1) = CONE
      KASE = 1
      REVCOM = 1
      JUMP = 3
      RETURN
C
C     ................ ENTRY   (JUMP = 3)
C     X HAS BEEN OVERWRITTEN BY A*X.
C
  300 CONTINUE
      CALL ZCOPY(N,X(1,3),1,X(1,1),1)
      ESTOLD = EST
      EST = DCSUM1(N,X(1,1),1)
  320 CONTINUE
C     TEST FOR CYCLING.                                C
      IF (EST .LE. ESTOLD) GOTO 410                   
      DO 330,I = 1,N
      IF (ABS(X(I,3)) .NE. ZERO) THEN
	X(I,3)=X(I,3)/DCMPLX(ABS(X(I,3)),0.D0)
      ELSE
	X(I,3)=CONE
      ENDIF
  330 CONTINUE
      KASE = 2
      REVCOM = 2
      JUMP = 4
      RETURN
C
C     ................ ENTRY   (JUMP = 4)
C     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
C
  400 CONTINUE
      JLAST = J
      J = ICMAX1(N,X(1,4),1)
      IF (   (  DREAL(X(JLAST,4)) .NE. ABS(DREAL(X(J,4)))  ) .AND.
     +       (ITER .LT. ITMAX)   ) THEN
         ITER = ITER + 1
         GOTO 220
      ENDIF
C
C     ITERATION COMPLETE.  FINAL STAGE. 
C
  410 CONTINUE
      ALTSGN = ONE
      DO 420,I = 1,N
         X(I,1) = DCMPLX(ALTSGN * (ONE + DBLE(I-1)/DBLE(N-1)),0.D0)
         ALTSGN = -ALTSGN
  420 CONTINUE
      KASE = 1
      REVCOM = 1
      JUMP = 5
      RETURN
C
C     ................ ENTRY   (JUMP = 5)
C     X HAS BEEN OVERWRITTEN BY A*X.
C
  500 CONTINUE
      TEMP = TWO*DCSUM1(N,X(1,3),1)/DBLE(3*N) 
      IF (TEMP. GT. EST) THEN 
         CALL ZCOPY(N,X(1,3),1,X(1,1),1)
         EST = TEMP 
      ENDIF
C
510   CONTINUE
      KASE = 3
      RETURN
C
      END




      INTEGER FUNCTION ICMAX1(N,CX,INCX)
C
C     FINDS THE INDEX OF ELEMENT WHOSE REAL PART HAS MAX.MODULE
C     BASED ON ICAMAX BY JACK DONGARRA, LINPACK, 3/11/78.
C     MODIFIED BY NICK HIGHAM, MAY 12, 1987.
C
      COMPLEX *16 CX(*)
      DOUBLE PRECISION DMAX
      INTEGER I,INCX,IX,N
      COMPLEX *16 ZDUM
      DOUBLE PRECISION DCABS1
C
C     ... NEXT LINE IS THE ONLY MODIFICATION.
C     DCABS1(ZDUM) = DABS(DREAL(ZDUM))
      DCABS1(ZDUM) = ABS(ZDUM)

C
      ICMAX1 = 0
      IF( N .LT. 1 ) RETURN
      ICMAX1 = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX = 1
      DMAX = DCABS1(CX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DCABS1(CX(IX)).LE.DMAX) GO TO 5
         ICMAX1 = I
         DMAX = DCABS1(CX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DMAX = DCABS1(CX(1))
      DO 30 I = 2,N
         IF(DCABS1(CX(I)).LE.DMAX) GO TO 30
         ICMAX1 = I
         DMAX = DCABS1(CX(I))
   30 CONTINUE
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION DCSUM1(N,CX,INCX)
C
C     TAKES THE SUM OF THE ABSOLUTE VALUES OF A COMPLEX VECTOR AND
C     RETURNS A DOUBLE PRECISION RESULT.
C     BASED ON SCASUM BY JACK DONGARRA, LINPACK, 3/11/78.
C     MODIFIED BY NICK HIGHAM, MAY 12, 1987.
C     THE CHANGE IS TO USE THE 'GENUINE' ABSOLUTE VALUE.
C
      COMPLEX *16 CX(*)
      DOUBLE PRECISION DTEMP
      INTEGER I,INCX,N,NINCX
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO = 0.D0)
C
      DCSUM1 = ZERO
      DTEMP = ZERO
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
C       ... NEXT LINE MODIFIED.
        DTEMP = DTEMP + CDABS(CX(I))
   10 CONTINUE
      DCSUM1 = DTEMP
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
   20 DO 30 I = 1,N
C       ... NEXT LINE MODIFIED.
        DTEMP = DTEMP + CDABS(CX(I))
   30 CONTINUE
      DCSUM1 = DTEMP
      RETURN
      END


