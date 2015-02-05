! calcul des harmoniques spheriques par fonctions elementales (fortran95)

module poly_legendre
   implicit none
   private ! by default everything is private. Don't pollute name space!
   public :: p ! list of exported routines
 contains
      elemental DOUBLE PRECISION FUNCTION P(L,M,Z)

        implicit none
        real*8,intent(in) :: z
        real*8 pq,q
        integer, intent(in) :: L,M
        integer k,i

      P=0.D0
      IF (M.GT.L) GOTO 10
        P=1d0
        DO K=1,M
        P=P*(2*K-1)
        END DO
        Q=0d0
        DO I=1,L-M
        PQ=P
        P=((2*(M+I)-1)*Z*P-(2*M+I-1)*Q)/I
        Q=PQ
        END DO
10    RETURN
      END function

end module

module mod_ylm
   use poly_legendre
   implicit none
   private        ! by default everything is private. Don't pollute name space!
   public :: ylm,dth_ylm,dphi_ylm,dtt_ylm,dtp_ylm ! list of exported routines
 contains
      elemental double precision function ylm(l,m,z) result(value)
        implicit none
        real*8,intent(in) :: z
        real*8 pi,g,en,z1
        integer, intent(in) :: l,m
        integer ma,k

      PI=3.141592653589793D0
      MA=IABS(M)
        G=1
        IF (M.EQ.0) GOTO 11
        DO K=L-MA+1,L+MA
        G=G*DSQRT(DFLOAT(K))
        END DO
11      EN=DSQRT(DFLOAT(2*L+1)/4D0/PI)/G
      Z1=DSQRT(1D0-Z*Z)
      IF (M.EQ.0) value=EN*P(L,MA,Z)
      IF (M.EQ.0) GOTO 12
      value=EN*Z1**MA*P(L,MA,Z)
12    RETURN
      end function

! dth_ylm calcule la derivee partielle d'un Y(l,m) par rapport a  theta
!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

      elemental DOUBLE PRECISION FUNCTION dth_ylm(L,M,Z)

      implicit none
      real*8, intent(in) :: z
      real*8 pi,g,z1,en
      integer, intent(in) :: m,l
      integer k,ma

      PI=3.141592653589793D0
      MA=IABS(M)
        G=1d0
        IF (MA.EQ.0) GOTO 13
        DO K=L-MA+1,L+MA
        G=G*K
        END DO
13      EN=DSQRT(DFLOAT(2*L+1)/4D0/PI/G)
      Z1=DSQRT(1D0-Z*Z)
      IF (MA.EQ.0) dth_ylm=-EN*Z1*P(L,MA+1,Z)
      IF (MA.EQ.1) dth_ylm=EN*(Z*P(L,MA,Z)-Z1**(MA+1)*P(L,MA+1,Z))
      IF (MA.LE.1) GO TO 14
      dth_ylm=EN*(MA*Z*Z1**(MA-1)*P(L,MA,Z)-Z1**(MA+1)*P(L,MA+1,Z))
14    RETURN
      END function

! la partie imaginaire de la derivee partielle de Y(l,m)
! par rapport a  phi divisee par sin(theta)
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      elemental double precision function dphi_ylm(L,M,Z)
      implicit none
      real*8, intent(in) :: z
      real*8 pi,g,z1,en
      integer, intent(in) :: m,l
      integer k,ma

      PI=3.141592653589793D0
      MA=IABS(M)
        G=1d0
        IF (MA.EQ.0) GOTO 15
        DO K=L-MA+1,L+MA
        G=G*K
        END DO
15      EN=DSQRT(DFLOAT(2*L+1)/4D0/PI/G)
      IF (MA.EQ.0) dphi_ylm=0D0
      IF (MA.EQ.1) dphi_ylm=EN*M*P(L,MA,Z)
      IF (MA.LE.1) GO TO 16
      dphi_ylm=EN*M*(DSQRT(1D0-Z*Z))**(MA-1)*P(L,MA,Z)
16    RETURN
      END function

! dtt_ylm calcule la derivee 2e partielle d'un Y(l,m) par rapport a  theta
!"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

      elemental DOUBLE PRECISION FUNCTION dtt_ylm(L,M,Z)

      implicit none
      real*8, intent(in) :: z
      real*8 pi,g,z1,en
      integer, intent(in) :: m,l
      integer k,ma

      PI=3.141592653589793D0
      MA = IABS(M)
        G=1d0
        IF (MA.EQ.0) GOTO 17
        DO K=L-MA+1,L+MA
        G=G*K
        END DO
17      EN=DSQRT(DFLOAT(2*L+1)/4D0/PI/G)
      Z1=DSQRT(1D0-Z*Z)
      IF (MA.EQ.0) dtt_ylm=EN*(-Z*P(L,MA+1,Z)+Z1**2*P(L,MA+2,Z))
      IF (MA.EQ.1) dtt_ylm=EN*(-Z1*P(L,MA,Z)-3d0*Z1*Z*P(L,MA+1,Z)&
                               +Z1**3*P(L,MA+2,Z))
      IF (MA.EQ.2) dtt_ylm=EN*((2d0-4d0*Z1**2)*P(L,MA,Z)&
                               -5d0*Z1**2*Z*P(L,MA+1,Z) &
                               +Z1**4*P(L,MA+2,Z))
      IF (MA.LE.2) GOTO 18
      dtt_ylm=EN*((DFLOAT(MA*(MA-1))*Z1**(MA-2)-MA*MA*Z1**MA)*P(L,MA,Z)&
                  -DFLOAT(2*MA+1)*Z1**MA*Z*P(L,MA+1,Z)                 &
                  +Z1**(MA+2)*P(L,MA+2,Z))
18    RETURN
      END function

! la derivee par rapport a theta de:
! la partie imaginaire de la derivee partielle de Y(l,m)
! par rapport a  phi divisee par sin(theta)
!""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
      elemental double precision function dtp_ylm(L,M,Z)
      implicit none
      real*8, intent(in) :: z
      real*8 pi,g,z1,en
      integer, intent(in) :: m,l
      integer k,ma

      PI=3.141592653589793D0
      MA=IABS(M)
        G=1d0
        IF (MA.EQ.0) GOTO 19
        DO K=L-MA+1,L+MA
        G=G*K
        END DO
        EN=DSQRT(DFLOAT(2*L+1)/4D0/PI/G)
      Z1=DSQRT(1D0-Z*Z)
19    IF (MA.EQ.0) dtp_ylm=0D0
      IF (MA.EQ.1) dtp_ylm=-EN*DFLOAT(M)*(Z1**MA*P(L,MA+1,Z))
      IF (MA.EQ.2) dtp_ylm=EN*DFLOAT(M)*(Z*P(L,MA,Z)    &
                                   -Z1**MA*P(L,MA+1,Z))
      IF (MA.LE.2) GO TO 20
      dtp_ylm=EN*DFLOAT(M)*(DFLOAT(MA-1)*Z1**(MA-2)*Z*P(L,MA,Z)&
                           -Z1**MA*P(L,MA+1,Z))
20    RETURN
      END function

end module
