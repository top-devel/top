      subroutine dcheby(n, deg, c, ea1, d, workc, flagcom,
     &                 vrcom, icycle, revcom)
c
c This routine performs the Chebyshev acceleration using the ellipse 
c parameters computed by dellip*.f.
c 
c parameter
c
c   n    : (input) integer dimension of the matrix a
c
c  deg   : (input) integer degre of the Chebychev polynomial
c
c   c    : (input) real the square of the focus distance
c
c   ea1  : (input) real the real semi axis
c
c    d   : (input) real the centre of the inverval
c
c  workc : (input) real array of size (n,2) working array
c
c flagcom: (input) integer used for reverse communication
c
c vrcom  : (input/output) real array of size(n,4) used for reverse
c        :                communication
c
c icycle : (input/output) integer Chebychev step
c
c revcom : (input/output) integer reverse communication parameter
c
c ARGUMENTS
c ...Scalar Arguments...
      integer n, deg, flagcom, icycle, revcom
      double precision c, ea1, d
c ...Array Arguments...
      double precision workc(n,*), vrcom(n,*)

c LOCAL DECLARATIONS
c ...Parameters...
      double precision zero, one, two
      parameter (zero = 0.d0, one = 1.d0, two = 2.d0)
c ...Local Scalars
      integer j, flagcom1
      double precision ec, sig, sigma, sigmam, t1, t2, t3, temp
c ...Save statements
      save ec, sig, sigma, sigmam, t1, t2, t3, temp, flagcom1, j
c ...Intrinsic functions
      double precision dsqrt, dabs
      intrinsic dsqrt, dabs
      
      if (icycle.eq.0) then
c
c  Set initial values and return asking for y<-Ax
c
        do j = 1, n
         workc(j,1) = zero
        enddo
        ec = dsqrt(c)/two
        sig = ea1/ec
        if (dabs(sig).eq.zero) sig = 1.d-3
        sigma = two/sig
        sigmam = zero
        icycle = 1
        flagcom1 = 0
        go to 1
      endif
c
      if (flagcom1.eq.4) then
        flagcom1 = 0
        temp = dsqrt(vrcom(1,3))
        call dcopy(n,vrcom(1,4),1,vrcom(1,1),1)
        call dcopy(n,vrcom(1,2),1,vrcom(1,3),1)
        temp = one/temp
        call dscal(n, temp, vrcom(1,1), 1)
        call dscal(n, temp, vrcom(1,3), 1)
        call dscal(n, temp, workc(1,1), 1)
        goto 150
      endif
c
c  Three term recurence for Chebyshev polynomials
c
      t1 = sigma/ec
      if (icycle.eq.1) t1 = t1/two
      t2 = t1*d
      t3 = sigma*sigmam
       if (mod(icycle,10).eq.0) then
         call dcopy(n,vrcom(1,1),1,vrcom(1,4),1)
         call dcopy(n,vrcom(1,3),1,vrcom(1,1),1)
         call dcopy(n,vrcom(1,3),1,vrcom(1,2),1)
         flagcom1 = 4 
         revcom = 3
         return
       endif
150   do j = 1, n 
        workc(j,2) = t1*vrcom(j,3)-t2*vrcom(j,1)-t3*workc(j,1)
      enddo
      if (mod(icycle,10).eq.0) then
        call dscal(n, temp, workc(1,2), 1)
      endif
      call dcopy(n, vrcom(1,1), 1, workc(1,1), 1)
      call dcopy(n, workc(1,2), 1, vrcom(1,1), 1)
      sigmam = sigma
      sigma = one/(sig-sigma)
      icycle = icycle + 1
      if (icycle.le.deg) then
        goto 1
      else
        flagcom = 9 
      endif
1     revcom = 1
      return 
      end
