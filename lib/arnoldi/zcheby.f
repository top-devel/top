      subroutine zcheby(n, deg, c, ea1, d, workc, flagcom,
     &                 vrcom, icycle, revcom)
c
c This routine performs the Chebyshev acceleration using the ellipse 
c parameters computed by cellip*.f.
c 
c parameter
c
c   n    : (input) integer dimension of the matrix a
c
c   deg   : (input) integer degre of the Chebychev polynomial
c
c   c    : (input) real the square of the focus distance
c
c   ea1  : (input) complex the complex semi axis
c
c    d   : (input) complex the centre of the ellipse
c
c  workc : (input) complex array of size (n,2) working array
c
c flagcom: (input) integer used for reverse communication
c
c vrcom  : (input/output) complex array of size(n,4) used for reverse
c        :                communication
c
c icycle : (input/output) integer Chebychev step
c
c revcom : (input/output) integer reverse communication parameter
c
c ARGUMENTS
c ...Scalar Arguments...
      integer n, deg, flagcom, icycle, revcom
      double precision c
      complex *16 ea1, d
c ...Array Arguments...
       complex *16 workc(n,*), vrcom(n,*)

c LOCAL DECLARATIONS
c ...Parameters...
      double precision zero, one, two
      parameter (zero = 0.d0, one = 1.d0, two = 2.d0)      
      complex *16 czero, cone, ctwo
      parameter (czero = (0.d0,0.d0) ,cone = (1.d0,0.d0),
     &           ctwo = (2.d0,0.d0))
c ...Scalar Arguments...
      integer j, flagcom1
      double precision ec
      complex *16 ec2 ,sig ,sigma, sigmam, t1, t2, t3, temp
c ...Save Statements
      save ec, sig, sigma, sigmam, flagcom1, t1, t2, t3
c ...Functions
      complex *16 dcmplx
      intrinsic dcmplx

      if (icycle.eq.0) then
c
c  Set initial values and return asking for y<-Ax
c
      	do j = 1, n
		workc(j,1) = czero
      	enddo
      	ec = dsqrt(c)/two
        ec2=dcmplx(ec,0.d0)
      	sig = ea1/ec2
      	if (abs(sig).eq.zero) sig = dcmplx(1.d-3,0.d0)
      	sigma = ctwo/sig
      	sigmam =czero
      	icycle = 1
        flagcom1 = 0
      	go to 1
      endif

      if (flagcom1.eq.4) then
        flagcom1 = 0
        temp=sqrt(vrcom(1,3))
        call zcopy(n,vrcom(1,4),1,vrcom(1,1),1)
        call zcopy(n,vrcom(1,2),1,vrcom(1,3),1)
        temp = one/temp
        call zscal(n, temp, vrcom(1,1), 1)
        call zscal(n, temp, vrcom(1,3), 1)
        call zscal(n, temp, workc(1,1), 1)
	goto 150
       endif



c
c  Three term recurence for Chebyshev polynomials
c
      t1 = sigma/ec2
      if (icycle.eq.1) t1 = t1/ctwo
      t2 = t1*d
      t3 = sigma*sigmam
       if (mod(icycle,10).eq.0) then
        call zcopy(n,vrcom(1,1),1,vrcom(1,4),1)
        call zcopy(n,vrcom(1,3),1,vrcom(1,1),1)
	call zcopy(n,vrcom(1,3),1,vrcom(1,2),1)
        flagcom1 = 4
        revcom=3
        return
       endif
150     continue
      do j = 1, n 
	workc(j,2) = t1*vrcom(j,3)-t2*vrcom(j,1)-t3*workc(j,1)
      enddo
      if (mod(icycle,10).eq.0) then
	call zscal(n, temp, workc(1,2), 1)
      endif
      call zcopy(n, vrcom(1,1), 1, workc(1,1), 1)
      call zcopy(n, workc(1,2), 1, vrcom(1,1), 1)
      sigmam = sigma
      sigma = cone/(sig-sigma)
      icycle = icycle + 1
      if (icycle.le.deg) then
	go to 1
      else
	flagcom = 9
      endif
1     revcom = 1
      return 
      end
