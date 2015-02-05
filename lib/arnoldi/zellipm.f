       subroutine zellipm(iarn, w, nev, d, c, ea1, n, iord, 
     &                    deg, degmax, che, info)
c
c This routine computes the Chebyshev ellipse enclosing the unwanted 
c part of the spectrum when we want to compute the eigenvalues having
c the largest modulii.
c
c parameters:
c
c  iarn : (input) integer number of eigenvalues
c
c  w    : (input) complex array of size iarn+1 contains the eigenvalues
c
c  nev  : (input) integer number of wanted eigenvalues
c
c   d   : (output) complex centre of the ellipse
c
c   c   : (output) real square of the focus distance
c
c  ea1  : (output) complex parameter used for the Chebychev 
c         acceleration
c
c iord  : (output) integer array of size n working array
c
c  deg  : (input/output) integer degree of the Chebyshev polynomial
c
c degmax: (input) integer maximum degree of the Chebyshev polynomial
c
c che   : (input) character specifies the way deg is set
c
c info  : (output) integer error flag: info = 0 normal exit
c       :                            : info = -18 no possible ellipse
c
       integer nev, iarn, info, n, iord(*), deg,
     &        degmax
       double precision c
       complex *16 w(*), d, ea1
       character*2 che
c
c Local variables
c
      integer i, ivp, j, icur, test, jmax, test1, test2, deg1
      double precision zero, two, bb, a, b, tb, db, inside,
     &                 aa, outside, eps2, dlamch, k1,
     &                 kmax, ktemp, tau, one, minIm, maxIm,
     &                 minRe, maxRe, xd, yd
      complex *16 czero, cone, ctwo 

      parameter (zero = 0.d0, one = 1.d0, two = 2.d0)
      parameter (czero = (0.d0, 0.d0), cone = (1.d0, 0.d0),
     &           ctwo = (2.d0, 0.d0))
c ...Functions...
      double precision dreal, dimag, dabs, dsqrt, zdel
      complex *16 dcmplx
      external zdel
      intrinsic dreal, dimag, dabs, dsqrt, dcmplx

c
      eps2 = one/dlamch('e')
      info = 0
      do i = 1, nev
	iord(i) = 0
      enddo
      do i = nev+1, iarn
	iord(i) = 1
      enddo
      icur = 0
      test = 0
      do while (test.eq.0.and.icur.le.3)
c
c Compute the centre d and the semi-axis a
c
	ivp = iarn-nev 
        minRe=10d16
        maxRe=-10d16
        minIm=10d16
        maxIm=-10d16
        do i=nev+1, iarn
           if (dimag(w(i)).lt.minIm) minIm=dimag(w(i))
           if (dimag(w(i)).gt.maxIm) maxIm=dimag(w(i))
           if (dreal(w(i)).lt.minRe) minRe=dreal(w(i))
           if (dreal(w(i)).gt.maxRe) maxRe=dreal(w(i))
        enddo
        xd=(minRe+maxRe)/2
        yd=(minIm+maxIm)/2
	d = dcmplx(xd, yd)
	k1 = dabs(dreal(w(nev)-xd))+
     &       dabs(dimag(w(nev)-yd))
	tau = zero
	a = zero
	jmax = 0
	do j = nev+1, iarn
		if (iord(j).eq.1) then
			aa = dabs(dreal(w(j))-xd)
			if (aa.gt.a) then
				a = aa
				jmax = j
			endif
		endif
	enddo
	if ((dimag(w(jmax))-yd).ne.zero) then
		a = a+((dsqrt((dreal(w(nev+1)-w(nev)))**2+
     &                 (dimag(w(nev+1)-w(nev)))**2)))/two
	endif
c
c Compute the semi-axis b
c
	b = zero
	jmax = 0
        kmax = 1.d2
	do j = nev+1, iarn
		if (iord(j).eq.1) then
			tb = dabs(a*(dimag(w(j))-yd))
			db = dsqrt(dabs(a*a)-(dreal(w(j))-xd)**2)
			if (db.ne.zero) then
				bb = tb/db
				if (bb.gt.b) then
					b = bb
					ktemp = (a+b)
					if (ktemp.lt.kmax) 
     &                                    	kmax = ktemp
					jmax = j
				endif
			endif
		endif
	enddo
c
c Check if the Ellipse E(d, a, b) encloses the unwanted part 
c of the spectrum and if wanted part of the spectrum is not 
c enclosed by the ellipse
c
	test1 = 0
	do i = 1, nev
		inside = zdel(w(i), d, a, b)
		if (inside.lt.1.d-10) then
			test1 = 1
			goto 1
		endif
	enddo
	test2 = 0
	do i = nev+1, iarn
		if (iord(i).eq.1) then
			outside = zdel(w(i), d, a, b)
			if (outside.gt.1.d-10) then
				test2 = 1
				goto 1
			endif
		endif
	enddo
1  	if (test1.eq.1.or.test2.eq.1) then
		icur = icur+1
		iord(jmax) = -2
	else
		test = 1
	endif	
      enddo
      if (test.ne.1) then
	info = -18
        go to 30
      endif
      ea1 = w(nev)-d
      c = dabs(a*a-b*b)
      if (che.eq.'DY') then
        deg1 = 0
        tau = k1/kmax
        if (tau.ne.zero) then
                deg1 = int((one/two)*(one+dabs(
     &              (dlog(eps2)/dlog(tau)))))
        endif
        if (deg1.le.0.or.deg1.gt.degmax) then
                deg = degmax
        else
                deg = deg1
        endif
        write(9,26)deg1
26      format('Computed degree of Chebyshev polynomial: ',i4)
      endif
      write(9,27)deg
27    format('Degree of Chebyshev polynomial: ',i4)
30    return 
      end
