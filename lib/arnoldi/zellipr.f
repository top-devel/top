       subroutine zellipr(iarn, w, nev, d, c, ea1, work, 
     &                   deg, degmax, che, info)
c
c This routine computes the Chebyshev ellipse enclosing the unwanted 
c part of the spectrum when we want to compute the eigenvalues having
c the largest or smallest real part.
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
c  ea1  : (output) complex semi-axis of the ellipse passing through l(nev+1)
c
c  work : (output) complex array of size (iarn+1) working array
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
       integer nev, iarn, info, deg, degmax
       double precision  c
       complex *16 w(*),work(*),d, ea1
       character*2 che
c
c Local variables
c
      integer i, ivp, deg1,etape, j
      double precision bb, a, test, b, rr, tb, db,
     &                 eps2, k1, dlamch, kmax, ktemp, tau,
     &                 xd,yd,minIm,maxIm,
     &                 zero,one,two
      parameter ( zero= 0.d0, one= 1.d0, two= 2.d0 )
      logical inside
c Functions
      double precision dreal, dimag, dabs, dsqrt, zdel
      complex *16 dcmplx
      external zdel
      intrinsic dreal, dimag, dabs, dsqrt, dcmplx

c
      eps2=one/dlamch('e')
      info = 0
      ivp = 0

      do j=nev+1,iarn
        ivp=ivp+1
        work(ivp)=w(j)
      enddo

c
c    Initialize the centre d and the real semi-axis a
c
      xd = (dreal(work(1))+dreal(work(ivp)))/two
      a = dreal(work(1)) -xd
c
c  calcul de yd
c

      minIm=10d16
      maxIm=-10d16
      do i=1,ivp
	if (dimag(work(i)).lt.minIm) minIm=dimag(work(i))
        if (dimag(work(i)).gt.maxIm) maxIm=dimag(work(i))
      enddo
      yd=(minIm+maxIm)/2 
      d=dcmplx(xd,yd)
      k1=dabs(dreal(w(nev))-xd)+dabs(dimag(w(nev)))
      tau = zero
c
c    Check for degenerated ellipse: wr(nev+1:iarn) = cste
c                                 : a = 0, b = max(Im(lambda(nev+1:iarn)))
c
      if (a.eq.zero) then
	b=maxIm-yd
        if (b.ne.zero) then
		tau = k1/b
	endif
	go to 20
      endif
c
c    Compute the imaginary semi-axis b
c
      etape=0
5     etape=etape+1
      b = zero
      kmax = 1.d20
      do i = 1,ivp
	rr = dreal(work(i))-xd
	tb = a*(dimag(work(i))-yd)
	db = dsqrt(dabs(a*a-rr*rr))
	if (db.gt.0.d0) then
		bb = dabs(tb/db)
		if (bb.gt.b) then
			b = bb
			ktemp = (a+b)
			if (ktemp.lt.kmax) kmax = ktemp
		endif
	endif
      enddo
c
c    Check for degenerated ellipse b=0 ,a=a+(real(w(nev))-real(work(1)))/2
c
c    Check if the ellipse E(d,a,b) encloses the unwanted part of the 
c    spectrum.
c

c      if (b.eq.zero.and.etape.eq.1) then 
c	a=a+dabs(dreal(work(1))-dreal(w(nev)))/2.d0
c	go to 5
c      endif

      inside = .false.
      test = zdel(work(1), d, a, b)
      if (test.gt.1.d-8) then
                inside = .true.
		if (etape.eq.1) then 
			a=a+dabs(dreal(work(1))-dreal(w(nev)))/2.d0
			go to 5
		endif
                go to 10
      endif

      test = zdel(work(ivp), d, a, b)
      if (test.gt.1d-8) then
                inside = .true.
                if (etape.eq.1) then
                        a=a+dabs(dreal(work(1))-dreal(w(nev)))/2.d0
                        go to 5
                endif
                go to 10
      endif

      do i = 2, ivp-1
	test = zdel(work(i), d, a, b)
	if (test.gt.1.d-8) then
		inside = .true.
		go to 10
	endif
      enddo
 10   if (inside) then
	info = -18
        go to 30
      endif
 20   ea1 = w(nev)-d
      c = dabs(a*a-b*b)
      if (che.eq.'DY') then
	deg1 = 0
	if (a.ne.zero) tau = k1/kmax
	if (tau.ne.zero) then
		deg1 = int((one/two)*(one+dabs(
     &            (dlog(eps2)/dlog(tau)))))
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

