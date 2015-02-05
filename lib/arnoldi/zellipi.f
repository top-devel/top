      subroutine zellipi(iarn, w, nev, d, c, ea1, work,
     &                   deg, degmax, che, info)
c
c This routine computes the Chebyshev ellipse enclosing the unwanted
c part of the spectrum when we want to compute the eigenvalues having
c the largest imaginary part.
c
c parameters:
c
c  iarn : (input) integer number of eigenvalues
c
c  w    : (input) complex array of size iarn+1 the eigenvalues
c
c  nev  : (input) integer number of wanted eigenvalues
c
c   d   : (output) complex centre of the ellipse
c
c   c   : (output) real square of the focus distance
c
c  ea1  : (output) real real semi-axis of the ellipse passing through l(nev+1)
c
c  work : (output) complex array of size (iarn+1) working array
c
c  deg  : (input/output) integer degree of the Chebyshev polynomial
c
c degmax: (input) integer maximum degree of the Chebyshev polynomial
c
c  info : (output) integer error flag: info = 0 normal exit
c       :                            : info = -18 no possible ellipse
c
       integer nev, iarn, info, deg, degmax
       double precision  c
       complex *16 w(*), work(*), d,ea1
       character*2 che
c
c Local variables
c
      integer ivp, i, deg1, etape
      double precision zero, two, b, a, aa, test, rr, db, tb,
     &                 eps2, dlamch, k1, kmax, ktemp, tau, one,
     &                 xd,yd,minRe,maxRe,minIm,maxIm
      parameter (zero = 0.d0, one = 1.d0, two = 2.d0)
      logical inside
c ...Functions...
      double precision dreal, dimag, dabs, dsqrt, zdel
      complex *16 dcmplx
      external zdel
      intrinsic dreal, dimag, dabs, dsqrt, dcmplx      
c
      eps2 = one/dlamch('e')  
      info = 0
      ivp = 0

      do i = nev+1, iarn
	ivp = ivp +1
	work(ivp) = w(i)
      enddo		
c
c    Initialize the centre d and the imaginary semi-axis b
c
      minRe=10d16
      maxRe=-10d16
      minIm=10d16
      maxIm=-10d16
      do i=1,ivp 
        if (dimag(work(i)).lt.minIm) minIm=dimag(work(i))
        if (dimag(work(i)).gt.maxIm) maxIm=dimag(work(i))
        if (dreal(work(i)).lt.minRe) minRe=dreal(work(i))
        if (dreal(work(i)).gt.maxRe) maxRe=dreal(work(i))
      enddo
      xd=(minRe+maxRe)/2
      yd=(minIm+maxIm)/2
      d=dcmplx(xd,yd)
      b=dimag(work(1))-yd
      k1 = dabs(dreal(w(nev))-xd)+dabs(dimag(w(nev)))
      tau = zero
c
c    Check for degenerated ellipse: dimag(w(nev+1:iarn)) = 0
c                                 : a=max(Re(lambda(nev+1:iarn)-xd), b = 0
c
      if (dimag(work(1)).eq.zero) then
        a = zero
	do i = 1, ivp
		aa = dabs(dreal(work(i))-xd)
		if (aa.gt.a) a = aa
	enddo
	if (a.ne.zero) then
		tau = k1/a
	endif
        go to 20
      endif
c
c    Compute the real semi-axis a
c
      etape =0
5     etape=etape+1 
      a = zero
      kmax = 1.d20
      do i = 1,ivp
	rr = dreal(work(i))-xd
	tb = dimag(work(i))-yd
	db = dsqrt(dabs(b*b-tb*tb))
	tb = dabs(rr*b)
	if (db.gt.zero) then
		aa = dabs(tb/db)
		if (aa.gt.a) then
			a = aa
			ktemp = (a+b)
			if (ktemp.lt.kmax) kmax = ktemp
		endif
	endif
      enddo

c
c    Check for degenerated ellipse if a=0
c    b is increased
c
      if (a.eq.zero) then
	b=b+dabs(dimag(w(nev))-dimag(work(1)))/two
	go to 5 
      endif
c
c    Check if the ellipse E(d,a,b) encloses the unwanted part of the
c    spectrum.
c
      inside = .false.

      test = zdel(work(1) , d, a, b)
      if (test.gt.1.d-8) then
		inside=.true.
		if (etape.eq.1) then
			b=b+dabs(dimag(w(nev))-dimag(work(1)))/two
			go to 5
		endif		
                go to 10
      endif 
      test = zdel(work(ivp) , d, a, b)
      if (test.gt.1.d-8) then
                inside=.true.
                if (etape.eq.1) then
                        b=b+dabs(dimag(w(nev))-dimag(work(1)))/two
                        go to 5
                endif
                go to 10
      endif
 
      do i = 2, ivp-1
         test = zdel(work(i) , d, a, b)
         if (test.gt.1.d-8) then
		inside = .true.
		go to 10
	endif
      enddo
10    if (inside) then
       	info = -18
      	go to 30
      endif
20    ea1 = w(nev)-d
      c = dabs(a*a-b*b)
      if (che.eq.'DY') then
        deg1 = 0
        if (dimag(work(1)).ne.zero) tau = k1/kmax
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
