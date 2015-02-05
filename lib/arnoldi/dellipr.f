       subroutine dellipr(iarn, wr, wi, nev, d, c, ea1, work, 
     &                   deg, degmax, che, info)
c
c Purpose
c =======
c This routine computes the Chebyshev ellipse enclosing the unwanted 
c part of the spectrum when we want to compute the eigenvalues having
c the largest or smallest real part.
c
c parameters:
c ==========
c
c  iarn : (input) integer number of eigenvalues
c
c  wr   : (input) real array of size iarn+1 contains the real part of the
c       :	                         eigenvalues
c
c  wi   : (input) real array of size iarn+1 contains the imaginary part of the
c       :	                         eigenvalues
c
c  nev  : (input) integer number of wanted eigenvalues
c
c   d   : (output) real centre of the ellipse
c
c   c   : (output) real square of the focus distance
c
c  ea1  : (output) real real semi-axis of the ellipse passing through l(nev+1)
c
c  work : (output) real array of size (iarn+1,2) working array
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
       double precision wr(*), wi(*), work(iarn+1,*), d, c, ea1
       character*2 che
c
c Local variables
c
      integer i, ivp, j, deg1
      double precision zero, two, bb, a, test, b, rr, tb, db,
     &                 eps2, dlamch, k1, kmax, ktemp, tau, one
      parameter (zero = 0.d0, one = 1.d0, two = 2.d0)
      logical inside
c
c Functions
      double precision dabs, ddel
      intrinsic dabs
      external ddel

      eps2 = one/dlamch('e')
      info = 0
      ivp = 0
c
c    Since the spectrum is symmetric about the x-axis the work 
c    is done only on half eigenvalues.
c
      do j = nev+1, iarn
     	if (wi(j).ge.zero) then
		ivp = ivp +1
		work(ivp,1) = wr(j)
		work(ivp,2) = wi(j)
	endif
      enddo
c
c    Initialize the centre d and the real semi-axis a
c
      d = (work(1,1)+work(ivp,1))/two
      a = work(1,1) -d
      k1 = dabs(wr(nev)-d)+dabs(wi(nev))
      tau = zero
      if (work(1,2).ne.zero.or.work(ivp,2).ne.zero) then
        a = a+dabs((work(1,1)-wr(nev))/two)
      endif
c
c    Check for degenerated ellipse: wr(nev+1:iarn) = cste
c                                 : a = 0, b = max(Im(lambda(nev+1:iarn)))
c
      if (a.eq.zero) then
	b = zero
	do i = 1, ivp
		if (work(i,2).gt.b) b = work(i,2)
	enddo
        if (b.ne.zero) then
		tau = k1/b
	endif
	go to 20
      endif
c
c    Compute the imaginary semi-axis b
c
      b = zero
      kmax = 1.d20
      do i = 1, ivp
	rr = work(i,1)-d
	tb = work(i,2)*a
	db = dsqrt(dabs(a*a-rr*rr))
	if (db.gt.zero) then
		bb = dabs(tb/db)
		if (bb.gt.b) then
			b = bb
			ktemp = (a+b)
			if (ktemp.lt.kmax) kmax = ktemp
		endif
	endif
      enddo
c
c    Check if the ellipse E(d,a,b) encloses the unwanted part of the 
c    spectrum.
c
      inside = .false.
      do i = 1, ivp
	test = ddel(work(i,1), work(i,2), d, a, b)
	if (test.gt.1.d-8) then
		inside = .true.
		go to 10
	endif
      enddo
 10   if (inside) then
	info = -18
        go to 30
      endif
 20   ea1 = wr(nev)-d
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
 30   return
      end
