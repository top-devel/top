      subroutine dellipi(iarn, wr, wi, nev, d, c, ea1, work,
     &                   deg, degmax, che, info)
c
c Purpose
c =======
c This routine computes the Chebyshev ellipse enclosing the unwanted
c part of the spectrum when we want to compute the eigenvalues having
c the largest imaginary part.
c
c parameters:
c ==========
c
c  iarn : (input) integer number of eigenvalues
c
c  wr   : (input) real array of size iarn+1 contains the real part of the
c       :                                eigenvalues
c
c  wi   : (input) real array of size iarn+1 contains the imaginary part of the
c       :                                eigenvalues
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
c  info : (output) integer error flag: info = 0 normal exit
c       :                            : info = -18 no possible ellipse
c
       integer nev, iarn, info, deg, degmax
       double precision wi(*), wr(*), work(iarn+1,*), d, ea1, c
       character*2 che
c
c Local variables
c
      integer ivp, i, deg1
      double precision zero, two, b, a, aa, test, rr, db, tb,
     &                 eps2, k1, kmax, ktemp, tau, one, dlamch
      parameter (zero = 0.d0, one = 1.d0, two = 2.d0)
      logical inside
c
c Functions
      double precision dabs
      intrinsic dabs
      double precision ddel
      external ddel

      eps2 = one/dlamch('e')
      info = 0
      ivp = 0
c
c    Since the spectrum is symmetric about the x-axis the work
c    is done only on half eigenvalues.
c
      do i = nev+1, iarn
	if (wi(i).ge.zero) then
		ivp = ivp +1
		work(ivp,1) = wr(i)
		work(ivp,2) = wi(i)
	endif
      enddo		
c
c    Initialize the centre d and the imaginary semi-axis b
c
      d = work(1,1)
      b = work(1,2)
      k1 = dabs(wr(nev)-d)+dabs(wi(nev))
      tau = zero
      if (work(1,2).eq.work(2,2)) then
	b = b+dabs((dabs(work(1,2))-dabs(wi(nev)))/two)
      endif	
c
c    Check for degenerated ellipse: wi(nev+1:iarn) = 0
c                                 : a=max(Re(lambda(nev+1:iarn)-d), b = 0
c
      if (work(1,2).eq.zero) then
        a = zero
	do i = 1, ivp
		aa = dabs(work(i,1)-d)
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
      a = zero
      kmax = 1.d20
      do i = 2, ivp
	rr = work(i,1)-d
	tb = work(i,2)
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
10    if (inside) then
       	info = -18
      	go to 30
      endif
20    ea1 = wr(nev)-d
      c = dabs(a*a-b*b)
      if (che.eq.'DY') then
        deg1 = 0
        if (work(1,2).ne.zero) tau = k1/kmax
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
