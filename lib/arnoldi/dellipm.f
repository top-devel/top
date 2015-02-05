       subroutine dellipm(iarn, wr, wi, nev1, d, c, ea1, n, iord, 
     &                    deg, degmax, che, info)
c
c Purpose
c =======
c This routine computes the Chebyshev ellipse enclosing the unwanted 
c part of the spectrum when we want to compute the eigenvalues having
c the largest modulii.
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
c  nev1 : (input) integer number of wanted eigenvalues
c
c   d   : (output) real centre of the ellipse
c
c   c   : (output) real square of the focus distance
c
c  ea1  : (output) real real semi-axis of the ellipse passing through l(nev+1)
c
c iord  : (output) integer array of size n working array
c
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
       integer nev1, iarn, info, n, iord(*), deg,
     &        degmax
       double precision wr(*), wi(*), d, c, ea1
       character*2 che
c
c Local variables
c
      integer i, ivp, j, icur, test, jmax, test1, test2, deg1
      double precision zero, two, bb, a, b, tb, db, inside,
     &                 dmin, dmax, aa, outside, eps2, k1,
     &                 kmax, ktemp, tau, one, dlamch
      parameter (zero = 0.d0, one = 1.d0, two = 2.d0)
c
c Functions
      double precision ddel, dabs
      intrinsic dabs
      external ddel

      eps2 = one/dlamch('e')
      info = 0
c
c Selection of the wanted eigenvalues
c iord(i)=1 for a wanted eigenvalue
c Since the spectrum is symmetric about the x-axis
c the work is done only on half eigenvalues
c
      do i = 1, nev1
	iord(i) = 0
      enddo
      do i = nev1+1, iarn
	iord(i) = 1
      enddo
      icur = 0
      test = 0
      do while (test.eq.0.and.icur.le.3)
	ivp = 0
	do j = nev1+1, iarn
		if (iord(j).eq.1) then
			if (wi(j).ge.zero) then
				ivp = ivp+1
			else
				iord(j) = -1
			endif
		endif
	enddo
c
c Compute the centre d and the real semi-axis a
c
	dmin = 1.d20
	dmax = -dmin
	do j = nev1+1, iarn
		if (iord(j).eq.1) then
			if (wr(j).lt.dmin) dmin = wr(j)
			if (wr(j).gt.dmax) dmax = wr(j)
		endif
	enddo
	d = (dmin+dmax)/two
	k1 = dabs(wr(nev1)-d)+dabs(wi(nev1))
	tau = zero
	a = zero
	jmax = 0
	do j = nev1+1, iarn
		if (iord(j).eq.1) then
			aa = dabs(wr(j)-d)
			if (aa.gt.a) then
				a = aa
				jmax = j
			endif
		endif
	enddo
	if (wi(jmax).ne.zero) then
		a = a+(dsqrt((wr(nev1+1)-wr(nev1))**2+
     &                 (wi(nev1+1)-wi(nev1))**2))
	endif
c
c Compute the imaginary semi-axis b
c
	b = zero
	jmax = 0
        kmax = 1.d2
	do j = nev1+1, iarn
		if (iord(j).eq.1) then
			tb = dabs(a*wi(j))
			db = dsqrt(dabs(a*a)-(wr(j)-d)**2)			
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
c Check if the ellipse E(d,a,b) encloses the unwanted part of
c the spectrum and if the wanted part of the spectrum is not
c enclosed by the ellipse 
c 
	test1 = 0
	do i = 1, nev1
		inside = ddel(wr(i), wi(i), d, a, b)
		if (inside.lt.1.d-10) then
			test1 = 1
			goto 1
		endif
	enddo
	test2 = 0
	do i = nev1+1, iarn
		if (iord(i).eq.1) then
			outside = ddel(wr(i), wi(i), d, a, b)
			if (outside.gt.1.d-10) then
				test2 = 1
				goto 1
			endif
		endif
	enddo
1  	if (test1.eq.1.or.test2.eq.1) then
		icur = icur+1
		iord(jmax) = -2
		if (wi(jmax).ne.zero) iord(jmax-1) = -2
	else
		test = 1
	endif	
      enddo
      if (test.ne.1) then
	info = -18
        go to 30
      endif
      ea1 = wr(nev1)-d
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
 30   return
      end
