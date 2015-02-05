       subroutine dordon(n, iord, wr, wi, iarn, y, eig)
c
c Purpose
c =======
c This routine sorts the computed Ritz values following the type 
c of eigen-computation defined by the user and coded by the 
c character eig.
c
c parameter:
c =========
c
c  n    : (input) integer leading dimension of the working arrays
c
c  iord : (output) integer array of size iarn+1 contains the new reordering
c
c  wr   : (input/output) real array of size iarn+1 contains the real parts of 
c       :		 the Ritz values
c
c  wi   : (input/output) real array of size iarn+1 contains the imaginary parts
c       :		 of the Ritz values
c
c iarn  : (input) integer the size of the projected problem
c
c  y    : (output) real array of size 2*(iarn+1) working array
c
c eig   : (input) character determines the type of eigen-computation
c
       integer n, iarn, iord(*)
       double precision wr(*), wi(*), y(*)
       character*2 eig
c
c Local variables
c
       double precision tb, xm, x
       integer j, im, k, i
c Functions
       double precision dabs, dsqrt
       intrinsic dabs, dsqrt
c
       tb = -1.d20
       if (eig.eq.'LR'.or.eig.eq.'SR') then 
       	call dcopy(iarn, wr, 1, y, 1)
       elseif (eig.eq.'LI') then
	do i = 1, iarn
		y(i) = dabs(wi(i))
	enddo
       elseif (eig.eq.'LM'.or.eig.eq.'SH') then
	do i = 1, iarn
		y(i) = dsqrt(wr(i)**2+wi(i)**2)
	enddo
       endif
c
c  Sort the Ritz values following the value of eig
c
       do i = 1, iarn
	xm = tb
	do j = 1, iarn
		x = y(j)
		if (x.gt.xm) then
			xm = x
			im = j
		endif
	enddo
	iord(i) = im
	y(im) = tb
       enddo
c
c  iord now contains the new ordering: resort wr and wi
c
       j = 1
       do i = 1, iarn
	k = iord(i)
	y(j) = wr(k)
	y(j+1) = wi(k)
	j = j+2       
       enddo
       j = 1
       do i = 1, iarn
	wr(i) = y(j)
	wi(i) = y(j+1)
	j = j+2
       enddo
       return
       end
