       subroutine zordon(n, iord, w, iarn, yr, yc, eig)
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
c  w    : (input/output) complex array of size iarn+1 contains the Ritz values 
c
c iarn  : (input) integer the size of the projected problem
c
c  yc   : (output) complex array of size iarn+1 working array
c
c  yr   : (output) real array of size iarn+1 working array
c
c eig   : (input) character determines the type of eigen-computation
c
       integer n, iarn, iord(*)
       double precision yr(*)
       complex *16 w(*), yc(*) 
       character*2 eig
c
c Local variables
c
       double precision tb, xm, x
       integer j, im, k, i
c ...Functions...
       double precision dreal, dimag
       intrinsic dreal, dimag
c
       tb = -1.d20
       if (eig.eq.'LR'.or.eig.eq.'SR') then 
	do i=1,iarn
		yr(i) = dreal(w(i))
	enddo
       elseif (eig.eq.'LI') then
	do i = 1, iarn
		yr(i) =dimag(w(i))
	enddo
       elseif (eig.eq.'LM'.or.eig.eq.'SH') then
	do i = 1, iarn
		yr(i) = abs(w(i))
	enddo
       endif
c
c  Sort the Ritz values following the value of eig
c
       do i = 1, iarn
	xm = tb
	do j = 1, iarn
		x = yr(j)
		if (x.gt.xm) then
			xm = x
			im = j
		endif
	enddo
	iord(i) = im
	yr(im) = tb
       enddo
c
c  iord now contains the new ordering: resort wr and wi
c
       j = 1
       do i = 1, iarn
	k = iord(i)
	yc(j) = w(k)
	j = j+1       
       enddo
       do i = 1, iarn
	w(i) = yc(i)
       enddo
       return
       end
