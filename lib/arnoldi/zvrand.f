      subroutine zvrand(n, v)
c
c This routine builds a random complex vector. (Qmrpack routine
c written by N. Nachtigal)
c 
c parameters
c
c  n  : (input) integer dimension of the vector
c
c  v  : (output) output random vector
c
      integer n
      complex *16 v(n)
c
c    Local variables
c
      integer im, i, j, seed,imax 
      double precision dmax,realv,vr,vi
c ...Functions...
      double precision dble
      complex *16 dcmplx 
      intrinsic dble, dcmplx

      im = 1
      j = 0
      do i = 1, 31
	j = j+1
	if (2*im.le.im) go to 1
	im = 2*im 
      enddo
1     imax = (im-1)*2+1
      dmax = dble(imax)
      do i = 1, mod(j,3)
	j = j-1
	im = im/2
      enddo
      im = im+5
      seed = iabs(mod(im*30107,imax))
      seed = (seed/2)*2+1
      do i = 1,n
	realv=dble(seed)/dmax
	seed=iabs(mod(im*seed,imax))
        vr=realv
        realv=dble(seed)/dmax
        seed=iabs(mod(im*seed,imax))
	vi=realv
	v(i)=dcmplx(vr,vi)
      enddo
      return
      end
