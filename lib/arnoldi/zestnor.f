       subroutine zestnor(n, nor, normA, iord, flagcom, vrcom, 
     &                   revcom)
c
c This routine computes an estimate of the norm of the matrix.
c
c parameters:
c 
c  n     : (input) integer size of the matrix
c
c nor    : (input) character if nor = 'GI' ||A|| is given
c        :                   if nor = 'PM' ||A|| is computed by n matrix-vector
c        :                                 products
c        :                   if nor = 'HH' ||A|| is computed using Higham's 
c        :                                 modification of Hager'algorithm
c  
c normA  : (input/output) real norm estimate
c
c iord   : (output) integer array of size n working array
c
c flagcom: (input/output) integer reverse comunication flag
c
c vrcom  : (input/output) complex array of size (n,4) use for reverse 
c        :                                           comunication
c
c revcom : (input/output) integer reverse comunication parameter
c
      integer flagcom, n, iord(*), revcom
      double precision normA
      complex *16 vrcom(n,4)
      character*2 nor

c
c Local variables
c
      integer icol, j
      save icol
      double precision zero, one
      complex *16 czero, cone
      parameter (zero = 0.d0, one = 1.d0)
      parameter (czero = (0.d0,0.d0) ,cone = (1.d0,0.d0))
c
c ...Functions...
      double precision dsqrt
      intrinsic dsqrt

      if (nor.eq.'GI') then 
c
c   ||A|| is given: update flagcom for Arnoldi step
c
	flagcom = 3
      elseif (nor.eq.'PM') then
c
c   ||A|| is computed by v<-A ei = A(:,i)  and ||A|| = sqrt(sum(norm(A(:,i),2)^2))
c 
	if (flagcom.eq.0) then
		icol = 1
		normA = zero
 10	continue
		do j = 1, n
			vrcom(j,1) =czero
		enddo
		vrcom(icol,1) = cone
		flagcom = 1
                revcom = 1
c
c    Ask the user to perform 
c           vrcom(:,3) = A vrcom(:,1)
c
		return
	elseif (flagcom.eq.1) then
		do j = 1, n
			normA = normA+dreal(vrcom(j,3))**2
     &                               +dimag(vrcom(j,3))**2     
		enddo
		icol = icol+1
		if (icol.le.n) then
		do j = 1, n
                        vrcom(j,1) =czero
                enddo
                vrcom(icol,1) = cone
                flagcom = 1
                revcom = 1
c
c    Ask the user to perform
c           vrcom(:,3) = A vrcom(:,1)
c
                return
		endif
	endif
	normA = dsqrt(normA)
	flagcom = 3
      elseif (nor.eq.'HH') then
c
c   ||A|| is computed using Higham's modification of Hager's algorithm
c   the routines y<-Ax and y<-A^* x are required
c
        call zhager (n, vrcom, iord, normA, flagcom, revcom)
      endif
      return
      end
