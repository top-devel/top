      subroutine zhessen(n, iarn, v, h, u, hmm,
     &                   flagcom, vrcom, revcom, ierr)
c
c  Purpose
c  =======
c  This routine computes the orthonormal basis V of the 
c  Krylov subspace K(u,A)= (u, Au, A^2u, ..., A^m-1 u) 
c  and the corresponding upper Hessenberg matrix H such 
c  that:
c                                     T
c       A V = V H + h(m+1,m) v(m+1) em
c
c Parameters:
c ==========
c
c  n     : (input) integer size of the matrix.
c
c iarn   : (input) integer size of the projection.
c
c  v     : (output) complex array of size(n,iarn+1) contains the basis of the 
c        :          Krylov subspace.
c
c  h     : (output) complex array of size(iarn+1,iarn+1) contains the upper Hessenberg
c        :          matrix.
c
c  u     : (input) complex array of size n contains the starting vector.
c
c hmm    : (output) complex contains h(iarn+1,iarn) used for the stopping 
c        :          criterion
c
c flagcom: (input/output) integer reverse communication flag
c
c vrcom  : (input/output) complex array of size(n,4) vrcom(:,1) contains the 
c        :                input vector and vrcom(:,3) contains the output 
c        :                vector when the user is asked to perform 
c        :                vrcom(:,3)<--A vrcom(:,1)
c
c revcom : (input/output) integer reverse communication parameter
c
c ierr   : (output) integer error flag
c
c ARGUMENTS
c ...Scalar Arguments...
      integer n, iarn, flagcom, revcom, ierr
      complex *16 hmm
c ...Array Arguments...
      complex *16 v(n,*), h(iarn+1,iarn+1), u(*), vrcom(n,4)

c LOCAL DECLARATIONS
c ...Parameters...
      complex *16 czero, cone, cfour
      parameter (czero = (0.d0,0.d0), cone = (1.d0,0.d0),
     &           cfour= (4.d0,0.d0))
c ...Local Scalars...
      integer i, j, istep, it
      complex *16 test, t, hinorm
c ...Save Statements...
      save istep, j, it
      save test, hinorm, t

c ...Functions
c     complex *16 zsqrt
c     intrinsic zsqrt 
c
c
c Drives the reverse communication
c
      if (flagcom.eq.3) then
          goto 11
      else if (flagcom.eq.4) then
          goto 12
      else if (flagcom.eq.5) then
          goto 13
      else if (flagcom.eq.6) then
          goto 14
      else if (flagcom.eq.7) then
          goto 15
      endif

11    continue
	do j = 1, iarn+1
		do i = 1, iarn+1
			h(i,j) = czero
		enddo
	enddo
c
c     Normalize the starting vector u
c
c copy u in vrcom(:,1) and vrcom(:,2) and ask the user 
c to perform vrcom(:,3) <-- vrcom(:,1) vrcom(:,2)
c

        call zcopy(n,u,1,vrcom(1,1),1)
        call zcopy(n,u,1,vrcom(1,2),1)
        flagcom=5
        revcom=3
        return
c
13      continue
	t = sqrt(vrcom(1,3))
	t = cone/t
	call zscal(n, t, u, 1)
	call zcopy(n, u, 1, v(1,1), 1)
	test = cfour
	istep = 1
c
c     Copy v(:,istep) in vrcom(:,1) and ask the user to perform 
c              vrcom(:,2)<-- A vrcom(:,1)
c
 10 	call zcopy(n, v(1,istep), 1, vrcom(1,1), 1)
	flagcom = 4
        revcom = 1
	return
c
12      continue
	call zcopy(n, vrcom(1,3), 1, v(1,istep+1), 1)
	it = 0
 20	hinorm = czero 
	it = it+1
c
c     Compute H(1:istep,istep)
c
        j = 1
18      continue
        call zcopy(n,v(1,j),1,vrcom(1,1),1)
        call zcopy(n,v(1,istep+1),1,vrcom(1,2),1)
        flagcom = 6
        revcom = 3
        return
14      continue
        t = vrcom(1,3)
        hinorm = hinorm+t**2
        h(j,istep) = h(j,istep)+t
        call zaxpy(n, -t, v(1,j), 1, v(1,istep+1), 1)
        j = j+1
        if (j.le.istep) goto 18
        
c
c     Compute h(istep+1,istep) and test if reorthogonalization is needed
c
        call zcopy(n,v(1,istep+1),1,vrcom(1,1),1)
        call zcopy(n,v(1,istep+1),1,vrcom(1,2),1)
        flagcom = 7
        revcom = 3
        return
15      continue
        t = vrcom(1,3)
	if ((abs(t)*abs(test)).le.abs(hinorm)
     &                    .and.it.le.3) go to 20	
	t = sqrt(t)
        if (t.eq.czero) then
		ierr = -16
		write(*,*)'Breakdown in the Arnoldi process'
		go to 100
	endif
	h(istep+1,istep) = t
c
c     Normalize v(:,istep+1)
c
	t = cone/t
	call zscal(n, t, v(1,istep+1), 1)
	istep = istep+1
        if (istep.le.iarn) go to 10
	flagcom = 8
c
c     Save h(iarn+1,iarn) for the stopping criterion
c
      hmm = h(iarn+1,iarn)
100   return
      end
