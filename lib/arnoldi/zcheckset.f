      subroutine zcheckset(n, iarn, deg, degmax, iparam,
     &                     tol, normA, shift, u, ierr, eig, 
     &                     nor, che, stdout, revcom)
c
c The routine checks the entries of the computation. It sets 
c the default values depending of the array iparam.
c
c parameters:
c
c   n   : (input) integer the size of the matrix.
c
c iarn  : (input/output) integer the size of the projection.
c
c deg   : (input/output) integer the degree pf the Chebyshev polynomial.
c
c degmax: (output) integer the maximal allowed values of deg.
c
c iparam: (input) integer array of size 10 contains the options set by the 
c       :         user:
c       :         iparam(1): the type of eigenproblem to be solve
c       :                    0-> eigenvalues of largest real part
c       :                    1-> eigenvalues of smallest real part
c       :                    2-> eigenvalues of largest imaginary part
c       :                    3-> eigenvalues of largest modulii
c       :                    4-> eigenvalues around a given shift
c       :         iparam(2): how the computations of ||A|| is done
c       :                    0-> ||A|| is computed by n matrix-vector
c       :                        products
c       :                    1-> ||A|| is given
c       :                    2-> ||A|| is computed using Higham's 
c       :                        modification of Hager's algorithm 
c       :                        (requires y<--Ax and y<--A^*x)
c       :         iparam(3): the degree of the Chebyshev polynomial
c       :                    0-> deg is dynamically set by the code
c       :                    1-> deg is given and fixed during all the 
c       :                        computation
c       :         iparam(4): the tolerance for the computation
c       :                    0-> tol is set to the n*machine precision
c       :                    1-> tol is given
c       :         iparam(5): the starting vector
c       :                    0-> the starting vector is a random vector
c       :                    1-> the starting vector is given
c
c       :         iparam(6): the number of wanted eigenvalues
c       :         iparam(10):Number maximum of iterations
c
c tol   : (input/output) real the tolerance (not referenced if iparam(4)=1)
c
c normA : (input) real the estimate of ||A|| (not referenced if iparam(2)=1,2)
c
c shift : (input) complex the given shift if iparam(1) = 4
c
c  u    : (input/output) complex array of size n contains the starting vector
c       :                (not referenced if iparam(5)=1)
c
c ierr  : (output) integer the error flag
c       :                  -1-> n<0
c       :                  -2-> iparam(6)<0
c       :                  -3-> iparam(1) different than 0 1 2 3 or 4
c       :                  -4-> iparam(2) different than 0 1 or 2
c       :                  -5-> iparam(3) different than 0 or 1
c       :                  -6-> iparam(4) different than 0 or 1
c       :                  -7-> iparam(5) different than 0 or 1
c       :                  -8-> ||A||<0
c       :                  -9-> iarn<0
c       :                 -10-> deg<0
c       :                 -11-> tol<0
c       :                 -12-> iarn>n
c       :                 -13-> iparam(6)>iarn
c
c eig   : (output) character define the type of eigenproblem to be solved
c
c nor   : (output) character define the way the norm estimate is done
c
c che   : (output) character define the way the degree of the degree of 
c       :          the Chebyshev polynomial is set
c
c stdout: (output) Integer defines the standard out for the error messages
c                  and for the computation progress
c
c revcom: (output) integer reverse comunication flag
c       :                on exit normal exit if revcom = 0
c       :                        if ierr<0 revcom = -1
c
      integer n, iarn, deg, ierr, iparam(10), degmax, 
     &        revcom,i,stdout
      double precision tol, normA, dlamch, tol1
      complex *16 u(*),shift
      character*2 eig, nor, che
      character*23 string
      double precision zero
      parameter (zero = 0.d0)
c
      ierr = 0
      if (n.le.0)                                  ierr = -1
      if (iparam(6).le.0)                                ierr = -2
      if (iparam(1).ne.0.and.iparam(1).ne.1.and.iparam(1).ne.
     &    2.and.iparam(1).ne.3.and.iparam(1).ne.4) ierr = -3
      if (iparam(2).ne.0.and.iparam(2).ne.1.and.iparam(2).ne.
     &    2)                                       ierr = -4
      if (iparam(3).ne.0.and.iparam(3).ne.1)       ierr = -5
      if (iparam(4).ne.0.and.iparam(4).ne.1)       ierr = -6
      if (iparam(5).ne.0.and.iparam(5).ne.1)       ierr = -7
      if (iparam(2).eq.1) then
	if(normA.le.zero)        ierr = -8
      endif
      if (iarn.le.0)                               ierr = -9
      if (iparam(3).eq.1) then
	if (deg.le.0)             ierr = -10
      endif
      if (iparam(4).eq.1) then
	if (tol.le.zero)          ierr = -11
      endif
      if (iarn.gt.n)                               ierr = -12
      if (iparam(6).gt.iarn)                             ierr = -13
      if (ierr.lt.0) go to 100
c
c    No error occurs, set parameters and defaults if needed
c
c
c 1) Set the standard out
c
      if (iparam(7).eq.0) then
                    open(11, file = 'Exec',status='unknown')
                    stdout=11
      else
                               stdout = 6
      endif
 
c
c 2) Set eigenproblem to be solved:
c
      if (iparam(1).eq.0) then
	eig = 'LR'
      elseif (iparam(1).eq.1) then
	eig = 'SR'
      elseif (iparam(1).eq.2) then
	eig = 'LI'
      elseif (iparam(1).eq.3) then
	eig = 'LM'
      elseif (iparam(1).eq.4) then
	eig = 'SH'
      endif
c
c 3) Set which computations of ||A|| is done:
c
      if (iparam(2).eq.0) then
	nor = 'PM'
      elseif (iparam(2).eq.1) then
	nor = 'GI'
      elseif (iparam(2).eq.2) then
	nor = 'HH'
      endif
c
c 4) Set the degree of the Chebyshev polynomial:
c
      degmax = min(200, int(n/2))
      if (iparam(3).eq.0) then
	che = 'DY'
      elseif (iparam(3).eq.1) then
	che = 'GI'
        if (deg.gt.degmax) then
		write(stdout,*)'Warning'
		write(stdout,*)'-------'
		write(stdout,*)'Chebyshev polynomial degree is set to',degmax
		write(stdout,*)'-------'
		deg = degmax
	endif
      endif
c
c 5) Set the tolerance for the computation:
c
      tol1 = dble(n)*dlamch('e')
      if (iparam(4).eq.0) then
	tol = tol1
      elseif (iparam(4).eq.1.and.tol.lt.tol1) then
	tol = tol1
      endif
c
c 6) Set the starting vector:
c
      if (iparam(5).eq.0) call zvrand(n, u)
      if (iparam(5).eq.1) then
c           do i=1,n
c               u(i)=dcmplx(dble(1),dble(0))
c           enddo
      print*,' initial vector in zcheckset:'
      print*,(u(i),i=1,10)
      endif

c
c    Display parameters and options of the computation
c
      write(stdout,*)'-----------------------------------------------'
c     write(7,*)'-----------------------------------------------'
      write(stdout,*)'Size of the eigenproblem: ',n
c     write(7,*)'Size of the eigenproblem: ',n
      write(stdout,*)'Size of the projection given to: ',iarn
c     write(7,*)'Size of the projection given to: ',iarn
      write(stdout,*)'Number of wanted eigenvalues: ',iparam(6)
c     write(7,*)'Number of wanted eigenvalues: ',iparam(6)
      if (eig.eq.'LR') then
        string = 'Largest real part'
      elseif (eig.eq.'SR') then
        string = 'Smallest real part'
      elseif (eig.eq.'LI') then
        string = 'Largest imaginary part'
      elseif (eig.eq.'LM') then
	string = 'Largest modulii'
      endif
      if (iparam(1).lt.4) then
      	write(stdout,*)'Type of resolution:  ',string
c     	write(7,*)'Type of resolution:  ',string
      else
      	write(stdout,44)shift
c     	write(7,44)shift
      endif
c 44   format('Type of resolution: Around the shift',d12.3)
 44   format('Type of resolution: Around the shift',2e12.3)
      if (nor.eq.'GI') then
	write(stdout,*)'||A|| is equal to: ', normA
c       write(7,*)'||A|| is equal to: ', normA
      elseif (nor.eq.'PM') then
	write(stdout,*)'||A|| estimated by n matrix-vector products'
c       write(7,*)'||A|| estimated by n matrix-vector products'
      elseif (nor.eq.'HH') then 
	write(stdout,*)'||A|| estimated by Hager''s algorithm'
c       write(7,*)'||A|| estimated by Hager''s algorithm'
      endif
      if (che.eq.'GI') then
	write(stdout,*)'Deg. of Chebyshev poly. is equal to: ',deg
c        write(7,*)'Deg. of Chebyshev poly. is equal to: ',deg
      elseif (che.eq.'DY') then
	write(stdout,*)'Deg. of Chebyshev poly. is dynamically set'
c       write(7,*)'Deg. of Chebyshev poly. is dynamically set'
      endif
      write(stdout,45)tol
c     write(7,45)tol
 45   format('Tolerance is equal to:',d12.3)
      if (iparam(5).eq.1) then
	write(stdout,*)'Starting vector is given'
c       write(7,*)'Starting vector is given'
      elseif (iparam(5).eq.0) then
	write(stdout,*)'Starting vector is set to a random one'
c       write(7,*)'Starting vector is set to a random one'
      endif
c     write(7,*)'Number maximum of iterations'
      write(stdout,*)'Number maximum of iterations',iparam(10)
      write(stdout,*)'-----------------------------------------------'
c     write(7,*)'-----------------------------------------------'
      revcom = 1
      go to 200
100   revcom = -1
200   return
      end
