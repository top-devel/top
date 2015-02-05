      subroutine zarncheb(n, iarn, deg, iparam, tol, shift,
     &                   normA, u, vrcom, iord, v, h, z, 
     &                   w, work, work1, workc, yr, yc, revcom, info)
c
c Purpose
c =======
c
c These routine compute the few eigenvalues of alarge
c complex non-symetric matrix A using a reverse communication
c strategy for the matrix-vector products
c
c Parameters
c ==========
c n	-Integer (INPUT)
c	n specifies the number of rows and columns
c       of the matrix A
c
c iarn  -Integer (INPUT)
c       On entry iarn specifies the size of the projection
c       for the krylov method
c
c deg   -Integer (INPUT/OUTPUT)
c       deg specifies the degree of the chebychev polynomial
c       According to the choice of the user(see iparam),
c       It could be an INPUT parameter or it could be set
c       dynamically by the program(INPUT/OUTPUT)
c
c iparam -Integer array of size 10 (INPUT)
c        This array specifies the informations and the choices
c        of the user.(See users guide or README file for more
c        information)
c
c tol 	-Double precision (INPUT/OUTPUT)
c        tol specifies the tolerance for the computation
c        According to the choice of the user(see iparam),
c        It could be an INPUT parameter or it could be set
c        dynamically by the program (INPUT/OUTPUT)
c
c shift -Complex (INPUT)
c       shift is the parameter used to compute the eigenvalue
c       closed to a given value(shift and invert problem)
c       (Only referenced if iparam(1)=4)
c
c normA -Double precision (INPUT/OUTPUT)
c        This is the estimated norm of the matrix A
c        iparam(2) indicates how the computation of normA
c        is done.
c
c u	-Complex Array of size (n) (INPUT/OUTPUT)
c       Only referenced if iparam(5)=1
c
c vrcom -Double precision array of size (n,4) (INPUT/OUTPUT)
c       This array contains three vectors and is used for the
c       reverse communication strategy.
c       If revcom = 1 the user is asked to perform :
c                 vrcom(:,3) <- A * vrcom(:,1)
c       If revcom = 2 the user is asked to perform :
c                 vrcom(:,4) <- A' * vrcom(:,3)
c       If revcom = 3 the user is asked to perform :
c                 vrcom(:,3) <- vrcom(:,1) * vrcom(:,2)
c
c
c iord 	-Integer Array of size n
c      	working array
c
c v    	-Complex array of size(n,iarn+1) (OUTPUT)
c      	v is an orthonormal system such as
c      	A = V' * H * V
c      	where H is a hessenberg matrix of size (iarn+1,iarn+1)
c
c h    	-Complex Array of size (iarn+1,iarn+1) (OUTPUT)
c      	H is a hessenberg matrix.
c      	it is the projection of A onto the Krylov subspace
c
c z    	-Complex Array of size (iarn+1,iarn+1) (OUTPUT)
c      	On exit z contains the schur vectors of the hessenberg
c      	matrix h
c
c w   	-Complex Array of size (iarn+1) (OUTPUT)
c      	On exit wr contains the real parts of the computed
c      	Eigenvalues
c
c work 	-Complex Array of size (iarn+1)
c      	Working Array
c
c work1 -Complex Array of size (iarn+1)
c       Working Array
c
c workc -Complex Array of size (n,2)
c       Working Array
c
c yr	-Double precision Array of size (iarn+1)
c	Working Array
c
c yc    -Complex Array of size (iarn+1)
c       Working Array
c
c revcom -Integer (INPUT/OUTPUT)
c        revcom is the reverse communication parameter.
c        see the users guide and the files README
c
c info  -Integer (OUTPUT)
c       integer error flag. info = 0 on a normal exit
c                           info < 0 if there is an error 
c
 
c
c ARGUMENTS
c
c ...Scalar Arguments...
      integer n, iarn, deg, revcom, info      
      double precision normA, tol
      complex *16 shift
c ...Array Arguments
      integer iparam(10), iord(n)
      double precision yr(iarn+1)
      complex *16 u(n), vrcom(n,4),
     &                 v(n,iarn+1), h(iarn+1,iarn+1),
     &                 z(iarn+1,iarn+1), w(iarn+1),
     &                 work(iarn+1), yc(iarn+1), work1(iarn+1),
     &                 workc(n,2)


c LOCAL DECLARATION
c ...Parameters...
      complex *16 cone, czero
      parameter (czero=(0.d0,0.d0),cone=(1.d0,0.d0))
c ...Local Scalars...
      integer degmax, flagcom, i, j , iconv, iter, icycle,
     &        ierrell, stdout
      double precision resmax, c
      complex *16 hmm, d, ea1
      character*2 eig, nor, che
c ...Save Statements...
      DATA flagcom /-1/
      save degmax, flagcom, iter, icycle, stdout
      save hmm, c, d, ea1
      save eig, nor, che
      
c Functions
      double precision zabs
      intrinsic zabs      

c
c degmax -Integer
c        Maximum degree for the Chebychev polynomial
c
c flagcom -Integer
c         Used for the reverse communication strategy
c         It allows the program to know where the execution
c         has been stopped
c
c iter -Integer
c      iter specifies the number of iterations of the Arnoldi
c      Algorithm that have been executed
c
c icycle -Integer
c      icycle specifies the number cycles that have been
c      executed in the Chebychev acceleration
c
c
c ierrell -Integer 
c         Error flag for the routines that compute the ellipse
c
c stdout -Integer
c        It defines the standard out for the error messages
c        and the computation progress
c        6 is the screen
c        11 is the file Exec
c
c hmm   -Complex
c       hmm specifies the coefficient of the Hessenberg matrix
c       used to compute the backward error
c
c c     -Double precision
c       Focal Distance of the chebychev ellipse
c
c ea1   -Complex
c       Parameter of the chebychev acceleration
c
c d     -Complex
c       Center of the chebychev ellipse

      if (flagcom.lt.3) then
c
c  if flagcom = -1: Initial call and error checking and set default if necessary
c  if 0 < flagcom <3 : Compute ||A||
c
      	if (flagcom.eq.-1) then 
   		info = 0
                open (7, file = 'result',status='unknown')
 		open (9, file = 'parameter',status='unknown')
		call zcheckset(n, iarn, deg,
     &                    degmax, iparam, tol, normA, shift, u,
     &                    info, eig, nor, che, stdout, revcom)
       	 	if (info.lt.0) goto 100
		flagcom = 0
        endif
c
c  Computation of the estimate of ||A||
c
      	call zestnor(n, nor, normA, iord, flagcom, 
     &                       vrcom, revcom)
        if (flagcom.lt.3) return
c
c The computation starts here
c
        iter = 0
600     iter = iter+1
        if (iter.gt.iparam(10)) then
		info=-17
		goto 100
        endif

	if (iter.eq.1) write(9,27)normA
27   	format('The estimated of ||A|| is equal to: ',d12.3,/)
        go to 700
      elseif (flagcom.ge.3.and.flagcom.le.7) then
c
c  Computation of the basis of the Krylov subspace v and of the
c  upper Hessenberg matrix h
c
c  If the user wants to compute the eigenvalues having the
c  smaller real parts then we work with -A
c
 700    if (flagcom.eq.4.and.eig.eq.'SR')
     &                call zscal(n, -cone, vrcom(1,3), 1)
        call zhessen(n, iarn, v, h, u, hmm, flagcom, vrcom,
     &               revcom, info)
        if (info.lt.0) then 
			info = -14
                      	goto 100
        endif
        if (flagcom.lt.8) return
c
c  Computation of the Schur vectors and the Schur form of the matrix h
c  using the Lapack routine zhseqr.f
c
	do j = 1, iarn+1
                do i = 1, iarn+1
                        z(i,j) = czero
                enddo
                z(j,j) = cone
        enddo
        i = 1
        j = iarn
        call zhseqr( 'S', 'I', iarn, i, j, h, iarn+1, w, z,
     &             iarn+1, work, iarn+1, info)
        if (info.ne.0) then
		info = -15	
                goto 100
        endif
c
c  Reordering of the Schur basis and the Schur form by swapping the
c  1x1 or 2x2 blocks of h. The work is done using the Lapack routine
c  dtrexc.f
c
        call zschur(iarn, n, iparam(6), h, z, w,
     &              work, yr, yc, eig, iord, info)
	if (info.ne.0) then 
		info = -16
		goto 100
	endif
c
c  Compute the corresponding residual based on the backward error and
c  Test for convergence
c
        call zresid(iarn, iparam(6), w, z, hmm, normA, work)
	iconv = 0
        resmax = zabs(work(1))
        do i = 1, iparam(6)
                if (zabs(work(i)).lt.tol) iconv = iconv+1
                if (zabs(work(i)).gt.resmax) resmax = zabs(work(i))
        enddo
c
c  Write history of the computation on the screen
c
        write(stdout,20)iter, resmax
20      format('Maximum residual at iteration No:',i4,' -> ',d12.3)
        write(stdout,21)iconv
21      format('Total number of converged eigenvalues : ',i4)
c
c  Get back to the original eigenvalues if the user wants to
c  compute either the eigenvalues having the smaller real part
c  (work done using -A) or the eigenvalues around the given
c  shift (lambda = shift+1/mu).
c
        if (eig.eq.'SR') then
                do i = 1, iparam(6)
                        w(i) = -w(i)
                enddo
        elseif (eig.eq.'SH') then
                do i = 1, iparam(6)
			w(i)=shift+cone/w(i)
                enddo
        endif
c
c   Write eigenvalues and residuals in file result
c
        write(7,22)iter
22      format(' ',' Iteration No',i4,/)
        write(7,23)(dreal(w(j)),dimag(w(j))
     &              ,zabs(work(j)),j = 1,iparam(6))
23      format(' ','    re(lambda) ',
     &     '         *       im(lambda)      *      res',/,(/3e24.16))
        write(7,24)
24      format(' ',//)
c
c   Convergence for the wanted eigenvalues is achieved, then Exit
c
        if (iconv.ge.iparam(6)) then
	   if (iparam(8).ne.0) then
		open(44, file = 'eigenval',status='unknown')
		write(44,*)'Dimension and number of wanted eigenvalues'
		write(44,*)n, iparam(6)
		write(44,*)'||A||'
		write(44,30)normA
30              format(e24.16)
                write(44,*)'Computed eigenvalues'
		i = 0
45              i = i+1
		write(44,28) dreal(w(i)),dimag(w(i))
		if (i.lt.iparam(6)) go to 45	
        	close(44, status = 'keep')
28      	format(2e24.16)
	   endif 
		goto 100
	endif
c
c   Convergence is not achieved:
c   Recompute the current eigenvalues if eig = 'SR' or eig = 'SH' in
c   order to compute the Chebyshev ellipse
c
        if (eig.eq.'SR') then
                do i = 1, iparam(6)
                        w(i) = -w(i)
                enddo
        elseif (eig.eq.'SH') then
                do i = 1, iparam(6)
			w(i)=cone/(w(i)-shift)
                enddo
        endif
c
c   Compute the ellipse enclosing the unwanted part of the spectrum
c   and the degree of Chebyshev polynomial is iparam(3) = 0
c
        if (eig.eq.'LR'.or.eig.eq.'SR') then
                call zellipr(iarn, w, iparam(6), d, c, ea1,
     &                       work1, deg, degmax, che, ierrell)
        elseif (eig.eq.'LI') then
                call zellipi(iarn, w, iparam(6), d, c, ea1,
     &                       work1, deg, degmax, che,ierrell)
        elseif (eig.eq.'LM'.or.eig.eq.'SH') then
		call zellipm(iarn, w, iparam(6), d, c, ea1,
     &                  n, iord, deg, degmax, che, ierrell)
        endif
c
c   Write parameters of Chebyshev acceleration in file parameter
c
        write(9,22)iter
        if (ierrell.lt.0) then
                write(9,*)'Restart with Arnoldi vectors'
        else
                write(9,25)dreal(d), dimag(d), c, 
     &                     dreal(ea1),dimag(ea1)
        endif
25      format(' Centre(real part): ',d12.4,
     &         ' Centre(imaginary part): ',d12.4,
     &         ' Focus distance: ',d12.4,
     &         ' ea1(real part): ',d12.4,
     &         ' ea1(imaginary part)',d12.4)
c
c  Compute starting vector either for Chebyshev accereration or for 
c  restarting using Arnoldi vectors
c
	call zscal(n, czero, u, 1)
        j = 0
300     j = j+1
        call zgemv('N', n, iarn, cone, v, n, z(1,j), 1,
     &                     czero, vrcom(1,1), 1)
        call zaxpy(n, work(j), vrcom(1,1), 1, u, 1)
        if (j.lt.iparam(6)) goto 300
c
c   if ierrell<0 then no ellipse can be built. Restart using a weighted 
c   sum of the Arnoldi vectors as starting vectors for the next 
c   Arnoldi step
c
	if (ierrell.lt.0) then
		ierrell = 0
		flagcom = 3
		go to 600
	endif
c
c   else Perform Chebysev acceleration
c
	call zcopy(n, u, 1, vrcom(1,1), 1)
        icycle = 0
        call zcheby(n, deg, c, ea1, d, workc, flagcom, vrcom,
     &              icycle, revcom)
        return
      elseif (flagcom.eq.8) then
        if (eig.eq.'SR'.and.revcom.ne.3) call zscal(n, -cone,
     &                                        vrcom(1,3), 1)
        call zcheby(n, deg, c, ea1, d, workc, flagcom, vrcom,
     &              icycle, revcom)
        if (flagcom.eq.8) return
c
c   Compute the next starting vector and restart the Arnoldi computation
c
        call zcopy(n, workc(1,2), 1, u, 1)
        flagcom = 3
        goto 600
      endif
100   if (info.ne.0) then
	revcom = -1
        flagcom = -1
        if (info.eq.-1) then
        write(stdout,*)'the size n of the matrix must be non negative' 
        write(stdout,*)'change the parameters and restart computation'
        elseif (info.eq.-2) then
        write(stdout,*)'iparam(6) must be nonnegative'
        write(stdout,*)'change the parameters and restart computation'
        elseif (info.eq.-3) then
        write(stdout,*)'iparam(1) must be 0 1 2 3 or 4'
        write(stdout,*)'change the parameters and restart computation'
        elseif (info.eq.-4)  then
        write(stdout,*)'iparam(2) must be 0 1 or 2'
        write(stdout,*)'change the parameters and restart computation'
        elseif (info.eq.-5) then 
        write(stdout,*)'iparam(3) must be 0 or 1'
        write(stdout,*)'change the parameters and restart computation'
        elseif (info.eq.-6) then
        write(stdout,*)'iparam(4) must be 0 or 1'
        write(stdout,*)'change the parameters and restart computation'
        elseif (info.eq.-7) then
        write(stdout,*)'iparam(5) must be 0 or 1'
        write(stdout,*)'change the parameters and restart computation'
        elseif (info.eq.-8) then
        write(stdout,*)'the given value of ||A|| must be non negative'
        write(stdout,*)'change the parameters and restart computation'
        elseif (info.eq.-9) then
        write(stdout,*)'the projection size iarn must be non negative'
        write(stdout,*)'change the parameters and restart computation'
        elseif (info.eq.-10) then
	write(stdout,*)'the Chebyshev degree deg must be non negative' 
        write(stdout,*)'change the parameters and restart computation'
        elseif (info.eq.-11) then 
	write(stdout,*)'the tolerance tol must be non negative'
        write(stdout,*)'change the parameters and restart computation'
	elseif (info.eq.-12) then
        write(stdout,*)'iarn must be less than n'
        write(stdout,*)'change the parameters and restart computation'
        elseif (info.eq.-13) then
        write(stdout,*)'iparam(6) must be less than iarn'
        write(stdout,*)'change the parameters and restart computation'
        elseif (info.eq.-14) then 
        write(stdout,*)'Breakdown in the Arnoldi process'
        elseif (info.eq.-15) then
        write(stdout,*)'Problem using the Lapack routine zhseqr'
        elseif (info.eq.-16) then
        write(stdout,*)'Problem swapping Schur form'
        elseif (info.eq.-17) then
	write(stdout,*)'Number maximum of iterations reached'
        endif
      else
	revcom = -99
        flagcom = -1
        info = 0
        write(stdout,*)'Normal exit'
        write(stdout,*)'The eigenvalues are in the file eigenval'
        write(stdout,*)'The computation parameters are in'//
     &' the file parameter'
        write(stdout,*)'The computation history is in the file result'
      endif
      close(7, status = 'keep')
      close(9, status = 'keep')
      return
      end
