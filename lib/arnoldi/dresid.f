      subroutine dresid(iarn, nev, wr, wi, s, hmm, normA, res)
c
c Purpose
c =======
c This routine computes the residuals associated with the computed 
c Ritz values:
c
c eta(i) = h(m+1,m) | s_i(m) | /||A||.
c
c The stopping criteria are based on the Arnoldi backward errors.
c
c parameters:
c ==========
c
c  iarn : (input) integer size of the projected problem
c
c  nev  : (input) integer number of wanted eigenvalues
c
c  wr   : (input) real array of size iarn+1 contains the real parts of the Rizt
c       :                                values
c
c  wi   : (input) real array of size iarn+1 contains the imaginary parts of the
c       :                                Ritz values
c
c  s    : (input) real array of size(iarn+1,iarn+1) contains the Schur vectors of 
c       :                                     the upper Hessenberg matrix H
c
c  hmm  : (input) real the coefficient h(m+1,m)
c
c normA : (input) real an estimate of the norm of the matrix A
c
c  res  : (output) real array of size iarn+1 contains the computed residuals
c
      integer iarn, nev
      double precision wr(*), wi(*), s(iarn+1,*), hmm, res(*),
     &                 normA
c
c  Local variables
c
      integer i
      double precision zero
      parameter (zero = 0.d0)
c Functions
      double precision dabs, dsqrt
      intrinsic dabs, dsqrt
      i = 1 
1     if (wi(i).eq.zero) then
c
c  The Ritz value is real
c
	res(i) = hmm*dabs(s(iarn,i))/normA
	i = i+1
	go to 2
      else
c
c  The Ritz value is complex
c
	res(i) = hmm*dsqrt(s(iarn,i)**2+s(iarn,i+1)**2)/normA
	res(i+1) = res(i)
	i = i+2
	go to 2
      endif
2     if (i.le.nev) go to 1
      return
      end	
