      subroutine zresid( iarn, nev, w, s, hmm, normA, res)
c
c This routine computes the residuals associated with the computed 
c Ritz values:
c
c eta(i) = h(m+1,m) | s_i(m) | /||A||.
c
c The stopping criteria are based on the Arnoldi backward errors.
c
c parameters:
c
c  iarn : (input) integer size of the projected problem
c
c  nev  : (input) integer number of wanted eigenvalues
c
c  w   : (input) complex array of size iarn+1 contains the Ritz values
c
c  s    : (input) complex array of size(iarn+1,iarn+1) contains the Schur vectors of 
c       :                                     the upper Hessenberg matrix H
c
c  hmm  : (input) complex the coefficient h(m+1,m)
c
c  normA : (input) real an estimate of the norm of the matrix A
c
c  res  : (output) complex array of size iarn+1 contains the computed residuals
c
      integer iarn, nev
      double precision normA
      complex *16 w(*), s(iarn+1,*), hmm, res(*)
c
c  Local variables
c
      integer i
      double precision zero, residu
      parameter (zero = 0.d0)
c ...Functions...
      complex *16 dcmplx
      intrinsic dcmplx
c
      i=0
1     i=i+1
      residu=(abs(hmm)*abs(s(iarn,i)))/normA
      res(i)=dcmplx(residu,0.d0)
2     if (i.lt.nev) go to 1
      return
      end	
