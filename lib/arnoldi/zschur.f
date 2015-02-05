      subroutine zschur(iarn, n, nev, h, z, w,
     &                  work, yr, yc, eig, iord, info)
c
c Purpose
c =======
c This routine computes the Ritz values using the lapack
c routine ztrexc and reordered them using the routine
c cordon
c
c Parameters
c ==========
c
c iarn : (input) integer The size of the projection
c
c n : (input) integer the leading dimension of the matrix A 
c
c nev : (input) integer the number of wanted eigenvalues
c
c h : (input/output) complex array of size (iarn+1,iarn+1)
c     On entry a matrix on its schur form
c     On exit the reordered schur form
c
c z : (input/output) the matrix of the schur vectors
c
c w : (input/output) complex array of size (n)
c     contains the Ritz Values
c  
c work : (input/output) complex array of size (iarn+1)
c                       working array
c
c yr : (input/output) real array of size (iarn+1)
c                     working array
c
c yc : (input/output) complex array of size (iarn+1)
c                     working array
c
c eig :(input) character determines the type of eigen-computation 
c
c iord : (output) integer array of size iarn+1 contains the new reordering
c
c info : (output) integer the error flag
c
      integer nev
      integer iarn, n, i, j, info, iord(*), ifst, ilst
      double precision yr(*)
      complex *16 h(iarn+1,*), z(iarn+1,*), w(*),
     &            work(*), yc(*) 
      complex *16 zero
      parameter (zero = (0.d0,0.d0))
      character*2 eig

c ...Functions
      double precision dimag
      intrinsic dimag

      i = 0
 11   i = i+1
      call zordon(n, iord, w, iarn, yr, yc, eig)
      j = iord(i)
      if (dimag(w(i)).eq.zero.and.eig.eq.'LI') go to 18
      if (i.eq.j) goto 13
      ifst = j
      ilst = i
      call ztrexc( 'V', iarn, h, iarn+1, z, iarn+1, ifst,
     &              ilst, work, info) 
      if (info.ne.0) then
		write(*,*)'Increase iarn and restart the computation'
      		goto 10
      endif
 13   j = 0
 12   j = j+1
      w(j) = h(j,j)
      if (j.lt.iarn) goto 12
      if (i.lt.nev) goto 11
 18   call zordon(n, iord, w, iarn, yr, yc, eig)
 10   return
      end
