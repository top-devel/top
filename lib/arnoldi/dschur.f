      subroutine dschur(iarn, n, nev, nevfirst, h, z, wr, wi,
     &                  work, y, eig, iord, info)
c
c Purpose
c =======
c
c dschur reorders the Schur basis and the schur form by swapping
c the 1x1 or 2x2 blocks of a.The work is done using the Lapack
c routine dtrexc.f
c
c Parameters:
c ==========
c 
c iarn	-Integer (INPUT)
c       On entry iarn specifies the dimension of the projection
c
c n     -Integer (INPUT)
c       On entry n specifies the dimension of A
c       (The matrix of the eigenvalue problem)
c 
c nev   -Integer (INPUT/OUTPUT)
c       On entry nev specifies the number of wanted eigenvalues.
c       On exit nev specifies the number of computed eigenvalues.
c       (It could be different because of the symmetry of the spectrum)
c
c nevfirst -Integer (INPUT)
c       It is the user's input value of iparam(6)
c
c h     -Double Precision Array of size (iarn+1,iarn+1) (INPUT/OUTPUT)
c	a specifies a matrix in its schur form
c       On exit, the schur form of a is reordered
c
c z     -Double precision Array of size (iarn+1,iarn+1) (INPUT/OUTPUT)
c       On entry, z contains the schur vectors of a
c       On exit z is reordered
c
c wr    -Double precision Array of size (iarn+1) (INPUT/OUTPUT)
c       wr contains the real parts of the eigenvalues of a
c
c wi    -Double precision Array of size (iarn+1) (INPUT/OUTPUT)
c       wi contains the imaginary parts of the eigenvalues of a
c
c work  -Double Precision Array of size (iarn+1)
c       working array
c
c y     -Double Precision Array of size (2*(iarn+1))
c       working array
c
c eig   -character (INPUT)
c       eig='LR' eigenvalues of Largest Real part
c       eig='SR' eigenvalues of Smallest Real part
c       eig='LI' eigenvalues of Largest Imaginary part
c       eig='LM' eigenvalues of Largest Modulii
c       eig='SH' eigenvalues around a given shift
c
c iord  -Integer Array of size (n)
c       iord specifies the new order of the eigenvalues
c
c info -Integer
c      error flag
c  
 
      integer nev, nevfirst, iarn
      integer n, i, j, info, iord(n), ifst, ilst, i1
      double precision h(iarn+1,*), z(iarn+1,*), wr(*),
     &                wi(*), work(*), y(*)
      double precision zero
      parameter (zero = 0.d0)
      character*2 eig
c Functions
      double precision dsqrt
      intrinsic dsqrt
      i = 0
 11   i = i+1
      i1 = i
c
c The Ritz values are sorted following the type of eigen-computation
c defined by the character eig
c
      call dordon(n, iord, wr, wi, iarn, y, eig)
      if (i.eq.1) then
        if ((wi(nevfirst).gt.zero).and.(nev.eq.nevfirst)) then
				nev=nev+1
        elseif ((wi(nevfirst).le.zero).and.((nev-1).eq.nevfirst)) then
				nev=nev-1
	endif 
      endif
      j = iord(i)
      if (wi(i).gt.zero) i = i+1
      if (wi(i).eq.zero.and.eig.eq.'LI') goto 18
      if (i1.eq.j) goto 13
      ifst = j
      ilst = i1
c
c dtrexc swaps two blocks of size 1x1 or 2x2
c (the blocks corresponding to the eigenvalues W(i) et W(iord(i))
c
      call dtrexc('V',iarn,h,iarn+1,z,iarn+1,ifst,ilst,work,info) 
      if (info.ne.0) then
		write(*,*)'Increase iarn and restart the computation'
      		goto 10
      endif
 13   j = 0
 12   j = j+1
      if (h(j+1,j).eq.zero) then
	wr(j) = h(j,j)
        wi(j) = zero
      else
	wr(j) = h(j,j)
	wr(j+1) = h(j,j)
	wi(j) = dsqrt(-h(j+1,j)*h(j,j+1))
	wi(j+1) = -wi(j)
	j = j+1
      endif
      if (j.lt.iarn) goto 12
      if (i.lt.nev) goto 11
 18   call dordon(n, iord, wr, wi, iarn, y, eig)
 10   return
      end
