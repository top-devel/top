      subroutine zeigvec(n, iarn, h, z, v, vec, w,
     &                    work2, work2real, select, z1, z2, iparam)
c
c Purpose
c =======
c
c This routine computes the right eigenvectors associated with the wanted 
c eigenvalues.
c
c Parameters
c ==========
c
c n	-Integer (INPUT)
c	n specifies the dimension of the problem
c	(Dimension of the matrix A)
c
c iarn	-Integer (INPUT)
c      iarn specifies the size of the projection
c
c h 	-Complex Array of size (iarn+1,iarn+1) (INPUT)
c      The schur form of the hessenberg matrix H
c
c z    -Complex Array of size (iarn+1, iarn+1) (INPUT)
c      z contains the schur vectors of the hessenberg matrix h
c
c v    -Complex Array of size (n,iarn+1) (INPUT)
c      v contains the Arnoldi basis(such as A = V' * H *V)
c
c vec  -Complex Array of size (n,iarn+1) (OUTPUT)
c      On exit, vec contains the eigenvectors of A corresponding
c      to the eigenvalues computed
c    
c      
c w    -Complex Array of size (iarn+1)   (INPUT)
c      On entry w specifies the computed eigenvalues
c
c work2 -Complex Array of size (3*(iarn+1))
c       Working Array
c
c work2real -Double Precission Array of size (iarn+1)
c	    Working Array
c
c select -Logical Array of size (iarn+1)
c
c z1     -Complex Array of size (iarn+1, iarn+1)
c        Working Array
c
c z2     -Complex Array of size (iarn+1, iarn+1)
c
c iparam -Integer Array of size (10)
c        Contains informations on the choices of the user
c        for the computation
c

c
c Local Variables
c
       integer iarn, n, iparam(10)
       complex *16  h(iarn+1,iarn+1), z(iarn+1,iarn+1),
     &             work2(3*(iarn+1)), v(n,iarn+1),
     &             vec(n,iarn+1), w(iarn+1), z1(iarn+1,iarn+1),
     &             z2(iarn+1,iarn+1)
       double precision work2real(iarn+1)
       logical select(iarn+1)
c
c    Local variables
c
       integer i, dim, j, ierr
       double precision zero, one, vl(1)
       parameter (zero = 0.d0, one = 1.d0)
       complex *16 czero, cone
       parameter (czero = (0.d0,0.d0), cone = (1.d0,0.d0))
c
c   Select the right eigenvectors to be computed: 1:iparam(6)
c
       do j = 1, iarn+1
	do i = 1, iarn+1
		z1(i,j) = czero
	enddo
	z1(j,j) = cone
       enddo
 
       do i= 1, iarn
        select(i)= .true.
       enddo

c
c   Compute the right eigenvectors of the submatrix T(1:iparam(6),1:iparam(6)) 
c   where T is the quasi-trigular matrix resulting from the Arnoldi
c   Chebyshev iteration
c
       dim = 1
       if (iparam(1).eq.1) then
	do j = 1, iarn
		call zscal(iarn, -cone, h(1,j), 1)
	enddo
       endif
       call ztrevc( 'R', 'S', select, iarn, h, iarn+1, vl, dim, z1,
     &             iarn+1 , iarn, dim, work2, work2real, ierr)
       if (ierr.lt.0) then
	write(*,*)'Problem computing the right eigenvectors '
	goto 100
       endif
c
c   Perform z2 <-- z*z1
c
       call zgemm('N', 'N', iarn, iparam(6), iarn, cone, z, iarn+1, z1, 
     &            iarn+1, czero, z2, iarn+1)
c
c   Compute the right eigenvectors corresponding to the wanted 
c   eigenvalues of the original matrix A: y<-v(:,iparam(6)) *z2
c
       call zgemm('N', 'N', n, iparam(6), iarn, cone, v, n, z2, iarn+1, 
     &             czero, vec, n)
c
c   If wanted
c   save the right eigenvectors in the file eigenvec
c   else exit
c
       if (iparam(9).eq.0) goto 100

       open(33, file = 'eigenvec',status='unknown')
       i = 0
1      i = i+1
       do j = 1, n
          write(33,10)vec(j,i)
       enddo
       if (i.lt.iparam(6)) go to 1
10     format( 2e24.16)
11     format( e24.16)
12     format( A8)
       close(33, status = 'keep')
       if (iparam(7).ne.0) write(*,16)
16     format('The eigenvectors are in the file eigenvec')
100    return
       end
