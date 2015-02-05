       subroutine deigvec(n, iarn, h, z, v, vec, wi,
     &                    work2, select, z1, z2, iparam)
c
c Purpose 
c =======
c
c This routine computes the right eigenvectors associated with the computed 
c eigenvalues.
c
c Parameters
c ==========
c
c n   -Integer (INPUT)
c     n specifies the dimension of the problem
c     (Dimension of the matrix A)  
c
c iarn -Integer (INPUT)
c      iarn specifies the dimension of the Krylov projection 
c
c h    -Double precision Array of size (iarn+1,iarn+1) (INPUT)
c      The schur form of the hessenberg matrix H
c
c z    -Double precision Array of size (iarn+1,iarn+1) (INPUT)
c      z contains the schur vectors of the hessenberg matrix h
c
c v    -Double precision Array of size (n,iarn+1) (INPUT)	
c      v contains the Arnoldi basis(such as A = V' * H *V)
c
c vec  -Double precision Array of size (n,iarn+1) (OUTPUT)
c      On exit, vec contains the eigenvectors of A corresponding
c      to the eigenvalues computed
c      If w1 is a complex eigenvalue of A then its conjugate w2 too
c      In this case we only need one eigenvector for w1 and w2
c      (Because the two corresponding eigenvectors are conjugated)
c      So,the real and the imaginary parts of this vector will
c      be stocked in two consecutive vectors of vec
c
c wi   -Double precision Array of size (iarn+1) (INPUT)
c      On entry wi specifies the imaginary parts of the
c      computed eigenvalues
c
c work2 -Double precision Array of size (3*(iarn+1))
c       Working array
c
c select Logical Array of size (iarn+1) 
c
c z1   -Double precision Array of size (iarn+1,iarn+1)  
c      Working Array
c
c z2   -Double precision Array of size (iarn+1,iarn+1)
c      Working Array
c
c iparam -Integer Array of size (10)
c        Contains informations on the choices of the user
c        for the computation


      
      integer iarn, n, iparam(10)
      double precision h(iarn+1,iarn+1), z(iarn+1,iarn+1),
     &                 work2(3*(iarn+1)),
     &                  v(n,iarn+1), vec(n,iarn+1), wi(iarn+1),
     &                  z1(iarn+1,iarn+1),
     &                  z2(iarn+1,iarn+1)
      logical select(iarn+1)
c
c    Local variables
c
       integer i, dim, j, ierr, stdout
       double precision zero, one, vl(1)
       parameter (zero = 0.d0, one = 1.d0)
c
c   Select the right eigenvectors to be computed: 1:iparam(6)
c
       do i = 1, iarn
	select(i) = .false.
       enddo
       do i = 1, iparam(6)
	if (wi(i).ge.zero) select(i) = .true.
       enddo
       do j = 1, iarn+1
	do i = 1, iarn+1
		z1(i,j) = zero
	enddo
	z1(j,j) = one
       enddo
c
c   Compute the right eigenvectors of the submatrix T(1:iparam(6),1:iparam(6)) 
c   where T is the quasi-trigular matrix resulting from the Arnoldi
c   Chebyshev iteration
c
       dim = 1
       if (iparam(1).eq.1) then
	do j = 1, iarn
		call dscal(iarn, -one, h(1,j), 1)
	enddo
       endif
       call dtrevc( 'R', 'S', select, iarn, h, iarn+1, vl, dim, z1,
     &                   iarn+1, iarn, dim, work2, ierr)
       if (ierr.lt.0) then
	write(*,*)'Problem computing the right eigenvectors '
	goto 100
       endif
c
c   Perform z2 <-- z*z1
c
       call dgemm('N', 'N', iarn, iparam(6), iarn, one, z, iarn+1,
     &             z1, iarn+1, zero, z2, iarn+1)
c
c   Compute the right eigenvectors corresponding to the computed 
c   eigenvalues of the original matrix A: y<-v*z2
c
       call dgemm('N', 'N', n, iparam(6), iarn, one, v, n, z2, iarn+1, 
     &             zero, vec, n)

c
c   If wanted 
c   save the right eigenvectors in the file eigenvec
c   Else exit
c
       if (iparam(9).eq.0) goto 100

       open(33, file = 'eigenvec')
       i = 0
1      i = i+1
       if (wi(i).lt.zero) then
	do j = 1, n
		write(33,10)vec(j,i), -vec(j,i+1)
	enddo
	i = i+1
       elseif (wi(i).gt.zero) then
	do j = 1, n
                write(33,10)vec(j,i), vec(j,i+1)
        enddo
        i = i+1
	elseif (wi(i).eq.zero) then
	do j = 1, n
                write(33,11)vec(j,i)
        enddo
       endif
       if (i.lt.iparam(6)) go to 1
10     format( 2e24.16)
11     format( e24.16)
       close(33, status = 'keep')
       if (iparam(7).ne.0) write(*,16)
16     format('The eigenvectors are in the file eigenvec')
100    return
       end
