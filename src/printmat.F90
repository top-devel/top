#include "config.h"
       subroutine printmat(a,lda,nz,kl,ku,filename)
       implicit none
       integer nz,kl,ku,i,j,lda
       double precision a(lda,nz)
       character*(*) filename

       open(unit = 17,file=filename,status='unknown')

       write(17,*) nz,nz
       do j=1,nz
         do i=max(1,j-ku),min(nz,j+kl)
           if (a(ku+kl+1+i-j,j).ne.0d0) then
             write(17,102) i,j,a(ku+kl+1+i-j,j)
           endif
         enddo
       enddo
102    format(I5,2X,I5,2(2X,1pe22.15))

       close(17)
       end
!--------------------*-*-*-*-*---------------*-*-*-*-*------------*-*-*--
       subroutine printmat_full(a,M,N,filename)
       implicit none
       integer i,j,M,N
       double precision a(M,N)
       character*(*) filename

       open(unit = 17,file=filename,status='unknown')

       write(17,*) M,N
       do i=1,M
         do j=1,N
           if (a(i,j).ne.0d0) write(17,102) i,j,a(i,j)
           !write(17,102) i,j,a(i,j)
         enddo
       enddo
102    format(I5,2X,I5,1(2X,1pe22.15))
103    format(I5,2X,I5,1(2X,I3))

       close(17)
       end
