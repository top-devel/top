!--------------------------------------------------------------
! Routine numerical recipes computing the weights and points for
! a Gauss-Legendre integration.
!--------------------------------------------------------------

       subroutine gauleg(x1,x2,x,w,n)

       implicit none
       integer, intent(in) :: n
       double precision, intent(in) :: x1, x2
       double precision, intent(inout) :: x(n), w(n)
       double precision, parameter :: eps=3.d-15
       double precision, parameter :: pi=3.14159265358979d0
       integer m, i, j
       double precision xm, xl, p1, p2, p3, pp, z, z1

       m=(n+1)/2
       xm=0.5d0*(x2+x1)
       xl=0.5d0*(x2-x1)
       do i=1,m
         z=cos(pi*(i-0.25d0)/(n+0.5d0))
         do
           p1=1d0
           p2=0d0
           do j=1,n
             p3=p2
             p2=p1
             p1=((2d0*dble(j)-1d0)*z*p2-(dble(j)-1d0)*p3)/dble(j)
           enddo
           pp=n*(z*p1-p2)/(z*z-1d0)
           z1=z
           z=z1-p1/pp
           if (abs(z-z1).le.eps) exit
         enddo
         x(i)=xm-xl*z
         x(n+1-i)=xm+xl*z
         w(i)=2d0*xl/((1d0-z*z)*pp*pp)
         w(n+1-i)=w(i)
       enddo
       return
       end subroutine gauleg
