#include "config.h"
       module Chebyshev

       integer, save :: nr_cheb, nr_cheb1, nr_cheb2
       double precision, allocatable, dimension(:,:), save :: &
            real2spec,spec2real,derive_spec,interpole
       double precision, allocatable, save :: derive(:,:,:)
       double precision, allocatable, save :: r_cheb(:)

contains

!---------------------------------------------------------------
       subroutine init_cheb(nr)

       implicit none
       integer nr 

       nr_cheb = nr
       if (allocated(real2spec))   deallocate(real2spec)
       if (allocated(spec2real))   deallocate(spec2real)
       if (allocated(derive_spec)) deallocate(derive_spec)
       if (allocated(derive))      deallocate(derive)
       if (allocated(r_cheb))      deallocate(r_cheb)

       call init_cheb_grid()
       call init_cheb_transform()
       end subroutine

!---------------------------------------------------------------
       subroutine init_cheb_grid()

       implicit none
       double precision pi
       integer i
       pi = dacos(-1d0)

       allocate(r_cheb(nr_cheb))

       do i=1,nr_cheb
         r_cheb(i) = 0.5d0*(1d0-dcos(pi*dble(i-1)/dble(nr_cheb-1)))
       enddo
       end subroutine
!---------------------------------------------------------------
      subroutine init_cheb_transform()

      implicit none
      double precision pi,aux
      integer i, j
      double precision, allocatable :: temp(:)

      pi = dacos(-1d0)
      allocate(real2spec(nr_cheb,nr_cheb),&
               spec2real(nr_cheb,nr_cheb),&
               temp(0:2*nr_cheb-2))

      do i=0,2*nr_cheb-2
        temp(i) = dcos(pi*dble(i)/dble(nr_cheb-1))
      enddo

      aux = 2d0/dble(nr_cheb-1)
      do i=1,nr_cheb
        do j=1,nr_cheb
          real2spec(i,j) = aux*temp(modulo((i-1)*(j-1),2*nr_cheb-2))
          spec2real(i,j) = temp(modulo((i-1)*(j-1),2*nr_cheb-2))
        enddo
      enddo
      deallocate(temp)

      do i=1,nr_cheb
        real2spec(1,i) = 0.5d0*real2spec(1,i)
        real2spec(nr_cheb,i) = 0.5d0*real2spec(nr_cheb,i)
        real2spec(i,1) = 0.5d0*real2spec(i,1)
        real2spec(i,nr_cheb) = 0.5d0*real2spec(i,nr_cheb)
      enddo
      end subroutine
!---------------------------------------------------------------
      subroutine init_cheb_derive(der_max)

      implicit none
      integer i,j,k,l,der_max
      double precision, allocatable :: temp(:,:)

      if (der_max.lt.1) &
         stop 'der_max_input must be >= 1 in init_cheb_derive'
      
      allocate(derive(nr_cheb,nr_cheb,der_max),&
               derive_spec(nr_cheb,nr_cheb),&
               temp(nr_cheb,nr_cheb))

      derive_spec(:,:) = 0d0
      derive(:,:,:) = 0d0
      temp(:,:) = 0d0

      do i=1,nr_cheb
        do j=i+1,nr_cheb,2
          derive_spec(i,j) = -dble(4*(j-1))
        enddo
      enddo

      do i=1,nr_cheb
        derive_spec(1,i) = 0.5d0*derive_spec(1,i)
      enddo

      call DGEMM('N','N',nr_cheb,nr_cheb,nr_cheb,1d0,derive_spec, &
                 nr_cheb,real2spec,nr_cheb,0d0,temp,nr_cheb)

      call DGEMM('N','N',nr_cheb,nr_cheb,nr_cheb,1d0,spec2real, &
                 nr_cheb,temp,nr_cheb,0d0,derive(:,:,1),nr_cheb)
      deallocate(temp)

      do l=2,der_max
        call DGEMM('N','N',nr_cheb,nr_cheb,nr_cheb,1d0,derive(:,:,1), &
                   nr_cheb,derive(:,:,l-1),nr_cheb, 0d0,derive(:,:,l), &
                   nr_cheb)
      enddo


      end subroutine

!---------------------------------------------------------------
       subroutine init_interpole(nr1,nr2)

       implicit none
       integer nr1,nr2,i,j
       double precision, allocatable :: temp1(:,:),temp2(:,:),aux(:)
       double precision pi

       if (allocated(interpole)) deallocate(interpole)
       nr_cheb1 = nr1
       nr_cheb2 = nr2

       pi = dacos(-1d0)

       allocate(aux(0:2*max(nr1,nr2)-2),temp1(nr1,nr1),&
                temp2(nr2,nr1),interpole(nr2,nr1))
       aux = 0d0
       temp1 = 0d0
       temp2 = 0d0
       interpole = 0d0

       do i=0,2*nr1-2
         aux(i) = 2d0*dcos(pi*dble(i)/dble(nr1-1))/dble(nr1-1)
       enddo

       do i=1,nr1
         do j=1,nr1
           temp1(i,j) = aux(modulo((i-1)*(j-1),2*nr1-2))
         enddo
       enddo
       do i=1,nr1
         temp1(1,i) = 0.5d0 * temp1(1,i) 
         temp1(nr1,i) = 0.5d0 * temp1(nr1,i) 
         temp1(i,1) = 0.5d0 * temp1(i,1) 
         temp1(i,nr1) = 0.5d0 * temp1(i,nr1) 
       enddo

       aux(:) = 0d0

       do i=0,2*nr2-2
         aux(i) = cos(pi*dble(i)/dble(nr2-1))
       enddo

       temp2(:,:) = 0d0
       do i=1,nr2
         do j=1,min(nr1,nr2)
           temp2(i,j) = aux(modulo((i-1)*(j-1),2*nr2-2))
         enddo
       enddo
      
       call DGEMM('N','N',nr2,nr1,nr1,1d0,temp2,nr2,temp1,nr1,0d0,&
                  interpole,nr2)
       deallocate(temp1,temp2,aux)
       end subroutine
!---------------------------------------------------------------
       subroutine project_function (x,f,fspec,n,nspec)

       implicit none
       integer n, nspec
       double precision x(n), f(n), fspec(nspec)

       double precision xmin, xmax
       double precision, allocatable :: theta(:), weights(:)
       integer i, j
       double precision, parameter :: pi = 3.14159265358979d0

       allocate(theta(1:n),weights(1:n))

       ! We're assuming the first point x(1) corresponds to theta = 0
       ! and the last point x(n) corresponds to theta = pi.

       xmin = x(1)
       xmax = x(n)
       do i=1,n
         theta(i) = dacos((xmax+xmin-2d0*x(i))/(xmax-xmin))
       enddo


       !!!! Warning: fspec(i) = a(i-1) where a(i) is the
       !!!! coefficient in the chebyshev expansion.

       do i=1,nspec
         fspec(i) = 0d0
         call find_weights(theta,weights,n,i-1)
         do j=1,n
           fspec(i) = fspec(i)+f(j)*weights(j)
         enddo
         fspec(i) = fspec(i)*2d0/pi
         if (i.eq.1) fspec(i) = fspec(i)/2d0
       enddo
       deallocate(theta,weights)
       end subroutine

!------------------------------------------------------------------------------ 
!  This subroutine sets up integration weights for an arbitrary grid.  It does 
!  this by integrating Lagrange interpolation polynomials calculated over a
!  sliding window which spans [i-order, i+order+1].  The function to be
!  integrated is then assumed to be multiplied by cos(nt grid(i)).
!------------------------------------------------------------------------------ 
! description of variables:
!
! grid(1:ngrid)    = grid on which the integration is carried out.
! weights(1:ngrid) = contains the integration weights.
! ngrid            = number of points in the grid
! nt               = the integration takes into account the weight function
!                    cos(nx)
! order            = positive integer which gives half the size of the window.
!------------------------------------------------------------------------------ 
      subroutine find_weights(grid,weights,ngrid,nt)

      implicit none
      integer ngrid, nt
      double precision grid(1:ngrid), weights(1:ngrid)
      integer i, j, k, l, start, finish
      double precision, allocatable :: mu(:), a(:)
      double precision prdct, my_sum
      integer, parameter :: order = 2
      
      ! check parameters:
      if (order.lt.0) stop "Please set order>= 0 in find_weights"

      allocate(mu(-order:(order+1)),a(-1:2*(order+1)))
      weights = 0d0

      do i=1,ngrid-1

        start  = -order
        finish = order+1

        ! special treatment do the endpoints
        ! Note: this leads to results which are less precise on the edges
        if ((i+start).lt.1) start = 1 - i
        if ((i+finish).gt.ngrid) finish = ngrid - i

        do j=start,finish
          mu(j) = grid(j+i) - grid(i)
        enddo

        do j=start,finish

          ! Initialise Lagrange polynomial
          a(-1) = 0d0
          a(0) = 1d0
          do l=1,2*(order+1)
            a(l) = 0d0
          enddo
          
          prdct = 1d0
          do k=start,j-1
            prdct = prdct*(mu(j)-mu(k))
          enddo
          do k=j+1,finish
            prdct = prdct*(mu(j)-mu(k))
          enddo
          a(0) = a(0)/prdct
          

          ! Calculate Lagrange polynomial, by calculating product of (x-mu(k))
          do k=start,j-1
            do l=2*(order+1),0,-1
              a(l) = -mu(k)*a(l) + a(l-1)
            enddo
          enddo
          do k=j+1,finish
            do l=2*(order+1),0,-1
              a(l) = -mu(k)*a(l) + a(l-1)
            enddo
          enddo

          ! Integrate the Lagrange polynomials, and add the results
          ! as weights at the appropriate points.

          if (nt.eq.0) then

            do l=0,2*(order+1)
              weights(i+j) = weights(i+j) + a(l)*mu(1)**(l+1)/dble(l+1)
            enddo

          else

            my_sum = (dsin(dble(nt)*grid(i+1))-dsin(dble(nt)*grid(i)))/dble(nt)
            do l=0,2*(order+1),2
              weights(i+j) = weights(i+j) + a(l)*my_sum
              my_sum =   mu(1)**(l+2)*dsin(dble(nt)*grid(i+1))/dble(nt)           &
                    + dble(l+2)*mu(1)**(l+1)*dcos(dble(nt)*grid(i+1))/dble(nt*nt) &
                    - dble((l+2)*(l+1))*my_sum/dble(nt*nt)
            enddo
            my_sum =  mu(1)*dsin(dble(nt)*grid(i+1))/dble(nt)                     &
                 + (dcos(dble(nt)*grid(i+1))-dcos(dble(nt)*grid(i)))/dble(nt*nt)
            do l=1,2*(order+1),2
              weights(i+j) = weights(i+j) + a(l)*my_sum
              my_sum =   mu(1)**(l+2)*dsin(dble(nt)*grid(i+1))/dble(nt)           &
                    + dble(l+2)*mu(1)**(l+1)*dcos(dble(nt)*grid(i+1))/dble(nt*nt) &
                    - dble((l+2)*(l+1))*my_sum/dble(nt*nt)
            enddo


          endif

        enddo
      enddo

      deallocate(a,mu)

      end subroutine

!------------------------------------------------------------------------------ 
!  This subroutine interpolate a function onto the Chebyshev grid of size
!  nr_cheb, using Lagrange interpolation, on the npoints nearest  points.
!------------------------------------------------------------------------------ 
! description of variables:
!
! r(1:n)       = original grid
! f(1:n)       = function to be interpolated
! n            = number of grid points in the original grid
! fout(1:nout) = function after interpolation
! nout         = number of grid points after interpolation
! npoints      = number of points used to calculate the Lagrange interpolation
!                polynomial
!------------------------------------------------------------------------------ 
      subroutine interpole_cheb(r,f,n,fout,nout)

      implicit none
      integer n, nout
      double precision r(1:n), f(1:n), fout(1:nout)
      double precision, parameter :: pi = 3.14159265358979d0
      double precision x, term
      integer start, finish, index, i, j, k
      integer, parameter :: npoints = 6

      ! check parameters:
      if (npoints.lt.2) stop "Please set npoints >= 2 in interpolate_cheb"

      fout = 0d0

      index = 1
      do i=1,nout

        x = r(1) + (r(n)-r(1))*(1d0-dcos(dble(i-1)*pi/dble(nout-1)))/2d0

        ! Find nearest point
        if (i.eq.nout) x = r(n) ! Just in case, to avoid numerical errors
                                ! messing up the next test.

        do while (r(index+1).lt.x)
          index = index + 1
        enddo

        if ((x-r(index)).lt.(r(index+1)-x)) then
          start = index
        else
          start = index+1
        endif
        finish = start

        ! Find an interval of npoint nearest points
        do j=2,npoints

          if (start.eq.1) then
            finish = npoints
            exit
          endif

          if (finish.eq.n) then
            start = n-npoints+1
            exit
          endif

          if ((x-r(start-1)).lt.(r(finish+1)-x)) then
            start = start - 1
          else
            finish = finish + 1
          endif
        enddo
  
        ! Lagrange interpolation
        do j=start,finish
          term = f(j)
          do k=start,j-1
            term = term*(x-r(k))/(r(j)-r(k))
          enddo
          do k=j+1,finish
            term = term*(x-r(k))/(r(j)-r(k))
          enddo
          fout(i) = fout(i) + term
        enddo
      enddo

      end subroutine
!------------------------------------------------------------------------------ 
      end module
