!------------------------------------------------------------------------------ 
!  This module sets up a derivation matrix "derive" in general form.
!  This matrix contains derivatives of order der_min to der_max.
!
!  Different representations of the dervation operator are implemented:
!    -> finite differences: applies to an (almost) arbitrary grid.  Different
!                  versions are implemented.
!    -> Chebyshev: assumes the function is expressed on the Gauss-Lobatto
!                  Chebyshev collocation grid.
!    -> b-splines: applies to an (almost) arbitrary grid.
!------------------------------------------------------------------------------ 

      module derivative

      type DERMAT
        integer :: der_max, der_min, der_id, ngrid
        integer, allocatable :: ubder(:), lbder(:)
        double precision, allocatable :: derive(:,:,:)
      end type DERMAT

!------------------------------------------------------------------------------ 
! Variables:
!
! der_min= mininal derivative order.  Is usually 0.  Sometimes a negative
!          value is used, as storage space for other matrices.
! der_max= highest order derivative
! der_id = derivative order, or matrix index which can be used to obtain
!          the original function (usually is equal to 0 but can be -1)
! ngrid  = number of grid points
! lbder  = number of lower bands
! ubder  = number of upper bands
!          Note: the matrix is not stored in banded form even if it may
!                banded
! derive = contains derivation matrices in general form.
!          derive(1:ngrid,1:ngrid,1:der_max) (where ngrid = number of grid points
!          and der_max = highest order derivative)
!------------------------------------------------------------------------------ 

contains
!------------------------------------------------------------------------------ 
!  This subroutine calls the correct subroutine to initialise the derivation
!  matrix, based on the value of dertype.
!------------------------------------------------------------------------------ 
      subroutine init_derive(dmat,grid,ngrid,der_max,der_min,order,&
                             dertype)

      implicit none
      type(DERMAT), intent(out)    :: dmat 
      integer, intent(in)          :: ngrid, der_max, der_min, order
      double precision, intent(in) :: grid(ngrid)
      character*(4), intent(in)    :: dertype

      if (trim(dertype).eq.'FD') then
        print*,"Derivative type: finite differences"
        call init_derive_FD(dmat,grid,ngrid,der_max,der_min,order)
      elseif (trim(dertype).eq.'FD2') then
        print*,"Derivative type: finite differences (type 2)"
        call init_derive_FD2(dmat,grid,ngrid,der_max,der_min,order)
      elseif (trim(dertype).eq.'CHEB') then
        print*,"Derivative type: spectral, using chebyshev polynomials"
        call init_derive_cheb(dmat,grid,ngrid,der_max,der_min)
      elseif (trim(dertype).eq.'BSPL') then
        print*,"Derivative type: based on B-splines"
        call init_derive_b_splines(dmat,grid,ngrid,der_max,der_min,order)
      elseif (trim(dertype).eq.'SMPL') then
        print*,"Derivative type: 1st order centered FD derivatives"
        call init_derive_FD_simple(dmat,grid,ngrid,der_max,der_min)
      elseif (trim(dertype).eq.'IFD') then
        print*,"Derivative type: improved finite differences"
        call init_derive_IFD(dmat,grid,ngrid,der_max,der_min,order)
      elseif (trim(dertype).eq.'MFD') then
        print*,"Derivative type: mid-point finite differences"
        call init_derive_MFD(dmat,grid,ngrid,der_max,der_min,order)
      elseif (trim(dertype).eq.'STAG') then
        print*,"Derivative type: staggered grid finite "// &
               "differences (experimental)"
        call init_derive_STAG(dmat,grid,ngrid,der_max,der_min)
      elseif (trim(dertype).eq.'MIX') then
        print*,"Derivative type: mixed finite differences "// &
               "(experimental)"
        call init_derive_MIX(dmat,grid,ngrid,der_max,der_min,order)
      else
        print*,"dertype = ",dertype
        stop 'dertype has a faulty value in init_derive'
      endif

      return
      end subroutine init_derive
!------------------------------------------------------------------------------ 
!  This subroutine clears a derivation matrix (so as to avoid memory leaks).
!------------------------------------------------------------------------------ 
      subroutine clear_derive(dmat)

      implicit none
      type(DERMAT), intent(out)    :: dmat 

      if (allocated(dmat%derive)) deallocate(dmat%derive)
      if (allocated(dmat%ubder))  deallocate(dmat%ubder)
      if (allocated(dmat%lbder))  deallocate(dmat%lbder)

      end subroutine clear_derive
!------------------------------------------------------------------------------ 
!  This subroutine sets up a banded derivation matrices (which is stored in
!  general form) using nth order finite difference schemes.  The
!  derivatives are calculated with respect to an (almost) arbitrary grid.
!  
!  It calculates the derivation coefficients based on Lagrange interpolation
!  polynomials (cf. Fornberg, 1988).
!------------------------------------------------------------------------------ 
! Variables:
!
! grid        = underlying on which the derivatives are calculated (grid(1:ngrid))
!            Warning: this grid must not contain two identical points.
! ngrid       = number of grid points (1:ngrid)
! order    = order of the scheme, divided by 2.  Note:  (order >= der_max)
!------------------------------------------------------------------------------ 
      subroutine init_derive_FD(dmat,grid,ngrid,der_max,der_min,order)

      implicit none
      type(DERMAT), intent(out)    :: dmat
      integer, intent(in)          :: ngrid, der_max, der_min, order
      double precision, intent(in) :: grid(ngrid)
      integer i, start, finish, j, k, l
      double precision prdct
      double precision, allocatable :: mu(:), a(:)

      ! check parameters:
      if (2*order.lt.der_max) then
       stop "Please set 2*order >= der_max in init_derive_FD"
      endif
      if (order.lt.der_max) then
       print*,"Warning: order < der_max in init_derive_FD"
      endif
      if (order.le.0) then
       stop "Please set order > 0 in init_derive_FD"
      endif
      if (der_max.lt.1) then
       stop "Please set der_max > 0 in init_derive_FD"
      endif
      if (der_min.ne.0) then
       stop "Please set der_min = 0 in init_derive_FD"
      endif

      ! initialise der_min, der_max and ngrid in dmat:
      dmat%der_min = der_min
      dmat%der_max = der_max
      dmat%der_id  = 0
      dmat%ngrid   = ngrid

      allocate(dmat%lbder(0:der_max),dmat%ubder(0:der_max))
      dmat%ubder(0) = 0
      dmat%lbder(0) = 0
      dmat%ubder(1:der_max) = order
      dmat%lbder(1:der_max) = order

      allocate(dmat%derive(1:ngrid,1:ngrid,0:der_max))
      dmat%derive(:,:,:) = 0d0

      allocate(a(-1:der_max),mu(-order:order))

      do i=1,ngrid
        dmat%derive(i,i,0) = 1d0
      enddo

      do i=1,ngrid

        start  = -order
        finish = order

        ! special treatment for the endpoints
        ! Note: this leads to results which are less precise on the
        ! edges, but preserves the banded form of the derivation matrix.
        if ((i+start).lt.1) start = 1 - i
        if ((i+finish).gt.ngrid) finish = ngrid - i

        do j=start,finish
          mu(j) = grid(j+i) - grid(i)
        enddo

        do j=start,finish

          ! Initialise Lagrange polynomial
          a(-1) = 0d0
          a(0) = 1d0
          do l=1,der_max
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
            do l=der_max,0,-1
              a(l) = -mu(k)*a(l) + a(l-1)
            enddo
          enddo
          do k=j+1,finish
            do l=der_max,0,-1
              a(l) = -mu(k)*a(l) + a(l-1)
            enddo
          enddo

          ! The coeffecients a(l) of the Lagrange polynomial are the
          ! derivation coefficients divided by l!.
          prdct = 1d0
          do l=1,der_max
            prdct = prdct*dble(l)
            dmat%derive(i,i+j,l) = a(l)*prdct
          enddo
        enddo
      enddo

      deallocate(a,mu)

      return
      end subroutine init_derive_FD

!------------------------------------------------------------------------------ 
!  This subroutine sets up a banded derivation matrices (which is stored in
!  general form) using nth order finite difference schemes.  The
!  derivatives are calculated with respect to an (almost) arbitrary grid.
!  
!  It calculates the first order derivation coefficients based on Lagrange
!  interpolation polynomials (cf. Fornberg, 1988).  Higher order derivatives
!  are obtained by multiplying calculating powers of the first order
!  derivation matrix.
!------------------------------------------------------------------------------ 
! Variables:
!
! grid        = underlying on which the derivatives are calculated (grid(1:ngrid))
!            Warning: this grid must not contain two identical points.
! ngrid       = number of grid points (1:ngrid)
! order    = order of the scheme, divided by 2.  Note:  (order >= der_max)
!------------------------------------------------------------------------------ 
      subroutine init_derive_FD2(dmat,grid,ngrid,der_max,der_min,order)

      implicit none
      type(DERMAT), intent(out)    :: dmat
      integer, intent(in)          :: ngrid, der_max, der_min, order
      double precision, intent(in) :: grid(ngrid)
      integer i, start, finish, j, k, l
      double precision prdct
      double precision, allocatable :: mu(:), a(:)

      ! check parameters:
      if (order.le.0) then
       stop "Please set order > 0 in init_derive_FD2"
      endif
      if (der_max.lt.1) then
       stop "Please set der_max > 0 in init_derive_FD2"
      endif
      if (der_min.ne.0) then
       stop "Please set der_min = 0 in init_derive_FD2"
      endif

      ! initialise der_min, der_max and ngrid in dmat:
      dmat%der_min = der_min
      dmat%der_max = der_max
      dmat%der_id  = 0
      dmat%ngrid   = ngrid

      allocate(dmat%lbder(0:der_max),dmat%ubder(0:der_max))
      dmat%ubder(0) = 0
      dmat%lbder(0) = 0
      dmat%ubder(1) = order
      dmat%lbder(1) = order

      allocate(dmat%derive(1:ngrid,1:ngrid,0:der_max))
      dmat%derive(:,:,:) = 0d0

      ! Calculate 1st order derivative.  Same as in init_derive_FD, except
      ! that we stop at 1 rather than der_max

      allocate(a(0:1),mu(-order:order))

      do i=1,ngrid
        dmat%derive(i,i,0) = 1d0
      enddo

      do i=1,ngrid

        start  = -order
        finish = order

        ! special treatment for the endpoints
        ! Note: this leads to results which are less precise on the
        ! edges, but preserves the banded form of the derivation matrix.
        if ((i+start).lt.1) start = 1 - i
        if ((i+finish).gt.ngrid) finish = ngrid - i

        do j=start,finish
          mu(j) = grid(j+i) - grid(i)
        enddo

        do j=start,finish

          ! Initialise Lagrange polynomial
          a(0) = 1d0
          a(1) = 0d0
          
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
            a(1) = -mu(k)*a(1) + a(0)
            a(0) = -mu(k)*a(0)
          enddo
          do k=j+1,finish
            a(1) = -mu(k)*a(1) + a(0)
            a(0) = -mu(k)*a(0)
          enddo

          dmat%derive(i,i+j,1) = a(1)
        enddo
      enddo

      do l=2,der_max
        call dgemm('N','N',ngrid,ngrid,ngrid,1d0,              &
                   dmat%derive(1:ngrid,1:ngrid,l-1),ngrid,     &
                   dmat%derive(1:ngrid,1:ngrid,1),ngrid,0d0,   &
                   dmat%derive(1:ngrid,1:ngrid,l),ngrid)

        ! When multiplying two band matrices, the result is a band matrix
        ! in which the number of upper (lower) bands is equal to the sum
        ! of the number upper (lower, resp.) bands of the two previous
        ! matrices  (where the diagonal is not counted as a band).
        dmat%lbder(l) = min(dmat%lbder(1)+dmat%lbder(l-1),ngrid)
        dmat%ubder(l) = min(dmat%ubder(1)+dmat%ubder(l-1),ngrid)
      enddo
      deallocate(a,mu)

      return
      end subroutine init_derive_FD2

!------------------------------------------------------------------------------ 
!  This module sets up a derivation matrices for functions expressed on a
!  Gauss-Lobatto Chebyshev collocation grid (cf. Canuto et al., 1988).
!------------------------------------------------------------------------------ 
! Variables:
!
! ngrid       = number of grid points (1:ngrid)
! der_max  = maximal order of derivative
!------------------------------------------------------------------------------ 
      subroutine init_derive_cheb(dmat,grid,ngrid,der_max,der_min)

      implicit none
      type(DERMAT), intent(out)    :: dmat
      integer, intent(in)          :: ngrid, der_max, der_min
      double precision, intent(in) :: grid(ngrid)
      double precision rspan
      integer i,j,k,l
      double precision, allocatable :: temp(:)
      double precision, dimension(:,:), allocatable :: &
             aux, real2spec, spec2real, derive_spec
      double precision, parameter :: pi = 3.14159265358979d0
      double precision cnst

      rspan = grid(ngrid)-grid(1)

      if (der_max.lt.1) then
       stop 'der_max_input must be >= 1 in init_derive_cheb'
      endif
      if (der_min.ne.0) then
       stop "Please set der_min = 0 in init_derive_cheb"
      endif

      ! initialise der_min, der_max and ngrid in dmat:
      dmat%der_min = der_min
      dmat%der_max = der_max
      dmat%der_id  = 0
      dmat%ngrid   = ngrid

      allocate(dmat%lbder(0:der_max),dmat%ubder(0:der_max))
      dmat%ubder(0) = 0
      dmat%lbder(0) = 0
      dmat%ubder(1:der_max) = ngrid-1
      dmat%lbder(1:der_max) = ngrid-1

      allocate(dmat%derive(1:ngrid,1:ngrid,0:der_max))
      dmat%derive(:,:,:) = 0d0

      allocate(derive_spec(ngrid,ngrid),temp(0:2*ngrid-2),&
               aux(ngrid,ngrid),real2spec(ngrid,ngrid),&
               spec2real(ngrid,ngrid))

      derive_spec(:,:) = 0d0
      real2spec(:,:) = 0d0
      spec2real(:,:) = 0d0
      aux(:,:) = 0d0
      temp(:) = 0d0

      do i=1,ngrid
        dmat%derive(i,i,0) = 1d0
      enddo

      ! Create matrices which transform a vector from real to
      ! spectral representation and vice-versa.
      do i=0,2*ngrid-2
        temp(i) = dcos(pi*dble(i)/dble(ngrid-1))
      enddo

      cnst = 2d0/dble(ngrid-1)
      do i=1,ngrid
        do j=1,ngrid
          real2spec(i,j) = cnst*temp(modulo((i-1)*(j-1),2*ngrid-2))
          spec2real(i,j) = temp(modulo((i-1)*(j-1),2*ngrid-2))
        enddo
      enddo
      deallocate(temp)

      do i=1,ngrid
        real2spec(1,i) = 0.5d0*real2spec(1,i)
        real2spec(ngrid,i) = 0.5d0*real2spec(ngrid,i)
        real2spec(i,1) = 0.5d0*real2spec(i,1)
        real2spec(i,ngrid) = 0.5d0*real2spec(i,ngrid)
      enddo

      ! Create a spectral representation of the derivation matrix
      do i=1,ngrid
        do j=i+1,ngrid,2
          derive_spec(i,j) = -dble(4*(j-1))
        enddo
      enddo

      do i=1,ngrid
        derive_spec(1,i) = 0.5d0*derive_spec(1,i)
      enddo

      ! Transform the derivation matrix from spectral to real
      ! representation
      call DGEMM('N','N',ngrid,ngrid,ngrid,1d0,derive_spec, &
                 ngrid,real2spec,ngrid,0d0,aux,ngrid)

      call DGEMM('N','N',ngrid,ngrid,ngrid,1d0,spec2real, &
                 ngrid,aux,ngrid,0d0,dmat%derive(:,:,1),ngrid)
      deallocate(aux, real2spec, spec2real, derive_spec)

      ! This is for grids which span a different interval than [0,1].
      dmat%derive(:,:,1) = dmat%derive(:,:,1)/rspan

      ! Multiply the derivation matrix by itself in order to obtain
      ! higher order derivatives
      do l=2,der_max
        call DGEMM('N','N',ngrid,ngrid,ngrid,1d0,dmat%derive(:,:,1), &
                   ngrid,dmat%derive(:,:,l-1),ngrid, 0d0,            &
                   dmat%derive(:,:,l),ngrid)
      enddo

      return
      end subroutine init_derive_cheb

!------------------------------------------------------------------------------ 
!  This sets up a simple 2nd order finite difference 1st derivative matrix.
!  The der = 0 matrix corresponds to simple interpolation on mid-points and
!  the der = -1 matrix to the identity matrix (which is useful for boundary
!  conditions and equations without any derivatives).
!------------------------------------------------------------------------------ 
      subroutine init_derive_FD_simple(dmat,grid,ngrid,der_max,der_min)

      implicit none
      type(DERMAT), intent(out)    :: dmat
      integer, intent(in)          :: ngrid, der_max, der_min
      double precision, intent(in) :: grid(ngrid)
      integer i

      ! check parameters:
      if (der_max.ne.1) then
       stop "Please set der_max = 1 in init_derive_FD_simple"
      endif
      if (der_min.ne.-1) then
       stop "Please set der_min = -1 in init_derive_FD_simple"
      endif

      ! initialise der_min, der_max and ngrid in dmat:
      dmat%der_min = der_min
      dmat%der_max = der_max
      dmat%der_id  = -1
      dmat%ngrid   = ngrid

      allocate(dmat%lbder(der_min:der_max),dmat%ubder(der_min:der_max))
      dmat%ubder(-1) = 0
      dmat%lbder(-1) = 0
      dmat%ubder(0) = 1
      dmat%lbder(0) = 0
      dmat%ubder(1) = 1
      dmat%lbder(1) = 0

      allocate(dmat%derive(1:ngrid,1:ngrid,der_min:der_max))
      dmat%derive(:,:,:) = 0d0

      do i=1,ngrid-1
        dmat%derive(i,i+1,0) = 0.5d0
        dmat%derive(i,i,0)   = 0.5d0
        dmat%derive(i,i+1,1) = 1d0/(grid(i+1)-grid(i))
        dmat%derive(i,i,1)   = -dmat%derive(i,i+1,1)
        dmat%derive(i,i,-1)  = 1d0
      enddo
      dmat%derive(ngrid,ngrid,-1)  = 1d0

      return
      end subroutine init_derive_FD_simple

!------------------------------------------------------------------------------ 
!  This subroutine sets up a banded derivation matrices (which is stored in
!  general form) using nth order finite difference schemes.  The
!  derivatives are calculated with respect to an (almost) arbitrary grid.
!  The order of the scheme is 2*order, even though the number of points
!  in the window is also 2*order.  This can be achieved for 1st order
!  differential systems by carefully defining a new grid which is
!  approximately (or exactly, when order = 1) located at the grid mid-
!  points.  This method hopefully minimises the effects of mesh drift.
!  
!  The derivation coefficients are based on Lagrange interpolation
!  polynomials (cf. Fornberg, 1988).
!------------------------------------------------------------------------------ 
! Variables:
!
! grid     = underlying on which the derivatives are calculated (grid(1:ngrid))
!            Warning: this grid must not contain two identical points.
! ngrid    = number of grid points (1:ngrid)
! order    = order of the scheme, divided by 2.
!------------------------------------------------------------------------------ 
! NOTE: the last line in the derive matrices is full of zeros (except
!       for derive(:,:,-1)).  This typically where one is expected to
!       place boundary conditions.  If one needs a boundary condition at
!       the first point of the domain, while keeping a banded structure
!       for the matrix, this can be achieved by setting up init_order
!       carefully.
!------------------------------------------------------------------------------ 
      subroutine init_derive_IFD(dmat,grid,ngrid,der_max,der_min,order)

      implicit none
      type(DERMAT), intent(out)    :: dmat
      integer, intent(in)          :: ngrid, der_max, der_min, order
      double precision, intent(in) :: grid(ngrid)
      integer i, j, k, l, start, finish, taille
      double precision prdct, num, den, xx, dx, acoef, bcoef
      double precision, allocatable :: mu(:), a(:), aa(:)
      double precision, allocatable :: new_grid(:)
      double precision, parameter :: eps = 1d-14
      integer, parameter :: n_iter_max = 100

      ! check parameters:
      if (order.lt.1) then
       stop "Please set order >= 1 in init_derive_IFD"
      endif
      if (der_max.lt.0) then
       stop "Please set der_max >= 0 in init_derive_IFD"
      endif
      if (der_min.ne.-1) then
       stop "Please set der_min = -1 in init_derive_IFD"
      endif

      ! initialise der_min, der_max and ngrid in dmat:
      dmat%der_min = der_min
      dmat%der_max = der_max
      dmat%der_id  = -1
      dmat%ngrid   = ngrid

      allocate(dmat%lbder(der_min:der_max),dmat%ubder(der_min:der_max))

      dmat%ubder(-1) = 0
      dmat%lbder(-1) = 0
      dmat%ubder(0) = order
      dmat%lbder(0) = order-1
      dmat%ubder(1) = order
      dmat%lbder(1) = order-1

      allocate(dmat%derive(1:ngrid,1:ngrid,der_min:der_max))
      dmat%derive(:,:,:) = 0d0

      ! set up identity matrix
      do i=1,ngrid
        dmat%derive(i,i,-1) = 1d0
      enddo

      ! define new grid
      allocate(new_grid(ngrid-1),aa(-1:2*order))
      do i=1,ngrid-1

        start  = -order+1
        finish = order

        ! special treatment for the endpoints
        ! Note: this leads to results which are less precise on the
        ! edges, but preserves the banded form of the derivation matrix.
        if ((i+start).lt.1) start = 1 - i
        if ((i+finish).gt.ngrid) finish = ngrid - i

        ! scale polynomial
        acoef = 2d0/(grid(i+1)-grid(i))
        bcoef = (grid(i)+grid(i+1))/(grid(i)-grid(i+1))

        ! find polynomial coefficients:
        aa(:) = 0d0
        aa(0) = 1d0
        do j=start,finish
          do k=2*order,0,-1
            aa(k) = aa(k-1)-(acoef*grid(i+j)+bcoef)*aa(k)
          enddo
        enddo

        ! find a root of the polynomial: this will become the new grid
        ! point
        xx = 0d0
        dx = 1d0

        ! Newton type iteration
        j = 0
        do while ((abs(dx).gt.eps).and.(j.lt.n_iter_max))
          num = aa(2*order)*dble(2*order)
          do k=2*order-1,1,-1 
            num = xx*num+dble(k)*aa(k)
          enddo
          den = dble(2*order*(2*order-1))*aa(2*order)
          do k=2*order-1,2,-1 
            den = xx*den+dble(k*(k-1))*aa(k)
          enddo
          dx = num/den
          xx = xx - dx
          j  = j + 1
        enddo
        xx = (xx-bcoef)/acoef

        ! sanity check (just in case)
        if ((xx.le.grid(i)).or.(xx.ge.grid(i+1)))  then
          write(*,*) (grid(i+j),j=start,finish)
          stop "Error in init_derive_IFD: stray grid point"
        endif

        ! assign new grid point
        new_grid(i) = xx
      enddo
      deallocate(aa)

      allocate(a(-1:der_max),mu(-order:order))

      do i=1,ngrid-1

        start  = -order+1
        finish = order

        ! special treatment for the endpoints
        ! Note: this leads to results which are less precise on the
        ! edges, but preserves the banded form of the derivation matrix.
        if ((i+start).lt.1) start = 1 - i
        if ((i+finish).gt.ngrid) finish = ngrid - i

        do j=start,finish
          mu(j) = grid(j+i) - new_grid(i)
        enddo

        do j=start,finish

          ! Initialise Lagrange polynomial
          a(-1) = 0d0
          a(0) = 1d0
          do l=1,der_max
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
            do l=der_max,0,-1
              a(l) = -mu(k)*a(l) + a(l-1)
            enddo
          enddo
          do k=j+1,finish
            do l=der_max,0,-1
              a(l) = -mu(k)*a(l) + a(l-1)
            enddo
          enddo

          ! The coeffecients a(l) of the Lagrange polynomial are the
          ! derivation coefficients divided by l! (where 0! = 1)
          prdct = 1d0
          do l=0,der_max
            prdct = prdct*dble(max(1,l))
            dmat%derive(i,i+j,l) = a(l)*prdct
          enddo
        enddo
      enddo

      deallocate(a,mu)

      return
      end subroutine init_derive_IFD

!------------------------------------------------------------------------------ 
!  This subroutine sets up a banded derivation matrices (which is stored in
!  general form) using nth order finite difference schemes.  The
!  derivatives are calculated with respect to an (almost) arbitrary grid.
!  The order of the scheme is 2*order-1 (in most cases) and uses windows
!  of 2*order grid points.  The differential systems is enforced at
!  the grid mid-points.  This reduces the effects of mesh drift, but
!  does not take advantage of boosted accuracy as in the
!  init_derive_IFD.
!  
!  The derivation coefficients are based on Lagrange interpolation
!  polynomials (cf. Fornberg, 1988).
!------------------------------------------------------------------------------ 
! Variables:
!
! grid     = underlying on which the derivatives are calculated (grid(1:ngrid))
!            Warning: this grid must not contain two identical points.
! ngrid    = number of grid points (1:ngrid)
! order    = order of the scheme, divided by 2.
!------------------------------------------------------------------------------ 
! NOTE: the last line in the derive matrices is full of zeros (except
!       for derive(:,:,-1)).  This typically where one is expected to
!       place boundary conditions.  If one needs a boundary condition at
!       the first point of the domain, while keeping a banded structure
!       for the matrix, this can be achieved by setting up init_order
!       carefully.
!------------------------------------------------------------------------------ 
      subroutine init_derive_MFD(dmat,grid,ngrid,der_max,der_min,order)

      implicit none
      type(DERMAT), intent(out)    :: dmat
      integer, intent(in)          :: ngrid, der_max, der_min, order
      double precision, intent(in) :: grid(ngrid)
      integer i, j, k, l, start, finish, taille
      double precision prdct, num, den, xx, dx, acoef, bcoef
      double precision, allocatable :: mu(:), a(:)
      double precision, allocatable :: new_grid(:)
      double precision, parameter :: eps = 1d-14
      integer, parameter :: n_iter_max = 100

      ! check parameters:
      if (order.lt.1) then
       stop "Please set order >= 1 in init_derive_MFD"
      endif
      if (der_max.lt.0) then
       stop "Please set der_max >= 0 in init_derive_MFD"
      endif
      if (der_min.ne.-1) then
       stop "Please set der_min = -1 in init_derive_MFD"
      endif

      ! initialise der_min, der_max and ngrid in dmat:
      dmat%der_min = der_min
      dmat%der_max = der_max
      dmat%der_id  = -1
      dmat%ngrid   = ngrid

      allocate(dmat%lbder(der_min:der_max),dmat%ubder(der_min:der_max))
      dmat%ubder(-1) = 0
      dmat%lbder(-1) = 0
      dmat%ubder(0) = order
      dmat%lbder(0) = order-1
      dmat%ubder(1) = order
      dmat%lbder(1) = order-1

      allocate(dmat%derive(1:ngrid,1:ngrid,der_min:der_max))
      dmat%derive(:,:,:) = 0d0

      ! set up identity matrix
      do i=1,ngrid
        dmat%derive(i,i,-1) = 1d0
      enddo

      ! define new grid
      allocate(new_grid(ngrid-1))
      do i=1,ngrid-1
        new_grid(i) = (grid(i) + grid(i+1))/2d0
      enddo

      allocate(a(-1:der_max),mu(-order:order))

      do i=1,ngrid-1

        start  = -order+1
        finish = order

        ! special treatment for the endpoints
        ! Note: this leads to results which are less precise on the
        ! edges, but preserves the banded form of the derivation matrix.
        if ((i+start).lt.1) start = 1 - i
        if ((i+finish).gt.ngrid) finish = ngrid - i

        do j=start,finish
          mu(j) = grid(j+i) - new_grid(i)
        enddo

        do j=start,finish

          ! Initialise Lagrange polynomial
          a(-1) = 0d0
          a(0) = 1d0
          do l=1,der_max
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
            do l=der_max,0,-1
              a(l) = -mu(k)*a(l) + a(l-1)
            enddo
          enddo
          do k=j+1,finish
            do l=der_max,0,-1
              a(l) = -mu(k)*a(l) + a(l-1)
            enddo
          enddo

          ! The coeffecients a(l) of the Lagrange polynomial are the
          ! derivation coefficients divided by l! (where 0! = 1)
          prdct = 1d0
          do l=0,der_max
            prdct = prdct*dble(max(1,l))
            dmat%derive(i,i+j,l) = a(l)*prdct
          enddo
        enddo
      enddo

      deallocate(a,mu)

      return
      end subroutine init_derive_MFD

!------------------------------------------------------------------------------ 
!  This sets up a simple 2nd order finite differences using the
!  staggered grid approach.  This is only for experimental purposes and
!  not meant for producing accurate results.
!------------------------------------------------------------------------------ 
      subroutine init_derive_STAG(dmat,grid,ngrid,der_max,der_min)

      implicit none
      type(DERMAT), intent(out)    :: dmat
      integer, intent(in)          :: ngrid, der_max, der_min
      double precision, intent(in) :: grid(ngrid)
      integer i

      ! check parameters:
      if (der_max.ne.1) then
       stop "Please set der_max = 1 in init_derive_STAG"
      endif
      if (der_min.ne.-1) then
       stop "Please set der_min = -1 in init_derive_STAG"
      endif

      ! initialise der_min, der_max and ngrid in dmat:
      dmat%der_min = der_min
      dmat%der_max = der_max
      dmat%der_id  = 0
      dmat%ngrid   = ngrid

      allocate(dmat%lbder(der_min:der_max),dmat%ubder(der_min:der_max))
      dmat%ubder(-1) = 1
      dmat%lbder(-1) = 0
      dmat%ubder(0) = 0
      dmat%lbder(0) = 0
      dmat%ubder(1) = 1
      dmat%lbder(1) = 0

      allocate(dmat%derive(1:ngrid,1:ngrid,der_min:der_max))
      dmat%derive(:,:,:) = 0d0

      do i=1,ngrid-1
        dmat%derive(i,i,0)   = 1d0
        dmat%derive(i,i+1,1) = 1d0/(grid(i+1)-grid(i))
        dmat%derive(i,i,1)   = -dmat%derive(i,i+1,1)
        dmat%derive(i,i+1,-1)= 1d0
      enddo
      dmat%derive(ngrid,ngrid,0)  = 1d0

      return
      end subroutine init_derive_STAG

!------------------------------------------------------------------------------ 
!  This subroutine sets up a banded derivation matrices (which is stored in
!  general form) using nth order finite difference schemes.  The
!  zeroth and first derivatives are calculated with respect to the
!  original grid, whereas the second derivative is calculated with
!  respect to the midpoints
!  
!  It calculates the derivation coefficients based on Lagrange interpolation
!  polynomials (cf. Fornberg, 1988).
!------------------------------------------------------------------------------ 
! Variables:
!
! grid        = underlying on which the derivatives are calculated (grid(1:ngrid))
!            Warning: this grid must not contain two identical points.
! ngrid       = number of grid points (1:ngrid)
! order    = order of the scheme, divided by 2.  Note:  (order >= der_max)
!------------------------------------------------------------------------------ 
      subroutine init_derive_MIX(dmat,grid,ngrid,der_max,der_min,order)

      implicit none
      type(DERMAT), intent(out)    :: dmat
      integer, intent(in)          :: ngrid, der_max, der_min, order
      double precision, intent(in) :: grid(ngrid)
      integer i, start, finish, j, k, l
      double precision prdct
      double precision, allocatable :: mu(:), a(:)

      ! check parameters:
      if (2*order.lt.der_max) then
       print*,"Warning: 2*order < der_max in init_derive_MIX"
      endif
      if (order.lt.der_max) then
       print*,"Warning: order < der_max in init_derive_MIX"
      endif
      if (order.le.0) then
       stop "Please set order > 0 in init_derive_MIX"
      endif
      if (der_max.gt.3) then
       stop "Please set der_max <= 3 in init_derive_MIX"
      endif
      if (der_max.lt.2) then
       stop "Please set der_max >= 2 in init_derive_MIX"
      endif
      if (der_min.ne.0) then
       stop "Please set der_min = 0 in init_derive_MIX"
      endif

      ! initialise der_min, der_max and ngrid in dmat:
      dmat%der_min = der_min
      dmat%der_max = der_max
      dmat%der_id  = 0
      dmat%ngrid   = ngrid

      allocate(dmat%lbder(0:der_max),dmat%ubder(0:der_max))
      dmat%ubder(0) = 0
      dmat%lbder(0) = 0
      dmat%ubder(1:der_max) = order
      dmat%lbder(1:der_max) = order

      allocate(dmat%derive(1:ngrid,1:ngrid,0:der_max))
      dmat%derive(:,:,:) = 0d0

      allocate(a(-1:der_max),mu(-order:order+1))

      do i=1,ngrid
        dmat%derive(i,i,0) = 1d0
      enddo

      do i=1,ngrid

        start  = -order
        finish = order

        ! special treatment for the endpoints
        ! Note: this leads to results which are less precise on the
        ! edges, but preserves the banded form of the derivation matrix.
        if ((i+start).lt.1) start = 1 - i
        if ((i+finish).gt.ngrid) finish = ngrid - i

        do j=start,finish
          mu(j) = grid(j+i) - grid(i)
        enddo

        do j=start,finish

          ! Initialise Lagrange polynomial
          a(-1) = 0d0
          a(0) = 1d0
          a(1) = 0d0
          
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
            do l=1,0,-1
              a(l) = -mu(k)*a(l) + a(l-1)
            enddo
          enddo
          do k=j+1,finish
            do l=1,0,-1
              a(l) = -mu(k)*a(l) + a(l-1)
            enddo
          enddo

          ! The coeffecients a(l) of the Lagrange polynomial are the
          ! derivation coefficients divided by l!.
          dmat%derive(i,i+j,1) = a(1)
        enddo
      enddo

      do i=1,ngrid-1

        start  = -order
        finish = order+1

        ! special treatment for the endpoints
        ! Note: this leads to results which are less precise on the
        ! edges, but preserves the banded form of the derivation matrix.
        if ((i+start).lt.1) start = 1 - i
        if ((i+finish).gt.ngrid) finish = ngrid - i

        do j=start,finish
          mu(j) = grid(j+i) - (grid(i)+grid(i+1))/2d0
        enddo

        do j=start,finish

          ! Initialise Lagrange polynomial
          a(-1) = 0d0
          a(0) = 1d0
          do l=1,2
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
            do l=2,0,-1
              a(l) = -mu(k)*a(l) + a(l-1)
            enddo
          enddo
          do k=j+1,finish
            do l=2,0,-1
              a(l) = -mu(k)*a(l) + a(l-1)
            enddo
          enddo

          ! The coeffecients a(l) of the Lagrange polynomial are the
          ! derivation coefficients divided by l! (where 0! = 1)
          dmat%derive(i,i+j,2) = a(2)*2d0
          if (der_max.eq.3) dmat%derive(i+1,i+j,3) = a(2)*2d0
        enddo
      enddo
      deallocate(a,mu)

      return
      end subroutine init_derive_MIX

!------------------------------------------------------------------------------ 
!  This module sets up a derivation matrices based on B-splines.
!  The unknowns in the equations will be the coefficients of the B-splines
!  rather than the values at different grid points.
!------------------------------------------------------------------------------ 
! Variables:
!
! ngrid    = number of grid points (1:ngrid)
! der_max  = maximal order of derivative
! order    = order of the B-splines (order >= der_max + 2)
!------------------------------------------------------------------------------ 
      subroutine init_derive_b_splines(dmat,grid,ngrid,der_max, &
                                       der_min,order)

      implicit none
      type(DERMAT), intent(out)    :: dmat
      integer, intent(in)          :: ngrid, der_max, der_min, order
      double precision, intent(in) :: grid(ngrid)
      integer i, j, der
      double precision, allocatable :: y(:,:), knots(:)

      ! check parameters:
      if (order.lt.der_max+2) then
        stop "Please increase order of B-splines"
      endif
      if (der_min.ne.0) then
       stop "Please set der_min = 0 in init_derive_b_splines"
      endif

      ! initialise some variables:
      allocate(y(order,0:der_max))
      y(:,:) = 0d0

      ! prepare knots for B-splines
      allocate(knots(ngrid+order))
      call find_knots(grid,ngrid,order,knots)

      ! initialise der_min, der_max and ngrid in dmat:
      dmat%der_min = der_min
      dmat%der_max = der_max
      dmat%der_id  = 0
      dmat%ngrid   = ngrid

      allocate(dmat%lbder(0:der_max),dmat%ubder(0:der_max))
      dmat%ubder(0) = 0
      dmat%lbder(0) = 0

      allocate(dmat%derive(1:ngrid,1:ngrid,0:der_max))
      dmat%derive(:,:,:) = 0d0

      ! create derivation-interpolation matrix:
      do i=1,ngrid

        ! We will consider groups of order+1 knots
        ! starting at j.  j will also index the
        ! B-splines in the derivation-interpolation
        ! matrix.
        do j=max(1,i-order),min(ngrid,i+order)

          if ((grid(i).ge.knots(j)).and.(grid(i).le.knots(j+order))) then
            ! calculate B-spline and derivatives:
            call find_b_spline(grid(i),knots(j:(j+order)), &
                               order,der_max,y,grid(1))
            do der = 0,der_max
              dmat%derive(i,j,der) = y(1,der)
            enddo
            if ((i-j).gt.dmat%lbder(0)) dmat%lbder(0) = i-j
            if ((j-i).gt.dmat%ubder(0)) dmat%ubder(0) = j-i
          endif
        enddo
      enddo

      ! fill in arrays lbder and ubder
      do der=1,der_max
        dmat%lbder(der) = dmat%lbder(0)
        dmat%ubder(der) = dmat%ubder(0)
      enddo

      deallocate(y)
      return
      end subroutine init_derive_b_splines

!---------------------------------------------------------------------
! This finds an adequate grid for setting up a a basis of B-splines,
! based on a formula by P. Morel (+ a minor correction).
!---------------------------------------------------------------------
! variables:
!
!   grid  = original grid
!   ngrid = number of points in grid
!   order = order of the B-splines (= degree + 1)
!   knots = output grid which will be used in setting up B-splines
!---------------------------------------------------------------------
        subroutine find_knots(grid,ngrid,order,knots)

        implicit none
        integer order, ngrid
        double precision grid(ngrid), knots(ngrid+order)
        integer i,j

        knots(1:order) = grid(1)

        do i=order+1,ngrid
          knots(i) = 0d0
          do j=1,order-1
            knots(i) = knots(i) + grid(i+j-order)
          enddo
          knots(i) = knots(i)/dble(order-1)
        enddo

        knots((ngrid+1):(ngrid+order)) = grid(ngrid)

        end subroutine

!---------------------------------------------------------------------
! This evaluates a B-spline and its derivatives at a point x.
!---------------------------------------------------------------------
! variables
!
!  x     = point at which to evaluate the B-spline and its derivatives
!  knots = points used in constructing the B-spline
!  order = order of the B-spline (= degree + 1)
!  dermax= maximum derivative (should obey dermax <= order-2)
!  y     = work space, and array in which the result is outputted:
!          y(1,0) = interpolation
!          y(1,1) = first derivative
!          y(1,2) = second derivative
!          ...
! left_boundary = left boundary of the domain.  This intervenes when
!                 calculating B-splines on that boundary.
!
!  Note: if x lies outside the grid knots, then y(*) = 0
!---------------------------------------------------------------------
       subroutine find_b_spline(x,knots,order,dermax,y,left_boundary)

       implicit none
       integer order, dermax
       double precision x, left_boundary
       double precision knots(order+1)
       double precision y(order,0:dermax) ! a work and solution space
       integer i,j,der,my_dermax

       ! check for easy exit:
       if ((x.lt.knots(1)).or.(x.gt.knots(order+1))) then
         y(:,:) = 0d0
         return
       endif

       ! check parameters:
       my_dermax = dermax
       if (dermax.gt.(order-2)) then
         print*,"Warning: dermax too large in find_b_splines."
         print*,"         Try increasing order of splines."
         my_dermax = order-2
       endif

       ! setting up initial conditions:
       y(:,:) = 0d0
       do i=1,order

         ! The strict inequality in the left-hand condition
         ! eliminates spurious peaks around the different
         ! knots.  However, B-splines on the left boundary
         ! need to be treated differently:

         if (knots(i).eq.left_boundary) then
           if ((knots(i).le.x).and.(x.le.knots(i+1))) then
             y(i,0) = 1d0
           else
             y(i,0) = 0d0
           endif
         else
           if ((knots(i).lt.x).and.(x.le.knots(i+1))) then
             y(i,0) = 1d0
           else
             y(i,0) = 0d0
           endif
         endif
       enddo

       ! Use recursion formula for B-splines.
       ! For derivatives: combine recursion with Leibniz' formula.

       do i=2,order
         do j=1,order+1-i
           if (knots(j).ne.knots(j+i-1)) then
             do der = my_dermax,1,-1
               y(j,der) = (y(j,der)*(x-knots(j))+dble(der)*y(j,der-1)) &
                        / (knots(j+i-1)-knots(j))
             enddo
             y(j,0) = y(j,0)*(x-knots(j))/(knots(j+i-1)-knots(j))
           else
             y(j,0:my_dermax) = 0d0
           endif
           if (knots(j+1).ne.knots(j+i)) then
             do der = my_dermax,1,-1
               y(j,der) = y(j,der) + (y(j+1,der)*(knots(j+i)-x) &
                        - dble(der)*y(j+1,der-1)) &
                        / (knots(j+i)-knots(j+1))
             enddo
             y(j,0) = y(j,0) + y(j+1,0)*(knots(j+i)-x)/(knots(j+i)-knots(j+1))
           endif
         enddo
       enddo

       return
       end subroutine

!------------------------------------------------------------------------------ 
      end module
