!-----------------------------------------------------------------------
! This module contains various subroutines for calculating matrices
! involved in calculating coupling coefficients, or projection
! matrices.  It uses efficient recursion formulas for calculating 
! the terms.
!-----------------------------------------------------------------------
      module fast_pylm

      double precision, private, parameter :: pi = 3.141592653589793d0

contains

!-----------------------------------------------------------------------
! This subroutine produces a 2D matrix which contains successive
! spherical harmonics, multiplied by a constant, on a user defined
! grid:
!    pylm(j,lndx(l)) = alpha.Y(l,mm,cth(j))
!
! alpha  = multiplicative constant
! pylm   = 2D matrix with successive
! nn_out = number of points on the user defined grid
! nn     = number of different spherical harmonics
! cth    = cos(theta(:)) where theta is the user-defined grid
! sth    = sin(theta(:)) where theta is the user-defined grid
! lndx   = l index array. lndx(l) says where to store the harmonic
!          of degree l in the pylm array.  lndx(l) = -1 means that
!          the given harmonic should not be stored in pylm.
! llmax  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine mat_ylm(alpha,pylm,nn_out,nn,cth,sth,lndx,llmax,mm)
 
      implicit none
      integer, intent(in) :: nn_out, nn, llmax, mm
      integer, intent(in) :: lndx(0:llmax)
      double precision, intent(in) :: alpha, cth(nn_out), sth(nn_out)
      double precision, intent(out) :: pylm(nn_out,nn)
      double precision, allocatable, dimension(:) :: plm1, plm2, plm3
      double precision cnst, coeff, coeff2
      integer j, l, ma

      ! initialisation
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out))

      ! preliminary settings
      ma = abs(mm)

      cnst = alpha
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma
      plm2 = plm1*dble(2*ma+1)*cth

      ! cnst is recycled, so no need to worry about alpha
      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (lndx(ma).ne.-1) then
        do j=1,nn_out
          pylm(j,lndx(ma)) = cnst*plm1(j)
        enddo
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (lndx(ma+1).ne.-1) then
        do j=1,nn_out
          pylm(j,lndx(ma+1)) = cnst*plm2(j)
        enddo
      endif

      do l = ma+2,llmax
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (lndx(l).ne.-1) then
          do j=1,nn_out
            pylm(j,lndx(l)) = cnst*plm3(j)
          enddo
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3)
      end subroutine mat_ylm

!-----------------------------------------------------------------------
! This subroutine produces a 2D matrix which contains successive
! spherical harmonics, multiplied by a constant, on a user defined
! grid:
!    pylm(j,lndx(l)) = alpha.dY(l,mm,cth(j))/dtheta
!
! alpha  = multiplicative constant
! pylm   = 2D matrix with successive
! nn_out = number of points on the user defined grid
! nn     = number of different spherical harmonics
! cth    = cos(theta(:)) where theta is the user-defined grid
! sth    = sin(theta(:)) where theta is the user-defined grid
! lndx   = l index array. lndx(l) says where to store the harmonic
!          of degree l in the pylm array.  lndx(l) = -1 means that
!          the given harmonic should not be stored in pylm.
! llmax  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine mat_dth_ylm(alpha,pylm,nn_out,nn,cth,sth,lndx,llmax,mm)

      implicit none
      integer, intent(in) :: nn_out, nn, llmax, mm
      integer, intent(in) :: lndx(0:llmax)
      double precision, intent(in) :: alpha, cth(nn_out), sth(nn_out)
      double precision, intent(out) :: pylm(nn_out,nn)
      double precision, allocatable, dimension(:) :: &
        plm1, plm2, plm3, plmp1, plmp2, plmp3
      double precision cnst, coeff, coeff2
      integer j, l, ma

      ! workspace and initialisation:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),plmp1(nn_out), &
               plmp2(nn_out),plmp3(nn_out))

      ! preliminary settings
      ma = abs(mm)

      cnst = alpha
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma
      plm2 = plm1*cth*dble(2*ma+1)

      if (ma.eq.0) then
        plmp1 = 0d0
      else
        plmp1 = cnst*dble(ma)*sth**(ma-1)*cth
      endif
      plmp2 = (plmp1*cth-plm1*sth)*dble(2*ma+1)

      ! cnst is recycled, so no need to worry about alpha
      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (lndx(ma).ne.-1) then
        do j=1,nn_out
          pylm(j,lndx(ma)) = cnst*plmp1(j)
        enddo
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (lndx(ma+1).ne.-1) then
        do j=1,nn_out
          pylm(j,lndx(ma+1)) = cnst*plmp2(j)
        enddo
      endif

      do l = ma+2,llmax
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (lndx(l).ne.-1) then
          do j=1,nn_out
            pylm(j,lndx(l)) = cnst*plmp3(j)
          enddo
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3)
      end subroutine mat_dth_ylm

!-----------------------------------------------------------------------
! This subroutine produces a 2D matrix which contains successive
! spherical harmonics, multiplied by a constant, on a user defined
! grid:
!    pylm(j,lndx(l)) = alpha.(dY(l,mm,cth(j))/dphi)/sin(theta)
!
! alpha  = multiplicative constant
! pylm   = 2D matrix with successive
! nn_out = number of points on the user defined grid
! nn     = number of different spherical harmonics
! cth    = cos(theta(:)) where theta is the user-defined grid
! sth    = sin(theta(:)) where theta is the user-defined grid
! lndx   = l index array. lndx(l) says where to store the harmonic
!          of degree l in the pylm array.  lndx(l) = -1 means that
!          the given harmonic should not be stored in pylm.
! llmax  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine mat_Dphi_ylm(alpha,pylm,nn_out,nn,cth,sth,lndx,llmax,mm)

      implicit none
      integer, intent(in) :: nn_out, nn, llmax, mm
      integer, intent(in) :: lndx(0:llmax)
      double precision, intent(in) :: alpha, cth(nn_out), sth(nn_out)
      double precision, intent(out) :: pylm(nn_out,nn)
      double precision, allocatable, dimension(:) :: plm1,plm2,plm3
      double precision cnst, coeff, coeff2
      integer j, l, ma

      ! easy exit, if possible:
      if (mm.eq.0) then
        pylm = 0d0
        return
      endif

      ! get workspace:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out))

      ! preliminary settings
      ma = abs(mm)

      cnst = alpha
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = dble(mm)*cnst*sth**(ma-1)
      plm2 = plm1*dble(2*ma+1)*cth

      ! cnst is recycled, so no need to worry about alpha
      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (lndx(ma).ne.-1) then
        do j=1,nn_out
          pylm(j,lndx(ma)) = cnst*plm1(j)
        enddo
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (lndx(ma+1).ne.-1) then
        do j=1,nn_out
          pylm(j,lndx(ma+1)) = cnst*plm2(j)
        enddo
      endif

      do l = ma+2,llmax
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (lndx(l).ne.-1) then
          do j=1,nn_out
            pylm(j,lndx(l)) = cnst*plm3(j)
          enddo
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3)
      end subroutine mat_Dphi_ylm

!-----------------------------------------------------------------------
! This subroutine produces a 2D matrix which contains successive
! spherical harmonics, multiplied by a constant, on a user defined
! grid:
!    pylm(j,lndx(l)) = alpha.d^2 Y(l,mm,cth(j))/d(theta)^2
!
! alpha  = multiplicative constant
! pylm   = 2D matrix with successive
! nn_out = number of points on the user defined grid
! nn     = number of different spherical harmonics
! cth    = cos(theta(:)) where theta is the user-defined grid
! sth    = sin(theta(:)) where theta is the user-defined grid
! lndx   = l index array. lndx(l) says where to store the harmonic
!          of degree l in the pylm array.  lndx(l) = -1 means that
!          the given harmonic should not be stored in pylm.
! llmax  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine mat_dtt_ylm(alpha,pylm,nn_out,nn,cth,sth,lndx,llmax,mm)

      implicit none
      integer, intent(in) :: nn_out, nn, llmax, mm
      integer, intent(in) :: lndx(0:llmax)
      double precision, intent(in) :: alpha, cth(nn_out), sth(nn_out)
      double precision, intent(out) :: pylm(nn_out,nn)
      double precision, allocatable, dimension(:) :: &
        plm1, plm2, plm3, plmp1, plmp2, plmp3, plmpp1, plmpp2, plmpp3
      double precision cnst, coeff, coeff2
      integer j, l, ma

      ! workspace and initialisation:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),plmp1(nn_out), &
               plmp2(nn_out),plmp3(nn_out),plmpp1(nn_out),           &
               plmpp2(nn_out),plmpp3(nn_out))

      ! preliminary settings
      ma = abs(mm)

      cnst = alpha
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma
      if (ma.eq.0) then
        plmp1 = 0d0
        plmpp1 = 0d0
      else
        plmp1 = cnst*dble(ma)*sth**(ma-1)*cth
        plmpp1 = cnst*(dble(ma*(ma-1))*sth**(ma-2)-dble(ma*ma)*sth**ma)
      endif

      plm2 = plm1*cth*dble(2*ma+1)
      plmp2 = (plmp1*cth-plm1*sth)*dble(2*ma+1)
      plmpp2 = ((plmpp1-plm1)*cth-2d0*plmp1*sth)*dble(2*ma+1)

      ! cnst is recycled, so no need to worry about alpha
      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (lndx(ma).ne.-1) then
        do j=1,nn_out
          pylm(j,lndx(ma)) = cnst*plmpp1(j)
        enddo
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (lndx(ma+1).ne.-1) then
        do j=1,nn_out
          pylm(j,lndx(ma+1)) = cnst*plmpp2(j)
        enddo
      endif

      do l = ma+2,llmax
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        plmpp3= coeff2*(cth*(plmpp2-plm2)-2d0*sth*plmp2) - coeff*plmpp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (lndx(l).ne.-1) then
          do j=1,nn_out
            pylm(j,lndx(l)) = cnst*plmpp3(j)
          enddo
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
        plmpp1 = plmpp2
        plmpp2 = plmpp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3,plmpp1,plmpp2,plmpp3)
      end subroutine mat_dtt_ylm

!-----------------------------------------------------------------------
! This subroutine produces a 2D matrix which contains successive
! spherical harmonics, multiplied by a constant, on a user defined
! grid:
!    pylm(j,lndx(l)) = alpha.d{[dY(l,mm,cth(j))/dphi]/sin(theta)}/d(theta)
!
! alpha  = multiplicative constant
! pylm   = 2D matrix with successive
! nn_out = number of points on the user defined grid
! nn     = number of different spherical harmonics
! cth    = cos(theta(:)) where theta is the user-defined grid
! sth    = sin(theta(:)) where theta is the user-defined grid
! lndx   = l index array. lndx(l) says where to store the harmonic
!          of degree l in the pylm array.  lndx(l) = -1 means that
!          the given harmonic should not be stored in pylm.
! llmax  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine mat_dtp_ylm(alpha,pylm,nn_out,nn,cth,sth,lndx,llmax,mm)

      implicit none
      integer, intent(in) :: nn_out, nn, llmax, mm
      integer, intent(in) :: lndx(0:llmax)
      double precision, intent(in) :: alpha, cth(nn_out), sth(nn_out)
      double precision, intent(out) :: pylm(nn_out,nn)
      double precision, allocatable, dimension(:) :: plm1,plm2,plm3,&
        plmp1,plmp2,plmp3
      double precision cnst, coeff, coeff2
      integer j, l, ma

      ! easy exit, if possible:
      if (mm.eq.0) then
        pylm = 0d0
        return
      endif

      ! get workspace:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),plmp1(nn_out), &
               plmp2(nn_out),plmp3(nn_out))

      ! preliminary settings
      ma = abs(mm)

      cnst = alpha
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*dble(mm)*sth**(ma-1)
      if (ma.le.1) then
        plmp1 = 0d0
      else
        plmp1 = cnst*dble(mm*(ma-1))*sth**(ma-2)*cth 
      endif
      plm2 = dble(2*ma+1)*plm1*cth
      plmp2= dble(2*ma+1)*(plmp1*cth-plm1*sth)

      ! cnst is recycled, so no need to worry about alpha
      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (lndx(ma).ne.-1) then
        do j=1,nn_out
          pylm(j,lndx(ma)) = cnst*plmp1(j)
        enddo
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (lndx(ma+1).ne.-1) then
        do j=1,nn_out
          pylm(j,lndx(ma+1)) = cnst*plmp2(j)
        enddo
      endif

      do l = ma+2,llmax
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (lndx(l).ne.-1) then
          do j=1,nn_out
            pylm(j,lndx(l)) = cnst*plmp3(j)
          enddo
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3)
      end subroutine mat_dtp_ylm
!-----------------------------------------------------------------------
      end module fast_pylm
