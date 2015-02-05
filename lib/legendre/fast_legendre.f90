!-----------------------------------------------------------------------
! This module contains various subroutines for projecting functions
! onto the harmonic basis or inversely.  It uses efficient recursion
! formulas for calculating the spherical harmonics as it's doing the
! projection or evaluation.
!-----------------------------------------------------------------------
! VERY IMPORTANT !!!
! ALL of these subroutines assume that the result is ADDED to the output
! variable, rather than written over the original content.  This allows
! users to make calculations with several subroutines in a row and add
! the results together in a natural way.
!-----------------------------------------------------------------------
      module fast_legendre

      double precision, private, parameter :: pi = 3.141592653589793d0

      interface eval_ylm
        module procedure eval_ylm_1D
        module procedure eval_ylm_2D
        module procedure eval_ylm_1D_z
        module procedure eval_ylm_2D_z
      end interface

      interface eval_dth_ylm
        module procedure eval_dth_ylm_1D
        module procedure eval_dth_ylm_2D
        module procedure eval_dth_ylm_1D_z
        module procedure eval_dth_ylm_2D_z
      end interface

      interface eval_Dphi_ylm
        module procedure eval_Dphi_ylm_1D
        module procedure eval_Dphi_ylm_2D
        module procedure eval_Dphi_ylm_1D_z
        module procedure eval_Dphi_ylm_2D_z
      end interface

      interface eval_ddth_ylm
        module procedure eval_ddth_ylm_1D
        module procedure eval_ddth_ylm_2D
        module procedure eval_ddth_ylm_1D_z
        module procedure eval_ddth_ylm_2D_z
      end interface

      interface project_ylm
        module procedure project_ylm_1D
        module procedure project_ylm_2D
        module procedure project_ylm_1D_z
        module procedure project_ylm_2D_z
      end interface

      interface project_dth_ylm
        module procedure project_dth_ylm_1D
        module procedure project_dth_ylm_2D
        module procedure project_dth_ylm_1D_z
        module procedure project_dth_ylm_2D_z
      end interface

      interface project_Dphi_ylm
        module procedure project_Dphi_ylm_1D
        module procedure project_Dphi_ylm_2D
        module procedure project_Dphi_ylm_1D_z
        module procedure project_Dphi_ylm_2D_z
      end interface

contains
!-----------------------------------------------------------------------
! This subroutine evaluates a function decomposed over the
! spherical harmonic basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l Ylm
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_ylm_1D(alpha,fs,f,cth,nn,nn_out,lstart,lincr,mm)
 
      implicit none
      integer, intent(in) ::  nn, nn_out, lstart, lincr, mm
      double precision, intent(in)   :: alpha, fs(nn), cth(nn_out)
      double precision, intent(inout):: f(nn_out)
      double precision, allocatable, dimension(:) :: plm1, plm2, plm3,&
                                                     sth
      double precision cnst, coeff, coeff2
      integer l, ma, ll

      ! initialisation
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),sth(nn_out)) 
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          coeff = cnst*alpha*fs(ll)
          f = f + coeff*plm1
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          coeff = cnst*alpha*fs(ll)
          f = f + coeff*plm2
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            coeff = cnst*alpha*fs(ll)
            f = f + coeff*plm3
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth) 
      end subroutine eval_ylm_1D

!-----------------------------------------------------------------------
! This subroutine evaluates the theta derivative of a function 
! decomposed over the spherical harmonic basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(theta)
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_dth_ylm_1D(alpha,fs,f,cth,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nn, nn_out, lstart, lincr, mm
      double precision, intent(in)   :: alpha,fs(nn),cth(nn_out)
      double precision, intent(inout):: f(nn_out)
      double precision, allocatable, dimension(:) :: &
        plm1, plm2, plm3, plmp1, plmp2, plmp3, sth
      double precision cnst, coeff, coeff2
      integer l, ma, ll

      ! workspace and initialisation:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),plmp1(nn_out), &
               plmp2(nn_out),plmp3(nn_out),sth(nn_out))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
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

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if ((ma.gt.0).and.(modulo(ma-lstart,lincr).eq.0)) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          coeff = cnst*alpha*fs(ll)
          f = f + coeff*plmp1
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          coeff = cnst*alpha*fs(ll)
          f = f + coeff*plmp2
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            coeff = cnst*alpha*fs(ll)
            f = f + coeff*plmp3
          endif
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3,sth)
      end subroutine eval_dth_ylm_1D

!-----------------------------------------------------------------------
! This subroutine evaluates the imaginary part of the phi derivative of
! a function decomposed over the spherical harmonic basis, divided by 
! sin(theta):
!    f(theta) =  alpha.Im { sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(phi)/sin(theta) }
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_Dphi_ylm_1D(alpha,fs,f,cth,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nn, nn_out, lstart, lincr, mm
      double precision, intent(in)   :: alpha,fs(nn),cth(nn_out)
      double precision, intent(inout):: f(nn_out)
      double precision, allocatable, dimension(:) :: plm1,plm2,plm3,sth
      double precision cnst, coeff, coeff2
      integer l, ma, ll

      ! easy exit, if possible:
      if (mm.eq.0) return

      ! get workspace:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),sth(nn_out))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = dble(mm)*cnst*sth**(ma-1)
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          coeff = cnst*alpha*fs(ll)
          f = f + coeff*plm1
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          coeff = cnst*alpha*fs(ll)
          f = f + coeff*plm2
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            coeff = cnst*alpha*fs(ll)
            f = f + coeff*plm3
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth)
      end subroutine eval_Dphi_ylm_1D

!-----------------------------------------------------------------------
! This subroutine evaluates the theta derivative of a function 
! decomposed over the spherical harmonic basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(theta)
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_ddth_ylm_1D(alpha,fs,f,cth,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nn, nn_out, lstart, lincr, mm
      double precision, intent(in)   :: alpha,fs(nn),cth(nn_out)
      double precision, intent(inout):: f(nn_out)
      double precision, allocatable, dimension(:) :: sth, &
        plm1, plm2, plm3, plmp1, plmp2, plmp3, plmpp1, plmpp2, plmpp3
      double precision cnst, coeff, coeff2
      integer l, ma, ll

      ! workspace and initialisation:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),plmp1(nn_out), &
               plmp2(nn_out),plmp3(nn_out),plmpp1(nn_out),           &
               plmpp2(nn_out),plmpp3(nn_out),sth(nn_out))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma
      plm2 = plm1*cth*dble(2*ma+1)

      if (ma.eq.0) then
        plmp1 = 0d0
        plmpp1 = 0d0
      else
        plmp1 = cnst*dble(ma)*sth**(ma-1)*cth
        if (ma.eq.1) then
          plmpp1 = -cnst*sth
        else
          plmpp1 =  cnst*ma*((ma-1)*sth**(ma-2)-ma*sth**ma)
        endif
      endif
      plmp2 = (plmp1*cth-plm1*sth)*dble(2*ma+1)
      plmpp2 = (plmpp1*cth-2d0*plmp1*sth-plm1*cth)*dble(2*ma+1)

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if ((ma.gt.0).and.(modulo(ma-lstart,lincr).eq.0)) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          coeff = cnst*alpha*fs(ll)
          f = f + coeff*plmpp1
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          coeff = cnst*alpha*fs(ll)
          f = f + coeff*plmpp2
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        plmpp3= coeff2*(cth*plmpp2-2d0*sth*plmp2-cth*plm2) - coeff*plmpp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            coeff = cnst*alpha*fs(ll)
            f = f + coeff*plmpp3
          endif
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
        plmpp1 = plmpp2
        plmpp2 = plmpp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3,plmpp1,plmpp2,plmpp3,sth)
      end subroutine eval_ddth_ylm_1D

!-----------------------------------------------------------------------
! This subroutine projects a function onto the spherical harmonic basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn_out-1)} fs_l Ylm
!
! alpha  = multiplicative constant
! f      = function in real space
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! w      = integration weights associated with the grid theta(:)
! fs     = spectral representation of the function (= OUTPUT)
! nn     = angular resolution in real space
! nn_out = number of spherical harmonics in output
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine project_ylm_1D(alpha,f,cth,w,fs,nn,nn_out,lstart,lincr,mm)
 
      implicit none
      integer, intent(in) ::  nn, nn_out, lstart, lincr, mm
      double precision, intent(in)   :: alpha,f(nn),cth(nn),w(nn)
      double precision, intent(inout):: fs(nn_out)
      double precision, allocatable, dimension(:) :: plm1,plm2,plm3,sth
      double precision cnst, coeff, coeff2, aux
      integer l, ma, ll, j

      ! initialisation
      allocate(plm1(nn),plm2(nn),plm3(nn),sth(nn))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma*w*f*(2d0*pi)
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = 0d0
          do j=1,nn
            aux = aux + plm1(j)
          enddo
          fs(ll) = fs(ll) + alpha*aux*cnst
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = 0d0
          do j=1,nn
            aux = aux + plm2(j)
          enddo
          fs(ll) = fs(ll) + alpha*aux*cnst
        endif
      endif

      do l = ma+2,lstart+lincr*(nn_out-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.gt.0) then
            aux = 0d0
            do j=1,nn
              aux = aux + plm3(j)
            enddo
            fs(ll) = fs(ll) + alpha*aux*cnst
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth) 
      end subroutine project_ylm_1D

!-----------------------------------------------------------------------
! This subroutine projects a function onto the theta derivatives of the
! the spherical harmonics:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(theta)
!
! alpha  = multiplicative constant
! f      = function in real space
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! w      = integration weights associated with the grid theta(:)
! fs     = spectral representation of the function (= OUTPUT)
! nn     = angular resolution in real space
! nn_out = number of spherical harmonics in output
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine project_dth_ylm_1D(alpha,f,cth,w,fs,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nn, nn_out, lstart, lincr, mm
      double precision, intent(in)   :: alpha,f(nn),cth(nn),w(nn)
      double precision, intent(inout):: fs(nn_out)
      double precision, allocatable, dimension(:) :: &
        plm1, plm2, plm3, plmp1, plmp2, plmp3, sth
      double precision cnst, coeff, coeff2, aux
      integer l, ma, ll, j

      ! get workspace:
      allocate(plm1(nn),plm2(nn),plm3(nn),plmp1(nn), &
               plmp2(nn),plmp3(nn),sth(nn))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo
 
      plm1 = (cnst*2d0*pi)*sth**ma*f*w
      plm2 = plm1*cth*dble(2*ma+1)

      if (ma.eq.0) then
        plmp1 = 0d0
      else
        plmp1 = (cnst*dble(ma)*2d0*pi)*sth**(ma-1)*cth*f*w
      endif
      plmp2 = (plmp1*cth-plm1*sth)*dble(2*ma+1)

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if ((ma.gt.0).and.(modulo(ma-lstart,lincr).eq.0)) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = 0d0
          do j=1,nn
            aux = aux + plmp1(j)
          enddo
          fs(ll) = fs(ll) + alpha*aux*cnst/dble(ma*(ma+1))
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = 0d0
          do j=1,nn
            aux = aux + plmp2(j)
          enddo
          fs(ll) = fs(ll) + alpha*aux*cnst/dble((ma+1)*(ma+2))
        endif
      endif

      do l = ma+2,lstart+lincr*(nn_out-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            aux = 0d0
            do j=1,nn
              aux = aux + plmp3(j)
            enddo
            fs(ll) = fs(ll) + alpha*aux*cnst/dble(l*(l+1))
          endif
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3,sth)
      end subroutine project_dth_ylm_1D

!-----------------------------------------------------------------------
! This subroutine projects a function onto the imaginary part of the 
! phi derivatives of the spherical harmonics, divided by sin(theta):
!    f(theta) =  alpha.Im { sum_{l=lstart}^{lstart+2(nn-1)} fs_l [d(Ylm)/d(phi)]/sin(theta) }
!
! alpha  = multiplicative constant
! f      = function in real space
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! w      = integration weights associated with the grid theta(:)
! fs     = spectral representation of the function (= OUTPUT)
! nn     = angular resolution in real space
! nn_out = number of spherical harmonics in output
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine project_Dphi_ylm_1D(alpha,f,cth,w,fs,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nn, nn_out, lstart, lincr, mm
      double precision, intent(in)   :: alpha,f(nn),cth(nn),w(nn)
      double precision, intent(inout):: fs(nn_out)
      double precision, allocatable, dimension(:) :: plm1,plm2,plm3,sth
      double precision cnst, coeff, coeff2, aux
      integer l, ma, ll, j

      ! easy exit, if possible:
      if (mm.eq.0) return

      ! get workspace:
      allocate(plm1(nn),plm2(nn),plm3(nn),sth(nn))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo
 
      plm1 = dble(mm)*cnst*sth**(ma-1)*f*w*(2d0*pi)
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = (4d0*pi)/dble(2*ma+1)
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(1d0/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = 0d0
          do j=1,nn
            aux = aux + plm1(j)
          enddo
          fs(ll) = fs(ll) + alpha*aux*cnst/dble(ma*(ma+1))
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = 0d0
          do j=1,nn
            aux = aux + plm2(j)
          enddo
          fs(ll) = fs(ll) + alpha*aux*cnst/dble((ma+1)*(ma+2))
        endif
      endif

      do l = ma+2,lstart+lincr*(nn_out-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.gt.0) then
            aux = 0d0
            do j=1,nn
              aux = aux + plm3(j)
            enddo
            fs(ll) = fs(ll) + alpha*aux*cnst/dble(l*(l+1))
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth)
      end subroutine project_Dphi_ylm_1D

!-----------------------------------------------------------------------
! This subroutine evaluates a 2D function decomposed over the
! spherical harmonic basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l Ylm
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nnr    = number of grid points in extra dimension
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_ylm_2D(alpha,fs,f,cth,nnr,nn,nn_out,lstart,lincr,mm)
 
      implicit none
      integer, intent(in) ::  nnr, nn, nn_out, lstart, lincr, mm
      double precision, intent(in)   :: alpha,fs(nnr,nn), cth(nn_out)
      double precision, intent(inout):: f(nnr,nn_out)
      double precision, allocatable, dimension(:) :: plm1, plm2, plm3,&
                                                     sth
      double precision cnst, coeff, coeff2
      integer i, l, ma, ll

      ! initialisation
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),sth(nn_out)) 
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            coeff = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + coeff*plm1
          enddo
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            coeff = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + coeff*plm2
          enddo
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            do i=1,nnr
              coeff = cnst*alpha*fs(i,ll)
              f(i,:) = f(i,:) + coeff*plm3
            enddo
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth) 
      end subroutine eval_ylm_2D

!-----------------------------------------------------------------------
! This subroutine evaluates the theta derivative of a 2D function 
! decomposed over the spherical harmonic basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(theta)
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nnr    = number of grid points in extra dimension
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_dth_ylm_2D(alpha,fs,f,cth,nnr,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nnr, nn, nn_out, lstart, lincr, mm
      double precision, intent(in)   :: alpha,fs(nnr,nn),cth(nn_out)
      double precision, intent(inout):: f(nnr,nn_out)
      double precision, allocatable, dimension(:) :: &
        plm1, plm2, plm3, plmp1, plmp2, plmp3, sth
      double precision cnst, coeff, coeff2
      integer i, l, ma, ll

      ! workspace and initialisation:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),plmp1(nn_out), &
               plmp2(nn_out),plmp3(nn_out),sth(nn_out))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
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

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if ((ma.gt.0).and.(modulo(ma-lstart,lincr).eq.0)) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            coeff = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + coeff*plmp1
          enddo
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            coeff = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + coeff*plmp2
          enddo
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            do i=1,nnr
              coeff = cnst*alpha*fs(i,ll)
              f(i,:) = f(i,:) + coeff*plmp3
            enddo
          endif
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3,sth)
      end subroutine eval_dth_ylm_2D

!-----------------------------------------------------------------------
! This subroutine evaluates the imaginary part of the phi derivative of
! a 2D function decomposed over the spherical harmonic basis, divided by 
! sin(theta):
!    f(theta) =  alpha.Im { sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(phi)/sin(theta) }
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nnr    = number of grid points in extra dimension
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_Dphi_ylm_2D(alpha,fs,f,cth,nnr,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nnr, nn, nn_out, lstart, lincr, mm
      double precision, intent(in)   :: alpha,fs(nnr,nn),cth(nn_out)
      double precision, intent(inout):: f(nnr,nn_out)
      double precision, allocatable, dimension(:) :: plm1,plm2,plm3,sth
      double precision cnst, coeff, coeff2
      integer i, l, ma, ll

      ! easy exit, if possible:
      if (mm.eq.0) return

      ! get workspace:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),sth(nn_out))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = dble(mm)*cnst*sth**(ma-1)
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            coeff = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + coeff*plm1
          enddo
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            coeff = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + coeff*plm2
          enddo
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            do i=1,nnr
              coeff = cnst*alpha*fs(i,ll)
              f(i,:) = f(i,:) + coeff*plm3
            enddo
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth)
      end subroutine eval_Dphi_ylm_2D

!-----------------------------------------------------------------------
! This subroutine evaluates the theta derivative of a 2D function 
! decomposed over the spherical harmonic basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(theta)
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nnr    = number of grid points in extra dimension
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_ddth_ylm_2D(alpha,fs,f,cth,nnr,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nnr, nn, nn_out, lstart, lincr, mm
      double precision, intent(in)   :: alpha,fs(nnr,nn),cth(nn_out)
      double precision, intent(inout):: f(nnr,nn_out)
      double precision, allocatable, dimension(:) :: sth, &
        plm1, plm2, plm3, plmp1, plmp2, plmp3, plmpp1, plmpp2, plmpp3
      double precision cnst, coeff, coeff2
      integer i, l, ma, ll

      ! workspace and initialisation:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),plmp1(nn_out), &
               plmp2(nn_out),plmp3(nn_out),plmpp1(nn_out),           &
               plmpp2(nn_out),plmpp3(nn_out),sth(nn_out))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma
      plm2 = plm1*cth*dble(2*ma+1)

      if (ma.eq.0) then
        plmp1 = 0d0
        plmpp1 = 0d0
      else
        plmp1 = cnst*dble(ma)*sth**(ma-1)*cth
        if (ma.eq.1) then
          plmpp1 = -cnst*sth
        else
          plmpp1 =  cnst*ma*((ma-1)*sth**(ma-2)-ma*sth**ma)
        endif
      endif
      plmp2 = (plmp1*cth-plm1*sth)*dble(2*ma+1)
      plmpp2 = (plmpp1*cth-2d0*plmp1*sth-plm1*cth)*dble(2*ma+1)

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if ((ma.gt.0).and.(modulo(ma-lstart,lincr).eq.0)) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            coeff = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + coeff*plmpp1
          enddo
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            coeff = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + coeff*plmpp2
          enddo
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        plmpp3= coeff2*(cth*plmpp2-2d0*sth*plmp2-cth*plm2) - coeff*plmpp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            do i=1,nnr
              coeff = cnst*alpha*fs(i,ll)
              f(i,:) = f(i,:) + coeff*plmpp3
            enddo
          endif
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
        plmpp1 = plmpp2
        plmpp2 = plmpp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3,plmpp1,plmpp2,plmpp3,sth)
      end subroutine eval_ddth_ylm_2D

!-----------------------------------------------------------------------
! This subroutine projects a 2D function onto the spherical harmonic
! basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn_out-1)} fs_l Ylm
!
! alpha  = multiplicative constant
! f      = function in real space
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! w      = integration weights associated with the grid theta(:)
! fs     = spectral representation of the function (= OUTPUT)
! nnr    = number of grid points in extra dimension
! nn     = angular resolution in real space
! nn_out = number of spherical harmonics in output
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine project_ylm_2D(alpha,f,cth,w,fs,nnr,nn,nn_out,lstart,lincr,mm)
 
      implicit none
      integer, intent(in) ::  nnr, nn, nn_out, lstart, lincr, mm
      double precision, intent(in)   :: alpha,f(nnr,nn),cth(nn),w(nn)
      double precision, intent(inout):: fs(nnr,nn_out)
      double precision, allocatable, dimension(:) :: plm1,plm2,plm3,sth
      double precision cnst, coeff, coeff2, aux
      integer l, ma, ll, i, j

      ! initialisation
      allocate(plm1(nn),plm2(nn),plm3(nn),sth(nn))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma*w*(2d0*pi)
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = 0d0
            do j=1,nn
              aux = aux + plm1(j)*f(i,j)
            enddo
            fs(i,ll) = fs(i,ll) + alpha*aux*cnst
          enddo
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = 0d0
            do j=1,nn
              aux = aux + plm2(j)*f(i,j)
            enddo
            fs(i,ll) = fs(i,ll) + alpha*aux*cnst
          enddo
        endif
      endif

      do l = ma+2,lstart+lincr*(nn_out-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.gt.0) then
            do i=1,nnr
              aux = 0d0
              do j=1,nn
                aux = aux + plm3(j)*f(i,j)
              enddo
              fs(i,ll) = fs(i,ll) + alpha*aux*cnst
            enddo
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth) 
      end subroutine project_ylm_2D

!-----------------------------------------------------------------------
! This subroutine projects a function onto the theta derivatives of the
! the spherical harmonics:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(theta)
!
! alpha  = multiplicative constant
! f      = function in real space
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! w      = integration weights associated with the grid theta(:)
! fs     = spectral representation of the function (= OUTPUT)
! nnr    = number of grid points in extra dimension
! nn     = angular resolution in real space
! nn_out = number of spherical harmonics in output
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine project_dth_ylm_2D(alpha,f,cth,w,fs,nnr,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nnr, nn, nn_out, lstart, lincr, mm
      double precision, intent(in)   :: alpha,f(nnr,nn),cth(nn),w(nn)
      double precision, intent(inout):: fs(nnr,nn_out)
      double precision, allocatable, dimension(:) :: &
        plm1, plm2, plm3, plmp1, plmp2, plmp3, sth
      double precision cnst, coeff, coeff2, aux
      integer l, ma, ll, i, j

      ! get workspace:
      allocate(plm1(nn),plm2(nn),plm3(nn),plmp1(nn), &
               plmp2(nn),plmp3(nn),sth(nn))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo
 
      plm1 = (cnst*2d0*pi)*sth**ma*w
      plm2 = plm1*cth*dble(2*ma+1)

      if (ma.eq.0) then
        plmp1 = 0d0
      else
        plmp1 = (cnst*dble(ma)*2d0*pi)*sth**(ma-1)*cth*w
      endif
      plmp2 = (plmp1*cth-plm1*sth)*dble(2*ma+1)

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if ((ma.gt.0).and.(modulo(ma-lstart,lincr).eq.0)) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = 0d0
            do j=1,nn
              aux = aux + plmp1(j)*f(i,j)
            enddo
            fs(i,ll) = fs(i,ll) + alpha*aux*cnst/dble(ma*(ma+1))
          enddo
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = 0d0
            do j=1,nn
              aux = aux + plmp2(j)*f(i,j)
            enddo
            fs(i,ll) = fs(i,ll) + alpha*aux*cnst/dble((ma+1)*(ma+2))
          enddo
        endif
      endif

      do l = ma+2,lstart+lincr*(nn_out-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            do i=1,nnr
              aux = 0d0
              do j=1,nn
                aux = aux + plmp3(j)*f(i,j)
              enddo
              fs(i,ll) = fs(i,ll) + alpha*aux*cnst/dble(l*(l+1))
            enddo
          endif
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3,sth)
      end subroutine project_dth_ylm_2D

!-----------------------------------------------------------------------
! This subroutine projects a function onto the imaginary part of the 
! phi derivatives of the spherical harmonics, divided by sin(theta):
!    f(theta) =  alpha.Im { sum_{l=lstart}^{lstart+2(nn-1)} fs_l [d(Ylm)/d(phi)]/sin(theta) }
!
! alpha  = multiplicative constant
! f      = function in real space
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! w      = integration weights associated with the grid theta(:)
! fs     = spectral representation of the function (= OUTPUT)
! nnr    = number of grid points in extra dimension
! nn     = angular resolution in real space
! nn_out = number of spherical harmonics in output
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine project_Dphi_ylm_2D(alpha,f,cth,w,fs,nnr,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nnr, nn, nn_out, lstart, lincr, mm
      double precision, intent(in)   :: alpha,f(nnr,nn),cth(nn),w(nn)
      double precision, intent(inout):: fs(nnr,nn_out)
      double precision, allocatable, dimension(:) :: plm1,plm2,plm3,sth
      double precision cnst, coeff, coeff2, aux
      integer l, ma, ll, i, j

      ! easy exit, if possible:
      if (mm.eq.0) return

      ! get workspace:
      allocate(plm1(nn),plm2(nn),plm3(nn),sth(nn))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo
 
      plm1 = dble(mm)*cnst*sth**(ma-1)*w*(2d0*pi)
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = (4d0*pi)/dble(2*ma+1)
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(1d0/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = 0d0
            do j=1,nn
              aux = aux + plm1(j)*f(i,j)
            enddo
            fs(i,ll) = fs(i,ll) + alpha*aux*cnst/dble(ma*(ma+1))
          enddo
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = 0d0
            do j=1,nn
              aux = aux + plm2(j)*f(i,j)
            enddo
            fs(i,ll) = fs(i,ll) + alpha*aux*cnst/dble((ma+1)*(ma+2))
          enddo
        endif
      endif

      do l = ma+2,lstart+lincr*(nn_out-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.gt.0) then
            do i=1,nnr
              aux = 0d0
              do j=1,nn
                aux = aux + plm3(j)*f(i,j)
              enddo
              fs(i,ll) = fs(i,ll) + alpha*aux*cnst/dble(l*(l+1))
            enddo
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth)
      end subroutine project_Dphi_ylm_2D

!-----------------------------------------------------------------------
! This subroutine evaluates a function decomposed over the
! spherical harmonic basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l Ylm
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_ylm_1D_z(alpha,fs,f,cth,nn,nn_out,lstart,lincr,mm)
 
      implicit none
      integer, intent(in) ::  nn, nn_out, lstart, lincr, mm
      double precision, intent(in) :: cth(nn_out)
      double complex, intent(in)   :: alpha,fs(nn)
      double complex, intent(inout):: f(nn_out)
      double precision, allocatable, dimension(:) :: plm1, plm2, plm3,&
                                                     sth
      double precision cnst, coeff, coeff2
      double complex aux
      integer l, ma, ll

      ! initialisation
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),sth(nn_out)) 
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = cnst*alpha*fs(ll)
          f = f + aux*plm1
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = cnst*alpha*fs(ll)
          f = f + aux*plm2
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            aux = cnst*alpha*fs(ll)
            f = f + aux*plm3
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth) 
      end subroutine eval_ylm_1D_z

!-----------------------------------------------------------------------
! This subroutine evaluates the theta derivative of a function 
! decomposed over the spherical harmonic basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(theta)
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_dth_ylm_1D_z(alpha,fs,f,cth,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nn, nn_out, lstart, lincr, mm
      double precision, intent(in) :: cth(nn_out)
      double complex, intent(in)   :: alpha,fs(nn)
      double complex, intent(inout):: f(nn_out)
      double precision, allocatable, dimension(:) :: &
        plm1, plm2, plm3, plmp1, plmp2, plmp3, sth
      double precision cnst, coeff, coeff2
      double complex aux
      integer l, ma, ll

      ! workspace and initialisation:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),plmp1(nn_out), &
               plmp2(nn_out),plmp3(nn_out),sth(nn_out))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
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

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if ((ma.gt.0).and.(modulo(ma-lstart,lincr).eq.0)) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = cnst*alpha*fs(ll)
          f = f + aux*plmp1
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = cnst*alpha*fs(ll)
          f = f + aux*plmp2
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            aux = cnst*alpha*fs(ll)
            f = f + aux*plmp3
          endif
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3,sth)
      end subroutine eval_dth_ylm_1D_z

!-----------------------------------------------------------------------
! This subroutine evaluates the imaginary part of the phi derivative of
! a function decomposed over the spherical harmonic basis, divided by 
! sin(theta):
!    f(theta) =  alpha.Im { sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(phi)/sin(theta) }
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_Dphi_ylm_1D_z(alpha,fs,f,cth,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nn, nn_out, lstart, lincr, mm
      double precision, intent(in) :: cth(nn_out)
      double complex, intent(in)   :: alpha,fs(nn)
      double complex, intent(inout):: f(nn_out)
      double precision, allocatable, dimension(:) :: plm1,plm2,plm3,sth
      double precision cnst, coeff, coeff2
      double complex aux
      integer l, ma, ll

      ! easy exit, if possible:
      if (mm.eq.0) return

      ! get workspace:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),sth(nn_out))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = dble(mm)*cnst*sth**(ma-1)
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = cnst*alpha*fs(ll)
          f = f + aux*plm1
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = cnst*alpha*fs(ll)
          f = f + aux*plm2
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            aux = cnst*alpha*fs(ll)
            f = f + aux*plm3
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth)
      end subroutine eval_Dphi_ylm_1D_z

!-----------------------------------------------------------------------
! This subroutine evaluates the theta derivative of a function 
! decomposed over the spherical harmonic basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(theta)
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_ddth_ylm_1D_z(alpha,fs,f,cth,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nn, nn_out, lstart, lincr, mm
      double precision, intent(in) :: cth(nn_out)
      double complex, intent(in)   :: alpha,fs(nn)
      double complex, intent(inout):: f(nn_out)
      double precision, allocatable, dimension(:) :: sth, &
        plm1, plm2, plm3, plmp1, plmp2, plmp3, plmpp1, plmpp2, plmpp3
      double precision cnst, coeff, coeff2
      double complex aux
      integer l, ma, ll

      ! workspace and initialisation:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),plmp1(nn_out), &
               plmp2(nn_out),plmp3(nn_out),plmpp1(nn_out),           &
               plmpp2(nn_out),plmpp3(nn_out),sth(nn_out))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma
      plm2 = plm1*cth*dble(2*ma+1)

      if (ma.eq.0) then
        plmp1 = 0d0
        plmpp1 = 0d0
      else
        plmp1 = cnst*dble(ma)*sth**(ma-1)*cth
        if (ma.eq.1) then
          plmpp1 = -cnst*sth
        else
          plmpp1 =  cnst*ma*((ma-1)*sth**(ma-2)-ma*sth**ma)
        endif
      endif
      plmp2 = (plmp1*cth-plm1*sth)*dble(2*ma+1)
      plmpp2 = (plmpp1*cth-2d0*plmp1*sth-plm1*cth)*dble(2*ma+1)

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if ((ma.gt.0).and.(modulo(ma-lstart,lincr).eq.0)) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = cnst*alpha*fs(ll)
          f = f + aux*plmpp1
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = cnst*alpha*fs(ll)
          f = f + aux*plmpp2
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        plmpp3= coeff2*(cth*plmpp2-2d0*sth*plmp2-cth*plm2) - coeff*plmpp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            aux = cnst*alpha*fs(ll)
            f = f + aux*plmpp3
          endif
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
        plmpp1 = plmpp2
        plmpp2 = plmpp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3,plmpp1,plmpp2,plmpp3,sth)
      end subroutine eval_ddth_ylm_1D_z

!-----------------------------------------------------------------------
! This subroutine projects a function onto the spherical harmonic basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn_out-1)} fs_l Ylm
!
! alpha  = multiplicative constant
! f      = function in real space
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! w      = integration weights associated with the grid theta(:)
! fs     = spectral representation of the function (= OUTPUT)
! nn     = angular resolution in real space
! nn_out = number of spherical harmonics in output
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine project_ylm_1D_z(alpha,f,cth,w,fs,nn,nn_out,lstart,lincr,mm)
 
      implicit none
      integer, intent(in) ::  nn, nn_out, lstart, lincr, mm
      double precision, intent(in) :: cth(nn),w(nn)
      double complex, intent(in)   :: alpha,f(nn)
      double complex, intent(inout):: fs(nn_out)
      double complex, allocatable, dimension(:) :: plm1,plm2,plm3
      double precision, allocatable, dimension(:) :: sth
      double precision cnst, coeff, coeff2
      double complex aux
      integer l, ma, ll, j

      ! initialisation
      allocate(plm1(nn),plm2(nn),plm3(nn),sth(nn))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma*w*f*(2d0*pi)
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = 0d0
          do j=1,nn
            aux = aux + plm1(j)
          enddo
          fs(ll) = fs(ll) + alpha*aux*cnst
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = 0d0
          do j=1,nn
            aux = aux + plm2(j)
          enddo
          fs(ll) = fs(ll) + alpha*aux*cnst
        endif
      endif

      do l = ma+2,lstart+lincr*(nn_out-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.gt.0) then
            aux = 0d0
            do j=1,nn
              aux = aux + plm3(j)
            enddo
            fs(ll) = fs(ll) + alpha*aux*cnst
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth) 
      end subroutine project_ylm_1D_z

!-----------------------------------------------------------------------
! This subroutine projects a function onto the theta derivatives of the
! the spherical harmonics:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(theta)
!
! alpha  = multiplicative constant
! f      = function in real space
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! w      = integration weights associated with the grid theta(:)
! fs     = spectral representation of the function (= OUTPUT)
! nn     = angular resolution in real space
! nn_out = number of spherical harmonics in output
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine project_dth_ylm_1D_z(alpha,f,cth,w,fs,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nn, nn_out, lstart, lincr, mm
      double precision, intent(in) :: cth(nn),w(nn)
      double complex, intent(in)   :: alpha,f(nn)
      double complex, intent(inout):: fs(nn_out)
      double complex, allocatable, dimension(:) :: &
        plm1, plm2, plm3, plmp1, plmp2, plmp3
      double precision, allocatable, dimension(:) :: sth
      double precision cnst, coeff, coeff2
      double complex aux
      integer l, ma, ll, j

      ! get workspace:
      allocate(plm1(nn),plm2(nn),plm3(nn),plmp1(nn), &
               plmp2(nn),plmp3(nn),sth(nn))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo
 
      plm1 = (cnst*2d0*pi)*sth**ma*f*w
      plm2 = plm1*cth*dble(2*ma+1)

      if (ma.eq.0) then
        plmp1 = 0d0
      else
        plmp1 = (cnst*dble(ma)*2d0*pi)*sth**(ma-1)*cth*f*w
      endif
      plmp2 = (plmp1*cth-plm1*sth)*dble(2*ma+1)

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if ((ma.gt.0).and.(modulo(ma-lstart,lincr).eq.0)) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = 0d0
          do j=1,nn
            aux = aux + plmp1(j)
          enddo
          fs(ll) = fs(ll) + alpha*aux*cnst/dble(ma*(ma+1))
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = 0d0
          do j=1,nn
            aux = aux + plmp2(j)
          enddo
          fs(ll) = fs(ll) + alpha*aux*cnst/dble((ma+1)*(ma+2))
        endif
      endif

      do l = ma+2,lstart+lincr*(nn_out-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            aux = 0d0
            do j=1,nn
              aux = aux + plmp3(j)
            enddo
            fs(ll) = fs(ll) + alpha*aux*cnst/dble(l*(l+1))
          endif
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3,sth)
      end subroutine project_dth_ylm_1D_z

!-----------------------------------------------------------------------
! This subroutine projects a function onto the imaginary part of the 
! phi derivatives of the spherical harmonics, divided by sin(theta):
!    f(theta) =  alpha.Im { sum_{l=lstart}^{lstart+2(nn-1)} fs_l [d(Ylm)/d(phi)]/sin(theta) }
!
! alpha  = multiplicative constant
! f      = function in real space
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! w      = integration weights associated with the grid theta(:)
! fs     = spectral representation of the function (= OUTPUT)
! nn     = angular resolution in real space
! nn_out = number of spherical harmonics in output
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine project_Dphi_ylm_1D_z(alpha,f,cth,w,fs,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nn, nn_out, lstart, lincr, mm
      double precision, intent(in) :: cth(nn),w(nn)
      double complex, intent(in)   :: alpha,f(nn)
      double complex, intent(inout):: fs(nn_out)
      double complex, allocatable, dimension(:) :: plm1,plm2,plm3
      double precision, allocatable, dimension(:) :: sth
      double precision cnst, coeff, coeff2
      double complex aux
      integer l, ma, ll, j

      ! easy exit, if possible:
      if (mm.eq.0) return

      ! get workspace:
      allocate(plm1(nn),plm2(nn),plm3(nn),sth(nn))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo
 
      plm1 = dble(mm)*cnst*sth**(ma-1)*f*w*(2d0*pi)
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = (4d0*pi)/dble(2*ma+1)
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(1d0/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = 0d0
          do j=1,nn
            aux = aux + plm1(j)
          enddo
          fs(ll) = fs(ll) + alpha*aux*cnst/dble(ma*(ma+1))
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          aux = 0d0
          do j=1,nn
            aux = aux + plm2(j)
          enddo
          fs(ll) = fs(ll) + alpha*aux*cnst/dble((ma+1)*(ma+2))
        endif
      endif

      do l = ma+2,lstart+lincr*(nn_out-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.gt.0) then
            aux = 0d0
            do j=1,nn
              aux = aux + plm3(j)
            enddo
            fs(ll) = fs(ll) + alpha*aux*cnst/dble(l*(l+1))
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth)
      end subroutine project_Dphi_ylm_1D_z

!-----------------------------------------------------------------------
! This subroutine evaluates a 2D function decomposed over the
! spherical harmonic basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l Ylm
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nnr    = number of grid points in extra dimension
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_ylm_2D_z(alpha,fs,f,cth,nnr,nn,nn_out,lstart,lincr,mm)
 
      implicit none
      integer, intent(in) ::  nnr, nn, nn_out, lstart, lincr, mm
      double precision, intent(in) :: cth(nn_out)
      double complex, intent(in)   :: alpha,fs(nnr,nn)
      double complex, intent(inout):: f(nnr,nn_out)
      double precision, allocatable, dimension(:) :: plm1, plm2, plm3,&
                                                     sth
      double precision cnst, coeff, coeff2
      double complex aux
      integer i, l, ma, ll

      ! initialisation
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),sth(nn_out)) 
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + aux*plm1
          enddo
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + aux*plm2
          enddo
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            do i=1,nnr
              aux = cnst*alpha*fs(i,ll)
              f(i,:) = f(i,:) + aux*plm3
            enddo
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth) 
      end subroutine eval_ylm_2D_z

!-----------------------------------------------------------------------
! This subroutine evaluates the theta derivative of a 2D function 
! decomposed over the spherical harmonic basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(theta)
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nnr    = number of grid points in extra dimension
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_dth_ylm_2D_z(alpha,fs,f,cth,nnr,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nnr, nn, nn_out, lstart, lincr, mm
      double precision, intent(in) :: cth(nn_out)
      double complex, intent(in)   :: alpha,fs(nnr,nn)
      double complex, intent(inout):: f(nnr,nn_out)
      double precision, allocatable, dimension(:) :: &
        plm1, plm2, plm3, plmp1, plmp2, plmp3, sth
      double precision cnst, coeff, coeff2
      double complex aux
      integer i, l, ma, ll

      ! workspace and initialisation:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),plmp1(nn_out), &
               plmp2(nn_out),plmp3(nn_out),sth(nn_out))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
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

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if ((ma.gt.0).and.(modulo(ma-lstart,lincr).eq.0)) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + aux*plmp1
          enddo
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + aux*plmp2
          enddo
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            do i=1,nnr
              aux = cnst*alpha*fs(i,ll)
              f(i,:) = f(i,:) + aux*plmp3
            enddo
          endif
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3,sth)
      end subroutine eval_dth_ylm_2D_z

!-----------------------------------------------------------------------
! This subroutine evaluates the imaginary part of the phi derivative of
! a 2D function decomposed over the spherical harmonic basis, divided by 
! sin(theta):
!    f(theta) =  alpha.Im { sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(phi)/sin(theta) }
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nnr    = number of grid points in extra dimension
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_Dphi_ylm_2D_z(alpha,fs,f,cth,nnr,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nnr, nn, nn_out, lstart, lincr, mm
      double precision, intent(in) :: cth(nn_out)
      double complex, intent(in)   :: alpha,fs(nnr,nn)
      double complex, intent(inout):: f(nnr,nn_out)
      double precision, allocatable, dimension(:) :: plm1,plm2,plm3,sth
      double precision cnst, coeff, coeff2
      double complex aux
      integer i, l, ma, ll

      ! easy exit, if possible:
      if (mm.eq.0) return

      ! get workspace:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),sth(nn_out))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = dble(mm)*cnst*sth**(ma-1)
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + aux*plm1
          enddo
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + aux*plm2
          enddo
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            do i=1,nnr
              aux = cnst*alpha*fs(i,ll)
              f(i,:) = f(i,:) + aux*plm3
            enddo
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth)
      end subroutine eval_Dphi_ylm_2D_z

!-----------------------------------------------------------------------
! This subroutine evaluates the theta derivative of a 2D function 
! decomposed over the spherical harmonic basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(theta)
!
! alpha  = multiplicative constant
! fs     = spectral representation of the function to be evaluated
! f      = function in real space (= OUTPUT)
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! nnr    = number of grid points in extra dimension
! nn     = number of terms in spectral representation
! nn_out = number of points in real space
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine eval_ddth_ylm_2D_z(alpha,fs,f,cth,nnr,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nnr, nn, nn_out, lstart, lincr, mm
      double precision, intent(in) :: cth(nn_out)
      double complex, intent(in)   :: alpha,fs(nnr,nn)
      double complex, intent(inout):: f(nnr,nn_out)
      double precision, allocatable, dimension(:) :: sth, &
        plm1, plm2, plm3, plmp1, plmp2, plmp3, plmpp1, plmpp2, plmpp3
      double precision cnst, coeff, coeff2
      double complex aux
      integer i, l, ma, ll

      ! workspace and initialisation:
      allocate(plm1(nn_out),plm2(nn_out),plm3(nn_out),plmp1(nn_out), &
               plmp2(nn_out),plmp3(nn_out),plmpp1(nn_out),           &
               plmpp2(nn_out),plmpp3(nn_out),sth(nn_out))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma
      plm2 = plm1*cth*dble(2*ma+1)

      if (ma.eq.0) then
        plmp1 = 0d0
        plmpp1 = 0d0
      else
        plmp1 = cnst*dble(ma)*sth**(ma-1)*cth
        if (ma.eq.1) then
          plmpp1 = -cnst*sth
        else
          plmpp1 =  cnst*ma*((ma-1)*sth**(ma-2)-ma*sth**ma)
        endif
      endif
      plmp2 = (plmp1*cth-plm1*sth)*dble(2*ma+1)
      plmpp2 = (plmpp1*cth-2d0*plmp1*sth-plm1*cth)*dble(2*ma+1)

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if ((ma.gt.0).and.(modulo(ma-lstart,lincr).eq.0)) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + aux*plmpp1
          enddo
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = cnst*alpha*fs(i,ll)
            f(i,:) = f(i,:) + aux*plmpp2
          enddo
        endif
      endif

      do l = ma+2,lstart+lincr*(nn-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        plmpp3= coeff2*(cth*plmpp2-2d0*sth*plmp2-cth*plm2) - coeff*plmpp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            do i=1,nnr
              aux = cnst*alpha*fs(i,ll)
              f(i,:) = f(i,:) + aux*plmpp3
            enddo
          endif
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
        plmpp1 = plmpp2
        plmpp2 = plmpp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3,plmpp1,plmpp2,plmpp3,sth)
      end subroutine eval_ddth_ylm_2D_z

!-----------------------------------------------------------------------
! This subroutine projects a 2D function onto the spherical harmonic
! basis:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn_out-1)} fs_l Ylm
!
! alpha  = multiplicative constant
! f      = function in real space
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! w      = integration weights associated with the grid theta(:)
! fs     = spectral representation of the function (= OUTPUT)
! nnr    = number of grid points in extra dimension
! nn     = angular resolution in real space
! nn_out = number of spherical harmonics in output
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine project_ylm_2D_z(alpha,f,cth,w,fs,nnr,nn,nn_out,lstart,lincr,mm)
 
      implicit none
      integer, intent(in) ::  nnr, nn, nn_out, lstart, lincr, mm
      double precision, intent(in) :: cth(nn),w(nn)
      double complex, intent(in)   :: alpha,f(nnr,nn)
      double complex, intent(inout):: fs(nnr,nn_out)
      double precision, allocatable, dimension(:) :: plm1,plm2,plm3,sth
      double precision cnst, coeff, coeff2
      double complex aux
      integer l, ma, ll, i, j

      ! initialisation
      allocate(plm1(nn),plm2(nn),plm3(nn),sth(nn))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo

      plm1 = cnst*sth**ma*w*(2d0*pi)
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = 0d0
            do j=1,nn
              aux = aux + plm1(j)*f(i,j)
            enddo
            fs(i,ll) = fs(i,ll) + alpha*aux*cnst
          enddo
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = 0d0
            do j=1,nn
              aux = aux + plm2(j)*f(i,j)
            enddo
            fs(i,ll) = fs(i,ll) + alpha*aux*cnst
          enddo
        endif
      endif

      do l = ma+2,lstart+lincr*(nn_out-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.gt.0) then
            do i=1,nnr
              aux = 0d0
              do j=1,nn
                aux = aux + plm3(j)*f(i,j)
              enddo
              fs(i,ll) = fs(i,ll) + alpha*aux*cnst
            enddo
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth) 
      end subroutine project_ylm_2D_z

!-----------------------------------------------------------------------
! This subroutine projects a function onto the theta derivatives of the
! the spherical harmonics:
!    f(theta) =  alpha.sum_{l=lstart}^{lstart+2(nn-1)} fs_l d(Ylm)/d(theta)
!
! alpha  = multiplicative constant
! f      = function in real space
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! w      = integration weights associated with the grid theta(:)
! fs     = spectral representation of the function (= OUTPUT)
! nnr    = number of grid points in extra dimension
! nn     = angular resolution in real space
! nn_out = number of spherical harmonics in output
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine project_dth_ylm_2D_z(alpha,f,cth,w,fs,nnr,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nnr, nn, nn_out, lstart, lincr, mm
      double precision, intent(in) :: cth(nn),w(nn)
      double complex, intent(in)   :: alpha,f(nnr,nn)
      double complex, intent(inout):: fs(nnr,nn_out)
      double precision, allocatable, dimension(:) :: &
        plm1, plm2, plm3, plmp1, plmp2, plmp3, sth
      double precision cnst, coeff, coeff2
      double complex aux
      integer l, ma, ll, i, j

      ! get workspace:
      allocate(plm1(nn),plm2(nn),plm3(nn),plmp1(nn), &
               plmp2(nn),plmp3(nn),sth(nn))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo
 
      plm1 = (cnst*2d0*pi)*sth**ma*w
      plm2 = plm1*cth*dble(2*ma+1)

      if (ma.eq.0) then
        plmp1 = 0d0
      else
        plmp1 = (cnst*dble(ma)*2d0*pi)*sth**(ma-1)*cth*w
      endif
      plmp2 = (plmp1*cth-plm1*sth)*dble(2*ma+1)

      cnst = 4d0*pi
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(dble(2*ma+1)/cnst)

      if ((ma.gt.0).and.(modulo(ma-lstart,lincr).eq.0)) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = 0d0
            do j=1,nn
              aux = aux + plmp1(j)*f(i,j)
            enddo
            fs(i,ll) = fs(i,ll) + alpha*aux*cnst/dble(ma*(ma+1))
          enddo
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = 0d0
            do j=1,nn
              aux = aux + plmp2(j)*f(i,j)
            enddo
            fs(i,ll) = fs(i,ll) + alpha*aux*cnst/dble((ma+1)*(ma+2))
          enddo
        endif
      endif

      do l = ma+2,lstart+lincr*(nn_out-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        plmp3 = coeff2*(cth*plmp2-sth*plm2) - coeff*plmp1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.ge.1) then
            do i=1,nnr
              aux = 0d0
              do j=1,nn
                aux = aux + plmp3(j)*f(i,j)
              enddo
              fs(i,ll) = fs(i,ll) + alpha*aux*cnst/dble(l*(l+1))
            enddo
          endif
        endif
        plm1 = plm2
        plm2 = plm3
        plmp1 = plmp2
        plmp2 = plmp3
      enddo

      deallocate(plm1,plm2,plm3,plmp1,plmp2,plmp3,sth)
      end subroutine project_dth_ylm_2D_z

!-----------------------------------------------------------------------
! This subroutine projects a function onto the imaginary part of the 
! phi derivatives of the spherical harmonics, divided by sin(theta):
!    f(theta) =  alpha.Im { sum_{l=lstart}^{lstart+2(nn-1)} fs_l [d(Ylm)/d(phi)]/sin(theta) }
!
! alpha  = multiplicative constant
! f      = function in real space
! cth    = cos(theta(:)) where theta is the grid on which f is defined
! w      = integration weights associated with the grid theta(:)
! fs     = spectral representation of the function (= OUTPUT)
! nnr    = number of grid points in extra dimension
! nn     = angular resolution in real space
! nn_out = number of spherical harmonics in output
! lstart = harmonic degree of first term
! lincr  = increment on harmonic degree
! mm     = azimuthal order
!-----------------------------------------------------------------------
      subroutine project_Dphi_ylm_2D_z(alpha,f,cth,w,fs,nnr,nn,nn_out,lstart,lincr,mm)

      implicit none
      integer, intent(in) ::  nnr, nn, nn_out, lstart, lincr, mm
      double precision, intent(in) :: cth(nn),w(nn)
      double complex, intent(in)   :: alpha,f(nnr,nn)
      double complex, intent(inout):: fs(nnr,nn_out)
      double precision, allocatable, dimension(:) :: plm1,plm2,plm3,sth
      double precision cnst, coeff, coeff2
      double complex aux
      integer l, ma, ll, i, j

      ! easy exit, if possible:
      if (mm.eq.0) return

      ! get workspace:
      allocate(plm1(nn),plm2(nn),plm3(nn),sth(nn))
      sth = sqrt(1d0-cth**2) 

      ! preliminary settings
      ma = abs(mm)

      cnst = 1d0
      do l = 1,2*ma-1,2
        cnst = cnst * dble(l)
      enddo
 
      plm1 = dble(mm)*cnst*sth**(ma-1)*w*(2d0*pi)
      plm2 = plm1*dble(2*ma+1)*cth

      cnst = (4d0*pi)/dble(2*ma+1)
      do l = 1,2*ma
        cnst = cnst*dble(l)
      enddo
      cnst = sqrt(1d0/cnst)

      if (modulo(ma-lstart,lincr).eq.0) then
        ll = 1+(ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = 0d0
            do j=1,nn
              aux = aux + plm1(j)*f(i,j)
            enddo
            fs(i,ll) = fs(i,ll) + alpha*aux*cnst/dble(ma*(ma+1))
          enddo
        endif
      endif

      cnst = cnst * sqrt(dble(2*ma+3))/dble(2*ma+1)
      if (modulo(1+ma-lstart,lincr).eq.0) then
        ll = 1+(1+ma-lstart)/lincr
        if ((ll.ge.1).and.(ll.le.nn)) then
          do i=1,nnr
            aux = 0d0
            do j=1,nn
              aux = aux + plm2(j)*f(i,j)
            enddo
            fs(i,ll) = fs(i,ll) + alpha*aux*cnst/dble((ma+1)*(ma+2))
          enddo
        endif
      endif

      do l = ma+2,lstart+lincr*(nn_out-1)
        coeff = dble(l-1+ma)/dble(l-ma)
        coeff2= dble(2*l-1)/dble(l-ma)
        plm3  = coeff2*cth*plm2 - coeff*plm1
        cnst  = cnst * sqrt(dble((2*l+1)*(l-ma))/dble((2*l-1)*(l+ma)))

        if (modulo(l-lstart,lincr).eq.0) then
          ll = 1+(l-lstart)/lincr
          if (ll.gt.0) then
            do i=1,nnr
              aux = 0d0
              do j=1,nn
                aux = aux + plm3(j)*f(i,j)
              enddo
              fs(i,ll) = fs(i,ll) + alpha*aux*cnst/dble(l*(l+1))
            enddo
          endif
        endif
        plm1 = plm2
        plm2 = plm3
      enddo

      deallocate(plm1,plm2,plm3,sth)
      end subroutine project_Dphi_ylm_2D_z
!-----------------------------------------------------------------------
      end module fast_legendre
