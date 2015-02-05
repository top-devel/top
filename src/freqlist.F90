#include "config.h"
      module freqlist

      integer, parameter :: max_data=1000
      integer, save, dimension(max_data) :: ll, mm
      double precision, save, dimension(max_data) :: fr, fi
      integer, save :: ndata 

!****************************************************************************
! ll = parity
! mm = azimuthal order m
! ff = frequency
!
! ndata = number of frequencies
!****************************************************************************

contains

!****************************************************************************
! Program to read in frequencies as calculated in JCD perturbative code,
! for use in non-perturbative code.
!****************************************************************************

      subroutine readfreq() 

      use inputs, only: freqlist
      IMPLICIT NONE

      OPEN(5,FILE=trim(freqlist))

      DO ndata=1,max_data
        READ(5,*,end=10) ll(ndata), mm(ndata), fr(ndata), fi(ndata)
      END DO

 10   CLOSE(5)

      ndata = ndata - 1

      END subroutine

!****************************************************************************
      end module
