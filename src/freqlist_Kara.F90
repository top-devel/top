#include "config.h"
      module freqlist

      integer, parameter :: max_data=1000
      integer, save, dimension(max_data) :: ll, nn, mm
      double precision, save, dimension(max_data) :: ww, ff
      integer, save :: ndata

!****************************************************************************
! ll = degree l,
! nn = radial order n
! mm = azimuthal order m
! ww = unperturbed, nonrot. freq.
! ff = fully corrected frequency
!
! ndata = number of frequencies
!****************************************************************************

contains

!****************************************************************************
! Program to read in frequencies as calculated in JCD perturbative code,
! for use in non-perturbative code.
!****************************************************************************

      subroutine readfreq()

      use inputs, only: splitfile
      IMPLICIT NONE
      double precision nnaux, llaux, mmaux

      OPEN(5,FILE=trim(splitfile))

      DO ndata=1,max_data

      READ(5,*,end=10) llaux, nnaux, mmaux, ww(ndata), ff(ndata)
      ll(ndata) = floor(llaux+0.5)
      nn(ndata) = floor(nnaux+0.5)
      mm(ndata) = floor(mmaux+0.5)

      END DO

 10   CLOSE(5)

      ndata = ndata - 1

      END subroutine

!****************************************************************************
      end module
