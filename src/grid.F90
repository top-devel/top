#include "config.h"
      module mod_grid
      use iso_c_binding

          type GRID
              character(len=4) :: mattype
              character(len=4) :: dertype
              integer(kind=c_int) :: order
              integer(kind=c_int) :: nr
              double precision, allocatable :: r(:)
          end type GRID

          integer(kind=c_int), save, bind(c) :: ndomains
          type(GRID), allocatable, save, target :: grd(:)

          integer, save :: nt

contains

          subroutine init_grid(ndom)
              implicit none
              integer, intent(in) :: ndom

              ndomains = ndom
              if (allocated(grd)) deallocate(grd)
              allocate(grd(ndom))
          end subroutine

          subroutine init_radial_grid()
              implicit none
              integer id

              do id=1, ndomains
                  if (allocated(grd(id)%r)) deallocate(grd(id)%r)
                  allocate(grd(id)%r(1:grd(id)%nr))
              end do

          end subroutine

      end module
