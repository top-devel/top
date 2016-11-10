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
          type(GRID), allocatable, save, target :: grd(:)

#ifndef USE_MULTI
          integer, pointer :: nr
#endif
#ifdef USE_1D
          double precision, pointer :: r(:)
#endif

          integer(kind=c_int), save, bind(c) :: ndomains

          integer, save :: nt

contains

          subroutine init_grid(ndom)
              implicit none
              integer, intent(in) :: ndom
              integer :: nd

              ndomains = ndom
              if (allocated(grd)) deallocate(grd)
              allocate(grd(ndom))

#ifndef USE_MULTI
              nr => grd(1)%nr
#endif
#ifdef USE_1D
              nt = 1
#else
              do nd = 1, ndom
                  grd(nd)%nr = 0
                  grd(nd)%order = 2
                  grd(nd)%dertype = 'CHEB'
                  grd(nd)%mattype = 'FULL'
              end do
#endif

          end subroutine init_grid

          subroutine init_radial_grid()
              implicit none
              integer id

              do id=1, ndomains
                  if (allocated(grd(id)%r)) deallocate(grd(id)%r)
                  allocate(grd(id)%r(1:grd(id)%nr))
              end do
#ifdef USE_1D
              nt = 1
              r => grd(1)%r
#endif

          end subroutine

      end module
