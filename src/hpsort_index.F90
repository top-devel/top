#include "config.h"
!------------------------------------------------------------------------------
! This program sorts an array based on the heapsort method.  It is copied from:
!
! "Numerical Recipes in Fortran: The Art of Scientific Computing"
! Second Edition.
! William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery
! Cambridge University Press, 1992.
!------------------------------------------------------------------------------
! This routine sorts creates an index array "index(.)" such that
! the array "ra(index(.))" is in ascending numerical order.
!------------------------------------------------------------------------------

      subroutine hpsort_index(n,ra,index)

          implicit none
          integer n
          double precision ra(1:n)
          integer index(1:n)
          integer i, ir, j, l
          integer iindex

          ! Initialise the index array
          do i=1,n
              index(i) = i
          enddo

          if (n.lt.2) return

          ! The index l will be decremented from its initial value down to 1
          ! during the "hiring" (heap creation) phase.  Once it reaches 1, the
          ! index ir will be decremented from its initial value down to 1 during
          ! the "retirement-and-promotion" (heap selection) phase.

          l = n/2+1
          ir = n
          do
              if (l.gt.1) then                ! Still in hiring phase
                  l=l-1
                  iindex = index(l)
              else                            ! In retiremnet-and-promotion phase.
                  iindex = index(ir)            ! Clear a space at end of array.
                  index(ir) = index(1)          ! Retire the top of the heap into it
                  ir=ir-1                       ! Decrease the size of the corporation.
                  if(ir.eq.1)then               ! Done with the last promotion.
                      index(1) = iindex           ! The least competent worker of all!
                      return
                  endif
              endif
              i=l                             ! Whether in the hiring phase or
              ! promotion phase, we here set up to
              ! sift down element rra to its proper
              ! level.
              j=l+l
              do while (j.le.ir)
                  if (j.lt.ir) then
                      if (ra(index(j)).lt.ra(index(j+1))) j=j+1
                      ! Compare to the better underling
                  endif
                  if (ra(iindex).lt.ra(index(j))) then
                      index(i)=index(j)           ! Demote ra(iindex)
                      i=j
                      j=j+j
                  else                          ! This is rra's level.  Set j to
                      ! terminate the sift-down.
                      j=ir+1
                  endif
              enddo
              index(i) = iindex               ! Put rra into its slot.
          enddo
      end subroutine
