!------------------------------------------------------------------------------ 
! This subroutine initialises the arrays ivar and ieq which govern the order
! in which the variables and the equations appear in the matrices.
!------------------------------------------------------------------------------ 

      subroutine init_order()
      implicit none
      integer i,var,eq,eqn
      integer, parameter :: eq_straight(3) = (/2, 3, 5/)
      integer, parameter :: eq_shifted(2) = (/1, 4/)

      if (allocated(dm(1)%ivar)) deallocate(dm(1)%ivar)
      if (allocated(dm(1)%ieq)) deallocate(dm(1)%ieq)
      allocate(dm(1)%ivar(dm(1)%nvar,nr),dm(1)%ieq(dm(1)%nvar,nr))

      do var=1,dm(1)%nvar
        do i=1,nr
          dm(1)%ivar(var,i) = var + (i-1)*dm(1)%nvar
          !dm(1)%ivar(var,i) = i + (var-1)*nr
        enddo
      enddo

      do eqn = 1, 3
        eq = eq_straight(eqn)
        do i=1,nr
          dm(1)%ieq(eq,i) = eq + (i-1)*dm(1)%nvar
        enddo
      enddo

      do eqn = 1, 2
        eq = eq_shifted(eqn)
        do i=1,nr-1
          dm(1)%ieq(eq,i) = eq + i*dm(1)%nvar
        enddo
        dm(1)%ieq(eq,nr) = eq
      enddo

      end subroutine
