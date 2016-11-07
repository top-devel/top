!------------------------------------------------------------------------------ 
! This subroutine initialises the arrays ivar and ieq which govern the order
! in which the variables and the equations appear in the matrices.
!------------------------------------------------------------------------------ 

      subroutine init_order()
      implicit none
      integer i,var,eq

      if (allocated(dm(1)%ivar)) deallocate(dm(1)%ivar)
      if (allocated(dm(1)%ieq)) deallocate(dm(1)%ieq)
      allocate(dm(1)%ivar(dm(1)%nvar,nr),dm(1)%ieq(dm(1)%nvar,nr))

      do var=1,dm(1)%nvar
        do i=1,nr
          dm(1)%ivar(var,i) = var + (i-1)*dm(1)%nvar
          !dm(1)%ivar(var,i) = i + (var-1)*nr
        enddo
      enddo

      do eq = 1,dm(1)%nvar
        do i=1,nr
          dm(1)%ieq(eq,i) = eq + (i-1)*dm(1)%nvar
        enddo
      enddo

      end subroutine
