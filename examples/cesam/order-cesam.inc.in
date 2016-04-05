!------------------------------------------------------------------------------ 
! This subroutine initialises the arrays ivar and ieq which govern the order
! in which the variables and the equations appear in the matrices.
!------------------------------------------------------------------------------ 
      subroutine init_order()
      implicit none
      integer, parameter :: eq_straight(4)  = (/2,3,5,6/)
      integer, parameter :: eq_shifted(2)   = (/1,4/)
      integer i,j,var,eq,eqn,offset

      print*, "init_order: nr, nt", nr, nt, dm(1)%nvar

      if (allocated(dm(1)%ivar)) deallocate(dm(1)%ivar)
      if (allocated(dm(1)%ieq)) deallocate(dm(1)%ieq)
      allocate(dm(1)%ivar(dm(1)%nvar,nr,nt),dm(1)%ieq(dm(1)%nvar,nr,nt))

      do var=1,6
        do i=1,nr
          do j=1,nt
            dm(1)%ivar(var,i,j) = j + nt*(var-1 + (dm(1)%nvar-2)*(i-1))
          enddo
        enddo
      enddo

      do eqn=1,4
        eq = eq_straight(eqn)
        do i=1,nr
          do j=1,nt
            dm(1)%ieq(eq,i,j) = j + nt*(eq-1 + (dm(1)%nvar-2)*(i-1))
          enddo
        enddo
      enddo

      do eqn=1,2
        eq = eq_shifted(eqn)
        do j=1,nt
          do i=1,nr-1
            dm(1)%ieq(eq,i,j) = j + nt*(eq-1 + (dm(1)%nvar-2)*i)
          enddo
          dm(1)%ieq(eq,nr,j) = j + nt*(eq-1)
        enddo
      enddo

      ! some delicate fine-tuning ...
      offset = 6*nr*nt
      do eq=7,7
        do j=1,nt
          do i=1,nr
            dm(1)%ieq(eq,i,j)  = offset + j + nt*(eq-7 + 2*(nr-i))
            dm(1)%ivar(eq,i,j) = offset + j + nt*(eq-7 + 2*(nr-i))
          enddo
        enddo
      enddo

      do eq=8,8
        do j=1,nt
          do i=1,nr
            dm(1)%ivar(eq,i,j) = offset + j + nt*(eq-7 + 2*(nr-i))
          enddo
          do i=1,nr-1
            dm(1)%ieq(eq,i,j) = offset + j + nt*(eq-7 + 2*(nr-i-1))
          enddo
          dm(1)%ieq(eq,nr,j) = offset + j + nt*(eq-7 + 2*(nr-1))
        enddo
      enddo

      end subroutine
