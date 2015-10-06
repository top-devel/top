!------------------------------------------------------------------------------ 
! This subroutine initialises the arrays ivar and ieq which govern the order
! in which the variables and the equations appear in the matrices.
! The following is a generic form:
!------------------------------------------------------------------------------ 

      subroutine init_order()
      implicit none
      integer i,j,id,var,eq

      do id=1,ndomains
        allocate(dm(id)%ivar(dm(id)%nvar,grd(id)%nr,nt))
        allocate(dm(id)%ieq(dm(id)%nvar,grd(id)%nr,nt))

        do var=1,dm(id)%nvar
          do i=1,grd(id)%nr
            do j=1,nt
              dm(id)%ivar(var,i,j) = j + nt*(var-1 + dm(id)%nvar*(i-1))
            enddo
          enddo
        enddo

        do eq=1,dm(id)%nvar
          do i=1,grd(id)%nr
            do j=1,nt
              dm(id)%ieq(eq,i,j) = j + nt*(eq-1 + dm(id)%nvar*(i-1))
            enddo
          enddo
        enddo
      enddo

      end subroutine init_order
