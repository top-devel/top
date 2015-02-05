      subroutine nwton(xin,yin,nin,xout,yout,nout)

      implicit none

      ! input and output parameters
      integer, intent(in) :: nin, nout
      double precision, intent(in), dimension(nin)  :: xin, yin
      double precision, intent(in), dimension(nout) :: xout
      double precision, intent(out), dimension(nout):: yout

      ! internal variables:
      double precision, allocatable, dimension(:) :: tab
      integer i, j

      allocate(tab(nin))
      do i=1,nin
        tab(i) = yin(i)
      enddo

      ! find polynomial coefficients in Newton form:
      do j=2,nin
        do i=nin,j,-1
          tab(i) = (tab(i)-tab(i-1))/(xin(i)-xin(i-j+1))
        enddo
      enddo

      ! apply Horner's method to evaluate the polynomial on the new grid
      do i=1,nout
        yout(i) = tab(nin)
      enddo

      do i=nin-1,1,-1
        do j=1,nout
          yout(j) = yout(j)*(xout(j)-xin(i)) + tab(i)
        enddo
      enddo

      deallocate(tab)
      return
      end subroutine
