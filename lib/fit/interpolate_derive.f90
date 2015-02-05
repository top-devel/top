        subroutine interpolate_derive(xin,yin,nin,xout,yout,nout)

        implicit none

        ! input and output parameters
        integer, intent(in) :: nin, nout
        double precision, intent(in), dimension(nin)  :: xin, yin
        double precision, intent(in), dimension(nout) :: xout
        double precision, intent(out), dimension(nout):: yout

        ! internal variables
        integer, parameter :: iopt = 0, k = 3
        double precision, parameter :: s = 0d0
        !integer, parameter :: iopt = 0, k = 5
        !double precision, parameter :: s = 1d-6
        integer n, ier, nest, lwrk
        double precision xlower, xupper, fp, d(k+1)
        double precision, allocatable, dimension(:)  :: c, t, wrk, w
        integer, allocatable, dimension(:)  :: iwrk
        integer i

        ! a few preliminary settings, and allocation of space:
        nest = nin + k + 1
        lwrk = nin*(k+1)+nest*(7+3*k)
        allocate(c(nest),t(nest),wrk(lwrk),iwrk(nest),w(nin))

        do i=1,nin
          w(i) = 1d0
        enddo

        n = nest
        xlower = min(xin(1),xout(1))
        xupper = max(xin(nin),xout(nout))

        !------------------------------------------------------------------
        ! nin       = number of data points
        ! xin(nin)  = array containing the values of the independant 
        !             variable, xin
        ! yin(nin)  = array containing the values of the dependant
        !             variable, yin
        ! w(nin)    = weights for the different points (if smoothing)
        ! xlower,   = bounds of the interval over which the smoothing/
        ! xupper      interpolation is calculated
        ! k         = order of the b-splines (k = 3 corresponds to cubic - a
        !             different convention is used in CESAM)
        ! s         = smoothing factor (s=0 corresponds to interpolation)
        ! nest      = maximum number of knots (nest=nin+k+1 for s=0)
        ! n         = will contain actual number of knots (but must be
        !             specified on entry
        ! t(nest)   = will contain the positions of the knots 
        ! c(nest)   = will contain the coefficients of the b-spline 
        !             representation
        ! fp        = will contain the weighted sum of squared residuals
        ! wrk(lwrk) = real array used as work space
        ! lwrk      = dimension of work space (>= nin*(k+1)+nest*(7+3*k)) 
        ! iwrk(nest)= integer array used as work space
        ! ier       = error flag (<= 0 on successful exit)
        !------------------------------------------------------------------
      
        ! This finds the coefficients to the b-spline representation  

        call curfit(iopt,nin,xin,yin,w,xlower,xupper,k,s,nest,n,t,c,fp, &
                    wrk,lwrk,iwrk,ier)

        if (ier.gt.0) then
          print*,"Failure in interpolate_derive.  ier = ",ier
          stop
        endif

        !------------------------------------------------------------------
        ! xout(nout)= grid on which to calculate (interpolate) function
        ! yout(nout)= values of interpolated function on grid xout
        !------------------------------------------------------------------

        ! This does the calculates the values of the function on the
        ! new grid based on the coefficients given by the previous
        ! subroutine.

        do i=1,nout
          call spalde(t,n,c,k+1,xout(i),d,ier)
          yout(i) = d(2)
        enddo

        deallocate(c,t,wrk,iwrk,w)
        return
        end subroutine
