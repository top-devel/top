        program testfit

        implicit none
        integer, parameter :: nin = 200, nout = 1000
        double precision, dimension(nin)  :: xin, yin
        double precision, dimension(nout) :: xout, yout
        double precision, parameter :: pi = 3.141592653589793d0
        integer i

        do i=1,nin
          xin(i) = dble(i-1)/dble(nin-1)
          yin(i) = dcos(5d0*xin(i)**2)
        enddo

        do i=1,nout
          xout(i) = dsin(dble(i-1)*pi/dble(nout*2-2))
        enddo

        call interpolate_derive(xin,yin,nin,xout,yout,nout)

        open(unit=2,file="delme",status="unknown")
        do i=1,nin
          write(2,*) xin(i), yin(i)
        enddo
        write(2,*) "&"
        do i=1,nout
          write(2,*) xout(i), yout(i)
        enddo
        write(2,*) "&"
        close(2)

        end program testfit
