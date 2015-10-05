#include "config.h"
      program TOP

      use model
      use matrices
      use eigensolve
      use inputs
      use postproc
      use udirectory

      implicit none
      character*(128) string
      character*(10) sparity
      integer i, diff

      call read_inputs()
      call init_model()

      if (iparity.eq.0) then
        sparity = "even modes"
      else
        sparity = "odd modes"
      endif

      write(string,101) mass, eta, alpha
      call system(trim(string))
      do m = m_start,m_end,m_incr
        write(string,102) m, sparity
        call system(trim(string))
        shift_real = 2d0*shift2 - shift1
        call init_a()
        call init_order()
        call init_bc_flag()
        call init_bc_range()

        call init_dir()
        call run_arncheb(dcmplx(shift_real,shift_imag))
        call write_output()
        print*,shift_real
        shift1 = shift2
        shift2 = dreal(omega(1))
        diff = abs(ldom(1)-abs(m)+1)
        do i=2,nsol_out
          if (abs(ldom(i)-abs(m)).lt.diff) then
            shift2 = dreal(omega(i))
            diff = abs(ldom(i)-abs(m)+1)
          endif
        enddo
      enddo

101   format("echo '# mass=",0pf5.2,", eta=",0pf4.2, &
             ", alpha=",0pf4.2,"' >> valps")
102   format("echo '# m=",I2,", ",a10,"' >> valps")
      end program TOP
