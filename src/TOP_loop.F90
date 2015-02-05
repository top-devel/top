#include "config.h"
      program TOP

      use model
      use matrices
      use eigensolve
      use inputs
      use postproc
      use udirectory

      implicit none

      call read_inputs()
      call init_model()

      ! Normally, I shouldn't have to do this, but on
      ! some computers, this could help.
      do m = m_start,m_end,m_incr
        shift_real = shift_start
        do while(shift_real.le.shift_end)
          call init_dir()
          shift_real = shift_real + shift_incr
        enddo
      enddo

      do m = m_start,m_end,m_incr
        call init_a()
        call init_order()
        call init_bc_flag()
        call init_bc_range()

        shift_real = shift_start
        do while(shift_real.le.shift_end)
          ! I still need to put the correct dir_path in
          ! memory:
          call init_dir()
          call run_arncheb(dcmplx(shift_real,shift_imag))
          call write_output()
          print*,shift_real
          shift_real = shift_real + shift_incr
        enddo
      enddo

      end program TOP
