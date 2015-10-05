#include "config.h"
      program TOP

          use model
          use mod_grid
          use matrices
          use eigensolve
          use inputs
          use postproc
          use udirectory

          implicit none
          character*(128) string
          character*(6) sparity
          character(len=255) arg

          call get_command_argument(1, arg)
          if (len(trim(arg)) == 0) then
              call read_inputs('/dev/stdin') ! yeah, this one is dirty
          else
              call read_inputs(trim(arg))
          endif
          call init_model()
          call init_dir()
          call init_a()
          call init_order()
          call init_bc_flag()
#ifdef USE_COMPLEX
          call run_arncheb(dcmplx(shift_real,shift_imag))
#else
          call run_arncheb(shift)
#endif
          call write_output()

      end program TOP
