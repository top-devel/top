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
      character*(6) sparity

      call read_inputs()
      call init_model()
      call init_dir()
      call init_a()
      call init_order()
      call init_bc_flag()
      call init_bc_range()
      call run_arncheb(dcmplx(shift_real,shift_imag))
      call write_output()

      end program TOP
