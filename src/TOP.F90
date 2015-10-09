      program TOP

      use model
      use matrices
      use eigensolve
      use inputs
      implicit none

      call read_inputs('/dev/stdin', 10)
      call init_model(dirmodel)
      call init_a()
      call init_order()
      call init_bc_flag()
      call run_arncheb(shift)

      end program TOP
