#include "config.h"

module toppy

    use eigensolve, only: run_arncheb, vec, omega, nsol_out
    use mod_grid, only: ndomains, grd, nt
    use abstract_model_mod, only: abstract_model, model_ptr
    use matrices, only: a_dim, init_order, init_a, init_bc_flag
    use model, only: init_model
    use inputs, only: read_inputs
    use postproc, only: write_output, get_sol

    implicit none

contains

    subroutine dtype(ret)
        character(len=4), intent(out) :: ret

#ifdef USE_COMPLEX
        ret = "cplx"
#else
        ret = "real"
#endif

    end subroutine dtype

    subroutine get_nsol_out(ret)
        integer, intent(out) :: ret

        ret = nsol_out
    end subroutine

    subroutine get_adim(ret)
        integer, intent(out) :: ret

        ret = a_dim
    end subroutine

    subroutine get_valps_cplx(ret, n)
        integer, intent(in) :: n
        complex(kind=8), intent(out) :: ret(n)

        ret = omega

    end subroutine get_valps_cplx

    subroutine get_valps_real(ret, n)
        integer, intent(in) :: n
        real(kind=8), intent(out) :: ret(n)

        ret = omega

    end subroutine get_valps_real

    subroutine get_vecps_cplx(vecps, n, adim)
        integer, intent(in) :: n, adim
        complex(kind=8), intent(out) :: vecps(adim, n)

        vecps = vec
    end subroutine get_vecps_cplx

    subroutine get_vecps_real(vecps, n, adim)
        integer, intent(in) :: n, adim
        real(kind=8), intent(out) :: vecps(adim, n)

        vecps = vec
    end subroutine get_vecps_real

    subroutine get_nr(n_r)
        integer, intent(out) :: n_r

        integer :: id

        n_r = 0
        do id=1, ndomains
            n_r = n_r + grd(id)%nr
        enddo

    end subroutine

    subroutine get_grid(nr, grid)
        integer, intent(in) :: nr
        real(kind=8), intent(out) :: grid(nr)

        integer :: ir, id

        ir = 1
        do id=1, ndomains
            grid(ir:ir+grd(id)%nr-1) = grd(id)%r(:)
            ir = ir + grd(id)%nr
        enddo

    end subroutine

    subroutine get_version(v)
        character(len=*), intent(out) :: v

        v = VERSION
    end subroutine

    subroutine py_init_model(modelfile)
        character(len=*), intent(in) :: modelfile

        call init_model(modelfile)
    end subroutine

    subroutine read_dati(dati)
        character(len=*), intent(in) :: dati

        call read_inputs(dati)
    end subroutine read_dati

    subroutine init_arncheb()

        call init_a()
        call init_order()
        call init_bc_flag()

    end subroutine init_arncheb

    subroutine py_run_arncheb(shift)
#ifdef USE_COMPLEX
        complex(kind=8), intent(in) :: shift
#else
        real(kind=8), intent(in) :: shift
#endif

        call run_arncheb(shift)
    end subroutine py_run_arncheb

    subroutine pywrite_output(dir)
        character(len=*), intent(in) :: dir

        call write_output(dir)
    end subroutine pywrite_output

    subroutine get_solsize(idom, nr, nth)
        integer, intent(in) :: idom
        integer, intent(out) :: nr, nth

        nr = grd(idom)%nr
        nth = nt

    end subroutine get_solsize

    subroutine get_sol_real(idom, isol, var, valp, vecp, nr, nt)
        integer, intent(in) :: idom, isol
        integer, intent(in) :: nr, nt
        character(len=*), intent(in) :: var
        real(kind=8), intent(out) :: valp, vecp(nr, nt)

        call get_sol(idom, isol, var, valp, vecp)
    end subroutine get_sol_real

end module toppy

module modelpy
    use abstract_model_mod
    implicit none

contains
    subroutine get_field_size(fname, n1, n2)
        character(len=*), intent(in) :: fname
        integer, intent(out) :: n1, n2

        real(kind=8), allocatable :: tmp(:, :)

        call model_ptr%get_field(fname, tmp)
        n1 = size(tmp, 1)
        n2 = size(tmp, 2)
        deallocate(tmp)
    end subroutine get_field_size

    subroutine get_field(fname, field, n1, n2)
        character(len=*), intent(in) :: fname
        real(kind=8), intent(out) :: field(n1, n2)
        integer, intent(in) :: n1, n2

        real(kind=8), allocatable :: tmp(:, :)

        call model_ptr%get_field(fname, tmp)

        field = tmp
        deallocate(tmp)

    end subroutine get_field

end module modelpy

