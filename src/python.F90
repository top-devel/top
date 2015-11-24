#include "config.h"
module toppy

    use eigensolve
    use inputs
    use model
    use matrices, only: a_dim

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
        integer :: i
        complex(kind=8), intent(out) :: ret(n)

        ret = omega

    end subroutine get_valps_cplx

    subroutine get_valps_real(ret, n)
        integer, intent(in) :: n
        integer :: i
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
        use mod_grid, only: ndomains, grd
        implicit none
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

    subroutine py_run_arncheb(shift_in)
#ifdef USE_COMPLEX
        complex(kind=8), intent(in) :: shift_in
#else
        real(kind=8), intent(in) :: shift_in
#endif

        call run_arncheb(shift)
    end subroutine py_run_arncheb

end module toppy
