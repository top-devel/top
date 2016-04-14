#include "config.h"

module toppy

    use eigensolve, only: run_arncheb, vec, omega, nsol_out
    use mod_grid, only: ndomains, grd, nt
    use abstract_model_mod, only: abstract_model, model_ptr
    use matrices, only: a_dim, init_order, init_a, init_bc_flag, dm
    use inputs, only: read_inputs, lres, init_default
    use postproc, only: write_output, get_sol, get_lvar_size, get_lvar

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
        do id=1, ndomains-1
            n_r = n_r + grd(id)%nr
        enddo

    end subroutine

    subroutine get_nt(n_t)
        integer, intent(out) :: n_t

        n_t = nt

    end subroutine

    subroutine get_lres(l_res)
        integer, intent(out) :: l_res

        l_res = lres

    end subroutine

    subroutine get_dom_nr(id, n_r)
        integer, intent(in) :: id
        integer, intent(out) :: n_r

        n_r = grd(id)%nr

    end subroutine

    subroutine get_nth(n_th)
        integer, intent(out) :: n_th

        n_th = nt

    end subroutine

    subroutine set_nr(n_r)
        integer, intent(in) :: n_r

        integer id

        do id=1, ndomains
            grd(id)%nr = n_r
        end do

    end subroutine
    subroutine set_nt(n_t)
        integer, intent(in) :: n_t

        nt = n_t

    end subroutine

    subroutine get_grid_size(nr, nt)
        integer, intent(out) :: nr, nt

        if (associated(model_ptr)) then
            call model_ptr%get_grid_size(nr, nt)
        else
            print*, "warning undefined model"
            nr = 1
            nt = 1
        endif

    end subroutine

    subroutine get_grid(grid, th, nr, nt)
        integer, intent(in) :: nr, nt
        real(kind=8), intent(out) :: grid(nr, nt), th(nt)

        if (associated(model_ptr)) then
            call model_ptr%get_grid(grid, th, nr, nt)
        else
            print*, "warning undefined model"
            grid = 0
            th = 0
        endif

    end subroutine

    subroutine get_zeta(zeta, nr)
        integer, intent(in) :: nr
        real(kind=8), intent(out) :: zeta(nr)

        integer :: id, skip, npts

        skip = 1
        do id=1, ndomains-1
            npts = grd(id)%nr
            zeta(skip:skip+npts) = grd(id)%r
            skip = skip + npts
        enddo

    end subroutine

    subroutine get_version(v)
        character(len=*), intent(out) :: v

        v = VERSION
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

    subroutine get_var_name(idom, ivar, var_name)
        integer, intent(in) :: idom, ivar
        character(len=64), intent(out) :: var_name

        var_name = dm(idom)%var_name(ivar)

    end subroutine get_var_name

    subroutine get_nvars(idom, nvars)
        integer, intent(in) :: idom
        integer, intent(out) :: nvars

        nvars = dm(idom)%nvar_keep
    end subroutine get_nvars

    subroutine get_ndom(ndom)
        integer, intent(out) :: ndom

        ndom = ndomains
    end subroutine get_ndom

    subroutine pyget_lvar_size(idom, var, lsize)
        integer, intent(in) :: idom
        character(len=*), intent(in) :: var
        integer, intent(out) :: lsize

        call get_lvar_size(idom, var, lsize)
    end subroutine pyget_lvar_size

    subroutine pyget_lvar(idom, var, lm, l)
        integer, intent(in) :: idom, lm
        character(len=*), intent(in) :: var
        integer, intent(out) :: l(lm)

        call get_lvar(idom, var, l)
    end subroutine pyget_lvar

    subroutine init_dati()
        call init_default()
    end subroutine init_dati
end module toppy

module modelpy
    use abstract_model_mod
    use model, only: init_model
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

    subroutine py_init_model(modelfile)
        character(len=*), intent(in) :: modelfile

        call init_model(modelfile)
    end subroutine

end module modelpy
