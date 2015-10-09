#include "config.h"
module python
    use iso_c_binding
    use string

contains

    subroutine get_valps(valps) bind(c)
        use eigensolve, only: omega, nsol_out
        use inputs, only: nsol
        implicit none
        real(kind=c_double), intent(out) :: valps(nsol_out)

        valps = omega
    end subroutine

    subroutine get_vecps(vecps) bind(c)
        use eigensolve, only: vec, a_dim, nsol_out
        use inputs, only: nsol
        implicit none
        real(kind=c_double), intent(out) :: vecps(a_dim, nsol_out)

        vecps = vec
    end subroutine

    subroutine get_nr(nr) bind(c)
        use mod_grid, only: ndomains, grd
        implicit none
        integer(kind=c_int), intent(out) :: nr

        integer :: id

        nr = 0
        do id=1, ndomains
            nr = nr + grd(id)%nr
        enddo

    end subroutine

    subroutine get_grid(nr, grid) bind(c)
        use mod_grid, only: ndomains, grd
        implicit none
        integer(kind=c_int), intent(in) :: nr
        real(kind=c_double), intent(out) :: grid(nr)

        integer :: ir, id

        ir = 1
        do id=1, ndomains
            grid(ir:ir+grd(id)%nr-1) = grd(id)%r(:)
            ir = ir + grd(id)%nr
        enddo

    end subroutine

    subroutine get_version(v, n) bind(c)
        implicit none
        integer(kind=c_int), intent(in) :: n
        character(kind=c_char), intent(out) :: v(n)

        v = transfer(VERSION, " ", size=len_trim(VERSION))
    end subroutine

    subroutine python_init_model(filename, n) bind(c)
        use model
        implicit none
        integer(kind=c_int), intent(in) :: n
        character(kind=c_char), intent(in) :: filename(n)

        print*, "model: ", cpy_str(filename)
        call init_model(cpy_str(filename))
    end subroutine

end module python
