module python
contains
    subroutine get_valps(valps)
        use eigensolve, only: omega, nsol_out
        use inputs, only: nsol
        implicit none
        double precision, intent(out) :: valps(nsol_out)

        valps = omega
    end subroutine

    subroutine get_vecps(vecps)
        use eigensolve, only: vec, a_dim, nsol_out
        use inputs, only: nsol
        implicit none
        double precision, intent(out) :: vecps(a_dim, nsol_out)

        vecps = vec
    end subroutine

    subroutine get_nr(nr)
        use mod_grid, only: ndomains, grd
        implicit none
        integer, intent(out) :: nr

        integer :: id

        nr = 0
        do id=1, ndomains
            nr = nr + grd(id)%nr
        enddo

    end subroutine

    subroutine get_grid(nr, grid)
        use mod_grid, only: ndomains, grd
        implicit none
        integer, intent(in) :: nr
        double precision, intent(out) :: grid(nr)

        integer :: ir, id

        ir = 1
        do id=1, ndomains
            print*, "dom: ", id, "size: ", grd(id)%nr
            print*, "write: ", ir, ir+grd(id)%nr-1
            grid(ir:ir+grd(id)%nr-1) = grd(id)%r(:)
            ir = ir + grd(id)%nr
        enddo

    end subroutine

end module python
