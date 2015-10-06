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

end module python
