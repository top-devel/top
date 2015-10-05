module python
contains
    subroutine get_valps(valps)
        use eigensolve, only: omega
        use inputs, only: nsol
        implicit none
        double precision, intent(out) :: valps(nsol)

        valps = omega
    end subroutine

    subroutine get_vecps(vecps)
        use eigensolve, only: vec, a_dim
        use inputs, only: nsol
        implicit none
        double precision, intent(out) :: vecps(a_dim, nsol)

        vecps = vec
    end subroutine

end module python
