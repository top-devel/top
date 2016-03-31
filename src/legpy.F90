module legpy
    use fast_legendre	, only: eval_ylm, eval_dth_ylm
    implicit none

contains

    subroutine eval2d(fs, cost, l, linc, m, f, nr, ns, nt)

        real(kind=8), intent(in) :: fs(:, :)
        real(kind=8), intent(in) :: cost(:)
        integer, intent(in) :: nr, ns, nt, m, l, linc
        real(kind=8), intent(out) :: f(nr, nt)
!f2py   integer intent(hide), depend(fs) :: nr=shape(fs, 0), ns=shape(fs, 1)
!f2py   integer intent(hide), depend(cost) :: nt=shape(cost, 0)

        f = 0.0
        call eval_ylm(1d0, fs, f, cost, nr, ns, nt, l, linc, m)
    end subroutine eval2d

    subroutine eval2d_dth(fs, cost, l, linc, m, f, nr, ns, nt)

        real(kind=8), intent(in) :: fs(:, :)
        real(kind=8), intent(in) :: cost(:)
        integer, intent(in) :: nr, ns, nt, m, l, linc
        real(kind=8), intent(out) :: f(nr, nt)
!f2py   integer intent(hide), depend(fs) :: nr=shape(fs, 0), ns=shape(fs, 1)
!f2py   integer intent(hide), depend(cost) :: nt=shape(cost, 0)

        f = 0.0
        call eval_dth_ylm(1d0, fs, f, cost, nr, ns, nt, l, linc, m)
    end subroutine eval2d_dth

    ! subroutine grid_to_spec()
    ! end subroutine grid_to_spec

end module legpy
