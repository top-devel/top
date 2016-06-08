module legpy
    use fast_legendre, only: eval_ylm, eval_dth_ylm
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

    subroutine eval2d_z(fs, cost, l, linc, m, f, nr, ns, nt)

        complex(kind=8), intent(in) :: fs(:, :)
        real(kind=8), intent(in) :: cost(:)
        integer, intent(in) :: nr, ns, nt, m, l, linc
        complex(kind=8), intent(out) :: f(nr, nt)
!f2py   integer intent(hide), depend(fs) :: nr=shape(fs, 0), ns=shape(fs, 1)
!f2py   integer intent(hide), depend(cost) :: nt=shape(cost, 0)

        f = dcmplx(0.d0, 0.d0)
        call eval_ylm(dcmplx(1.d0, 0.d0), fs, f, cost, nr, ns, nt, l, linc, m)
    end subroutine eval2d_z

    subroutine eval2d_dth_z(fs, cost, l, linc, m, f, nr, ns, nt)

        complex(kind=8), intent(in) :: fs(:, :)
        real(kind=8), intent(in) :: cost(:)
        integer, intent(in) :: nr, ns, nt, m, l, linc
        complex(kind=8), intent(out) :: f(nr, nt)
!f2py   integer intent(hide), depend(fs) :: nr=shape(fs, 0), ns=shape(fs, 1)
!f2py   integer intent(hide), depend(cost) :: nt=shape(cost, 0)

        f = dcmplx(0.d0, 0.d0)
        call eval_dth_ylm(dcmplx(1.d0, 1.d0), fs, f, cost, nr, ns, nt, l, linc, m)
    end subroutine eval2d_dth_z

end module legpy
