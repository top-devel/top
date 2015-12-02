#include "config.h"

module postproc

    use mod_grid
    use eigensolve, only: nsol_out, omega, vec
    use matrices
    use inputs, only: nsol, write_stamp, write_inputs

    implicit none

    integer, allocatable, save :: Ndex(:)
    integer, allocatable, save :: ldom(:)

contains
    !-----------------------------------------------------
    subroutine write_output(dir)

        character(len=*), intent(in) :: dir

        call init_index()
        call write_vecp(dir)
        call write_valp(dir)
        call write_grid(dir)
        call find_lmax()
    end subroutine
    !-----------------------------------------------------
    subroutine init_index()

        if (allocated(Ndex)) deallocate(Ndex)
        allocate(Ndex(nsol_out))
#ifdef USE_COMPLEX
        call hpsort_index(nsol_out, dreal(omega), Ndex)
#else
        call hpsort_index(nsol_out, omega, Ndex)
#endif
        return
    end subroutine
    !-----------------------------------------------------
    subroutine get_sol(idom, nsol, var, valp, vecp)

        integer, intent(in) :: idom, nsol
        character(len=*), intent(in) :: var
#ifdef USE_COMPLEX
        complex(kind=8), intent(out) :: vecp(:, :), valp
        complex(kind=8) :: vec_value
#else
        real(kind=8), intent(out) :: vecp(:, :), valp
        real(kind=8) :: vec_value
#endif
        integer :: iv, is, ivar, j, i, ii, der_id, lbder, ubder
        integer :: isol

        if (idom > ndomains) then
            print'(A, I2, A)', "error only have ", ndomains, " domains"
            return
        endif

        if (nsol > nsol_out) then
            print'(A, I2, A)', "error only have ", nsol_out, " solutions"
            return
        endif

        call init_index()
        isol = Ndex(nsol)

        ivar = 0
        do iv=1, dm(idom)%nvar_keep
            if (trim(var) == trim(dm(idom)%var_name(iv))) then
                ivar = iv
            endif
        enddo
        if (ivar == 0) then
            print*, "no variable named "//trim(var)
            return
        endif

        do j=1, nt
            do i=1, grd(idom)%nr
#ifdef USE_COMPLEX
                vec_value = (0d0, 0d0)
#else
                vec_value = 0d0
#endif
                der_id = dmat(idom)%der_id
                lbder = dmat(idom)%lbder(der_id)
                ubder = dmat(idom)%ubder(der_id)
                do ii=max(1, i-lbder), min(grd(idom)%nr, i+ubder)
                    vec_value = vec_value &
                        + dmat(idom)%derive(i, ii, der_id) &
                        * vec(dm(idom)%offset+dm(idom)%ivar(ivar, ii, j), isol)
                enddo
                vecp(i, j) = vec_value
            enddo
        enddo

        valp = omega(isol)

    end subroutine get_sol
    !-----------------------------------------------------
    subroutine write_vecp(dir)

        ! use inputs, only: write_inputs, write_stamp
        ! use eigensolve
        ! use mod_grid
        ! use derivative
        ! use matrices

        character(len=*), intent(in) :: dir

        integer isol, iisol, id, i, ii, j, var, vvar, der_id, lbder, ubder
        double complex vec_value
        character*(64) str

        open(unit=2, file=trim(dir)//"vecp", status="unknown")
        call write_stamp(2)
        write(str, 97) ndomains
        write(2, str) (trim(grd(id)%dertype), id=1, ndomains)
        write(2, *) ndomains
        write(str, 98) ndomains
        write(2, str) (grd(id)%nr, id=1, ndomains)
        write(str, 99) 2*ndomains
        write(2, str) (grd(id)%r(1), grd(id)%r(grd(id)%nr), id=1, ndomains)
97      format("('dertype = ', ", I2, "(X, A))")
98      format("(", I2, "(X, I4))")
99      format("(", I3, "(X, f7.5))")
        call write_inputs(2)
        write(str, 98) ndomains+1
        write(2, str) nsol_out, (dm(id)%nvar_keep , id=1, ndomains)
        do iisol=1, nsol_out
            isol = Ndex(iisol)
            write(2, 100) omega(isol)
            do id=1, ndomains
                der_id = dmat(id)%der_id
                lbder = dmat(id)%lbder(der_id)
                ubder = dmat(id)%ubder(der_id)
                do vvar=1, dm(id)%nvar_keep
                    var = dm(id)%var_list(vvar)
                    write(2, '(a, 2X, I2)') trim(dm(id)%var_name(var)), id
                    do j=1, nt
                        write(2, 102) j, dm(id)%lvar(j, var)
                        do i=1, grd(id)%nr
                            vec_value = (0d0, 0d0)
                            do ii=max(1, i-lbder), min(grd(id)%nr, i+ubder)
                                vec_value = vec_value &
                                    + dmat(id)%derive(i, ii, der_id) &
                                    * vec(dm(id)%offset+dm(id)%ivar(var, ii, j), isol)
                            enddo
                            write(2, 101) vec_value
                        enddo
                    enddo
                enddo
            enddo
        enddo
        close(2)
        call system("bzip2 --force "//trim(dir)//"vecp")
100     format(2X, 1pe22.15, 2X, 1pe22.15)
101     format(1pe22.15, 2X, 1pe22.15)
102     format(2X, I3, 2X, I3)

    end subroutine
    !-----------------------------------------------------
    subroutine write_grid(dir)

        character(len=*), intent(in) :: dir
        integer i, id

        open(unit=2, file=trim(dir)//"grid", status="unknown")
        write(2, *) ndomains
        write(2, 101) (grd(id)%nr, id=1, ndomains)
        do id=1, ndomains
            write(2, 102) (grd(id)%r(i), i=1, grd(id)%nr)
        enddo
        close(2)

101     format(10(X, I4))
102     format(1pe22.15)

    end subroutine
    !-----------------------------------------------------
    subroutine write_valp(dir)

        character(len=*), intent(in) :: dir
        integer i

        open(unit=2, file=trim(dir)//"valp", status="unknown")
        call write_stamp(2)
        call write_inputs(2)
        do i=1, nsol_out
            write(2, 100) omega(Ndex(i))
        enddo
        close(2)

100     format(2X, 1pe22.15, 2X, 1pe22.15)

    end subroutine
    !-----------------------------------------------------
    subroutine find_lmax()

        integer isol, iisol, id, i, j, var, lmax
        integer, allocatable :: iu(:)
        double precision my_max, my_sum
        character*(2) :: str

        if (allocated(ldom)) deallocate(ldom)
        allocate(ldom(nsol_out), iu(ndomains-1))

        do id=1, ndomains-1
        write(str, '(I2)') id
        iu(id) = 0
        do var = 1, dm(id)%nvar
        if (  (trim(dm(id)%var_name(var)).eq.'u' //trim(adjustl(str)))   &
            .or.(trim(dm(id)%var_name(var)).eq.'Er'//trim(adjustl(str)))) then
            iu(id) = var
            exit
        endif
        enddo
        if (iu(id).eq.0) stop "Couldn't find variable u"
        enddo

        do iisol=1, nsol_out
        isol = Ndex(iisol)
        lmax = -1
        my_max = 0d0
        do j=1, nt
        my_sum = 0d0
        do id=1, ndomains-1
        do i=1, grd(id)%nr
        my_sum = my_sum + &
            abs(vec(dm(id)%offset+dm(id)%ivar(iu(id), i, j), isol))**2
        enddo
        enddo
        if (my_sum.gt.my_max) then
            my_max = my_sum
            lmax = dm(1)%lvar(j, iu(1))
        endif
        enddo
        ldom(iisol) = lmax
        enddo

        deallocate(iu)

    end subroutine

    subroutine get_valp(i, val)

        integer, intent(in) :: i
        double precision, intent(out) :: val

        call init_index()
        val = omega(Ndex(i))
    end subroutine

    subroutine get_valps(vals)

        double precision, intent(out) :: vals(nsol)

        call init_index()
        vals = omega
    end subroutine

    !-----------------------------------------------------
end module
