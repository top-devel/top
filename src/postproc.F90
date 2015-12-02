#include "config.h"

module postproc

      use mod_grid
      integer, allocatable, save :: Ndex(:)
      integer, allocatable, save :: ldom(:)

contains
!-----------------------------------------------------
      subroutine write_output(dir)

          implicit none
          character(len=*), intent(in) :: dir

          print*, "writting to: " // trim(dir)
          call init_index()
          call write_vecp(dir)
          call write_valp(dir)
          call write_grid(dir)
          call find_lmax()
      end subroutine
!-----------------------------------------------------
      subroutine init_index()

          use eigensolve, only: nsol_out, omega
          implicit none

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
      subroutine write_vecp(dir)

          use inputs, only: write_inputs, write_stamp
          use eigensolve
          use mod_grid
          use derivative
          use matrices

          implicit none

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
97        format("('dertype = ', ", I2, "(X, A))")
98        format("(", I2, "(X, I4))")
99        format("(", I3, "(X, f7.5))")
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
100       format(2X, 1pe22.15, 2X, 1pe22.15)
101       format(1pe22.15, 2X, 1pe22.15)
102       format(2X, I3, 2X, I3)

      end subroutine
!-----------------------------------------------------
      subroutine write_grid(dir)

          use mod_grid

          implicit none
          character(len=*), intent(in) :: dir
          integer i, id

          open(unit=2, file=trim(dir)//"grid", status="unknown")
          write(2, *) ndomains
          write(2, 101) (grd(id)%nr, id=1, ndomains)
          do id=1, ndomains
              write(2, 102) (grd(id)%r(i), i=1, grd(id)%nr)
          enddo
          close(2)

101       format(10(X, I4))
102       format(1pe22.15)

      end subroutine
!-----------------------------------------------------
      subroutine write_valp(dir)

          use inputs, only: write_inputs, write_stamp
          use eigensolve

          implicit none
          character(len=*), intent(in) :: dir
          integer i

          open(unit=2, file=trim(dir)//"valp", status="unknown")
          call write_stamp(2)
          call write_inputs(2)
          do i=1, nsol_out
              write(2, 100) omega(Ndex(i))
          enddo
          close(2)

100       format(2X, 1pe22.15, 2X, 1pe22.15)

      end subroutine
!-----------------------------------------------------
      subroutine find_lmax()

          use eigensolve
          use mod_grid
          use matrices

          implicit none

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
          use eigensolve, only: omega
          implicit none

          integer, intent(in) :: i
          double precision, intent(out) :: val

          val = omega(i)
      end subroutine

      subroutine get_valps(vals)
          use eigensolve, only: omega, nsol
          implicit none

          double precision, intent(out) :: vals(nsol)

          vals = omega
      end subroutine

!-----------------------------------------------------
  end module
