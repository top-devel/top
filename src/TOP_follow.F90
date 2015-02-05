#include "config.h"
      program TOP

      use model
      use matrices
      use eigensolve
      use inputs
      use postproc
      use udirectory
      use script

      implicit none
      character*(128) string
      character*(10) sparity
      integer ijob, i, imass, diff, ltarget, j, lbest
      double precision freq1, freq2, freq3
      double precision mss1, mss2, mss3

      call read_inputs()  ! This reads the dati file.
      call read_script()  ! This reads a script file which says which
                          ! frequencies to follow as a function of mass.

      if (index(stamp,"even modes").ne.0) then
        sparity = "even modes"
      else
        sparity = "odd modes"
      endif

      do ijob=1,njobs

        do i=1,jl(ijob)%list_size

          ! Initialise the two first frequencies:
          call set_filename(jl(ijob)%mass_list(1),jl(ijob)%eta, &
                            jl(ijob)%alpha)
          call read_model()

          ! Normalise frequency by Keplerian break-up rotation rate:
          freq1 = jl(ijob)%w1(i)*sqrt(p1D(1)*radius2/(rho1D(1)*G*mass2))
          mss1 = jl(ijob)%mass_list(1)

          call set_filename(jl(ijob)%mass_list(2),jl(ijob)%eta, &
                            jl(ijob)%alpha)
          call read_model()

          ! Normalise frequency by Keplerian break-up rotation rate:
          freq2 = jl(ijob)%w2(i)*sqrt(p1D(1)*radius2/(rho1D(1)*G*mass2))
          mss2 = jl(ijob)%mass_list(2)

          ! set the azimuthal order:
          m = jl(ijob)%m(i)

          ! set the target l value:
          ltarget = jl(ijob)%ltarget(i)

          call system("echo >> run.log")
          ! follow the frequency as a function of mass:
          do imass=3,jl(ijob)%nmass
            mss3 = jl(ijob)%mass_list(imass)
            freq3 = freq1+(mss3-mss1)*(freq2-freq1)/(mss2-mss1)
            call set_filename(jl(ijob)%mass_list(imass), &
                              jl(ijob)%eta,jl(ijob)%alpha)
            call init_model()

            ! set the shift:
            shift_real = freq3*sqrt(rho1D(1)*G*mass2/(p1D(1)*radius2))
            shift_imag = 0d0

            call init_a()
            call init_order()
            call init_bc_flag()
            call init_bc_range()
            call init_dir()
            call run_arncheb(dcmplx(shift_real,shift_imag))
            call write_output()

            shift_real = dreal(omega(1))
            lbest = ldom(1)
            do j=2,nsol_out
              if (abs(ldom(j)-ltarget).lt.abs(lbest-ltarget)) then
                shift_real = dreal(omega(j))
                lbest = ldom(j)
              endif
            enddo

            ! output results in run.log:
            write (string,101) jl(ijob)%mass_list(imass), &
                               jl(ijob)%eta,              &
                               jl(ijob)%alpha,            &
                               m, lbest, shift_real
            call system(string)

            ! prepare for next iteration:
            freq1 = freq2
            freq2 = shift_real*sqrt(p1D(1)*radius2/(rho1D(1)*G*mass2))
            mss1 = mss2
            mss2 = mss3
          enddo
        enddo
      enddo

101   format("echo '",3(0pf5.2,X),I3,X,I3,X,1pe22.15,"' >> run.log")
      end program TOP
