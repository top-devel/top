#include "config.h"
      program TOP

      use model
      use matrices
      use eigensolve
      use inputs
      use postproc
      use udirectory
      use freqlist

      implicit none
      integer i, isol, iisol, lbest
      double complex best_freq
      character*(256) str

      call read_inputs()
      call init_model()
      call readfreq()

      do i=1,ndata
        iparity = modulo(ll(i)+abs(mm(i)),2)
        shift_real = fr(i)
        shift_imag = fi(i)
        m = mm(i)
        call init_dir()
      enddo

      do i=1,ndata
        iparity = modulo(ll(i)+abs(mm(i)),2)
        shift_real = fr(i)
        shift_imag = fi(i)
        m = mm(i)
        call init_a()
        call init_order()
        call init_bc_flag()
        call init_bc_range()
        call run_arncheb(dcmplx(shift_real,shift_imag))
        call init_dir()
        call write_output()

        ! find solution closest to target l and frequency
        lbest = 100000
        best_freq = 0d0
        do iisol=1,nsol_out
          isol = Ndex(iisol)
          if (abs(lbest-ll(i)).eq.abs(ldom(iisol)-ll(i))) then
            if (abs(best_freq-dcmplx(fr(i),fi(i))).gt. &
                abs(omega(isol)-dcmplx(fr(i),fi(i)))) then
              best_freq = omega(isol)
              lbest = ldom(iisol)
            endif
          endif
          if (abs(lbest-ll(i)).gt.abs(ldom(iisol)-ll(i))) then
            best_freq = omega(isol)
            lbest = ldom(iisol)
          endif
        enddo
        write(str,101) ll(i), mm(i), fr(i), fi(i), lbest, best_freq
        call system(str)
      enddo

101   format('echo "',2(I3,X),2(0pf9.6,X),I2,2(X,0pf18.15),'" >> valps')
      end program TOP
