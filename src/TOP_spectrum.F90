#include "config.h"
      program TOP

      use model
      use matrices
      use eigensolve
      use inputs
      use postproc
      use udirectory
      use spectrum

      implicit none
      character*(600) st
      character*(1) svalid
      integer i, j, isol, ltarget
      double precision freq

      call read_inputs()
      call init_model()
      call read_input_list()
      call complete_list()
      call sort_list(dcompare)
      call init_mat()
    
      call system("echo '-' >> vecp_list.unsorted") 
      do i = nfreq_input+1,nfreq

        ! set some initial values
        shift_real   = find_freq(i)
        m       = mode_list(i)%m
        iparity = modulo(mode_list(i)%n,2)
        ltarget = 2*mode_list(i)%l+abs(m)+iparity

        call init_dir()
        call init_a()
        call init_order()
        call init_bc_flag()
        call init_bc_range()
        call run_arncheb(dcmplx(shift_real,shift_imag))
        call write_output()

        isol = 1 
        do j=2,nsol_out
          if (abs(ldom(isol)-ltarget).gt.abs(ldom(j)-ltarget)) isol = j
          if (abs(ldom(isol)-ltarget).eq.abs(ldom(j)-ltarget)) then
            if (abs(dreal(omega(Ndex(isol)))-shift_real).gt.  &
                abs(dreal(omega(Ndex(j)))-shift_nreal)) isol = j
          endif
        enddo
        mode_list(i)%freq = dreal(omega(Ndex(isol)))
        mode_list(i)%valid = (abs(ldom(isol)-ltarget).le.ltol)
        mode_list(i)%path = trim(dir)//'vecp.gz'
        mode_list(i)%isol = isol
        mode_list(i)%ldom = ldom(isol)
        write(st,101) mode_list(i)%freq, &
                      mode_list(i)%n,    &
                      mode_list(i)%l,    &
                      mode_list(i)%m,    &
                      ldom(isol)
        call system(trim(st))
        write(st,102) isol,trim(dir)//'vecp.gz'
        call system(trim(st))
      enddo
      call system("echo >> run.log.unsorted")
      call sort_list(lcompare)

      open(unit=22,file='run.log',status='unknown')
      do i=1,nfreq
        if (mode_list(i)%is_input) cycle
        if (mode_list(i)%valid) then
          svalid = 'T'
        else
          svalid = 'F'
        endif
        write(22,103) mode_list(i)%freq, &
                      mode_list(i)%n,    &
                      mode_list(i)%l,    &
                      mode_list(i)%m,    &
                      mode_list(i)%ldom
      enddo
      close(22)

      open(unit=22,file='vecp_list',status='unknown')
      write(22,'(a)') "-"
      do i=1,nfreq
        if (mode_list(i)%is_input) cycle
        write(22,104) mode_list(i)%isol, &
                      trim(mode_list(i)%path)
      enddo
      close(22)

101   format("echo '",0pf19.15,2(X,I2),X,I3," | ",I3, &
             "' >> run.log.unsorted")
103   format(0pf19.15,2(X,I2),X,I3," | ",I3)
102   format("echo '",I2,X,A,"' >> vecp_list.unsorted")
104   format(I2,X,A)
      end program TOP
