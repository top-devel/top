#include "config.h"
      program TOP

      use model
      use matrices
      use eigensolve
      use inputs
      use postproc
      use udirectory

      implicit none
      character*(128) string
      character*(6) sparity
      double precision, allocatable :: shift_table(:)
      integer nshifts, i

      call read_inputs()
      call init_model()

      if (index(stamp,"even").eq.0) then
        sparity = "l odd "
      else
        sparity = "l even"
      endif

      nshifts = 0
      open(unit=3,file='shifts',status='old')
      do
        read(3,*,end=99)
        nshifts = nshifts + 1
      enddo 
99    rewind(3)
      allocate(shift_table(nshifts))
      do i=1,nshifts
        read(3,*) shift_table(i)
      enddo
      close(3)

      ! Normally, I shouldn't have to do this, but on
      ! some computers, this could help.
      do m = m_start,m_end,m_incr
        do i=1,nshifts
          shift_real = shift_table(i)
          call init_dir()
        enddo
      enddo

      write(string,101) mass, eta, alpha
      call system(trim(string))
      do m = m_start,m_end,m_incr
        call init_a()
        call init_order()
        call init_bc_flag()
        call init_bc_range()

        write(string,102) m, sparity
        call system(trim(string))
        do i=1,nshifts
          shift_real = shift_table(i)
          ! I still need to put the correct dir_path in
          ! memory:
          call init_dir()
          call run_arncheb(dcmplx(shift_real,shift_imag))
          call write_output()
          print*,shift_real
        enddo
      enddo

101   format("echo '# mass=",0pf5.2,", eta=",0pf4.2, &
             ", alpha=",0pf4.2,"' >> valps") 
102   format("echo '# m=",I2,", ",a6,"' >> valps")
      end program TOP
