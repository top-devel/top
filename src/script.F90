#include "config.h"
      module script

      integer, parameter :: list_size_max = 500
      type JOBLIST
        double precision, pointer :: mass_list(:)
        integer :: nmass
        double precision :: eta
        double precision :: alpha
        double precision :: w1(list_size_max)
        double precision :: w2(list_size_max)
        integer :: ltarget(list_size_max)
        integer :: m(list_size_max)
        integer :: list_size
      end type

      type(JOBLIST), allocatable, save :: jl(:)
      integer, save :: njobs


contains
!----------------------------------------------------------------------
! This reads the script file with the different jobs that need to
! be done.
!----------------------------------------------------------------------
! File name: script_TOP
!----------------------------------------------------------------------
! File format:
!
! model eta mass1 mass2 ...
! follow freq1 freq2 ltarget m
!----------------------------------------------------------------------
! For each mass command, the program loads a new list of model masses
! to be used when tracking the frequencies.
!   IMPORTANT: the first two masses are associated with freq1 and freq2.
!              Therefore the program will start calculating frequencies
!              from mass3 and onwards.
! 
! The follow command gives two frequencies, a target l value and an
! m value with which to track frequencies.
!----------------------------------------------------------------------
      subroutine read_script()

      implicit none
      integer i, ijob, istart
      character*(1024) oneline

      open(unit=9, file="script_TOP", status="unknown")
      njobs = 0
      do
        read(9,'(a1024)',end=9) oneline
        if (oneline(1:1).eq."#") cycle
        if (index(oneline,'mass').ne.0) njobs = njobs + 1
      enddo
9     rewind(9)
      if (allocated(jl)) deallocate(jl)
      allocate(jl(njobs))
      ijob = 0
      do
        read(9,'(a1024)',end=11) oneline

        ! Ignore comments
        if (oneline(1:1).eq."#") cycle

        ! Read a line which contains the keyword 'mass'
        istart = index(oneline,'mass') ! ignore what comes before 'mass'
        if (istart.ne.0) then
          ijob = ijob + 1
          jl(ijob)%nmass = 0
          jl(ijob)%list_size = 0
          do i=istart,len(oneline)
            ! count the number of decimal points to determine the number
            ! of different masses:
            if (oneline(i:i).eq.".") jl(ijob)%nmass = jl(ijob)%nmass+1
          enddo
          jl(ijob)%nmass = jl(ijob)%nmass-2
          allocate(jl(ijob)%mass_list(jl(ijob)%nmass))
          read(oneline(istart+4:),*) jl(ijob)%eta, jl(ijob)%alpha, &
                        (jl(ijob)%mass_list(i),i=1,jl(ijob)%nmass)
          cycle
        endif

        ! Read a line which contains the keyword 'follow'
        istart = index(oneline,'follow')
        if (istart.ne.0) then
          if (jl(ijob)%list_size.eq.list_size_max) then
           print*,"WARNING: some of the jobs have been cancelled."
           print*,"         Please increase list_size_max in script.f90"
           cycle
          endif
          jl(ijob)%list_size = jl(ijob)%list_size + 1
          i = jl(ijob)%list_size
          read(oneline(istart+6:),*) jl(ijob)%w1(i), jl(ijob)%w2(i), &
                                     jl(ijob)%ltarget(i), jl(ijob)%m(i)
          cycle
        endif
      enddo
11    close(9)
      end subroutine read_script
!----------------------------------------------------------------------
      subroutine set_filename(mass_value,eta_value,alpha_value)

      use inputs
      implicit none
      double precision, intent(in) :: mass_value, eta_value, alpha_value
      character*(5) s_mass,s_eta,s_alpha

      write(s_mass, 101) mass_value
      write(s_eta,  101) eta_value
      write(s_alpha,101) alpha_value
101   format(0pf5.2)

      filename = trim(adjustl(modeldir))//"m_"// &
                 trim(adjustl(s_mass))  //"_" // &
                 trim(adjustl(s_eta))   //"_" // &
                 trim(adjustl(s_alpha))
      end subroutine set_filename
!----------------------------------------------------------------------
      subroutine dump_script()

      implicit none
      integer ijob, i
      open(unit=11,file='script_dumped',status='unknown')
      write(11,*) 'njobs: ',njobs
      do ijob=1,njobs
        write(11,*) '***********************'
        write(11,*) 'nmass: ',jl(ijob)%nmass
        do i=1,jl(ijob)%nmass
          write(11,*) '   ',jl(ijob)%mass_list(i)
        enddo
        write(11,*) 'eta:   ',jl(ijob)%eta
        write(11,*) 'alpha: ',jl(ijob)%alpha
        write(11,*) 'list_size: ',jl(ijob)%list_size
        do i=1,jl(ijob)%list_size
          write(11,*) jl(ijob)%w1(i), jl(ijob)%w2(i), &
                      jl(ijob)%ltarget(i), jl(ijob)%m(i)
        enddo
      enddo
      close(11)

      end subroutine dump_script
!----------------------------------------------------------------------
      end module script
