#include "config.h"
      module spectrum

      type MODE_ID
        integer :: n
        integer :: l
        integer :: m
        double precision :: freq
        double precision :: distance
        logical :: valid
        logical :: is_input
        integer :: ldom
        integer :: isol
        character*(512) :: path
      end type MODE_ID

      type(MODE_ID), allocatable, save :: mode_list(:)
      integer, save :: nfreq, nfreq_input
      double precision, save :: mat(5,5), v(5), w(5)

contains

!-----------------------------------------------------------------------------
! This reads a list of already identified modes from a file called 'input_list'
!-----------------------------------------------------------------------------
      subroutine read_input_list()

      use inputs
      implicit none
      character*(80) st
      double precision auxr
      integer aux1, aux2, aux3
      logical bool
      integer i

      open(unit=9,file='input_list',status='old')
      nfreq_input = 0
      do
        read(9,'(a)',end=99) st
        call isblank(st,bool)
        if (.not.bool) nfreq_input = nfreq_input + 1
      enddo
99    rewind(9)
      print*,"Number of entries: ", nfreq_input
      nfreq = nfreq_input + (1+(n_end-n_start)/n_incr) &
                          * (1+(l_end-l_start)/l_incr) &
                          * (1+(m_end-m_start)/m_incr)
      if (allocated(mode_list)) deallocate(mode_list)
      print*,"nfreq",nfreq
      allocate(mode_list(nfreq))
      i = 0
      do
        read(9,'(a)',end=98) st
        call isblank(st,bool)
        if (.not.bool) then
          i = i + 1
          read(st,*) mode_list(i)%freq, mode_list(i)%n, &
                     mode_list(i)%l, mode_list(i)%m
          mode_list(i)%distance = 0d0
          mode_list(i)%valid    = .true.
          mode_list(i)%is_input = .true.
          mode_list(i)%ldom     = -1
          mode_list(i)%isol     = 0
          mode_list(i)%path     = ""
        endif
      enddo
98    close(9)
      end subroutine read_input_list

!-----------------------------------------------------------------------------
! This completes mode_list by adding those modes which still need to be
! calculated so as to comply with the requirements in the dati file.
!-----------------------------------------------------------------------------
      subroutine complete_list()

      use inputs

      implicit none
      integer tn, tl, tm, i, j
      double precision term, total
      logical is_new

      i = nfreq_input
      do tn = n_start, n_end, n_incr
        do tl = l_start, l_end, l_incr
          do tm = m_start, m_end, m_incr
            j = 1
            is_new = .true.
            total = 0d0
            do while(is_new.and.(j.le.nfreq_input))
              term = dble((mode_list(j)%n-tn)**2) &
                   + dble((mode_list(j)%l-tl)**2) &
                   + dble((mode_list(j)%m-tm)**2)
              if (abs(term).lt.0.5d0) then
                is_new = .false.
              else
                total = total + term
              endif
              j = j + 1
            enddo
            if (is_new) then
              i = i + 1
              mode_list(i)%n = tn
              mode_list(i)%l = tl
              mode_list(i)%m = tm
              mode_list(i)%distance = total
              mode_list(i)%valid    = .false.
              mode_list(i)%is_input = .false.
            endif
          enddo
        enddo
      enddo
      nfreq = i
      end subroutine complete_list

!------------------------------------------------------------------------------
!  This initialises the matrices involved in the least-squares calculation
!  of the coefficients in the asymptotic formula.
!------------------------------------------------------------------------------
      subroutine init_mat()
      implicit none

      double precision aux(6)
      integer i, j, k

      mat = 0d0
      v = 0d0

      do i=1,nfreq_input
        aux(1) = mode_list(i)%n
        aux(2) = mode_list(i)%l
        aux(3) = abs(mode_list(i)%m)
        aux(4) = 1d0
        aux(5) =-mode_list(i)%m
        aux(6) = mode_list(i)%freq

        do j=1,5
          do k=j,5
            mat(j,k) = mat(j,k) +aux(j)*aux(k)
          enddo
        enddo

        do j=1,5
          v(j) = v(j) + aux(j)*aux(6)
        enddo
      enddo

      ! Make the matrix symmetric
      do j=2,5
        do k=1,j-1
          mat(j,k) = mat(k,j)
        enddo
      enddo

      call solve()
      end subroutine init_mat

!------------------------------------------------------------------------------
!  This updates the matrices involved in the least-squares calculation
!  of the coefficients in the asymptotic formula.
!------------------------------------------------------------------------------
      subroutine update_mat(ifreq)

      implicit none
      integer, intent(in) :: ifreq
      double precision aux(6)
      integer j, k

      ! if correct identification has failed, do not modify the matrices:
      if (.not.mode_list(ifreq)%valid) return

      aux(1) = mode_list(ifreq)%n
      aux(2) = mode_list(ifreq)%l
      aux(3) = abs(mode_list(ifreq)%m)
      aux(4) = 1d0
      aux(5) =-mode_list(ifreq)%m
      aux(6) = mode_list(ifreq)%freq

      do j=1,5
        do k=1,5
          mat(j,k) = mat(j,k) +aux(j)*aux(k)
        enddo
      enddo

      do j=1,5
        v(j) = v(j) + aux(j)*aux(6)
      enddo

      call solve()
      end subroutine update_mat

!------------------------------------------------------------------------------
! This uses the parameters stored in the w(:) array to calculate the
! frequency for the mode identification ifreq.
!------------------------------------------------------------------------------
      double precision function find_freq(ifreq)

      implicit none
      integer ifreq
      find_freq = w(1)*mode_list(ifreq)%n    &
                + w(2)*mode_list(ifreq)%l    &
                + w(3)*abs(mode_list(ifreq)%m) &
                - w(5)*mode_list(ifreq)%m    &
                + w(4)
      return
      end function find_freq
!------------------------------------------------------------------------------
! This program sorts an array based on the heapsort method.  It is copied from:
!
! "Numerical Recipes in Fortran: The Art of Scientific Computing"
! Second Edition.
! William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery
! Cambridge University Press, 1992.
!------------------------------------------------------------------------------
! This routine sorts mode_list such that
! compare(mode_list(i),mode_list(i+1)) is true throughout.
!------------------------------------------------------------------------------
       subroutine sort_list(compare)

       implicit none
       logical, external :: compare
       integer i, ir, j, l
       type(MODE_ID) :: a_mode

       if (nfreq.lt.2) return
       
         ! The index l will be decremented from its initial value down to 1
         ! during the "hiring" (heap creation) phase.  Once it reaches 1, the
         ! index ir will be decremented from its initial value down to 1 during
         ! the "retirement-and-promotion" (heap selection) phase.

       l = nfreq/2+1
       ir = nfreq
       do
         if (l.gt.1) then                ! Still in hiring phase
           l=l-1
           a_mode=mode_list(l)
         else                            ! In retiremnet-and-promotion phase.
           a_mode=mode_list(ir)          ! Clear a space at end of aa_modey.
           mode_list(ir)=mode_list(1)    ! Retire the top of the heap into it
           ir=ir-1                       ! Decrease the size of the corporation.
           if(ir.eq.1)then               ! Done with the last promotion.
             mode_list(1)=a_mode         ! The least competent worker of all!
             return
           endif
         endif
         i=l                             ! Whether in the hiring phase or
                                         ! promotion phase, we here set up to
                                         ! sift down element a_mode to its proper
                                         ! level.
         j=l+l
         do while (j.le.ir)
           if (j.lt.ir) then
             if (compare(mode_list(j),mode_list(j+1))) j=j+1 ! Compare to the better underling
           endif
           if (compare(a_mode,mode_list(j))) then         ! Demote a_mode
             mode_list(i)=mode_list(j)
             i=j
             j=j+j
           else                          ! This is a_mode's level.  Set j to
                                         ! terminate the sift-down.
             j=ir+1
           endif
         enddo
         mode_list(i)=a_mode             ! Put a_mode into its slot.
       enddo
       end subroutine

!-------------------------------------------------------------
!  This compares two MODE_ID type elements.  It returns
!  .true. if element1%distance < element2%distance
!-------------------------------------------------------------
       logical function dcompare(mode1,mode2)

       implicit none
       type(MODE_ID), intent(in) :: mode1, mode2
       dcompare = (mode1%distance.lt.mode2%distance)
       return
       end function dcompare

!-------------------------------------------------------------
!  This compares two MODE_ID type elements.  It returns
!  .true. if mode1(l,n,m) comes before mode2(l,n,m) in
!  terms of lexicographic order. 
!-------------------------------------------------------------
       logical function lcompare(mode1,mode2)

       implicit none
       type(MODE_ID), intent(in) :: mode1, mode2

       if (mode1%l.lt.mode2%l) then
         lcompare = .true.
         return
       endif

       if (mode1%l.gt.mode2%l) then
         lcompare = .false.
         return
       endif

       if (mode1%n.lt.mode2%n) then
         lcompare = .true.
         return
       endif

       if (mode1%n.gt.mode2%n) then
         lcompare = .false.
         return
       endif

       lcompare = (mode1%m.lt.mode2%m)
       return

       end function lcompare

!-------------------------------------------------------------
!  This subroutine determines whether a line is blank
!-------------------------------------------------------------
       subroutine isblank(string,bool)
       implicit none
       character*(*) string
       character*(*), parameter :: blanks = ' '
       logical bool
       integer i

       bool = .true.
       do i=1,len(string)
         if (index(blanks,string(i:i)) == 0) then
           bool = .false.
           return
         endif
       enddo
       return
       end subroutine
!------------------------------------------------------------------------------
! This swaps two elements.
!------------------------------------------------------------------------------

       subroutine swap (nombre1, nombre2)
       implicit none
       double precision nombre1, nombre2, aux

       aux = nombre1
       nombre1 = nombre2
       nombre2 = aux
       return
       end subroutine

!------------------------------------------------------------------------------
! This solves a system of equations.  If the matrix is non-inversible, it
! gives an error message and stops.
!------------------------------------------------------------------------------
       subroutine solve()
       implicit none
       integer, parameter ::  debut = 1, taille = 5
       double precision matrix(5,5), v_temp(5)
       double precision r
       integer icolumn, i, j

       ! copy the original matrices to avoid overwriting them
       do i=1,5
         do j=1,5
           matrix(i,j) = mat(i,j)
         enddo
         v_temp(i) = v(i)
       enddo

       do icolumn=debut,taille
         ! Search for first non-zero value
         do i=icolumn,taille
           if (matrix(i,icolumn).ne.0d0) exit
         enddo

         ! If the matrix is singular, then stop:
         if (i.gt.taille) then
           print*,"Singular matrix in spectrum.f90:"
           print*,"Try adding extra frequencies to input_list."
           stop
         endif

         ! If need be, swap the two rows:
         if (i.ne.icolumn) then
           do j=icolumn,taille
             call swap(matrix(i,j),matrix(icolumn,j))
           enddo
           call swap(v_temp(i),v_temp(icolumn))
         endif
         
         ! normalise row icolumn so that its pivot is 1:
         r = 1d0/matrix(icolumn,icolumn)
         do j=icolumn+1,taille
           matrix(icolumn,j) = matrix(icolumn,j)*r
         enddo
         v_temp(icolumn) = r*v_temp(icolumn)
         
         ! cancel out the elements beneath the pivot:
         do i=icolumn+1,taille
           r = -matrix(i,icolumn)
           do j=icolumn+1,taille
             matrix(i,j) = matrix(i,j) + r*matrix(icolumn,j)
           enddo
           v_temp(i) = v_temp(i) + r*v_temp(icolumn)
         enddo
       enddo
       
       ! By this point, the matrix should be in upper triangular form.
       ! We finish solving the system:
       do i=taille,debut,-1
         w(i) = v_temp(i)
         do j=i+1,taille
           w(i) = w(i) - matrix(i,j)*w(j)
         enddo
       enddo

       return
       end subroutine

!-------------------------------------------------------------------------

      end module spectrum
