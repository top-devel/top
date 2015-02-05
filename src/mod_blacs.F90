#include "config.h"
       module mod_blacs

       integer, save :: ictxt, ictxt1D
       integer, save :: nprocs, iproc
       integer, save :: nrows, ncols, irow, icol

       character, save :: topo = " "
       character, save :: blacs_order = "R"

!------------------------------------------------------------------------------
! ictxt  = context number used by BLACS for the 2D grid
! ictxt1D= context number used by BLACS for the 1D grid
! nprocs = number of processes avalaible (from command line)
! iproc  = number of the current process
! nrows  = number of rows in the BLACS process grid
! ncols  = number of columns in the BLACS process grid
! irow   = row number of the current process in the BLACS process grid
! icol   = column number of the current process in the BLACS process grid
! topo   = topology of the parallel machine (" " = system dependent default)
! blacs_order = order in which to distribute the processes on the process grid.
!------------------------------------------------------------------------------

contains
!------------------------------------------------------------------------------
! Obtains the number of processes, the number of the current process and
! a BLACS context number.  Set up a first 1D process grid which will be
! used for communicating input parameters and for doing sparse matrix
! multiplications.
!------------------------------------------------------------------------------
       subroutine start_blacs()
#ifdef USE_MPI
       implicit none
       integer :: nrows1D
       nrows1D = 1
       call BLACS_PINFO(iproc,nprocs)
       call BLACS_GET(-1,0,ictxt1D)
       call BLACS_GRIDINIT(ictxt1D,blacs_order,nrows1D,nprocs)
#else
       iproc = 0
       nprocs = 1
#endif
       end subroutine start_blacs

!------------------------------------------------------------------------------
! Starts the BLACS process grid, gets relevant information and checks
! the number of processes.
!------------------------------------------------------------------------------
       subroutine start_blacs_grid()
       implicit none

#ifdef USE_MPI
       call BLACS_GET(-1,0,ictxt)
#endif
       ! check to see that the number of processes match the size of the
       ! process grid (as defined in the input file)
       if (nrows*ncols.ne.nprocs) then
         print*,"ERROR: Number of processes does not match"
         print*,"       size of process grid.  Exiting."
         print*,"ncols = ", ncols
         print*,"nrows = ", nrows
         print*,"nprocs= ", nprocs
#ifdef USE_MPI
         call BLACS_ABORT(ictxt,1)
#else
         stop
#endif
       endif
#ifdef USE_MPI
       call BLACS_GRIDINIT(ictxt,blacs_order,nrows,ncols)
       call BLACS_GRIDINFO(ictxt,nrows,ncols,irow,icol)
#endif
       end subroutine start_blacs_grid

!------------------------------------------------------------------------------
! This subroutine ends the BLACS process grids and should only be called
! at the end of the program.
!------------------------------------------------------------------------------
       subroutine stop_blacs()
#ifdef USE_MPI
       call BLACS_GRIDEXIT(ictxt)
       call BLACS_GRIDEXIT(ictxt1D)
       call BLACS_EXIT(0)
#endif
       end subroutine stop_blacs
!--------------------------------------------------------------------
        subroutine string2iarray(str,iarray,narray)

        implicit none
        character*(*), intent(in) :: str
        integer, intent(in) :: narray
        integer, intent(inout) :: iarray(narray)
        integer i

        do i=1,min(narray,len(str))
          iarray(i)=ichar(str(i:i))
        enddo
        do i=min(narray,len(str))+1,narray
          iarray(i)=32 ! the ascii value of a blank space
        enddo

        end subroutine string2iarray
!--------------------------------------------------------------------
        subroutine iarray2string(str,iarray,narray)

        implicit none
        character*(*), intent(inout) :: str
        integer, intent(in) :: narray
        integer, intent(in) :: iarray(narray)
        integer i

        do i=1,min(narray,len(str))
          str(i:i)=char(iarray(i))
        enddo
        do i=min(narray,len(str))+1,len(str)
          str(i:i) = " "
        enddo

        end subroutine iarray2string
!------------------------------------------------------------------------------
       end module mod_blacs
