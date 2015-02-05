       module itemlist

       integer, private,  parameter :: item_length = 256
       character*(item_length), public, allocatable, save :: &
                                      item(:), value(:)
       integer, public, save :: nitem

       ! Important: when using fetch, always send it into a
       !            variable.  Doing "print*,fetch('var',2)"
       !            causes an error (at least with gfortran).

       public :: fetch
       interface fetch
         module procedure dfetch
         module procedure ifetch
         module procedure rfetch
         module procedure sfetch
       end interface
       
contains

!---------------------------------------------------------------------
!  This fetches a double precision value corresponding to
!  a specified item in the item list.
!---------------------------------------------------------------------
       double precision function dfetch(item_name, default_value)

       implicit none
       character*(*) item_name
       double precision default_value
       integer i

       dfetch = default_value
       do i=1,nitem
         if (trim(item_name).eq.trim(item(i))) then
           read(value(i),*,err=101) dfetch
           return
         endif
       enddo
       return
101    write(*,102) trim(item_name)
102    format('Warning: item "',a,'" does not contain a double ',&
              'precision value.')
       return
       end function

!---------------------------------------------------------------------
!  This fetches an integer value corresponding to
!  a specified item in the item list.
!---------------------------------------------------------------------
       integer function ifetch(item_name, default_value)

       implicit none
       character*(*) item_name
       integer default_value
       integer i

       ifetch = default_value
       do i=1,nitem
         if (trim(item_name).eq.trim(item(i))) then
           read(value(i),*,err=101) ifetch
           return
         endif
       enddo
       return
101    write(*,102) trim(item_name)
102    format('Warning: item "',a,'" does not contain an integer value.')
       return
       end function

!---------------------------------------------------------------------
!  This fetches a real (or single precision) value corresponding to
!  a specified item in the item list.
!---------------------------------------------------------------------
       real function rfetch(item_name, default_value)

       implicit none
       character*(*) item_name
       real default_value
       integer i

       rfetch = default_value
       do i=1,nitem
         if (trim(item_name).eq.trim(item(i))) then
           read(value(i),*,err=101) rfetch
           return
         endif
       enddo
       return
101    write(*,102) trim(item_name)
102    format('Warning: item "',a,'" does not contain a real ',&
              'value.')
       return
       end function

!---------------------------------------------------------------------
!  This fetches a string corresponding to
!  a specified item in the item list.
!---------------------------------------------------------------------
       character*(item_length) function sfetch(item_name, default_value)

       implicit none
       character*(*) item_name
       character*(*) default_value
       integer i

       sfetch = default_value
       do i=1,nitem
         if (trim(item_name).eq.trim(item(i))) then
           sfetch = value(i)
           return
         endif
       enddo
       return
       end function

!---------------------------------------------------------------------
!  This subroutine parses one line of text ("oneline") of the form
!           item1 = value1   item2 = value2 ...
!  into a series of item/value pairs which are both stored as strings
!  of length item_length.  It does two sweeps: the first determines
!  the number of item/value pairs, and the second reads the
!  names and values of the different pairs and stores them.
!---------------------------------------------------------------------

       subroutine readItemlist(oneline)

       implicit none
       integer length, i, istart, ifinish, mode, nn
       character*(*) oneline
       character*(item_length) padding

       ! This is a string full of blank spaces - it
       ! serves as a padding to ensure that the different
       ! items and values have trailing blank spaces
       ! rather than something else.
       do i=1,item_length
         padding(i:i) = " "
       enddo

       length = len(oneline)
       nitem = 0
       do i=1,length
         if (oneline(i:i).eq."=") nitem = nitem + 1
       enddo

       if (allocated(item)) deallocate(item)
       if (allocated(value)) deallocate(value)
       allocate(item(nitem),value(nitem))
       do nn=1,nitem
         item(nn)  = padding
         value(nn) = padding
       enddo

       ! As the subroutine reads "oneline" one character at a time,
       ! it can be in different modes, depending on what it is
       ! looking for:
       !
       !  mode = 1: looking for item name
       !  mode = 2: reading item name
       !  mode = 3: finished reading item name, looking for "="
       !  mode = 4: found "=", looking for item value
       !  mode = 5: reading value
       !
       !  nn      = number of item name or value
       !  istart  = start of item name or value
       !  ifinish = end of item name or value

       istart = 0
       ifinish = 0
       nn = 1
       mode = 1

       do i=1,length

         select case(mode)

           case (1)

             if (oneline(i:i).eq."=") then
               if (nn.gt.1) then
                 item(nn) = value(nn-1)
                 value(nn-1) = padding
               else
                 item(nn) = padding
               endif
               mode = 4
             elseif (oneline(i:i).ne." ") then
               istart = i
               mode = 2
             endif

           case (2)

             if (oneline(i:i).eq." ") mode = 3
             if (oneline(i:i).eq."=") mode = 4
             if (mode.ne.2) then 
               ifinish = i-1
               item(nn) = oneline(istart:ifinish)//padding
             endif

           case (3)

             if (oneline(i:i).eq."=") then
               mode = 4
             elseif (oneline(i:i).ne." ") then
               istart = i
               mode = 2
             endif

           case (4)

             if (oneline(i:i).eq."=") then
               value(nn) = padding
               nn = nn + 1
               item(nn) = padding
               mode = 4
             elseif (oneline(i:i).ne." ") then
               istart = i
               mode = 5
             endif

           case (5)

             if (oneline(i:i).eq."=") then
               ifinish = i-1
               value(nn) = padding
               nn = nn + 1
               item(nn) = oneline(istart:ifinish)//padding
               mode = 4
             elseif (oneline(i:i).eq." ") then
               ifinish = i-1
               value(nn) = oneline(istart:ifinish)//padding
               nn = nn + 1
               mode = 1
             endif
             if (nn.gt.nitem) exit

         end select

       enddo

       if (mode.eq.5) then
         value(nn) = oneline(istart:length)//padding
       endif

       return

       end subroutine readItemlist

!-------------------------------------------------------------------
!  This displays the content of the itemlist for the purposes of
!  debbuging the program.
!-------------------------------------------------------------------
       subroutine displayItemlist()

       implicit none
       integer i

       do i=1,nitem
         write(*,101) i, trim(item(i)), trim(value(i))
       enddo

101    format(I3,X,A,'=',A)
       end subroutine displayItemlist

!-------------------------------------------------------------------
!  Another useful method.  This parses a string str using separators
!  and stores the results in two arrays: num_list and str_list.
!-------------------------------------------------------------------

       subroutine parse(str,separators,num_list,str_list,taille)

       ! cette routine prend les characteres de "separators"
       ! pour diviser la chaine "str".  Le resultat est
       ! mis dans num_list et str_list.

       implicit none
       integer taille
       double precision num_list(taille)
       character*(*) str_list(taille)
       character*(*) str, separators
       integer i,istart,i_list
       logical issep,issepn

       if (taille.eq.0) return
       issepn = .true.
       i_list = 0
       do i=1,len(str)
         issep = issepn
         issepn = (index(separators,str(i:i)) /= 0)
         if (issep.and.(.not.issepn)) istart = i
         if ((.not.issep).and.issepn) then
           i_list = i_list + 1
           str_list(i_list) = str(istart:i-1)
           num_list(i_list) = 0d0
           read(str_list(i_list),*,err=10,end=10) num_list(i_list)
10         if (i_list.eq.taille) exit
         endif
       enddo
       if (i_list.eq.taille) return
       issep = issepn
       issepn = .true.
       if ((.not.issep).and.issepn) then
         i_list = i_list + 1
         str_list(i_list) = str(istart:i-1)
         num_list(i_list) = 0d0
         read(str_list(i_list),*,err=11,end=11) num_list(i_list)
11     endif
       taille = i_list
       return
       end subroutine
!---------------------------------------------------------------------
       end module
