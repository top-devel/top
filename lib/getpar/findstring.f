c     $Author: raus $
c     $Revision: 1.2 $
c     $Log: GP_find.f,v $
c Revision 1.2  1993/04/21  16:53:35  raus
c size of char strings up to 1024
c
c Revision 1.1  1993/04/20  17:43:17  raus
c 1st save; find_string to findstring
c
c***********************************************************************
      INTEGER FUNCTION FINDSTRING(name_variable)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      CHARACTER*(*) name_variable
      INTEGER i, count, l, l1
      LOGICAL flag_ok
      CHARACTER*1 char
      CHARACTER*1024 name
c***
      l1 = LEN(name_variable)
      FINDSTRING = 0
      if (nlines .lt. 1) return
      do count=1,nlines
         l = name_var_len(count)
         if (l1 .ne. l) then
            flag_ok = .false.
         else 
            name = lines_of_text(count)
            flag_ok = .true.
            do i=1,l1
               char = name_variable(i:i)
               if (name(i:i) .ne. char) flag_ok = .false.
            enddo
         endif 
         if (flag_ok) FINDSTRING = count
      enddo
c***
      return
      end
