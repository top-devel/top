c $Author: raus $
c $Revision: 1.2 $
c $Log: GP_term.f,v $
c Revision 1.2 1993/04/21 16:58:05 raus
c fixed a bug (i counter)
c
c Revision 1.1 1993/04/20 17:48:27 raus
c 1st save
c
c***********************************************************************
      SUBROUTINE term_string(string)
c***


      INCLUDE 'undefined.h'
      INTEGER len_string, i, ibegin
      CHARACTER*(*) string
      CHARACTER*1 end_char, blank, char
      DATA end_char/' '/, blank/' '/
c***
      len_string = LEN(string)
      ibegin = 1
      i = 1
      char = string(1:1)
      do while (char .ne. blank .and. i .le. len_string)
         ibegin = ibegin + 1
         char = string(ibegin:ibegin)
      enddo
c***
c--- print*,'>>>',string,'<<<'
c***
      do i=ibegin,len_string
         string(i:i) = end_char
      enddo


c***
c--- print*,'>>>',string,'<<<'
c***
      return
      end
