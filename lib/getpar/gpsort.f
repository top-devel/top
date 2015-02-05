c     $Author: raus $
c     $Revision: 1.2 $
c     $Log: GP_sort.f,v $
c Revision 1.2  1993/04/21  16:57:23  raus
c string size up to 1024 char
c
c Revision 1.1  1993/04/20  17:47:06  raus
c 1st save
c
c*****************************************************************
      SUBROUTINE GPsort
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
c***
      call sortwithkey(nlines,name_var_len,lines_of_text)
c***
      return
      end
