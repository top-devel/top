c     $Author: raus $
c     $Revision: 1.2 $
c     $Log: GP_fetch.f,v $
c Revision 1.2  1993/04/21  16:52:26  raus
c size of the char function up to 1024
c
c Revision 1.1  1993/04/20  17:41:37  raus
c 1st save; find_string to findstring
c
c***********************************************************************
      CHARACTER*1024 FUNCTION SFETCH(name_variable,default)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      CHARACTER*(*) name_variable, default
      CHARACTER*1024 GETSTRING
      INTEGER line, FINDSTRING
      EXTERNAL FINDSTRING, GETSTRING
c***
      line = FINDSTRING(name_variable)
      if (line .eq. 0) then
         SFETCH = default
      else 
         SFETCH = GETSTRING(line)
      endif 
c***
c----      write (*,*) SFETCH
c***
      return 
      end
