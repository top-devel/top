c***********************************************************************
      REAL FUNCTION RFETCH(name_variable,default)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      REAL  default, GETREAL
      CHARACTER*(*) name_variable
      INTEGER line, FINDSTRING
      EXTERNAL GETREAL, FINDSTRING
c***
      line = FINDSTRING(name_variable)
      if (line .eq. 0) then
         RFETCH = default
       else
         RFETCH = GETREAL(line)
      endif
c***
c---       write (*,*) RFETCH
c***
      return
      end
