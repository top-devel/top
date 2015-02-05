c***********************************************************************
      INTEGER FUNCTION IFETCH(name_variable,default)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      INTEGER  default, GETINTEGER
      CHARACTER*(*) name_variable
      INTEGER line, FINDSTRING
      EXTERNAL GETINTEGER, FINDSTRING
c***
      line = FINDSTRING(name_variable)
      if (line .eq. 0) then
         IFETCH = default
       else
         IFETCH = GETINTEGER(line)
      endif
c***
c---       write (*,*) IFETCH
c***
      return
      end
