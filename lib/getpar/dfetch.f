c***********************************************************************
      DOUBLE PRECISION FUNCTION DFETCH(name_variable,default)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      DOUBLE PRECISION  default, GETDOUBLE
      CHARACTER*(*) name_variable
      INTEGER line, FINDSTRING
      EXTERNAL GETDOUBLE, FINDSTRING
c***
      line = FINDSTRING(name_variable)
      if (line .eq. 0) then
         DFETCH = default
       else
         DFETCH = GETDOUBLE(line)
      endif

      return
      end
