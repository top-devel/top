      CHARACTER*1024 FUNCTION GETSTRING(i_line)
c***
      INCLUDE 'undefined.h'
      INCLUDE 'GP_common.h'
      INTEGER i_line
      INTEGER ibegin, i, iend
      CHARACTER*1024 line1,line2
      CHARACTER*1 eq_char, blank
      DATA eq_char/'='/, blank/' '/
      INTRINSIC INDEX




c***
      line1 = lines_of_text(i_line)
      ibegin = 1 + INDEX(line1,eq_char)
      do i = 1, 1024
         line2(i:i) = blank
      enddo
      line2(1:1024 - ibegin + 1) = line1(ibegin:1024)



      read(line2,'(A)') GETSTRING

c***
c***
      return
      end
