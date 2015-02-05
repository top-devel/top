      function readc2i(line2)
      integer readc2i
      character*(*) line2
c use fort.99:
      rewind 99
      write(99,*)line2
      rewind 99
      read(99,*)readc2i
      rewind 99
      return
      end
