      function readc2r(line2)
      real readc2r
      character*(*) line2
c use fort.99:
      rewind 99
      write(99,*)line2
      rewind 99
      read(99,*)readc2r
      rewind 99
      return
      end
