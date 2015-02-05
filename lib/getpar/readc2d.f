      function readc2d(line2)
      double precision readc2d
      character*(*) line2
c use fort.99:
      rewind 99
      write(99,*)line2
      rewind 99
      read(99,*)readc2d
      rewind 99
      return
      end
