      function readc2s(line2)
      character*(*) readc2s
      character*(*) line2
c use fort.99:
      rewind 99
      write(99,*)'''',line2,''''
      rewind 99
      read(99,*)readc2s
      rewind 99
      return
      end
