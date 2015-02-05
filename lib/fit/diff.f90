       program diff

       implicit none
       integer i
       double precision x, y1, y2
       character*(3) st

       open(unit=2,file="delme_pcspadr",status="old")
       open(unit=3,file="delme_biruni",status="old")
       open(unit=23,file="delme_diff",status="unknown")
       do i=1,251
         read(2,*) x, y1
         read(3,*) x, y2
         write(23,*) x,y1-y2
       enddo
       read(2,*) st
       read(3,*) st
       write(23,*) st
       do i=1,1001
         read(2,*) x, y1
         read(3,*) x, y2
         write(23,*) x,y1-y2
       enddo
       read(2,*) st
       read(3,*) st
       write(23,*) st
       close(23)
       close(3)
       close(2)

       end program
