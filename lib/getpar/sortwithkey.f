c*****************************************************************
      SUBROUTINE sortwithkey(n,ra,key)
c***
c     Sort array ra of length n into decreasing size order
c     using the Heapsort algorithm, and sort key as well
c***
      INTEGER n
      INTEGER ra(n)
      INTEGER rra
      CHARACTER*1024 key(n)
      CHARACTER*1024 kkey
      INTEGER i,j,k,ir
c***
c     The sort algorithm has a run-time error if n .eq. 1 
c***
      if (n .le. 1) then
         return
      endif
c***
c     k is decremented during heap creation
c***
      k = (n/2) + 1
c***
c     ir is decremented during heap selection
c***
      ir = n
10    continue
      if (k .gt. 1) then
c***  
c     in heap creation phase
c***  
         k = k - 1
         rra = ra(k)
         kkey = key(k)
      else
c***  
c     in heap selection phase
c***  
         rra = ra(ir)
         kkey = key(ir)
         ra(ir) = ra(1)
         key(ir) = key(1)
         ir = ir - 1
         if (ir .eq. 1) then
            ra(1) = rra
            key(1) = kkey
            return
         endif
      endif
      i = k
      j = k + k
      do while (j .le. ir)
         if (j .lt. ir) then
            if (ra(j) .gt. ra(j+1)) then
               j = j + 1
            endif
         endif
         if (rra .gt. ra(j)) then
            ra(i) = ra(j)
            key(i) = key(j)
            i = j
            j = j + j
         else
            j = ir + 1
         endif
      enddo
      ra(i) = rra
      key(i) = kkey
      go to 10
c***  
      end
