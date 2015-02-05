       program test_legendre

       use fast_legendre

       implicit none
       integer, parameter :: n = 100, ns = 50
       integer, parameter :: lstart = 0, lincr = 2, m = 0
       double precision, parameter :: alpha = 1d0
       double precision, dimension(n) :: cth, sth, w, f, f2, df, df2
       double precision, dimension(ns):: fs
       integer i

       call gauleg(-1d0,1d0,cth,w,n)
       sth = sqrt(1d0-cth*cth)
       do i=1,n
         f(i) = 1d0
       enddo
       fs = 0d0
       call project_ylm(alpha,f,cth,w,fs,n,ns,lstart,lincr,m)

       open(unit=3,file="plotme",status="unknown")
       do i=1,ns
         write(3,*) 2*(i-1),fs(i)
       enddo
       close(3)

       end program test_legendre
