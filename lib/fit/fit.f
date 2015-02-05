c project the star field y from 'star' on the 
c Chebyshev grid. Derive it if ider=1.

        subroutine fit(np,x,y,sf,k,s,ider,r,yinter,dyinter)

	implicit double precision(a-h,o-z)

        common/dim/ nz

	parameter (npm=2000,km=5,nest=npm+km+1)
	dimension x(npm),y(npm),w(npm),
     &            t(nest),c(nest),d(km+1),
     &            wrk(npm*(km+1)+nest*(7+3*km))
	integer iwrk(nest)
	dimension r(nz+1),yinter(nz+1),dyinter(nz+1)
	character*6 sf

c we set up the weights of the data points
	do i=1,npm
	 w(i)=1.d0
	enddo

c we set up the boundaries of the approximation interval
	xb=x(1)
	xe=x(np)
	
c we set up the dimension information
 	lwrk=npm*(km+1)+nest*(7+3*km)
	iopt=0

	call curfit(iopt,np,x,y,w,xb,xe,k,s,nest,n,
     &              t,c,fp,wrk,lwrk,iwrk,ier)

	rap=float(n)/float(np)
 	if ( ier .le. 0 ) write(6,5) n,sf,rap
   5	format('   # ',i4,' knots for ',a6,' then knots/ncouche= ',f5.2)

c spline evaluation of y on the Chebyschev grid
        call splev(t,n,c,k,r,yinter,nz+1,ier)

c derive the field with spalde if ider=1
	if ( ider .eq. 1 ) then
	 write(6,6) sf
   6	 format('   # derive ',a6)
         k1=k+1
         do i=1,nz+1
          x0=r(i)
          call spalde(t,n,c,k1,x0,d,ier)
          dyinter(i)=d(2)
         enddo
	endif

	return
	end
