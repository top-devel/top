      subroutine drzarncheb(n,iarn,deg,iparam,tol,
     &                     shift,normA,irc,info,
     &                     intwork,nintwork,rwork,nrwork,
     &                     zwork, nzwork)
 
 
      integer n,iarn,deg,iparam(10),
     &        nintwork,nzwork,nrwork,
     &        irc(6),info 
      integer intwork(nintwork)
      double precision tol,normA,rwork(nrwork)
      complex *16 shift,zwork(nzwork)
 
c
c Local variables and pointers on the workspace
c
 
      integer revcom, nzworkmin
      integer pu,pv,ph,pz,pw,pwork,pwork1,pworkc,pyc
 
      pw=1
      irc(1)=pw+(iarn+1)
      irc(2)=irc(1)+n
      irc(3)=irc(2)+n
      irc(4)=irc(3)+n
      pu=irc(4)+n
      pv=pu+n
      ph=pv+n*(iarn+1)
      pz=ph+(iarn+1)*(iarn+1)
      pwork=pz+(iarn+1)*(iarn+1)
      pwork1=pwork+(iarn+1)
      pworkc=pwork1+(iarn+1)
      pyc=pworkc+2*n
 
      nzworkmin=6*n+2*(iarn+1)*(iarn+1)+n*(iarn+1)+4*(iarn+1)
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Tests on the size of the workspace
c
      if (nintwork.lt.n) then
                         irc(5)=-101
                         write(*,*)'Increase intwork'
		 	 info = -17	
                         goto 10
      endif
      if (nrwork.lt.iarn) then
                          irc(5)=-101
                          write(*,*)'Increase rwork'
	                  info = -17
                          goto 10
      endif
      if (nzwork.lt.nzworkmin) then
                          irc(5)=-101
                          write(*,*)'Increase zwork'
			  info = -17
                          goto 10
      endif
 
      revcom=irc(5)
      call zarncheb(n,iarn,deg,iparam,tol,shift,
     &                   normA,zwork(pu),zwork(irc(1)),intwork,
     &                   zwork(pv),
     &                   zwork(ph),zwork(pz),zwork(pw),
     &                   zwork(pwork),zwork(pwork1),zwork(pworkc),
     &                   rwork,zwork(pyc),revcom,info)
 
10    irc(5)=revcom
      irc(6)=1
 
      return
      end

