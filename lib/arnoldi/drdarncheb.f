      subroutine drdarncheb(n,iarn,deg,iparam,tol,
     &                      shift,normA,irc,info,
     &                      intwork,nintwork,rwork,nrwork)

      integer n,iarn,deg,iparam(10),
     &              nintwork,nrwork,
     &              irc(7),info
      integer intwork(nintwork)
      double precision tol,shift,normA,rwork(nrwork)

c
c Local Variables and pointers on the workspace
c

      integer revcom,nrworkmin
      integer pwr, pwi, pu,pv,ph,pz,pwork,py,pwork1,pworkc

c
c Set up the pointers on the workspace
c
      pwr=1
      pwi=pwr+(iarn+1)
      irc(1)=pwi+(iarn+1)
      irc(2)=irc(1)+n
      irc(3)=irc(2)+n
      irc(4)=irc(3)+n
      pu=irc(4)+n
      pv=pu+n
      ph=pv+n*(iarn+1)
      pz=ph+(iarn+1)*(iarn+1)
      pwork=pz+(iarn+1)*(iarn+1)
      py=pwork+(iarn+1)
      pwork1=py+2*(iarn+1)
      pworkc=pwork1+2*(iarn+1)

      nrworkmin=7*n+(7+n)*(iarn+1)+2*(iarn+1)*(iarn+1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Tests on the sizes of the workspaces
c
      if (nintwork.lt.n) then
                         revcom=-101
                         write(*,*)'Increase intwork'
			 info = -17
                         goto 10
      endif
      if (nrwork.lt.nrworkmin) then
                          revcom=-101
                          write(*,*)'Increase rwork'
			  info = -17	
                          goto 10
      endif


      revcom=irc(5)
      call darncheb(n,iarn,deg,iparam,tol,shift,
     &                   normA,rwork(pu),rwork(irc(1)),intwork,
     &                   rwork(pv),rwork(ph),rwork(pz),
     &                   rwork(pwr),rwork(pwi),rwork(pwork),
     &                   rwork(pwork1),rwork(pworkc),
     &                   rwork(py),revcom,info)
      irc(5)=revcom
      irc(6)=pwr
      irc(7)=pwi
      
 10   return
      end
