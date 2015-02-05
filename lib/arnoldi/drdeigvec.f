      subroutine drdeigvec(n,iarn,rwork,nrwork,
     &                     rworkvec, nrworkvec,
     &                     select,iparam)

      integer n, iarn, iparam(10), nrwork,
     &        nrworkvec
      logical select(iarn+1)
      double precision rwork(nrwork)
      double precision rworkvec(nrworkvec)
c
c Local Variables and pointers on the workspace
c

      integer pu,pv,ph,pz,pwr,pwi,pwork,py,pwork1,pworkc,
     &        pvec, pwork2, pz1, pz2,
     &        col1, col2, col3, col4, nrworkvecmin

      nrworkvecmin=2*iarn**2+n*iarn+5*iarn+4
      if (nrworkvec.lt.nrworkvecmin) then
       write(*,*)'*****************************************************'
       write(*,*)'There is not enough space to compute the eigenvectors'
       write(*,*)'*****************************************************'
       goto 10
      endif
     

c
c Set up the pointers on the workspace
c

      pwr=1
      pwi=pwr+(iarn+1)
      col1=pwi+(iarn+1)
      col2=col1+n
      col3=col2+n
      col4=col3+n
      pu=col4+n
      pv=pu+n
      ph=pv+n*(iarn+1)
      pz=ph+(iarn+1)*(iarn+1)
      pwork=pwi+(iarn+1)*(iarn+1)
      py=pwork+(iarn+1)
      pwork1=py+2*(iarn+1)
      pworkc=pwork1+2*(iarn+1)

      pvec=1
      pwork2=pvec+n*(iarn+1)
      pz1=pwork2+3*(iarn+1)
      pz2=pz1+(iarn+1)*(iarn+1)

      call deigvec(n, iarn, rwork(ph), rwork(pz),
     &             rwork(pv), rworkvec(pvec),
     &             rwork(pwi), rworkvec(pwork2),
     &             select, rworkvec(pz1), rworkvec(pz2), iparam)

10    return
      end 
       
