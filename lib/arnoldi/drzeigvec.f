      subroutine drzeigvec(n,iarn,zwork,nzwork,
     &                     zworkvec,nzworkvec,
     &                     rwork,nrwork,select,iparam)
 
      integer n, iarn, iparam(10), nzwork,
     &        nrwork, nzworkvec
      logical select(iarn+1)
      double precision rwork(nrwork)
      complex *16 zwork(nzwork), zworkvec(nzworkvec)
 
c
c Local variables and pointers on the workspace
c
 
      integer col1, col2, col3, col4
      integer pu,pv,ph,pz,pw,pwork,pwork1,pworkc,pyc
      integer pvec,pwork2,pz1,pz2
      integer nzworkvecmin
 
      pw=1
      col1=pw+(iarn+1)
      col2=col1+n
      col3=col2+n
      col4=col3+n
      pu=col4+n
      pv=pu+n
      ph=pv+n*(iarn+1)
      pz=ph+(iarn+1)*(iarn+1)
      pwork=pw+(iarn+1)*(iarn+1)
      pwork1=pwork+(iarn+1)
      pworkc=pwork1+(iarn+1)
      pyc=pworkc+2*n
 
      pvec=1
      pwork2=pvec+n*(iarn+1)
      pz1=pwork2+3*(iarn+1)
      pz2=pz1+(iarn+1)*(iarn+1)
 
      nzworkvecmin=n*(iarn+1)+3*(iarn+1)+2*(iarn+1)*(iarn+1)
      if (nzworkvec.lt.nzworkvecmin) then
       write(*,*)'There is not enough space to compute the eigenvectors'
       goto 10
      endif
 
      call zeigvec(n,iarn,zwork(ph),zwork(pz),
     &             zwork(pv),zworkvec(pvec),zwork(pw),
     &             zworkvec(pwork2),rwork,select,zworkvec(pz1),
     &             zworkvec(pz2),iparam)
 
10    return
      end

