!-------------------------------------------------------------------
! This module provides subroutines for calculating coupling
! integrals.  The subroutines have been optimised using BLAS, and
! easy exits when possible.
!-------------------------------------------------------------------
! List of integrals and associated meanings:
!
!               | {Ylm}* | {dt(Ylm)}* | {Dp(Ylm)}*
!  -----------------------------------------------------------------
!  Yl'm         | Illm   | Jllmc      | i.Kllmc
!  dt(Yl'm)     | Jllm   | Lllm       | i.Mllmc
!  Dp(Yl'm)     | i.Kllm | i.Mllm     | Nllm
!  dt(dt(Yl'm)) | Ollm   | Tllm       | i.Xllm
!  dt(Dp(Yl'm)) | i.Qllm | i.Ullm     | Yllm
!
!  where:
!    Ylm    = Y_l^m spherical harmonic
!    Yl'm   = Y_{l'}^m spherical harmonic
!    l      = harmonic degree (1st index in results from integrals)
!    l'     = harmonic degree (2nd index in results from integrals)
!    m      = azimuthal order
!    dt(f)  = df/d(theta)
!    Dp(f)  = [df/d(phi)]/sin(theta) = im.f/sin(theta)
!    {f}*   = complex conjugate of f
!    i      = sqrt(-1)
!    dOmega = sin(theta)d(theta)d(phi)
!
!  Examples:
!    * G(r,theta) is assumed to be real
!    *   Jllm(G)  = integral_(4pi) G.dt(Yl'm).{Ylm}* dOmega
!    * i.Kllm(G)  = integral_(4pi) G.dP(Yl'm).{Ylm}* dOmega
!    * i.Mllmc(G) = integral_(4pi) G.dt(Yl'm).{Dp(Ylm)}* dOmega
!
!  Note: the last example contains {Dp(Ylm)}*.  Since, Dp(Ylm) is
!        purely imaginary, this introduces a minus sign, since
!        {Dp(Ylm)}* is the complex conjugate. Hence, the subroutines
!        Kllmc, Mllmc and Xllm have a minus sign which intervenes.
!-------------------------------------------------------------------
      module integrales

      use fast_pylm
    
      implicit none
      integer, private, save :: mynt,mynr,mym,mylres,mylmax
      double precision,allocatable,dimension(:,:),save,private :: &
                paux, pylm, pylm2, ptmp, faux, avg_mat
      double precision,allocatable,dimension(:),save,private :: z,w,snt
      integer, allocatable, save, private :: lndx(:)
      integer, private, save :: avg_lb, avg_ub
      public Illm,Jllm,Kllm,Lllm,Mllm,Nllm,Jllmc,Kllmc
      public Mllmc,Ollm,Qllm,Tllm,Ullm,Xllm,Yllm
      public Illmbc,Jllmbc,Kllmbc,Lllmbc,Mllmbc,Nllmbc,Jllmcbc,Kllmcbc
      public Mllmcbc,Ollmbc,Qllmbc,Tllmbc,Ullmbc,Xllmbc,Yllmbc
      public Illm_a,Jllm_a,Kllm_a,Lllm_a,Mllm_a,Nllm_a,Jllmc_a,Kllmc_a
      public Mllmc_a,Ollm_a,Qllm_a,Tllm_a,Ullm_a,Xllm_a,Yllm_a
      public initialisation
contains       

!-------------------------------------------------------------------
!                           SUBROUTINE initialization
!-------------------------------------------------------------------
      subroutine initialisation(nt,nr,m,lres,lmax)
      
      implicit none
      integer, intent(in) :: nt,nr,m,lres,lmax

      mynt   = nt
      mynr   = nr
      mym    = m
      mylres = lres
      mylmax = lmax

      if (allocated(lndx))  deallocate(lndx)
      if (allocated(pylm2)) deallocate(pylm2)
      if (allocated(pylm))  deallocate(pylm)
      if (allocated(ptmp))  deallocate(ptmp)
      if (allocated(paux))  deallocate(paux)
      if (allocated(snt))   deallocate(snt)
      if (allocated(z))     deallocate(z)
      if (allocated(w))     deallocate(w)
 
      allocate(pylm(mylres,mynt),pylm2(mylres,mynt),z(mylres), &
               w(mylres),snt(mylres),lndx(0:mylmax),           &
               ptmp(mylres,mynt),paux(mynt,mynt))
      call gauleg(-1d0,1d0,z,w,mylres)
      w = w * (3.141592653589793d0*2d0)
      snt = sqrt(1d0-z*z)

      end subroutine

!-------------------------------------------------------------------
!                           SUBROUTINE init_avg
!-------------------------------------------------------------------
      subroutine init_avg(mat,lb,ub)

      implicit none
      double precision, intent(in) :: mat(mynr,mynr)
      integer, intent(in) :: lb, ub
      integer i, j

      avg_lb = lb
      avg_ub = ub
      if (allocated(faux))    deallocate(faux)
      if (allocated(avg_mat)) deallocate(avg_mat)
      allocate(faux(mynr,mylres),avg_mat(mynr,mynr))
      do i=1,mynr
        do j=1,mynr
          avg_mat(i,j) = mat(i,j)
        enddo
      enddo
      end subroutine init_avg

!-------------------------------------------------------------------
!                           SUBROUTINE integrales_clear
!-------------------------------------------------------------------
! This subroutine clears auxiliary arrays to free up memory.
!-------------------------------------------------------------------
      subroutine integrales_clear()
      implicit none
      if (allocated(avg_mat)) deallocate(avg_mat)
      if (allocated(pylm2))   deallocate(pylm2)
      if (allocated(pylm))    deallocate(pylm)
      if (allocated(ptmp))    deallocate(ptmp)
      if (allocated(paux))    deallocate(paux)
      if (allocated(faux))    deallocate(faux)
      if (allocated(lndx))    deallocate(lndx)
      if (allocated(snt))     deallocate(snt)
      if (allocated(z))       deallocate(z)
      if (allocated(w))       deallocate(w)
      end subroutine integrales_clear

!-------------------------------------------------------------------
!                           SUBROUTINE avg
!-------------------------------------------------------------------
! This subroutine takes an input array fin(mynr,mylres) and 
! calculates an output array fout(mynr,mylres) which corresponds to fin
! interpolated onto a new grid.  In some cases, this reduces to a simple
! average.
!-------------------------------------------------------------------
! NOTE:  In general, fout(mynr,mylres) = 0
!        (thanks to the definition of avg_mat)
!-------------------------------------------------------------------
      subroutine avg(fin,fout)

      implicit none
      double precision, intent(in)  :: fin(mynr,mylres)
      double precision, intent(out) :: fout(mynr,mylres)
      integer i,ii,j

      fout = 0d0
      do j=1,mylres
        do i=1,mynr
          do ii=max(i-avg_lb,1),min(i+avg_ub,mynr)
            fout(i,j)  = fout(i,j) + avg_mat(i,ii)*fin(ii,j)
          enddo
        enddo
      enddo

      end subroutine avg

!-------------------------------------------------------------------
!                           SUBROUTINE find_coef
!-------------------------------------------------------------------
      subroutine find_coeff(f,tf)

      implicit none
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      double precision temp
      integer i,l,ll,k

      do i=1,mynr
        do k=1,mylres
          temp = f(i,k)*w(k)
          do l=1,mynt
            ptmp(k,l) = pylm2(k,l)*temp
          enddo
        enddo
        call DGEMM("T","N",mynt,mynt,mylres,1d0,ptmp,mylres, &
                   pylm,mylres,0d0,paux,mynt)
        tf(i,:,:) = paux
      enddo

      end subroutine find_coeff

!-------------------------------------------------------------------
!                           SUBROUTINE find_coef_bc
!-------------------------------------------------------------------
      subroutine find_coeff_bc(f,tf)

      implicit none
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      double precision temp, temp2
      integer l,ll,k

      do k=1,mylres
        temp = f(k)*w(k)
        do l=1,mynt
          ptmp(k,l) = pylm2(k,l)*temp
        enddo
      enddo
      call DGEMM("T","N",mynt,mynt,mylres,1d0,ptmp,mylres, &
                 pylm,mylres,0d0,tf,mynt)

      end subroutine find_coeff_bc

!-------------------------------------------------------------------
!                           SUBROUTINE find_lndx
!-------------------------------------------------------------------
      subroutine find_lndx(larray)

      implicit none
      integer, intent(in) :: larray(mynt)
      integer j

      lndx = -1 ! initialisation
      do j=1,mynt
        lndx(larray(j)) = j
      enddo
      end subroutine find_lndx

!-------------------------------------------------------------------
!                           SUBROUTINE Illm
!-------------------------------------------------------------------
      subroutine Illm(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      call find_lndx(lvar)
      call mat_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Illm

!-------------------------------------------------------------------
!                           SUBROUTINE Jllm
!-------------------------------------------------------------------
      subroutine Jllm(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      call find_lndx(lvar)
      call mat_dth_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Jllm

!-------------------------------------------------------------------
!                           SUBROUTINE Kllm
!-------------------------------------------------------------------
      subroutine Kllm(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dphi_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Kllm

!-------------------------------------------------------------------
!                           SUBROUTINE Lllm
!-------------------------------------------------------------------
      subroutine Lllm(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      call find_lndx(lvar)
      call mat_dth_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dth_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Lllm

!-------------------------------------------------------------------
!                           SUBROUTINE Mllm
!-------------------------------------------------------------------
      subroutine Mllm(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dphi_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dth_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Mllm

!-------------------------------------------------------------------
!                           SUBROUTINE Nllm
!-------------------------------------------------------------------
      subroutine Nllm(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dphi_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dphi_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Nllm
        
!-------------------------------------------------------------------
!                           SUBROUTINE Jllmc
!-------------------------------------------------------------------
      subroutine Jllmc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      call find_lndx(lvar)
      call mat_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dth_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Jllmc

!-------------------------------------------------------------------
!                           SUBROUTINE Kllmc
!-------------------------------------------------------------------
      subroutine Kllmc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dphi_ylm(-1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Kllmc

!-------------------------------------------------------------------
!                           SUBROUTINE Mllmc
!-------------------------------------------------------------------
      subroutine Mllmc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dth_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dphi_ylm(-1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Mllmc

!-------------------------------------------------------------------
!                           SUBROUTINE Ollm
!-------------------------------------------------------------------
      subroutine Ollm(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      call find_lndx(lvar)
      call mat_dtt_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Ollm

!-------------------------------------------------------------------
!                           SUBROUTINE Qllm
!-------------------------------------------------------------------
      subroutine Qllm(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dtp_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Qllm

!-------------------------------------------------------------------
!                           SUBROUTINE Tllm
!-------------------------------------------------------------------
      subroutine Tllm(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      call find_lndx(lvar)
      call mat_dtt_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dth_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Tllm

!-------------------------------------------------------------------
!                           SUBROUTINE Ullm
!-------------------------------------------------------------------
      subroutine Ullm(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dtp_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dth_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Ullm

!-------------------------------------------------------------------
!                           SUBROUTINE Xllm
!-------------------------------------------------------------------
      subroutine Xllm(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dtt_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dphi_ylm(-1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Xllm

!-------------------------------------------------------------------
!                           SUBROUTINE Yllm
!-------------------------------------------------------------------
      subroutine Yllm(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dtp_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dphi_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff(f,tf)
      end subroutine Yllm

!-------------------------------------------------------------------
!                           SUBROUTINE Illmbc
!-------------------------------------------------------------------
      subroutine Illmbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      call find_lndx(lvar)
      call mat_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Illmbc

!-------------------------------------------------------------------
!                           SUBROUTINE Jllmbc
!-------------------------------------------------------------------
      subroutine Jllmbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      call find_lndx(lvar)
      call mat_dth_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Jllmbc

!-------------------------------------------------------------------
!                           SUBROUTINE Kllmbc
!-------------------------------------------------------------------
      subroutine Kllmbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dphi_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Kllmbc

!-------------------------------------------------------------------
!                           SUBROUTINE Lllmbc
!-------------------------------------------------------------------
      subroutine Lllmbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      call find_lndx(lvar)
      call mat_dth_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dth_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Lllmbc

!-------------------------------------------------------------------
!                           SUBROUTINE Mllmbc
!-------------------------------------------------------------------
      subroutine Mllmbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dphi_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dth_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Mllmbc

!-------------------------------------------------------------------
!                           SUBROUTINE Nllmbc
!-------------------------------------------------------------------
      subroutine Nllmbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dphi_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dphi_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Nllmbc
        
!-------------------------------------------------------------------
!                           SUBROUTINE Jllmcbc
!-------------------------------------------------------------------
      subroutine Jllmcbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      call find_lndx(lvar)
      call mat_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dth_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Jllmcbc

!-------------------------------------------------------------------
!                           SUBROUTINE Kllmcbc
!-------------------------------------------------------------------
      subroutine Kllmcbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dphi_ylm(-1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Kllmcbc

!-------------------------------------------------------------------
!                           SUBROUTINE Mllmcbc
!-------------------------------------------------------------------
      subroutine Mllmcbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dth_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dphi_ylm(-1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Mllmcbc

!-------------------------------------------------------------------
!                           SUBROUTINE Ollmbc
!-------------------------------------------------------------------
      subroutine Ollmbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      call find_lndx(lvar)
      call mat_dtt_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Ollmbc

!-------------------------------------------------------------------
!                           SUBROUTINE Qllmbc
!-------------------------------------------------------------------
      subroutine Qllmbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dtp_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Qllmbc

!-------------------------------------------------------------------
!                           SUBROUTINE Tllmbc
!-------------------------------------------------------------------
      subroutine Tllmbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      call find_lndx(lvar)
      call mat_dtt_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dth_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Tllmbc

!-------------------------------------------------------------------
!                           SUBROUTINE Ullmbc
!-------------------------------------------------------------------
      subroutine Ullmbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dtp_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dth_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Ullmbc

!-------------------------------------------------------------------
!                           SUBROUTINE Xllmbc
!-------------------------------------------------------------------
      subroutine Xllmbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dtt_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dphi_ylm(-1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Xllmbc

!-------------------------------------------------------------------
!                           SUBROUTINE Yllmbc
!-------------------------------------------------------------------
      subroutine Yllmbc(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mylres)
      double precision, intent(out) :: tf(mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call find_lndx(lvar)
      call mat_dtp_ylm(1d0,pylm,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_lndx(leq)
      call mat_dphi_ylm(1d0,pylm2,mylres,mynt,z,snt,lndx,mylmax,mym)
      call find_coeff_bc(f,tf)
      end subroutine Yllmbc

!-------------------------------------------------------------------
!                           SUBROUTINE Illm_a
!-------------------------------------------------------------------
      subroutine Illm_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      call avg(f,faux)
      call Illm(faux,tf,leq,lvar)
      end subroutine Illm_a
!-------------------------------------------------------------------
!                           SUBROUTINE Jllm_a
!-------------------------------------------------------------------
      subroutine Jllm_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      call avg(f,faux)
      call Jllm(faux,tf,leq,lvar)
      end subroutine Jllm_a
!-------------------------------------------------------------------
!                           SUBROUTINE Kllm_a
!-------------------------------------------------------------------
      subroutine Kllm_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call avg(f,faux)
      call Kllm(faux,tf,leq,lvar)
      end subroutine Kllm_a
!-------------------------------------------------------------------
!                           SUBROUTINE Lllm_a
!-------------------------------------------------------------------
      subroutine Lllm_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      call avg(f,faux)
      call Lllm(faux,tf,leq,lvar)
      end subroutine Lllm_a
!-------------------------------------------------------------------
!                           SUBROUTINE Mllm_a
!-------------------------------------------------------------------
      subroutine Mllm_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call avg(f,faux)
      call Mllm(faux,tf,leq,lvar)
      end subroutine Mllm_a
!-------------------------------------------------------------------
!                           SUBROUTINE Nllm_a
!-------------------------------------------------------------------
      subroutine Nllm_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call avg(f,faux)
      call Nllm(faux,tf,leq,lvar)
      end subroutine Nllm_a
!-------------------------------------------------------------------
!                           SUBROUTINE Jllmc_a
!-------------------------------------------------------------------
      subroutine Jllmc_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      call avg(f,faux)
      call Jllmc(faux,tf,leq,lvar)
      end subroutine Jllmc_a
!-------------------------------------------------------------------
!                           SUBROUTINE Kllmc_a
!-------------------------------------------------------------------
      subroutine Kllmc_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call avg(f,faux)
      call Kllmc(faux,tf,leq,lvar)
      end subroutine Kllmc_a
!-------------------------------------------------------------------
!                           SUBROUTINE Mllmc_a
!-------------------------------------------------------------------
      subroutine Mllmc_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call avg(f,faux)
      call Mllmc(faux,tf,leq,lvar)
      end subroutine Mllmc_a
!-------------------------------------------------------------------
!                           SUBROUTINE Ollm_a
!-------------------------------------------------------------------
      subroutine Ollm_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      call avg(f,faux)
      call Ollm(faux,tf,leq,lvar)
      end subroutine Ollm_a
!-------------------------------------------------------------------
!                           SUBROUTINE Qllm_a
!-------------------------------------------------------------------
      subroutine Qllm_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call avg(f,faux)
      call Qllm(faux,tf,leq,lvar)
      end subroutine Qllm_a
!-------------------------------------------------------------------
!                           SUBROUTINE Tllm_a
!-------------------------------------------------------------------
      subroutine Tllm_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      call avg(f,faux)
      call Tllm(faux,tf,leq,lvar)
      end subroutine Tllm_a
!-------------------------------------------------------------------
!                           SUBROUTINE Ullm_a
!-------------------------------------------------------------------
      subroutine Ullm_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call avg(f,faux)
      call Ullm(faux,tf,leq,lvar)
      end subroutine Ullm_a
!-------------------------------------------------------------------
!                           SUBROUTINE Xllm_a
!-------------------------------------------------------------------
      subroutine Xllm_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call avg(f,faux)
      call Xllm(faux,tf,leq,lvar)
      end subroutine Xllm_a
!-------------------------------------------------------------------
!                           SUBROUTINE Yllm_a
!-------------------------------------------------------------------
      subroutine Yllm_a(f,tf,leq,lvar)
      implicit none
      integer, intent(in) :: leq(mynt), lvar(mynt)
      double precision, intent(in)  :: f(mynr,mylres)
      double precision, intent(out) :: tf(mynr,mynt,mynt)
      if (mym.eq.0) then
        tf = 0d0
        return
      endif
      call avg(f,faux)
      call Yllm(faux,tf,leq,lvar)
      end subroutine Yllm_a
!-------------------------------------------------------------------
      end module integrales
