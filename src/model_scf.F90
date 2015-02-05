#include "config.h"
       module model
       use grid
       use inputs
       ! works with MacGregor et al. models (2D)

       ! raw variables
       integer, save :: nr_mod, ntheta_mod
       double precision, save :: mass1, luminosity1, radius1, eta1,alfa1
       double precision, save :: mass2, luminosity2, radius2, x, y
       double precision, allocatable, dimension(:), save :: xx, qq, bb, &
              T1D, rho1D, p1D, beta, log_op, n_plus_one
       character*(4), allocatable, save :: conv(:)
       double precision, save :: luminosity3, radius3, tc, rhoc, pc
       double precision, save :: mass3, Log_L, Log_Te, Log_g, Re,   &
                        Log_Pc, Log_Tc, Log_rho, &
                        T_over_W, Ve, xcc, qcc, xce, qce
       double precision, allocatable, dimension(:), save :: theta, gravity, temp
       double precision, save :: mass4, alfa2, eta2, Log_L2, Req,   &
                        Rp_over_Req, Z_over_Req, Teff, Teq, Tp, Veq
       double precision, allocatable, save ::  xk(:,:)

       ! variables needed for pulsation calculations
       double precision, allocatable, dimension(:,:), save :: rhom, &
                   rhom_z, rhom_t, pm, pm_z, pm_t, c2, Gamma1, &
                   Rota, Rota_prime
       double precision, save :: Lambda, Rota_0

       ! geometric terms
       double precision,allocatable,dimension(:,:),save::r_t,r_z,r_map,&
                        re_t,re_z,re_map,r_zz,r_zt,r_tt,re_zz,re_zt,&
                        re_tt,zeta,cost,sint,cott,r_aux
       double precision,allocatable,dimension(:),save::cth,sth,r
       double precision,save::K,epsilon
       integer, save :: nr
        

       ! different physical constants
       double precision, parameter :: solar_mass = 1.9891d33    !g
       double precision, parameter :: solar_radius = 6.95508d10 !cm
       double precision, parameter :: G = 6.672d-8  !cm^3.g^-1.s^-2

!----------------------------------------------------------------------------------
! A few comments:
!
!   mass1       = mass (in M_odot)
!   radius1     = radius (in R_odot)
!   mass2       = mass (in g)
!   radius2     = radius (in cm)
!   xx          = tilde(r)/Re
!   qq          = M(xx)/M_total
!   bb          = fractional luminosity
!   T           = temperature in K
!   rho         = density (in g/cm^3)
!   p           = pressure (in g/(cm.s^2) = dyne/cm^2)
!   beta        = Rosseland mean opacity (in cm^2/g)
!   log op      = ???
!   n+1         = effective polytropic index (essentially dlnT/dlnp) ?
!----------------------------------------------------------------------------------


contains

!----------------------------------------------------------------------------------
!  This subroutine calls all the necessary subroutines to set up the a model
!  from MacGregor et al.
!----------------------------------------------------------------------------------
       subroutine init_model()

       call read_model()
       call make_mapping()
       call init_fields()
       call init_Rota()
       !call write_fields() ! This is only for debugging.
       !                    ! It will make the program stop.

       end subroutine

!----------------------------------------------------------------------------------
!  This subroutine reads the raw variables.
!----------------------------------------------------------------------------------
       subroutine read_model()
       implicit none

       character*(512) oneline
       integer i, j

       nr_mod = 251
       ntheta_mod = 35

       if (allocated(xx))         deallocate(xx)
       if (allocated(qq))         deallocate(qq)
       if (allocated(bb))         deallocate(bb)
       if (allocated(T1D))        deallocate(T1D)
       if (allocated(rho1D))      deallocate(rho1D)
       if (allocated(p1D))        deallocate(p1D)
       if (allocated(beta))       deallocate(beta)
       if (allocated(log_op))     deallocate(log_op)
       if (allocated(n_plus_one)) deallocate(n_plus_one)
       if (allocated(conv))       deallocate(conv)
       if (allocated(theta))      deallocate(theta)
       if (allocated(gravity))    deallocate(gravity)
       if (allocated(temp))       deallocate(temp)
       if (allocated(xk))         deallocate(xk)

       allocate(xx(nr_mod), qq(nr_mod), bb(nr_mod), T1D(nr_mod), rho1D(nr_mod), p1D(nr_mod), &
                beta(nr_mod), log_op(nr_mod), n_plus_one(nr_mod), conv(nr_mod))
       allocate(theta(0:ntheta_mod), gravity(0:ntheta_mod), temp(0:ntheta_mod))
       allocate(xk(nr_mod,2*ntheta_mod-1))

      
       open(unit=2,file=trim(filename),status="old")

       call read_next_line(2,oneline)
       call make_numeric(oneline)
       read(oneline,*) mass1,luminosity1,radius1,eta1,alfa1

       call read_next_line(2,oneline)
       call make_numeric(oneline)
       read(oneline,*) mass2,luminosity2,radius2,x,y

       call read_next_line(2,oneline)
       read(2,*) (xx(i),qq(i),bb(i),T1D(i),rho1D(i),p1D(i), &
                  beta(i),log_op(i),n_plus_one(i),conv(i),i=1,nr_mod)

       call read_next_line(2,oneline)
       call read_next_line(2,oneline)
       call make_numeric(oneline)
       read(oneline,*) luminosity3, radius3, tc, rhoc, pc

       call read_next_line(2,oneline)
       call read_next_line(2,oneline)
       read(oneline,*) mass3, Log_L, Log_Te, Log_g, Re, Log_Pc, Log_Tc, Log_rho, &
                       T_over_W, Ve, xcc, qcc, xce, qce

       call read_next_line(2,oneline)
       read(2,*) (theta(i), gravity(i), temp(i),i=(1+ntheta_mod)/2,0,-1)
       
       call read_next_line(2,oneline)
       call read_next_line(2,oneline)
       read(2,*) mass4, alfa2, eta2, Log_L2, Req, Rp_over_Req, Z_over_Req, &
                 Teff, Teq, Tp, Veq

       call read_next_line(2,oneline)
       read(2,*) ((xk(i,j),j=(1+ntheta_mod)/2,1,-1),i=1,nr_mod)
       
       close(2)
       
       ! This constructs the rest of the arrays theta(:), gravity(:), temp(:)
       ! and xk(:,:) by equatorial symmetry.
       do j=1,ntheta_mod/2
         theta(ntheta_mod-j+1)   = 180d0-theta(j)
         gravity(ntheta_mod-j+1) = gravity(j)
         temp(ntheta_mod-j+1)    = temp(j)
         do i = 1,nr_mod
           xk(i,ntheta_mod-j+1)  = xk(i,j)
         enddo
       enddo

       ! initialise some contants
       mass  = mass1
       eta   = eta1
       alpha = alfa1

       end subroutine

!--------------------------------------------------------------------------
!  This prepares the mapping r(i,j) and other geometrical functions.
!  It is necessary to interpolate in the angular direction - this is done
!  by spectral methods involving Legendre polynomials.
!--------------------------------------------------------------------------
        subroutine make_mapping()
        
        use mod_legendre

        implicit none
        double precision, allocatable, dimension(:) :: &
          rs_spec, rs_ll, rs, rsp, rss, w
        double precision, allocatable :: rs_aux(:,:)
        double precision xi
        integer i,j,l
        
        if (lres.lt.ntheta_mod)  &
                stop 'Case lres < ntheta_mod not implemented'

        nr = grd(1)%nr

        if (allocated(zeta))         deallocate(zeta)
        if (allocated(sint))         deallocate(sint)
        if (allocated(cost))         deallocate(cost)
        if (allocated(cott))         deallocate(cott)
        if (allocated(cth))          deallocate(cth)
        if (allocated(sth))          deallocate(sth)
        if (allocated(r_map))        deallocate(r_map)
        if (allocated(r_t))          deallocate(r_t)
        if (allocated(r_z))          deallocate(r_z)
        if (allocated(r_zz))         deallocate(r_zz)
        if (allocated(r_zt))         deallocate(r_zt)
        if (allocated(r_tt))         deallocate(r_tt)
        if (allocated(re_map))       deallocate(re_map)
        if (allocated(re_t))         deallocate(re_t)
        if (allocated(re_z))         deallocate(re_z)
        if (allocated(re_zz))        deallocate(re_zz)
        if (allocated(re_zt))        deallocate(re_zt)
        if (allocated(re_tt))        deallocate(re_tt)
        if (allocated(r_aux))        deallocate(r_aux)

        allocate(rs_spec(lres),rs(lres),rsp(lres),  &
                 rss(lres),w(lres),zeta(nr,lres),   &
                 sint(nr,lres),cost(nr,lres),       &
                 cott(nr,lres),cth(lres),sth(lres), &
                 r_map(nr,lres),r_t(nr,lres),       &
                 r_z(nr,lres),re_map(nr,lres),      &
                 re_t(nr,lres),re_z(nr,lres),       &
                 r_zz(nr,lres),r_zt(nr,lres),       &
                 r_tt(nr,lres),re_zz(nr,lres),      &
                 re_zt(nr,lres),re_tt(nr,lres),     &
                 rs_ll(lres),r_aux(nr_mod,lres),    &
                 rs_aux(nr_mod,lres))

        rs_spec = 0d0;rs = 0d0;rsp = 0d0;
        rss = 0d0;w = 0d0;zeta = 0d0;
        sint = 0d0;cost = 0d0;
        cott = 0d0;cth = 0d0;sth = 0d0;
        r_map = 0d0;r_t = 0d0;
        r_z = 0d0;re_map = 0d0;
        re_t = 0d0;re_z = 0d0;
        r_zz = 0d0;r_zt = 0d0;
        r_tt = 0d0;re_zz = 0d0;
        re_zt = 0d0;re_tt = 0d0;
        rs_ll = 0d0;rs_aux = 0d0;

        grd(2)%nr = grd(1)%nr
        allocate(grd(1)%r(grd(1)%nr))
        allocate(grd(2)%r(grd(2)%nr))

        allocate(r(nr))
        do i=1,nr
          r(i) = dble(i-1)/dble(nr-1)
          zeta(i,:) = r(i)
          grd(1)%r(i) = r(i)
          grd(2)%r(i) = r(i)
        enddo

        call gauleg(-1d0,1d0,cth,w,lres)
        sth = sqrt(1-cth**2)
        do i=1,nr
          sint(i,:) = sth(:)
          cost(i,:) = cth(:)
          cott(i,:) = cth(:)/sth(:)
        enddo

        ! real-spectral transformation
        call legendre(xk(nr_mod,1:ntheta_mod),rs_spec(1:ntheta_mod), &
                      ntheta_mod,0,1)
        do l = 0,ntheta_mod-1
          rs_ll(l+1) = -dble(l*(l+1))*rs_spec(l+1)
        enddo

        ! spectral-real transformation
        call legendre(rs_spec(1:lres),rs(1:lres),lres,0,-1) 
        call legendrep(rs_spec(1:lres),rsp(1:lres),lres,0)
        call legendre(rs_ll(1:lres),rss(1:lres),lres,0,-1)
        rss = rss - cott(1,:)*rsp

        ! Find the inner domain:
        K = eval_ylm(rs_spec(1:ntheta_mod),ntheta_mod,dcos(0d0))
        epsilon = 1d0-K
        
        do i=1,nr
          do j=1,lres
            r_map(i,j) = K*r(i) + 0.5d0*(5d0*r(i)**3-3d0*r(i)**5)*(rs(j)-K)
            r_z  (i,j) = K    + 0.5d0*(15d0*r(i)**2-15d0*r(i)**4)*(rs(j)-K)
            r_zz (i,j) =        0.5d0*(30d0*r(i)   -60d0*r(i)**3)*(rs(j)-K)
            r_t  (i,j) =        0.5d0*(5d0*r(i)**3-3d0*r(i)**5)*rsp(j)
            r_zt (i,j) =        0.5d0*(15d0*r(i)**2-15d0*r(i)**4)*rsp(j)
            r_tt (i,j) =        0.5d0*(5d0*r(i)**3-3d0*r(i)**5)*rss(j)
          enddo
        enddo

        ! real-spectral transformation
        call legendre(xk(1:nr_mod,1:ntheta_mod),     &
                      rs_aux(1:nr_mod,1:ntheta_mod), &
                      ntheta_mod,nr_mod,0,1)
        call legendre(rs_aux(1:nr_mod,1:lres),       &
                      r_aux(1:nr_mod,1:lres),        &
                      lres,nr_mod,0,-1)

        ! Find the outer domain (we fold this domain onto the inner domain):
        ! Note: we're derive with respect to r(i), not xi, hence the sign change for
        ! re_z and re_zt.

        do i=1,nr
          xi = 2d0-r(i)
          do j=1,lres
            re_map(i,j) = 2d0*epsilon + K*xi + (2d0*xi**3-9d0*xi**2+12d0*xi-4d0)*(rs(j)-1d0-epsilon)
            re_z(i,j)   =             - K    - (6d0*xi**2-18d0*xi+12d0)*(rs(j)-1d0-epsilon)
            re_zz(i,j)  =                      (12d0*xi-18d0)*(rs(j)-1d0-epsilon)
            re_t(i,j)   =                      (2d0*xi**3-9d0*xi**2+12d0*xi-4d0)*rsp(j)
            re_zt(i,j)  =                    - (6d0*xi**2-18d0*xi+12d0)*rsp(j)
            re_tt(i,j)  =                      (2d0*xi**3-9d0*xi**2+12d0*xi-4d0)*rss(j)
          enddo
        enddo

        deallocate(rs_spec,rs_ll,rs,rsp,rss,w,rs_aux)

        end subroutine

!--------------------------------------------------------------
! This function evaluates a function decomposed over the
! spherical harmonic basis at a given z = cos(theta)
!--------------------------------------------------------------
        double precision function eval_ylm(f,nn,z)
        
        implicit none
        integer nn
        double precision f(1:nn), z
        double precision, parameter :: pi = 3.14159265358979d0
        double precision yl1, yl2, yl3
        integer l
        
        yl1 = sqrt(0.25d0/pi)
        yl2 = sqrt(0.75d0/pi)*z
        eval_ylm = yl1*f(1)+yl2*f(2)
        ! Beware: f(l+1) corresponds to f at l
        do l = 2,nn-1
          yl3 = (sqrt(dble(4*l*l-1))/dble(l))*z*yl2  &
              + (1d0/dble(l)-1d0)*sqrt(dble(2*l+1)/(2*l-3))*yl1
          eval_ylm = eval_ylm + yl3*f(l+1)
          yl1 = yl2
          yl2 = yl3
        enddo
        return
        end function
!--------------------------------------------------------------
!  This subroutine prepares different non-dimensional variables
!  for the pulsation equations.
!--------------------------------------------------------------
         subroutine init_fields()

         implicit none
         integer i
         double precision, parameter :: gamma = 5d0/3d0
         double precision, parameter :: pi = 3.14159265358979d0
         double precision, allocatable, dimension(:) :: &
                rho_aux, p_aux, c2_aux, Gamma1_aux
         double precision, allocatable :: aux(:,:)

         if (allocated(rhom))   deallocate(rhom)
         if (allocated(rhom_z)) deallocate(rhom_z)
         if (allocated(rhom_t)) deallocate(rhom_t)
         if (allocated(pm))     deallocate(pm)
         if (allocated(pm_z))   deallocate(pm_z)
         if (allocated(pm_t))   deallocate(pm_t)
         if (allocated(c2))     deallocate(c2)
         if (allocated(Gamma1)) deallocate(Gamma1)

         allocate(rhom(nr,lres), rhom_z(nr,lres), rhom_t(nr,lres), &
                  pm(nr,lres), pm_z(nr,lres), pm_t(nr,lres),       &
                  c2(nr,lres), Gamma1(nr,lres))
         
         allocate(rho_aux(nr_mod), p_aux(nr_mod),     &
                  c2_aux(nr_mod), Gamma1_aux(nr_mod))

         Lambda = 4d0*pi*G*radius2**2*rho1D(1)**2/p1D(1)
         do i = 1, nr_mod
           rho_aux(i)    = rho1D(i)/rho1D(1)
           p_aux(i)      = p1D(i)/p1D(1)
           Gamma1_aux(i) = 1d0 + (4d0-3d0*beta(i))*(gamma-1d0)/    &
                           (beta(i)**2 + 3d0*(gamma-1d0)*          &
                           (1d0-beta(i))*(4d0+beta(i)))
           c2_aux(i)     = Gamma1_aux(i)* p1D(i)/rho1D(i)*rho1D(1)/p1D(1)
         enddo
         
         allocate(aux(nr,lres))
         call map2D(Gamma1_aux, Gamma1)
         call map2D(c2_aux, c2)
         call map2D_der(rho_aux, rhom, rhom_t, rhom_z, aux)
         call map2D_der(p_aux, pm, pm_t, pm_z, aux)
         
         deallocate(rho_aux, p_aux, Gamma1_aux, c2_aux, aux)

         end subroutine

!--------------------------------------------------------------------------
!  This subroutine maps a field from the iso-potentials grid to a
!  new mapping, which introduces dependance of theta.
!--------------------------------------------------------------------------
         subroutine map2D(f1D,f2D)

         use mod_legendre

         implicit none
         double precision f1D(nr_mod), f2D(nr,lres)
         integer i, j

         ! make use of equatorial symmetry
         do j=1,(lres+1)/2
           call interpolate(r_aux(1:nr_mod,j),f1D,nr_mod,r_map(1:nr,j),f2D(1:nr,j),nr)
           f2D(1:nr,lres+1-j) = f2D(1:nr,j)
         enddo
         
         end subroutine

!--------------------------------------------------------------------------
!  This subroutine maps a field from the iso-potentials grid to a
!  new mapping, which introduces dependance of theta.  It also calculates
!  the theta and zeta derivatives on the new mapping.
!--------------------------------------------------------------------------
         subroutine map2D_der(f1D,f2D,f2D_t,f2D_z,aux)

         use mod_legendre

         implicit none
         double precision f1D(nr_mod), aux(nr,lres)
         double precision f2D(nr,lres), f2D_t(nr,lres), f2D_z(nr,lres)
         integer i, j

         ! make use of equatorial symmetry
         do j=1,(lres+1)/2
           call interpolate(r_aux(1:nr_mod,j),f1D,nr_mod,r_map(1:nr,j),f2D(1:nr,j),nr)
           call interpolate_derive(r_aux(1:nr_mod,j),f1D,nr_mod,r_map(1:nr,j),f2D_z(1:nr,j),nr)
           f2D(1:nr,lres+1-j)   = f2D(1:nr,j)
           f2D_z(1:nr,lres+1-j) = f2D_z(1:nr,j)
         enddo

         aux = 0d0
         call legendre(f2D(1:nr,1:lres),aux(1:nr,1:lres),lres,nr,0,1)
         call legendrep(aux(1:nr,1:lres),f2D_t(1:nr,1:lres),lres,nr,0)

         end subroutine

!--------------------------------------------------------------------------
!  This prepares the rotation profile.  This subroutine must be called
!  after the mapping has been initialised.
!--------------------------------------------------------------------------
       subroutine init_Rota()

       implicit none
       integer i,j

       if (alfa1.ge.0d0) then
         Rota_0 = (1d0 + alfa1**2)*(Ve*1d5)*sqrt(rho1D(1)/p1D(1))
       else
         Rota_0 = (Ve*1d5)*sqrt(rho1D(1)/p1D(1))/(1d0+alfa1**2)
       endif

       if (allocated(Rota))       deallocate(Rota)
       if (allocated(Rota_prime)) deallocate(Rota_prime)

       allocate (Rota(nr,lres), Rota_prime(nr,lres))

       if (alfa1.ge.0d0) then
         do i=1,nr
           do j=1,lres
             Rota(i,j) = Rota_0/(1d0+(alfa1*r_map(i,j)*sint(i,j))**2)
             Rota_prime(i,j)= -2d0*Rota_0*alfa1**2*r_map(i,j)*sint(i,j)/ &
                              (1d0+(alfa1*r_map(i,j)*sint(i,j))**2)**2
           enddo
         enddo
       else
         do i=1,nr
           do j=1,lres
             Rota(i,j) = Rota_0*(1d0+(alfa1*r_map(i,j)*sint(i,j))**2)
             Rota_prime(i,j)= 2d0*Rota_0*alfa1**2*r_map(i,j)*sint(i,j)
           enddo
         enddo
       endif

       end subroutine
!-------------------------------------------------------------
!  This subroutine prints different fields so as to help debug
!  the program...
!-------------------------------------------------------------
       subroutine write_fields()
       implicit none
       integer i,j

       open(unit=999,file="delme_grid",status="unknown")
       write(999,*) nr,lres
       write(999,'(a)') "# x"
       write(999,101) ((r_map(i,j)*sint(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# y"
       write(999,101) ((r_map(i,j)*cost(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# xe"
       write(999,101) ((re_map(i,j)*sint(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# ye"
       write(999,101) ((re_map(i,j)*cost(i,j),i=1,nr),j=1,lres)
       close(999)

       open(unit=999,file="delme_fields",status="unknown")
       write(999,*) nr,lres
       write(999,'(a)') "# rhom"
       write(999,101) ((rhom(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# rhom_z"
       write(999,101) ((rhom_z(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# rhom_t"
       write(999,101) ((rhom_t(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# pm"
       write(999,101) ((pm(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# pm_z"
       write(999,101) ((pm_z(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# pm_t"
       write(999,101) ((pm_t(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# c2"
       write(999,101) ((c2(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# Gamma1"
       write(999,101) ((Gamma1(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# Rota"
       write(999,101) ((Rota(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# Rota_prime"
       write(999,101) ((Rota_prime(i,j),i=1,nr),j=1,lres)
       close(999)

       open(unit=999,file="delme_geometry",status="unknown")
       write(999,*) nr,lres
       write(999,'(a)') "# r_map"
       write(999,101) ((r_map(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# r_z"
       write(999,101) ((r_z(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# r_t"
       write(999,101) ((r_t(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# r_zz"
       write(999,101) ((r_zz(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# r_zt"
       write(999,101) ((r_zt(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# r_tt"
       write(999,101) ((r_tt(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# re_map"
       write(999,101) ((re_map(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# re_z"
       write(999,101) ((re_z(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# re_t"
       write(999,101) ((re_t(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# re_zz"
       write(999,101) ((re_zz(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# re_zt"
       write(999,101) ((re_zt(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# re_tt"
       write(999,101) ((re_tt(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# zeta"
       write(999,101) ((zeta(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# cost"
       write(999,101) ((cost(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# sint"
       write(999,101) ((sint(i,j),i=1,nr),j=1,lres)
       write(999,'(a)') "# cott"
       write(999,101) ((cott(i,j),i=1,nr),j=1,lres)
       close(999)

       open(unit=999,file="delme_constants",status="unknown")
       write(999,'(a)') "# Lambda"
       write(999,'(1pe22.15)') Lambda
       close(999)

101    format(3(1pe22.15,2X))
       stop
       end subroutine

!-------------------------------------------------------------
!  This subroutine determines whether a line is blank
!-------------------------------------------------------------
       subroutine isblank(string,bool)
       implicit none
       character*(*) string
       character*(*), parameter :: blanks = ' '
       logical bool
       integer i
       
       bool = .true.
       do i=1,len(string)
         if (index(blanks,string(i:i)) == 0) then
           bool = .false.
           return
         endif
       enddo
       return
       end subroutine

!-------------------------------------------------------------
!  This subroutine reads the next non-blank, non-commentary
!  line.
!-------------------------------------------------------------
       subroutine read_next_line(iu,oneline)
       implicit none
       integer iu
       character*(512) oneline
       logical bool

       do
         read(iu,'(a512)') oneline
         if (oneline(1:1) .eq. "=") cycle
         call isblank(oneline,bool)
         if (bool) cycle
         return
       enddo
       end subroutine
!-------------------------------------------------------------------
!  This subroutine filters a string so as to only leave
!  numeric characters and 'D' or 'E' (in the middle of a string). 
!-------------------------------------------------------------------
       subroutine make_numeric(string)
       implicit none
       character*(*) string
       character*(*), parameter :: numerics = '0123456789+-.DE '
       integer i
       
       do i=1,len(string)
         if (index(numerics,string(i:i)) == 0) string(i:i) = ' '
       enddo
       do i=1,len(string)-1
         if ((string(i:i+1) == ' D').or.    &
             (string(i:i+1) == ' E').or.    &
             (string(i:i+1) == 'D ').or.    &
             (string(i:i+1) == 'E ')) then
             string(i:i+1) = "  "
         endif
       enddo
       return
       end subroutine

       end module
