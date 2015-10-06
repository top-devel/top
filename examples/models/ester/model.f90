      module model

      use inputs
      use mod_grid


!-------------------------------------------------------------------------
! Model file structure
! 
! nrm: Number of radial points
! nthm: Number of angular points
! ndom: Number of domains
! nconv: Number of convective domains from the center of the star
! npts_i: Number of radial points in the i-th domain
! R_i(th): External boundary of the i-th domain
! 
! Xc gives the Hydrogen abundance in the convective core X(core)=Xc*X
! The units are cgs
! Rotation is in radians/sec
! 
!-------------------------------------------------------------------------

      ! domain dependant variables
      type STAR_DOMAIN

        ! outer surface of domain
        double precision, pointer, dimension(:)   :: Rs

        ! raw variables (of dimension grd(id)%nr x nthm)
        double precision, pointer, dimension(:,:) ::                &
                       rr, p_raw, rho_raw, Gamma1_raw, Rota_raw,    &
                       p_z_raw, rho_z_raw, Rota_z_raw

        ! final variables (of dimension grd(id)%nr x lres)
        double precision, pointer, dimension(:,:) ::                &
                       pm, rhom, Gamma1, pm_z, pm_t, rhom_z, rhom_t,&
                       c2, NNtoz, Rota, Rota_z, Rota_t, pm_zz,      &
                       pm_zt, pm_tz, pm_tt, pm_ez, grd_pe_z,        &
                       grd_pe_t, grd_pe_zz, grd_pe_zt, grd_pe_tz,   &
                       grd_pe_tt, det_rhom_pm
        double precision, pointer, dimension(:,:) :: r_t,r_z,r_map, &
                       r_zz,r_zt,r_tt,zeta,cost,sint,cott,rrt,roz

      end type STAR_DOMAIN

      integer, save :: nrm, nthm, ndom, nconv, npts_max
      double precision, save :: mass_input, radius, luminosity, rotation
      double precision, save :: X, Z, Xc, Req, Rp
      double precision, save :: rho_ref, p_ref, t_ref
      double precision, allocatable, save :: theta(:),zeta(:),sth(:),cth(:)
      type(STAR_DOMAIN), allocatable, save :: s(:)
      logical, save :: first_model = .true.

      ! different physical constants
      double precision, parameter :: pi = 3.141592653589793d0
      double precision, parameter :: Lambda = 4d0*pi
      double precision, parameter :: solar_mass   = 1.9891d33  !g
      double precision, parameter :: solar_radius = 6.95508d10 !cm
      double precision, parameter :: G = 6.672d-8 !cm^3.g^-1.s^-2

contains

!-------------------------------------------------------------------------
! This suboutine calls all of the relevant subroutines for initialising
! a model.
!-------------------------------------------------------------------------
      subroutine init_model()

      implicit none
      integer id,j

      call read_model()
      call interpolate_model()
      call make_mapping()
      call find_NNtoz()
      call find_grd_pm()
      call find_grd_pe()
      !call write_fields() ! this stops the program

      end subroutine init_model

!-------------------------------------------------------------------------
! This suboutine reads the model
!-------------------------------------------------------------------------
      subroutine read_model()

      use mod_legendre

      implicit none
      integer i, j, id
      double precision, allocatable :: r_spec(:), r_temp(:)
      double precision :: aux

      if (.not.first_model) call clear_model()
      first_model = .false.

      open(unit=37,file=trim(filename),status="old")
      read(37,*) nrm, nthm, ndom, nconv
      ndom = ndom + 1 ! add an extra outer domain
      if (ndom.ne.ndomains) stop "Incorrect number of domains"
      allocate(s(ndom), theta(nthm), zeta(nrm))
      read(37,*) (grd(id)%nr,id=1,ndom-1)
      npts_max = 0
      do id=1,ndom
        if (npts_max.lt.grd(id)%nr) npts_max = grd(id)%nr
      enddo
      call allocate_model()
      read(37,*)
      read(37,*) mass_input, radius, luminosity, rotation, X, Z, Xc
      read(37,*)
      read(37,*) (theta(j),j=1,nthm)
      read(37,*)
      read(37,*) (zeta(i),i=1,nrm)
      read(37,*)
      do id=1,ndom-1
        read(37,*) (s(id)%Rs(j),j=1,nthm)
      enddo
      read(37,*)
      do id=1,ndom-1
        read(37,*) ((s(id)%rr(i,j),j=1,nthm),i=1,grd(id)%nr)
      enddo
      read(37,*)
      do id=1,ndom-1
        read(37,*) ((s(id)%p_raw(i,j),j=1,nthm),i=1,grd(id)%nr)
      enddo
      read(37,*)
      do id=1,ndom-1
        read(37,*) ((s(id)%rho_raw(i,j),j=1,nthm),i=1,grd(id)%nr)
      enddo
      read(37,*)
      do id=1,ndom-1
        read(37,*) ((s(id)%Gamma1_raw(i,j),j=1,nthm),i=1,grd(id)%nr)
      enddo
      read(37,*)
      do id=1,ndom-1
        read(37,*) ((s(id)%Rota_raw(i,j),j=1,nthm),i=1,grd(id)%nr)
      enddo
      close(37)

      ! Find equatorial and polar radii
      allocate(r_temp(2*nthm),r_spec(2*nthm))
      do j=1,nthm
        r_temp(nthm+j)   = s(ndom-1)%Rs(j)
        r_temp(nthm+1-j) = s(ndom-1)%Rs(j)
      enddo
      r_spec = 0d0
      call legendre(r_temp,r_spec,2*nthm,0,1)
      Rp  = eval_ylm(r_spec,2*nthm,1d0)
      Req = eval_ylm(r_spec,2*nthm,0d0)
      deallocate(r_temp,r_spec)

      ! sanity check on mapping
      do id=2,ndom-1
        aux = s(id)%Rs(nthm)-s(id-1)%Rs(nthm)
        do j=1,nthm-1
          if ((s(id)%Rs(j)-s(id-1)%Rs(j)).lt.aux) then
            print*,"Problem with make_mapping."
            stop
          endif
        enddo
      enddo

      ! Find refence quantities
      rho_ref = mass_input/Req**3 
      p_ref   = G*mass_input**2/Req**4 
      t_ref   = sqrt(Req**3/(G*mass_input))
      mass    = mass_input/solar_mass
      rota_avg= rotation*t_ref
      !print*,rotation*Req*1d-5 ! = veq

      end subroutine read_model

!-------------------------------------------------------------------------
! This allocates all the necessary arrays for the model
!-------------------------------------------------------------------------
      subroutine allocate_model()

      implicit none
      integer id

      do id=1,ndom
        allocate(grd(id)%r(grd(id)%nr))
      enddo
      do id=1,ndom-1
        allocate(s(id)%p_raw(grd(id)%nr,nthm))
        allocate(s(id)%p_z_raw(grd(id)%nr,nthm))
        allocate(s(id)%rho_raw(grd(id)%nr,nthm))
        allocate(s(id)%rho_z_raw(grd(id)%nr,nthm))
        allocate(s(id)%Gamma1_raw(grd(id)%nr,nthm))
        allocate(s(id)%Rota_raw(grd(id)%nr,nthm))
        allocate(s(id)%Rota_z_raw(grd(id)%nr,nthm))
        allocate(s(id)%pm(grd(id)%nr,lres))
        allocate(s(id)%rhom(grd(id)%nr,lres))
        allocate(s(id)%Gamma1(grd(id)%nr,lres))
        allocate(s(id)%Rota(grd(id)%nr,lres))
        allocate(s(id)%c2(grd(id)%nr,lres))
        allocate(s(id)%pm_z(grd(id)%nr,lres))
        allocate(s(id)%pm_t(grd(id)%nr,lres))
        allocate(s(id)%rhom_z(grd(id)%nr,lres))
        allocate(s(id)%rhom_t(grd(id)%nr,lres))
        allocate(s(id)%Rota_z(grd(id)%nr,lres))
        allocate(s(id)%Rota_t(grd(id)%nr,lres))
        allocate(s(id)%NNtoz(grd(id)%nr,lres))
        allocate(s(id)%pm_zz(grd(id)%nr,lres))
        allocate(s(id)%pm_tz(grd(id)%nr,lres))
        allocate(s(id)%pm_zt(grd(id)%nr,lres))
        allocate(s(id)%pm_tt(grd(id)%nr,lres))
        allocate(s(id)%pm_ez(grd(id)%nr,lres))
        allocate(s(id)%grd_pe_z(grd(id)%nr,lres))
        allocate(s(id)%grd_pe_t(grd(id)%nr,lres))
        allocate(s(id)%grd_pe_zz(grd(id)%nr,lres))
        allocate(s(id)%grd_pe_zt(grd(id)%nr,lres))
        allocate(s(id)%grd_pe_tz(grd(id)%nr,lres))
        allocate(s(id)%grd_pe_tt(grd(id)%nr,lres))
        allocate(s(id)%det_rhom_pm(grd(id)%nr,lres))
        s(id)%NNtoz = 0d0
        s(id)%pm_zz = 0d0
        s(id)%pm_zt = 0d0
        s(id)%pm_tz = 0d0
        s(id)%pm_tt = 0d0
        s(id)%pm_ez = 0d0
        s(id)%grd_pe_z = 0d0
        s(id)%grd_pe_t = 0d0
        s(id)%grd_pe_zz = 0d0
        s(id)%grd_pe_zt = 0d0
        s(id)%grd_pe_tz = 0d0
        s(id)%grd_pe_tt = 0d0
        s(id)%det_rhom_pm = 0d0
      enddo
      do id=1,ndom
        allocate(s(id)%Rs(nthm))
        allocate(s(id)%rr(grd(id)%nr,nthm))
        allocate(s(id)%r_map(grd(id)%nr,lres))
        allocate(s(id)%r_z(grd(id)%nr,lres))
        allocate(s(id)%r_t(grd(id)%nr,lres))
        allocate(s(id)%r_zz(grd(id)%nr,lres))
        allocate(s(id)%r_zt(grd(id)%nr,lres))
        allocate(s(id)%r_tt(grd(id)%nr,lres))
        allocate(s(id)%zeta(grd(id)%nr,lres))
        allocate(s(id)%cost(grd(id)%nr,lres))
        allocate(s(id)%sint(grd(id)%nr,lres))
        allocate(s(id)%cott(grd(id)%nr,lres))
        allocate(s(id)%roz(grd(id)%nr,lres))
        allocate(s(id)%rrt(grd(id)%nr,lres))
        s(id)%roz = 0d0
        s(id)%rrt = 0d0
      enddo
      end subroutine allocate_model

!-------------------------------------------------------------------------
! This deallocates the model to make space for a new model
!-------------------------------------------------------------------------
      subroutine clear_model()

      implicit none
      integer id

      do id=1,ndom
        deallocate(grd(id)%r)
      enddo

      do id=1,ndom-1
        deallocate(s(id)%p_raw)
        deallocate(s(id)%p_z_raw)
        deallocate(s(id)%rho_raw)
        deallocate(s(id)%rho_z_raw)
        deallocate(s(id)%Gamma1_raw)
        deallocate(s(id)%Rota_raw)
        deallocate(s(id)%Rota_z_raw)
        deallocate(s(id)%pm)
        deallocate(s(id)%rhom)
        deallocate(s(id)%Gamma1)
        deallocate(s(id)%Rota)
        deallocate(s(id)%c2)
        deallocate(s(id)%pm_z)
        deallocate(s(id)%pm_t)
        deallocate(s(id)%rhom_z)
        deallocate(s(id)%rhom_t)
        deallocate(s(id)%Rota_z)
        deallocate(s(id)%Rota_t)
        deallocate(s(id)%NNtoz)
        deallocate(s(id)%pm_zz)
        deallocate(s(id)%pm_tz)
        deallocate(s(id)%pm_zt)
        deallocate(s(id)%pm_tt)
        deallocate(s(id)%pm_ez)
        deallocate(s(id)%grd_pe_z)
        deallocate(s(id)%grd_pe_t)
        deallocate(s(id)%grd_pe_zz)
        deallocate(s(id)%grd_pe_zt)
        deallocate(s(id)%grd_pe_tz)
        deallocate(s(id)%grd_pe_tt)
        deallocate(s(id)%det_rhom_pm)
      enddo

      do id=1,ndom
        deallocate(s(id)%Rs)
        deallocate(s(id)%rr)
        deallocate(s(id)%r_map)
        deallocate(s(id)%r_z)
        deallocate(s(id)%r_t)
        deallocate(s(id)%r_zz)
        deallocate(s(id)%r_zt)
        deallocate(s(id)%r_tt)
        deallocate(s(id)%zeta)
        deallocate(s(id)%cost)
        deallocate(s(id)%sint)
        deallocate(s(id)%cott)
        deallocate(s(id)%roz)
        deallocate(s(id)%rrt)
      enddo
      deallocate(s, theta, zeta)
      end subroutine clear_model

!-------------------------------------------------------------------------
!  This interpolates the model to a higher angular resolution
!-------------------------------------------------------------------------
      subroutine interpolate_model()

      use mod_legendre
      use derivative
      
      implicit none
      integer i, j, id, npid, nrs, nrf
      double precision, allocatable :: mat_raw(:,:), mat(:,:)
      type(DERMAT) :: dmat

      ! Don't forget to initialise the radial derivative of p, rho and
      ! Rota.
      nrs = 1
      nrf = grd(1)%nr
      do id=1,ndomains-1
        call init_derive_cheb(dmat,zeta(nrs),grd(id)%nr,1,0)
        call dgemm("N","N",grd(id)%nr,nthm,grd(id)%nr,1d0, &
                   dmat%derive(:,:,1),grd(id)%nr,          &
                   s(id)%p_raw,grd(id)%nr,0d0,             &
                   s(id)%p_z_raw,grd(id)%nr)
        call dgemm("N","N",grd(id)%nr,nthm,grd(id)%nr,1d0, &
                   dmat%derive(:,:,1),grd(id)%nr,          &
                   s(id)%rho_raw,grd(id)%nr,0d0,           &
                   s(id)%rho_z_raw,grd(id)%nr)
        call dgemm("N","N",grd(id)%nr,nthm,grd(id)%nr,1d0, &
                   dmat%derive(:,:,1),grd(id)%nr,          &
                   s(id)%Rota_raw,grd(id)%nr,0d0,           &
                   s(id)%Rota_z_raw,grd(id)%nr)
        nrs = nrs + grd(id)%nr
        nrf = nrs + grd(id+1)%nr - 1
        call clear_derive(dmat)
      enddo

      allocate(mat_raw(npts_max,nthm*2),mat(npts_max,lres))
      
      do id=1,ndom-1
        npid = grd(id)%nr
        do j=1,nthm
          mat_raw(1:npid,nthm+j)   = s(id)%p_raw(:,j)/p_ref
          mat_raw(1:npid,nthm+1-j) = s(id)%p_raw(:,j)/p_ref
        enddo
        mat = 0d0
        call legendre(mat_raw(1:npid,:),mat(1:grd(id)%nr,1:(2*nthm)), &
                     (2*nthm),grd(id)%nr,0,1)
        call legendre(mat(1:grd(id)%nr,1:lres),s(id)%pm,lres,         &
                      grd(id)%nr,0,-1)
        call legendrep(mat(1:grd(id)%nr,1:lres),s(id)%pm_t,lres,      &
                      grd(id)%nr,0)

        do j=1,nthm
          mat_raw(1:npid,nthm+j)   = s(id)%rho_raw(:,j)/rho_ref
          mat_raw(1:npid,nthm+1-j) = s(id)%rho_raw(:,j)/rho_ref
        enddo
        mat = 0d0
        call legendre(mat_raw(1:npid,:),mat(1:grd(id)%nr,1:(2*nthm)), &
                     (2*nthm),grd(id)%nr,0,1)
        call legendre(mat(1:grd(id)%nr,1:lres),s(id)%rhom,lres,       &
                      grd(id)%nr,0,-1)
        call legendrep(mat(1:grd(id)%nr,1:lres),s(id)%rhom_t,lres,    &
                      grd(id)%nr,0)

        do j=1,nthm
          mat_raw(1:npid,nthm+j)   = s(id)%Rota_raw(:,j)*t_ref
          mat_raw(1:npid,nthm+1-j) = s(id)%Rota_raw(:,j)*t_ref
        enddo
        mat = 0d0
        call legendre(mat_raw(1:npid,:),mat(1:grd(id)%nr,1:(2*nthm)), &
                     (2*nthm),grd(id)%nr,0,1)
        call legendre(mat(1:grd(id)%nr,1:lres),s(id)%Rota,lres,       &
                      grd(id)%nr,0,-1)
        call legendrep(mat(1:grd(id)%nr,1:lres),s(id)%Rota_t,lres,    &
                      grd(id)%nr,0)

        do j=1,nthm
          mat_raw(1:npid,nthm+j)   = s(id)%Gamma1_raw(:,j)
          mat_raw(1:npid,nthm+1-j) = s(id)%Gamma1_raw(:,j)
        enddo
        mat = 0d0
        call legendre(mat_raw(1:npid,:),mat(1:grd(id)%nr,1:(2*nthm)), &
                     (2*nthm),grd(id)%nr,0,1)
        call legendre(mat(1:grd(id)%nr,1:lres),s(id)%Gamma1,lres, &
                      grd(id)%nr,0,-1)

        s(id)%c2 = s(id)%Gamma1*s(id)%pm/s(id)%rhom

        do j=1,nthm
          mat_raw(1:npid,nthm+j)   = s(id)%p_z_raw(:,j)/p_ref
          mat_raw(1:npid,nthm+1-j) = s(id)%p_z_raw(:,j)/p_ref
        enddo
        mat = 0d0
        call legendre(mat_raw(1:npid,:),mat(1:grd(id)%nr,1:(2*nthm)), &
                     (2*nthm),grd(id)%nr,0,1)
        call legendre(mat(1:grd(id)%nr,1:lres),s(id)%pm_z,lres, &
                      grd(id)%nr,0,-1)

        do j=1,nthm
          mat_raw(1:npid,nthm+j)   = s(id)%rho_z_raw(:,j)/rho_ref
          mat_raw(1:npid,nthm+1-j) = s(id)%rho_z_raw(:,j)/rho_ref
        enddo
        mat = 0d0
        call legendre(mat_raw(1:npid,:),mat(1:grd(id)%nr,1:(2*nthm)), &
                     (2*nthm),grd(id)%nr,0,1)
        call legendre(mat(1:grd(id)%nr,1:lres),s(id)%rhom_z,lres, &
                      grd(id)%nr,0,-1)

        do j=1,nthm
          mat_raw(1:npid,nthm+j)   = s(id)%Rota_z_raw(:,j)*t_ref
          mat_raw(1:npid,nthm+1-j) = s(id)%Rota_z_raw(:,j)*t_ref
        enddo
        mat = 0d0
        call legendre(mat_raw(1:npid,:),mat(1:grd(id)%nr,1:(2*nthm)), &
                     (2*nthm),grd(id)%nr,0,1)
        call legendre(mat(1:grd(id)%nr,1:lres),s(id)%Rota_z,lres, &
                      grd(id)%nr,0,-1)

      enddo
      deallocate(mat, mat_raw)
      end subroutine interpolate_model


!-------------------------------------------------------------------------
!  This calculates the array NNtoz: it needs to be called after the
!  mapping is created.
!-------------------------------------------------------------------------
      subroutine find_NNtoz()

      implicit none
      integer id, i, j
      
      id = 1
      do j=1,lres
        do i=2,grd(id)%nr
          s(id)%NNtoz(i,j) = (-s(id)%c2(i,j)*s(id)%rhom_t(i,j)  &
                              +s(id)%pm_t(i,j))/s(id)%zeta(i,j)
        enddo
        s(id)%NNtoz(1,j) = 0d0
      enddo
      end subroutine find_NNtoz

!-------------------------------------------------------------------------
! This finds extra gradients of the pressure, used for an alternate
! boundary condition.
!-------------------------------------------------------------------------
      subroutine find_grd_pm()

      use mod_legendre
      use derivative
      
      implicit none
      integer i, j, id, npid, nrs, nrf
      double precision, allocatable :: mat(:,:), smat(:,:)
      type(DERMAT) :: dmat

      ! find various indices
      nrs = 1
      do id=1,ndomains-2
        nrs = nrs + grd(id)%nr
      enddo
      id = ndomains - 1
      nrf = nrs + grd(id)%nr - 1

      ! initialise mat and dmat
      call init_derive_cheb(dmat,zeta(nrs),grd(id)%nr,1,0)
      allocate(mat(grd(id)%nr,lres))
      allocate(smat(grd(id)%nr,lres))

      ! find pm_zz
      mat = s(id)%zeta**2*s(id)%pm_z/(s(id)%r_map**2*s(id)%r_z)
      call dgemm("N","N",grd(id)%nr,lres,grd(id)%nr,1d0, &
                 dmat%derive(:,:,1),grd(id)%nr,          &
                 mat,grd(id)%nr,0d0,             &
                 s(id)%pm_zz,grd(id)%nr)

      ! find pm_tz
      call legendre(mat,smat,lres,grd(id)%nr,0,1)
      call legendrep(smat,s(id)%pm_tz,lres,grd(id)%nr,0)

      ! find pm_zt
      mat = s(id)%zeta*s(id)%pm_t/(s(id)%r_map**2*s(id)%r_z)
      call dgemm("N","N",grd(id)%nr,lres,grd(id)%nr,1d0, &
                 dmat%derive(:,:,1),grd(id)%nr,          &
                 mat,grd(id)%nr,0d0,             &
                 s(id)%pm_zt,grd(id)%nr)

      ! find pm_tt: (quick and ugly solution for d(sin t P(cos t))/dt)
      mat = mat/s(id)%sint
      s(id)%pm_tt = s(id)%cost*mat
      call legendre(mat,smat,lres,grd(id)%nr,0,1)
      call legendrep(smat,mat,lres,grd(id)%nr,0)
      s(id)%pm_tt = s(id)%pm_tt + (s(id)%cost**2-1d0)*mat

      ! find pm_ez (just a convenient quantity):
      s(id)%pm_ez = (s(id)%r_map**2+s(id)%r_t**2)*s(id)%pm_z  &
                  / (s(id)%r_map**2*s(id)%r_z**2)             &
                  -  s(id)%r_t*s(id)%pm_t                     &
                  / (s(id)%r_map**2*s(id)%r_z)

      deallocate(mat,smat)
      call clear_derive(dmat)

      end subroutine find_grd_pm

!--------------------------------------------------------------------------
! This subroutine initialises all of the grd_pe arrays in such a way
! that there are no singularities in the center of the star.
!--------------------------------------------------------------------------
       subroutine find_grd_pe()

       use derivative
       use mod_legendre

       implicit none
       type(DERMAT) :: dm
       double precision, allocatable :: aux(:,:)
       integer id, j, i, ii


       do id=1,ndomains-1
         call init_derive_cheb(dm,grd(id)%r,grd(id)%nr,1,0)
       
         s(id)%det_rhom_pm = s(id)%zeta*(s(id)%rhom_z*s(id)%pm_t &
                        - s(id)%pm_z*s(id)%rhom_t)               &
                        / (s(id)%r_map**2*s(id)%r_z*s(id)%rhom**2)
         s(id)%grd_pe_z = -s(id)%zeta**2*s(id)%pm_z &
                        / (s(id)%r_map**2*s(id)%r_z*s(id)%rhom)
         s(id)%grd_pe_t = -s(id)%zeta*s(id)%pm_t    &
                        / (s(id)%r_map**2*s(id)%r_z*s(id)%rhom)

         ! This has to be done prior to calculating numerical
         ! derivatives:
         if (id.eq.1) then
           s(id)%det_rhom_pm(1,:) = 0d0
           s(id)%grd_pe_z(1,:) =-s(id)%pm_z(1,:)/(s(id)%roz(1,:)**2 &
                               * s(id)%r_z(1,:)*s(id)%rhom(1,:))
           s(id)%grd_pe_t(1,:) = 0d0
         endif

         ! Numerically calculate zeta derivatives: this avoids various
         ! difficulties with singularities in the center of the star:
         call dgemm("N","N",grd(id)%nr,lres,grd(id)%nr,1d0, &
                   dm%derive(:,:,1),grd(id)%nr,             &
                   s(id)%grd_pe_z,grd(id)%nr,0d0,           &
                   s(id)%grd_pe_zz,grd(id)%nr)
         call dgemm("N","N",grd(id)%nr,lres,grd(id)%nr,1d0, &
                   dm%derive(:,:,1),grd(id)%nr,             &
                   s(id)%grd_pe_t,grd(id)%nr,0d0,           &
                   s(id)%grd_pe_zt,grd(id)%nr)

         ! Numerically calculate theta derivatives
         allocate(aux(grd(id)%nr,lres))
         call legendre(s(id)%grd_pe_z,aux,lres,grd(id)%nr,0,1)
         call legendrep(aux,s(id)%grd_pe_tz,lres,grd(id)%nr,0)

         ! This is not the most elegant solution but the legendre grid
         ! does not go to theta = 0, so this should work (at reasonable
         ! resolutions)
         call legendre(s(id)%grd_pe_t/s(id)%sint,aux,lres,grd(id)%nr,0,1)
         call legendrep(aux,s(id)%grd_pe_tt,lres,grd(id)%nr,0)
         s(id)%grd_pe_tt = s(id)%sint*s(id)%grd_pe_tt+s(id)%cott*s(id)%grd_pe_t

         call clear_derive(dm)
         deallocate(aux)

       enddo

       end subroutine find_grd_pe

!-------------------------------------------------------------------------
!  This sets up the mapping for the model
!-------------------------------------------------------------------------
      subroutine make_mapping()

      use mod_legendre

      implicit none
      integer i, j, l, id, npts_sum
      double precision :: xi, zz, eta_p, eta_m, aa, bb
      double precision, allocatable, dimension(:) :: a,ap,as,r_spec,  &
             r_temp,rs,rsp,rss,rs_prev,rsp_prev,rss_prev,w,llr_spec

      if (allocated(cth)) deallocate(cth)
      if (allocated(sth)) deallocate(sth)

      allocate(a(npts_max),ap(npts_max),as(npts_max),r_spec(lres), &
               r_temp(2*nthm), rs(lres),rsp(lres),rss(lres),       &
               rs_prev(lres),rsp_prev(lres),rss_prev(lres),w(lres),&
               llr_spec(lres),cth(lres),sth(lres))

      ! calculation of sin(theta), cos(theta)
      call gauleg(-1d0,1d0,cth,w,lres)
      sth = sqrt(1-cth**2)

      npts_sum = 0
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! inner domain
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      rs_prev = 0d0
      rsp_prev= 0d0
      rss_prev= 0d0
      do j=1,nthm
        r_temp(nthm+j)   = s(1)%Rs(j)/Req
        r_temp(nthm+1-j) = s(1)%Rs(j)/Req
      enddo
      r_spec = 0d0
      call legendre(r_temp,r_spec(1:(2*nthm)),2*nthm,0,1)
      call legendre(r_spec,rs,lres,0,-1)
      call legendrep(r_spec,rsp,lres,0)
      do l=1,lres ! Beware: r_spec(l) is r_spec of (ell-1)
        llr_spec(l)=-dble(l*(l-1))*r_spec(l)
      enddo
      call legendre(llr_spec(1:lres),rss(1:lres),lres,0,-1)
      rss = rss - cth*rsp/sth

      ! find radial variables
      eta_m = zeta(1)
      eta_p = zeta(grd(1)%nr)
      do i=1,grd(1)%nr
        zz = zeta(i)
        grd(1)%r(i) = zz
        s(1)%zeta(i,1:lres) = zz
        xi    = zz/eta_p
        a(i)  = 0.5d0*( 5d0*xi**3- 3d0*xi**5)
        ap(i) = 0.5d0*(15d0*xi**2-15d0*xi**4)/eta_p
        as(i) = 0.5d0*(30d0*xi   -60d0*xi**3)/eta_p**2
        !a(i)  =   6d0*xi**5- 15d0*xi**4+10d0*xi**3
        !ap(i) = (30d0*xi**4- 60d0*xi**3+30d0*xi**2)/eta_p
        !as(i) =(120d0*xi**3-180d0*xi**2+60d0*xi   )/eta_p**2
      enddo

      ! find r and it's derivatives
      aa = Rp/Req
      do i=1,grd(1)%nr
        zz = zeta(i)
        s(1)%r_map(i,:)= aa*zz + a(i)*( rs(:)-aa*eta_p)
        s(1)%r_z(i,:)  = aa    +ap(i)*( rs(:)-aa*eta_p)
        s(1)%r_t(i,:)  =         a(i)* rsp(:)
        s(1)%r_zz(i,:) =        as(i)*( rs(:)-aa*eta_p)
        s(1)%r_zt(i,:) =        ap(i)* rsp(:)
        s(1)%r_tt(i,:) =         a(i)* rss(:)
      enddo
      
      ! find extra geometrical terms
      do i=1,grd(1)%nr
        s(1)%cost(i,:) = cth(:)
        s(1)%sint(i,:) = sth(:)
        s(1)%cott(i,:) = cth(:)/sth(:)
      enddo

      ! define roz and rrt:
      do i=2,grd(1)%nr
        s(1)%roz(i,:) = s(1)%r_map(i,:)/s(1)%zeta(i,:)
        s(1)%rrt(i,:) = s(1)%r_t(i,:)/s(1)%r_map(i,:)
      enddo
      s(1)%roz(1,:) = s(1)%r_z(1,:)
      s(1)%rrt(1,:) = s(1)%r_zt(1,:)/s(1)%r_z(1,:)

      npts_sum = grd(1)%nr
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! intermediate domains
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do id=2,ndom-1

        ! inner bound = previous outer bound
        rs_prev  = rs
        rsp_prev = rsp
        rss_prev = rss

        ! interpolate outer bound
        do j=1,nthm
          r_temp(nthm+j)   = s(id)%Rs(j)/Req
          r_temp(nthm+1-j) = s(id)%Rs(j)/Req
        enddo
        r_spec = 0d0
        call legendre(r_temp,r_spec(1:(2*nthm)),2*nthm,0,1)
        call legendre(r_spec,rs,lres,0,-1)
        call legendrep(r_spec,rsp,lres,0)
        do l=1,lres ! Beware: r_spec(l) is r_spec of (ell-1)
          llr_spec(l)=-dble(l*(l-1))*r_spec(l)
        enddo
        call legendre(llr_spec(1:lres),rss(1:lres),lres,0,-1)
        rss = rss - cth*rsp/sth ! Please don't forget this line!

        ! find radial variables
        eta_m = zeta(npts_sum+1)
        eta_p = zeta(npts_sum+grd(id)%nr)
        do i=1,grd(id)%nr
          zz = zeta(npts_sum+i)
          grd(id)%r(i) = zz
          s(id)%zeta(i,1:lres) = zz
          xi    = (2d0*zz-eta_p-eta_m)/(eta_p-eta_m)
          a(i)  = 0.25d0*(    -xi**3+3d0*xi+2d0)
          ap(i) = 0.25d0*(-3d0*xi**2+3d0)*2d0/(eta_p-eta_m)
          as(i) = 0.25d0*(-6d0*xi)*4d0/(eta_p-eta_m)**2
          !a(i)  = 0.0625d0*( 3d0*xi**5-10d0*xi**3+15d0*xi+8d0)
          !ap(i) = 0.0625d0*(15d0*xi**4-30d0*xi**2+15d0)/(eta_p-eta_m)
          !as(i) = 0.0625d0*(60d0*xi**3-60d0*xi)/(eta_p-eta_m)**2
        enddo

        ! find r and it's derivatives
        aa = Rp/Req
        do i=1,grd(id)%nr
          zz = zeta(npts_sum+i)
          s(id)%r_map(i,:)= aa*zz+ a(i)*(rs(:)-aa*eta_p)+(1d0-a(i))*(rs_prev(:)-aa*eta_m)
          s(id)%r_z(i,:)  = aa   +ap(i)*(rs(:)-aa*eta_p)-    ap(i) *(rs_prev(:)-aa*eta_m)
          s(id)%r_t(i,:)  =        a(i)*rsp(:)          +(1d0-a(i))*rsp_prev(:)
          s(id)%r_zz(i,:) =       as(i)*(rs(:)-aa*eta_p)-    as(i) *(rs_prev(:)-aa*eta_m)
          s(id)%r_zt(i,:) =       ap(i)*rsp(:)          -    ap(i) *rsp_prev(:)
          s(id)%r_tt(i,:) =        a(i)*rss(:)          +(1d0-a(i))*rss_prev(:)
        enddo
        
        ! find extra geometrical terms
        do i=1,grd(id)%nr
          s(id)%cost(i,:) = cth(:)
          s(id)%sint(i,:) = sth(:)
          s(id)%cott(i,:) = cth(:)/sth(:)
        enddo

        npts_sum = npts_sum + grd(id)%nr

      enddo
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! outer domain
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      do i=1,grd(ndom)%nr
        zz = 1d0+0.5d0*(1d0-dcos(dble(i-1)*pi/dble(grd(ndom)%nr-1)))
        grd(ndom)%r(i) = zz
        s(ndom)%zeta(i,1:lres) = zz
        a(i)  =  2d0*zz**3 -  9d0*zz**2 + 12d0*zz - 4d0
        ap(i) =  6d0*zz**2 - 18d0*zz    + 12d0
        as(i) = 12d0*zz    - 18d0
      enddo

      aa = Rp/Req
      bb = 2d0*(1d0-aa)
      do i=1,grd(ndom)%nr
        zz = s(ndom)%zeta(i,1)
        s(ndom)%r_map(i,:)= aa*zz+bb+ a(i)*( rs(:)+aa-2d0)
        s(ndom)%r_z(i,:)  = aa      +ap(i)*( rs(:)+aa-2d0)
        s(ndom)%r_t(i,:)  =           a(i)* rsp(:)
        s(ndom)%r_zz(i,:) =          as(i)*( rs(:)+aa-2d0)
        s(ndom)%r_zt(i,:) =          ap(i)* rsp(:)
        s(ndom)%r_tt(i,:) =           a(i)* rss(:)
      enddo

      do i=1,grd(ndom)%nr
        s(ndom)%cost(i,:) = cth(:)
        s(ndom)%sint(i,:) = sth(:)
        s(ndom)%cott(i,:) = cth(:)/sth(:)
      enddo

      deallocate(a,ap,as,r_spec,r_temp,rs,rsp,rss,rs_prev, &
                 rsp_prev,rss_prev,w,llr_spec)
      end subroutine make_mapping

!--------------------------------------------------------------
! This checks the mapping to make sure it agrees with the
! contents of the file with the model.
!
! IMPORTANT: this test is only valid when lres = 2*nthm
!--------------------------------------------------------------
      subroutine check_mapping()

      implicit none
      integer id, i, j
      double precision erreur

      erreur = 0d0
      do id=1,ndom-1
        do i=1,grd(id)%nr
          do j=1,nthm
            erreur = erreur + (s(id)%r_map(i,j+nthm)-s(id)%rr(i,j)/Req)**2
          enddo
        enddo
      enddo

      erreur = sqrt(erreur/(nthm*nrm))
      print*,"Erreur = ",erreur
      end subroutine check_mapping

!--------------------------------------------------------------
! This function evaluates a function decomposed over the
! spherical harmonic basis at a given z = cos(theta)
!--------------------------------------------------------------
      double precision function eval_ylm(f,nn,z)
        
      implicit none
      integer nn
      double precision f(1:nn), z
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
      end function eval_ylm

!--------------------------------------------------------------
! This subroutine prints different fields so as to help debug
! the program
!--------------------------------------------------------------
      subroutine write_fields()

      implicit none
      integer i, j, id
      character*(2), parameter :: dir="./"

      open(unit=999,file=dir//"delme_grid",status="unknown")
      write(999,*) nrm,lres
      write(999,'(a)') "# x"
      write(999,101) (((s(id)%r_map(i,j)*s(id)%sint(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# y"
      write(999,101) (((s(id)%r_map(i,j)*s(id)%cost(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,*) grd(ndom)%nr,lres
      write(999,'(a)') "# xe"
      write(999,101) ((s(ndom)%r_map(i,j)*s(ndom)%sint(i,j),i=1,grd(ndom)%nr),j=1,lres)
      write(999,'(a)') "# ye"
      write(999,101) ((s(ndom)%r_map(i,j)*s(ndom)%cost(i,j),i=1,grd(ndom)%nr),j=1,lres)
      close(999)

      open(unit=999,file=dir//"delme_fields",status="unknown")
      write(999,*) nrm,lres
      write(999,'(a)') "# rhom"
      write(999,101) (((s(id)%rhom(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# rhom_z"
      write(999,101) (((s(id)%rhom_z(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# rhom_t"
      write(999,101) (((s(id)%rhom_t(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# pm"
      write(999,101) (((s(id)%pm(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# pm_z"
      write(999,101) (((s(id)%pm_z(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# pm_t"
      write(999,101) (((s(id)%pm_t(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# pm_zz"
      write(999,101) (((s(id)%pm_zz(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# pm_tz"
      write(999,101) (((s(id)%pm_tz(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# pm_zt"
      write(999,101) (((s(id)%pm_zt(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# pm_tt"
      write(999,101) (((s(id)%pm_tt(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# pm_ez"
      write(999,101) (((s(id)%pm_ez(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# Gamma1"
      write(999,101) (((s(id)%Gamma1(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# c2"
      write(999,101) (((s(id)%c2(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# NNtoz"
      write(999,101) (((s(id)%NNtoz(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# Rota"
      write(999,101) (((s(id)%Rota(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# Rota_z"
      write(999,101) (((s(id)%Rota_z(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# Rota_t"
      write(999,101) (((s(id)%Rota_t(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# grd_pe_z"
      write(999,101) (((s(id)%grd_pe_z(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# grd_pe_t"
      write(999,101) (((s(id)%grd_pe_t(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# grd_pe_zz"
      write(999,101) (((s(id)%grd_pe_zz(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# grd_pe_zt"
      write(999,101) (((s(id)%grd_pe_zt(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# grd_pe_tz"
      write(999,101) (((s(id)%grd_pe_tz(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# grd_pe_tt"
      write(999,101) (((s(id)%grd_pe_tt(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# det_rhom_pm"
      write(999,101) (((s(id)%det_rhom_pm(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      close(999)

      open(unit=999,file=dir//"delme_geometry",status="unknown")
      write(999,*) nrm,lres
      write(999,'(a)') "# r_map"
      write(999,101) (((s(id)%r_map(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# r_z"
      write(999,101) (((s(id)%r_z(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# r_t"
      write(999,101) (((s(id)%r_t(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# r_zz"
      write(999,101) (((s(id)%r_zz(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# r_zt"
      write(999,101) (((s(id)%r_zt(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# r_tt"
      write(999,101) (((s(id)%r_tt(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# zeta"
      write(999,101) (((s(id)%zeta(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# cost"
      write(999,101) (((s(id)%cost(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# sint"
      write(999,101) (((s(id)%sint(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# cott"
      write(999,101) (((s(id)%cott(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# roz"
      write(999,101) (((s(id)%roz(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      write(999,'(a)') "# rrt"
      write(999,101) (((s(id)%rrt(i,j),i=1,grd(id)%nr),id=1,ndom-1),j=1,lres)
      close(999)

      open(unit=999,file=dir//"delme_constants",status="unknown")
      write(999,'(a)') "# Req"
      write(999,'(1pe22.15)') Req
      write(999,'(a)') "# Rp"
      write(999,'(1pe22.15)') Rp
      close(999)

      open(unit=999,file=dir//"plotme",status="unknown")
      write(999,'(a)') "# grd_pe_zz"
      j = lres/4
      do id=1,ndom-1
        do i=1,grd(id)%nr
          write(999,*)  s(id)%r_map(i,j),s(id)%grd_pe_zz(i,j)
        enddo
      enddo
      close(999)

101   format(3(1pe22.15,2X))
      stop
      end subroutine write_fields
!-------------------------------------------------------------------------
      end module
