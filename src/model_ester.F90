#include "config.h"

module model
    use abstract_model_mod
    use mod_grid, only: ndomains, grd
    use inputs, only: lres, mass, rota_avg
    use iso_c_binding
    use mod_legendre
    use derivative, only: init_derive_cheb

    implicit none

    type, extends(abstract_model) :: ester_model
    contains
        procedure :: init => init_ester_model
        procedure :: get_field => getField
    end type ester_model

    type STAR_DOMAIN

        ! outer surface of domain
        real(kind=c_double), pointer, dimension(:)   :: Rs

        ! raw variables (of dimension grd(id)%nr x nthm)
        real(kind=c_double), pointer, dimension(:, :) ::   &
            rr, p_raw, rho_raw, Gamma1_raw, Rota_raw,   &
            p_z_raw, rho_z_raw, Rota_z_raw

        ! final variables (of dimension grd(id)%nr x lres)
        real(kind=c_double), pointer, dimension(:, :) ::       &
            pm, rhom, Gamma1, pm_z, pm_t, rhom_z, rhom_t,   &
            c2, NNtoz, Rota, Rota_z, Rota_t, pm_zz,         &
            pm_zt, pm_tz, pm_tt, pm_ez, grd_pe_z,           &
            grd_pe_t, grd_pe_zz, grd_pe_zt, grd_pe_tz,      &
            grd_pe_tt, det_rhom_pm
        real(kind=c_double), pointer, dimension(:, :) :: r_t, r_z, r_map, &
            r_zz, r_zt, r_tt, zeta, cost, sint, cott, rrt, roz

    end type STAR_DOMAIN

    integer, save :: nrm, nthm, ndom, nconv, npts_max
    real(kind=c_double), save :: mass_input, radius, luminosity, rotation
    real(kind=c_double), save :: X, Z, Xc, Req, Rp
    real(kind=c_double), save :: rho_ref, p_ref, t_ref
    real(kind=c_double), allocatable, save :: theta(:), zeta(:), sth(:), cth(:)
    type(STAR_DOMAIN), allocatable, save :: s(:)
    logical, save :: first_model = .true.

    ! different physical constants
    real(kind=c_double), parameter :: pi = 3.141592653589793d0
    real(kind=c_double), parameter :: Lambda = 4d0*pi
    real(kind=c_double), parameter :: solar_mass   = 1.9891d33  !g
    real(kind=c_double), parameter :: solar_radius = 6.95508d10 !cm
    real(kind=c_double), parameter :: G = 6.672d-8 !cm^3.g^-1.s^-2

contains

    !------------------------------------------------------------------------
    subroutine init_ester_model(this, filename)

        class(ester_model), target :: this
        character(len=*), intent(in) :: filename

        call read_model(filename)
        call interpolate_model()
        call make_mapping()
        call find_NNtoz()
        call find_grd_pm()
        call find_grd_pe()

        model_ptr => this

    end subroutine init_ester_model


    subroutine read_model(filename)

        character(len=*), intent(in) :: filename
        integer(c_int), allocatable :: npts(:)

        integer :: i, j
        double precision, allocatable :: r_spec(:), r_temp(:)
        double precision :: aux

#ifdef USE_LIBESTER

        call cpp_read_ester_model(filename, nrm, nthm, ndom)

        ndom = ndom+1
        if (ndomains /= ndom) then
            print*, "Incorrect number of domains"
            return
        endif

        if (allocated(s)) call clear_model()

        allocate(s(ndom), theta(nthm), zeta(nrm))

        allocate(npts(ndom-1))

        call get_npts(npts)
        do i=1, ndom-1
            grd(i)%nr = npts(i)
        enddo
        call get_nex(grd(ndom)%nr)
        deallocate(npts)

        do i=1, ndom
            if (npts_max < grd(i)%nr) npts_max = grd(i)%nr
        enddo
        call allocate_model()

        call get_mass(mass_input)
        call get_radius(radius)
        call get_lum(luminosity)
        call get_omega(rotation)
        call get_X(X)
        call get_Z(Z)
        call get_Xc(Xc)

        call get_theta(theta)
        call get_zeta(zeta)

        do i=1, ndom-1
            call get_Rs(i,  s(i)%Rs)
            call get_rr(i,  s(i)%rr)
            call get_p(i,   s(i)%p_raw)
            call get_rho(i, s(i)%rho_raw)
            call get_G1(i,  s(i)%Gamma1_raw)
            call get_w(i,   s(i)%Rota_raw)
        enddo

        ! Find equatorial and polar radii
        allocate(r_temp(2*nthm), r_spec(2*nthm))
        do j=1, nthm
            r_temp(nthm+j)   = s(ndom-1)%Rs(j)
            r_temp(nthm+1-j) = s(ndom-1)%Rs(j)
        enddo

        r_spec = 0d0
        call legendre(r_temp, r_spec, 2*nthm, 0, 1)
        Rp  = eval_ylm(r_spec, 2*nthm, 1d0)
        Req = eval_ylm(r_spec, 2*nthm, 0d0)
        deallocate(r_temp, r_spec)

        ! sanity check on mapping
        do i=2, ndom-1
            aux = s(i)%Rs(nthm)-s(i-1)%Rs(nthm)
            do j=1, nthm-1
                if ((s(i)%Rs(j)-s(i-1)%Rs(j)).lt.aux) then
                    print*, "Problem with make_mapping."
                    stop
                endif
            enddo
        enddo

        rho_ref = mass_input/Req**3
        p_ref   = G*mass_input**2/Req**4
        t_ref   = sqrt(Req**3/(G*mass_input))
        mass    = mass_input/solar_mass
        rota_avg= rotation*t_ref

#else
        print*, "TOP was compiled without libester support..."
        print*, "Cannot read Ester models"
#endif

    end subroutine read_model
    !------------------------------------------------------------------------
    subroutine getField(this, fname, field)

        class(ester_model) :: this
        character(len=*), intent(in) :: fname
        real(kind=8), allocatable, intent(out) :: field(:, :)

        if (fname == 'zeta') then
            allocate(field(nrm, 1))
            field(:, 1) = zeta
        elseif (fname == 'theta') then
            allocate(field(1, nthm))
            field(1, :) = theta
        else
            print*, 'Unknown field:', fname
            allocate(field(1, 1))
            field = 0.0
        endif

    end subroutine getField
    !------------------------------------------------------------------------
    subroutine init_model(filename)

        implicit none
        character(len=*), intent(in) :: filename
        class(ester_model), allocatable :: model

        allocate(model)

        call model%init(filename)

    end subroutine init_model
    !------------------------------------------------------------------------
    subroutine allocate_model()

        implicit none
        integer id

        do id=1, ndom
            allocate(grd(id)%r(grd(id)%nr))
        enddo

        do id=1, ndom-1
            allocate(s(id)%p_raw(grd(id)%nr, nthm))
            allocate(s(id)%p_z_raw(grd(id)%nr, nthm))
            allocate(s(id)%rho_raw(grd(id)%nr, nthm))
            allocate(s(id)%rho_z_raw(grd(id)%nr, nthm))
            allocate(s(id)%Gamma1_raw(grd(id)%nr, nthm))
            allocate(s(id)%Rota_raw(grd(id)%nr, nthm))
            allocate(s(id)%Rota_z_raw(grd(id)%nr, nthm))
            allocate(s(id)%pm(grd(id)%nr, lres))
            allocate(s(id)%rhom(grd(id)%nr, lres))
            allocate(s(id)%Gamma1(grd(id)%nr, lres))
            allocate(s(id)%Rota(grd(id)%nr, lres))
            allocate(s(id)%c2(grd(id)%nr, lres))
            allocate(s(id)%pm_z(grd(id)%nr, lres))
            allocate(s(id)%pm_t(grd(id)%nr, lres))
            allocate(s(id)%rhom_z(grd(id)%nr, lres))
            allocate(s(id)%rhom_t(grd(id)%nr, lres))
            allocate(s(id)%Rota_z(grd(id)%nr, lres))
            allocate(s(id)%Rota_t(grd(id)%nr, lres))
            allocate(s(id)%NNtoz(grd(id)%nr, lres))
            allocate(s(id)%pm_zz(grd(id)%nr, lres))
            allocate(s(id)%pm_tz(grd(id)%nr, lres))
            allocate(s(id)%pm_zt(grd(id)%nr, lres))
            allocate(s(id)%pm_tt(grd(id)%nr, lres))
            allocate(s(id)%pm_ez(grd(id)%nr, lres))
            allocate(s(id)%grd_pe_z(grd(id)%nr, lres))
            allocate(s(id)%grd_pe_t(grd(id)%nr, lres))
            allocate(s(id)%grd_pe_zz(grd(id)%nr, lres))
            allocate(s(id)%grd_pe_zt(grd(id)%nr, lres))
            allocate(s(id)%grd_pe_tz(grd(id)%nr, lres))
            allocate(s(id)%grd_pe_tt(grd(id)%nr, lres))
            allocate(s(id)%det_rhom_pm(grd(id)%nr, lres))
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

        do id=1, ndom
            allocate(s(id)%Rs(nthm))
            allocate(s(id)%rr(grd(id)%nr, nthm))
            allocate(s(id)%r_map(grd(id)%nr, lres))
            allocate(s(id)%r_z(grd(id)%nr, lres))
            allocate(s(id)%r_t(grd(id)%nr, lres))
            allocate(s(id)%r_zz(grd(id)%nr, lres))
            allocate(s(id)%r_zt(grd(id)%nr, lres))
            allocate(s(id)%r_tt(grd(id)%nr, lres))
            allocate(s(id)%zeta(grd(id)%nr, lres))
            allocate(s(id)%cost(grd(id)%nr, lres))
            allocate(s(id)%sint(grd(id)%nr, lres))
            allocate(s(id)%cott(grd(id)%nr, lres))
            allocate(s(id)%roz(grd(id)%nr, lres))
            allocate(s(id)%rrt(grd(id)%nr, lres))
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

        do id=1, ndom
            deallocate(grd(id)%r)
        enddo

        do id=1, ndom-1
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

        do id=1, ndom
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

    double precision function eval_ylm(f, nn, z)

        implicit none
        integer nn
        double precision f(1:nn), z
        double precision yl1, yl2, yl3
        integer l

        yl1 = sqrt(0.25d0/pi)
        yl2 = sqrt(0.75d0/pi)*z
        eval_ylm = yl1*f(1)+yl2*f(2)
        ! Beware: f(l+1) corresponds to f at l
        do l = 2, nn-1
            yl3 = (sqrt(dble(4*l*l-1))/dble(l))*z*yl2  &
                + (1d0/dble(l)-1d0)*sqrt(dble(2*l+1)/(2*l-3))*yl1
            eval_ylm = eval_ylm + yl3*f(l+1)
            yl1 = yl2
            yl2 = yl3
        enddo
        return
    end function eval_ylm
!-------------------------------------------------------------------------
!  This calculates the array NNtoz: it needs to be called after the
!  mapping is created.
!-------------------------------------------------------------------------
      subroutine find_NNtoz()

          implicit none
          integer id, i, j

          id = 1
          do j=1, lres
              do i=2, grd(id)%nr
                  s(id)%NNtoz(i, j) = (-s(id)%c2(i, j)*s(id)%rhom_t(i, j)  &
                      +s(id)%pm_t(i, j))/s(id)%zeta(i, j)
              enddo
              s(id)%NNtoz(1, j) = 0d0
          enddo
      end subroutine find_NNtoz

!-------------------------------------------------------------------------
!  This interpolates the model to a higher angular resolution
!-------------------------------------------------------------------------
      subroutine interpolate_model()

          use mod_legendre
          use derivative

          implicit none
          integer i, j, id, npid, nrs, nrf
          double precision, allocatable :: mat_raw(:, :), mat(:, :)
          type(DERMAT) :: dmat

          ! Don't forget to initialise the radial derivative of p, rho and
          ! Rota.
          nrs = 1
          nrf = grd(1)%nr
          do id=1, ndomains-1
          call init_derive_cheb(dmat, zeta(nrs), grd(id)%nr, 1, 0)
          call dgemm("N", "N", grd(id)%nr, nthm, grd(id)%nr, 1d0, &
              dmat%derive(:, :, 1), grd(id)%nr,          &
              s(id)%p_raw, grd(id)%nr, 0d0,             &
              s(id)%p_z_raw, grd(id)%nr)
          call dgemm("N", "N", grd(id)%nr, nthm, grd(id)%nr, 1d0, &
              dmat%derive(:, :, 1), grd(id)%nr,          &
              s(id)%rho_raw, grd(id)%nr, 0d0,           &
              s(id)%rho_z_raw, grd(id)%nr)
          call dgemm("N", "N", grd(id)%nr, nthm, grd(id)%nr, 1d0, &
              dmat%derive(:, :, 1), grd(id)%nr,          &
              s(id)%Rota_raw, grd(id)%nr, 0d0,           &
              s(id)%Rota_z_raw, grd(id)%nr)
          nrs = nrs + grd(id)%nr
          nrf = nrs + grd(id+1)%nr - 1
          call clear_derive(dmat)
          enddo

          allocate(mat_raw(npts_max, nthm*2), mat(npts_max, lres))

          do id=1, ndom-1
          npid = grd(id)%nr
          do j=1, nthm
          mat_raw(1:npid, nthm+j)   = s(id)%p_raw(:, j)/p_ref
          mat_raw(1:npid, nthm+1-j) = s(id)%p_raw(:, j)/p_ref
          enddo
          mat = 0d0
          call legendre(mat_raw(1:npid, :), mat(1:grd(id)%nr, 1:(2*nthm)), &
              (2*nthm), grd(id)%nr, 0, 1)
          call legendre(mat(1:grd(id)%nr, 1:lres), s(id)%pm, lres,         &
              grd(id)%nr, 0, -1)
          call legendrep(mat(1:grd(id)%nr, 1:lres), s(id)%pm_t, lres,      &
              grd(id)%nr, 0)

          do j=1, nthm
          mat_raw(1:npid, nthm+j)   = s(id)%rho_raw(:, j)/rho_ref
          mat_raw(1:npid, nthm+1-j) = s(id)%rho_raw(:, j)/rho_ref
          enddo
          mat = 0d0
          call legendre(mat_raw(1:npid, :), mat(1:grd(id)%nr, 1:(2*nthm)), &
              (2*nthm), grd(id)%nr, 0, 1)
          call legendre(mat(1:grd(id)%nr, 1:lres), s(id)%rhom, lres,       &
              grd(id)%nr, 0, -1)
          call legendrep(mat(1:grd(id)%nr, 1:lres), s(id)%rhom_t, lres,    &
              grd(id)%nr, 0)

          do j=1, nthm
          mat_raw(1:npid, nthm+j)   = s(id)%Rota_raw(:, j)*t_ref
          mat_raw(1:npid, nthm+1-j) = s(id)%Rota_raw(:, j)*t_ref
          enddo
          mat = 0d0
          call legendre(mat_raw(1:npid, :), mat(1:grd(id)%nr, 1:(2*nthm)), &
              (2*nthm), grd(id)%nr, 0, 1)
          call legendre(mat(1:grd(id)%nr, 1:lres), s(id)%Rota, lres,       &
              grd(id)%nr, 0, -1)
          call legendrep(mat(1:grd(id)%nr, 1:lres), s(id)%Rota_t, lres,    &
              grd(id)%nr, 0)

          do j=1, nthm
          mat_raw(1:npid, nthm+j)   = s(id)%Gamma1_raw(:, j)
          mat_raw(1:npid, nthm+1-j) = s(id)%Gamma1_raw(:, j)
          enddo
          mat = 0d0
          call legendre(mat_raw(1:npid, :), mat(1:grd(id)%nr, 1:(2*nthm)), &
              (2*nthm), grd(id)%nr, 0, 1)
          call legendre(mat(1:grd(id)%nr, 1:lres), s(id)%Gamma1, lres, &
              grd(id)%nr, 0, -1)

          s(id)%c2 = s(id)%Gamma1*s(id)%pm/s(id)%rhom

          do j=1, nthm
          mat_raw(1:npid, nthm+j)   = s(id)%p_z_raw(:, j)/p_ref
          mat_raw(1:npid, nthm+1-j) = s(id)%p_z_raw(:, j)/p_ref
          enddo
          mat = 0d0
          call legendre(mat_raw(1:npid, :), mat(1:grd(id)%nr, 1:(2*nthm)), &
              (2*nthm), grd(id)%nr, 0, 1)
          call legendre(mat(1:grd(id)%nr, 1:lres), s(id)%pm_z, lres, &
              grd(id)%nr, 0, -1)

          do j=1, nthm
          mat_raw(1:npid, nthm+j)   = s(id)%rho_z_raw(:, j)/rho_ref
          mat_raw(1:npid, nthm+1-j) = s(id)%rho_z_raw(:, j)/rho_ref
          enddo
          mat = 0d0
          call legendre(mat_raw(1:npid, :), mat(1:grd(id)%nr, 1:(2*nthm)), &
              (2*nthm), grd(id)%nr, 0, 1)
          call legendre(mat(1:grd(id)%nr, 1:lres), s(id)%rhom_z, lres, &
              grd(id)%nr, 0, -1)

          do j=1, nthm
          mat_raw(1:npid, nthm+j)   = s(id)%Rota_z_raw(:, j)*t_ref
          mat_raw(1:npid, nthm+1-j) = s(id)%Rota_z_raw(:, j)*t_ref
          enddo
          mat = 0d0
          call legendre(mat_raw(1:npid, :), mat(1:grd(id)%nr, 1:(2*nthm)), &
              (2*nthm), grd(id)%nr, 0, 1)
          call legendre(mat(1:grd(id)%nr, 1:lres), s(id)%Rota_z, lres, &
              grd(id)%nr, 0, -1)

          enddo
          deallocate(mat, mat_raw)
      end subroutine interpolate_model

!-------------------------------------------------------------------------
!  This sets up the mapping for the model
!-------------------------------------------------------------------------
      subroutine make_mapping()

          use mod_legendre

          implicit none
          integer i, j, l, id, npts_sum
          double precision :: xi, zz, eta_p, eta_m, aa, bb
          double precision, allocatable, dimension(:) :: a, ap, as, r_spec,  &
              r_temp, rs, rsp, rss, rs_prev, rsp_prev, rss_prev, w, llr_spec

          if (allocated(cth)) deallocate(cth)
          if (allocated(sth)) deallocate(sth)

          allocate(a(npts_max), ap(npts_max), as(npts_max), r_spec(lres), &
              r_temp(2*nthm), rs(lres), rsp(lres), rss(lres),       &
              rs_prev(lres), rsp_prev(lres), rss_prev(lres), w(lres), &
              llr_spec(lres), cth(lres), sth(lres))

          ! calculation of sin(theta), cos(theta)
          call gauleg(-1d0, 1d0, cth, w, lres)
          sth = sqrt(1-cth**2)

          npts_sum = 0

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! inner domain
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          rs_prev = 0d0
          rsp_prev= 0d0
          rss_prev= 0d0
          do j=1, nthm
          r_temp(nthm+j)   = s(1)%Rs(j)/Req
          r_temp(nthm+1-j) = s(1)%Rs(j)/Req
          enddo
          r_spec = 0d0
          call legendre(r_temp, r_spec(1:(2*nthm)), 2*nthm, 0, 1)
          call legendre(r_spec, rs, lres, 0, -1)
          call legendrep(r_spec, rsp, lres, 0)
          do l=1, lres ! Beware: r_spec(l) is r_spec of (ell-1)
          llr_spec(l)=-dble(l*(l-1))*r_spec(l)
          enddo
          call legendre(llr_spec(1:lres), rss(1:lres), lres, 0, -1)
          rss = rss - cth*rsp/sth

          ! find radial variables
          eta_m = zeta(1)
          eta_p = zeta(grd(1)%nr)
          do i=1, grd(1)%nr
          zz = zeta(i)
          grd(1)%r(i) = zz
          s(1)%zeta(i, 1:lres) = zz
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
          do i=1, grd(1)%nr
          zz = zeta(i)
          s(1)%r_map(i, :)= aa*zz + a(i)*( rs(:)-aa*eta_p)
          s(1)%r_z(i, :)  = aa    +ap(i)*( rs(:)-aa*eta_p)
          s(1)%r_t(i, :)  =         a(i)* rsp(:)
          s(1)%r_zz(i, :) =        as(i)*( rs(:)-aa*eta_p)
          s(1)%r_zt(i, :) =        ap(i)* rsp(:)
          s(1)%r_tt(i, :) =         a(i)* rss(:)
          enddo

          ! find extra geometrical terms
          do i=1, grd(1)%nr
          s(1)%cost(i, :) = cth(:)
          s(1)%sint(i, :) = sth(:)
          s(1)%cott(i, :) = cth(:)/sth(:)
          enddo

          ! define roz and rrt:
          do i=2, grd(1)%nr
          s(1)%roz(i, :) = s(1)%r_map(i, :)/s(1)%zeta(i, :)
          s(1)%rrt(i, :) = s(1)%r_t(i, :)/s(1)%r_map(i, :)
          enddo
          s(1)%roz(1, :) = s(1)%r_z(1, :)
          s(1)%rrt(1, :) = s(1)%r_zt(1, :)/s(1)%r_z(1, :)

          npts_sum = grd(1)%nr

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! intermediate domains
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          do id=2, ndom-1

          ! inner bound = previous outer bound
          rs_prev  = rs
          rsp_prev = rsp
          rss_prev = rss

          ! interpolate outer bound
          do j=1, nthm
          r_temp(nthm+j)   = s(id)%Rs(j)/Req
          r_temp(nthm+1-j) = s(id)%Rs(j)/Req
          enddo
          r_spec = 0d0
          call legendre(r_temp, r_spec(1:(2*nthm)), 2*nthm, 0, 1)
          call legendre(r_spec, rs, lres, 0, -1)
          call legendrep(r_spec, rsp, lres, 0)
          do l=1, lres ! Beware: r_spec(l) is r_spec of (ell-1)
          llr_spec(l)=-dble(l*(l-1))*r_spec(l)
          enddo
          call legendre(llr_spec(1:lres), rss(1:lres), lres, 0, -1)
          rss = rss - cth*rsp/sth ! Please don't forget this line!

          ! find radial variables
          eta_m = zeta(npts_sum+1)
          eta_p = zeta(npts_sum+grd(id)%nr)
          do i=1, grd(id)%nr
          zz = zeta(npts_sum+i)
          grd(id)%r(i) = zz
          s(id)%zeta(i, 1:lres) = zz
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
          do i=1, grd(id)%nr
          zz = zeta(npts_sum+i)
          s(id)%r_map(i, :)= aa*zz+ a(i)*(rs(:)-aa*eta_p)+(1d0-a(i))*(rs_prev(:)-aa*eta_m)
          s(id)%r_z(i, :)  = aa   +ap(i)*(rs(:)-aa*eta_p)-    ap(i) *(rs_prev(:)-aa*eta_m)
          s(id)%r_t(i, :)  =        a(i)*rsp(:)          +(1d0-a(i))*rsp_prev(:)
          s(id)%r_zz(i, :) =       as(i)*(rs(:)-aa*eta_p)-    as(i) *(rs_prev(:)-aa*eta_m)
          s(id)%r_zt(i, :) =       ap(i)*rsp(:)          -    ap(i) *rsp_prev(:)
          s(id)%r_tt(i, :) =        a(i)*rss(:)          +(1d0-a(i))*rss_prev(:)
          enddo

          ! find extra geometrical terms
          do i=1, grd(id)%nr
          s(id)%cost(i, :) = cth(:)
          s(id)%sint(i, :) = sth(:)
          s(id)%cott(i, :) = cth(:)/sth(:)
          enddo

          npts_sum = npts_sum + grd(id)%nr

          enddo

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! outer domain
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          do i=1, grd(ndom)%nr
          zz = 1d0+0.5d0*(1d0-dcos(dble(i-1)*pi/dble(grd(ndom)%nr-1)))
          grd(ndom)%r(i) = zz
          s(ndom)%zeta(i, 1:lres) = zz
          a(i)  =  2d0*zz**3 -  9d0*zz**2 + 12d0*zz - 4d0
          ap(i) =  6d0*zz**2 - 18d0*zz    + 12d0
          as(i) = 12d0*zz    - 18d0
          enddo

          aa = Rp/Req
          bb = 2d0*(1d0-aa)
          do i=1, grd(ndom)%nr
          zz = s(ndom)%zeta(i, 1)
          s(ndom)%r_map(i, :)= aa*zz+bb+ a(i)*( rs(:)+aa-2d0)
          s(ndom)%r_z(i, :)  = aa      +ap(i)*( rs(:)+aa-2d0)
          s(ndom)%r_t(i, :)  =           a(i)* rsp(:)
          s(ndom)%r_zz(i, :) =          as(i)*( rs(:)+aa-2d0)
          s(ndom)%r_zt(i, :) =          ap(i)* rsp(:)
          s(ndom)%r_tt(i, :) =           a(i)* rss(:)
          enddo

          do i=1, grd(ndom)%nr
          s(ndom)%cost(i, :) = cth(:)
          s(ndom)%sint(i, :) = sth(:)
          s(ndom)%cott(i, :) = cth(:)/sth(:)
          enddo

          deallocate(a, ap, as, r_spec, r_temp, rs, rsp, rss, rs_prev, &
              rsp_prev, rss_prev, w, llr_spec)
      end subroutine make_mapping

!-------------------------------------------------------------------------
! This finds extra gradients of the pressure, used for an alternate
! boundary condition.
!-------------------------------------------------------------------------
      subroutine find_grd_pm()

          use mod_legendre
          use derivative

          implicit none
          integer i, j, id, npid, nrs, nrf
          double precision, allocatable :: mat(:, :), smat(:, :)
          type(DERMAT) :: dmat

          ! find various indices
          nrs = 1
          do id=1, ndomains-2
          nrs = nrs + grd(id)%nr
          enddo
          id = ndomains - 1
          nrf = nrs + grd(id)%nr - 1

          ! initialise mat and dmat
          call init_derive_cheb(dmat, zeta(nrs), grd(id)%nr, 1, 0)
          allocate(mat(grd(id)%nr, lres))
          allocate(smat(grd(id)%nr, lres))

          ! find pm_zz
          mat = s(id)%zeta**2*s(id)%pm_z/(s(id)%r_map**2*s(id)%r_z)
          call dgemm("N", "N", grd(id)%nr, lres, grd(id)%nr, 1d0, &
              dmat%derive(:, :, 1), grd(id)%nr,          &
              mat, grd(id)%nr, 0d0,             &
              s(id)%pm_zz, grd(id)%nr)

          ! find pm_tz
          call legendre(mat, smat, lres, grd(id)%nr, 0, 1)
          call legendrep(smat, s(id)%pm_tz, lres, grd(id)%nr, 0)

          ! find pm_zt
          mat = s(id)%zeta*s(id)%pm_t/(s(id)%r_map**2*s(id)%r_z)
          call dgemm("N", "N", grd(id)%nr, lres, grd(id)%nr, 1d0, &
              dmat%derive(:, :, 1), grd(id)%nr,          &
              mat, grd(id)%nr, 0d0,             &
              s(id)%pm_zt, grd(id)%nr)

          ! find pm_tt: (quick and ugly solution for d(sin t P(cos t))/dt)
          mat = mat/s(id)%sint
          s(id)%pm_tt = s(id)%cost*mat
          call legendre(mat, smat, lres, grd(id)%nr, 0, 1)
          call legendrep(smat, mat, lres, grd(id)%nr, 0)
          s(id)%pm_tt = s(id)%pm_tt + (s(id)%cost**2-1d0)*mat

          ! find pm_ez (just a convenient quantity):
          s(id)%pm_ez = (s(id)%r_map**2+s(id)%r_t**2)*s(id)%pm_z  &
              / (s(id)%r_map**2*s(id)%r_z**2)             &
              -  s(id)%r_t*s(id)%pm_t                     &
              / (s(id)%r_map**2*s(id)%r_z)

          deallocate(mat, smat)
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
          double precision, allocatable :: aux(:, :)
          integer id, j, i, ii


          do id=1, ndomains-1
          call init_derive_cheb(dm, grd(id)%r, grd(id)%nr, 1, 0)

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
              s(id)%det_rhom_pm(1, :) = 0d0
              s(id)%grd_pe_z(1, :) =-s(id)%pm_z(1, :)/(s(id)%roz(1, :)**2 &
                  * s(id)%r_z(1, :)*s(id)%rhom(1, :))
              s(id)%grd_pe_t(1, :) = 0d0
          endif

          ! Numerically calculate zeta derivatives: this avoids various
          ! difficulties with singularities in the center of the star:
          call dgemm("N", "N", grd(id)%nr, lres, grd(id)%nr, 1d0, &
              dm%derive(:, :, 1), grd(id)%nr,             &
              s(id)%grd_pe_z, grd(id)%nr, 0d0,           &
              s(id)%grd_pe_zz, grd(id)%nr)
          call dgemm("N", "N", grd(id)%nr, lres, grd(id)%nr, 1d0, &
              dm%derive(:, :, 1), grd(id)%nr,             &
              s(id)%grd_pe_t, grd(id)%nr, 0d0,           &
              s(id)%grd_pe_zt, grd(id)%nr)

          ! Numerically calculate theta derivatives
          allocate(aux(grd(id)%nr, lres))
          call legendre(s(id)%grd_pe_z, aux, lres, grd(id)%nr, 0, 1)
          call legendrep(aux, s(id)%grd_pe_tz, lres, grd(id)%nr, 0)

          ! This is not the most elegant solution but the legendre grid
          ! does not go to theta = 0, so this should work (at reasonable
          ! resolutions)
          call legendre(s(id)%grd_pe_t/s(id)%sint, aux, lres, grd(id)%nr, 0, 1)
          call legendrep(aux, s(id)%grd_pe_tt, lres, grd(id)%nr, 0)
          s(id)%grd_pe_tt = s(id)%sint*s(id)%grd_pe_tt+s(id)%cott*s(id)%grd_pe_t

          call clear_derive(dm)
          deallocate(aux)

          enddo

      end subroutine find_grd_pe

end module model
