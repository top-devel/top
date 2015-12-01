#include "config.h"

module model
    use abstract_model_mod
    use mod_grid, only: ndomains, grd
    use inputs, only: lres
    use iso_c_binding

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
    real(kind=c_double), pointer, save :: theta(:), zeta(:), sth(:), cth(:)
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
        type(c_ptr) :: ester_ptr
        integer(c_int), allocatable :: npts(:)

        integer :: i

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
        call get_nex(grd(ndom)%nr)
        do i=1, ndom-1
            grd(i)%nr = npts(i)
        enddo
        deallocate(npts)

        do i=1, ndom
            if (npts_max < grd(i)%nr) npts_max = grd(i)%nr
        enddo
        call allocate_model()

        call get_theta(theta)
        call get_zeta(zeta)

        do i=1, ndom-1
            call get_rho(i, s(i)%rho_raw)
        enddo

        model_ptr => this
#else
        print*, "TOP was compiled without libester support..."
        print*, "Cannot read Ester models"
#endif

    end subroutine init_ester_model
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

end module model
