module model
    use mod_grid
    use derivative
    use mod_legendre
    use inputs, only: lres, rota, mass, dertype, orderFD, pert_model, grid_type
    use abstract_model_mod, only: model_ptr, abstract_model

    implicit none

    ! raw variables:
    double precision, allocatable, save :: var(:,:),glob(:)
    integer,save :: iconst,ivar,ivers

    ! some 1D fields:
    integer, save :: nrmod
    double precision, allocatable, save :: Gamma1_1D(:),uhx(:),&
        rho1D(:),p1D(:),&
        r_model(:),pe1D(:), &
        mass1D(:), NN1D(:)

    ! stellar parameters:
    double precision, save :: radius ! in cm
    double precision, save :: p_ref          ! in g.cm^-1.s^-2
    double precision, save :: rho_ref        ! in g.cm^-3
    double precision, save :: phi_ref        ! in cm.s^-2

    ! workspace:
    double precision, allocatable, save :: ws1(:), ws2(:)
    double precision :: rota_pert

    ! variables needed for pulsation calculations
#ifdef USE_1D
    double precision, allocatable, dimension(:), save :: &
        rhom, rhom_z, pm, pm_z, c2, g_m, dg_m, Gamma1, NN, gm_r, &
        rF, diffrho_r2
#else
    double precision, allocatable, dimension(:,:), save :: rhom, NNt,&
        rhom_z, rhom_t, pm, pm_z, pm_t, Gamma1, c2,pe,pe_z, &
        pe_t, grd_pe_z, grd_pe_t, grd_pe_zz, grd_pe_zt,NNr, &
        grd_pe_tz, grd_pe_tt, NN
#endif
    ! /!\ NN is never perturbated
    double precision, save :: Lambda
#ifdef USE_1D
    double precision, save :: CC, rhom_zz
#endif

    ! geometric terms
#ifdef USE_1D
    double precision, allocatable, dimension(:), save :: r_map
#else
    double precision, allocatable, dimension(:,:), save :: r_t,r_z,r_map,&
        re_t,re_z,re_map,r_zz,r_zt,r_tt,re_zz,re_zt,&
        re_tt,zeta,cost,sint,cott,roz,rrt,r_aux
    double precision, allocatable, dimension(:), save :: cth, sth
#endif
    double precision, save :: K

    ! different physical constants
    double precision, parameter :: solar_mass = 1.98919d33 ! g
    double precision, parameter :: solar_radius = 6.959933134279d10 ! cm
    double precision, parameter :: G = 6.67259d-8          ! cm^3.g^-1.s^-2
    double precision, parameter :: pi = 3.141592653589793d0

    type, extends(abstract_model) :: cesam_model
    contains
        procedure :: init => init_cesam_model
        procedure :: get_field => cesam_get_field
        procedure :: get_grid => cesam_get_grid
        procedure :: get_grid_size => cesam_get_grid_size
    end type cesam_model

contains

    subroutine init_cesam_model(this, filename)

        class(cesam_model), target :: this
        character(len=*), intent(in) :: filename

#ifdef USE_1D
        call read_model(filename)
        call init_radial_grid()
        call init_fields()
#else
        call read_model(filename)
        ! call init_radial_grid_file()
        call init_radial_grid_g_modes()
        if (pert_model == 1) then
            rota_pert = rota
            print*, "Pert model: rota=", rota_pert
            call find_pert()
        else
            rota_pert = 0.d0
            print*, "No model pert"
            call no_pert() ! spherical case: but watch out because Rota is != 0
        endif
        call make_mapping()
        !call find_pe1D() ! integrate Poisson's equation
        call find_pe1D_alt() ! integrate gravity
        call init_fields()
        call init_grd_pe()
        ! call write_radial_grid() ! this stops the program
        ! call write_amdl()        ! this stops the program
        ! call write_fields()      ! this stops the program
#endif

        model_ptr => this

    end subroutine init_cesam_model
  
    subroutine cesam_get_field(this, fname, field)

        class(cesam_model) :: this
        character(len=*), intent(in) :: fname
        real(kind=8), allocatable, intent(out) :: field(:, :)

#ifdef USE_1D
        if (fname == 'rhom') then
            allocate(field(nr, 1))
            field(:, 1) = rhom
        elseif (fname == 'rhom_z') then
            allocate(field(nr, 1))
            field(:, 1) = rhom_z
        ! elseif (fname == 'rhom_t') then
        !     allocate(field(nr, 1))
        !     field(:, 1) = rhom_t
        elseif (fname == 'pm') then
            allocate(field(nr, 1))
            field(:, 1) = pm
        elseif (fname == 'pm_z') then
            allocate(field(nr, 1))
            field(:, 1) = pm_z
        ! elseif (fname == 'pm_t') then
        !     allocate(field(nr, 1))
        !     field(:, 1) = pm_t
        elseif (fname == 'Gamma1') then
            allocate(field(nr, 1))
            field(:, 1) = Gamma1
        elseif (fname == 'NN') then
            allocate(field(nr, 1))
            field(:, 1) = NN
        ! elseif (fname == 'NNr') then
        !     allocate(field(nr))
        !     field = NNr
        ! elseif (fname == 'NNt') then
        !     allocate(field(nr))
        !     field = NNt
        ! elseif (fname == 'pe') then
        !     allocate(field(nr))
        !     field = pe
        ! elseif (fname == 'pe_z') then
        !     allocate(field(nr))
        !     field = pe_z
        ! elseif (fname == 'pe_t') then
        !     allocate(field(nr))
        !     field = pe_t
        ! elseif (fname == 'grd_pe_z') then
        !     allocate(field(nr))
        !     field = grd_pe_z
        ! elseif (fname == 'grd_pe_t') then
        !     allocate(field(nr))
        !     field = grd_pe_t
        ! elseif (fname == 'grd_pe_zz') then
        !     allocate(field(nr))
        !     field = grd_pe_zz
        ! elseif (fname == 'grd_pe_zt') then
        !     allocate(field(nr))
        !     field = grd_pe_zt
        ! elseif (fname == 'grd_pe_tz') then
        !     allocate(field(nr))
        !     field = grd_pe_tz
        ! elseif (fname == 'grd_pe_tt') then
        !     allocate(field(nr))
        !     field = grd_pe_tt
        else
            allocate(field(1, 1))
            field = 0.d0
        endif
#else
        if (fname == 'rhom') then
            allocate(field(nr, lres))
            field = rhom
        elseif (fname == 'rhom_z') then
            allocate(field(nr, lres))
            field = rhom_z
        elseif (fname == 'rhom_t') then
            allocate(field(nr, lres))
            field = rhom_t
        elseif (fname == 'pm') then
            allocate(field(nr, lres))
            field = pm
        elseif (fname == 'pm_z') then
            allocate(field(nr, lres))
            field = pm_z
        elseif (fname == 'pm_t') then
            allocate(field(nr, lres))
            field = pm_t
        elseif (fname == 'Gamma1') then
            allocate(field(nr, lres))
            field = Gamma1
        elseif (fname == 'NN') then
            allocate(field(nr, lres))
            field(:, :) = NN
        elseif (fname == 'NNr') then
            allocate(field(nr, lres))
            field = NNr
        elseif (fname == 'NNt') then
            allocate(field(nr, lres))
            field = NNt
        elseif (fname == 'pe') then
            allocate(field(nr, lres))
            field = pe
        elseif (fname == 'pe_z') then
            allocate(field(nr, lres))
            field = pe_z
        elseif (fname == 'pe_t') then
            allocate(field(nr, lres))
            field = pe_t
        elseif (fname == 'grd_pe_z') then
            allocate(field(nr, lres))
            field = grd_pe_z
        elseif (fname == 'grd_pe_t') then
            allocate(field(nr, lres))
            field = grd_pe_t
        elseif (fname == 'grd_pe_zz') then
            allocate(field(nr, lres))
            field = grd_pe_zz
        elseif (fname == 'grd_pe_zt') then
            allocate(field(nr, lres))
            field = grd_pe_zt
        elseif (fname == 'grd_pe_tz') then
            allocate(field(nr, lres))
            field = grd_pe_tz
        elseif (fname == 'grd_pe_tt') then
            allocate(field(nr, lres))
            field = grd_pe_tt
        else
            allocate(field(1, 1))
            field = 0.d0
        endif
#endif

    end subroutine cesam_get_field

    subroutine cesam_get_grid(this, r, th, nr, nt)

        class(cesam_model) :: this
        real(kind=8), intent(out) :: r(nr, nt), th(nt)
        integer, intent(in) :: nr, nt

#ifdef USE_1D
        r(:, 1) = r_map
        th = 0
#else
        r = r_map
        th = acos(cth)
#endif

    end subroutine cesam_get_grid

    subroutine cesam_get_grid_size(this, n_r, n_t)

        class(cesam_model) :: this
        integer, intent(out) :: n_r, n_t
#ifdef USE_1D
        n_r = nr
        n_t = 1
#else
        n_r = grd(1)%nr
        n_t = lres
#endif

    end subroutine cesam_get_grid_size

    !--------------------------------------------------------------------------
    ! This subroutine calls all of the necessary subroutines to initialise
    ! the model
    !--------------------------------------------------------------------------
    subroutine init_model(modelfile)

        character(len=*), intent(in) :: modelfile
        class(cesam_model), allocatable :: model

        allocate(model)

        call model%init(modelfile)

    end subroutine
    !--------------------------------------------------------------------------
    ! This subroutine reads the model
    !--------------------------------------------------------------------------
    subroutine read_model(modelfile)

        integer i, j, nbelem
        double precision, allocatable :: aux(:)
        character(len=*), intent(in) :: modelfile

        if (allocated(var))  deallocate(var)
        if (allocated(glob)) deallocate(glob)

        ! Lecture du fichier CESAM:
        open(unit=37, form='formatted', status='old', file=trim(modelfile))

        ! Lecture des lignes sans interets 
        read(37, *)
        read(37, *)
        read(37, *)
        read(37, *)
        read(37, *)

140    format(7i10)
143    format(1p5d19.12)

        read(37,140) nrmod,iconst,ivar,nbelem
        nr = nrmod
        allocate(glob(iconst),var(ivar+nbelem,nrmod))

        read(37,143) (glob(i),i=1,iconst)
        ! Beware: CESAM models go in reverse order
        do j=1,nrmod
            read(37,143) (var(i,nrmod-j+1),i=1, ivar+nbelem)
        enddo    

        close(37)


        ! setting some global variables:
!        mass    = glob(1)*dexp(var(2, nrmod)) ! in g
!        radius  = var(1, nrmod) ! in cm
         mass    = glob(1) ! in g
         radius  = glob(2) ! in cm
        ! age     = glob(13) ! in some unknown unit
        p_ref   = G * mass**2 / radius**4
        rho_ref = mass / radius**3
        phi_ref = G * mass / radius
        Lambda = 4d0*3.141592653589793d0
        write(*,*) nrmod,iconst,ivar,nbelem
        write(*,*) 'mass:', mass
        write(*,*) 'radius:',radius
        ! write(*,*) 'age:',age
        write(*,*) 'glob:',glob


        ! convert mass into solar units:
        mass = mass/solar_mass

        ! create radial grid
        if (allocated(r_model)) deallocate(r_model)
        allocate(r_model(nrmod))
        do j=1, nrmod
            r_model(j) = var(1, j)/radius
        enddo
        ! calculate some 1D fields:
        if (allocated(Gamma1_1D)) deallocate(Gamma1_1D)
        if (allocated(rho1D))     deallocate(rho1D)
        if (allocated(p1D))       deallocate(p1D)
        if (allocated(mass1D))    deallocate(mass1D)
        if (allocated(NN1D))      deallocate(NN1D)

        allocate(Gamma1_1D(nrmod), rho1D(nrmod), p1D(nrmod),mass1D(nrmod),NN1D(nrmod))

        ! open(111,file='mon_model',form='formatted',status="unknown")
        do i=1,nrmod
            rho1D(i)     = var(5, i)/rho_ref
            p1D(i)       = var(4, i)/p_ref
            Gamma1_1D(i) = var(10, i)
            !mass1D(i)    = dexp(var(2, i)-var(2,nrmod))
            mass1D(i)    = dexp(var(2, i))
            ! write(111,'(5e15.7)') r_model(i),rho1D(i),p1D(i),Gamma1_1D(i),mass1D(i)
        enddo
        NN1D(1) = 0d0
        do i=2,nrmod
           NN1D(i)      = var(15,i)*mass1D(i)/r_model(i)**3
        enddo
        ! close(111)

    end subroutine

    !--------------------------------------------------------------------------
    ! This does a simple initialisation for the radial grid
    !--------------------------------------------------------------------------
    subroutine init_radial_grid_simple()

        integer i

#ifdef USE_1D
        nr = nrmod
        call init_radial_grid()
        do i=1, nr
            r(i) = var(1, i) / var(1, nr)
        enddo
#else
        grd(1)%nr = nrmod
        call init_radial_grid()
        do i=1, grd(1)%nr
            grd(1)%r(i) = var(1, i) / var(1, grd(1)%nr)
        enddo
#endif

    end subroutine init_radial_grid_simple

    !--------------------------------------------------------------------------
    ! This does a simple initialisation for the radial grid
    !--------------------------------------------------------------------------
    subroutine init_radial_grid_file()

        use inputs, only: gridfile
        integer i, j

        open(unit=37,file=trim(gridfile),status="old")
#ifdef USE_1D
        nr = 0
#else
        grd(1)%nr = 0
#endif
        do
            read(37,*,end=10)
#ifdef USE_1D
            nr = nr + 1
#else
            grd(1)%nr = grd(1)%nr + 1
#endif
        enddo
        10     rewind(37)
        print*,"New grid resolution:", grd(1)%nr
        grd(1)%nr = grd(1)%nr
        call init_radial_grid()
        do i=1, grd(1)%nr
            read(37,*) j, grd(1)%r(i)
        enddo
        close(37)

        ! very important
        do i=1,grd(1)%nr
            grd(1)%r(i) = grd(1)%r(i)/grd(1)%r(grd(1)%nr)
        enddo

    end subroutine init_radial_grid_file

    !--------------------------------------------------------------------------
    ! This subroutine calculates a radial grid appropriate for g-modes:
    !--------------------------------------------------------------------------
    subroutine init_radial_grid_g_modes()

        use inputs, C0_dati => C0, C1_dati => C1, C2_dati => C2, C3_dati => C3

        double precision, allocatable ::  PP(:), CC1(:), CC2(:), CC3(:), V_son(:)
        double precision :: P_total, P_target, mu, C1, C2, C3, C0
        integer i, j

        allocate(PP(nrmod), CC1(nrmod), CC2(nrmod), CC3(nrmod), V_son(nrmod))

        print*, "Grid type: ", trim(grid_type)
        if (trim(grid_type) == 'grid_g') then
            C0 = 1d0
            C1 = 2.5d-2
            C2 = 1d-1
            C3 = 1d-4
        else if (trim(grid_type) == 'grid_p') then
            C0 = 1d0
            C1 = 10d0
            C2 = 1d-2
            C3 = 1.5d-2
        else if (trim(grid_type) == 'grid_coef') then
            C0 = C0_dati
            C1 = C1_dati
            C2 = C2_dati
            C3 = C3_dati
            print*, "C0, C1, C2, C3 = ", C0, C1, C2, C3
        else
            stop "Unsupported grid in dati"
        end if

        ! find scaled Brunt-Vaisala frequency: N/r
        do i = 2, nrmod

            CC1(i) = rho1D(i)/Gamma1_1D(i)/p1D(i) ! Corresp au terme en C1
            CC2(i) = abs(NN1D(i))/r_model(i)**2
            CC3(i) = (mass1D(i)*rho1D(i)/(r_model(i)**2*p1D(i)))**2

        enddo

        CC1(1)=rho1D(1)/Gamma1_1D(1)/p1D(1)
        CC2(1)=0d0
        CC3(1)=0d0

        ! approximate time travel integral (from center), with
        ! possible saturation near surface.
        PP(1) = 0d0
        do i= 2, nrmod
            PP(i)=PP(i-1)+(r_model(i)-r_model(i-1))* &
                sqrt(C1*((CC1(i)+CC1(i-1))/2d0)+C2*((CC2(i)+CC2(i-1))/2d0)+ &
                C3*((CC3(i)+CC3(i-1))/2d0)+C0) ! 09/03
            ! print*, "PP: ", PP(i)
        enddo

        ! find total travel time from center to surface
        P_total = PP(nrmod)-PP(1)

        ! define new grid
        call init_radial_grid()
        grd(1)%r(1) = 0d0
        grd(1)%r(grd(1)%nr) = r_model(nrmod)

        j = 1
        do i=2, grd(1)%nr-1
            P_target = P_total*dble(i-1)/dble(grd(1)%nr-1)
            do while (PP(j).lt.P_target)
                j = j + 1
            enddo

            mu   = (PP(j)-P_target)/(PP(j)-PP(j-1))

            grd(1)%r(i) = r_model(j-1)*mu + r_model(j)*(1d0-mu)
        enddo

        do i=1, grd(1)%nr
            grd(1)%r(i) = grd(1)%r(i)/r_model(nrmod)
        enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! open(unit=3, file="info_diagram", status="unknown")
        do i=1, nrmod
            V_son(i)=Gamma1_1D(i)*p1D(i)/rho1D(i)
            ! write(3, *) r_model(i) , NN1D(i)*r_model(i)*r_model(i), V_son(i) !, sqrt(2*Gamma1_1D(i)*p1D(i)/rho1D(i)/(r_model(i)**2))
        enddo

        ! close(3)

        !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        deallocate(PP)
        deallocate(CC1)
        deallocate(CC3)
        deallocate(V_son)

        ! print*, NN1D(1)
    end subroutine init_radial_grid_g_modes

    !--------------------------------------------------------------------------
    ! This subroutine calculates a radial grid appropriate for p-modes:
    !--------------------------------------------------------------------------
#if 0
    subroutine init_radial_grid_p_modes()

        double precision, allocatable :: c_aux(:), TT(:)
        double precision :: c_max, T_total, T_target, mu
        integer i, j

        allocate(c_aux(nrmod),TT(nrmod))

        ! find sound velocity
        c_max = 0d0
        do i = 1, nrmod
            c_aux(i) = sqrt(Gamma1_1D(i)* p1D(i)/rho1D(i))
            if (c_aux(i).gt.c_max) c_max = c_aux(i)
        enddo

        ! approximate time travel integral (from center), with
        ! possible saturation near surface.
        TT(1) = 0d0
        do i= 2, nrmod
            if (c_aux(i).lt.(1d-4*c_max)) then
                c_aux(i) = 1d4/c_max
                TT(i) = TT(i-1) + (r_model(i)-r_model(i-1))*1d4/c_max
            else
                TT(i) = TT(i-1) + (r_model(i)-r_model(i-1))/c_aux(i)
            endif
        enddo

        ! find total travel time from center to surface
        T_total = TT(nrmod)-TT(1)

        ! define new grid
        call init_radial_grid()
        grd(1)%r(1) = 0d0
        grd(1)%r(grd(1)%nr) = r_model(nrmod)
        j = 1
        do i=2, grd(1)%nr-1
            T_target = T_total*dble(i-1)/dble(grd(1)%nr-1)
            do while (TT(j).lt.T_target)
                j = j + 1
            enddo
            mu   = (TT(j)-T_target)/(TT(j)-TT(j-1))
            grd(1)%r(i) = r_model(j-1)*mu + r_model(j)*(1d0-mu)
        enddo

        do i=1, grd(1)%nr
            grd(1)%r(i) = grd(1)%r(i)/r_model(nrmod)
        enddo

        deallocate(c_aux, TT)

    end subroutine init_radial_grid_p_modes
#endif

    !--------------------------------------------------------------------------
    ! This writes the radial grid as a function of index and stops the
    ! program.
    !--------------------------------------------------------------------------
    subroutine write_radial_grid()

        integer i

        open(unit=3,file="delme_grid",status="unknown")
        do i=1,grd(1)%nr
            write(3,*) i, grd(1)%r(i)
        enddo
        close(3)
        stop

    end subroutine write_radial_grid

    !--------------------------------------------------------------------------
    !  This finds rotational perturbations, assuming that the spherical
    !  symmetric component of the perturbation to the density profile is
    !  zero.
    !--------------------------------------------------------------------------
#ifndef USE_1D
    subroutine find_pert()

        type(DERMAT) :: dm
        integer, allocatable :: ipiv(:)
        double precision, allocatable :: mat(:,:), v(:), aux(:)
        integer kl, ku, ldmat, offset, info
        integer i, ii

        if (allocated(uhx)) deallocate(uhx)
        allocate(uhx(nrmod), aux(nrmod))
        call init_derive(dm, r_model,nrmod, 1, -1, 2,'IFD ')

        ! Find auxiliary function:
        do i=2, nrmod
            aux(i) = 2d0*Lambda*rho1D(i)*r_model(i)**3/mass1D(i)
        enddo
        aux(1) = 6d0

        ! initialise workspace and associated constants
        kl = 4 ; ku = 4; ldmat = 2*kl+ku+1; offset = kl+ku+1
        allocate(mat(ldmat, 2*nrmod), v(2*nrmod),ipiv(2*nrmod))

        ! Find uhx, where uhx is defined as follows:
        !       r(x,theta) = x - uhx(x) P_2(cos theta)
        mat = 0d0; v = 0d0;
        do i=1,nrmod-1
            do ii=max(1,i-dm%lbder(1)),min(nrmod,i+dm%ubder(1))
                mat(offset+2*(i-ii)+1,2*ii-1) = dm%derive(i,ii,1)*r_model(ii)**2 &
                    + dm%derive(i,ii,0)       &
                    * r_model(ii)*(aux(ii)-2d0)
                mat(offset+2*(i-ii),2*ii)     = dm%derive(i,ii,0)*(aux(ii)-6d0)
                mat(offset+2*(i-ii)+2,2*ii-1) =-dm%derive(i,ii,0)
                mat(offset+2*(i-ii)+1,2*ii)   = dm%derive(i,ii,1)
            enddo
        enddo
        ! Boundary condition (center):
        mat(offset, 1) = 1d0
        v(1) = 0d0
        ! Boundary condition (surface):
        mat(offset+1, 2*nrmod-1) = r_model(nrmod)
        mat(offset, 2*nrmod) = 2d0
        v(2*nrmod) = 5d0*r_model(nrmod)**3/(3d0*mass1D(nrmod))
        call DGBSV(2*nrmod, kl, ku, 1,mat, ldmat, ipiv, v, 2*nrmod, info)
        if (info.ne.0) then
            print*,"Problem with finding total potential. Info =", info
            stop
        endif
        do i=1,nrmod
            uhx(i) = v(2*i)*r_model(i)*rota_pert**2
        enddo

        call clear_derive(dm)
        deallocate(aux)

    end subroutine find_pert
#endif

    !--------------------------------------------------------------------------
    !  Set rotational perturbations to zero - i.e. no centrifugal deformation.
    !--------------------------------------------------------------------------
    subroutine no_pert()
        if (allocated(uhx))  deallocate(uhx)
        allocate(uhx(nrmod))
        uhx = 0d0
    end subroutine no_pert

    !--------------------------------------------------------------------------
    !  This prepares the mapping r(i,j) and other geometrical functions.
    !  It is necessary to interpolate in the angular direction - this is done
    !  by spectral methods involving Legendre polynomials.
    !--------------------------------------------------------------------------
#ifndef USE_1D
    subroutine make_mapping()

        double precision xi, epsilon, uprime
        double precision, allocatable, dimension(:) :: rs,rsp,rss,w
        integer i, ii, der, j, itemp

        ! this does a basic check on uhx to see if it
        ! will produce a singular mapping:
        do i=2,nrmod
            uprime = (uhx(i)-uhx(i-1))/(r_model(i)-r_model(i-1))
            if (uprime.ge.1) then
                print*,i
                stop "d(uhx)/dx >= 1"
            endif
            if (uprime.le.-2) then
                print*,i
                stop "d(uhx)/dx <= -2"
            endif
        enddo

        if (allocated(r_map))  deallocate(r_map)
        if (allocated(r_z))    deallocate(r_z)
        if (allocated(r_t))    deallocate(r_t)
        if (allocated(r_zz))   deallocate(r_zz)
        if (allocated(r_zt))   deallocate(r_zt)
        if (allocated(r_tt))   deallocate(r_tt)
        if (allocated(re_map)) deallocate(re_map)
        if (allocated(re_z))   deallocate(re_z)
        if (allocated(re_t))   deallocate(re_t)
        if (allocated(re_zz))  deallocate(re_zz)
        if (allocated(re_zt))  deallocate(re_zt)
        if (allocated(re_tt))  deallocate(re_tt)
        if (allocated(cth))    deallocate(cth)
        if (allocated(sth))    deallocate(sth)
        if (allocated(cost))   deallocate(cost)
        if (allocated(sint))   deallocate(sint)
        if (allocated(cott))   deallocate(cott)
        if (allocated(zeta))   deallocate(zeta)
        if (allocated(roz))    deallocate(roz)
        if (allocated(rrt))    deallocate(rrt)
        if (allocated(r_aux))  deallocate(r_aux)

        allocate(rs(lres),rsp(lres),rss(lres),  &
            r_map(grd(1)%nr,lres),r_z(grd(1)%nr,lres),   &
            r_t(grd(1)%nr,lres),r_zz(grd(1)%nr,lres),    &
            r_zt(grd(1)%nr,lres),r_tt(grd(1)%nr,lres),   &
            re_map(grd(1)%nr,lres),re_z(grd(1)%nr,lres), &
            re_t(grd(1)%nr,lres),re_zz(grd(1)%nr,lres),  & 
            re_zt(grd(1)%nr,lres),re_tt(grd(1)%nr,lres), &
            cth(lres), sth(lres), w(lres), &
            cost(grd(1)%nr,lres), sint(grd(1)%nr,lres),  &
            cott(grd(1)%nr,lres), zeta(grd(1)%nr,lres),  &
            r_aux(nrmod,lres),             &
            roz(grd(1)%nr,lres),rrt(grd(1)%nr,lres))
        rs = 0d0; rsp = 0d0; rss = 0d0;
        r_map = 0d0; r_z = 0d0; r_t = 0d0;
        r_zz = 0d0; r_zt = 0d0; r_tt = 0d0;
        re_map = 0d0; re_z = 0d0; re_t = 0d0;
        re_zz = 0d0; re_zt = 0d0; re_tt = 0d0;
        cth = 0d0; sth = 0d0; w = 0d0;
        cost = 0d0; sint = 0d0; cott = 0d0;
        zeta = 0d0; r_aux = 0d0; roz = 0d0;
        rrt = 0d0;

        do i=1,grd(1)%nr
            zeta(i,:) = grd(1)%r(i)
        enddo

        call gauleg(-1d0,1d0,cth,w,lres)
        sth = sqrt(1-cth**2)
        do i=1,grd(1)%nr
            sint(i,:) = sth(:)
            cost(i,:) = cth(:)
            cott(i,:) = cth(:)/sth(:)
        enddo

        do j=1,lres
            rs(j)   = r_model(nrmod) - (1.5d0*cth(j)**2-0.5d0)*uhx(nrmod)
            rsp(j)  = 3.0d0*sth(j)*cth(j)*uhx(nrmod)
            rss(j)  = (6.0d0*cth(j)**2-3d0)*uhx(nrmod)
        enddo

        K = r_model(nrmod) - uhx(nrmod)
        epsilon = 1d0-K

        ! Find the inner domain:
        do i=1,grd(1)%nr
            do j=1,lres
                r_map(i,j) = K*grd(1)%r(i) + 0.5d0*(5d0*grd(1)%r(i)**3-3d0*grd(1)%r(i)**5)*(rs(j)-K)
                r_z  (i,j) = K    + 0.5d0*(15d0*grd(1)%r(i)**2-15d0*grd(1)%r(i)**4)*(rs(j)-K)
                r_zz (i,j) =        0.5d0*(30d0*grd(1)%r(i)   -60d0*grd(1)%r(i)**3)*(rs(j)-K)
                r_t  (i,j) =        0.5d0*(5d0*grd(1)%r(i)**3-3d0*grd(1)%r(i)**5)*rsp(j)
                r_zt (i,j) =        0.5d0*(15d0*grd(1)%r(i)**2-15d0*grd(1)%r(i)**4)*rsp(j)
                r_tt (i,j) =        0.5d0*(5d0*grd(1)%r(i)**3-3d0*grd(1)%r(i)**5)*rss(j)
            enddo
        enddo

        ! Find the outer domain (we fold this domain onto the inner domain):
        ! Note: we're derive with respect to r(i), not xi, hence the sign change for
        ! re_z and re_zt.

        do i=1,grd(1)%nr
            xi = 2d0-grd(1)%r(i)
            do j=1,lres
                re_map(i,j) = 2d0*epsilon + K*xi + (2d0*xi**3-9d0*xi**2+12d0*xi-4d0)*(rs(j)-1d0-epsilon)
                re_z(i,j)   =             - K    - (6d0*xi**2-18d0*xi+12d0)*(rs(j)-1d0-epsilon)
                re_zz(i,j)  =                      (12d0*xi-18d0)*(rs(j)-1d0-epsilon)
                re_t(i,j)   =                      (2d0*xi**3-9d0*xi**2+12d0*xi-4d0)*rsp(j)
                re_zt(i,j)  =                    - (6d0*xi**2-18d0*xi+12d0)*rsp(j)
                re_tt(i,j)  =                      (2d0*xi**3-9d0*xi**2+12d0*xi-4d0)*rss(j)
            enddo
        enddo

        do i=2,nrmod-1
            do j=1,lres
                r_aux(i,j) = r_model(i) - (1.5d0*cth(j)**2-0.5d0)*uhx(i)
            enddo
        enddo

        do j=1,lres
            r_aux(1,j)     = r_map(1,j)
            r_aux(nrmod,j) = r_map(grd(1)%nr,j)
        enddo

        ! Supplementary geometric terms, which require special attention
        ! in the center:
        do j=1,lres
            do i=2,grd(1)%nr
                rrt(i,j) = r_t(i,j)/r_map(i,j)
                roz(i,j) = r_map(i,j)/zeta(i,j)
            enddo
            rrt(1,j) = r_zt(1,j)/r_z(1,j)
            roz(1,j) = r_z(1,j)
        enddo

        deallocate(rs, rsp, rss, w)
    end subroutine make_mapping
#endif

    !--------------------------------------------------------------
    !  This subroutine prepares different non-dimensional variables
    !  for the pulsation equations.
    !--------------------------------------------------------------
    subroutine init_fields()

        double precision, allocatable :: aux(:,:), p_aux(:), p1D_bis(:)
        integer i, j
#ifdef USE_1D
        type(DERMAT) :: dm
        integer :: order, der_max, der_min
        integer :: ii
        character*(4) :: dertype

       if (allocated(r_map)) deallocate(r_map)
       allocate(r_map(nr))
       call init_radial_grid()
       
       do j=1,nr
         r(j) = var(1,j)/glob(2)
         r_map(j) = r(j)
       enddo

       order = 1
       der_max = 1
       der_min = 0
       dertype = "FD  "
       call init_derive(dm,r,nr,der_max,der_min,order,dertype)

       if (allocated(rhom))       deallocate(rhom)
       if (allocated(rhom_z))     deallocate(rhom_z)
       if (allocated(diffrho_r2)) deallocate(diffrho_r2)
       if (allocated(pm))         deallocate(pm)
       if (allocated(pm_z))       deallocate(pm_z)
       if (allocated(c2))         deallocate(c2)
       if (allocated(g_m))        deallocate(g_m)
       if (allocated(dg_m))       deallocate(dg_m)
       if (allocated(gm_r))       deallocate(gm_r)
       if (allocated(Gamma1))     deallocate(Gamma1)
       if (allocated(NN))         deallocate(NN)
       if (allocated(rF))         deallocate(rF)

       ! calculate the different fields
       allocate(rhom(nr),rhom_z(nr),pm(nr),pm_z(nr),c2(nr),g_m(nr), &
                dg_m(nr),Gamma1(nr),NN(nr),gm_r(nr),rF(nr),         &
                diffrho_r2(nr))
       do i=1,nr
           rhom(i)   = var(5,i)/rho_ref
           pm(i)     = var(4,i)/p_ref
           c2(i)     = var(10,i)*pm(i)/rhom(i)
           Gamma1(i) = var(10,i)
       enddo
       do i=2,nr
           g_m(i)        = exp(var(2,i))/r_map(i)**2
           gm_r(i)       = exp(var(2,i))/r_map(i)**3
           diffrho_r2(i) = (1d0 - 4d0*pi*r_map(i)**3*rhom(i) &
                         / (3d0*exp(var(2,i))))/(r_map(i)**2)
       enddo
       g_m(1) = 0d0
       gm_r(1) = 4d0*pi*rhom(1)/3d0

       ! finding the radial derivative of various fields:
       do i=1,nr
         rhom_z(i) = 0d0
         pm_z(i)   = 0d0
         dg_m(i)    = 0d0
         do ii=max(1,i-dm%lbder(1)),min(nr,i+dm%ubder(1))
           rhom_z(i) = rhom_z(i) + dm%derive(i,ii,1)*rhom(ii)
           pm_z(i)   = pm_z(i)   + dm%derive(i,ii,1)*pm(ii)
           dg_m(i)   = dg_m(i)   + dm%derive(i,ii,1)*g_m(ii)
         enddo
       enddo

       rhom_z(1) = 0d0 ! it is extremely important to do this before
                       ! calculating the 2nd derivative
       rhom_zz = 0d0
       do ii=1,min(nr,1+dm%ubder(1))
         rhom_zz = rhom_zz + dm%derive(1,ii,1)*rhom_z(ii)
       enddo
       diffrho_r2(1) = -rhom_zz/(5d0*rhom(1))
       CC = 20d0*pi*rhom(1)**2/(3d0*c2(1)*rhom_zz)

       ! find the quantity NN which is related to the
       ! Brunt-Vaisala frequency:
       do i=2,nr
         NN(i) = (pm_z(i)-c2(i)*rhom_z(i))/r(i)
         rF(i) = 4d0*pi*rhom(i)*r_map(i)**3/(3d0*exp(var(2,i))) &
               - g_m(i)/(r_map(i)*diffrho_r2(i)*c2(i))          &
               - r_map(i)*rhom_z(i)/(2d0*rhom(i))
       enddo
       !NN(708:) = 0d0

       ! the center needs special treatment
       NN(1) = 2d0*(pm(2)-pm(1)-c2(1)*(rhom(2)-rhom(1)))/r(2)**2
       rF(1) = 1 + CC
#else
        if (allocated(rhom))    deallocate(rhom)
        if (allocated(rhom_z))  deallocate(rhom_z)
        if (allocated(rhom_t))  deallocate(rhom_t)
        if (allocated(pm))      deallocate(pm)
        if (allocated(pm_z))    deallocate(pm_z)
        if (allocated(pm_t))    deallocate(pm_t)
        if (allocated(pe))      deallocate(pe)
        if (allocated(pe_z))    deallocate(pe_z)
        if (allocated(pe_t))    deallocate(pe_t)
        if (allocated(Gamma1))  deallocate(Gamma1)
        if (allocated(c2))      deallocate(c2)
        if (allocated(NNt))     deallocate(NNt)
        if (allocated(NNr))     deallocate(NNr)
        if (allocated(NN))     deallocate(NN)

        allocate(rhom(grd(1)%nr,lres), rhom_z(grd(1)%nr,lres), rhom_t(grd(1)%nr,lres), &
            pm(grd(1)%nr,lres), pm_z(grd(1)%nr,lres), pm_t(grd(1)%nr,lres),       &
            Gamma1(grd(1)%nr,lres), aux(grd(1)%nr,lres), NNt(grd(1)%nr,lres),     &
            pe(grd(1)%nr,lres), pe_z(grd(1)%nr,lres), pe_t(grd(1)%nr,lres),       &
            c2(grd(1)%nr,lres),ws1(nrmod),ws2(grd(1)%nr), NNr(grd(1)%nr,lres),    &
            p_aux(nrmod), p1D_bis(nrmod), NN(grd(1)%nr, lres))

        ! The spherically averaged pressure is modifed as follows:
        do i=1,nrmod
        p_aux(i) = 2d0*r_model(i)*rho1D(i)/3d0
        enddo
        call antiderivative_down(r_model,p_aux,p1D_bis,nrmod,1)
        do i=1,nrmod
        p1D_bis(i) = rota_pert**2*(p1D_bis(i)-p1D_bis(nrmod))+p1D(i)
        enddo

        call map2D_der_p_bis(rho1D, rhom, rhom_t, rhom_z, aux)
        call map2D_der_p_bis(p1D_bis, pm, pm_t, pm_z, aux)
        call map2D(Gamma1_1D, Gamma1)
        call map2D(NN1D, NN)
        c2 = Gamma1*pm/rhom
        call map2D_der_bis(pe1D, pe, pe_t, pe_z, aux)

        ws1 = var(15,:) 
        call interpolate(r_aux(1:nrmod,1)**2,ws1,nrmod,r_map(1:grd(1)%nr,1)**2,ws2,grd(1)%nr)
        do i=2,grd(1)%nr
        do j=1,lres
        NNr(i,j) = Gamma1(i,j)*pm(i,j)*ws2(i)/r_map(i,j)
        enddo
        enddo
        NNr(1,:) = 0d0

        do j=1,lres
        do i=2,grd(1)%nr
        NNt(i,j) = zeta(i,j)*(pm_t(i,j)-Gamma1(i,j)*rhom_t(i,j))/(r_map(i,j)**2*r_z(i,j))
        enddo
        NNt(1,j) = 0d0
        enddo

        deallocate(aux,ws1,ws2,p_aux,p1D_bis)
#endif

    end subroutine

    !--------------------------------------------------------------------------
    ! This subroutine initialises all of the grd_pe arrays in such a way
    ! that there are no singularities in the center of the star.
    !--------------------------------------------------------------------------
#ifndef USE_1D
    subroutine init_grd_pe()

        type(DERMAT) :: dm
        double precision, allocatable :: aux(:,:)
        integer der_min, i, ii, j

        if ((trim(dertype).eq."IFD").or.(trim(dertype).eq."SMPL")) then
            der_min = -1
        else
            der_min = 0
        endif
        call init_derive(dm,grd(1)%r,grd(1)%nr,1,der_min,orderFD,dertype)

        if (allocated(grd_pe_z))   deallocate(grd_pe_z)
        if (allocated(grd_pe_t))   deallocate(grd_pe_t)
        if (allocated(grd_pe_zz))  deallocate(grd_pe_zz)
        if (allocated(grd_pe_zt))  deallocate(grd_pe_zt)
        if (allocated(grd_pe_tz))  deallocate(grd_pe_tz)
        if (allocated(grd_pe_tt))  deallocate(grd_pe_tt)
        allocate(grd_pe_z(grd(1)%nr,lres), grd_pe_t(grd(1)%nr,lres),   &
            grd_pe_zz(grd(1)%nr,lres),grd_pe_zt(grd(1)%nr,lres),  &
            grd_pe_tz(grd(1)%nr,lres),grd_pe_tt(grd(1)%nr,lres),  &
            aux(grd(1)%nr,lres))

        grd_pe_z = pe_z/(roz**2*r_z)
        do j=1,lres
        do i=2,grd(1)%nr
        grd_pe_t(i,j) = zeta(i,j)*pe_t(i,j)/(r_map(i,j)**2*r_z(i,j))
        enddo
        grd_pe_t(1,j) = 0d0
        enddo

        ! Numerically calculate zeta derivatives: this avoids various
        ! difficulties with singularities in the center of the star:
        ! IMPORTANT: if dertype = IFD or SMPL, this produces derivatives
        ! at mid (or nearly mid) grid points rather than on the original
        ! grid.
        grd_pe_zz = 0d0
        grd_pe_zt = 0d0
        do j=1,lres
        do i=1,grd(1)%nr
        do ii=max(1,i-dm%lbder(1)),min(grd(1)%nr,i+dm%ubder(1))
        grd_pe_zz(i,j) = grd_pe_zz(i,j) + dm%derive(i,ii,1)*grd_pe_z(ii,j)
        grd_pe_zt(i,j) = grd_pe_zt(i,j) + dm%derive(i,ii,1)*grd_pe_t(ii,j)
        enddo
        enddo
        enddo

        ! Numerically calculate theta derivatives
        aux = 0d0
        call legendre(grd_pe_z,aux,lres,grd(1)%nr,0,1)
        call legendrep(aux,grd_pe_tz,lres,grd(1)%nr,0)
        ! This is not the most elegant solution but the legendre grid
        ! does not go to theta = 0, so this should work (at reasonable
        ! resolutions)
        aux = 0d0
        call legendre(grd_pe_t/sint,aux,lres,grd(1)%nr,0,1)
        call legendrep(aux,grd_pe_tt,lres,grd(1)%nr,0)
        grd_pe_tt = sint*grd_pe_tt+cott*grd_pe_t

        call clear_derive(dm)
        deallocate(aux)

        if ((trim(dertype).ne."IFD").and. &
            (trim(dertype).ne."SMPL")) return

        ! the following sets the radial derivatives grd_pe_zz and
        ! grd_pe_zt at the last point (nr):
        der_min = 0
        call init_derive(dm,grd(1)%r,grd(1)%nr,1,der_min,orderFD,"FD  ")
        do j=1,lres
        grd_pe_zz(grd(1)%nr,j) = 0d0
        grd_pe_zt(grd(1)%nr,j) = 0d0
        do ii=max(1,grd(1)%nr-dm%lbder(1)),min(grd(1)%nr,i+dm%ubder(1))
        grd_pe_zz(grd(1)%nr,j) = grd_pe_zz(grd(1)%nr,j) + dm%derive(grd(1)%nr,ii,1)*grd_pe_z(ii,j)
        grd_pe_zt(grd(1)%nr,j) = grd_pe_zt(grd(1)%nr,j) + dm%derive(grd(1)%nr,ii,1)*grd_pe_t(ii,j)
        enddo
        enddo

        call clear_derive(dm)

    end subroutine init_grd_pe
#endif

    !--------------------------------------------------------------------------
    ! This writes an amdl which can be used by ADIPLS.
    !--------------------------------------------------------------------------
#ifndef USE_1D
    subroutine write_amdl()

        double precision dta(8), mu
        double precision, allocatable :: aa(:,:)
        integer i, j

        allocate(aa(5,grd(1)%nr))

        ! Since we're using 2D variables to write a 1D model, a specific
        ! latitude must be chose (as expressed by j):
        j = 1

        ! write binary model for adipls
        dta(1) = mass*solar_mass
        dta(2) = radius
        dta(3) = pm(1,j)*p_ref
        dta(4) = rhom(1,j)*rho_ref
        dta(5) = -glob(9)/Gamma1(1,j)
        dta(6) = -glob(10)
        dta(7) = -1d0
        dta(8) = 0d0

        do i=2,grd(1)%nr
        !aa(1,i) = -pm_z(i,j)/(r_map(i,j)*rhom(i,j))
        !aa(2,i) = -r_map(i,j)*pm_z(i,j)/(Gamma1(i,j)*pm(i,j))
        !aa(3,i) = Gamma1(i,j)
        !aa(4,i) = -aa(2,i)-r_map(i,j)*rhom_z(i,j)/rhom(i,j)
        !aa(5,i) = -4d0*pi*rhom(i,j)**2*r_map(i,j)/pm_z(i,j)
        aa(1,i) = pe_z(i,j)/r_map(i,j)
        aa(2,i) = r_map(i,j)*pe_z(i,j)*rhom(i,j)/(Gamma1(i,j)*pm(i,j))
        aa(3,i) = Gamma1(i,j)
        aa(4,i) = -aa(2,i)-r_map(i,j)*rhom_z(i,j)/rhom(i,j)
        aa(5,i) = 4d0*pi*rhom(i,j)*r_map(i,j)/pe_z(i,j)
        enddo
        !call interpolate(r_aux(1:nrmod,j)**2,var(15,1:nrmod),nrmod, &
        !                 r_map(1:grd(1)%nr,j)**2,aa(4,1:grd(1)%nr),grd(1)%nr)

        ! central values:
        aa(1,1) = 4d0*pi*rhom(1,j)/3d0
        aa(2,1) = 0d0
        aa(3,1) = Gamma1(1,j)
        aa(4,1) = 0d0
        aa(5,1) = 3d0

        open(unit=3,file="amdl",status="unknown",form="unformatted")
        write(3) 0,grd(1)%nr,(dta(j),j=1,8),(r_map(i,1), (aa(j,i),j=1,5),i=1,grd(1)%nr)
        close(3)

        deallocate(aa)
        stop

    end subroutine write_amdl
#endif

    !--------------------------------------------------------------------------
    ! This subroutine calculates pe1D using Poisson's equation.
    !--------------------------------------------------------------------------
    subroutine find_pe1D()

        double precision, allocatable, dimension(:)  :: f, ff
        integer i

        if (allocated(pe1D)) deallocate(pe1D)
        allocate(f(nrmod),ff(nrmod),pe1D(nrmod))

        do i=1,nrmod
        f(i) = -4d0*pi*r_model(i)**2*rho1D(i)
        enddo

        call antiderivative_down(r_model,f,ff,nrmod,1)

        pe1D(1) = 0d0
        do i=2,nrmod
        pe1D(i) = (ff(i)-ff(1))/r_model(i)
        enddo 

        do i=1,nrmod
        f(i) = -4d0*pi*r_model(i)*rho1D(i)
        enddo

        call antiderivative_down(r_model,f,ff,nrmod,1)

        do i=1,nrmod
        pe1D(i) = pe1D(i) + (ff(nrmod) - ff(i)) &
            - (rota_pert*r_model(i))**2/3d0
        enddo

        deallocate(f,ff)
    end subroutine find_pe1D

    !--------------------------------------------------------------------------
    ! This subroutine calculates pe1D by integrating the local gravity.
    !--------------------------------------------------------------------------

    subroutine find_pe1D_alt()

        double precision, allocatable, dimension(:)  :: aux
        integer i

        if (allocated(pe1D)) deallocate(pe1D)
        allocate(aux(nrmod),pe1D(nrmod))

        do i=2,nrmod
        aux(i) = mass1D(i)/r_model(i)**2
        enddo
        aux(1) = 0d0
        call antiderivative_up(r_model,aux,pe1D,nrmod,1)
        do i=1,nrmod
        pe1D(i) = pe1D(i) - pe1D(nrmod) - mass1D(nrmod)/r_model(nrmod) &
            - (rota_pert*r_model(i))**2/3d0
        enddo

        deallocate(aux)

    end subroutine find_pe1D_alt

    !--------------------------------------------------------------------------
    !  This subroutine maps a field from the iso-potentials grid to a
    !  new mapping, which introduces dependance of theta.
    !--------------------------------------------------------------------------
#ifndef USE_1D
    subroutine map2D(f1D,f2D)

        double precision f1D(nrmod), f2D(grd(1)%nr,lres)
        integer i, j

        ! make use of equatorial symmetry
        do j=1,(lres+1)/2
        call interpolate(r_aux(1:nrmod,j)**2,f1D,nrmod,r_map(1:grd(1)%nr,j)**2,f2D(1:grd(1)%nr,j),grd(1)%nr)
        f2D(1:grd(1)%nr,lres+1-j) = f2D(1:grd(1)%nr,j)
        enddo

    end subroutine map2D
#endif

    !--------------------------------------------------------------------------
    !  This subroutine maps a field from the iso-potentials grid to a
    !  new mapping, which introduces dependance of theta.  It also calculates
    !  the theta and zeta derivatives on the new mapping.
    !--------------------------------------------------------------------------
#ifndef USE_1D
    subroutine map2D_der(f1D,f2D,f2D_t,f2D_z,aux)

        double precision f1D(nrmod), aux(grd(1)%nr,lres)
        double precision, dimension(grd(1)%nr,lres) :: f2D, f2D_t, f2D_z
        integer i, j

        ! make use of equatorial symmetry
        do j=1,(lres+1)/2
        call interpolate(r_aux(1:nrmod,j)**2,f1D,nrmod,r_map(1:grd(1)%nr,j)**2,f2D(1:grd(1)%nr,j),grd(1)%nr)
        call interpolate_derive(r_aux(1:nrmod,j)**2,f1D,nrmod,r_map(1:grd(1)%nr,j)**2,f2D_z(1:grd(1)%nr,j),grd(1)%nr)
        f2D(1:grd(1)%nr,lres+1-j)   = f2D(1:grd(1)%nr,j)
        f2D_z(1:grd(1)%nr,lres+1-j) = f2D_z(1:grd(1)%nr,j)
        enddo

        aux = 0d0
        call legendre(f2D(1:grd(1)%nr,1:lres),aux(1:grd(1)%nr,1:lres),lres,grd(1)%nr,0,1)
        call legendrep(aux(1:grd(1)%nr,1:lres),f2D_t(1:grd(1)%nr,1:lres),lres,grd(1)%nr,0)

        f2D_z = 2d0*f2D_z*r_z*r_map

    end subroutine
#endif

    !--------------------------------------------------------------------------
    !  This subroutine maps a field from the iso-potentials grid to a new
    !  mapping, which introduces a dependance on theta.  It also calculates
    !  the theta and zeta derivatives on the new mapping.
    !--------------------------------------------------------------------------
#ifndef USE_1D
    subroutine map2D_der_bis(f1D,f2D,f2D_t,f2D_z,aux)

        double precision f1D(nrmod), aux(grd(1)%nr,lres)
        double precision, dimension(grd(1)%nr,lres) :: f2D, f2D_t, f2D_z
        integer i, ii, j, order
        type(DERMAT) :: dm

        order = 2
        call init_derive(dm,grd(1)%r,grd(1)%nr,1,0,order,"FD  ")

        ! make use of equatorial symmetry
        do j=1,(lres+1)/2
        call interpolate(r_aux(1:nrmod,j)**2,f1D,nrmod,r_map(1:grd(1)%nr,j)**2,f2D(1:grd(1)%nr,j),grd(1)%nr)
        do i=1,grd(1)%nr
        f2D_z(i,j) = 0d0
        do ii=max(1,i-dm%lbder(1)),min(grd(1)%nr,i+dm%ubder(1))
        f2D_z(i,j) = f2D_z(i,j) + dm%derive(i,ii,1)*f2D(ii,j)
        enddo
        enddo
        f2D(1:grd(1)%nr,lres+1-j)   = f2D(1:grd(1)%nr,j)
        f2D_z(1:grd(1)%nr,lres+1-j) = f2D_z(1:grd(1)%nr,j)
        enddo

        aux = 0d0
        call legendre(f2D(1:grd(1)%nr,1:lres),aux(1:grd(1)%nr,1:lres),lres,grd(1)%nr,0,1)
        call legendrep(aux(1:grd(1)%nr,1:lres),f2D_t(1:grd(1)%nr,1:lres),lres,grd(1)%nr,0)

        !f2D_z = 2d0*f2D_z*r_z*r_map
        call clear_derive(dm)

    end subroutine map2D_der_bis
#endif

    !--------------------------------------------------------------------------
    !  This subroutine maps a field from the iso-potentials grid to a
    !  new mapping, which introduces dependance of theta.
    !--------------------------------------------------------------------------
#ifndef USE_1D
    subroutine map2D_p(f1D,f2D)

        use mod_legendre

        double precision f1D(nrmod), f2D(grd(1)%nr,lres)
        integer i, j

        ws1 = log(f1D)

        ! make use of equatorial symmetry
        do j=1,(lres+1)/2
        call interpolate(r_aux(1:nrmod,j)**2,ws1,nrmod,r_map(1:grd(1)%nr,j)**2,ws2,grd(1)%nr)
        f2D(1:grd(1)%nr,j) = exp(ws2)
        f2D(1:grd(1)%nr,lres+1-j) = f2D(1:grd(1)%nr,j)
        enddo

    end subroutine
#endif

    !--------------------------------------------------------------------------
    !  This subroutine maps a field from the iso-potentials grid to a
    !  new mapping, which introduces dependance of theta.  It also calculates
    !  the theta and zeta derivatives on the new mapping.
    !--------------------------------------------------------------------------
#ifndef USE_1D
    subroutine map2D_der_p(f1D,f2D,f2D_t,f2D_z,aux)

        use mod_legendre

        double precision f1D(nrmod)
        double precision, dimension(grd(1)%nr,lres) :: f2D, f2D_t, f2D_z, aux
        integer i, j

        ws1 = log(f1D)

        ! make use of equatorial symmetry
        do j=1,(lres+1)/2
        call interpolate(r_aux(1:nrmod,j)**2,ws1,nrmod,r_map(1:grd(1)%nr,j)**2,ws2,grd(1)%nr)
        f2D(1:grd(1)%nr,j) = exp(ws2)
        f2D(1:grd(1)%nr,lres+1-j)   = f2D(1:grd(1)%nr,j)
        call interpolate_derive(r_aux(1:nrmod,j)**2,ws1,nrmod,r_map(1:grd(1)%nr,j)**2,ws2,grd(1)%nr)
        f2D_z(1:grd(1)%nr,j) = f2D(1:grd(1)%nr,j)*ws2(1:grd(1)%nr)
        f2D_z(1:grd(1)%nr,lres+1-j) = f2D_z(1:grd(1)%nr,j)
        enddo

        aux = 0d0
        call legendre(f2D(1:grd(1)%nr,1:lres),aux(1:grd(1)%nr,1:lres),lres,grd(1)%nr,0,1)
        call legendrep(aux(1:grd(1)%nr,1:lres),f2D_t(1:grd(1)%nr,1:lres),lres,grd(1)%nr,0)

        f2D_z = 2d0*f2D_z*r_z*r_map

    end subroutine
#endif

    !--------------------------------------------------------------------------
    !  This subroutine maps a field from the iso-potentials grid to a
    !  new mapping, which introduces dependance of theta.  It also calculates
    !  the theta and zeta derivatives on the new mapping.
    !--------------------------------------------------------------------------
#ifndef USE_1D
    subroutine map2D_der_p_bis(f1D,f2D,f2D_t,f2D_z,aux)

        use mod_legendre
        use derivative

        double precision f1D(nrmod)
        double precision, dimension(grd(1)%nr,lres) :: f2D, f2D_t, f2D_z, aux
        integer i, ii, j, order
        type(DERMAT) :: dm

        order = 2
        call init_derive(dm, grd(1)%r, grd(1)%nr, 1, 0, order, "FD  ")

        ws1 = log(f1D)

        ! make use of equatorial symmetry
        do j=1,(lres+1)/2
            call interpolate(r_aux(1:nrmod,j)**2,ws1,nrmod,r_map(1:grd(1)%nr,j)**2,ws2,grd(1)%nr)
            f2D(1:grd(1)%nr,j) = exp(ws2)
            f2D(1:grd(1)%nr,lres+1-j)   = f2D(1:grd(1)%nr,j)
            do i=1,grd(1)%nr
            f2D_z(i,j) = 0d0
            do ii=max(1,i-dm%lbder(1)),min(grd(1)%nr,i+dm%ubder(1))
                f2D_z(i,j) = f2D_z(i,j) + dm%derive(i,ii,1)*f2D(ii,j)
            enddo
        enddo
        f2D_z(1:grd(1)%nr,lres+1-j) = f2D_z(1:grd(1)%nr,j)
        enddo

        aux = 0d0
        call legendre(f2D(1:grd(1)%nr,1:lres),aux(1:grd(1)%nr,1:lres),lres,grd(1)%nr,0,1)
        call legendrep(aux(1:grd(1)%nr,1:lres),f2D_t(1:grd(1)%nr,1:lres),lres,grd(1)%nr,0)

        !f2D_z = 2d0*f2D_z*r_z*r_map
        call clear_derive(dm)

    end subroutine map2D_der_p_bis
#endif

    !------------------------------------------------------------------------------ 
    ! This subroutine calculates the anti-derivative of a function f defined on an
    ! arbitrary strictly monotonic grid.  It does this by integrating Lagrange
    ! interpolation polynomials calculated over a sliding window which spans
    ! [i-wnd, i+wnd+1].
    !------------------------------------------------------------------------------ 
    ! description of variables:
    !
    ! grid(1:ngrid) = grid on which the function is defined
    ! f(1:ngrid)    = input function
    ! ff(1:ngrid)   = output anti-derivative
    ! ngrid         = number of grid points
    ! wnd           = positive integer which gives half the size of the window.
    !------------------------------------------------------------------------------ 
    subroutine antiderivative_down(grid,f,ff,ngrid,wnd)

        integer, intent(in) :: ngrid, wnd
        double precision, intent(in) :: grid(ngrid), f(ngrid)
        double precision, intent(out):: ff(ngrid)
        integer i, j, k, l, start, finish
        double precision, allocatable :: mu(:), a(:)
        double precision prdct, my_sum

        ! check parameters:
        if (wnd.lt.0) stop "Please set wnd>= 0 in find_weights"

        ! initialise ff:
        ff(ngrid) = 0d0

        allocate(mu(-wnd:(wnd+1)),a(-1:2*(wnd+1)))

        do i=ngrid-1,1,-1

        ff(i) = ff(i+1)

        start  = -wnd
        finish = wnd+1

        ! special treatment do the endpoints
        ! Note: this leads to results which are less precise on the edges
        if ((i+start).lt.1) start = 1 - i
        if ((i+finish).gt.ngrid) finish = ngrid - i

        do j=start,finish
        mu(j) = grid(j+i) - grid(i)
        enddo

        do j=start,finish

        ! Initialise Lagrange polynomial
        a(-1) = 0d0
        a(0) = 1d0
        do l=1,2*(wnd+1)
        a(l) = 0d0
        enddo

        prdct = 1d0
        do k=start,j-1
        prdct = prdct*(mu(j)-mu(k))
        enddo
        do k=j+1,finish
        prdct = prdct*(mu(j)-mu(k))
        enddo
        a(0) = a(0)/prdct


        ! Calculate Lagrange polynomial, by calculating product of (x-mu(k))
        do k=start,j-1
        do l=2*(wnd+1),0,-1
        a(l) = -mu(k)*a(l) + a(l-1)
        enddo
        enddo
        do k=j+1,finish
        do l=2*(wnd+1),0,-1
        a(l) = -mu(k)*a(l) + a(l-1)
        enddo
        enddo

        ! Integrate the Lagrange polynomials, and add the results
        ! as weights at the appropriate points.

        do l=0,2*(wnd+1)
        ff(i) = ff(i) - f(i+j)*a(l)*mu(1)**(l+1)/dble(l+1)
        enddo

        enddo
        enddo

        deallocate(a,mu)

    end subroutine antiderivative_down

    !------------------------------------------------------------------------------ 
    ! This subroutine calculates the anti-derivative of a function f defined on an
    ! arbitrary strictly monotonic grid.  It does this by integrating Lagrange
    ! interpolation polynomials calculated over a sliding window which spans
    ! [i-wnd, i+wnd+1].
    !------------------------------------------------------------------------------ 
    ! description of variables:
    !
    ! grid(1:ngrid) = grid on which the function is defined
    ! f(1:ngrid)    = input function
    ! ff(1:ngrid)   = output anti-derivative
    ! ngrid         = number of grid points
    ! wnd           = positive integer which gives half the size of the window.
    !------------------------------------------------------------------------------ 
    subroutine antiderivative_up(grid,f,ff,ngrid,wnd)

        integer, intent(in) :: ngrid, wnd
        double precision, intent(in) :: grid(ngrid), f(ngrid)
        double precision, intent(out):: ff(ngrid)
        integer i, j, k, l, start, finish
        double precision, allocatable :: mu(:), a(:)
        double precision prdct, my_sum

        ! check parameters:
        if (wnd.lt.0) stop "Please set wnd>= 0 in find_weights"

        ! initialise ff:
        ff(1) = 0d0

        allocate(mu(-wnd:(wnd+1)),a(-1:2*(wnd+1)))

        do i=1,ngrid-1

        ff(i+1) = ff(i)

        start  = -wnd
        finish = wnd+1

        ! special treatment do the endpoints
        ! Note: this leads to results which are less precise on the edges
        if ((i+start).lt.1) start = 1 - i
        if ((i+finish).gt.ngrid) finish = ngrid - i

        do j=start,finish
        mu(j) = grid(j+i) - grid(i)
        enddo

        do j=start,finish

        ! Initialise Lagrange polynomial
        a(-1) = 0d0
        a(0) = 1d0
        do l=1,2*(wnd+1)
        a(l) = 0d0
        enddo

        prdct = 1d0
        do k=start,j-1
        prdct = prdct*(mu(j)-mu(k))
        enddo
        do k=j+1,finish
        prdct = prdct*(mu(j)-mu(k))
        enddo
        a(0) = a(0)/prdct


        ! Calculate Lagrange polynomial, by calculating product of (x-mu(k))
        do k=start,j-1
        do l=2*(wnd+1),0,-1
        a(l) = -mu(k)*a(l) + a(l-1)
        enddo
        enddo
        do k=j+1,finish
        do l=2*(wnd+1),0,-1
        a(l) = -mu(k)*a(l) + a(l-1)
        enddo
        enddo

        ! Integrate the Lagrange polynomials, and add the results
        ! as weights at the appropriate points.

        do l=0,2*(wnd+1)
        ff(i+1) = ff(i+1) + f(i+j)*a(l)*mu(1)**(l+1)/dble(l+1)
        enddo

        enddo
        enddo

        deallocate(a,mu)

    end subroutine antiderivative_up

    !--------------------------------------------------------------------------
    !  This subroutine removes half the points from the model so as to do
    !  Richardson extrapolation
    !--------------------------------------------------------------------------
#ifndef USE_1D
    subroutine reduce_model()

        integer i, nr_in, nr_out
        double precision, allocatable :: aux(:,:), aux1D(:)

        nr_in  = grd(1)%nr
        nr_out = 1 + floor(dble(grd(1)%nr)/2d0)
        grd(1)%nr = nr_out
        allocate(aux(nr_in,lres),aux1D(nr_in))

        aux = rhom   
        deallocate(rhom   )
        allocate(rhom   (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,rhom   ,nr_out)
        aux = rhom_z 
        deallocate(rhom_z )
        allocate(rhom_z (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,rhom_z ,nr_out)
        aux = rhom_t 
        deallocate(rhom_t )
        allocate(rhom_t (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,rhom_t ,nr_out)
        aux = pm     
        deallocate(pm     )
        allocate(pm     (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,pm     ,nr_out)
        aux = pm_z   
        deallocate(pm_z   )
        allocate(pm_z   (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,pm_z   ,nr_out)
        aux = pm_t   
        deallocate(pm_t   )
        allocate(pm_t   (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,pm_t   ,nr_out)
        aux = c2     
        deallocate(c2     )
        allocate(c2     (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,c2     ,nr_out)
        aux = Gamma1     
        deallocate(Gamma1     )
        allocate(Gamma1     (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,Gamma1     ,nr_out)
        aux = NNt    
        deallocate(NNt    )
        allocate(NNt    (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,NNt    ,nr_out)
        aux = r_map  
        deallocate(r_map  )
        allocate(r_map  (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,r_map  ,nr_out)
        aux = r_z    
        deallocate(r_z    )
        allocate(r_z    (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,r_z    ,nr_out)
        aux = r_t    
        deallocate(r_t    )
        allocate(r_t    (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,r_t    ,nr_out)
        aux = r_zz   
        deallocate(r_zz   )
        allocate(r_zz   (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,r_zz   ,nr_out)
        aux = r_zt   
        deallocate(r_zt   )
        allocate(r_zt   (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,r_zt   ,nr_out)
        aux = r_tt   
        deallocate(r_tt   )
        allocate(r_tt   (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,r_tt   ,nr_out)
        aux = re_map 
        deallocate(re_map )
        allocate(re_map (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,re_map ,nr_out)
        aux = re_z   
        deallocate(re_z   )
        allocate(re_z   (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,re_z   ,nr_out)
        aux = re_t   
        deallocate(re_t   )
        allocate(re_t   (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,re_t   ,nr_out)
        aux = re_zz  
        deallocate(re_zz  )
        allocate(re_zz  (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,re_zz  ,nr_out)
        aux = re_zt  
        deallocate(re_zt  )
        allocate(re_zt  (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,re_zt  ,nr_out)
        aux = re_tt  
        deallocate(re_tt  )
        allocate(re_tt  (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,re_tt  ,nr_out)
        aux = zeta   
        deallocate(zeta   )
        allocate(zeta   (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,zeta   ,nr_out)
        aux = cost   
        deallocate(cost   )
        allocate(cost   (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,cost   ,nr_out)
        aux = sint   
        deallocate(sint   )
        allocate(sint   (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,sint   ,nr_out)
        aux = cott   
        deallocate(cott   )
        allocate(cott   (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,cott   ,nr_out)
        aux = roz    
        deallocate(roz    )
        allocate(roz    (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,roz    ,nr_out)
        aux = rrt    
        deallocate(rrt    )
        allocate(rrt    (nr_out,lres))
        call reduce_2Darray(aux,nr_in,lres,rrt    ,nr_out)
        aux1D = r_model(1:nr_in)
        deallocate(r_model)
        allocate(r_model(nr_out))
        call reduce_1Darray(aux1D,nr_in,r_model,nr_out)
        aux1D = grd(1)%r(1:nr_in)
        deallocate(grd(1)%r)
        call init_radial_grid()
        call reduce_1Darray(aux1D,nr_in,grd(1)%r,nr_out)
        aux1D = rho1D  
        deallocate(rho1D  )
        allocate(rho1D  (nr_out))
        call reduce_1Darray(aux1D,nr_in,rho1D  ,nr_out)
        aux1D = p1D    
        deallocate(p1D    )
        allocate(p1D    (nr_out))
        call reduce_1Darray(aux1D,nr_in,p1D    ,nr_out)
        aux1D = Gamma1_1D  
        deallocate(Gamma1_1D  )
        allocate(Gamma1_1D  (nr_out))
        call reduce_1Darray(aux1D,nr_in,Gamma1_1D  ,nr_out)
        aux1D = uhx    
        deallocate(uhx    )
        allocate(uhx    (nr_out))
        call reduce_1Darray(aux1D,nr_in,uhx    ,nr_out)

        deallocate(aux,aux1D)

    end subroutine reduce_model
    !--------------------------------------------------------------------------
    !  This copies half of 2D array f_in into f_out.
    !--------------------------------------------------------------------------
    subroutine reduce_2Darray(f_in,nr_in,llres,f_out,nr_out)

        integer, intent(in) ::  nr_in, nr_out, llres
        double precision, intent(in) :: f_in(nr_in,llres)
        double precision, intent(out) :: f_out(nr_out,llres)
        integer i, ii, j

        if (nr_out.ne.(1 + floor(dble(nr_in)/2d0))) then
            stop "nr_out has a faulty value in reduce_2Darray"
        endif

        ii = 1
        do i=1,nr_in,2
        do j=1,llres
        f_out(ii,j) = f_in(i,j)
        enddo
        ii = ii + 1
        enddo
        do j=1,llres
        f_out(nr_out,j) = f_in(nr_in,j)
        enddo

    end subroutine
    !--------------------------------------------------------------------------
    !  This copies half of 1D array f_in into f_out.
    !--------------------------------------------------------------------------
    subroutine reduce_1Darray(f_in,nr_in,f_out,nr_out)

        integer, intent(in) ::  nr_in, nr_out
        double precision, intent(in) :: f_in(nr_in)
        double precision, intent(out) :: f_out(nr_out)
        integer i, ii

        if (nr_out.ne.(1 + floor(dble(nr_in)/2d0))) then
            stop "nr_out has a faulty value in reduce_2Darray"
        endif

        ii = 1
        do i=1,nr_in,2
            f_out(ii) = f_in(i)
            ii = ii + 1
        enddo
        f_out(nr_out) = f_in(nr_in)

    end subroutine
#endif
    !-------------------------------------------------------------
    !  This subroutine prints different fields so as to help debug
    !  the program...
    !-------------------------------------------------------------
#ifndef USE_1D
    subroutine write_fields()
        integer i,j

        open(unit=999,file="delme_grid",status="unknown")
        write(999,*) grd(1)%nr,lres
        write(999,'(a)') "# x"
        write(999,101) ((r_map(i,j)*sint(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# y"
        write(999,101) ((r_map(i,j)*cost(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# xe"
        write(999,101) ((re_map(i,j)*sint(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# ye"
        write(999,101) ((re_map(i,j)*cost(i,j),i=1,grd(1)%nr),j=1,lres)
        close(999)

        open(unit=999,file="delme_fields",status="unknown")
        write(999,*) grd(1)%nr,lres
        write(999,'(a)') "# rhom"
        write(999,101) ((rhom(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# rhom_z"
        write(999,101) ((rhom_z(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# rhom_t"
        write(999,101) ((rhom_t(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# pm"
        write(999,101) ((pm(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# pm_z"
        write(999,101) ((pm_z(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# pm_t"
        write(999,101) ((pm_t(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# Gamma1"
        write(999,101) ((Gamma1(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# c2"
        write(999,101) ((c2(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# NNt"
        write(999,101) ((NNt(i,j),i=1,grd(1)%nr),j=1,lres)
        close(999)

        open(unit=999,file="delme_geometry",status="unknown")
        write(999,*) grd(1)%nr,lres
        write(999,'(a)') "# r_map"
        write(999,101) ((r_map(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# r_z"
        write(999,101) ((r_z(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# r_t"
        write(999,101) ((r_t(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# r_zz"
        write(999,101) ((r_zz(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# r_zt"
        write(999,101) ((r_zt(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# r_tt"
        write(999,101) ((r_tt(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# re_map"
        write(999,101) ((re_map(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# re_z"
        write(999,101) ((re_z(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# re_t"
        write(999,101) ((re_t(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# re_zz"
        write(999,101) ((re_zz(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# re_zt"
        write(999,101) ((re_zt(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# re_tt"
        write(999,101) ((re_tt(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# roz"
        write(999,101) ((roz(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# rrt"
        write(999,101) ((rrt(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# zeta"
        write(999,101) ((zeta(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# cost"
        write(999,101) ((cost(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# sint"
        write(999,101) ((sint(i,j),i=1,grd(1)%nr),j=1,lres)
        write(999,'(a)') "# cott"
        write(999,101) ((cott(i,j),i=1,grd(1)%nr),j=1,lres)
        close(999)

        open(unit=999,file="delme_constants",status="unknown")
        write(999,'(a)') "# Lambda"
        write(999,'(1pe22.15)') Lambda
        write(999,'(a)') "# Rota"
        write(999,'(1pe22.15)') rota_pert
        close(999)

        101    format(3(1pe22.15,2X))
        stop
    end subroutine write_fields
#endif
    !--------------------------------------------------------------------------
end module
