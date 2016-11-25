#include "config.h"
!------------------------------------------------------------------------------
!  This module contains all the matrices relevant to the oscillations problem.
!  It is based on the idea of additive operators, so as to allow efficient
!  storage of coupling coefficients.
!------------------------------------------------------------------------------
      module matrices

#ifndef USE_MULTI
      use mod_grid, only: grd, nr, nt, ndomains
#else
      use mod_grid, only: grd, nt, ndomains
#endif
      use inputs
      use derivative
#ifndef USE_1D
      use integrales
#endif

      use cfg

      implicit none

      type DOMAIN
        ! these quantities are the same for all processes (for the
        ! moment)
        integer :: nvar, d_dim, ku, kl, nvar_keep, offset
#ifdef USE_1D
        integer :: nas, nar, nasbc
        double precision, allocatable :: as(:), ar(:,:), asbc(:)
        integer, allocatable :: asi(:,:), ari(:,:), asbci(:,:)
#else
        integer :: nas, nart, nartt
        double precision, allocatable :: as(:), art(:,:,:),artt(:,:,:,:)
        integer, allocatable :: asi(:,:), arti(:,:), artti(:,:)
#endif
#ifdef USE_1D
        integer, allocatable :: ivar(:,:), ieq(:,:)
#else
        integer, allocatable :: ivar(:,:,:), ieq(:,:,:)
#endif
        integer, allocatable :: lvar(:,:), leq(:,:)
        integer, allocatable :: var_der_max(:), var_der_min(:)
        integer, allocatable :: var_list(:)
        logical, allocatable :: bc_flag(:), var_keep(:)
#ifdef USE_COMPLEX
        double complex, allocatable :: vect_der(:,:)
#else
        double precision, allocatable :: vect_der(:,:)
#endif
        character*(10), allocatable :: var_name(:), eq_name(:)
      end type DOMAIN

      type INTERDOMAIN
        integer :: natbc, nattbc
        double precision, allocatable :: atbc(:,:),attbc(:,:,:)
        integer, allocatable :: atbci(:,:), attbci(:,:)
        integer, allocatable :: v_bc_range(:), h_bc_range(:)
        integer, allocatable :: inv_v_bc_range(:), inv_h_bc_range(:)
        integer :: n_v_bc, n_h_bc
        integer :: v_bc_min, v_bc_max, h_bc_min, h_bc_max
      end type INTERDOMAIN

      logical, save :: first_run = .true.
      type(DERMAT),      allocatable, save :: dmat(:)
      type(DOMAIN),      allocatable, save :: dm(:)
      type(INTERDOMAIN), allocatable, save :: idm(:,:)

      integer, save :: a_dim, llmax, power_max
#ifdef USE_MULTI
      integer, save :: d_dim_max
#else
      logical, save :: first_dmat = .true.
#endif

!------------------------------------------------------------------------------
! List of variables:
!
! dm(:)        = these contain a sparse representation of the equations
!                within a same domain
! idm(:)       = these contain a sparse representation of the boundary
!                conditions and the interface conditions between different
!                domains
! dmat(:)      = these contain the derivation matrices for each domain
! power_max    = highest power of the eigenvalue
! a_dim        = total dimension of problem (in this case sum of d_dim)
! d_dim_max    = maximum value for d_dim
! first_run    = this keeps track of whether this is a first run or not so
!                as to know whether the different arrays in domain,
!                interdomain and dmat can be freed or not
!------------------------------------------------------------------------------
! List of variables in DOMAIN:
!
! as, art, artt, asi, arti, artti, nas, nart, nartt: contains sparse matrices
!                and relevent information.  See below for details.
! nvar         = number of variables
! nvar_keep    = number of variables in output file
! d_dim        = dimension of domain (i.e. nvar*nr*nt)
! offset       = offset of the domain with respect to the entire problem
! ku           = number of upper bands (useful for band storage)
! kl           = number of lower bands (useful for band storage)
! ivar(var,i,j)= index which corresponds to variable var at grid point i, j.
! ieq(eq,i,j)  = index which corresponds to equation eq at grid point i, j.
!                These determine the order of variables and equations, and
!                and can be used to obtain an efficient matrix storage.
!                They are set in init_order()
! lvar(j,var)  = these give the values of l for each variable var at
!                each "point" j.
! leq(j,var)   = these give the values of l for each equation eq at
!                each "point" j.
! var_list     = list of variables to print in output file (size nvar_keep)
! var_keep     = indicates which variables to print in output file (size nvar)
! var_der_max()= array which gives the highest derivative for each var
! var_der_min()= array which gives the lowest derivative for each var
! bc_flag      = This keeps track of which equations are replaced by boundary
!                conditions.  When a boundary condition is specified on a
!                certain line of an array, then that line cannot contain
!                any information from the equations, but must be used
!                exclusively on boundary condition(s).
! vect_der     = workspace which contains various derivatives of vect.
!                First index: ivar(var,i).  Second index: derivative order.
! var_name()   = names of the variables (it is used when writing vecp.gz)
! eq_name()    = names of the equations
!------------------------------------------------------------------------------
! List of variables in INTERDOMAIN:
!
! atbc, attbc, atbci, attbci, natbc, nattbc: contains boundary and interface
!                            conditions, and relevent information.  See below
!                            for details.
! n_v_bc                   = number of rows with b.c. or i.c.
! n_h_bc                   = number of columns with b.c. or i.c.
! v_bc_min                 = highest row with b.c. or i.c.
! v_bc_max                 = lowest row with b.c. or i.c.
! h_bc_min                 = left-most column with b.c. or i.c.
! h_bc_max                 = right-most column with b.c. or i.c.
! v_bc_range(index)        = position
! h_bc_range(index)        = position
! inv_v_bc_range(position) = index
! inv_h_bc_range(position) = index
!
! where
!
!   v = vertical, h = horizontal
!   1 <= index <= number of lines/columns with b.c. or i.c.
!   1 <= position <= number of lines/columns in matrix
!------------------------------------------------------------------------------
! The remaining variables are named through the following scheme:
! [{},n]a[{},bc][s,r][{},i] = coupling matrices and relevant information
!                       [{},n]:  {} = matrices themselves (or relevant
!                                     to specific matrices)
!                                n  = number of matrices
!                       [{},bc]: {} = operators, bc = boundary conditions
!                       [s,rt,rtt,t,tt]:
!                                s  = scalar
!                                rt = diagonal matrice in radial and angular directions
!                                rtt= diagonal matrice in radial direction,
!                                     full matrix in angular direction
!                                t  = diagonal matrix in angular direction, scalar in
!                                     radial direction
!                                tt = full matrix in angular direction, scalar in
!                                     radial direction
!                       [{},i]:  {} = matrix or scalar
!                                i  = information on matrix or scalar
! The last dimension to all these arrays is the index.
!
! A specific example:
! artt(1:nr,1:nt,1:nt,i): First coordinate: position in radial grid
!                         Second coordinate: angular coordinate
!                         Third coordinate: angular coordinate
!                         Fourth coordinate: the index
! artti(info,i): First coordinate (in [1,6]): tells what information is being given
!                  info = 1: power of the eigenvalue (e.g. w^2.a(2).x + w.a(1).x + a(0).x = 0)
!                  info = 2: derivative order
!                  info = 3: number of equation
!                  info = 4: number of variable
!                 (info = 5: equation location, only for boundary conditions (typically, 1 or nr))
!                 (info = 6: variable location, only for boundary conditions (typically, 1 or nr))
!                  info = 7: 0 for real, 1 for imaginary
!
!------------------------------------------------------------------------------

contains

      subroutine dump_asigma_matrix(asigma, filename, ierr)

          character(len=*), intent(in) :: filename
          character(256) :: fname
          integer, intent(out) :: ierr
#ifdef USE_COMPLEX
          double complex, intent(in) :: asigma(:, :)
#else
          double precision, intent(in) :: asigma(:, :)
#endif
          integer :: i, j, eq

          fname = filename // trim(tag)
          write(*, "(A, A)") "dump asigma to file: ", trim(fname)
          open(unit=42, file=trim(fname))

          do i = lbound(asigma, 1), ubound(asigma, 1)
              do j = lbound(asigma, 2), ubound(asigma, 2)
                  if (abs(asigma(i, j)) >= 1d-32) then
                      if (grd(1)%mattype == 'BAND') then
                          write(42, *)                              &
                              i+j-dm(1)%nvar+dm(1)%kl-dm(1)%ku-1,   &
                              j,                                    &
                              asigma(i, j)
                      else
                          write(42, *)      &
                              i,            &
                              j,            &
                              asigma(i, j)
                      endif
                  endif
              enddo
          enddo
          close(42)
          if (stop_after_dump) ierr = 2
      end subroutine dump_asigma_matrix

      subroutine dump_aterms(ierr)
          integer, intent(out) :: ierr

#ifdef USE_1D
          integer :: i, j
          character(256) :: fname

          fname = "asi" // trim(tag)
          write(*, "(A, A)") "dump asi to ", trim(fname)
          open(unit=42, file=trim(fname))
          do i = lbound(dm(1)%asi, 1), ubound(dm(1)%asi, 1)
              do j = lbound(dm(1)%asi, 2), ubound(dm(1)%asi, 2)
                  write(42, *) i, j, dm(1)%asi(i, j)
              enddo
          enddo
          close(42)

          fname = "as" // trim(tag)
          write(*, "(A, A)") "dump as to ", trim(fname)
          open(unit=42, file=trim(fname))
          do i = lbound(dm(1)%as, 1), ubound(dm(1)%as, 1)
              write(42, *) i, dm(1)%as(i)
          enddo
          close(42)

          fname = "ari" // trim(tag)
          write(*, "(A, A)") "dump ari to ", trim(fname)
          open(unit=42, file=trim(fname))
          do i = lbound(dm(1)%ari, 1), ubound(dm(1)%ari, 1)
              do j = lbound(dm(1)%ari, 2), ubound(dm(1)%ari, 2)
                  write(42, *) i, j, dm(1)%ari(i, j)
              enddo
          enddo
          close(42)

          fname = "ar" // trim(tag)
          write(*, "(A, A)") "dump ar to ", trim(fname)
          open(unit=42, file=trim(fname))
          do j = lbound(dm(1)%ar, 2), ubound(dm(1)%ar, 2)
              do i = lbound(dm(1)%ar, 1), ubound(dm(1)%ar, 1)
                  write(42, *) i, j, dm(1)%ar(i, j)
              enddo
          enddo
          close(42)

          fname = "asbci" // trim(tag)
          write(*, "(A, A)") "dump asbci to ", trim(fname)
          open(unit=42, file=trim(fname))
          do i = lbound(dm(1)%asbci, 1), ubound(dm(1)%asbci, 1)
              do j = lbound(dm(1)%asbci, 2), ubound(dm(1)%asbci, 2)
                  write(42, *) i, j, dm(1)%asbci(i, j)
              enddo
          enddo
          close(42)

          fname = "asbc" // trim(tag)
          write(*, "(A, A)") "dump asbc to ", trim(fname)
          open(unit=42, file=trim(fname))
          do i = lbound(dm(1)%asbc, 1), ubound(dm(1)%asbc, 1)
              write(42, *) i, dm(1)%asbc(i)
          enddo
          close(42)

          if (stop_after_dump) ierr = 2
#endif

          if (stop_after_dump) ierr = 2
      end subroutine dump_aterms


#ifdef USE_1D
      subroutine avg(vin, vout)

          double precision, intent(in) :: vin(nr)
          double precision, intent(out) :: vout(nr)
          integer i, ii

          do i=1, nr
              vout(i) = 0d0
              do ii=max(1,i-dmat(1)%lbder(0)),min(nr,i+dmat(1)%ubder(0))
                  vout(i) = vout(i) + dmat(1)%derive(i,ii,0)*vin(ii)
              enddo
          enddo
      end subroutine
#endif

!------------------------------------------------------------------------------
! This subroutine initialises the arrays dm, idm and dmat.  The arrays dm and
! idm, which contain the equations and boundary/interface conditions are
! initialised in "matrices.inc" which is written by readeq.
!------------------------------------------------------------------------------
          include "matrices.inc"

!------------------------------------------------------------------------------
! This subroutine finds the maximum l value in the lvar and leq arrays.
!------------------------------------------------------------------------------
#ifndef USE_1D
      subroutine find_llmax()

      integer id,j,var

      llmax = 0
      do id=1,ndomains
      do var = 1, dm(id)%nvar
          do j=1,nt
          if (dm(id)%lvar(j,var).gt.llmax) llmax = dm(id)%lvar(j,var)
          if (dm(id)%leq(j,var).gt.llmax)  llmax = dm(id)%leq(j,var)
          enddo
        enddo
      enddo

      end subroutine find_llmax
#endif

!------------------------------------------------------------------------------
! This subroutine initialises var_list, the list of variables to be
! printed in the output file, from the logical array var_keep.  It also
! sets the appropriate value for nvar_keep.
!------------------------------------------------------------------------------
! #ifndef USE_1D
      subroutine init_var_list()

      integer id, var, var_out

      do id=1,ndomains
          dm(id)%nvar_keep = 0
          do var = 1, dm(id)%nvar
              if (dm(id)%var_keep(var)) &
                & dm(id)%nvar_keep=dm(id)%nvar_keep+1
          enddo

          allocate(dm(id)%var_list(dm(id)%nvar_keep))

          var_out = 1
          do var = 1, dm(id)%nvar
              if (dm(id)%var_keep(var)) then
                  dm(id)%var_list(var_out) = var
                  var_out = var_out + 1
              endif
          enddo
      enddo


      end subroutine init_var_list
! #endif

!------------------------------------------------------------------------------
! This subroutine clears dm, idm and dmat so as to avoid memory leaks.
!------------------------------------------------------------------------------
      subroutine clear_all()
#ifdef USE_1D

          if (allocated(dmat)) deallocate(dmat)
          if (allocated(idm)) deallocate(idm)
          if (allocated(dm)) deallocate(dm)
#else

      use derivative
      integer id, id2

      ! clear derivation matrices
#ifdef USE_MULTI
      do id=1, ndomains
        call clear_derive(dmat(id))
      enddo
#endif
      if (allocated(dmat)) &
          deallocate(dmat)

      ! clear domain matrices
      do id=1, ndomains
          if (allocated(dm(id)%as)) &
              deallocate(dm(id)%as)
          if (allocated(dm(id)%art)) &
              deallocate(dm(id)%art)
          if (allocated(dm(id)%artt)) &
              deallocate(dm(id)%artt)

          if (allocated(dm(id)%asi)) &
              deallocate(dm(id)%asi)
          if (allocated(dm(id)%arti)) &
              deallocate(dm(id)%arti)
          if (allocated(dm(id)%artti)) &
              deallocate(dm(id)%artti)

          if (allocated(dm(id)%ivar)) &
              deallocate(dm(id)%ivar)
          if (allocated(dm(id)%ieq)) &
              deallocate(dm(id)%ieq)

          if (allocated(dm(id)%lvar)) &
              deallocate(dm(id)%lvar)
          if (allocated(dm(id)%leq)) &
              deallocate(dm(id)%leq)

          if (allocated(dm(id)%var_der_max)) &
              deallocate(dm(id)%var_der_max)
          if (allocated(dm(id)%var_der_min)) &
              deallocate(dm(id)%var_der_min)

          if (allocated(dm(id)%bc_flag)) &
              deallocate(dm(id)%bc_flag)
          if (allocated(dm(id)%vect_der)) &
              deallocate(dm(id)%vect_der)

          if (allocated(dm(id)%var_name)) &
              deallocate(dm(id)%var_name)
          if (allocated(dm(id)%eq_name)) &
              deallocate(dm(id)%eq_name)

          if (allocated(dm(id)%var_keep)) &
              deallocate(dm(id)%var_keep)
          if (allocated(dm(id)%var_list)) &
              deallocate(dm(id)%var_list)

      enddo

      if (allocated(dm)) &
          deallocate(dm)

      ! clear interdomain matrices
      do id=1, ndomains
          do id2=1, ndomains
              if (allocated(idm(id, id2)%atbc)) &
                  deallocate(idm(id, id2)%atbc)
              if (allocated(idm(id, id2)%attbc)) &
                  deallocate(idm(id, id2)%attbc)

              if (allocated(idm(id, id2)%atbci)) &
                  deallocate(idm(id, id2)%atbci)
              if (allocated(idm(id, id2)%attbci)) &
                  deallocate(idm(id, id2)%attbci)

              if (allocated(idm(id, id2)%v_bc_range)) &
                  deallocate(idm(id, id2)%v_bc_range)
              if (allocated(idm(id, id2)%h_bc_range)) &
                  deallocate(idm(id, id2)%h_bc_range)

              if (allocated(idm(id, id2)%inv_v_bc_range)) &
                  deallocate(idm(id, id2)%inv_v_bc_range)
              if (allocated(idm(id, id2)%inv_h_bc_range)) &
                  deallocate(idm(id, id2)%inv_h_bc_range)

          enddo
      enddo

      if (allocated(idm)) &
          deallocate(idm)

#endif
      end subroutine clear_all

!------------------------------------------------------------------------------
! This subroutine allocates most of the arrays in dm and idm
!------------------------------------------------------------------------------
      subroutine allocate_all()
#ifdef USE_1D

        allocate(dm(1)%as(dm(1)%nas))
        allocate(dm(1)%ar(1:nr, dm(1)%nar))
        allocate(dm(1)%asbc(dm(1)%nasbc))

        dm(1)%as = 0d0
        dm(1)%ar = 0d0
        dm(1)%asbc = 0d0
#else

      integer, parameter :: nai   = 7
      integer, parameter :: nabci = 7
      integer id, id2

      ! allocate and initialise domain matrices
      do id=1,ndomains
        allocate(dm(id)%as(dm(id)%nas))
        allocate(dm(id)%art(grd(id)%nr,nt,dm(id)%nart))
        allocate(dm(id)%artt(grd(id)%nr,nt,nt,dm(id)%nartt))
        allocate(dm(id)%asi(nai,dm(id)%nas))
        allocate(dm(id)%arti(nai,dm(id)%nart))
        allocate(dm(id)%artti(nai,dm(id)%nartt))
        allocate(dm(id)%lvar(nt,dm(id)%nvar))
        allocate(dm(id)%leq(nt,dm(id)%nvar))
        allocate(dm(id)%var_name(dm(id)%nvar))
        allocate(dm(id)%eq_name(dm(id)%nvar))
        allocate(dm(id)%var_keep(dm(id)%nvar))
        dm(id)%as   = 0d0
        dm(id)%art  = 0d0
        dm(id)%artt = 0d0
        dm(id)%lvar = 0
        dm(id)%leq  = 0
      enddo

      ! allocate and initialise interdomain matrices
      do id=1,ndomains
        do id2=1,ndomains
          allocate(idm(id,id2)%atbc(nt,idm(id,id2)%natbc))
          allocate(idm(id,id2)%attbc(nt,nt,idm(id,id2)%nattbc))
          allocate(idm(id,id2)%atbci(nabci,idm(id,id2)%natbc))
          allocate(idm(id,id2)%attbci(nabci,idm(id,id2)%nattbc))
          idm(id,id2)%atbc  = 0d0
          idm(id,id2)%attbc = 0d0
        enddo
      enddo
#endif

      end subroutine allocate_all

!------------------------------------------------------------------------------
! This subroutine allows you to do an ad-hoc modification to a coupling matrix
! so as to make the second velocity component, v, artificially appear in the
! equations for l = 0 (only applies when m = 0).
!------------------------------------------------------------------------------
! Variables:
!
! nr                = radial resolution
! mm                = azimuthal order (or some other appropriate quantity)
! f(1:nr,1:nt,1:nt) = coupling matrix which is being modified
!------------------------------------------------------------------------------
#ifdef USE_MULTI
      subroutine modify_l0(f,nr,mm)

      integer, intent(in)             :: nr, mm
      double precision, intent(inout) :: f(1:nr,1:nt,1:nt)
      integer i

      if (mm.eq.0) then
        do i=1,nr
          f(i,1,1) = 1d0
        enddo
      endif

      end subroutine modify_l0
#else
      subroutine modify_l0(f,m)

      double precision f(1:grd(1)%nr,1:nt,1:nt)
      integer m, i

      if (m.eq.0) then
        do i=1,grd(1)%nr
          f(i,1,1) = 1d0
        enddo
      endif

      end subroutine modify_l0
#endif

!------------------------------------------------------------------------------
! This subroutine interpolates background equilibrium quantities so that
! improved (alternate) finite differences approach can be applied.
!
! IMPORTANT: this only applies to the 1D case.  For the 2D case, use
!            subroutines such as Illm_a, Jllm_a, etc.
!------------------------------------------------------------------------------
! Variables:
!
! vin(1:grd(id)%nr)  = variable before interpolation
! vout(1:grd(id)%nr) = variable after interpolation
!------------------------------------------------------------------------------
#ifdef USE_MULTI
      subroutine avg1D(vin,vout,id)

      integer, intent(in) :: id
      double precision, intent(in) :: vin(grd(id)%nr)
      double precision, intent(out) :: vout(grd(id)%nr)
      integer i, ii

      do i=1,grd(id)%nr
        vout(i) = 0d0
        do ii=max(1,i-dmat(id)%lbder(0)),min(grd(id)%nr,i+dmat(id)%ubder(0))
          vout(i) = vout(i) + dmat(id)%derive(i,ii,0)*vin(ii)
        enddo
      enddo
      end subroutine
#endif

!------------------------------------------------------------------------------
! This subroutine assigns a 1D array to a 2D array, by repeating the
! values on the second coordinate.  This is useful for doing 1D
! simulations in this 2D context.
!------------------------------------------------------------------------------
! Variables:
!
! v1D(1:grd(id)%nr)    = 1D variable before assignement
! v2D(1:grd(id)%nr,nt) = 2D variable after assignement
!------------------------------------------------------------------------------
      subroutine assign1D(v1D,v2D,id)

      integer, intent(in) :: id
      double precision, intent(in) :: v1D(grd(id)%nr)
      double precision, intent(out) :: v2D(grd(id)%nr,nt)
      integer i, j

      do i=1,grd(id)%nr
        do j=1,nt
          v2D(i,j) = v1D(i)
        enddo
      enddo
      end subroutine

!------------------------------------------------------------------------------
! This subroutine finds the min and max derivative domain id and for
! each variable in the domain.  The global min and max derivatives
! are set in der_min and der_max, and the min and max derivatives for
! each variable are set in dm(id)%var_der_min and dm(id)%var_der_max.)
!------------------------------------------------------------------------------
! Variables:
!
! der_min = lowest derivative for domain id
! der_max = highest derivative for domain id
! id      = number of the domain
!------------------------------------------------------------------------------
#ifdef USE_MULTI
      subroutine find_der_range(der_min,der_max,id)
#else
      subroutine find_der_range(der_min,der_max)
#endif

#ifdef USE_MULTI
      integer, intent(in)  :: id
#else
      integer, parameter :: id = 1
#endif
      integer, intent(out) :: der_min, der_max
      integer n, der, var, i, id2

#ifdef USE_1D
      integer :: eq

      if (allocated(dm(1)%var_der_max)) deallocate(dm(1)%var_der_max)
      if (allocated(dm(1)%var_der_min)) deallocate(dm(1)%var_der_min)
      allocate(dm(1)%var_der_max(dm(1)%nvar), dm(1)%var_der_min(dm(1)%nvar))

      der_min = 0
      der_max = 0
      dm(1)%var_der_max(:) = 0
      dm(1)%var_der_min(:) = 0

      do n=1, dm(1)%nas
        der = dm(1)%asi(2,n)
        eq  = dm(1)%asi(3,n)
        var = dm(1)%asi(4,n)
        if (der.gt.dm(1)%var_der_max(var)) dm(1)%var_der_max(var) = der
        if (der.lt.dm(1)%var_der_min(var)) dm(1)%var_der_min(var) = der
        if (der.gt.der_max)          der_max = der
        if (der.lt.der_min)          der_min = der
      enddo

      do n=1, dm(1)%nar
        der = dm(1)%ari(2,n)
        eq  = dm(1)%ari(3,n)
        var = dm(1)%ari(4,n)
        if (der.gt.dm(1)%var_der_max(var)) dm(1)%var_der_max(var) = der
        if (der.lt.dm(1)%var_der_min(var)) dm(1)%var_der_min(var) = der
        if (der.gt.der_max)          der_max = der
        if (der.lt.der_min)          der_min = der
      enddo

      do n=1, dm(1)%nasbc
        der = dm(1)%asbci(2,n)
        eq  = dm(1)%asbci(3,n)
        var = dm(1)%asbci(4,n)
        if (der.gt.dm(1)%var_der_max(var)) dm(1)%var_der_max(var) = der
        if (der.lt.dm(1)%var_der_min(var)) dm(1)%var_der_min(var) = der
        if (der.gt.der_max)          der_max = der
        if (der.lt.der_min)          der_min = der
      enddo

#else

      allocate(dm(id)%var_der_max(dm(id)%nvar))
      allocate(dm(id)%var_der_min(dm(id)%nvar))

      der_min = 0
      der_max = 0
      dm(id)%var_der_max(:) = 0
      dm(id)%var_der_min(:) = 0

      do n=1, dm(id)%nas
        der = dm(id)%asi(2,n)
        var = dm(id)%asi(4,n)
        if (der.gt.dm(id)%var_der_max(var)) dm(id)%var_der_max(var) = der
        if (der.lt.dm(id)%var_der_min(var)) dm(id)%var_der_min(var) = der
        if (der.gt.der_max)          der_max = der
        if (der.lt.der_min)          der_min = der
      enddo

      do n=1,dm(id)%nart
        der = dm(id)%arti(2,n)
        var = dm(id)%arti(4,n)
        if (der.gt.dm(id)%var_der_max(var)) dm(id)%var_der_max(var) = der
        if (der.lt.dm(id)%var_der_min(var)) dm(id)%var_der_min(var) = der
        if (der.gt.der_max)          der_max = der
        if (der.lt.der_min)          der_min = der
      enddo

      do n=1,dm(id)%nartt
        der = dm(id)%artti(2,n)
        var = dm(id)%artti(4,n)
        if (der.gt.dm(id)%var_der_max(var)) dm(id)%var_der_max(var) = der
        if (der.lt.dm(id)%var_der_min(var)) dm(id)%var_der_min(var) = der
        if (der.gt.der_max)          der_max = der
        if (der.lt.der_min)          der_min = der
      enddo

      ! Boundary conditions:
      ! Reminder: variables go according to columns,
      !           equations go according to rows ...
      do id2=1,ndomains
        do n=1,idm(id2,id)%natbc
          der = idm(id2,id)%atbci(2,n)
          var = idm(id2,id)%atbci(4,n)
          if (der.gt.dm(id)%var_der_max(var)) dm(id)%var_der_max(var) = der
          if (der.lt.dm(id)%var_der_min(var)) dm(id)%var_der_min(var) = der
          if (der.gt.der_max)          der_max = der
          if (der.lt.der_min)          der_min = der
        enddo
      enddo

      do id2=1,ndomains
        do n=1,idm(id2,id)%nattbc
          der = idm(id2,id)%attbci(2,n)
          var = idm(id2,id)%attbci(4,n)
          if (der.gt.dm(id)%var_der_max(var)) dm(id)%var_der_max(var) = der
          if (der.lt.dm(id)%var_der_min(var)) dm(id)%var_der_min(var) = der
          if (der.gt.der_max)          der_max = der
          if (der.lt.der_min)          der_min = der
        enddo
      enddo
#endif

      end subroutine find_der_range

!------------------------------------------------------------------------------
! This subroutine initialises the arrays ivar and ieq which govern the order
! in which the variables and the equations appear in the matrices.
! The following is a generic form:
!------------------------------------------------------------------------------
      include "order.inc"
!------------------------------------------------------------------------------
! This subroutine checks the arrays ivar and ieq to make sure they have been
! correctly initialised.
!------------------------------------------------------------------------------

#ifndef USE_1D
      subroutine check_order(ierr)
      integer i,j,id,var,eq,max_value, min_value
      logical, allocatable :: eq_flag(:), var_flag(:)
      integer, intent(out) :: ierr

      do id=1,ndomains
        min_value = 1
        max_value = dm(id)%nvar*grd(id)%nr*nt
        do var=1,dm(id)%nvar
          do i=1,grd(id)%nr
            do j=1,nt
              if ((dm(id)%ieq(var,i,j).gt.max_value).or. &
                  (dm(id)%ieq(var,i,j).lt.min_value)) then
                print*, "ieq takes on values outside allowed range: ", &
                         1, max_value
                print*, "faulty value = ",dm(id)%ieq(var,i,j)
                print*, "should be: ", j + nt*(var-1 + (dm(1)%nvar-2)*(i-1))
                print*, "nvar:", dm(id)%nvar
                print*, "nr:", grd(id)%nr
                print*, "nt:", nt
                print*, "var = ",var
                print*, "i = ",i
                print*, "j = ",j
                print*, "id = ",id
                ierr = 1
                return
                ! stop ! the following tests won't function correctly
              endif
            enddo
          enddo
        enddo

        do var=1,dm(id)%nvar
          do i=1,grd(id)%nr
            do j=1,nt
              if ((dm(id)%ivar(var,i,j).gt.max_value).or. &
                  (dm(id)%ivar(var,i,j).lt.min_value)) then
                print*,"ivar takes on values outside allowed range"
                print*,"faulty value = ",dm(id)%ivar(var,i,j)
                print*,"var = ",var
                print*,"i = ",i
                print*,"j = ",j
                ierr = 1
                return
                ! stop ! the following tests won't function correctly
              endif
            enddo
          enddo
        enddo

        allocate(eq_flag(max_value),var_flag(max_value))

        do i=1,dm(id)%nvar*grd(id)%nr*nt
          eq_flag(i) = .false.
          var_flag(i) = .false.
        enddo

        do var=1,dm(id)%nvar
          do i=1,grd(id)%nr
            do j=1,nt
              eq_flag(dm(id)%ieq(var,i,j)) = .true.
              var_flag(dm(id)%ivar(var,i,j)) = .true.
            enddo
          enddo
        enddo

        do i=1,dm(id)%nvar*grd(id)%nr*nt
          if (.not.eq_flag(i)) then
            print*,"ieq not properly initialised"
            exit
          endif
        enddo

        do i=1,dm(id)%nvar*grd(id)%nr*nt
          if (.not.var_flag(i)) then
            print*,"ivar not properly initialised"
            exit
          endif
        enddo

        deallocate(var_flag,eq_flag)
      enddo

      print*,"ieq and ivar seem to be initialised correctly."
      print*,"Please remove 'check_order' and run program again."
      ierr = 1
      return
      ! stop

      end subroutine check_order
#endif

!------------------------------------------------------------------------------
! This subroutine initialises the array bc_flag which specifies which lines
! are taken up by boundary conditions.
!
! This must come after init_order.
!------------------------------------------------------------------------------

      subroutine init_bc_flag() bind(c)
#ifdef USE_1D
      integer n,eq,loc,l

      if (allocated(dm(1)%bc_flag)) deallocate(dm(1)%bc_flag)
      allocate(dm(1)%bc_flag(a_dim))

      ! The boundary conditions:
      dm(1)%bc_flag(:) = .false.

      do n=1, dm(1)%nasbc
        eq  = dm(1)%asbci(3,n)
        loc = dm(1)%asbci(5,n)
        l = dm(1)%ieq(eq,loc)
        dm(1)%bc_flag(l) = .true.
      enddo
#else
      integer id,id2,n,j,eq,loc

      do id=1,ndomains
        allocate(dm(id)%bc_flag(dm(id)%d_dim))

        ! The boundary conditions:
        dm(id)%bc_flag(:) = .false.

        do id2=1,ndomains
          do n=1,idm(id,id2)%natbc
            eq  = idm(id,id2)%atbci(3,n)
            loc = idm(id,id2)%atbci(5,n)
            do j=1,nt
              dm(id)%bc_flag(dm(id)%ieq(eq,loc,j)) = .true.
            enddo
          enddo

          do n=1,idm(id,id2)%nattbc
            eq  = idm(id,id2)%attbci(3,n)
            loc = idm(id,id2)%attbci(5,n)
            do j=1,nt
              dm(id)%bc_flag(dm(id)%ieq(eq,loc,j)) = .true.
            enddo
          enddo
        enddo
      enddo

#ifdef USE_MULTI
      call init_bc_range()
#endif
#endif
      end subroutine init_bc_flag

!------------------------------------------------------------------------------
! This subroutine initialises the arrays v_bc_range, h_bc_range,
! inv_v_bc_range and inv_h_bc_range which are used when calculating
! corrections to the matrices which operate within one domain (i.e.
! id = id2).  Their lower and upper bounds intervene in determining
! the number of lower and upper diagonal bands in these matrices,
! when working with band storage.
!
! This must come after init_order.
!------------------------------------------------------------------------------
! REMINDER:
!
! n_v_bc                   = number of rows with b.c. or i.c.
! n_h_bc                   = number of columns with b.c. or i.c.
! v_bc_min                 = highest row with b.c. or i.c.
! v_bc_max                 = lowest row with b.c. or i.c.
! h_bc_min                 = left-most column with b.c. or i.c.
! h_bc_max                 = right-most column with b.c. or i.c.
! v_bc_range(index)        = position
! h_bc_range(index)        = position
! inv_v_bc_range(position) = index
! inv_h_bc_range(position) = index
!
! where
!
!   v = vertical, h = horizontal
!   1 <= index <= number of lines/columns with b.c. or i.c.
!   1 <= position <= number of lines/columns in matrix
!------------------------------------------------------------------------------
#ifdef USE_MULTI
      subroutine init_bc_range()

      integer id,id2,n,i,j,der,eq,var,eqloc,varloc,l,c
      integer n_v_bc, n_h_bc
      logical, allocatable :: v_bc_flag(:), h_bc_flag(:)

      allocate(v_bc_flag(d_dim_max), h_bc_flag(d_dim_max))

      do id=1,ndomains
        do id2=1,ndomains
          v_bc_flag(1:dm(id)%d_dim)  = .false.
          h_bc_flag(1:dm(id2)%d_dim) = .false.

          ! The boundary conditions:
          do n=1,idm(id,id2)%natbc
            der = idm(id,id2)%atbci(2,n)
            eq  = idm(id,id2)%atbci(3,n)
            var = idm(id,id2)%atbci(4,n)
            eqloc  = idm(id,id2)%atbci(5,n)
            varloc = idm(id,id2)%atbci(6,n)
            do j=1,nt
              l = dm(id)%ieq(eq,eqloc,j)
              v_bc_flag(l) = .true.
              do i=max(1,varloc-dmat(id2)%lbder(der)), &
                   min(grd(id2)%nr,varloc+dmat(id2)%ubder(der))
                c = dm(id2)%ivar(var,i,j)
                h_bc_flag(c) = .true.
              enddo
            enddo
          enddo

          do n=1,idm(id,id2)%nattbc
            der = idm(id,id2)%attbci(2,n)
            eq  = idm(id,id2)%attbci(3,n)
            var = idm(id,id2)%attbci(4,n)
            eqloc  = idm(id,id2)%attbci(5,n)
            varloc = idm(id,id2)%attbci(6,n)
            do j=1,nt
              l = dm(id)%ieq(eq,eqloc,j)
              v_bc_flag(l) = .true.
              do i=max(1,varloc-dmat(id2)%lbder(der)), &
                   min(grd(id2)%nr,varloc+dmat(id2)%ubder(der))
                ! it is not necessary to do another loop with jj=1,nt
                c = dm(id2)%ivar(var,i,j)
                h_bc_flag(c) = .true.
              enddo
            enddo
          enddo

          ! the corresponding ranges
          n_v_bc = 0
          do i=1,dm(id)%d_dim
            if (v_bc_flag(i)) n_v_bc = n_v_bc + 1
          enddo

          idm(id,id2)%n_v_bc = n_v_bc
          allocate(idm(id,id2)%v_bc_range(n_v_bc))
          allocate(idm(id,id2)%inv_v_bc_range(dm(id)%d_dim))
          n_v_bc = 0
          do i=1,dm(id)%d_dim
            if (v_bc_flag(i)) then
              n_v_bc = n_v_bc + 1
              idm(id,id2)%v_bc_range(n_v_bc) = i
              idm(id,id2)%inv_v_bc_range(i) = n_v_bc
            else
              idm(id,id2)%inv_v_bc_range(i) = 0
            endif
          enddo

          n_h_bc = 0
          do i=1,dm(id2)%d_dim
            if (h_bc_flag(i)) n_h_bc = n_h_bc + 1
          enddo

          idm(id,id2)%n_h_bc = n_h_bc
          allocate(idm(id,id2)%h_bc_range(n_h_bc))
          allocate(idm(id,id2)%inv_h_bc_range(dm(id2)%d_dim))
          n_h_bc = 0
          do i=1,dm(id2)%d_dim
            if (h_bc_flag(i)) then
              n_h_bc = n_h_bc + 1
              idm(id,id2)%h_bc_range(n_h_bc) = i
              idm(id,id2)%inv_h_bc_range(i) = n_h_bc
            else
              idm(id,id2)%inv_h_bc_range(i) = 0
            endif
          enddo

          ! This condition is necessary to avoid accessing elements
          ! from an empty array.
          if (n_v_bc.gt.0) then
            idm(id,id2)%v_bc_min = idm(id,id2)%v_bc_range(1)
            idm(id,id2)%v_bc_max = idm(id,id2)%v_bc_range(n_v_bc)
          else
            ! these conditions ensure that ku and kl are not modified
            ! in increase_ku_kl_upward/downward
            idm(id,id2)%v_bc_min = 0
            idm(id,id2)%v_bc_max = 0
          endif

          ! This condition is necessary to avoid accessing elements
          ! from an empty array.
          if (n_h_bc.gt.0) then
            idm(id,id2)%h_bc_min = idm(id,id2)%h_bc_range(1)
            idm(id,id2)%h_bc_max = idm(id,id2)%h_bc_range(n_h_bc)
          else
            ! these conditions ensure that ku and kl are not modified
            ! in increase_ku_kl_upward/downward
            idm(id,id2)%h_bc_min = 0
            idm(id,id2)%h_bc_max = 0
          endif
        enddo
      enddo

      deallocate(v_bc_flag, h_bc_flag)
      end subroutine init_bc_range
#endif

!------------------------------------------------------------------------------
! For a given domain, this subroutine sets appropriate values for dm(id)%ku
! and dm(id)%kl, the number of upper and lower bands, when using band storage
! for the matrix asigma(id,id):
!      asigma = a(0) + sigma.a(1) + sigma^2.a(2) +  ...
!------------------------------------------------------------------------------
! Variables:
!
! id = number of the domain
!------------------------------------------------------------------------------
#ifdef USE_MULTI
      subroutine init_ku_kl(id)
#else
      subroutine init_ku_kl()
#endif

#ifdef USE_MULTI
      integer, intent(in) :: id
#else
      integer, parameter :: id = 1
#endif
      integer n, i, j, ii, jj, l, c, der, eq, var, eqloc, varloc

      dm(id)%ku = 0
      dm(id)%kl = 0

#ifdef USE_1D

      ! The equations:
      do n=1, dm(1)%nas
        der = dm(1)%asi(2,n)
        eq  = dm(1)%asi(3,n)
        var = dm(1)%asi(4,n)
        do i=1,nr
          l = dm(1)%ieq(eq,i)
          if (.not.dm(1)%bc_flag(l)) then
            do ii=max(1,i-dmat(1)%lbder(der)),min(nr,i+dmat(1)%ubder(der))
              c = dm(1)%ivar(var,ii)
              if (dm(1)%kl.lt.(l-c)) dm(1)%kl = l-c
              if (dm(1)%ku.lt.(c-l)) dm(1)%ku = c-l
            enddo
          endif
        enddo
      enddo

      do n=1, dm(1)%nar
        der = dm(1)%ari(2,n)
        eq  = dm(1)%ari(3,n)
        var = dm(1)%ari(4,n)
        do i=1,nr
          l = dm(1)%ieq(eq,i)
          if (.not.dm(1)%bc_flag(l)) then
            do ii=max(1,i-dmat(1)%lbder(der)),min(nr,i+dmat(1)%ubder(der))
              c = dm(1)%ivar(var,ii)
              if (dm(1)%kl.lt.(l-c)) dm(1)%kl = l-c
              if (dm(1)%ku.lt.(c-l)) dm(1)%ku = c-l
            enddo
          endif
        enddo
      enddo

      ! Boundary conditions:
      do n=1, dm(1)%nasbc
        der = dm(1)%asbci(2,n)
        eq  = dm(1)%asbci(3,n)
        var = dm(1)%asbci(4,n)
        eqloc  = dm(1)%asbci(5,n)
        varloc = dm(1)%asbci(6,n)
        l = dm(1)%ieq(eq,eqloc)
        do ii=max(1,varloc-dmat(1)%lbder(der)),min(nr,varloc+dmat(1)%ubder(der))
          c = dm(1)%ivar(var,ii)
          if (dm(1)%kl.lt.(l-c)) dm(1)%kl = l-c
          if (dm(1)%ku.lt.(c-l)) dm(1)%ku = c-l
        enddo
      enddo

      print*, "Number of lower bands: ", dm(1)%kl
      print*, "Number of upper bands: ", dm(1)%ku

#else
      ! The equations:
      do n=1,dm(id)%nas
        der = dm(id)%asi(2,n)
        eq  = dm(id)%asi(3,n)
        var = dm(id)%asi(4,n)
        do i=1,grd(id)%nr
          do j=1,nt
            l = dm(id)%ieq(eq,i,j)
            if (.not.dm(id)%bc_flag(l)) then
              do ii=max(1,i-dmat(id)%lbder(der)), &
                  & min(grd(id)%nr,i+dmat(id)%ubder(der))
                c = dm(id)%ivar(var,ii,j)
                if (dm(id)%kl.lt.(l-c)) dm(id)%kl = l-c
                if (dm(id)%ku.lt.(c-l)) dm(id)%ku = c-l
              enddo
            endif
          enddo
        enddo
      enddo

      do n=1,dm(id)%nart
        der = dm(id)%arti(2,n)
        eq  = dm(id)%arti(3,n)
        var = dm(id)%arti(4,n)
        do i=1,grd(id)%nr
          do j=1,nt
            l = dm(id)%ieq(eq,i,j)
            if (.not.dm(id)%bc_flag(l)) then
              do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                c = dm(id)%ivar(var,ii,j)
                if (dm(id)%kl.lt.(l-c)) dm(id)%kl = l-c
                if (dm(id)%ku.lt.(c-l)) dm(id)%ku = c-l
              enddo
            endif
          enddo
        enddo
      enddo

      do n=1,dm(id)%nartt
        der = dm(id)%artti(2,n)
        eq  = dm(id)%artti(3,n)
        var = dm(id)%artti(4,n)
        do i=1,grd(id)%nr
          do j=1,nt
            l = dm(id)%ieq(eq,i,j)
            if (.not.dm(id)%bc_flag(l)) then
              do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                do jj=1,nt
                  c = dm(id)%ivar(var,ii,jj)
                  if (dm(id)%kl.lt.(l-c)) dm(id)%kl = l-c
                  if (dm(id)%ku.lt.(c-l)) dm(id)%ku = c-l
                enddo
              enddo
            endif
          enddo
        enddo
      enddo

      ! Boundary conditions:
      do n=1,idm(id,id)%natbc
        der = idm(id,id)%atbci(2,n)
        eq  = idm(id,id)%atbci(3,n)
        var = idm(id,id)%atbci(4,n)
        eqloc  = idm(id,id)%atbci(5,n)
        varloc = idm(id,id)%atbci(6,n)
        do j=1,nt
          l = dm(id)%ieq(eq,eqloc,j)
          do ii=max(1,varloc-dmat(id)%lbder(der)),min(grd(id)%nr,varloc+dmat(id)%ubder(der))
            c = dm(id)%ivar(var,ii,j)
            if (dm(id)%kl.lt.(l-c)) dm(id)%kl = l-c
            if (dm(id)%ku.lt.(c-l)) dm(id)%ku = c-l
          enddo
        enddo
      enddo

      do n=1,idm(id,id)%nattbc
        der = idm(id,id)%attbci(2,n)
        eq  = idm(id,id)%attbci(3,n)
        var = idm(id,id)%attbci(4,n)
        eqloc  = idm(id,id)%attbci(5,n)
        varloc = idm(id,id)%attbci(6,n)
        do j=1,nt
          l = dm(id)%ieq(eq,eqloc,j)
          do ii=max(1,varloc-dmat(id)%lbder(der)),min(grd(id)%nr,varloc+dmat(id)%ubder(der))
            do jj=1,nt
              c = dm(id)%ivar(var,ii,jj)
              if (dm(id)%kl.lt.(l-c)) dm(id)%kl = l-c
              if (dm(id)%ku.lt.(c-l)) dm(id)%ku = c-l
            enddo
          enddo
        enddo
      enddo

      print*,"**************************************"
#ifdef USE_MULTI
      print*,"DOMAIN: ",id
#endif
      print*,"Number of lower bands (before): ",dm(id)%kl
      print*,"Number of upper bands (before): ",dm(id)%ku
#endif

      end subroutine init_ku_kl

!------------------------------------------------------------------------------
! This adjusts the values of dm(id)%ku and dm(id)%kl so that the system can
! be solved using a downward sweep followed by an upward sweep.
!
! This MUST come after init_ku_kl and init_bc_range.
!------------------------------------------------------------------------------
! Variables:
!
! id = number of the domain
!------------------------------------------------------------------------------
#ifdef USE_MULTI
      subroutine increase_ku_kl_downward (id)

      integer, intent(in) :: id
      integer l, c

      if (id.eq.1) return

      l = idm(id,id-1)%v_bc_min
      c = idm(id-1,id)%h_bc_max
      if (dm(id)%ku.lt.(c-l)) dm(id)%ku = c-l
      l = idm(id,id-1)%v_bc_max
      c = idm(id-1,id)%h_bc_min
      if (dm(id)%kl.lt.(l-c)) dm(id)%kl = l-c

      print*,"Number of lower bands (after):  ",dm(id)%kl
      print*,"Number of upper bands (after):  ",dm(id)%ku

      end subroutine increase_ku_kl_downward
#endif

!------------------------------------------------------------------------------
! This adjusts the values of dm(id)%ku and dm(id)%kl so that the system can
! be solved using a upward sweep followed by an downward sweep.
!
! This must come after init_ku_kl
!------------------------------------------------------------------------------
! Variables:
!
! id = number of the domain
!------------------------------------------------------------------------------
#ifndef USE_1D
      subroutine increase_ku_kl_upward (id)

      integer, intent(in) :: id
      integer l, c

      if (id.eq.ndomains) return

      l = idm(id,id+1)%v_bc_min
      c = idm(id+1,id)%h_bc_max
      if (dm(id)%ku.lt.(c-l)) dm(id)%ku = c-l
      l = idm(id,id+1)%v_bc_max
      c = idm(id+1,id)%h_bc_min
      if (dm(id)%kl.lt.(l-c)) dm(id)%kl = l-c

      print*, "Number of lower bands (after):  ", dm(id)%kl
      print*, "Number of upper bands (after):  ", dm(id)%ku

      end subroutine increase_ku_kl_upward
#endif

!------------------------------------------------------------------------------
! This subroutine performs the matrix vector product a(power).x for the whole
! system, where a(power) commes from the eigenvalue problem:
!      w^2.a(2).x + w.a(1).x + a(0).x = 0.
!------------------------------------------------------------------------------
! List of variables:
!
! vect_in(1:a_dim)  = vecteur undergoing the matrix product
! vect_out(1:a_dim) = vecteur receiving the result from the matrix-vector
!                     product
! power             = this tells which matrix to use (it corresponds to the
!                     power of the eigenvalue in front)
!------------------------------------------------------------------------------

#ifdef USE_MULTI
      subroutine a_product_total(vect_in, vect_out, power)

#ifdef USE_COMPLEX
      double complex, intent(in) :: vect_in(a_dim)
      double complex, intent(out):: vect_out(a_dim)
      double complex :: zero = (0d0, 0d0)
      double complex cfactor
#else
      double precision, intent(in) :: vect_in(a_dim)
      double precision, intent(out):: vect_out(a_dim)
      double precision :: zero = 0d0
      double precision :: cfactor = 1d0
#endif
      integer, intent(in) :: power
      integer id, id2, i, j, ii, jj, n, der, eq, var, a_index, v_index
      integer eqloc, varloc

      vect_out(1:a_dim) = zero

      ! Preliminary calculations: derivative(s) of vect
      do id = 1, ndomains
        dm(id)%vect_der(1:dm(id)%d_dim,dmat(id)%der_min:dmat(id)%der_max) = zero
        do var=1,dm(id)%nvar
          do der=dm(id)%var_der_min(var),dm(id)%var_der_max(var)
            do i=1,grd(id)%nr
              do j=1,nt
                v_index = dm(id)%ivar(var,i,j)
                do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                  a_index = dm(id)%offset+dm(id)%ivar(var,ii,j)
                  dm(id)%vect_der(v_index,der) = dm(id)%vect_der(v_index,der) + &
                                 dmat(id)%derive(i,ii,der)*vect_in(a_index)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      ! The equations:
      do id = 1, ndomains
        do n=1,dm(id)%nas
          if (dm(id)%asi(1,n).eq.power) then
            der = dm(id)%asi(2,n)
            eq  = dm(id)%asi(3,n)
            var = dm(id)%asi(4,n)
#ifdef USE_COMPLEX
            if (dm(id)%asi(7,n).eq.0) then
              cfactor = (1d0,0d0)
            else
              cfactor = (0d0,1d0)
            endif
#endif
            do i=1,grd(id)%nr
              do j=1,nt
                v_index = dm(id)%ieq(eq,i,j)
                if (.not.dm(id)%bc_flag(v_index)) then
                  v_index = v_index + dm(id)%offset
                  a_index = dm(id)%ivar(var,i,j)
                  vect_out(v_index) = vect_out(v_index) +    &
                    cfactor*dm(id)%as(n)*dm(id)%vect_der(a_index,der)
                endif
              enddo
            enddo
          endif
        enddo

        do n=1,dm(id)%nart
          if (dm(id)%arti(1,n).eq.power) then
            der = dm(id)%arti(2,n)
            eq  = dm(id)%arti(3,n)
            var = dm(id)%arti(4,n)
#ifdef USE_COMPLEX
            if (dm(id)%arti(7,n).eq.0) then
              cfactor = (1d0,0d0)
            else
              cfactor = (0d0,1d0)
            endif
#endif
            do i=1,grd(id)%nr
              do j=1,nt
                v_index = dm(id)%ieq(eq,i,j)
                if (.not.dm(id)%bc_flag(v_index)) then
                  v_index = v_index + dm(id)%offset
                  a_index = dm(id)%ivar(var,i,j)
                  vect_out(v_index) = vect_out(v_index) + &
                    cfactor*dm(id)%art(i,j,n)*dm(id)%vect_der(a_index,der)
                endif
              enddo
            enddo
          endif
        enddo

        do n=1,dm(id)%nartt
          if (dm(id)%artti(1,n).eq.power) then
            der = dm(id)%artti(2,n)
            eq  = dm(id)%artti(3,n)
            var = dm(id)%artti(4,n)
#ifdef USE_COMPLEX
            if (dm(id)%artti(7,n).eq.0) then
              cfactor = (1d0,0d0)
            else
              cfactor = (0d0,1d0)
            endif
#endif
            do i=1,grd(id)%nr
              do j=1,nt
                v_index = dm(id)%ieq(eq,i,j)
                if (.not.dm(id)%bc_flag(v_index)) then
                  v_index = v_index + dm(id)%offset
                  do jj=1,nt
                    a_index = dm(id)%ivar(var,i,jj)
                    vect_out(v_index) = vect_out(v_index) + &
                      cfactor*dm(id)%artt(i,j,jj,n)*dm(id)%vect_der(a_index,der)
                  enddo
                endif
              enddo
            enddo
          endif
        enddo

        ! Boundary conditions:
        do id2=1,ndomains
          do n=1,idm(id,id2)%natbc
            if (idm(id,id2)%atbci(1,n).eq.power) then
              der = idm(id,id2)%atbci(2,n)
              eq  = idm(id,id2)%atbci(3,n)
              var = idm(id,id2)%atbci(4,n)
              eqloc  = idm(id,id2)%atbci(5,n)
              varloc = idm(id,id2)%atbci(6,n)
#ifdef USE_COMPLEX
              if (idm(id,id2)%atbci(7,n).eq.0) then
                cfactor = (1d0,0d0)
              else
                cfactor = (0d0,1d0)
              endif
#endif
              do j=1,nt
                v_index = dm(id)%offset+dm(id)%ieq(eq,eqloc,j)
                a_index = dm(id2)%ivar(var,varloc,j)
                vect_out(v_index) = vect_out(v_index) + &
                  cfactor*idm(id,id2)%atbc(j,n)*dm(id2)%vect_der(a_index,der)
              enddo
            endif
          enddo

          do n=1,idm(id,id2)%nattbc
            if (idm(id,id2)%attbci(1,n).eq.power) then
              der = idm(id,id2)%attbci(2,n)
              eq  = idm(id,id2)%attbci(3,n)
              var = idm(id,id2)%attbci(4,n)
              eqloc  = idm(id,id2)%attbci(5,n)
              varloc = idm(id,id2)%attbci(6,n)
#ifdef USE_COMPLEX
              if (idm(id,id2)%attbci(7,n).eq.0) then
                cfactor = (1d0,0d0)
              else
                cfactor = (0d0,1d0)
              endif
#endif
              do j=1,nt
                v_index = dm(id)%offset+dm(id)%ieq(eq,eqloc,j)
                do jj=1,nt
                  a_index = dm(id2)%ivar(var,varloc,jj)
                  vect_out(v_index) = vect_out(v_index) + &
                    cfactor*idm(id,id2)%attbc(j,jj,n)*dm(id2)%vect_der(a_index,der)
                enddo
              enddo
            endif
          enddo
        enddo
      enddo

      end subroutine a_product_total
#endif

!------------------------------------------------------------------------------
! This does:
!  vect_out(id) <-  vect_out(id) - asigma(id,id2)*vect_in(id2)
! using a sparse multiplication method.
!------------------------------------------------------------------------------
! Variables:
!
! vect_in  = input vector (dim = dm(id2)%d_dim)
! vect_out = output vector (dim = dm(id)%d_dim)
! id, id2  = domain coordinates of matrix, from matrix-vector product
! sigma    = the eigenvalue shift
!------------------------------------------------------------------------------

#ifdef USE_MULTI
      subroutine asigma_product_subtract_local(vect_in, vect_out, id, id2, sigma)

      integer, intent(in) :: id, id2
#ifdef USE_COMPLEX
      double complex, intent(in)    :: sigma
      double complex, intent(in)    :: vect_in(dm(id2)%d_dim)
      double complex, intent(inout) :: vect_out(dm(id)%d_dim)
      double complex sigmap, cfactor
      double complex :: zero = (0d0, 0d0)
#else
      double precision, intent(in)    :: sigma
      double precision, intent(in)    :: vect_in(dm(id2)%d_dim)
      double precision, intent(inout) :: vect_out(dm(id)%d_dim)
      double precision :: sigmap, cfactor = 1d0
      double precision :: zero = 0d0
#endif
      integer i, j, ii, jj, n, der, eq, var, power
      integer a_index, v_index, eqloc, varloc

      ! Preliminary calculations: derivative(s) of vect
      dm(id2)%vect_der(1:dm(id2)%d_dim,dmat(id2)%der_min:dmat(id2)%der_max) = zero
      do var=1,dm(id2)%nvar
        do der=dm(id2)%var_der_min(var),dm(id2)%var_der_max(var)
          do i=1,grd(id2)%nr
            do j=1,nt
              v_index = dm(id2)%ivar(var,i,j)
              do ii=max(1,i-dmat(id2)%lbder(der)),min(grd(id2)%nr,i+dmat(id2)%ubder(der))
                a_index = dm(id2)%ivar(var,ii,j)
                dm(id2)%vect_der(v_index,der) = dm(id2)%vect_der(v_index,der) + &
                               dmat(id2)%derive(i,ii,der)*vect_in(a_index)
              enddo
            enddo
          enddo
        enddo
      enddo

      ! Boundary conditions:
      do n=1,idm(id,id2)%natbc
        power = idm(id,id2)%atbci(1,n)
        der = idm(id,id2)%atbci(2,n)
        eq  = idm(id,id2)%atbci(3,n)
        var = idm(id,id2)%atbci(4,n)
        eqloc  = idm(id,id2)%atbci(5,n)
        varloc = idm(id,id2)%atbci(6,n)
#ifdef USE_COMPLEX
        if (idm(id,id2)%atbci(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do j=1,nt
          v_index = dm(id)%ieq(eq,eqloc,j)
          a_index = dm(id2)%ivar(var,varloc,j)
          vect_out(v_index) = vect_out(v_index) - &
            sigmap*idm(id,id2)%atbc(j,n)*dm(id2)%vect_der(a_index,der)
        enddo
      enddo

      do n=1,idm(id,id2)%nattbc
        power = idm(id,id2)%attbci(1,n)
        der = idm(id,id2)%attbci(2,n)
        eq  = idm(id,id2)%attbci(3,n)
        var = idm(id,id2)%attbci(4,n)
        eqloc  = idm(id,id2)%attbci(5,n)
        varloc = idm(id,id2)%attbci(6,n)
#ifdef USE_COMPLEX
        if (idm(id,id2)%attbci(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do j=1,nt
          v_index = dm(id)%ieq(eq,eqloc,j)
          do jj=1,nt
            a_index = dm(id2)%ivar(var,varloc,jj)
            vect_out(v_index) = vect_out(v_index) - &
              sigmap*idm(id,id2)%attbc(j,jj,n)*dm(id2)%vect_der(a_index,der)
          enddo
        enddo
      enddo

      if (id.ne.id2) return

      ! The equations:
      do n=1,dm(id)%nas
        power = dm(id)%asi(1,n)
        der = dm(id)%asi(2,n)
        eq  = dm(id)%asi(3,n)
        var = dm(id)%asi(4,n)
#ifdef USE_COMPLEX
        if (dm(id)%asi(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do i=1,grd(id)%nr
          do j=1,nt
            v_index = dm(id)%ieq(eq,i,j)
            if (.not.dm(id)%bc_flag(v_index)) then
              a_index = dm(id)%ivar(var,i,j)
              vect_out(v_index) = vect_out(v_index) -    &
                sigmap*dm(id)%as(n)*dm(id)%vect_der(a_index,der)
            endif
          enddo
        enddo
      enddo

      do n=1,dm(id)%nart
        power = dm(id)%arti(1,n)
        der = dm(id)%arti(2,n)
        eq  = dm(id)%arti(3,n)
        var = dm(id)%arti(4,n)
#ifdef USE_COMPLEX
        if (dm(id)%arti(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do i=1,grd(id)%nr
          do j=1,nt
            v_index = dm(id)%ieq(eq,i,j)
            if (.not.dm(id)%bc_flag(v_index)) then
              a_index = dm(id)%ivar(var,i,j)
              vect_out(v_index) = vect_out(v_index) - &
                sigmap*dm(id)%art(i,j,n)*dm(id)%vect_der(a_index,der)
            endif
          enddo
        enddo
      enddo

      do n=1,dm(id)%nartt
        power = dm(id)%artti(1,n)
        der = dm(id)%artti(2,n)
        eq  = dm(id)%artti(3,n)
        var = dm(id)%artti(4,n)
#ifdef USE_COMPLEX
        if (dm(id)%artti(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do i=1,grd(id)%nr
          do j=1,nt
            v_index = dm(id)%ieq(eq,i,j)
            if (.not.dm(id)%bc_flag(v_index)) then
              do jj=1,nt
                a_index = dm(id)%ivar(var,i,jj)
                vect_out(v_index) = vect_out(v_index) - &
                  sigmap*dm(id)%artt(i,j,jj,n)*dm(id)%vect_der(a_index,der)
              enddo
            endif
          enddo
        enddo
      enddo

      end subroutine asigma_product_subtract_local
#endif

!------------------------------------------------------------------------------
! This does:
!   mat_out <-  mat_out - asigma(id,id-1)*mat_in
! using a sparse multiplication method.
!------------------------------------------------------------------------------
! Variables:
!
! mat_in   = matrix which should correspond to:
!               asigma(id-1,id-1)^(-1)*asigma(id-1,id)
!            and is compressesd horizontally (only those columns given by
!            the h_bc_range(:) array are defined)
! mat_out  = matrix which receives the result (see above formula), and
!            should, in the end, correspond to tilde(asigma(id,id)).
!            mat_out is a FULL matrix
! id       = the number of the domain
! sigma    = the eigenvalue shift
!------------------------------------------------------------------------------
! IMPORTANT: id must be greater than 1 (otherwise the declarations will go out
!            of array bounds)
!------------------------------------------------------------------------------

#ifdef USE_MULTI
      subroutine asigma_full_product_subtract_local_hcomp(mat_in, &
                        mat_out, id, sigma)

      integer, intent(in)           :: id
#ifdef USE_COMPLEX
      double complex, intent(in)    :: mat_in(dm(id-1)%d_dim,idm(id-1,id)%n_h_bc)
      double complex, intent(inout) :: mat_out(dm(id)%d_dim,dm(id)%d_dim)
      double complex, intent(in)    :: sigma
      double complex sigmap, cfactor
#else
      double precision, intent(in)    :: mat_in(dm(id-1)%d_dim,idm(id-1,id)%n_h_bc)
      double precision, intent(inout) :: mat_out(dm(id)%d_dim,dm(id)%d_dim)
      double precision, intent(in)    :: sigma
      double precision :: sigmap, cfactor = 1d0
#endif
      integer i, j, ii, jj, n, k, kk, der, eq, var, a_index, v_index
      integer eqloc, varloc, power

      do k=1,idm(id-1,id)%n_h_bc
        ! Preliminary calculations: derivative(s) of vect
        dm(id-1)%vect_der(1:dm(id-1)%d_dim,dmat(id-1)%der_min:dmat(id-1)%der_max) = (0d0,0d0)
        do var=1,dm(id-1)%nvar
          do der=dm(id-1)%var_der_min(var),dm(id-1)%var_der_max(var)
            do i=1,grd(id-1)%nr
              do j=1,nt
                v_index = dm(id-1)%ivar(var,i,j)
                do ii=max(1,i-dmat(id-1)%lbder(der)),min(grd(id-1)%nr,i+dmat(id-1)%ubder(der))
                  a_index = dm(id-1)%ivar(var,ii,j)
                  dm(id-1)%vect_der(v_index,der) = dm(id-1)%vect_der(v_index,der) + &
                                 dmat(id-1)%derive(i,ii,der)*mat_in(a_index,k)
                enddo
              enddo
            enddo
          enddo
        enddo

        ! Boundary conditions:
        kk = idm(id-1,id)%h_bc_range(k)
        do n=1,idm(id,id-1)%natbc
          power  = idm(id,id-1)%atbci(1,n)
          der    = idm(id,id-1)%atbci(2,n)
          eq     = idm(id,id-1)%atbci(3,n)
          var    = idm(id,id-1)%atbci(4,n)
          eqloc  = idm(id,id-1)%atbci(5,n)
          varloc = idm(id,id-1)%atbci(6,n)
#ifdef USE_COMPLEX
          if (idm(id,id-1)%atbci(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          sigmap = cfactor*sigma**power
          do j=1,nt
            v_index = dm(id)%ieq(eq,eqloc,j)
            a_index = dm(id-1)%ivar(var,varloc,j)
            mat_out(v_index,kk) = mat_out(v_index,kk) - &
              sigmap*idm(id,id-1)%atbc(j,n)*dm(id-1)%vect_der(a_index,der)
          enddo
        enddo

        do n=1,idm(id,id-1)%nattbc
          power  = idm(id,id-1)%attbci(1,n)
          der    = idm(id,id-1)%attbci(2,n)
          eq     = idm(id,id-1)%attbci(3,n)
          var    = idm(id,id-1)%attbci(4,n)
          eqloc  = idm(id,id-1)%attbci(5,n)
          varloc = idm(id,id-1)%attbci(6,n)
#ifdef USE_COMPLEX
          if (idm(id,id-1)%attbci(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          sigmap = cfactor*sigma**power
          do j=1,nt
            v_index = dm(id)%ieq(eq,eqloc,j)
            do jj=1,nt
              a_index = dm(id-1)%ivar(var,varloc,jj)
              mat_out(v_index,kk) = mat_out(v_index,kk) - &
                sigmap*idm(id,id-1)%attbc(j,jj,n)*dm(id-1)%vect_der(a_index,der)
            enddo
          enddo
        enddo
      enddo

      end subroutine asigma_full_product_subtract_local_hcomp
#endif

!------------------------------------------------------------------------------
! This does:
!   mat_out <-  mat_out - asigma(id,id-1)*mat_in
! using a sparse multiplication method.
!------------------------------------------------------------------------------
! Variables:
!
! mat_in   = matrix which should correspond to:
!               asigma(id-1,id-1)^(-1)*asigma(id-1,id)
!            and is compressesd horizontally (only those columns given by
!            the h_bc_range(:) are defined)
! mat_out  = matrix which receives the result (see above formula), and
!            should, in the end, correspond to tilde(asigma(id,id))
!            mat_out is a BAND matrix
! id       = the number of the domain
!------------------------------------------------------------------------------
! IMPORTANT: id must be greater than 1 (otherwise the declarations will go out
!            of array bounds)
!------------------------------------------------------------------------------

#ifdef USE_MULTI
      subroutine asigma_band_product_subtract_local_hcomp(mat_in, &
                 mat_out, id, sigma)

      integer, intent(in)           :: id
#ifdef USE_COMPLEX
      double complex, intent(in)    :: mat_in(dm(id-1)%d_dim,idm(id-1,id)%n_h_bc)
      double complex, intent(inout) :: mat_out(2*dm(id)%kl+dm(id)%ku &
                                                 +1,dm(id)%d_dim)
      double complex, intent(in)    :: sigma
      double complex sigmap, cfactor
#else
      double precision, intent(in)    :: mat_in(dm(id-1)%d_dim,idm(id-1,id)%n_h_bc)
      double precision, intent(inout) :: mat_out(2*dm(id)%kl+dm(id)%ku &
                                                   +1,dm(id)%d_dim)
      double precision, intent(in)    :: sigma
      double precision :: sigmap, cfactor = 1d0
#endif
      integer i, j, ii, jj, n, k, kk, der, eq, var, a_index, v_index, pos
      integer eqloc, varloc, power

      pos = dm(id)%kl + dm(id)%ku + 1

      do k=1,idm(id-1,id)%n_h_bc
        ! Preliminary calculations: derivative(s) of vect
        dm(id-1)%vect_der(1:dm(id-1)%d_dim,dmat(id-1)%der_min:dmat(id-1)%der_max) = (0d0,0d0)
        do var=1,dm(id-1)%nvar
          do der=dm(id-1)%var_der_min(var),dm(id-1)%var_der_max(var)
            do i=1,grd(id-1)%nr
              do j=1,nt
                v_index = dm(id-1)%ivar(var,i,j)
                do ii=max(1,i-dmat(id-1)%lbder(der)),min(grd(id-1)%nr,i+dmat(id-1)%ubder(der))
                  a_index = dm(id-1)%ivar(var,ii,j)
                  dm(id-1)%vect_der(v_index,der) = dm(id-1)%vect_der(v_index,der) + &
                                 dmat(id-1)%derive(i,ii,der)*mat_in(a_index,k)
                enddo
              enddo
            enddo
          enddo
        enddo

        ! Boundary conditions:
        kk = idm(id-1,id)%h_bc_range(k)
        do n=1,idm(id,id-1)%natbc
          power  = idm(id,id-1)%atbci(1,n)
          der    = idm(id,id-1)%atbci(2,n)
          eq     = idm(id,id-1)%atbci(3,n)
          var    = idm(id,id-1)%atbci(4,n)
          eqloc  = idm(id,id-1)%atbci(5,n)
          varloc = idm(id,id-1)%atbci(6,n)
#ifdef USE_COMPLEX
          if (idm(id,id-1)%atbci(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          sigmap = cfactor*sigma**power
          do j=1,nt
            v_index = dm(id)%ieq(eq,eqloc,j)
            a_index = dm(id-1)%ivar(var,varloc,j)
            mat_out(pos+v_index-kk,kk) = mat_out(pos+v_index-kk,kk) - &
              sigmap*idm(id,id-1)%atbc(j,n)*dm(id-1)%vect_der(a_index,der)
          enddo
        enddo

        do n=1,idm(id,id-1)%nattbc
          power  = idm(id,id-1)%attbci(1,n)
          der    = idm(id,id-1)%attbci(2,n)
          eq     = idm(id,id-1)%attbci(3,n)
          var    = idm(id,id-1)%attbci(4,n)
          eqloc  = idm(id,id-1)%attbci(5,n)
          varloc = idm(id,id-1)%attbci(6,n)
#ifdef USE_COMPLEX
          if (idm(id,id-1)%attbci(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          sigmap = cfactor*sigma**power
          do j=1,nt
            v_index = dm(id)%ieq(eq,eqloc,j)
            do jj=1,nt
              a_index = dm(id-1)%ivar(var,varloc,jj)
              mat_out(pos+v_index-kk,kk) = mat_out(pos+v_index-kk,kk) - &
                sigmap*idm(id,id-1)%attbc(j,jj,n)*dm(id-1)%vect_der(a_index,der)
            enddo
          enddo
        enddo
      enddo

      end subroutine asigma_band_product_subtract_local_hcomp
#endif

!------------------------------------------------------------------------------
! This subroutine performs the matrix vector product
! conjg(transpose(a(power))).x for the whole system, where a(power) comes
! from the eigenvalue problem:
!      w^2.a(2).x + w.a(1).x + a(0).x = 0.
!------------------------------------------------------------------------------
! List of variables:
!
! vect_in(1:a_dim)  = vecteur undergoing the matrix product
! vect_out(1:a_dim) = vecteur receiving the result from the matrix-vector
!                     product
! power             = this tells which matrix to use (it corresponds to the
!                     power of the eigenvalue in front)
!------------------------------------------------------------------------------

#ifdef USE_MULTI
      subroutine a_product_total_transpose(vect_in,vect_out,power)

#ifdef USE_COMPLEX
      double complex, intent(in)  :: vect_in(a_dim)
      double complex, intent(out) :: vect_out(a_dim)
      double complex cfactor
      double complex :: zero = (0d0, 0d0)
#else
      double precision, intent(in)  :: vect_in(a_dim)
      double precision, intent(out) :: vect_out(a_dim)
      double precision :: cfactor = 1d0
      double precision :: zero = 0d0
#endif
      integer, intent(in)         :: power
      integer id,id2,i,j,ii,jj,n,der,eq,var,eqloc,varloc,a_index,v_index

      ! Very important: multiplications by derivatives only takes place
      ! at the end since the transpose operation swaps the order of
      ! matrix multiplications.

      ! initialisation:
      do id=1,ndomains
        dm(id)%vect_der(1:dm(id)%d_dim,dmat(id)%der_min:dmat(id)%der_max) = zero
      enddo

      ! The equations:
      do id=1,ndomains
        do n=1,dm(id)%nas
          if (dm(id)%asi(1,n).eq.power) then
            der = dm(id)%asi(2,n)
            eq  = dm(id)%asi(3,n)
            var = dm(id)%asi(4,n)
#ifdef USE_COMPLEX
            if (dm(id)%asi(7,n).eq.0) then
              cfactor = (1d0,0d0)
            else
              cfactor = (0d0,-1d0)
            endif
#endif
            do i=1,grd(id)%nr
              do j=1,nt
                a_index = dm(id)%ieq(eq,i,j)
                if (.not.dm(id)%bc_flag(a_index)) then
                  a_index = a_index + dm(id)%offset
                  v_index = dm(id)%ivar(var,i,j)
                  dm(id)%vect_der(v_index,der) =   &
                    dm(id)%vect_der(v_index,der) + &
                    cfactor*dm(id)%as(n)*vect_in(a_index)
                endif
              enddo
            enddo
          endif
        enddo

        do n=1,dm(id)%nart
          if (dm(id)%arti(1,n).eq.power) then
            der = dm(id)%arti(2,n)
            eq  = dm(id)%arti(3,n)
            var = dm(id)%arti(4,n)
#ifdef USE_COMPLEX
            if (dm(id)%arti(7,n).eq.0) then
              cfactor = (1d0,0d0)
            else
              cfactor = (0d0,-1d0)
            endif
#endif
            do i=1,grd(id)%nr
              do j=1,nt
                a_index = dm(id)%ieq(eq,i,j)
                if (.not.dm(id)%bc_flag(a_index)) then
                  a_index = a_index + dm(id)%offset
                  v_index = dm(id)%ivar(var,i,j)
                  dm(id)%vect_der(v_index,der) =   &
                    dm(id)%vect_der(v_index,der) + &
                    cfactor*dm(id)%art(i,j,n)*vect_in(a_index)
                endif
              enddo
            enddo
          endif
        enddo

        do n=1,dm(id)%nartt
          if (dm(id)%artti(1,n).eq.power) then
            der = dm(id)%artti(2,n)
            eq  = dm(id)%artti(3,n)
            var = dm(id)%artti(4,n)
#ifdef USE_COMPLEX
            if (dm(id)%artti(7,n).eq.0) then
              cfactor = (1d0,0d0)
            else
              cfactor = (0d0,-1d0)
            endif
#endif
            do i=1,grd(id)%nr
              do j=1,nt
                a_index = dm(id)%ieq(eq,i,j)
                if (.not.dm(id)%bc_flag(a_index)) then
                  a_index = a_index+dm(id)%offset
                  do jj=1,nt
                    v_index = dm(id)%ivar(var,i,jj)
                    dm(id)%vect_der(v_index,der) =   &
                      dm(id)%vect_der(v_index,der) + &
                      cfactor*dm(id)%artt(i,j,jj,n)*vect_in(a_index)
                  enddo
                endif
              enddo
            enddo
          endif
        enddo

        ! Boundary conditions:
        do id2 = 1,ndomains
          do n=1,idm(id,id2)%natbc
            if (idm(id,id2)%atbci(1,n).eq.power) then
              der = idm(id,id2)%atbci(2,n)
              eq  = idm(id,id2)%atbci(3,n)
              var = idm(id,id2)%atbci(4,n)
              eqloc  = idm(id,id2)%atbci(5,n)
              varloc = idm(id,id2)%atbci(6,n)
#ifdef USE_COMPLEX
              if (idm(id,id2)%atbci(7,n).eq.0) then
                cfactor = (1d0,0d0)
              else
                cfactor = (0d0,-1d0)
              endif
#endif
              do j=1,nt
                a_index = dm(id)%offset+dm(id)%ieq(eq,eqloc,j)
                v_index = dm(id2)%ivar(var,varloc,j)
                dm(id2)%vect_der(v_index,der) =   &
                  dm(id2)%vect_der(v_index,der) + &
                  cfactor*idm(id,id2)%atbc(j,n)*vect_in(a_index)
              enddo
            endif
          enddo

          do n=1,idm(id,id2)%nattbc
            if (idm(id,id2)%attbci(1,n).eq.power) then
              der = idm(id,id2)%attbci(2,n)
              eq  = idm(id,id2)%attbci(3,n)
              var = idm(id,id2)%attbci(4,n)
              eqloc  = idm(id,id2)%attbci(5,n)
              varloc = idm(id,id2)%attbci(6,n)
#ifdef USE_COMPLEX
              if (idm(id,id2)%attbci(7,n).eq.0) then
                cfactor = (1d0,0d0)
              else
                cfactor = (0d0,-1d0)
              endif
#endif
              do j=1,nt
                a_index = dm(id)%offset+dm(id)%ieq(eq,eqloc,j)
                do jj=1,nt
                  v_index = dm(id2)%ivar(var,varloc,jj)
                  dm(id2)%vect_der(v_index,der) =   &
                    dm(id2)%vect_der(v_index,der) + &
                    cfactor*idm(id,id2)%attbc(j,jj,n)*vect_in(a_index)
                enddo
              enddo
            endif
          enddo
        enddo
      enddo

      ! Posterior calculations: transpose(mat_der)*vect
      vect_out(1:a_dim) = zero
      do id=1,ndomains
        do var=1,dm(id)%nvar
          do der=dm(id)%var_der_min(var),dm(id)%var_der_max(var)
            do j=1,nt
              do i=1,grd(id)%nr
                v_index = dm(id)%offset+dm(id)%ivar(var,i,j)
                do ii=max(1,i-dmat(id)%ubder(der)),min(grd(id)%nr,i+dmat(id)%lbder(der))
                  a_index = dm(id)%ivar(var,ii,j)
                  ! To understand this, you need to draw a picture...
                  vect_out(v_index) = vect_out(v_index) + &
                    dmat(id)%derive(ii,i,der)*dm(id)%vect_der(a_index,der)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      end subroutine a_product_total_transpose
#endif

!------------------------------------------------------------------------------
! This does:
!   vect_out(id2) <- vect_out(id2) - CONJG(TRANSPOSE(asigma(id,id2)))*vect_in(id)
! using a sparse multiplication method.
!------------------------------------------------------------------------------
! Variables:
!
! vect_in  = input vector (dim = dm(id)%d_dim)
! vect_out = output vector (dim = dm(id2)%d_dim)
! id, id2  = domain coordinates of matrix, from matrix-vector product
! sigma    = the eigenvalue shift
!------------------------------------------------------------------------------
! NOTE: The eigenvalue shift will be conjugated within this subroutine, so the
!       subroutine should be called with the eigenvalue shift and NOT its
!       conjugate.
!------------------------------------------------------------------------------

#ifdef USE_MULTI
      subroutine asigma_product_subtract_local_transpose(vect_in, &
                                        vect_out, id, id2, sigma)

      integer, intent(in)           :: id, id2
#ifdef USE_COMPLEX
      double complex, intent(in)    :: vect_in(dm(id)%d_dim)
      double complex, intent(inout) :: vect_out(dm(id2)%d_dim)
      double complex, intent(in)    :: sigma
      double complex sigmap, cfactor
      double complex :: zero = (0d0, 0d0)
#else
      double precision, intent(in)    :: vect_in(dm(id)%d_dim)
      double precision, intent(inout) :: vect_out(dm(id2)%d_dim)
      double precision, intent(in)    :: sigma
      double precision sigmap
      double precision :: zero = 0d0
#endif
      integer i, j, ii, jj, n, der, eq, var, power
      integer a_index, v_index, eqloc, varloc

      ! Initialisation
      dm(id2)%vect_der(1:dm(id2)%d_dim,dmat(id2)%der_min:dmat(id2)%der_max) = zero

      ! Boundary conditions:
      do n=1,idm(id,id2)%natbc
        power  = idm(id,id2)%atbci(1,n)
        der    = idm(id,id2)%atbci(2,n)
        eq     = idm(id,id2)%atbci(3,n)
        var    = idm(id,id2)%atbci(4,n)
        eqloc  = idm(id,id2)%atbci(5,n)
        varloc = idm(id,id2)%atbci(6,n)
#ifdef USE_COMPLEX
        if (idm(id,id2)%atbci(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,-1d0)
        endif
        sigmap = cfactor*dconjg(sigma**power)
#else
        sigmap = sigma**power
#endif
        do j=1,nt
          a_index = dm(id)%ieq(eq,eqloc,j)
          v_index = dm(id2)%ivar(var,varloc,j)
          dm(id2)%vect_der(v_index,der) = dm(id2)%vect_der(v_index,der) + &
            sigmap*idm(id,id2)%atbc(j,n)*vect_in(a_index)
        enddo
      enddo

      do n=1,idm(id,id2)%nattbc
        power  = idm(id,id2)%attbci(1,n)
        der    = idm(id,id2)%attbci(2,n)
        eq     = idm(id,id2)%attbci(3,n)
        var    = idm(id,id2)%attbci(4,n)
        eqloc  = idm(id,id2)%attbci(5,n)
        varloc = idm(id,id2)%attbci(6,n)
#ifdef USE_COMPLEX
        if (idm(id,id2)%attbci(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,-1d0)
        endif
        sigmap = cfactor*dconjg(sigma**power)
#else
        sigmap = sigma**power
#endif
        do j=1,nt
          a_index = dm(id)%ieq(eq,eqloc,j)
          do jj=1,nt
            v_index = dm(id2)%ivar(var,varloc,jj)
            dm(id2)%vect_der(v_index,der) = dm(id2)%vect_der(v_index,der) + &
              sigmap*idm(id,id2)%attbc(j,jj,n)*vect_in(a_index)
          enddo
        enddo
      enddo

      if (id.eq.id2) then

        ! The equations:
        do n=1,dm(id)%nas
          power  = dm(id)%asi(1,n)
          der    = dm(id)%asi(2,n)
          eq     = dm(id)%asi(3,n)
          var    = dm(id)%asi(4,n)
#ifdef USE_COMPLEX
          if (dm(id)%asi(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,-1d0)
          endif
          sigmap = cfactor*dconjg(sigma**power)
#else
          sigmap = sigma**power
#endif
          do i=1,grd(id)%nr
            do j=1,nt
              a_index = dm(id)%ieq(eq,i,j)
              if (.not.dm(id)%bc_flag(a_index)) then
                v_index = dm(id)%ivar(var,i,j)
                dm(id)%vect_der(v_index,der) = dm(id)%vect_der(v_index,der) + &
                  sigmap*dm(id)%as(n)*vect_in(a_index)
              endif
            enddo
          enddo
        enddo

        do n=1,dm(id)%nart
          power  = dm(id)%arti(1,n)
          der    = dm(id)%arti(2,n)
          eq     = dm(id)%arti(3,n)
          var    = dm(id)%arti(4,n)
#ifdef USE_COMPLEX
          if (dm(id)%arti(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,-1d0)
          endif
          sigmap = cfactor*dconjg(sigma**power)
#else
          sigmap = sigma**power
#endif
          do i=1,grd(id)%nr
            do j=1,nt
              a_index = dm(id)%ieq(eq,i,j)
              if (.not.dm(id)%bc_flag(a_index)) then
                v_index = dm(id)%ivar(var,i,j)
                dm(id)%vect_der(v_index,der) = dm(id)%vect_der(v_index,der) + &
                  sigmap*dm(id)%art(i,j,n)*vect_in(a_index)
              endif
            enddo
          enddo
        enddo

        do n=1,dm(id)%nartt
          power  = dm(id)%artti(1,n)
          der    = dm(id)%artti(2,n)
          eq     = dm(id)%artti(3,n)
          var    = dm(id)%artti(4,n)
#ifdef USE_COMPLEX
          if (dm(id)%artti(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,-1d0)
          endif
          sigmap = cfactor*dconjg(sigma**power)
#else
          sigmap = sigma**power
#endif
          do i=1,grd(id)%nr
            do j=1,nt
              a_index = dm(id)%ieq(eq,i,j)
              if (.not.dm(id)%bc_flag(a_index)) then
                do jj=1,nt
                  v_index = dm(id)%ivar(var,i,jj)
                  dm(id)%vect_der(v_index,der) = dm(id)%vect_der(v_index,der) + &
                    sigmap*dm(id)%artt(i,j,jj,n)*vect_in(a_index)
                enddo
              endif
            enddo
          enddo
        enddo

      endif

      ! Finish the calculations: derivative(s) of vect
      do var=1,dm(id2)%nvar
        do der=dm(id2)%var_der_min(var),dm(id2)%var_der_max(var)
          do i=1,grd(id2)%nr
            do j=1,nt
              v_index = dm(id2)%ivar(var,i,j)
              do ii=max(1,i-dmat(id2)%ubder(der)),min(grd(id2)%nr,i+dmat(id2)%lbder(der))
                a_index = dm(id2)%ivar(var,ii,j)
                ! To understand this, you need to draw a picture...
                vect_out(v_index) = vect_out(v_index) - &
                  dmat(id2)%derive(ii,i,der)*dm(id2)%vect_der(a_index,der)
              enddo
            enddo
          enddo
        enddo
      enddo

      end subroutine asigma_product_subtract_local_transpose
#endif

!------------------------------------------------------------------------------
! This subroutine creates the matrix asigma using FULL storage:
!      asigma(id,id2) = a(power=0,id,id2) + sigma.a(power=1,id,id2)
!                     + sigma^2.a(power=2,id,id2) +  ...
!------------------------------------------------------------------------------
! List of variables:
!
! sigma  = value for sigma
! asigma = matrix which will contain asigma
! id,id2 = domain coordinates for asigma
!------------------------------------------------------------------------------
#ifdef USE_MULTI
      subroutine make_asigma_full_local(sigma, asigma, id, id2, ierr)

      integer, intent(in)         :: id, id2
#ifdef USE_COMPLEX
      double complex, intent(in)  :: sigma
      double complex, intent(out) :: asigma(dm(id)%d_dim,dm(id2)%d_dim)
      double complex sigmap, temp, cfactor
      double complex :: zero = (0d0, 0d0)
#else
      double precision, intent(in)  :: sigma
      double precision, intent(out) :: asigma(dm(id)%d_dim,dm(id2)%d_dim)
      double precision :: sigmap, temp, cfactor = 1d0
      double precision :: zero = 0d0
#endif
      integer n, i, j, ii, jj, l, c, power, der, eq, var, eqloc, varloc
      integer, intent(out) :: ierr

      asigma = zero

      ! Boundary conditions:
      do n=1,idm(id,id2)%natbc
        power  = idm(id,id2)%atbci(1,n)
        der    = idm(id,id2)%atbci(2,n)
        eq     = idm(id,id2)%atbci(3,n)
        var    = idm(id,id2)%atbci(4,n)
        eqloc  = idm(id,id2)%atbci(5,n)
        varloc = idm(id,id2)%atbci(6,n)
#ifdef USE_COMPLEX
        if (idm(id,id2)%atbci(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do ii=max(1,varloc-dmat(id2)%lbder(der)),min(grd(id2)%nr,varloc+dmat(id2)%ubder(der))
          temp = sigmap*dmat(id2)%derive(varloc,ii,der)
          do j=1,nt
            l = dm(id)%ieq(eq,eqloc,j)
            c = dm(id2)%ivar(var,ii,j)
            asigma(l,c) = asigma(l,c)+temp*idm(id,id2)%atbc(j,n)
          enddo
        enddo
      enddo

      do n=1,idm(id,id2)%nattbc
        power  = idm(id,id2)%attbci(1,n)
        der    = idm(id,id2)%attbci(2,n)
        eq     = idm(id,id2)%attbci(3,n)
        var    = idm(id,id2)%attbci(4,n)
        eqloc  = idm(id,id2)%attbci(5,n)
        varloc = idm(id,id2)%attbci(6,n)
#ifdef USE_COMPLEX
        if (idm(id,id2)%attbci(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do ii=max(1,varloc-dmat(id2)%lbder(der)),min(grd(id2)%nr,varloc+dmat(id2)%ubder(der))
          temp = sigmap*dmat(id2)%derive(varloc,ii,der)
          do j=1,nt
            l = dm(id)%ieq(eq,eqloc,j)
            do jj=1,nt
              c = dm(id2)%ivar(var,ii,jj)
              asigma(l,c) = asigma(l,c)+temp*idm(id,id2)%attbc(j,jj,n)
            enddo
          enddo
        enddo
      enddo

      if (id.ne.id2) return

      ! The equations:
      do n=1,dm(id)%nas
        power  = dm(id)%asi(1,n)
        der    = dm(id)%asi(2,n)
        eq     = dm(id)%asi(3,n)
        var    = dm(id)%asi(4,n)
#ifdef USE_COMPLEX
        if (dm(id)%asi(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        temp   = sigmap*dm(id)%as(n)
        do i=1,grd(id)%nr
          do j=1,nt
            l = dm(id)%ieq(eq,i,j)
            if (.not.dm(id)%bc_flag(l)) then
              do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                c = dm(id)%ivar(var,ii,j)
                asigma(l,c) = asigma(l,c)+temp*dmat(id)%derive(i,ii,der)
              enddo
            endif
          enddo
        enddo
      enddo

      do n=1,dm(id)%nart
        power  = dm(id)%arti(1,n)
        der    = dm(id)%arti(2,n)
        eq     = dm(id)%arti(3,n)
        var    = dm(id)%arti(4,n)
#ifdef USE_COMPLEX
        if (dm(id)%arti(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do i=1,grd(id)%nr
          do j=1,nt
            l = dm(id)%ieq(eq,i,j)
            if (.not.dm(id)%bc_flag(l)) then
              temp = sigmap*dm(id)%art(i,j,n)
              do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                c = dm(id)%ivar(var,ii,j)
                asigma(l,c) = asigma(l,c)+temp*dmat(id)%derive(i,ii,der)
              enddo
            endif
          enddo
        enddo
      enddo

      do n=1,dm(id)%nartt
        power  = dm(id)%artti(1,n)
        der    = dm(id)%artti(2,n)
        eq     = dm(id)%artti(3,n)
        var    = dm(id)%artti(4,n)
#ifdef USE_COMPLEX
        if (dm(id)%artti(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do i=1,grd(id)%nr
          do j=1,nt
            l = dm(id)%ieq(eq,i,j)
            if (.not.dm(id)%bc_flag(l)) then
              do jj=1,nt
                temp = sigmap*dm(id)%artt(i,j,jj,n)
                do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                  c = dm(id)%ivar(var,ii,jj)
                  asigma(l,c) = asigma(l,c)+temp*dmat(id)%derive(i,ii,der)
                enddo
              enddo
            endif
          enddo
        enddo
      enddo

      if (dump_asigma) then
          call dump_asigma_matrix(asigma, "asigma_full_local", ierr)
          if (ierr /= 0) return
      endif

      end subroutine make_asigma_full_local
#else
      subroutine make_asigma_full(sigma, asigma, ierr)

      double precision sigma
      double precision asigma(a_dim,a_dim)
      double precision sigmap, temp
      integer n, i, j, ii, jj, l, c, power, der, eq, var, eqloc, varloc
      integer, intent(out) :: ierr

      asigma = 0d0
#ifdef USE_1D
            ! The equations:
      do n=1, dm(1)%nas
        power = dm(1)%asi(1,n)
        der = dm(1)%asi(2,n)
        eq  = dm(1)%asi(3,n)
        var = dm(1)%asi(4,n)
        sigmap = sigma**power
        do i=1,nr
          l = dm(1)%ieq(eq,i)
          if (.not.dm(1)%bc_flag(l)) then
            do ii=max(1,i-dmat(1)%lbder(der)),min(nr,i+dmat(1)%ubder(der))
              c = dm(1)%ivar(var,ii)
              asigma(l,c) = asigma(l,c)+sigmap*dm(1)%as(n)*dmat(1)%derive(i,ii,der)
            enddo
          endif
        enddo
      enddo

      do n=1, dm(1)%nar
        power = dm(1)%ari(1,n)
        der = dm(1)%ari(2,n)
        eq  = dm(1)%ari(3,n)
        var = dm(1)%ari(4,n)
        sigmap = sigma**power
        do i=1,nr
          l = dm(1)%ieq(eq,i)
          if (.not.dm(1)%bc_flag(l)) then
            do ii=max(1,i-dmat(1)%lbder(der)),min(nr,i+dmat(1)%ubder(der))
              c = dm(1)%ivar(var,ii)
              asigma(l,c) = asigma(l,c)+sigmap*dm(1)%ar(i,n)*dmat(1)%derive(i,ii,der)
            enddo
          endif
        enddo
      enddo

      ! Boundary conditions:
      do n=1, dm(1)%nasbc
        power = dm(1)%asbci(1,n)
        der = dm(1)%asbci(2,n)
        eq  = dm(1)%asbci(3,n)
        var = dm(1)%asbci(4,n)
        eqloc  = dm(1)%asbci(5,n)
        varloc = dm(1)%asbci(6,n)
        sigmap = sigma**power
        l = dm(1)%ieq(eq,eqloc)
        do ii=max(1,varloc-dmat(1)%lbder(der)),min(nr,varloc+dmat(1)%ubder(der))
          c = dm(1)%ivar(var,ii)
          asigma(l,c) = asigma(l,c)+sigmap*dm(1)%asbc(n)*dmat(1)%derive(varloc,ii,der)
        enddo
      enddo

#else
      ! The equations:
      do n=1,dm(1)%nas
        power = dm(1)%asi(1,n)
        der = dm(1)%asi(2,n)
        eq  = dm(1)%asi(3,n)
        var = dm(1)%asi(4,n)
        sigmap = sigma**power
        temp = sigmap*dm(1)%as(n)
        do i=1,grd(1)%nr
          do j=1,nt
            l = dm(1)%ieq(eq,i,j)
            if (.not.dm(1)%bc_flag(l)) then
              do ii=max(1,i-dmat(1)%lbder(der)),min(grd(1)%nr,i+dmat(1)%ubder(der))
                c = dm(1)%ivar(var,ii,j)
                asigma(l,c) = asigma(l,c)+temp*dmat(1)%derive(i,ii,der)
              enddo
            endif
          enddo
        enddo
      enddo

      do n=1,dm(1)%nart
        power = dm(1)%arti(1,n)
        der = dm(1)%arti(2,n)
        eq  = dm(1)%arti(3,n)
        var = dm(1)%arti(4,n)
        sigmap = sigma**power
        do i=1,grd(1)%nr
          do j=1,nt
            l = dm(1)%ieq(eq,i,j)
            if (.not.dm(1)%bc_flag(l)) then
              temp = sigmap*dm(1)%art(i,j,n)
              do ii=max(1,i-dmat(1)%lbder(der)),min(grd(1)%nr,i+dmat(1)%ubder(der))
                c = dm(1)%ivar(var,ii,j)
                asigma(l,c) = asigma(l,c)+temp*dmat(1)%derive(i,ii,der)
              enddo
            endif
          enddo
        enddo
      enddo

      do n=1,dm(1)%nartt
        power = dm(1)%artti(1,n)
        der = dm(1)%artti(2,n)
        eq  = dm(1)%artti(3,n)
        var = dm(1)%artti(4,n)
        sigmap = sigma**power
        do i=1,grd(1)%nr
          do j=1,nt
            l = dm(1)%ieq(eq,i,j)
            if (.not.dm(1)%bc_flag(l)) then
              do jj=1,nt
                temp = sigmap*dm(1)%artt(i,j,jj,n)
                do ii=max(1,i-dmat(1)%lbder(der)),min(grd(1)%nr,i+dmat(1)%ubder(der))
                  c = dm(1)%ivar(var,ii,jj)
                  asigma(l,c) = asigma(l,c)+temp*dmat(1)%derive(i,ii,der)
                enddo
              enddo
            endif
          enddo
        enddo
      enddo

      ! Boundary conditions:
      do n=1,idm(1, 1)%natbc
        power = idm(1, 1)%atbci(1,n)
        der = idm(1, 1)%atbci(2,n)
        eq  = idm(1, 1)%atbci(3,n)
        var = idm(1, 1)%atbci(4,n)
        eqloc  = idm(1, 1)%atbci(5,n)
        varloc = idm(1, 1)%atbci(6,n)
        sigmap = sigma**power
        do ii=max(1,varloc-dmat(1)%lbder(der)),min(grd(1)%nr,varloc+dmat(1)%ubder(der))
          temp = sigmap*dmat(1)%derive(varloc,ii,der)
          do j=1,nt
            l = dm(1)%ieq(eq,eqloc,j)
            c = dm(1)%ivar(var,ii,j)
            asigma(l,c) = asigma(l,c)+temp*idm(1, 1)%atbc(j,n)
          enddo
        enddo
      enddo

      do n=1,idm(1, 1)%nattbc
        power = idm(1, 1)%attbci(1,n)
        der = idm(1, 1)%attbci(2,n)
        eq  = idm(1, 1)%attbci(3,n)
        var = idm(1, 1)%attbci(4,n)
        eqloc  = idm(1, 1)%attbci(5,n)
        varloc = idm(1, 1)%attbci(6,n)
        sigmap = sigma**power
        do ii=max(1,varloc-dmat(1)%lbder(der)),min(grd(1)%nr,varloc+dmat(1)%ubder(der))
          temp = sigmap*dmat(1)%derive(varloc,ii,der)
          do j=1,nt
            l = dm(1)%ieq(eq,eqloc,j)
            do jj=1,nt
              c = dm(1)%ivar(var,ii,jj)
              asigma(l,c) = asigma(l,c)+temp*idm(1, 1)%attbc(j,jj,n)
            enddo
          enddo
        enddo
      enddo
#endif

      if (dump_asigma) then
          call dump_asigma_matrix(asigma, "asigma_full", ierr)
          if (ierr /= 0) return
      endif

      end subroutine make_asigma_full
#endif

!------------------------------------------------------------------------------
! This subroutine creates the TOTAL matrix asigma using FULL storage:
!      asigma = a(power=0) + sigma.a(power=1) + sigma^2.a(power=2) +  ...
!------------------------------------------------------------------------------
! List of variables:
!
! sigma  = value for sigma
! asigma = matrix which will contain asigma
!------------------------------------------------------------------------------
#ifndef USE_1D
      subroutine make_asigma_full_total(sigma, asigma, ierr)

#ifdef USE_COMPLEX
      double complex, intent(in)  :: sigma
      double complex, intent(out) :: asigma(a_dim,a_dim)
      double complex sigmap, temp, cfactor
      double complex :: zero = (0d0, 0d0)
#else
      double precision, intent(in)  :: sigma
      double precision, intent(out) :: asigma(a_dim,a_dim)
      double precision :: sigmap, temp, cfactor = 1d0
      double precision :: zero = 0d0
#endif
      integer n, i, j, ii, jj, l, c, power, der, eq, var, eqloc, varloc
      integer id, id2
      integer, intent(out) :: ierr

      asigma = zero

      ! Boundary conditions:
      do id=1,ndomains
        do id2=1,ndomains
          do n=1,idm(id,id2)%natbc
            power  = idm(id,id2)%atbci(1,n)
            der    = idm(id,id2)%atbci(2,n)
            eq     = idm(id,id2)%atbci(3,n)
            var    = idm(id,id2)%atbci(4,n)
            eqloc  = idm(id,id2)%atbci(5,n)
            varloc = idm(id,id2)%atbci(6,n)
#ifdef USE_COMPLEX
            if (idm(id,id2)%atbci(7,n).eq.0) then
              cfactor = (1d0,0d0)
            else
              cfactor = (0d0,1d0)
            endif
#endif
            sigmap = cfactor*sigma**power
            do ii=max(1,varloc-dmat(id2)%lbder(der)),min(grd(id2)%nr,varloc+dmat(id2)%ubder(der))
              temp = sigmap*dmat(id2)%derive(varloc,ii,der)
              do j=1,nt
                l = dm(id)%ieq(eq,eqloc,j)+dm(id)%offset
                c = dm(id2)%ivar(var,ii,j)+dm(id2)%offset
                asigma(l,c) = asigma(l,c)+temp*idm(id,id2)%atbc(j,n)
              enddo
            enddo
          enddo

          do n=1,idm(id,id2)%nattbc
            power  = idm(id,id2)%attbci(1,n)
            der    = idm(id,id2)%attbci(2,n)
            eq     = idm(id,id2)%attbci(3,n)
            var    = idm(id,id2)%attbci(4,n)
            eqloc  = idm(id,id2)%attbci(5,n)
            varloc = idm(id,id2)%attbci(6,n)
#ifdef USE_COMPLEX
            if (idm(id,id2)%attbci(7,n).eq.0) then
              cfactor = (1d0,0d0)
            else
              cfactor = (0d0,1d0)
            endif
#endif
            sigmap = cfactor*sigma**power
            do ii=max(1,varloc-dmat(id2)%lbder(der)),min(grd(id2)%nr,varloc+dmat(id2)%ubder(der))
              temp = sigmap*dmat(id2)%derive(varloc,ii,der)
              do j=1,nt
                l = dm(id)%ieq(eq,eqloc,j)+dm(id)%offset
                do jj=1,nt
                  c = dm(id2)%ivar(var,ii,jj)+dm(id2)%offset
                  asigma(l,c) = asigma(l,c)+temp*idm(id,id2)%attbc(j,jj,n)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      ! The equations:
      do id=1,ndomains
        do n=1,dm(id)%nas
          power  = dm(id)%asi(1,n)
          der    = dm(id)%asi(2,n)
          eq     = dm(id)%asi(3,n)
          var    = dm(id)%asi(4,n)
#ifdef USE_COMPLEX
          if (dm(id)%asi(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          sigmap = cfactor*sigma**power
          temp   = sigmap*dm(id)%as(n)
          do i=1,grd(id)%nr
            do j=1,nt
              l = dm(id)%ieq(eq,i,j)
              if (.not.dm(id)%bc_flag(l)) then
                l = l + dm(id)%offset
                do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                  c = dm(id)%ivar(var,ii,j)+dm(id)%offset
                  asigma(l,c) = asigma(l,c)+temp*dmat(id)%derive(i,ii,der)
                enddo
              endif
            enddo
          enddo
        enddo

        do n=1,dm(id)%nart
          power  = dm(id)%arti(1,n)
          der    = dm(id)%arti(2,n)
          eq     = dm(id)%arti(3,n)
          var    = dm(id)%arti(4,n)
#ifdef USE_COMPLEX
          if (dm(id)%arti(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          sigmap = cfactor*sigma**power
          do i=1,grd(id)%nr
            do j=1,nt
              l = dm(id)%ieq(eq,i,j)
              if (.not.dm(id)%bc_flag(l)) then
                l = l + dm(id)%offset
                temp = sigmap*dm(id)%art(i,j,n)
                do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                  c = dm(id)%ivar(var,ii,j)+dm(id)%offset
                  asigma(l,c) = asigma(l,c)+temp*dmat(id)%derive(i,ii,der)
                enddo
              endif
            enddo
          enddo
        enddo

        do n=1,dm(id)%nartt
          power  = dm(id)%artti(1,n)
          der    = dm(id)%artti(2,n)
          eq     = dm(id)%artti(3,n)
          var    = dm(id)%artti(4,n)
#ifdef USE_COMPLEX
          if (dm(id)%artti(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          sigmap = cfactor*sigma**power
          do i=1,grd(id)%nr
            do j=1,nt
              l = dm(id)%ieq(eq,i,j)
              if (.not.dm(id)%bc_flag(l)) then
                l = l + dm(id)%offset
                do jj=1,nt
                  temp = sigmap*dm(id)%artt(i,j,jj,n)
                  do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                    c = dm(id)%ivar(var,ii,jj)+dm(id)%offset
                    asigma(l,c) = asigma(l,c)+temp*dmat(id)%derive(i,ii,der)
                  enddo
                enddo
              endif
            enddo
          enddo
        enddo
      enddo

      if (dump_asigma) then
          call dump_asigma_matrix(asigma, "asigma_full_total", ierr)
          if (ierr /= 0) return
      endif

      end subroutine make_asigma_full_total
#endif

!------------------------------------------------------------------------------
! This subroutine creates the matrix asigma using FULL storage, and compressed
! in the horizontal direction:
!      asigma(id,id2) = a(power=0,id,id2) + sigma.a(power=1,id,id2)
!                     + sigma^2.a(power=2,id,id2) +  ...
!
! Horizontal compression means that only non-zero columns of asigma are
! created.  Their indices are given by h_bc_range.
!------------------------------------------------------------------------------
! List of variables:
!
! sigma  = value for sigma
! asigma = matrix which will contain asigma
! id,id2 = domain coordinates for asigma
!------------------------------------------------------------------------------
! IMPORTANT: the case id == id2 is not allowed, because h_bc_range in this
!            case would not correspond to the non-zero columns.
!------------------------------------------------------------------------------

#ifdef USE_MULTI
      subroutine make_asigma_full_local_hcomp(sigma, asigma, id, id2, ierr)

      integer, intent(in)         :: id, id2
      integer, intent(out) :: ierr
#ifdef USE_COMPLEX
      double complex, intent(in)  :: sigma
      double complex, intent(out) :: asigma(dm(id)%d_dim,idm(id,id2)%n_h_bc)
      double complex sigmap, temp, cfactor
      double complex :: zero = (0d0, 0d0)
#else
      double precision, intent(in)  :: sigma
      double precision, intent(out) :: asigma(dm(id)%d_dim,idm(id,id2)%n_h_bc)
      double precision :: sigmap, temp, cfactor = 1d0
      double precision :: zero = 0d0
#endif
      integer n, i, j, ii, jj, l, c, power, der, eq, var, eqloc, varloc

      ! Just in case:
      if (id.eq.id2) then
        print*, "id == id2 not allowed in make_asigma_full_local_hcomp"
        ierr = 1
        return
      endif

      asigma = zero

      ! Boundary conditions:
      do n=1,idm(id,id2)%natbc
        power  = idm(id,id2)%atbci(1,n)
        der    = idm(id,id2)%atbci(2,n)
        eq     = idm(id,id2)%atbci(3,n)
        var    = idm(id,id2)%atbci(4,n)
        eqloc  = idm(id,id2)%atbci(5,n)
        varloc = idm(id,id2)%atbci(6,n)
#ifdef USE_COMPLEX
        if (idm(id,id2)%atbci(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do ii=max(1,varloc-dmat(id2)%lbder(der)),min(grd(id2)%nr,varloc+dmat(id2)%ubder(der))
          temp = sigmap*dmat(id2)%derive(varloc,ii,der)
          do j=1,nt
            l = dm(id)%ieq(eq,eqloc,j)
            c = idm(id,id2)%inv_h_bc_range(dm(id2)%ivar(var,ii,j))
            asigma(l,c) = asigma(l,c)+temp*idm(id,id2)%atbc(j,n)
          enddo
        enddo
      enddo

      do n=1,idm(id,id2)%nattbc
        power  = idm(id,id2)%attbci(1,n)
        der    = idm(id,id2)%attbci(2,n)
        eq     = idm(id,id2)%attbci(3,n)
        var    = idm(id,id2)%attbci(4,n)
        eqloc  = idm(id,id2)%attbci(5,n)
        varloc = idm(id,id2)%attbci(6,n)
#ifdef USE_COMPLEX
        if (idm(id,id2)%attbci(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do ii=max(1,varloc-dmat(id2)%lbder(der)),min(grd(id2)%nr,varloc+dmat(id2)%ubder(der))
          temp = sigmap*dmat(id2)%derive(varloc,ii,der)
          do j=1,nt
            l = dm(id)%ieq(eq,eqloc,j)
            do jj=1,nt
              c = idm(id,id2)%inv_h_bc_range(dm(id2)%ivar(var,ii,jj))
              asigma(l,c) = asigma(l,c)+temp*idm(id,id2)%attbc(j,jj,n)
            enddo
          enddo
        enddo
      enddo

      if (dump_asigma) then
          call dump_asigma_matrix(asigma, "asigma_full_local_hcomp", ierr)
          if (ierr /= 0) return
      endif

      end subroutine make_asigma_full_local_hcomp
#endif

!------------------------------------------------------------------------------
! This subroutine creates the matrix apower using FULL storage, where
! apower = a(0) when power = 0, apower = a(1) when power = 1, etc.
!------------------------------------------------------------------------------
! List of variables:
!
! power  = exponent of the eigenvalue which goes in front of the matrix in
!          the full eigenvalue problem
! apower = matrix into which the result is stored
! id,id2 = domain coordinates for apower
!------------------------------------------------------------------------------
#ifdef USE_MULTI
      subroutine make_a_full_local(power, apower, id, id2)

      integer, intent(in)         :: power, id, id2
#ifdef USE_COMPLEX
      double complex, intent(out) :: apower(dm(id)%d_dim,dm(id2)%d_dim)
      double complex cfactor
      double complex :: zero = (0d0, 0d0)
#else
      double precision, intent(out) :: apower(dm(id)%d_dim,dm(id2)%d_dim)
      double precision cfactor
      double precision :: zero = 0d0
#endif
      integer n, i, j, ii, jj, l, c, der, eq, var, eqloc, varloc

      apower = zero

      ! The boundary conditions:
      do n=1,idm(id,id2)%natbc
        if (power.eq.idm(id,id2)%atbci(1,n)) then
          der    = idm(id,id2)%atbci(2,n)
          eq     = idm(id,id2)%atbci(3,n)
          var    = idm(id,id2)%atbci(4,n)
          eqloc  = idm(id,id2)%atbci(5,n)
          varloc = idm(id,id2)%atbci(6,n)
#ifdef USE_COMPLEX
          if (idm(id,id2)%atbci(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          do j=1,nt
            l = dm(id)%ieq(eq,eqloc,j)
            do ii=max(1,varloc-dmat(id2)%lbder(der)),min(grd(id2)%nr,varloc+dmat(id2)%ubder(der))
              c = dm(id2)%ivar(var,ii,j)
              apower(l,c) = apower(l,c)+cfactor*idm(id,id2)%atbc(j,n)*dmat(id2)%derive(varloc,ii,der)
            enddo
          enddo
        endif
      enddo

      do n=1,idm(id,id2)%nattbc
        if (power.eq.idm(id,id2)%attbci(1,n)) then
          der    = idm(id,id2)%attbci(2,n)
          eq     = idm(id,id2)%attbci(3,n)
          var    = idm(id,id2)%attbci(4,n)
          eqloc  = idm(id,id2)%attbci(5,n)
          varloc = idm(id,id2)%attbci(6,n)
#ifdef USE_COMPLEX
          if (idm(id,id2)%attbci(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          do j=1,nt
            l = dm(id)%ieq(eq,eqloc,j)
            do ii=max(1,varloc-dmat(id2)%lbder(der)),min(grd(id2)%nr,varloc+dmat(id2)%ubder(der))
              do jj=1,nt
                c = dm(id2)%ivar(var,ii,jj)
                apower(l,c) = apower(l,c)+cfactor*idm(id,id2)%attbc(j,jj,n)*dmat(id2)%derive(varloc,ii,der)
              enddo
            enddo
          enddo
        endif
      enddo

      if (id.ne.id2) return

      ! The equations:
      do n=1,dm(id)%nas
        if (power.eq.dm(id)%asi(1,n)) then
          der = dm(id)%asi(2,n)
          eq  = dm(id)%asi(3,n)
          var = dm(id)%asi(4,n)
#ifdef USE_COMPLEX
          if (dm(id)%asi(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          do i=1,grd(id)%nr
            do j=1,nt
              l = dm(id)%ieq(eq,i,j)
              if (.not.dm(id)%bc_flag(l)) then
                do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                  c = dm(id)%ivar(var,ii,j)
                  apower(l,c) = apower(l,c)+cfactor*dm(id)%as(n)*dmat(id)%derive(i,ii,der)
                enddo
              endif
            enddo
          enddo
        endif
      enddo

      do n=1,dm(id)%nart
        if (power.eq.dm(id)%arti(1,n)) then
          der = dm(id)%arti(2,n)
          eq  = dm(id)%arti(3,n)
          var = dm(id)%arti(4,n)
#ifdef USE_COMPLEX
          if (dm(id)%arti(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          do i=1,grd(id)%nr
            do j=1,nt
              l = dm(id)%ieq(eq,i,j)
              if (.not.dm(id)%bc_flag(l)) then
                do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                  c = dm(id)%ivar(var,ii,j)
                  apower(l,c) = apower(l,c)+cfactor*dm(id)%art(i,j,n)*dmat(id)%derive(i,ii,der)
                enddo
              endif
            enddo
          enddo
        endif
      enddo

      do n=1,dm(id)%nartt
        if (power.eq.dm(id)%artti(1,n)) then
          der = dm(id)%artti(2,n)
          eq  = dm(id)%artti(3,n)
          var = dm(id)%artti(4,n)
#ifdef USE_COMPLEX
          if (dm(id)%artti(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          do i=1,grd(id)%nr
            do j=1,nt
              l = dm(id)%ieq(eq,i,j)
              if (.not.dm(id)%bc_flag(l)) then
                do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                  do jj=1,nt
                    c = dm(id)%ivar(var,ii,jj)
                    apower(l,c) = apower(l,c)+cfactor*dm(id)%artt(i,j,jj,n)*dmat(id)%derive(i,ii,der)
                  enddo
                enddo
              endif
            enddo
          enddo
        endif
      enddo

      end subroutine make_a_full_local
#endif

!------------------------------------------------------------------------------
! This subroutine creates the TOTAL matrix apower using FULL storage, where
! apower = a(0) when power = 0, apower = a(1) when power = 1, etc.
!------------------------------------------------------------------------------
! List of variables:
!
! power  = exponent of the eigenvalue which goes in front of the matrix in
!          the full eigenvalue problem
! apower = matrix into which the result is stored
!------------------------------------------------------------------------------
#ifndef USE_1D
      subroutine make_a_full_total(power, apower)

      integer, intent(in)         :: power
#ifdef USE_COMPLEX
      double complex, intent(out) :: apower(a_dim,a_dim)
      double complex cfactor
      double complex :: zero = (0d0, 0d0)
#else
      double precision, intent(out) :: apower(a_dim,a_dim)
      double precision cfactor
      double precision :: zero = 0d0
#endif
      integer n, i, j, ii, jj, l, c, der, eq, var, eqloc, varloc, id, id2

      apower = zero

      ! The boundary conditions:
      do id = 1, ndomains
        do id2 = 1, ndomains
          do n=1,idm(id,id2)%natbc
            if (power.eq.idm(id,id2)%atbci(1,n)) then
              der    = idm(id,id2)%atbci(2,n)
              eq     = idm(id,id2)%atbci(3,n)
              var    = idm(id,id2)%atbci(4,n)
              eqloc  = idm(id,id2)%atbci(5,n)
              varloc = idm(id,id2)%atbci(6,n)
#ifdef USE_COMPLEX
              if (idm(id,id2)%atbci(7,n).eq.0) then
                cfactor = (1d0,0d0)
              else
                cfactor = (0d0,1d0)
              endif
#endif
              do j=1,nt
                l = dm(id)%ieq(eq,eqloc,j)+dm(id)%offset
                do ii=max(1,varloc-dmat(id2)%lbder(der)),min(grd(id2)%nr,varloc+dmat(id2)%ubder(der))
                  c = dm(id2)%ivar(var,ii,j)+dm(id2)%offset
                  apower(l,c) = apower(l,c)+cfactor*idm(id,id2)%atbc(j,n)*dmat(id2)%derive(varloc,ii,der)
                enddo
              enddo
            endif
          enddo

          do n=1,idm(id,id2)%nattbc
            if (power.eq.idm(id,id2)%attbci(1,n)) then
              der    = idm(id,id2)%attbci(2,n)
              eq     = idm(id,id2)%attbci(3,n)
              var    = idm(id,id2)%attbci(4,n)
              eqloc  = idm(id,id2)%attbci(5,n)
              varloc = idm(id,id2)%attbci(6,n)
#ifdef USE_COMPLEX
              if (idm(id,id2)%attbci(7,n).eq.0) then
                cfactor = (1d0,0d0)
              else
                cfactor = (0d0,1d0)
              endif
#endif
              do j=1,nt
                l = dm(id)%ieq(eq,eqloc,j)+dm(id)%offset
                do ii=max(1,varloc-dmat(id2)%lbder(der)),min(grd(id2)%nr,varloc+dmat(id2)%ubder(der))
                  do jj=1,nt
                    c = dm(id2)%ivar(var,ii,jj)+dm(id2)%offset
                    apower(l,c) = apower(l,c)+cfactor*idm(id,id2)%attbc(j,jj,n)*dmat(id2)%derive(varloc,ii,der)
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
      enddo

      ! The equations:
      do id=1,ndomains
        do n=1,dm(id)%nas
          if (power.eq.dm(id)%asi(1,n)) then
            der = dm(id)%asi(2,n)
            eq  = dm(id)%asi(3,n)
            var = dm(id)%asi(4,n)
#ifdef USE_COMPLEX
            if (dm(id)%asi(7,n).eq.0) then
              cfactor = (1d0,0d0)
            else
              cfactor = (0d0,1d0)
            endif
#endif
            do i=1,grd(id)%nr
              do j=1,nt
                l = dm(id)%ieq(eq,i,j)
                if (.not.dm(id)%bc_flag(l)) then
                  l = l + dm(id)%offset
                  do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                    c = dm(id)%ivar(var,ii,j)+dm(id)%offset
                    apower(l,c) = apower(l,c)+cfactor*dm(id)%as(n)*dmat(id)%derive(i,ii,der)
                  enddo
                endif
              enddo
            enddo
          endif
        enddo

        do n=1,dm(id)%nart
          if (power.eq.dm(id)%arti(1,n)) then
            der = dm(id)%arti(2,n)
            eq  = dm(id)%arti(3,n)
            var = dm(id)%arti(4,n)
#ifdef USE_COMPLEX
            if (dm(id)%arti(7,n).eq.0) then
              cfactor = (1d0,0d0)
            else
              cfactor = (0d0,1d0)
            endif
#endif
            do i=1,grd(id)%nr
              do j=1,nt
                l = dm(id)%ieq(eq,i,j)
                if (.not.dm(id)%bc_flag(l)) then
                  l = l + dm(id)%offset
                  do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                    c = dm(id)%ivar(var,ii,j)+dm(id)%offset
                    apower(l,c) = apower(l,c)+cfactor*dm(id)%art(i,j,n)*dmat(id)%derive(i,ii,der)
                  enddo
                endif
              enddo
            enddo
          endif
        enddo

        do n=1,dm(id)%nartt
          if (power.eq.dm(id)%artti(1,n)) then
            der = dm(id)%artti(2,n)
            eq  = dm(id)%artti(3,n)
            var = dm(id)%artti(4,n)
#ifdef USE_COMPLEX
            if (dm(id)%artti(7,n).eq.0) then
              cfactor = (1d0,0d0)
            else
              cfactor = (0d0,1d0)
            endif
#endif
            do i=1,grd(id)%nr
              do j=1,nt
                l = dm(id)%ieq(eq,i,j)
                if (.not.dm(id)%bc_flag(l)) then
                  l = l + dm(id)%offset
                  do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                    do jj=1,nt
                      c = dm(id)%ivar(var,ii,jj)+dm(id)%offset
                      apower(l,c) = apower(l,c)+cfactor*dm(id)%artt(i,j,jj,n)*dmat(id)%derive(i,ii,der)
                    enddo
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

      end subroutine make_a_full_total
#endif

!------------------------------------------------------------------------------
! This subroutine creates the matrix asigma using BAND storage:
!      asigma(id,id) = a(0) + sigma.a(1) + sigma^2.a(2) +  ...
!------------------------------------------------------------------------------
! List of variables:
!
! sigma  = value for sigma
! asigma = matrix which will contain asigma
! id     = domain number
!------------------------------------------------------------------------------
! NOTE: a second domain argument, id2, is not given because band storage
!       only seems appropriate for square matrices, which usually only occurs
!       when id == id2 (and the corresponding ku, kl are only defined for
!       id == id2)
!------------------------------------------------------------------------------
#ifdef USE_MULTI
      subroutine make_asigma_band_local(sigma, asigma, id, ierr)

      integer, intent(in)         :: id
#ifdef USE_COMPLEX
      double complex, intent(in)  :: sigma
      double complex, intent(out) :: asigma(2*dm(id)%kl+dm(id)%ku+1,dm(id)%d_dim)
      double complex sigmap, temp, cfactor
#else
      double precision, intent(in)  :: sigma
      double precision, intent(out) :: asigma(2*dm(id)%kl+dm(id)%ku+1,dm(id)%d_dim)
      double precision sigmap, temp, cfactor
#endif
      integer n, i, j, ii, jj, l, c, power, der, eq, var, eqloc, varloc, pos
      integer, intent(out) :: ierr

      asigma = 0d0
      pos = dm(id)%kl + dm(id)%ku + 1

      ! The equations:
      do n=1,dm(id)%nas
        power  = dm(id)%asi(1,n)
        der    = dm(id)%asi(2,n)
        eq     = dm(id)%asi(3,n)
        var    = dm(id)%asi(4,n)
#ifdef USE_COMPLEX
        if (dm(id)%asi(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        temp   = sigmap*dm(id)%as(n)
        do i=1,grd(id)%nr
          do j=1,nt
            l = dm(id)%ieq(eq,i,j)
            if (.not.dm(id)%bc_flag(l)) then
              do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                c = dm(id)%ivar(var,ii,j)
                asigma(pos+l-c,c) = asigma(pos+l-c,c)+temp*dmat(id)%derive(i,ii,der)
              enddo
            endif
          enddo
        enddo
      enddo

      do n=1,dm(id)%nart
        power  = dm(id)%arti(1,n)
        der    = dm(id)%arti(2,n)
        eq     = dm(id)%arti(3,n)
        var    = dm(id)%arti(4,n)
#ifdef USE_COMPLEX
        if (dm(id)%arti(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do i=1,grd(id)%nr
          do j=1,nt
            l = dm(id)%ieq(eq,i,j)
            if (.not.dm(id)%bc_flag(l)) then
              temp = sigmap*dm(id)%art(i,j,n)
              do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                c = dm(id)%ivar(var,ii,j)
                asigma(pos+l-c,c) = asigma(pos+l-c,c)+temp*dmat(id)%derive(i,ii,der)
              enddo
            endif
          enddo
        enddo
      enddo

      do n=1,dm(id)%nartt
        power  = dm(id)%artti(1,n)
        der    = dm(id)%artti(2,n)
        eq     = dm(id)%artti(3,n)
        var    = dm(id)%artti(4,n)
#ifdef USE_COMPLEX
        if (dm(id)%artti(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do i=1,grd(id)%nr
          do j=1,nt
            l = dm(id)%ieq(eq,i,j)
            if (.not.dm(id)%bc_flag(l)) then
              do jj=1,nt
                temp = sigmap*dm(id)%artt(i,j,jj,n)
                do ii=max(1,i-dmat(id)%lbder(der)),min(grd(id)%nr,i+dmat(id)%ubder(der))
                  c = dm(id)%ivar(var,ii,jj)
                  asigma(pos+l-c,c) = asigma(pos+l-c,c)+temp*dmat(id)%derive(i,ii,der)
                enddo
              enddo
            endif
          enddo
        enddo
      enddo

      ! Boundary conditions:
      do n=1,idm(id,id)%natbc
        power  = idm(id,id)%atbci(1,n)
        der    = idm(id,id)%atbci(2,n)
        eq     = idm(id,id)%atbci(3,n)
        var    = idm(id,id)%atbci(4,n)
        eqloc  = idm(id,id)%atbci(5,n)
        varloc = idm(id,id)%atbci(6,n)
#ifdef USE_COMPLEX
        if (idm(id,id)%atbci(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do ii=max(1,varloc-dmat(id)%lbder(der)),min(grd(id)%nr,varloc+dmat(id)%ubder(der))
          temp = sigmap*dmat(id)%derive(varloc,ii,der)
          do j=1,nt
            l = dm(id)%ieq(eq,eqloc,j)
            c = dm(id)%ivar(var,ii,j)
            asigma(pos+l-c,c) = asigma(pos+l-c,c)+temp*idm(id,id)%atbc(j,n)
          enddo
        enddo
      enddo

      do n=1,idm(id,id)%nattbc
        power  = idm(id,id)%attbci(1,n)
        der    = idm(id,id)%attbci(2,n)
        eq     = idm(id,id)%attbci(3,n)
        var    = idm(id,id)%attbci(4,n)
        eqloc  = idm(id,id)%attbci(5,n)
        varloc = idm(id,id)%attbci(6,n)
#ifdef USE_COMPLEX
        if (idm(id,id)%attbci(7,n).eq.0) then
          cfactor = (1d0,0d0)
        else
          cfactor = (0d0,1d0)
        endif
#endif
        sigmap = cfactor*sigma**power
        do ii=max(1,varloc-dmat(id)%lbder(der)),min(grd(id)%nr,varloc+dmat(id)%ubder(der))
          temp = sigmap*dmat(id)%derive(varloc,ii,der)
          do j=1,nt
            l = dm(id)%ieq(eq,eqloc,j)
            do jj=1,nt
              c = dm(id)%ivar(var,ii,jj)
              asigma(pos+l-c,c) = asigma(pos+l-c,c)+temp*idm(id,id)%attbc(j,jj,n)
            enddo
          enddo
        enddo
      enddo

      if (dump_asigma) then
          call dump_asigma_matrix(asigma, "asigma_band_local", ierr)
          if (ierr /= 0) return
      endif

      end subroutine make_asigma_band_local
#endif

!------------------------------------------------------------------------------
! This subroutine creates the matrix apower using BAND storage, where
! apower = a(0) when power = 0, apower = a(1) when power = 1, etc.
!------------------------------------------------------------------------------
! List of variables:
!
! power  = exponent of the eigenvalue which goes in front of the matrix in
!          the full eigenvalue problem
! apower = matrix into which the result is stored
! id,id2 = domain coordinates for apower
!------------------------------------------------------------------------------
! NOTE: a second domain argument, id2, is not given because band storage
!       only seems appropriate for square matrices, which usually only occurs
!       when id == id2 (and the corresponding ku, kl are only defined for
!       id == id2)
!------------------------------------------------------------------------------
#ifdef USE_MULTI
      subroutine make_a_band_local(power, apower, id)

      integer, intent(in)         :: power, id
#ifdef USE_COMPLEX
      double complex, intent(out) :: apower(2*dm(id)%kl+dm(id)%ku+1,dm(id)%d_dim)
      double complex cfactor
      double complex :: zero = (0d0, 0d0)
#else
      double precision, intent(out) :: apower(2*dm(id)%kl+dm(id)%ku+1,dm(id)%d_dim)
      double precision cfactor
      double precision :: zero = 0d0
#endif
      integer n, i, j, ii, jj, l, c, der, eq, var, eqloc, varloc, pos

      apower = zero
      pos = dm(id)%kl + dm(id)%ku + 1

      ! The equations:
      do n=1,dm(id)%nas
        if (power.eq.dm(id)%asi(1,n)) then
          der = dm(id)%asi(2,n)
          eq  = dm(id)%asi(3,n)
          var = dm(id)%asi(4,n)
#ifdef USE_COMPLEX
          if (dm(id)%asi(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          do i=1,grd(id)%nr
            do j=1,nt
              l = dm(id)%ieq(eq,i,j)
              if (.not.dm(id)%bc_flag(l)) then
                do ii=max(1,i-dmat(id)%lbder(der)), &
                      min(grd(id)%nr,i+dmat(id)%ubder(der))
                  c = dm(id)%ivar(var,ii,j)
                  apower(pos+l-c,c) = apower(pos+l-c,c) + &
                    cfactor*dm(id)%as(n)*dmat(id)%derive(i,ii,der)
                enddo
              endif
            enddo
          enddo
        endif
      enddo

      do n=1,dm(id)%nart
        if (power.eq.dm(id)%arti(1,n)) then
          der = dm(id)%arti(2,n)
          eq  = dm(id)%arti(3,n)
          var = dm(id)%arti(4,n)
#ifdef USE_COMPLEX
          if (dm(id)%arti(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          do i=1,grd(id)%nr
            do j=1,nt
              l = dm(id)%ieq(eq,i,j)
              if (.not.dm(id)%bc_flag(l)) then
                do ii=max(1,i-dmat(id)%lbder(der)), &
                      min(grd(id)%nr,i+dmat(id)%ubder(der))
                  c = dm(id)%ivar(var,ii,j)
                  apower(pos+l-c,c) = apower(pos+l-c,c) + &
                    cfactor*dm(id)%art(i,j,n)*dmat(id)%derive(i,ii,der)
                enddo
              endif
            enddo
          enddo
        endif
      enddo

      do n=1,dm(id)%nartt
        if (power.eq.dm(id)%artti(1,n)) then
          der = dm(id)%artti(2,n)
          eq  = dm(id)%artti(3,n)
          var = dm(id)%artti(4,n)
#ifdef USE_COMPLEX
          if (dm(id)%artti(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          do i=1,grd(id)%nr
            do j=1,nt
              l = dm(id)%ieq(eq,i,j)
              if (.not.dm(id)%bc_flag(l)) then
                do ii=max(1,i-dmat(id)%lbder(der)), &
                      min(grd(id)%nr,i+dmat(id)%ubder(der))
                  do jj=1,nt
                    c = dm(id)%ivar(var,ii,jj)
                    apower(pos+l-c,c) = apower(pos+l-c,c) + &
                      cfactor*dm(id)%artt(i,j,jj,n)*dmat(id)%derive(i,ii,der)
                  enddo
                enddo
              endif
            enddo
          enddo
        endif
      enddo

      ! The boundary conditions:
      do n=1,idm(id,id)%natbc
        if (power.eq.idm(id,id)%atbci(1,n)) then
          der    = idm(id,id)%atbci(2,n)
          eq     = idm(id,id)%atbci(3,n)
          var    = idm(id,id)%atbci(4,n)
          eqloc  = idm(id,id)%atbci(5,n)
          varloc = idm(id,id)%atbci(6,n)
#ifdef USE_COMPLEX
          if (idm(id,id)%atbci(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          do j=1,nt
            l = dm(id)%ieq(eq,eqloc,j)
            do ii=max(1,varloc-dmat(id)%lbder(der)), &
                  min(grd(id)%nr,varloc+dmat(id)%ubder(der))
              c = dm(id)%ivar(var,ii,j)
              apower(pos+l-c,c) = apower(pos+l-c,c) + &
                cfactor*idm(id,id)%atbc(j,n)*dmat(id)%derive(varloc,ii,der)
            enddo
          enddo
        endif
      enddo

      do n=1,idm(id,id)%nattbc
        if (power.eq.idm(id,id)%attbci(1,n)) then
          der    = idm(id,id)%attbci(2,n)
          eq     = idm(id,id)%attbci(3,n)
          var    = idm(id,id)%attbci(4,n)
          eqloc  = idm(id,id)%attbci(5,n)
          varloc = idm(id,id)%attbci(6,n)
#ifdef USE_COMPLEX
          if (idm(id,id)%attbci(7,n).eq.0) then
            cfactor = (1d0,0d0)
          else
            cfactor = (0d0,1d0)
          endif
#endif
          do j=1,nt
            l = dm(id)%ieq(eq,eqloc,j)
            do ii=max(1,varloc-dmat(id)%lbder(der)), &
                  min(grd(id)%nr,varloc+dmat(id)%ubder(der))
              do jj=1,nt
                c = dm(id)%ivar(var,ii,jj)
                apower(pos+l-c,c) = apower(pos+l-c,c) + &
                  cfactor*idm(id,id)%attbc(j,jj,n)*dmat(id)%derive(varloc,ii,der)
              enddo
            enddo
          enddo
        endif
      enddo

      end subroutine make_a_band_local
#endif
!------------------------------------------------------------------------------
      subroutine a_product_transpose(vect_in,vect_out,power)

#ifdef USE_COMPLEX
      double complex, intent(in) :: vect_in(a_dim)
      double complex, intent(out):: vect_out(a_dim)
#else
      double precision, intent(in) :: vect_in(a_dim)
      double precision, intent(out):: vect_out(a_dim)
#endif
      integer, intent(in) :: power
      integer i,j,ii,jj,n,der,eq,var,eqloc,varloc,a_index,v_index

      ! Initialisation:
#ifdef USE_1D
            ! Initialisation:
      dm(1)%vect_der(1:a_dim,dmat(1)%der_min:dmat(1)%der_max) = 0d0

      ! The equations:
      do n=1, dm(1)%nas
        if (dm(1)%asi(1,n).eq.power) then
          der = dm(1)%asi(2,n)
          eq  = dm(1)%asi(3,n)
          var = dm(1)%asi(4,n)
          do i=1,nr
            a_index = dm(1)%ieq(eq,i)
            if (.not.dm(1)%bc_flag(a_index)) then
              v_index = dm(1)%ivar(var,i)
              dm(1)%vect_der(v_index,der) = dm(1)%vect_der(v_index,der) &
                    & + dm(1)%as(n)*vect_in(a_index)
            endif
          enddo
        endif
      enddo

      do n=1, dm(1)%nar
        if (dm(1)%ari(1,n).eq.power) then
          der = dm(1)%ari(2,n)
          eq  = dm(1)%ari(3,n)
          var = dm(1)%ari(4,n)
          do i=1,nr
            a_index = dm(1)%ieq(eq,i)
            if (.not.dm(1)%bc_flag(a_index)) then
              v_index = dm(1)%ivar(var,i)
              dm(1)%vect_der(v_index,der) = dm(1)%vect_der(v_index,der) &
                    & + dm(1)%ar(i,n)*vect_in(a_index)
            endif
          enddo
        endif
      enddo

      ! Boundary conditions:
      do n=1, dm(1)%nasbc
        if (dm(1)%asbci(1,n).eq.power) then
          der = dm(1)%asbci(2,n)
          eq  = dm(1)%asbci(3,n)
          var = dm(1)%asbci(4,n)
          eqloc  = dm(1)%asbci(5,n)
          varloc = dm(1)%asbci(6,n)
          a_index = dm(1)%ieq(eq,eqloc)
          v_index = dm(1)%ivar(var,varloc)
          dm(1)%vect_der(v_index,der) = dm(1)%vect_der(v_index,der) + &
            & dm(1)%asbc(n)*vect_in(a_index)
        endif
      enddo

      ! Finish the calculations: transpose(mat_der)*vect
      vect_out(1:a_dim) = 0d0
      do var=1, dm(1)%nvar
        do der=dm(1)%var_der_min(var),dm(1)%var_der_max(var)
          do i=1,nr
            a_index = dm(1)%ivar(var,i)
            do ii=max(1,i-dmat(1)%ubder(der)),min(nr,i+dmat(1)%lbder(der))
              v_index = dm(1)%ivar(var,ii)

              ! To understand this, you need to draw a picture...
              vect_out(a_index) = vect_out(a_index) + &
                      dmat(1)%derive(ii,i,der)*dm(1)%vect_der(v_index,der)
            enddo
          enddo
        enddo
      enddo
#else
      dm(1)%vect_der(1:a_dim,dmat(1)%der_min:dmat(1)%der_max) = 0d0

      ! The equations:
      do n=1,dm(1)%nas
        if (dm(1)%asi(1,n).eq.power) then
          der = dm(1)%asi(2,n)
          eq  = dm(1)%asi(3,n)
          var = dm(1)%asi(4,n)
          do i=1,grd(1)%nr
            do j=1,nt
              a_index = dm(1)%ieq(eq,i,j)
              if (.not.dm(1)%bc_flag(a_index)) then
                v_index = dm(1)%ivar(var,i,j)
                dm(1)%vect_der(v_index,der) = dm(1)%vect_der(v_index,der) + dm(1)%as(n)*vect_in(a_index)
              endif
            enddo
          enddo
        endif
      enddo

      do n=1,dm(1)%nart
        if (dm(1)%arti(1,n).eq.power) then
          der = dm(1)%arti(2,n)
          eq  = dm(1)%arti(3,n)
          var = dm(1)%arti(4,n)
          do i=1,grd(1)%nr
            do j=1,nt
              a_index = dm(1)%ieq(eq,i,j)
              if (.not.dm(1)%bc_flag(a_index)) then
                v_index = dm(1)%ivar(var,i,j)
                dm(1)%vect_der(v_index,der) = dm(1)%vect_der(v_index,der) + dm(1)%art(i,j,n)*vect_in(a_index)
              endif
            enddo
          enddo
        endif
      enddo

      do n=1,dm(1)%nartt
        if (dm(1)%artti(1,n).eq.power) then
          der = dm(1)%artti(2,n)
          eq  = dm(1)%artti(3,n)
          var = dm(1)%artti(4,n)
          do i=1,grd(1)%nr
            do j=1,nt
              a_index = dm(1)%ieq(eq,i,j)
              if (.not.dm(1)%bc_flag(a_index)) then
                do jj=1,nt
                  v_index = dm(1)%ivar(var,i,jj)
                  dm(1)%vect_der(v_index,der) = dm(1)%vect_der(v_index,der) + dm(1)%artt(i,j,jj,n)*vect_in(a_index)
                enddo
              endif
            enddo
          enddo
        endif
      enddo

      ! Boundary conditions:
      do n=1,idm(1, 1)%natbc
        if (idm(1, 1)%atbci(1,n).eq.power) then
          der = idm(1, 1)%atbci(2,n)
          eq  = idm(1, 1)%atbci(3,n)
          var = idm(1, 1)%atbci(4,n)
          eqloc  = idm(1, 1)%atbci(5,n)
          varloc = idm(1, 1)%atbci(6,n)
          do j=1,nt
            a_index = dm(1)%ieq(eq,eqloc,j)
            v_index = dm(1)%ivar(var,varloc,j)
            dm(1)%vect_der(v_index,der) = dm(1)%vect_der(v_index,der) + idm(1, 1)%atbc(j,n)*vect_in(a_index)
          enddo
        endif
      enddo

      do n=1,idm(1, 1)%nattbc
        if (idm(1, 1)%attbci(1,n).eq.power) then
          der = idm(1, 1)%attbci(2,n)
          eq  = idm(1, 1)%attbci(3,n)
          var = idm(1, 1)%attbci(4,n)
          eqloc  = idm(1, 1)%attbci(5,n)
          varloc = idm(1, 1)%attbci(6,n)
          do j=1,nt
            a_index = dm(1)%ieq(eq,eqloc,j)
            do jj=1,nt
              v_index = dm(1)%ivar(var,varloc,jj)
              dm(1)%vect_der(v_index,der) = dm(1)%vect_der(v_index,der) + idm(1, 1)%attbc(j,jj,n)*vect_in(a_index)
            enddo
          enddo
        endif
      enddo

      ! Finish the calculations: transpose(mat_der)*vect
      vect_out = 0d0 ! don't forget this
      do var=1,dm(1)%nvar
        do der=dm(1)%var_der_min(var),dm(1)%var_der_max(var)
          do j=1,nt
            do i=1,grd(1)%nr
              v_index = dm(1)%ivar(var,i,j)
              do ii=max(1,i-dmat(1)%ubder(der)),min(grd(1)%nr,i+dmat(1)%lbder(der))
                a_index = dm(1)%ivar(var,ii,j)
                ! To understand this, you need to draw a picture...
                vect_out(v_index) = vect_out(v_index) + &
                        dmat(1)%derive(ii,i,der)*dm(1)%vect_der(a_index,der)
              enddo
            enddo
          enddo
        enddo
      enddo

#endif
      end subroutine a_product_transpose
!------------------------------------------------------------------------------
      subroutine make_asigma_band(sigma, asigma, ierr)

#ifdef USE_COMPLEX
      double complex sigma
      double complex asigma(2*dm(1)%kl+dm(1)%ku+1,a_dim)
      double complex sigmap, temp
#else
      double precision sigma
      double precision asigma(2*dm(1)%kl+dm(1)%ku+1,a_dim)
      double precision sigmap, temp
#endif
      integer n, i, j, ii, jj, l, c, power, der, eq, var, eqloc, varloc
      integer pos
      integer, intent(out) :: ierr

      asigma = 0d0

#ifdef USE_1D

      ! The equations:
      do n=1, dm(1)%nas
        power = dm(1)%asi(1,n)
        der = dm(1)%asi(2,n)
        eq  = dm(1)%asi(3,n)
        var = dm(1)%asi(4,n)
        sigmap = sigma**power
        do i=1,nr
          l = dm(1)%ieq(eq,i)
          if (.not.dm(1)%bc_flag(l)) then
            do ii=max(1,i-dmat(1)%lbder(der)),min(nr,i+dmat(1)%ubder(der))
              c = dm(1)%ivar(var,ii)
              asigma(l-c+dm(1)%kl+dm(1)%ku+1,c) = &
                & asigma(l-c+dm(1)%kl+dm(1)%ku+1,c)+ &
                & sigmap*dm(1)%as(n)*dmat(1)%derive(i,ii,der)
            enddo
          endif
        enddo
      enddo

      do n=1, dm(1)%nar
        power = dm(1)%ari(1,n)
        der = dm(1)%ari(2,n)
        eq  = dm(1)%ari(3,n)
        var = dm(1)%ari(4,n)
        sigmap = sigma**power
        do i=1,nr
          l = dm(1)%ieq(eq,i)
          if (.not.dm(1)%bc_flag(l)) then
            do ii=max(1,i-dmat(1)%lbder(der)),min(nr,i+dmat(1)%ubder(der))
              c = dm(1)%ivar(var,ii)
              asigma(l-c+dm(1)%kl+dm(1)%ku+1,c) = &
                & asigma(l-c+dm(1)%kl+dm(1)%ku+1,c)+ &
                & sigmap*dm(1)%ar(i,n)*dmat(1)%derive(i,ii,der)
            enddo
          endif
        enddo
      enddo

      ! Boundary conditions:
      do n=1, dm(1)%nasbc
        power = dm(1)%asbci(1,n)
        der = dm(1)%asbci(2,n)
        eq  = dm(1)%asbci(3,n)
        var = dm(1)%asbci(4,n)
        eqloc  = dm(1)%asbci(5,n)
        varloc = dm(1)%asbci(6,n)
        sigmap = sigma**power
        l = dm(1)%ieq(eq,eqloc)
        do ii=max(1,varloc-dmat(1)%lbder(der)),min(nr,varloc+dmat(1)%ubder(der))
          c = dm(1)%ivar(var,ii)
          asigma(l-c+dm(1)%kl+dm(1)%ku+1,c) = &
            & asigma(l-c+dm(1)%kl+dm(1)%ku+1,c)+ &
            & sigmap*dm(1)%asbc(n)*dmat(1)%derive(varloc,ii,der)
        enddo
      enddo
#else

      pos = dm(1)%kl + dm(1)%ku + 1

      ! The equations:
      do n=1,dm(1)%nas
        power = dm(1)%asi(1,n)
        der = dm(1)%asi(2,n)
        eq  = dm(1)%asi(3,n)
        var = dm(1)%asi(4,n)
        sigmap = sigma**power
        temp = sigmap*dm(1)%as(n)
        do i=1,grd(1)%nr
          do j=1,nt
            l = dm(1)%ieq(eq,i,j)
            if (.not.dm(1)%bc_flag(l)) then
              do ii=max(1,i-dmat(1)%lbder(der)),min(grd(1)%nr,i+dmat(1)%ubder(der))
                c = dm(1)%ivar(var,ii,j)
                asigma(l-c+pos,c) = asigma(l-c+pos,c)+temp*dmat(1)%derive(i,ii,der)
              enddo
            endif
          enddo
        enddo
      enddo

      do n=1,dm(1)%nart
        power = dm(1)%arti(1,n)
        der = dm(1)%arti(2,n)
        eq  = dm(1)%arti(3,n)
        var = dm(1)%arti(4,n)
        sigmap = sigma**power
        do j=1,nt
          do i=1,grd(1)%nr
            l = dm(1)%ieq(eq,i,j)
            if (.not.dm(1)%bc_flag(l)) then
              temp = sigmap*dm(1)%art(i,j,n)
              do ii=max(1,i-dmat(1)%lbder(der)),min(grd(1)%nr,i+dmat(1)%ubder(der))
                c = dm(1)%ivar(var,ii,j)
                asigma(l-c+pos,c) = asigma(l-c+pos,c)+temp*dmat(1)%derive(i,ii,der)
              enddo
            endif
          enddo
        enddo
      enddo

      do n=1,dm(1)%nartt
        power = dm(1)%artti(1,n)
        der = dm(1)%artti(2,n)
        eq  = dm(1)%artti(3,n)
        var = dm(1)%artti(4,n)
        sigmap = sigma**power
        do i=1,grd(1)%nr
          do j=1,nt
            l = dm(1)%ieq(eq,i,j)
            if (.not.dm(1)%bc_flag(l)) then
              do jj=1,nt
                temp = sigmap*dm(1)%artt(i,j,jj,n)
                do ii=max(1,i-dmat(1)%lbder(der)),min(grd(1)%nr,i+dmat(1)%ubder(der))
                  c = dm(1)%ivar(var,ii,jj)
                  asigma(l-c+pos,c) = asigma(l-c+pos,c)+temp*dmat(1)%derive(i,ii,der)
                enddo
              enddo
            endif
          enddo
        enddo
      enddo

      ! Boundary conditions:
      do n=1,idm(1, 1)%natbc
        power = idm(1, 1)%atbci(1,n)
        der = idm(1, 1)%atbci(2,n)
        eq  = idm(1, 1)%atbci(3,n)
        var = idm(1, 1)%atbci(4,n)
        eqloc  = idm(1, 1)%atbci(5,n)
        varloc = idm(1, 1)%atbci(6,n)
        sigmap = sigma**power
        do ii=max(1,varloc-dmat(1)%lbder(der)),min(grd(1)%nr,varloc+dmat(1)%ubder(der))
          temp = sigmap*dmat(1)%derive(varloc,ii,der)
          do j=1,nt
            l = dm(1)%ieq(eq,eqloc,j)
            c = dm(1)%ivar(var,ii,j)
            asigma(l-c+pos,c) = asigma(l-c+pos,c)+temp*idm(1, 1)%atbc(j,n)
          enddo
        enddo
      enddo

      do n=1,idm(1, 1)%nattbc
        power = idm(1, 1)%attbci(1,n)
        der = idm(1, 1)%attbci(2,n)
        eq  = idm(1, 1)%attbci(3,n)
        var = idm(1, 1)%attbci(4,n)
        eqloc  = idm(1, 1)%attbci(5,n)
        varloc = idm(1, 1)%attbci(6,n)
        sigmap = sigma**power
        do ii=max(1,varloc-dmat(1)%lbder(der)),min(grd(1)%nr,varloc+dmat(1)%ubder(der))
          temp = sigmap*dmat(1)%derive(varloc,ii,der)
          do j=1,nt
            l = dm(1)%ieq(eq,eqloc,j)
            do jj=1,nt
              c = dm(1)%ivar(var,ii,jj)
              asigma(l-c+pos,c) = asigma(l-c+pos,c)+temp*idm(1, 1)%attbc(j,jj,n)
            enddo
          enddo
        enddo
      enddo

#endif

      if (dump_asigma) then
          call dump_asigma_matrix(asigma, "asigma_band", ierr)
          if (ierr /= 0) return
      endif

      end subroutine make_asigma_band
!------------------------------------------------------------------------------
#ifdef USE_1D
      subroutine a_product(vect_in, vect_out, power)

      double precision, intent(in) :: vect_in(a_dim)
      double precision, intent(out) :: vect_out(a_dim)
      integer, intent(in) :: power
      integer i, ii, n, der, eq, var, a_index, v_index
      integer eqloc, varloc

      ! Preliminary calculations: derivative(s) of vect
      dm(1)%vect_der(1:a_dim, dmat(1)%der_min:dmat(1)%der_max) = 0d0
      do var=1, dm(1)%nvar
        do der= dm(1)%var_der_min(var),  dm(1)%var_der_max(var)
          do i=1, nr
            a_index = dm(1)%ivar(var, i)
            do ii=max(1,i-dmat(1)%lbder(der)),min(nr,i+dmat(1)%ubder(der))
              v_index = dm(1)%ivar(var,ii)
              dm(1)%vect_der(a_index,der) = dm(1)%vect_der(a_index,der) + &
                               dmat(1)%derive(i,ii,der)*vect_in(v_index)
            enddo
          enddo
        enddo
      enddo

      vect_out(1:a_dim) = 0d0

      ! The equations:
      do n=1,  dm(1)%nas
        if (dm(1)%asi(1,n).eq.power) then
          der = dm(1)%asi(2,n)
          eq  = dm(1)%asi(3,n)
          var = dm(1)%asi(4,n)
          do i=1,nr
            v_index = dm(1)%ieq(eq,i)
            if (.not.dm(1)%bc_flag(v_index)) then
              a_index = dm(1)%ivar(var,i)
              vect_out(v_index) = vect_out(v_index) + &
                & dm(1)%as(n)*dm(1)%vect_der(a_index,der)
            endif
          enddo
        endif
      enddo

        do n=1, dm(1)%nar
        if (dm(1)%ari(1,n).eq.power) then
          der = dm(1)%ari(2,n)
          eq  = dm(1)%ari(3,n)
          var = dm(1)%ari(4,n)
          do i=1,nr
            v_index = dm(1)%ieq(eq,i)
            if (.not.dm(1)%bc_flag(v_index)) then
              a_index = dm(1)%ivar(var,i)
              vect_out(v_index) = vect_out(v_index) + &
                & dm(1)%ar(i,n)*dm(1)%vect_der(a_index,der)
            endif
          enddo
        endif
      enddo

      ! Boundary conditions:
      do n=1, dm(1)%nasbc
        if (dm(1)%asbci(1,n).eq.power) then
          der = dm(1)%asbci(2,n)
          eq  = dm(1)%asbci(3,n)
          var = dm(1)%asbci(4,n)
          eqloc  = dm(1)%asbci(5,n)
          varloc = dm(1)%asbci(6,n)
          a_index = dm(1)%ivar(var,varloc)
          v_index = dm(1)%ieq(eq,eqloc)
          vect_out(v_index) = vect_out(v_index) + &
            & dm(1)%asbc(n)*dm(1)%vect_der(a_index,der)
        endif
      enddo

      end subroutine a_product
#else
      subroutine a_product(vect_in, vect_out, power)

#ifdef USE_COMPLEX
      double complex, intent(in)  :: vect_in(a_dim)
      double complex, intent(out) :: vect_out(a_dim)
#else
      double precision, intent(in)  :: vect_in(a_dim)
      double precision, intent(out) :: vect_out(a_dim)
#endif
      integer, intent(in) :: power
      integer i, j, ii, jj, n, der, eq, var, a_index, v_index
      integer eqloc, varloc

      ! Preliminary calculations: derivative(s) of vect
      dm(1)%vect_der(1:a_dim,dmat(1)%der_min:dmat(1)%der_max) = 0d0
      do var=1,dm(1)%nvar
        do der=dm(1)%var_der_min(var),dm(1)%var_der_max(var)
          do i=1,grd(1)%nr
            do j=1,nt
              v_index = dm(1)%ivar(var,i,j)
              do ii=max(1,i-dmat(1)%lbder(der)),min(grd(1)%nr,i+dmat(1)%ubder(der))
                a_index = dm(1)%ivar(var,ii,j)
                dm(1)%vect_der(v_index,der) = dm(1)%vect_der(v_index,der) + &
                               dmat(1)%derive(i,ii,der)*vect_in(a_index)
              enddo
            enddo
          enddo
        enddo
      enddo

      vect_out(1:a_dim) = 0d0

      ! The equations:
      do n=1,dm(1)%nas
        if (dm(1)%asi(1,n).eq.power) then
          der = dm(1)%asi(2,n)
          eq  = dm(1)%asi(3,n)
          var = dm(1)%asi(4,n)
          do i=1,grd(1)%nr
            do j=1,nt
              v_index = dm(1)%ieq(eq,i,j)
              if (.not.dm(1)%bc_flag(v_index)) then
                a_index = dm(1)%ivar(var,i,j)
                vect_out(v_index) = vect_out(v_index) + dm(1)%as(n)*dm(1)%vect_der(a_index,der)
              endif
            enddo
          enddo
        endif
      enddo

      do n=1,dm(1)%nart
        if (dm(1)%arti(1,n).eq.power) then
          der = dm(1)%arti(2,n)
          eq  = dm(1)%arti(3,n)
          var = dm(1)%arti(4,n)
          do i=1,grd(1)%nr
            do j=1,nt
              v_index = dm(1)%ieq(eq,i,j)
              if (.not.dm(1)%bc_flag(v_index)) then
                a_index = dm(1)%ivar(var,i,j)
                vect_out(v_index) = vect_out(v_index) + &
                                  dm(1)%art(i,j,n)*dm(1)%vect_der(a_index,der)
              endif
            enddo
          enddo
        endif
      enddo

      do n=1,dm(1)%nartt
        if (dm(1)%artti(1,n).eq.power) then
          der = dm(1)%artti(2,n)
          eq  = dm(1)%artti(3,n)
          var = dm(1)%artti(4,n)
          do i=1,grd(1)%nr
            do j=1,nt
              v_index = dm(1)%ieq(eq,i,j)
              if (.not.dm(1)%bc_flag(v_index)) then
                do jj=1,nt
                  a_index = dm(1)%ivar(var,i,jj)
                  vect_out(v_index) = vect_out(v_index) + &
                                    dm(1)%artt(i,j,jj,n)*dm(1)%vect_der(a_index,der)
                enddo
              endif
            enddo
          enddo
        endif
      enddo

      ! Boundary conditions:
      do n=1,idm(1, 1)%natbc
        if (idm(1, 1)%atbci(1,n).eq.power) then
          der = idm(1, 1)%atbci(2,n)
          eq  = idm(1, 1)%atbci(3,n)
          var = idm(1, 1)%atbci(4,n)
          eqloc  = idm(1, 1)%atbci(5,n)
          varloc = idm(1, 1)%atbci(6,n)
          do j=1,nt
            v_index = dm(1)%ieq(eq,eqloc,j)
            a_index = dm(1)%ivar(var,varloc,j)
            vect_out(v_index) = vect_out(v_index) &
                              + idm(1, 1)%atbc(j,n)*dm(1)%vect_der(a_index,der)
          enddo
        endif
      enddo

      do n=1,idm(1, 1)%nattbc
        if (idm(1, 1)%attbci(1,n).eq.power) then
          der = idm(1, 1)%attbci(2,n)
          eq  = idm(1, 1)%attbci(3,n)
          var = idm(1, 1)%attbci(4,n)
          eqloc  = idm(1, 1)%attbci(5,n)
          varloc = idm(1, 1)%attbci(6,n)
          do j=1,nt
            v_index = dm(1)%ieq(eq,eqloc,j)
            do jj=1,nt
              a_index = dm(1)%ivar(var,varloc,jj)
              vect_out(v_index) = vect_out(v_index) &
                                + idm(1, 1)%attbc(j,jj,n)*dm(1)%vect_der(a_index,der)
            enddo
          enddo
        endif
      enddo

      end subroutine a_product
#endif
!------------------------------------------------------------------------------
#ifndef USE_MULTI
      subroutine dump_coef()
          use mod_grid, only: grd, nt, ndomains
          use model
          integer i
#ifdef USE_1D
          open(42, file="as.coef")
          do i=1, dm(1)%nas
            write(42, *) dm(1)%as(i)
          enddo
          close(42)

          open(42, file="ar.coef")
          do i=1, dm(1)%nar
            write(42, *) dm(1)%ar(1:nr, i)
          enddo
          close(42)

          open(42, file="asbc.coef")
          do i=1, dm(1)%nasbc
            write(42, *) dm(1)%asbc(i)
          enddo
          close(42)

#else
          open(42, file="as.coef")
          do i=1, dm(1)%nas
            write(42, *) dm(1)%as(i)
          enddo
          close(42)

          open(42, file="art.coef")
          do i=1, dm(1)%nart
            write(42, *) dm(1)%art(1:nr, 1:nt, i)
          enddo
          close(42)

          open(42, file="artt.coef")
          do i=1, dm(1)%nartt
            write(42, *) "i=", i, dm(1)%artt(1:nr, 1:nt, 1:nt, i)
          enddo
          close(42)

          open(42, file="atbc.coef")
          do i=1, idm(1, 1)%natbc
            write(42, *) "i=", i, idm(1, 1)%atbc(:, i)
          enddo
          close(42)

          open(42, file="attbc.coef")
          do i=1, idm(1, 1)%nattbc
            write(42, *) "i=", i, idm(1, 1)%attbc(:, :, i)
          enddo
          close(42)
#endif
      end subroutine dump_coef
#endif

      end module matrices
