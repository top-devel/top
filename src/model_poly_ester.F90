module model
      use mod_grid
      use abstract_model_mod
      use inputs

      implicit none

      double precision, allocatable, dimension(:, :), save :: hh, &
          hht, hhz, hhzz, hhzt, lnhht
      double precision, save :: lambda, alpha, aplat, omega_K, omga
      double precision, allocatable, dimension(:, :), save::r_t, r_z,  &
          r_map, re_t, re_z, re_map, r_zz, r_zt, r_tt, re_zz, re_zt,   &
          re_tt, zeta, cost, sint, cott
      double precision, allocatable, dimension(:), save::cth, sth

      type, extends(abstract_model) :: model_poly_ester
      contains
          procedure :: init => init_poly_ester_model
          procedure :: get_field => poly_ester_get_field
      end type model_poly_ester

contains

!------------------------------------------------------------------------
      subroutine poly_ester_get_field(this, fname, field)

          class(model_poly_ester) :: this
          character(len=*), intent(in) :: fname
          real(kind=8), allocatable, intent(out) :: field(:, :)

          if (fname == 'hh') then
              allocate(field(grd(1)%nr, lres))
              field = hh
              return
          endif

          if (fname == 'hht') then
              allocate(field(grd(1)%nr, lres))
              field = hht
              return
          endif

          print*, 'Unknown field:', fname
          allocate(field(1, 1))
          field = 0.0

      end subroutine poly_ester_get_field
!------------------------------------------------------------------------
      subroutine init_poly_ester_model(this, filename)

          class(model_poly_ester), target :: this
          character(len=*), intent(in) :: filename

          call read_model(filename)
          call read_poly()
          call init_radial_grid()
          call write_grid()
          call make_mapping(lres)

          model_ptr => this

      end subroutine init_poly_ester_model
!------------------------------------------------------------------------
      subroutine init_model(filename)

          character(len=*), intent(in) :: filename

          class(model_poly_ester), allocatable :: model

          allocate(model)

          call model%init(filename)

      end subroutine
!------------------------------------------------------------------------
      subroutine read_model(dirmodel)

          use inputs, only: rota, pindex
          use mod_grid

          implicit none
          character(len=*), intent(in) :: dirmodel

          integer :: nr_tmp
          double precision :: tmp, pindex1, pindex2, rota_tmp
          character(len=255) :: filename

          filename = "lambda_eps"
          write(*, *) trim(dirmodel)
          open(unit=2, file=trim(dirmodel)//filename, &
              status="old", action="read")

          read(2, *) lambda
          read(2, *) aplat
          read(2, *) alpha
          read(2, *) tmp
          read(2, *) pindex1
          if (pindex /= pindex1) &
              write(*, &
              '("Warning: new value for pindex: ", f5.3, " -> ", f5.3)') &
              pindex, pindex1
          pindex = pindex1
          read(2, *) pindex2
          read(2, *) rota_tmp
          if (rota /= rota_tmp) &
              write(*, &
              '("Warning: new value for rota: ", f6.4, " -> ", f6.4)') &
              rota, rota_tmp
          rota=rota_tmp
          read(2, *) omega_K
          omga = omega_K/sqrt(3d0*alpha)
          read(2, *) nr_tmp
          if (grd(1)%nr.ne.nr_tmp) &
              write(*, &
              '("Warning: new value for nr: ", i5, " -> ", i5)') &
              grd(1)%nr, nr_tmp
          grd(1)%nr=nr_tmp
          close(2)

      end subroutine read_model

!------------------------------------------------------------------------
      subroutine read_poly()

          use mod_legendre
          use inputs
          implicit none

          integer i, j, l, k, nr_temp, lmod_temp
          double precision, allocatable :: hh_spec(:, :), hhz_spec(:, :)
          double precision, allocatable :: hhzz_spec(:, :)
          double precision aux
          character*(3) str
          character*(512) filename

          filename = "enthalpy_et_der"
          open(unit=2, file=trim(dirmodel)//filename, &
              status="old", &
              action="read")
          read(2, *) nr_temp, lmod_temp  ! this is in case I use a preexisting model
!      if (nr.ne.nr_temp) print*, "Warning: new value for nr"
          write(*, *) nr_temp, lmod_temp
          if (lmod.ne.lmod_temp) &
              write(*, '("Warning: new value for lmod: ", i5, " -> ", i5)') &
              lmod, lmod_temp
          lmod   = lmod_temp

          lmod = lmod - 1

          if (lmod.ge.lres) stop 'lmod.ge.lres: case not implemented'

          if (allocated(hh)) deallocate(hh)
          if (allocated(hhz)) deallocate(hhz)
          if (allocated(hht)) deallocate(hht)
          if (allocated(hhzz)) deallocate(hhzz)
          if (allocated(hhzt)) deallocate(hhzt)
          if (allocated(lnhht)) deallocate(lnhht)
          allocate(hh_spec(grd(1)%nr, lres), hhz_spec(grd(1)%nr, lres), &
              hhzz_spec(grd(1)%nr, lres))
          allocate(hh(grd(1)%nr, lres), hhz(grd(1)%nr, lres), &
              hht(grd(1)%nr, lres), hhzz(grd(1)%nr, lres),  &
              lnhht(grd(1)%nr, lres), hhzt(grd(1)%nr, lres))

          hh_spec = 0d0 ; hhz_spec = 0d0 ; hhzz_spec = 0d0
          hh = 0d0; hhz = 0d0; hht = 0d0 ; hhzz = 0d0; hhzt = 0d0; lnhht = 0d0

          read(2, *) str
          do k=1, lmod/2+1
          read(2, *) i
          read(2, *) (aux, j=1, grd(1)%nr)
          enddo
          read(2, *) str
          do l=0, lmod, 2
          read(2, *) i
          ! verif
          if(i .ne. l) stop 'problem reading h'
          read(2, *) (hh_spec(j, i+1), j=1, grd(1)%nr)
          enddo

          read(2, *) str
          do l=0, lmod, 2
          read(2, *) i
          ! verif
          if(i .ne. l) stop 'problem reading hz'
          read(2, *) (hhz_spec(j, i+1), j=1, grd(1)%nr)
          enddo

          read(2, *) str
          do l=0, lmod, 2
          read(2, *) i
          ! verif
          if(i .ne. l) stop 'problem reading hzz'
          read(2, *) (hhzz_spec(j, i+1), j=1, grd(1)%nr)
          enddo
          close(2)

          call legendre(hh_spec(1:grd(1)%nr, 1:lres), hh(1:grd(1)%nr, 1:lres), &
              lres, grd(1)%nr, 0, -1)

          call legendrep(hh_spec(1:grd(1)%nr, 1:lres), hht(1:grd(1)%nr, 1:lres), &
              lres, grd(1)%nr, 0)

          call legendre(hhz_spec(1:grd(1)%nr, 1:lres), hhz(1:grd(1)%nr, 1:lres), &
              lres, grd(1)%nr, 0, -1)

          call legendrep(hhz_spec(1:grd(1)%nr, 1:lres), hhzt(1:grd(1)%nr, 1:lres), &
              lres, grd(1)%nr, 0)

          call legendre(hhzz_spec(1:grd(1)%nr, 1:lres), hhzz(1:grd(1)%nr, 1:lres), &
              lres, grd(1)%nr, 0, -1)

          deallocate(hh_spec, hhz_spec, hhzz_spec)
          where (hh < 0d0) hh = abs(hh)

              do i=1, grd(1)%nr-1
              do j=1, lres
              lnhht(i, j) = hht(i, j)/hh(i, j)
              enddo
              enddo

              do j=1, lres
              lnhht(grd(1)%nr, j) = hhzt(i, j)/hhz(i, j)
              enddo
          end subroutine

!-----------------------------------------------------------------------
          subroutine write_grid()
              implicit none
              integer i
              double precision, parameter :: pi = 3.14159265358979d0

              do i=1, grd(1)%nr
              grd(1)%r(i) = 0.5d0*(1d0-dcos(dble(i-1)*pi/dble(grd(1)%nr-1)))
              enddo
          end subroutine

!-----------------------------------------------------------------------
          subroutine make_mapping(lres)

              use mod_legendre
              use inputs, only: dirmodel, lmod

              implicit none
              integer, intent(in) :: lres
              integer i, lmod_temp, l
              double precision cnst, pi, xi, aux
              double precision, allocatable, dimension(:)::r_spec, llr_spec
              double precision, allocatable, dimension(:)::rs, rsp, a, ae, ap, &
                  aep, as, aes, w, rss

              double precision, pointer :: r(:)
              character*(512) filename

              ! double precision, pointer :: r(:)
              r => grd(1)%r

              filename = "domains_boundaries"
              open(unit=3, file=trim(dirmodel)//filename, status="old")
              read(3, *) lmod_temp
              if (lmod.ne.lmod_temp) &
                  stop "Error: lmod in domains_boundaries is not compatible"
              if (lmod.ge.lres) stop 'lmod.ge.lres: not implemented'
              if (allocated(r_map)) deallocate ( r_map, r_t,          &
                  r_z, re_map,         &
                  re_t, re_z,          &
                  r_zz, r_zt,          &
                  r_tt, re_zz,         &
                  re_zt, re_tt,        &
                  zeta, cth, sth, cost, sint, cott)

              allocate(r_spec(lres), rs(lres), rsp(lres))
              allocate(llr_spec(lres), rss(lres))
              allocate(r_map(grd(1)%nr, lres), r_t(grd(1)%nr, lres))
              allocate(r_z(grd(1)%nr, lres), re_map(grd(1)%nr, lres))
              allocate(re_t(grd(1)%nr, lres), re_z(grd(1)%nr, lres))
              allocate(r_zz(grd(1)%nr, lres), r_zt(grd(1)%nr, lres))
              allocate(r_tt(grd(1)%nr, lres), re_zz(grd(1)%nr, lres))
              allocate(re_zt(grd(1)%nr, lres), re_tt(grd(1)%nr, lres))
              allocate(zeta(grd(1)%nr, lres), a(grd(1)%nr), ap(grd(1)%nr), as(grd(1)%nr))
              allocate(ae(grd(1)%nr), aep(grd(1)%nr), aes(grd(1)%nr), cth(lres))
              allocate(sth(lres), w(lres), cost(grd(1)%nr, lres))
              allocate(sint(grd(1)%nr, lres), cott(grd(1)%nr, lres))

              r_spec=0d0;rs=0d0;rsp=0d0;r_map=0d0;r_t=0d0;r_z=0d0;re_map=0d0;
              re_t=0d0;re_z=0d0;zeta=0d0;a=0d0;ap=0d0;ae=0d0;aep=0d0;
              read(3, *) (aux, i=1, lmod+1)
              read(3, *) r_spec(1:lmod+1)
              close(3)

              ! calculation of sin(theta), cos(theta)
              call gauleg(-1d0, 1d0, cth, w, lres)
              sth = sqrt(1d0-cth**2)
              do i=1, grd(1)%nr
              sint(i, :) = sth(:)
              cost(i, :) = cth(:)
              cott(i, :) = cth(:)/sth(:)
              enddo

              call legendre(r_spec(1:lres), rs(1:lres), lres, 0, -1)
              call legendrep(r_spec(1:lres), rsp(1:lres), lres, 0)
              do l=1, lres ! Beware: r_spec(l) is r_spec of (l-1)
              llr_spec(l)=-l*(l-1)*r_spec(l)
              enddo
              call legendre(llr_spec(1:lres), rss(1:lres), lres, 0, -1)
              rss = rss - cth*rsp/sth

              pi = dacos(-1d0)
              do i=1, grd(1)%nr
              zeta(i, 1:lres) = r(i)
              enddo

              ! Les nouveaux modeles utilisent le nouveau mapping:
              do i=1, grd(1)%nr
              a(i)  = ( 5d0*r(i)**3 -   3d0*r(i)**5)/2d0
              ap(i) = (15d0*r(i)**2 -  15d0*r(i)**4)/2d0
              as(i) = (30d0*r(i)    -  60d0*r(i)**3)/2d0
              enddo

              ! On simule du bi-domaine, en reutilisant r(i)
              ! Dans le deuxième «domaine» on a:
              !       zeta = 2 - r(i), donc xi = 1-2*r(i)
              do i=1, grd(1)%nr
              xi = 1d0-2d0*r(i)
              ae(i)  = (    xi**3-3d0*xi+2d0)/4d0
              aep(i) =-(3d0*xi**2-3d0       )/2d0
              aes(i) = (6d0*xi              )
              enddo

              cnst = 1d0-aplat
              do i=1, grd(1)%nr
              r_map(i, :)= cnst*r(i)  + a(i) *(rs(:)-cnst)
              r_z(i, :)  = cnst       + ap(i)*(rs(:)-cnst)
              r_t(i, :)  =              a(i) * rsp(:)
              r_zz(i, :) =              as(i)*(rs(:)-cnst)
              r_zt(i, :) =              ap(i)* rsp(:)
              r_tt(i,  :) =              a(i) * rss(:)
              enddo

              ! On simule du bi-domaine, en reutilisant r(i)
              do i=1, grd(1)%nr
              re_map(i, :)=-cnst*r(i)  + 2d0 + ae(i) *(rs(:)-1d0-aplat)
              re_z(i, :)  =-cnst             + aep(i)*(rs(:)-1d0-aplat)
              re_t(i, :)  =                    ae(i) * rsp(:)
              re_zz(i, :) =                    aes(i)*(rs(:)-1d0-aplat)
              re_zt(i, :) =                    aep(i)* rsp(:)
              re_tt(i, :) =                    ae(i) * rss(:)
              enddo
              deallocate(rs, rsp, rss, r_spec, llr_spec, a, ae, ap, aep, as, aes, w)
          end subroutine
!-----------------------------------------------------------------------
      end module
