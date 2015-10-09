#include "config.h"
      module model
      use mod_grid

      double precision,allocatable,dimension(:,:),save::hh,hht,hhz,&
              hhzz, hhzt, lnhht, hhtoz
      double precision,save :: lambda, alpha, aplat, omega_K,omga
      double precision,allocatable,dimension(:,:),save::r_t,r_z,r_map,&
                       re_t,re_z,re_map,r_zz,r_zt,r_tt,re_zz,re_zt,&
                       re_tt,zeta,cost,sint,cott,cotte,roz,rrt
      double precision,allocatable,dimension(:),save::cth,sth
      integer, save :: nr

contains
!------------------------------------------------------------------------
      subroutine init_model(filename)
      use inputs, only: lres
      implicit none
      character(len=*), intent(in) :: filename

      !call run_bckg2D()
      call read_poly(filename)
      call read_lambda(filename)
      call write_grids()
      call make_mapping(lres)
      end subroutine
!------------------------------------------------------------------------
      subroutine run_bckg2D()

      use inputs

      implicit none
      character*100 string1,string2,string3,string4,string5,string6

      write(string1,'("rotation = ",1pd12.5)') rota
      write(string2,'("p_index = ",f6.3)') pindex
      write(string3,'("iz = ",i6)') nr*2-2
      write(string4,'("nzd(1) = ",i6)') nr
      write(string5,'("nzd(2) = ",i6)') nr-1
      write(string6,'("lmax =  ",i6)') lmod
      open(unit=1,file='param',status='unknown')
       write(1,*) string1
       write(1,*) string2
       write(1,*) string3
       write(1,*) 'ndomains = 2'
       write(1,*) 'zd(2) = 1.0'
       write(1,*) string4
       write(1,*) string5
       write(1,*) 'lmin = 0'
       write(1,*) string6
      close(1)

      print*,'fabrication du background'
      call system('cat param dati_bckg2D > dati_poisson2D')
      call system('bckg_full2D < dati_poisson2D > out')

      end subroutine
!------------------------------------------------------------------------
      subroutine read_poly(filename)

      use mod_legendre
      use inputs

      implicit none
      character(len=*), intent(in) :: filename

      integer i,j,l,k,lmod_temp
      double precision pindex_temp
      double precision, allocatable :: hh_spec(:,:),hhz_spec(:,:)
      double precision, allocatable :: hhzz_spec(:,:)
      double precision aux
      character(len=3) str

      open(unit=2,file=filename//"enthalpy",status="old")
      read(2,*) pindex_temp        ! this is in case I use a preexisting model
      if (pindex.ne.pindex_temp) print*,"Warning: new value for pindex"
      pindex = pindex_temp
      read(2,*) nr,lmod_temp  ! this is in case I use a preexisting model
      if (lmod.ne.lmod_temp) print*,"Warning: new value for lmod"
      lmod   = lmod_temp

      lmod = lmod - 1

      if (lmod.ge.lres) stop 'lmod.ge.lres: case not implemented'
      allocate(hh_spec(nr,lres),hhz_spec(nr,lres),hhzz_spec(nr,lres))
      allocate(hh(nr,lres),hhz(nr,lres),hht(nr,lres),hhzz(nr,lres),&
               lnhht(nr,lres),hhzt(nr,lres))

      hh_spec = 0d0 ; hhz_spec = 0d0 ; hhzz_spec = 0d0
      hh = 0d0; hhz = 0d0; hht = 0d0 ; hhzz = 0d0; hhzt = 0d0; lnhht = 0d0

      read(2,*) str
      do k=1,lmod/2+1
          read(2,*) i
          read(2,*) (aux,j=1,nr+2)
      enddo
      read(2,*) str
      do l=0,lmod,2
        read(2,*) i
        read(2,*) (hh_spec(j,i+1),j=1,nr),aux,aux
      enddo
      close(2)

      open(unit=2,file=filename//"enth_der_final",status="old")
      read(2,*) str
      read(2,*) str
      do l=0,lmod,2
        read(2,*) i
        read(2,*) (hhz_spec(j,i+1),j=1,nr),aux,aux
      enddo
      do l=0,lmod,2
        read(2,*) i
        read(2,*) (hhzz_spec(j,i+1),j=1,nr),aux,aux
      enddo
      close(2)

      call legendre(hh_spec(1:nr,1:lres),hh(1:nr,1:lres),&
                    lres,nr,0,-1)

      call legendrep(hh_spec(1:nr,1:lres),hht(1:nr,1:lres),&
                    lres,nr,0)

      call legendre(hhz_spec(1:nr,1:lres),hhz(1:nr,1:lres),&
                    lres,nr,0,-1)

      call legendrep(hhz_spec(1:nr,1:lres),hhzt(1:nr,1:lres),&
                    lres,nr,0)

      call legendre(hhzz_spec(1:nr,1:lres),hhzz(1:nr,1:lres),&
                    lres,nr,0,-1)

      deallocate(hh_spec,hhz_spec, hhzz_spec)
      where (hh < 0d0) hh = 0d0

      do i=1,nr-1
        do j=1,lres
          lnhht(i,j) = hht(i,j)/hh(i,j)
        enddo
      enddo

      do j=1,lres
        lnhht(nr,j) = hhzt(i,j)/hhz(i,j)
      enddo
      end subroutine
!-----------------------------------------------------------------------
      subroutine read_lambda(filename)

      use inputs, only: rota

      implicit none
      character(len=*), intent(in) :: filename

      integer i
      double precision aux
      double precision rota_temp

      open(unit=3,file=filename//"lambda_relax",status="old")
      read(3,*)
      do
        read(3,*,end=404) i, lambda, aplat, omega_K, alpha
      enddo
  404 close(3)

      rota_temp = omega_K*sqrt(lambda/3d0/alpha) ! this is in case we're using a prexisting model
      if (rota.ne.rota_temp) print*,"Warning: new value for rota"
      rota = rota_temp
      omga = omega_K/sqrt(3d0*alpha)

      end subroutine
!-----------------------------------------------------------------------
      subroutine write_grids()
      implicit none
      integer i
      double precision, parameter :: pi = 3.14159265358979d0

      grd(1)%nr = nr
      !grd(2)%nr = nr

      !if (allocated(grd(1)%r)) deallocate(grd(1)%r)
      !if (allocated(grd(2)%r)) deallocate(grd(2)%r)
      allocate(grd(1)%r(grd(1)%nr),grd(2)%r(grd(2)%nr))

      do i=1,grd(1)%nr
        grd(1)%r(i) = 0.5d0*(1d0-dcos(dble(i-1)*pi/dble(grd(1)%nr-1)))
      enddo

      do i=1,grd(2)%nr
        grd(2)%r(i) = 0.5d0*(3d0-dcos(dble(i-1)*pi/dble(grd(2)%nr-1)))
      enddo

      end subroutine
!-----------------------------------------------------------------------
      subroutine make_mapping(lres)

      use mod_legendre

      implicit none
      integer aux,nths,i,j,lres,lmod,l
      character*(3) str
      double precision cnst,pi
      double precision,allocatable,dimension(:)::r_spec,llr_spec
      double precision,allocatable,dimension(:)::rs,rsp,a,ae,ap,&
                                                 aep,as,aes,w,rss
      character*(512) filename

      filename = "surface_b"
      open(unit=3,file=filename,status="old")
      read(3,*) aux
      read(3,*) aux,lmod
      if (lmod.ge.lres) stop 'lmod.ge.lres: not implemented'
      allocate(r_spec(lres),rs(lres),rsp(lres),                    &
               llr_spec(lres),rss(lres),                           &
               r_map(grd(1)%nr,lres),r_t(grd(1)%nr,lres),          &
               r_z(grd(1)%nr,lres),re_map(grd(2)%nr,lres),         &
               re_t(grd(2)%nr,lres),re_z(grd(2)%nr,lres),          &
               r_zz(grd(1)%nr,lres),r_zt(grd(1)%nr,lres),          &
               r_tt(grd(1)%nr,lres),re_zz(grd(2)%nr,lres),         &
               re_zt(grd(2)%nr,lres),re_tt(grd(2)%nr,lres),        &
               zeta(grd(1)%nr,lres),a(grd(1)%nr),ap(grd(1)%nr),    &
               as(grd(1)%nr),ae(grd(2)%nr),aep(grd(2)%nr),         &
               aes(grd(2)%nr),cth(lres),cotte(grd(2)%nr,lres),     &
               sth(lres),w(lres),cost(grd(1)%nr,lres),             &
               sint(grd(1)%nr,lres),cott(grd(1)%nr,lres),          &
               roz(grd(1)%nr,lres),rrt(grd(1)%nr,lres),            &
               hhtoz(grd(1)%nr,lres))

      r_spec=0d0;rs=0d0;rsp=0d0;r_map=0d0;r_t=0d0;r_z=0d0;re_map=0d0;
      re_t=0d0;re_z=0d0;zeta=0d0;a=0d0;ap=0d0;ae=0d0;aep=0d0;cotte=0d0;
      roz=0d0;rrt=0d0;hhtoz=0d0;

      read(3,*) str
      read(3,'(1pe21.14)') r_spec(1:lmod+1:2)
      close(3)

      ! calculation of sin(theta), cos(theta)
      call gauleg(-1d0,1d0,cth,w,lres)
      sth = sqrt(1-cth**2)
      do i=1,grd(1)%nr
        sint(i,:) = sth(:)
        cost(i,:) = cth(:)
        cott(i,:) = cth(:)/sth(:)
      enddo

      do i=1,grd(2)%nr
        cotte(i,:) = cott(1,:)
      enddo

      call legendre(r_spec(1:lres),rs(1:lres),lres,0,-1)
      call legendrep(r_spec(1:lres),rsp(1:lres),lres,0)
      do l=1,lres ! Beware: r_spec(l) is r_spec of (l-1)
        llr_spec(l)=-l*(l-1)*r_spec(l)
      enddo
      call legendre(llr_spec(1:lres),rss(1:lres),lres,0,-1)
      rss = rss - cth*rsp/sth

      pi = dacos(-1d0)
      do i=1,grd(1)%nr
        zeta(i,1:lres) = grd(1)%r(i)
      enddo

      ! Les nouveaux modeles utilisent le nouveau mapping:
      do i=1,grd(1)%nr
        a(i)  = ( 5d0*grd(1)%r(i)**3 -   3d0*grd(1)%r(i)**5)/2d0
        ap(i) = (15d0*grd(1)%r(i)**2 -  15d0*grd(1)%r(i)**4)/2d0
        as(i) = (30d0*grd(1)%r(i)    -  60d0*grd(1)%r(i)**3)/2d0
      enddo

      cnst = 1d0-aplat
      do i=1,grd(1)%nr
        r_map(i,:)= cnst*grd(1)%r(i)  + a(i) *(rs(:)-cnst)
        r_z(i,:)  = cnst       + ap(i)*(rs(:)-cnst)
        r_t(i,:)  =              a(i) * rsp(:)
        r_zz(i,:) =              as(i)*(rs(:)-cnst)
        r_zt(i,:) =              ap(i)* rsp(:)
        r_tt(i,:) =              a(i) * rss(:)
      enddo

      ! Deuxième domaine:
      do i=1,grd(2)%nr
        ae(i)  =  2d0*grd(2)%r(i)**3 -  9d0*grd(2)%r(i)**2 + 12d0*grd(2)%r(i) - 4d0
        aep(i) =  6d0*grd(2)%r(i)**2 - 18d0*grd(2)%r(i)    + 12d0
        aes(i) = 12d0*grd(2)%r(i)    - 18d0
      enddo

      do i=1,grd(2)%nr
        re_map(i,:)= 2d0*aplat + cnst*grd(2)%r(i) + ae(i) *(rs(:)-1d0-aplat)
        re_z(i,:)  =             cnst             + aep(i)*(rs(:)-1d0-aplat)
        re_t(i,:)  =                                ae(i) * rsp(:)
        re_zz(i,:) =                                aes(i)*(rs(:)-1d0-aplat)
        re_zt(i,:) =                                aep(i)* rsp(:)
        re_tt(i,:) =                                ae(i) * rss(:)
      enddo
      deallocate(rs,rsp,rss,r_spec,llr_spec,a,ae,ap,aep,as,aes,w)

      ! roz(1:nr_poly,1:nth_poly)   = r_map/zeta
      ! rrt(1:nr_poly,1:nth_poly)   = r_t/r_map
      ! hhtoz(1:nr_poly,1:nth_poly) = hht/zeta
      do j=1,lres
        do i=2,grd(1)%nr
          roz(i,j)   = r_map(i,j)/zeta(i,j)
          rrt(i,j)   = r_t(i,j)/r_map(i,j)
          hhtoz(i,j) = hht(i,j)/zeta(i,j)
        enddo
        roz(1,j)   = r_z(1,j)
        rrt(1,j)   = r_zt(1,j)/r_z(1,j)
        hhtoz(1,j) = hhzt(1,j)
      enddo

      end subroutine
!-----------------------------------------------------------------------
      end module
