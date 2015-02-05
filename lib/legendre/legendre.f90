        module mod_legendre

        use fast_legendre

        implicit none
        private ! by default everything is private. Don't pollute name space!
        public :: legendre,legendrep,legendrepp  ! list of exported routines
        interface legendre
           module procedure ds1D_legendre
           module procedure ds_legendre
           module procedure zs_legendre
           module procedure vr_legendre
           module procedure v_legendre
        end interface
        interface legendrep
           module procedure ds1D_legendrep
           module procedure ds_legendrep
           module procedure zs_legendrep
        end interface
        interface legendrepp
           module procedure ds1D_legendrepp
           module procedure ds_legendrepp
           module procedure zs_legendrepp
        end interface

contains
!------------------------------------------------------------------------
! performs legendre transform on a 1D real scalar field
! f  is the field in real space
! tf its projection on spherical harmonics of order m
!------------------------------------------------------------------------

        subroutine ds1D_legendre(f,tf,nth,m,index)

        implicit none

        integer, intent(in) :: m, nth, index
        double precision, intent(in) :: f(nth)
        double precision, intent(out):: tf(nth)
        double precision, allocatable :: cth(:), w(:)

        tf = 0d0 ! VERY important
        allocate(cth(nth),w(nth))
        call gauleg(-1d0,1d0,cth,w,nth)
        if (m.ne.0) print*,'beware! your field is real and m#0'
        if (index.eq.1) then       ! real-spectral transform
          call project_ylm(1d0,f,cth,w,tf,nth,nth,0,1,m)
        else if (index.eq.-1) then ! spectral-real transform
          call eval_ylm(1d0,f,tf,cth,nth,nth,0,1,m)
        else
          stop 'There is a problem in the value of index'
        endif
        deallocate(cth,w)

        end subroutine ds1D_legendre

!------------------------------------------------------------------------
! performs legendre transform on a 2D real scalar field
! f  is the field in real space
! tf its projection on spherical harmonics of order m
!------------------------------------------------------------------------

        subroutine ds_legendre(f,tf,nth,nr,m,index)

        implicit none

        integer, intent(in) :: nth, nr, m, index
        double precision, intent(in) :: f(nr,nth)
        double precision, intent(out):: tf(nr,nth)
        double precision, allocatable :: cth(:), w(:)

        tf = 0d0
        allocate(cth(nth),w(nth))
        call gauleg(-1d0,1d0,cth,w,nth) 
        if (m.ne.0) print*,'beware! your field is real and m#0'
        if (index.eq.1) then       ! real-spectral transform
          call project_ylm(1d0,f,cth,w,tf,nr,nth,nth,0,1,m)
        else if (index.eq.-1) then ! spectral-real transform
          call eval_ylm(1d0,f,tf,cth,nr,nth,nth,0,1,m)
        else
          stop 'There is a problem in the value of index'
        endif
        deallocate(cth,w)

        end subroutine ds_legendre

!------------------------------------------------------------------------
! performs legendre transform on a 2D complex scalar field
! f  is the field in real space
! tf its projection on spherical harmonics of order m
!------------------------------------------------------------------------

        subroutine zs_legendre(f,tf,nth,nr,m,index)

        implicit none

        integer, intent(in) :: nth, nr, m, index
        double complex, intent(in) :: f(nr,nth)
        double complex, intent(out):: tf(nr,nth)
        double precision, allocatable :: cth(:), w(:)

        tf = (0d0,0d0)
        allocate(cth(nth),w(nth))
        call gauleg(-1d0,1d0,cth,w,nth)
        if (index.eq.1) then       ! real-spectral transform
          call project_ylm((1d0,0d0),f,cth,w,tf,nr,nth,nth,0,1,m)
        else if (index.eq.-1) then ! spectral-real transform
          call eval_ylm((1d0,0d0),f,tf,cth,nr,nth,nth,0,1,m)
        else
          stop 'There is a problem in the value of index'
        endif
        deallocate(cth,w)

        end subroutine zs_legendre

!------------------------------------------------------------------------
! performs legendre transform on a 2D vector field in spherical coordinates
! f  is the field in real space:
! f(:,:,1) = Vr
! f(:,:,2) = Vtheta
! i.f(:,:,3) = Vphi
!
! tf its projection on spherical harmonics of order m :
!
! tf(:,:,1) = u^l_m
! tf(:,:,2) = v^l_m
! -i.tf(:,:,3) = w^l_m
!------------------------------------------------------------------------
        subroutine vr_legendre(f,tf,nth,nr,m,index)

        implicit none

        integer, intent(in) :: nth, nr, m, index
        double precision, intent(in) :: f(nr,nth,3)
        double precision, intent(out):: tf(nr,nth,3)
        double precision, allocatable :: cth(:), w(:)

        tf = 0d0
        allocate(cth(nth),w(nth))
        call gauleg(-1d0,1d0,cth,w,nth)
        if (m.ne.0) print*,'beware! your field is real and m#0'
        if (index.eq.1) then       ! performs real-spectral transform
          call project_ylm(1d0,f(:,:,1),cth,w,tf(:,:,1),nr,nth,nth,0,1,m)
          call project_dth_ylm(1d0,f(:,:,2),cth,w,tf(:,:,2),nr,nth,nth,0,1,m)
          call project_Dphi_ylm(1d0,f(:,:,3),cth,w,tf(:,:,2),nr,nth,nth,0,1,m)
          call project_Dphi_ylm(1d0,f(:,:,2),cth,w,tf(:,:,3),nr,nth,nth,0,1,m)
          call project_dth_ylm(1d0,f(:,:,3),cth,w,tf(:,:,3),nr,nth,nth,0,1,m)
        else if (index.eq.-1) then ! performs spectral-real transform
          call eval_ylm(1d0,f(:,:,1),tf(:,:,1),cth,nr,nth,nth,0,1,m)
          call eval_dth_ylm(1d0,f(:,:,2),tf(:,:,2),cth,nr,nth,nth,0,1,m)
          call eval_Dphi_ylm(1d0,f(:,:,3),tf(:,:,2),cth,nr,nth,nth,0,1,m)
          call eval_Dphi_ylm(1d0,f(:,:,2),tf(:,:,3),cth,nr,nth,nth,0,1,m)
          call eval_dth_ylm(1d0,f(:,:,3),tf(:,:,3),cth,nr,nth,nth,0,1,m)
        else
          stop 'There is a problem in the value of index'
        endif
        deallocate(cth,w)

        end subroutine vr_legendre

!------------------------------------------------------------------------
! performs legendre transform on a 2D vector field in spherical coordinates
! f  is the field in real space:
! f(:,:,1) = Vr
! f(:,:,2) = Vtheta
! f(:,:,3) = Vphi
!
! tf its projection on spherical harmonics of order m :
!
! tf(:,:,1) = u^l_m
! tf(:,:,2) = v^l_m
! tf(:,:,3) = w^l_m
!------------------------------------------------------------------------
        subroutine v_legendre(f,tf,nth,nr,m,index)

        implicit none

        integer, intent(in) :: nth, nr, m, index
        double complex, intent(in) :: f(nr,nth,3)
        double complex, intent(out):: tf(nr,nth,3)
        double precision, allocatable :: cth(:), w(:)

        tf = (0d0,0d0)
        allocate(cth(nth),w(nth))
        call gauleg(-1d0,1d0,cth,w,nth)
        if (index.eq.1) then       ! performs real-spectral transform
          call project_ylm((1d0,0d0),f(:,:,1),cth,w,tf(:,:,1),nr,nth,nth,0,1,m)
          call project_dth_ylm((1d0,0d0),f(:,:,2),cth,w,tf(:,:,2),nr,nth,nth,0,1,m)
          call project_Dphi_ylm((0d0,-1d0),f(:,:,3),cth,w,tf(:,:,2),nr,nth,nth,0,1,m)
          call project_Dphi_ylm((0d0,-1d0),f(:,:,2),cth,w,tf(:,:,3),nr,nth,nth,0,1,m)
          call project_dth_ylm((-1d0,0d0),f(:,:,3),cth,w,tf(:,:,3),nr,nth,nth,0,1,m)
        else if (index.eq.-1) then ! performs spectral-real transform
          call eval_ylm((1d0,0d0),f(:,:,1),tf(:,:,1),cth,nr,nth,nth,0,1,m)
          call eval_dth_ylm((1d0,0d0),f(:,:,2),tf(:,:,2),cth,nr,nth,nth,0,1,m)
          call eval_Dphi_ylm((0d0,1d0),f(:,:,3),tf(:,:,2),cth,nr,nth,nth,0,1,m)
          call eval_Dphi_ylm((0d0,1d0),f(:,:,2),tf(:,:,3),cth,nr,nth,nth,0,1,m)
          call eval_dth_ylm((-1d0,0d0),f(:,:,3),tf(:,:,3),cth,nr,nth,nth,0,1,m)
        else
          stop 'There is a problem in the value of index'
        endif
        deallocate(cth,w)

        end subroutine v_legendre

!------------------------------------------------------------------------
! Performs a 'd-legendre' transform on a l-spectral 1D real scalar field
! It yields the derivative in theta of the 1D field
!
! f  is the field in spectral space
! tf is in real space
!------------------------------------------------------------------------

        subroutine ds1D_legendrep(f,tf,nth,m)

        implicit none

        integer, intent(in) :: nth, m
        double precision, intent(in) :: f(nth)
        double precision, intent(out):: tf(nth)
        double precision, allocatable :: cth(:), w(:)

        tf = 0d0
        allocate(cth(nth),w(nth))
        call gauleg(-1d0,1d0,cth,w,nth)
        if (m.ne.0) print*,'beware! your field is real and m#0'
        ! Always spectral-real transform
        call eval_dth_ylm(1d0,f,tf,cth,nth,nth,0,1,m)
        deallocate(cth,w)

        end subroutine ds1D_legendrep

!------------------------------------------------------------------------
! Performs a 'd-legendre' transform on a l-spectral 2D real scalar field
! It yields the derivative in theta of the 2D field
!
! f  is the field in spectral space
! tf is in real space
!------------------------------------------------------------------------

        subroutine ds_legendrep(f,tf,nth,nr,m)

        implicit none

        integer, intent(in) :: nth, nr, m
        double precision, intent(in) :: f(nr,nth)
        double precision, intent(out):: tf(nr,nth)
        double precision, allocatable :: cth(:), w(:)

        tf = 0d0
        allocate(cth(nth),w(nth))
        call gauleg(-1d0,1d0,cth,w,nth)
        if (m.ne.0) print*,'beware! your field is real and m#0'
        ! Always spectral-real transform
        call eval_dth_ylm(1d0,f,tf,cth,nr,nth,nth,0,1,m)
        deallocate(cth,w)

        end subroutine ds_legendrep

!------------------------------------------------------------------------
! Performs a 'd-legendre' transform on a l-spectral 1D real scalar field
! It yields the derivative in theta of the 1D field
!
! f  is the field in spectral space
! tf is in real space
!------------------------------------------------------------------------

        subroutine zs_legendrep(f,tf,nth,nr,m)

        implicit none

        integer, intent(in) :: nth, nr, m
        double complex, intent(in) :: f(nr,nth)
        double complex, intent(out):: tf(nr,nth)
        double precision, allocatable :: cth(:), w(:)

        tf = (0d0,0d0)
        allocate(cth(nth),w(nth))
        call gauleg(-1d0,1d0,cth,w,nth)
        ! Always spectral-real transform
        call eval_dth_ylm((1d0,0d0),f,tf,cth,nr,nth,nth,0,1,m)
        deallocate(cth,w)

        end subroutine zs_legendrep

!------------------------------------------------------------------------
! Performs a 'd2-legendre' transform on a l-spectral 1D real scalar field
! It yields the second derivative in theta of the 1D field
!
! f  is the field in spectral space
! tf is in real space
!------------------------------------------------------------------------

        subroutine ds1D_legendrepp(f,tf,nth,m)

        implicit none

        integer, intent(in) :: nth, m
        double precision, intent(in) :: f(nth)
        double precision, intent(out):: tf(nth)
        double precision, allocatable :: cth(:), w(:)

        tf = 0d0
        allocate(cth(nth),w(nth))
        call gauleg(-1d0,1d0,cth,w,nth)
        if (m.ne.0) print*,'beware! your field is real and m#0'
        ! Always spectral-real transform
        call eval_ddth_ylm(1d0,f,tf,cth,nth,nth,0,1,m)
        deallocate(cth,w)

        end subroutine ds1D_legendrepp

!------------------------------------------------------------------------
! Performs a 'd2-legendre' transform on a l-spectral 2D real scalar field
! It yields the second derivative in theta of the 2D field
!
! f  is the field in spectral space
! tf is in real space
!------------------------------------------------------------------------

        subroutine ds_legendrepp(f,tf,nth,nr,m)

        implicit none

        integer, intent(in) :: nth, nr, m
        double precision, intent(in) :: f(nr,nth)
        double precision, intent(out):: tf(nr,nth)
        double precision, allocatable :: cth(:), w(:)

        tf = 0d0
        allocate(cth(nth),w(nth))
        call gauleg(-1d0,1d0,cth,w,nth)
        if (m.ne.0) print*,'beware! your field is real and m#0'
        ! Always spectral-real transform
        call eval_ddth_ylm(1d0,f,tf,cth,nr,nth,nth,0,1,m)
        deallocate(cth,w)

        end subroutine ds_legendrepp

!------------------------------------------------------------------------
! Performs a 'd2-legendre' transform on a l-spectral 1D real scalar field
! It yields the second derivative in theta of the 1D field
!
! f  is the field in spectral space
! tf is in real space
!------------------------------------------------------------------------

        subroutine zs_legendrepp(f,tf,nth,nr,m)

        implicit none

        integer, intent(in) :: nth, nr, m
        double complex, intent(in) :: f(nr,nth)
        double complex, intent(out):: tf(nr,nth)
        double precision, allocatable :: cth(:), w(:)

        tf = (0d0,0d0)
        allocate(cth(nth),w(nth))
        call gauleg(-1d0,1d0,cth,w,nth)
        ! Always spectral-real transform
        call eval_ddth_ylm((1d0,0d0),f,tf,cth,nr,nth,nth,0,1,m)
        deallocate(cth,w)

        end subroutine zs_legendrepp

!------------------------------------------------------------------------
        end module mod_legendre
