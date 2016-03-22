#include "config.h"
module eigensolve

      use mod_grid
      use mod_blacs
      use matrices
      use inputs, only: nsol, lres
      use iso_c_binding

      type MULTI_MAT

#ifdef USE_COMPLEX
          double complex, pointer :: mat(:,:)
#else
          double precision, pointer :: mat(:,:)
#endif
          integer, pointer        :: ipiv(:)
#ifdef USE_MPI
          integer, pointer        :: desc(:)
#endif
      end type MULTI_MAT

#ifdef USE_COMPLEX
      double complex, save, allocatable :: omega(:), vec(:,:)
#else
      double precision, save, allocatable :: omega(:), vec(:,:)
#endif
      integer, save :: lda, t_dim, nsol_out
#ifdef USE_COMPLEX
      double complex, save, allocatable :: temp(:),tmp_s(:), &
          temp2(:,:)
#else
      double precision, save, allocatable :: temp(:),tmp_s(:), &
          temp2(:,:)
#endif

      type(MULTI_MAT), save, allocatable :: asigma(:)
      integer ku, kl

contains

!--------------------------------------------------------------
! This subroutine is a wrapper around the Arnoldi-Chebyshev
! program ARNCHEB, written by Braconnier.  It provides the
! necessary work space, the matrix-vector/vector-vector
! products needed to run the ARNCHEB program and the
! arrays which contain the output eigensolutions.
!--------------------------------------------------------------
! Variables:
!
! sigma = eigenvalue shift
!--------------------------------------------------------------

      subroutine run_arncheb(sigma) bind(c)

          implicit none

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Input variables:

#ifdef USE_COMPLEX
          double complex, intent(in) :: sigma
#else
          real(kind=c_double), intent(in) :: sigma
#endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Variables used in the LU factorization
!

          integer info_lapack

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Variables used in the arncheb code:
!

          ! variables in zarncheb.f
          integer iarn, deg, iparam(10), revcom, info
          integer, allocatable :: iord(:)
          double precision tol, normA
#ifdef USE_COMPLEX
          double complex shift
          double complex, allocatable, dimension (:) :: &
              u, w, work, work1, yc
#else
          double precision shift
          double precision, allocatable, dimension (:) :: &
              u, wr, wi, work, y
          double precision, allocatable, dimension (:,:) :: &
              work1
#endif
          double precision, allocatable, dimension (:) :: yr
#ifdef USE_COMPLEX
          double complex, allocatable, dimension (:,:) :: &
              vrcom, v, h, z, workc
#else
          double precision, allocatable, dimension (:,:) :: &
              vrcom, v, h, z, workc
#endif

          ! additional variables in zeigvec.f
#ifdef USE_COMPLEX
          double complex, allocatable, dimension(:,:) :: z1,z2,vec_arncheb
          double complex, allocatable, dimension(:)   :: work2
#else
          double precision, allocatable, dimension(:,:) :: z1,z2,vec_arncheb
          double precision, allocatable, dimension(:)   :: work2
#endif
          double precision, allocatable, dimension(:) :: work2real
          logical, allocatable :: select(:)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Other variables and functions:
!

#ifdef USE_COMPLEX
          double complex aux
#else
          double precision aux
#endif
#ifdef USE_COMPLEX
          double complex, external :: zdotc
#else
          double complex, external :: ddot
#endif
#ifdef USE_MPI
          integer, external :: NUMROC
#endif

          integer iteration, i, j, var, id, d_dim

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! iparam(1)
! 0:Largest Real Part
! 1:Smallest Real Part
! 2:Largest Imaginary Part
! 3:Largest Modulii
! 4:Around a given Shift
          iparam(1)=4

! iparam(2)
! 0:||A|| is computed by n matrix-vector products
! 1:||A|| is given
! 2:||A|| is computed using Higham's modifications of Hager's Algorithm
          iparam(2)=2

! iparam(3)
! 0:deg is dynamically set by the code
! 1:deg is given and fixed during all the computation
          iparam(3)=0

! iparam(4)
! 0:tol is set to the n*machine precision
! 1:tol is given
          iparam(4)=0

! 0:the starting vector is a random vector
! 1:the starting vector is given ([1 1 1 1 ..... 1 1])
          iparam(5)=0

! iparam(6)
! number of wanted eigenvalues
          iparam(6)=nsol

! standard out for error message and computation progress 1=screen
          iparam(7)=1

! standard out for eigenvalue information
          iparam(8)=0

! standard out for eigenvector information
          iparam(9)=0

! maximum number of iterations
          iparam(10)=100

! Please note:
! deg the degree of the chebyshev polynomial,only referenced if iparam(3)=1
! tol the tolerance, only referenced if iparam(4)=1
! normA only referenced if iparam(2)=1

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Create workspaces by allocating various arrays:
!
          if (.not. allocated(asigma)) &
              allocate(asigma(ndomains))

          if (allocated(omega))  deallocate(omega)
          if (allocated(vec))    deallocate(vec)
          if (allocated(temp))   deallocate(temp)
          if (allocated(temp2))  deallocate(temp2)
          if (allocated(tmp_s))  deallocate(tmp_s)

!  For the arncheb program:
          iarn = 40
          shift = sigma
          t_dim = a_dim * power_max

#ifdef USE_MPI
          if (iproc.eq.0) then
#endif
#ifdef USE_COMPLEX
              allocate(w(iarn+1), yr(iarn+1),yc(iarn+1))
              allocate(work1(iarn+1))
#else
              allocate(wr(iarn+1), wi(iarn+1), y(2*iarn+2))
              allocate(work1(iarn+1, 2))
#endif
              allocate(u(t_dim), vrcom(t_dim,4), iord(t_dim))
              allocate(v(t_dim, iarn+1), h(iarn+1, iarn+1))
              allocate(z(iarn+1,iarn+1))
              allocate(work(iarn+1), workc(t_dim,2))
#ifdef USE_MPI
          endif
#endif

!  Extra space needed for the matrix vector products
          allocate(temp(1:a_dim),temp2(1:a_dim,0:power_max-1))
#ifdef USE_MULTI
          allocate(tmp_s(1:d_dim_max))
#endif
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Factorize the matrix A -shift*B
!

#ifdef USE_MULTI
          do id=1,ndomains
          if (grd(id)%mattype.eq.'BAND') then
#ifdef USE_MPI
              if (iproc.eq.0) then
                  print*,"Domain: ",id
                  print*,"Banded matrices: currently unimplemented in MPI TOP"
              endif
              call BLACS_ABORT(ictxt,1)
#else
              call init_ku_kl(id)
              call increase_ku_kl_downward(id)
#endif
          elseif (grd(id)%mattype.ne.'FULL') then
#ifdef USE_MPI
              if (iproc.eq.0) then
#endif
                  print*,"Domain: ",id
                  print*,"Faulty mattype: ",grd(id)%mattype
#ifdef USE_MPI
              endif
              call BLACS_ABORT(ictxt,1)
#else
              stop "wrong mattype"
#endif
          endif
          enddo
#else
          if (grd(1)%mattype.eq.'BAND') call init_ku_kl()
#endif

          ! downward sweep method:
          do id=1,ndomains
          d_dim = dm(id)%d_dim
#ifdef USE_MPI
          dm(id)%vlsize = NUMROC(dm(id)%d_dim ,blk,irow,0,nrows)
          dm(id)%hlsize = NUMROC(dm(id)%d_dim ,blk,icol,0,ncols)
          allocate(asigma(id)%ipiv(dm(id)%vlsize +blk))
          allocate(asigma(id)%desc(9))
#else
          allocate(asigma(id)%ipiv(dm(id)%d_dim ))
#endif
          if (grd(id)%mattype.eq.'FULL') then
#ifdef USE_MPI
              allocate(asigma(id)%mat(dm(id)%vlsize ,dm(id)%hlsize ))
              call DESCINIT(asigma(id)%desc,d_dim,d_dim,blk,blk,0,0,ictxt, &
                  dm(id)%vlsize ,info_lapack)
              if (info_lapack.ne.0) then
                  print*,info_lapack,' info'
                  print*,'DESCINIT problem in domain ',id
                  call BLACS_ABORT(ictxt,1)
              endif
#else
#ifdef USE_MULTI
              allocate(asigma(id)%mat(d_dim, d_dim))
#else
              allocate(asigma(id)%mat(a_dim, a_dim))
#endif
#endif
#ifdef USE_MULTI
              call make_asigma_full_local(sigma,asigma(id)%mat, id, id)
              if (id.gt.1) call correct_asigma_downward(id,sigma)

#else
              call make_asigma_full_total(sigma,asigma(id)%mat)
#endif

#ifdef USE_COMPLEX
              call ZGETRF(d_dim,d_dim,asigma(id)%mat,d_dim,asigma(id)%ipiv,info_lapack)
#else
              call DGETRF(d_dim,d_dim,asigma(id)%mat,d_dim,asigma(id)%ipiv,info_lapack)
#endif
              if (info_lapack.ne.0) then
                  print*, info_lapack,' info'
                  print*, 'Factorisation (1) problem in domain ', id
                  stop
              endif
          elseif (grd(id)%mattype.eq.'BAND') then
              ku = dm(id)%ku
              kl = dm(id)%kl
              lda = 2*kl + ku + 1
              allocate(asigma(id)%mat(lda,d_dim))
#ifdef USE_MULTI
              call make_asigma_band_local(sigma,asigma(id)%mat,id)
              if (id.gt.1) call correct_asigma_downward(id,sigma)
#else
              call make_asigma_band(sigma,asigma(id)%mat)
#endif
          else
              print*, "mattype:", grd(id)%mattype
              stop 'faulty matrix type'

#if USE_MPI
#ifdef USE_COMPLEX
              call PZGETRF(d_dim,d_dim,asigma(id)%mat,1,1,asigma(id)%desc, &
                  asigma(id)%ipiv,info_lapack)
              if (info_lapack.ne.0) then
                  print*,info_lapack,' info'
                  print*,'Factorisation (2) problem in domain ',id
                  call BLACS_ABORT(ictxt,1)
              endif
#else
              call PDGETRF(d_dim,d_dim,asigma(id)%mat,1,1,asigma(id)%desc, &
                  asigma(id)%ipiv,info_lapack)
              if (info_lapack.ne.0) then
                  if (iproc.eq.0) then
                      print*,info_lapack,' info'
                      print*,'Factorisation (2) problem in domain ',id
                  endif
                  call BLACS_ABORT(ictxt,1)
              endif
#endif
#else
#ifdef USE_COMPLEX
              call ZGBTRF(d_dim,d_dim,kl,ku,asigma(id)%mat,lda,asigma(id)%ipiv,info_lapack)
#else
              call DGBTRF(d_dim,d_dim,kl,ku,asigma(id)%mat,lda,asigma(id)%ipiv,info_lapack)
#endif
              if (info_lapack.ne.0) then
                  print*,info_lapack,' info'
                  print*,'Factorisation (3) problem in domain ',id
                  stop
              endif
#endif
          endif
          enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Reverse communication interface
!

          iteration=1
          revcom=0
! big loop starts here ------------------------------------------
          do
#ifdef USE_MPI
          if (iproc.eq.0) then
#endif
#ifdef USE_COMPLEX
              call zarncheb(t_dim,iarn,deg,iparam,tol,shift,normA,u,  &
                  vrcom,iord,v,h,z,w,work,work1,workc,  &
                  yr,yc,revcom,info)
#else
              call darncheb(t_dim,iarn,deg,iparam,tol,shift,normA,u,  &
                  vrcom,iord,v,h,z,wr,wi,work,work1,workc,  &
                  y,revcom,info)
#endif
#ifdef USE_MPI
              call IGEBS2D(ictxt,'A',topo,1,1,revcom,1)
          else
              call IGEBR2D(ictxt,'A',topo,1,1,revcom,1,0,0)
          endif
#endif

          if (revcom.le.0) exit

          if (revcom.eq.1) then
              ! Compute col3 <-- B*col1
              call b_times(vrcom(1:t_dim,1), vrcom(1:t_dim,3))
              ! Compute col3 <-- inv(A-shift*B)*col3
              call solve_amsigmab(sigma, vrcom(1:t_dim,3))
          endif

          if (revcom.eq.2) then
              ! Compute col3 <-- inv(A-shift*B)^H*col3
              call solve_amsigmab_transpose(sigma,vrcom(1:t_dim,3))
              ! Compute col4 <-- B^H*col3
              call b_times_transpose(vrcom(1:t_dim,3), vrcom(1:t_dim,4))
          endif

          if (revcom.eq.3) then
#ifdef USE_MPI
              if (iproc.eq.0) then
#endif
#ifdef USE_COMPLEX
                  aux = zdotc(t_dim,vrcom(1:t_dim,1), 1, vrcom(1:t_dim,2), 1)
#else
                  aux = ddot(t_dim,vrcom(1:t_dim,1), 1, vrcom(1:t_dim,2), 1)
#endif
                  vrcom(1,3) = aux
#ifdef USE_MPI
              endif
#endif
          endif

          iteration=iteration+1

          enddo
! big loop ends here ------------------------------------------
#ifdef USE_MULTI
          do id=1,ndomains
          deallocate(asigma(id)%mat,asigma(id)%ipiv)
          enddo
          deallocate(vrcom,asigma)
#else
          deallocate(vrcom)
          deallocate(asigma(1)%mat)
#endif
#ifdef USE_MPI
          if (iproc.eq.0) then
#endif
#ifdef USE_COMPLEX
              deallocate(u,iord,work,work1,workc,yr,yc)
#else
              deallocate(u,iord,work,work1,workc,y)
#endif
#ifdef USE_MPI
          endif
#endif

          ! print and store the eigenvalues:
          nsol_out = iparam(6) ! the number solutions could change
          allocate(omega(nsol_out))
          do i=1,nsol_out
#ifdef USE_COMPLEX
          omega(i)=w(i)
#else
          omega(i)=dcmplx(wr(i), wi(i))
#endif
          print*,omega(i)
          enddo

!
!  Compute the corresponding eigenvectors
          if (info /= 0) then
              print*, "at iteration ", iteration
              stop "Failure in the Arnoldi-Chebyshev process"
          endif

#ifdef USE_COMPLEX
          allocate(work2(3*iarn+3),z1(iarn+1,iarn+1),z2(iarn+1,iarn+1), &
              work2real(iarn+1),select(iarn+1),vec_arncheb(t_dim,iarn+1))
#else
          allocate(work2(3*iarn+3),z1(iarn+1,iarn+1),z2(iarn+1,iarn+1), &
              vec_arncheb(t_dim,iarn+1),select(iarn+1))
#endif

          allocate(vec(a_dim,nsol_out))

#ifdef USE_COMPLEX
          call zeigvec(t_dim,iarn,h,z,v,vec_arncheb,w,work2,work2real,  &
              select,z1,z2,iparam)
#else
          call deigvec(t_dim,iarn,h,z,v,vec_arncheb,wi,work2,select, &
              z1,z2,iparam)
#endif
          do i=1,nsol_out
          vec(1:a_dim,i) = vec_arncheb(1:a_dim,i)
          enddo

#ifdef USE_COMPLEX
          deallocate(h,z,v,w,work2,work2real,select,z1,z2,vec_arncheb)
#else
          deallocate(h,z,v,vec_arncheb,wr,wi,work2,select,z1,z2)
#endif
          deallocate(temp,temp2)

          return
      end subroutine run_arncheb

!--------------------------------------------------------------
! This routine does:
!    vec2 = "B" * vec1
!--------------------------------------------------------------
! Variables:
!
! vec1 = input vector
! vec2 = output vector
!--------------------------------------------------------------
      subroutine b_times(vec1, vec2)

          implicit none
#ifdef USE_COMPLEX
          double complex, intent(in)  :: vec1(t_dim)
          double complex, intent(out) :: vec2(t_dim)
          double complex :: zero = (0d0, 0d0)
#else
          double precision, intent(in)  :: vec1(t_dim)
          double precision, intent(out) :: vec2(t_dim)
          double precision :: zero = 0d0
#endif
          integer i

#ifdef USE_MPI
          if (iproc.eq.0) then
              call ZGEBS2D(ictxt1D,"A",topo,t_dim,1,vec1,t_dim)
          else
              call ZGEBR2D(ictxt1D,"A",topo,t_dim,1,vec1,t_dim,0,0)
              if (iproc.eq.0) vec2(1:t_dim) = 0d0
          endif
#else
          vec2(1:t_dim) = zero
          temp(1:a_dim) = zero
#endif
          do i=1,power_max
#ifdef USE_MULTI
          call a_product_total(vec1(1+(i-1)*a_dim:i*a_dim), temp(1:a_dim), i)
#else
          call a_product(vec1(1+(i-1)*a_dim:i*a_dim), temp(1:a_dim), i)
#endif
#ifdef USE_MPI
          if (iproc.eq.0) then
#endif
              vec2(1:a_dim) = vec2(1:a_dim) - temp(1:a_dim)
#ifdef USE_MPI
          endif
#endif
          enddo
#ifdef USE_MPI
          if (iproc.eq.0) then
#endif
              do i=2,power_max
              vec2(1+(i-1)*a_dim:i*a_dim) = vec1(1+(i-2)*a_dim:(i-1)*a_dim)
              enddo
#ifdef USE_MPI
          endif
#endif

          return
      end subroutine b_times

!--------------------------------------------------------------
! This routine does:
!    vec2 = transpose("B") * vec1
!--------------------------------------------------------------
! Variables:
!
! vec1 = input vector
! vec2 = output vector
!--------------------------------------------------------------
      subroutine b_times_transpose(vec1, vec2)

          implicit none
#ifdef USE_COMPLEX
          double complex, intent(in)  :: vec1(t_dim)
          double complex, intent(out) :: vec2(t_dim)
          double complex :: zero = (0d0, 0d0)
#else
          double precision, intent(in)  :: vec1(t_dim)
          double precision, intent(out) :: vec2(t_dim)
          double precision :: zero = 0d0
#endif
          integer i

#ifdef USE_MPI
          if (iproc.eq.0) then
              call ZGEBS2D(ictxt1D,"A",topo,t_dim,1,vec1,t_dim)
          else
              call ZGEBR2D(ictxt1D,"A",topo,t_dim,1,vec1,t_dim,0,0)
          endif
#endif
          vec2(1:t_dim) = zero
          do i=1,power_max
#ifdef USE_MULTI
          call a_product_total_transpose(vec1(1:a_dim),&
              vec2(1+(i-1)*a_dim:i*a_dim),i)
#else
          call a_product_transpose(vec1(1:a_dim),&
              vec2(1+(i-1)*a_dim:i*a_dim),i)
#endif
          enddo
#ifdef USE_MPI
          if (iproc.eq.0) then
#endif
              vec2(1:t_dim) = -vec2(1:t_dim)
              do i=1,power_max-1
              vec2(1+(i-1)*a_dim:i*a_dim) = vec2(1+(i-1)*a_dim:i*a_dim) &
                  + vec1(1+i*a_dim:(i+1)*a_dim)
              enddo
#ifdef USE_MPI
          endif
#endif

          return
      end subroutine b_times_transpose

!--------------------------------------------------------------
! This routine does :
!   vect <- inv(A-sigma.B) vect
! where (A-sigma.B) corresponds to the total system, taking
! into account the powers of sigma
!--------------------------------------------------------------
! Variables:
!
! sigma = the eigenvalue shift
! vect  = input and output vector
!--------------------------------------------------------------
      subroutine solve_amsigmab(sigma, vect)

          implicit none
#ifdef USE_COMPLEX
          double complex, intent(in)    :: sigma
          double complex, intent(inout) :: vect(t_dim)
          double complex :: zero = (0d0, 0d0)
#else
          double precision, intent(in)    :: sigma
          double precision, intent(inout) :: vect(t_dim)
          double precision :: zero = 0d0
#endif
          integer i, info_lapack

          if (sigma.eq.zero) then
#ifdef USE_MULTI
              call solve_multi_domain(sigma,vect(1:a_dim))
#else
              if (grd(1)%mattype.eq.'FULL') then
                  call DGETRV('N', a_dim, asigma(1)%mat, a_dim, asigma(1)%ipiv, &
                      vect(1:a_dim), info_lapack)
              elseif (grd(1)%mattype.eq.'BAND') then
                  call DGBTRS('N', a_dim, kl, ku, 1, asigma(1)%mat, lda, asigma(1)%ipiv, &
                      vect(1:a_dim), a_dim, info_lapack)
              else
                  stop 'mattype has a faulty value in solve_amsigmb'
              endif
              if (info_lapack.ne.0) &
                  stop 'Problem solving linear system, in arncheb'
#endif
          else
#ifdef USE_MPI
              if (iproc.eq.0) then
#endif
                  if (power_max.gt.1) then
                      temp2(1:a_dim,1) = sigma*vect(1+a_dim:2*a_dim)
                  endif
#ifdef USE_MPI
              endif
              do i=2,power_max-1
              temp2(1:a_dim,i) = temp2(1:a_dim,i-1)+&
                  vect(1+i*a_dim:(i+1)*a_dim)
              temp2(1:a_dim,i) = temp2(1:a_dim,i)*sigma
              enddo
              call DGEBS2D(ictxt1D,"A",topo,a_dim,power_max,temp2,a_dim)
          else
              call DGEBR2D(ictxt1D,"A",topo,a_dim,power_max,temp2,a_dim,0,0)
          endif
#endif
          do i=2,power_max-1
          temp2(1:a_dim,i) = temp2(1:a_dim,i-1)+&
              vect(1+i*a_dim:(i+1)*a_dim)
          temp2(1:a_dim,i) = temp2(1:a_dim,i)*sigma
          enddo
          do i=1,power_max-1
#ifdef USE_MULTI
          call a_product_total(temp2(1:a_dim,i),temp(1:a_dim), i+1)
#else
          call a_product(temp2(1:a_dim,i),temp(1:a_dim), i+1)
#endif
#ifdef USE_MPI
          if (iproc.eq.0) then
#endif
              vect(1:a_dim) = vect(1:a_dim)-temp(1:a_dim)
#ifdef USE_MPI
          endif
#endif
          enddo
#ifdef USE_MULTI
          call solve_multi_domain(sigma,vect(1:a_dim))
#else
          if (grd(1)%mattype.eq.'FULL') then
              call DGETRV('N', a_dim, asigma(1)%mat, a_dim, asigma(1)%ipiv, &
                  vect(1:a_dim), info_lapack)
          elseif (grd(1)%mattype.eq.'BAND') then
              call DGBTRS('N', a_dim, kl, ku, 1, asigma(1)%mat, lda, asigma(1)%ipiv, &
                  vect(1:a_dim), a_dim, info_lapack)
          else
              stop 'mattype has a faulty value in solve_amsigmb'
          endif
          if (info_lapack.ne.0) &
              stop 'Problem solving linear system, in arncheb'
          do i=1,power_max-1
          vect(1+i*a_dim:(i+1)*a_dim) = vect(1+i*a_dim:(i+1)*a_dim)+&
              sigma*vect(1+(i-1)*a_dim:i*a_dim)
          enddo
#endif
#ifdef USE_MPI
          if (iproc.eq.0) then
#endif
              do i=1,power_max-1
              vect(1+i*a_dim:(i+1)*a_dim) = vect(1+i*a_dim:(i+1)*a_dim)+&
                  sigma*vect(1+(i-1)*a_dim:i*a_dim)
              enddo
#ifdef USE_MPI
          endif
#endif
      endif
      return
  end subroutine solve_amsigmab

!--------------------------------------------------------------
! This routine does:
!   vect <- inv(A-sigma.B)^H vect
! where (A-sigma.B) corresponds to the total system, taking
! into account the powers of sigma
!--------------------------------------------------------------
! Variables:
!
! sigma = the eigenvalue shift
! vect  = input and output vector
!--------------------------------------------------------------
  subroutine solve_amsigmab_transpose(sigma, vect)

  implicit none
#ifdef USE_COMPLEX
  double complex, intent(in)    :: sigma
  double complex, intent(inout) :: vect(t_dim)
  double complex invsigma
  double complex :: zero = (0d0, 0d0)
#else
  double precision, intent(in)    :: sigma
  double precision, intent(inout) :: vect(t_dim)
  double precision invsigma
  double precision :: zero = 0d0
#endif
  integer i, info_lapack

  if (sigma.eq.zero) then
#ifdef USE_MULTI
  call solve_multi_domain_transpose(sigma,vect(1:a_dim))
#endif
          else
              temp(1:a_dim) = vect(1+(power_max-1)*a_dim:power_max*a_dim)
              do i=power_max-2,0,-1
#ifdef USE_COMPLEX
              temp(1:a_dim) = dconjg(sigma)*temp(1:a_dim)+vect(1+i*a_dim:(i+1)*a_dim)
#else
              temp(1:a_dim) = sigma*temp(1:a_dim)+vect(1+i*a_dim:(i+1)*a_dim)
#endif
              enddo
#ifdef USE_MULTI
              call solve_multi_domain_transpose(sigma,temp(1:a_dim))
#endif
              do i=0,power_max-1
              temp2(1:a_dim,i) = vect(1+i*a_dim:(i+1)*a_dim)
              enddo
#ifdef USE_COMPLEX
              invsigma = (1d0,0d0)/dconjg(sigma)
#else
              invsigma = 1d0/sigma
#endif
              vect(1:a_dim) = temp(1:a_dim)
              if (power_max.gt.1) then
#ifdef USE_MULTI
                  call a_product_total_transpose(vect(1:a_dim),vect(1+a_dim:2*a_dim),0)
                  call a_product_total_transpose(vect(1:a_dim),temp(1:a_dim),1)
#else
                  call a_product_transpose(vect(1:a_dim),vect(1+a_dim:2*a_dim),0)
                  call a_product_transpose(vect(1:a_dim),temp(1:a_dim),1)
#endif
                  vect(1+a_dim:2*a_dim) = vect(1+a_dim:2*a_dim)          &
                      - temp2(1:a_dim,0)
                  vect(1+a_dim:2*a_dim) = invsigma*vect(1+a_dim:2*a_dim) &
                      + temp(1:a_dim)
              endif
              do i=2,power_max-1
#ifdef USE_MULTI
              call a_product_total_transpose(vect(1:a_dim),temp(1:a_dim),i)
#else
              call a_product_transpose(vect(1:a_dim),temp(1:a_dim),i)
#endif
              vect(1+i*a_dim:(i+1)*a_dim)=vect(1+(i-1)*a_dim:i*a_dim)   &
                  -temp2(1:a_dim,i-1)
              vect(1+i*a_dim:(i+1)*a_dim)=invsigma*vect(1+i*a_dim:(i+1)*a_dim) &
                  +temp(1:a_dim)
              enddo
          endif
          return
      end subroutine solve_amsigmab_transpose

!--------------------------------------------------------------
! Correction to matrices before factorisation, in a
! multi-domain context.
!--------------------------------------------------------------
! Variables:
!
! id    = domain number
! sigma = eigenvalue shift
!--------------------------------------------------------------
#ifdef USE_MULTI
      subroutine correct_asigma_downward(id,sigma)

          implicit none
          integer, intent(in)         :: id
#ifdef USE_COMPLEX
          double complex, intent(in)  :: sigma
          double complex, allocatable :: aiip(:,:)
#else
          double precision, intent(in)  :: sigma
          double precision, allocatable :: aiip(:,:)
#endif
          integer info_lapack, ku, kl, lda, d_dim, n_h_bc

          ! Just in case:
          if (id.eq.1) return

          ! Easy exit if possible
          if (idm(id-1,id)%n_v_bc == 0) return
          if (idm(id,id-1)%n_v_bc == 0) return

          d_dim = dm(id-1)%d_dim
          n_h_bc = idm(id-1,id)%n_h_bc

          allocate(aiip(d_dim,n_h_bc))
          call make_asigma_full_local_hcomp(sigma,aiip,id-1,id)
          if (grd(id-1)%mattype.eq.'FULL') then
#ifdef USE_COMPLEX
              call ZGETRS('N',d_dim,n_h_bc,asigma(id-1)%mat,d_dim,     &
                  asigma(id-1)%ipiv,aiip,d_dim,info_lapack)
#else
              call DGETRS('N',d_dim,n_h_bc,asigma(id-1)%mat,d_dim,     &
                  asigma(id-1)%ipiv,aiip,d_dim,info_lapack)
#endif
          else
              ku = dm(id-1)%ku
              kl = dm(id-1)%kl
              lda = 2*kl + ku + 1
#ifdef USE_COMPLEX
              call ZGBTRS('N',d_dim,kl,ku,n_h_bc,asigma(id-1)%mat,lda, &
                  asigma(id-1)%ipiv,aiip,d_dim,info_lapack)
#else
              call DGBTRS('N',d_dim,kl,ku,n_h_bc,asigma(id-1)%mat,lda, &
                  asigma(id-1)%ipiv,aiip,d_dim,info_lapack)
#endif
          endif

          if (grd(id)%mattype.eq.'FULL') then
              call asigma_full_product_subtract_local_hcomp(aiip, &
                  asigma(id)%mat,id,sigma)
          else
              call asigma_band_product_subtract_local_hcomp(aiip, &
                  asigma(id)%mat,id,sigma)
          endif
          deallocate(aiip)
      end subroutine correct_asigma_downward
#endif

!--------------------------------------------------------------
! Solve the system in a multi-domain context:
!    vect <- inv(asigma) vect
! where asigma = a(0) + sigma.a(1) + sigma^2.a(2) ...
!--------------------------------------------------------------
! Variables:
!
! sigma = the eigenvalue shift
! vect  = input and output vector
!--------------------------------------------------------------
#ifdef USE_MULTI
      subroutine solve_multi_domain(sigma,vect)

          implicit none
#ifdef USE_COMPLEX
          double complex, intent(in)    :: sigma
          double complex, intent(inout) :: vect(a_dim)
#else
          double precision, intent(in)    :: sigma
          double precision, intent(inout) :: vect(a_dim)
#endif
          integer i, id, d_dim, offset, info_lapack, ku, kl, lda
          integer d_dimn, offsetn

          ! Downward sweep
          d_dimn = dm(1)%d_dim
          offsetn= dm(1)%offset
          do id=2,ndomains
          d_dim = d_dimn
          offset = offsetn
          d_dimn = dm(id)%d_dim
          offsetn= dm(id)%offset

          ! this is to avoid unecessary calculations
          if (idm(id,id-1)%n_v_bc.eq.0) cycle

          ! The recursive formula uses tilde(y) rather than y.
          ! Therefore, tmp_s can only be set right before the
          ! calculations.
          tmp_s(1:d_dim) = vect((offset+1):(offset+d_dim))

          if (grd(id-1)%mattype.eq.'FULL') then
#ifdef USE_COMPLEX
              call ZGETRV('N',d_dim,asigma(id-1)%mat,d_dim,             &
                  asigma(id-1)%ipiv,tmp_s(1:d_dim),info_lapack)
#else
              call DGETRV('N',d_dim,asigma(id-1)%mat,d_dim,             &
                  asigma(id-1)%ipiv,tmp_s(1:d_dim),info_lapack)
#endif
          else
              ku = dm(id-1)%ku
              kl = dm(id-1)%kl
              lda = 2*kl + ku + 1
              call ZGBTRS('N',d_dim,kl,ku,1,asigma(id-1)%mat,lda,        &
                  asigma(id-1)%ipiv,tmp_s(1:d_dim),              &
                  d_dim,info_lapack)
          endif
          call asigma_product_subtract_local(tmp_s(1:d_dim),           &
              vect((offsetn+1):(offsetn+d_dimn)), id, id-1, sigma)
          enddo

          ! Upward sweep:
          if (grd(ndomains)%mattype.eq.'FULL') then
#ifdef USE_COMPLEX
              call ZGETRV('N',d_dimn,asigma(ndomains)%mat,d_dimn,   &
                  asigma(ndomains)%ipiv,                    &
                  vect((offsetn+1):(offsetn+d_dimn)),       &
                  info_lapack)
#else
              call DGETRV('N',d_dimn,asigma(ndomains)%mat,d_dimn,   &
                  asigma(ndomains)%ipiv,                    &
                  vect((offsetn+1):(offsetn+d_dimn)),       &
                  info_lapack)
#endif
          else
              ku = dm(ndomains)%ku
              kl = dm(ndomains)%kl
              lda = 2*kl + ku + 1
#ifdef USE_COMPLEX
              call ZGBTRS('N',d_dimn,kl,ku,1,asigma(ndomains)%mat,lda, &
                  asigma(ndomains)%ipiv,                       &
                  vect((offsetn+1):(offsetn+d_dimn)),          &
                  d_dimn,info_lapack)
#else
              call DGBTRS('N',d_dimn,kl,ku,1,asigma(ndomains)%mat,lda, &
                  asigma(ndomains)%ipiv,                       &
                  vect((offsetn+1):(offsetn+d_dimn)),          &
                  d_dimn,info_lapack)
#endif
          endif

          do id=ndomains-1,1,-1
          d_dim = d_dimn
          offset = offsetn
          d_dimn = dm(id)%d_dim
          offsetn= dm(id)%offset
          call asigma_product_subtract_local(vect((offset+1):(offset+d_dim)),&
              vect((offsetn+1):(offsetn+d_dimn)), id, id+1, sigma)
          if (grd(id)%mattype.eq.'FULL') then
#ifdef USE_COMPLEX
              call ZGETRV('N',d_dimn,asigma(id)%mat,d_dimn,   &
                  asigma(id)%ipiv,                    &
                  vect((offsetn+1):(offsetn+d_dimn)), &
                  info_lapack)
#else
              call DGETRV('N',d_dimn,asigma(id)%mat,d_dimn,   &
                  asigma(id)%ipiv,                    &
                  vect((offsetn+1):(offsetn+d_dimn)), &
                  info_lapack)
#endif
          else
              ku = dm(id)%ku
              kl = dm(id)%kl
              lda = 2*kl + ku + 1
#ifdef USE_COMPLEX
              call ZGBTRS('N',d_dimn,kl,ku,1,asigma(id)%mat,lda, &
                  asigma(id)%ipiv,                       &
                  vect((offsetn+1):(offsetn+d_dimn)),    &
                  d_dimn,info_lapack)
#else
              call DGBTRS('N',d_dimn,kl,ku,1,asigma(id)%mat,lda, &
                  asigma(id)%ipiv,                       &
                  vect((offsetn+1):(offsetn+d_dimn)),    &
                  d_dimn,info_lapack)
#endif
          endif
          enddo

      end subroutine solve_multi_domain
#endif

!--------------------------------------------------------------
! Solve the conjugate transposed system in a multi-domain context:
!    vect <- inv(asigma)^H vect
! where asigma = a(0) + sigma.a(1) + sigma^2.a(2) ...
!--------------------------------------------------------------
! Variables:
!
! sigma = the eigenvalue shift
! vect  = input and output vector
!--------------------------------------------------------------
#ifdef USE_MULTI
      subroutine solve_multi_domain_transpose(sigma,vect)

          implicit none
#ifdef USE_COMPLEX
          double complex, intent(in)    :: sigma
          double complex, intent(inout) :: vect(a_dim)
#else
          double precision, intent(in)    :: sigma
          double precision, intent(inout) :: vect(a_dim)
#endif
          integer i, id, d_dim, d_dimn, offset, offsetn, ku, kl, lda
          integer info_lapack

          ! Downward sweep:
          d_dimn = dm(1)%d_dim
          offsetn= dm(1)%offset
          do id=2,ndomains
          d_dim  = d_dimn
          offset = offsetn
          d_dimn = dm(id)%d_dim
          offsetn= dm(id)%offset

          ! This is to avoid unecessary calculations:
          if (idm(id-1,id)%n_v_bc.eq.0) cycle

          ! The recursive formula uses tilde(y) rather than y.
          ! Therefore, tmp_s can only be set right before the
          ! calculations.
          tmp_s(1:d_dim) = vect((offset+1):(offset+d_dim))

          if (grd(id-1)%mattype.eq.'FULL') then
#ifdef USE_COMPLEX
              call ZGETRV('C',d_dim,asigma(id-1)%mat,d_dim,   &
                  asigma(id-1)%ipiv,                  &
                  tmp_s(1:d_dim),info_lapack)
#else
              call DGETRV('T',d_dim,asigma(id-1)%mat,d_dim,   &
                  asigma(id-1)%ipiv,                  &
                  tmp_s(1:d_dim),info_lapack)
#endif
          else
              ku = dm(id-1)%ku
              kl = dm(id-1)%kl
              lda = 2*kl + ku + 1
#ifdef USE_COMPLEX
              call ZGBTRS('C',d_dim,kl,ku,1,asigma(id-1)%mat,lda, &
                  asigma(id-1)%ipiv,                      &
                  tmp_s(1:d_dim),d_dim,info_lapack)
#else
              call DGBTRS('T',d_dim,kl,ku,1,asigma(id-1)%mat,lda, &
                  asigma(id-1)%ipiv,                      &
                  tmp_s(1:d_dim),d_dim,info_lapack)
#endif
          endif
          call asigma_product_subtract_local_transpose(tmp_s(1:d_dim),&
              vect((offsetn+1):(offsetn+d_dimn)), id-1, id, sigma)
          enddo

          ! Upward sweep:
          if (grd(ndomains)%mattype.eq.'FULL') then
#ifdef USE_COMPLEX
              call ZGETRV('C',d_dimn,asigma(ndomains)%mat,d_dimn,   &
                  asigma(ndomains)%ipiv,                    &
                  vect((offsetn+1):(offsetn+d_dimn)),       &
                  info_lapack)
#else
              call DGETRV('T',d_dimn,asigma(ndomains)%mat,d_dimn,   &
                  asigma(ndomains)%ipiv,                    &
                  vect((offsetn+1):(offsetn+d_dimn)),       &
                  info_lapack)
#endif
          else
              ku = dm(ndomains)%ku
              kl = dm(ndomains)%kl
              lda = 2*kl + ku + 1
#ifdef USE_COMPLEX
              call ZGBTRS('C',d_dimn,kl,ku,1,asigma(ndomains)%mat,lda, &
                  asigma(ndomains)%ipiv,                       &
                  vect((offsetn+1):(offsetn+d_dimn)),          &
                  d_dimn,info_lapack)
#else
              call DGBTRS('T',d_dimn,kl,ku,1,asigma(ndomains)%mat,lda, &
                  asigma(ndomains)%ipiv,                       &
                  vect((offsetn+1):(offsetn+d_dimn)),          &
                  d_dimn,info_lapack)
#endif
          endif
          do id=ndomains-1,1,-1
          d_dim = d_dimn
          offset = offsetn
          d_dimn = dm(id)%d_dim
          offsetn= dm(id)%offset
          call asigma_product_subtract_local_transpose(vect((offset+1):(offset+d_dim)),&
              vect((offsetn+1):(offsetn+d_dimn)), id+1, id, sigma)
          if (grd(id)%mattype.eq.'FULL') then
#ifdef USE_COMPLEX
              call ZGETRV('C',d_dimn,asigma(id)%mat,d_dimn,   &
                  asigma(id)%ipiv,                    &
                  vect((offsetn+1):(offsetn+d_dimn)), &
                  info_lapack)
#else
              call DGETRV('T',d_dimn,asigma(id)%mat,d_dimn,   &
                  asigma(id)%ipiv,                    &
                  vect((offsetn+1):(offsetn+d_dimn)), &
                  info_lapack)
#endif
          else
              ku = dm(id)%ku
              kl = dm(id)%kl
              lda = 2*kl + ku + 1
#ifdef USE_COMPLEX
              call ZGBTRS('C',d_dimn,kl,ku,1,asigma(id)%mat,lda, &
                  asigma(id)%ipiv,                       &
                  vect((offsetn+1):(offsetn+d_dimn)),    &
                  d_dimn,info_lapack)
#else
              call DGBTRS('T',d_dimn,kl,ku,1,asigma(id)%mat,lda, &
                  asigma(id)%ipiv,                       &
                  vect((offsetn+1):(offsetn+d_dimn)),    &
                  d_dimn,info_lapack)
#endif
          endif
          enddo

      end subroutine solve_multi_domain_transpose
#endif
!--------------------------------------------------------------
  end module
