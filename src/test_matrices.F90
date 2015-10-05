#include "config.h"
      module eigensolve

      use mod_grid
      use matrices

      type MULTI_MAT
        double complex, pointer :: mat(:,:)
        integer, pointer        :: ipiv(:)
      end type MULTI_MAT

      integer, save :: t_dim, nsol_out
      double complex, save, allocatable :: omega(:), vec(:,:)
      double complex, save, allocatable :: temp0(:),temp1(:), &
                                           temp2(:),tmp_s(:)
      double complex, save, allocatable :: big_mat(:,:)
      type(MULTI_MAT), save, allocatable :: asigma(:)

contains

!--------------------------------------------------------------
! This subroutine is used to test matrices.f90, to make sure
! the different subroutines produce the correct results.
!--------------------------------------------------------------
! Variables:
!
! sigma = shift that would have been used in the Arnoldi-
!         Chebyshev algorithm
!--------------------------------------------------------------

      subroutine test_matrices(sigma)

      implicit none

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Input variables:

      double complex, intent(in) :: sigma

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Variables used in the LU factorization
!

      integer info_lapack
      integer, allocatable :: ipiv(:)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Other variables and functions:
!

      double precision aux, auxr, auxi, max_diff, max_diff_transpose
      integer i, j, k, id, d_dim, var, ku, kl, lda
      double complex, external :: zdotc

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Create workspaces by allocating various arrays:
!
      if (allocated(temp0))   deallocate(temp0)
      if (allocated(temp1))   deallocate(temp1)
      if (allocated(temp2))   deallocate(temp2)
      if (allocated(tmp_s))   deallocate(tmp_s)
      if (allocated(big_mat)) deallocate(big_mat)
      if (allocated(ipiv))    deallocate(ipiv)

      t_dim = a_dim * power_max

!  Extra space needed for the matrix vector products
      allocate(temp0(1:a_dim),temp1(1:a_dim), &
               temp2(1:a_dim),tmp_s(1:a_dim), &
               big_mat(1:a_dim,1:a_dim),ipiv(1:a_dim))

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Let the tests begin ...
!

      print*,"Let the test begin ..."
      temp0 = (0d0,0d0)
      temp1 = (0d0,0d0)
      temp2 = (0d0,0d0)
      max_diff = 0d0
      max_diff_transpose = 0d0
      do i=0,power_max
        call make_a_full_total(i,big_mat)
        do j=1,a_dim
          temp0(j) = (1d0,0d0)
          if (j.gt.1) temp0(j-1) = (0d0,0d0)
          call a_product_total(temp0,temp1,i)
          call ZGEMV('N',a_dim,a_dim,(1d0,0d0),big_mat,a_dim,temp0,1,(0d0,0d0),temp2,1)
          aux = 0d0
          do k = 1,a_dim
            auxr= dreal(temp1(k)) - dreal(temp2(k))
            auxi= dimag(temp1(k)) - dimag(temp2(k))
            aux = aux + auxr*auxr + auxi*auxi
          enddo
          aux = sqrt(aux/dble(a_dim))
          if (aux.gt.max_diff) max_diff = aux
        enddo

        temp0(a_dim) = (0d0,0d0)
        do j=1,a_dim
          temp0(j) = (1d0,0d0)
          if (j.gt.1) temp0(j-1) = 0d0
          call a_product_total_transpose(temp0,temp1,i)
          call ZGEMV('C',a_dim,a_dim,(1d0,0d0),big_mat,a_dim,temp0,1,(0d0,0d0),temp2,1)
          aux = 0d0
          do k = 1,a_dim
            auxr= dreal(temp1(k)) - dreal(temp2(k))
            auxi= dimag(temp1(k)) - dimag(temp2(k))
            aux = aux + auxr*auxr + auxi*auxi
          enddo
          aux = sqrt(aux/dble(a_dim))
          if (aux.gt.max_diff_transpose) max_diff_transpose = aux
        enddo

      enddo
      print*,"*****************************"
      print*,"Matrix products:"
      print*,"Max diff (normal) =",max_diff
      print*,"Max diff (trans.) =",max_diff_transpose

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Factorize the matrix A -shift*B
!

      do id=1,ndomains
        if (grd(id)%mattype.eq.'BAND') then
          call init_ku_kl(id)
          call increase_ku_kl_downward(id)
        elseif (grd(id)%mattype.ne.'FULL') then
          print*,"Domain: ",id
          print*,"Faulty mattype: ",grd(id)%mattype
          stop
        endif
      enddo

      allocate(asigma(ndomains))

      ! downward sweep method:
      do id=1,ndomains
        d_dim = dm(id)%d_dim
        allocate(asigma(id)%ipiv(d_dim))
        if (grd(id)%mattype.eq.'FULL') then
          allocate(asigma(id)%mat(d_dim,d_dim))
          call make_asigma_full_local(sigma,asigma(id)%mat,id,id)
          if (id.gt.1) call correct_asigma_downward(id,sigma)

          call ZGETRF(d_dim,d_dim,asigma(id)%mat,d_dim,asigma(id)%ipiv,info_lapack)
          if (info_lapack.ne.0) then
            print*,info_lapack,' info'
            print*,'Factorisation problem in domain ',id
            stop
          endif
        else
          ku = dm(id)%ku
          kl = dm(id)%kl
          lda = 2*kl + ku + 1
          allocate(asigma(id)%mat(lda,d_dim))
          call make_asigma_band_local(sigma,asigma(id)%mat,id)
          if (id.gt.1) call correct_asigma_downward(id,sigma)

          call ZGBTRF(d_dim,d_dim,kl,ku,asigma(id)%mat,lda,asigma(id)%ipiv,info_lapack)
          if (info_lapack.ne.0) then
            print*,info_lapack,' info'
            print*,'Factorisation problem in domain ',id
            stop
          endif
        endif
      enddo

      call make_asigma_full_total(sigma,big_mat)
      call ZGETRF(a_dim,a_dim,big_mat,a_dim,ipiv,info_lapack)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Further tests ...
!

      max_diff = 0d0
      max_diff_transpose = 0d0
      do i=1,a_dim
        temp1 = (0d0,0d0)
        temp2 = (0d0,0d0)
        temp1(i) = (1d0,0d0)
        temp2(i) = (1d0,0d0)

        call ZGETRS('N',a_dim,1,big_mat,a_dim,     &
                    ipiv,temp1,a_dim,info_lapack)
        call solve_multi_domain(sigma,temp2)

        aux = 0d0
        do j=1,a_dim
          auxr= dreal(temp1(k)) - dreal(temp2(k))
          auxi= dimag(temp1(k)) - dimag(temp2(k))
          aux = aux + auxr*auxr + auxi*auxi
        enddo
        aux = sqrt(aux/dble(a_dim))
        if (aux.gt.max_diff) max_diff = aux
      enddo

      do i=1,a_dim
        temp1 = (0d0,0d0)
        temp2 = (0d0,0d0)
        temp1(i) = (1d0,0d0)
        temp2(i) = (1d0,0d0)

        call ZGETRS('C',a_dim,1,big_mat,a_dim,     &
                    ipiv,temp1,a_dim,info_lapack)
        call solve_multi_domain_transpose(sigma,temp2)

        aux = 0d0
        do j=1,a_dim
          auxr= dreal(temp1(k)) - dreal(temp2(k))
          auxi= dimag(temp1(k)) - dimag(temp2(k))
          aux = aux + auxr*auxr + auxi*auxi
        enddo
        aux = sqrt(aux/dble(a_dim))
        if (aux.gt.max_diff_transpose) max_diff_transpose = aux
      enddo

      print*,"*****************************"
      print*,"Matrix inversions:"
      print*,"Max diff (normal) =",max_diff
      print*,"Max diff (trans.) =",max_diff_transpose

      end subroutine test_matrices

!--------------------------------------------------------------
! Correction to matrices before factorisation, in a
! multi-domain context.
!--------------------------------------------------------------
! Variables:
!
! id    = domain number
! sigma = eigenvalue shift
!--------------------------------------------------------------
      subroutine correct_asigma_downward(id,sigma)

      implicit none
      integer, intent(in)         :: id
      double complex, intent(in)  :: sigma
      double complex, allocatable :: aiip(:,:)
      integer info_lapack, ku, kl, lda, d_dim, n_h_bc

      ! Just in case:
      if (id.eq.1) return

      ! Easy exit if possible
      if (idm(id-1,id)%n_v_bc.eq.0) return
      if (idm(id,id-1)%n_v_bc.eq.0) return

      d_dim = dm(id-1)%d_dim
      n_h_bc = idm(id-1,id)%n_h_bc

      allocate(aiip(d_dim,n_h_bc))
      call make_asigma_full_local_hcomp(sigma,aiip,id-1,id)
      if (grd(id-1)%mattype.eq.'FULL') then
        call ZGETRS('N',d_dim,n_h_bc,asigma(id-1)%mat,d_dim,     &
                    asigma(id-1)%ipiv,aiip,d_dim,info_lapack)
      else
        ku = dm(id-1)%ku
        kl = dm(id-1)%kl
        lda = 2*kl + ku + 1
        call ZGBTRS('N',d_dim,kl,ku,n_h_bc,asigma(id-1)%mat,lda, &
                    asigma(id-1)%ipiv,aiip,d_dim,info_lapack)
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
      subroutine solve_multi_domain(sigma,vect)

      implicit none
      double complex, intent(in)    :: sigma
      double complex, intent(inout) :: vect(a_dim)
      integer i, id, d_dim, offset, info_lapack, ku, kl, lda
      integer d_dimn, offsetn

      ! make a copy of vect:
      do i = 1, a_dim
        tmp_s(i) = vect(i)
      enddo

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

        if (grd(id-1)%mattype.eq.'FULL') then
          call ZGETRS('N',d_dim,1,asigma(id-1)%mat,d_dim,                &
                      asigma(id-1)%ipiv,tmp_s((offset+1):(offset+d_dim)),&
                      d_dim,info_lapack)
        else
          ku = dm(id-1)%ku
          kl = dm(id-1)%kl
          lda = 2*kl + ku + 1
          call ZGBTRS('N',d_dim,kl,ku,1,asigma(id-1)%mat,lda,            &
                      asigma(id-1)%ipiv,tmp_s((offset+1):(offset+d_dim)),&
                      d_dim,info_lapack)
        endif
        call asigma_product_subtract_local(tmp_s((offset+1):(offset+d_dim)),&
               vect((offsetn+1):(offsetn+d_dimn)), id, id-1, sigma)
      enddo

      ! Upward sweep:
      if (grd(ndomains)%mattype.eq.'FULL') then
        call ZGETRS('N',d_dimn,1,asigma(ndomains)%mat,d_dimn, &
                    asigma(ndomains)%ipiv,                    &
                    vect((offsetn+1):(offsetn+d_dimn)),       &
                    d_dimn,info_lapack)
      else
        ku = dm(ndomains)%ku
        kl = dm(ndomains)%kl
        lda = 2*kl + ku + 1
        call ZGBTRS('N',d_dimn,kl,ku,1,asigma(ndomains)%mat,lda, &
                    asigma(ndomains)%ipiv,                       &
                    vect((offsetn+1):(offsetn+d_dimn)),          &
                    d_dimn,info_lapack)
      endif

      do id=ndomains-1,1,-1
        d_dim = d_dimn
        offset = offsetn
        d_dimn = dm(id)%d_dim
        offsetn= dm(id)%offset
        call asigma_product_subtract_local(vect((offset+1):(offset+d_dim)),&
               vect((offsetn+1):(offsetn+d_dimn)), id, id+1, sigma)
        if (grd(id)%mattype.eq.'FULL') then
          call ZGETRS('N',d_dimn,1,asigma(id)%mat,d_dimn, &
                      asigma(id)%ipiv,                    &
                      vect((offsetn+1):(offsetn+d_dimn)), &
                      d_dimn,info_lapack)
        else
          ku = dm(id)%ku
          kl = dm(id)%kl
          lda = 2*kl + ku + 1
          call ZGBTRS('N',d_dimn,kl,ku,1,asigma(id)%mat,lda, &
                      asigma(id)%ipiv,                       &
                      vect((offsetn+1):(offsetn+d_dimn)),    &
                      d_dimn,info_lapack)
        endif
      enddo

      end subroutine solve_multi_domain

!--------------------------------------------------------------
! Solve the transposed system in a multi-domain context:
!    vect <- inv(asigma)^T vect
! where asigma = a(0) + sigma.a(1) + sigma^2.a(2) ...
!--------------------------------------------------------------
! Variables:
!
! sigma = the eigenvalue shift
! vect  = input and output vector
!--------------------------------------------------------------
      subroutine solve_multi_domain_transpose(sigma,vect)

      implicit none
      double complex, intent(in)    :: sigma
      double complex, intent(inout) :: vect(a_dim)
      integer i, id, d_dim, d_dimn, offset, offsetn, ku, kl, lda
      integer info_lapack

      ! make a copy of vect:
      do i = 1, a_dim
        tmp_s(i) = vect(i)
      enddo

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

        if (grd(id-1)%mattype.eq.'FULL') then
          call ZGETRS('C',d_dim,1,asigma(id-1)%mat,d_dim, &
                      asigma(id-1)%ipiv,                  &
                      tmp_s((offset+1):(offset+d_dim)),   &
                      d_dim,info_lapack)
        else
          ku = dm(id-1)%ku
          kl = dm(id-1)%kl
          lda = 2*kl + ku + 1
          call ZGBTRS('C',d_dim,kl,ku,1,asigma(id-1)%mat,lda, &
                      asigma(id-1)%ipiv,                      &
                      vect((offset+1):(offset+d_dim)),        &
                      d_dim,info_lapack)
        endif
        call asigma_product_subtract_local_transpose(tmp_s((offset+1):(offset+d_dim)),&
               vect((offsetn+1):(offsetn+d_dimn)), id-1, id, sigma)
      enddo

      ! Upward sweep:
      if (grd(ndomains)%mattype.eq.'FULL') then
        call ZGETRS('C',d_dimn,1,asigma(ndomains)%mat,d_dimn, &
                    asigma(ndomains)%ipiv,                    &
                    vect((offsetn+1):(offsetn+d_dimn)),       &
                    d_dimn,info_lapack)
      else
        ku = dm(ndomains)%ku
        kl = dm(ndomains)%kl
        lda = 2*kl + ku + 1
        call ZGBTRS('C',d_dimn,kl,ku,1,asigma(ndomains)%mat,lda, &
                    asigma(ndomains)%ipiv,                       &
                    vect((offsetn+1):(offsetn+d_dimn)),          &
                    d_dimn,info_lapack)
      endif
      do id=ndomains-1,1,-1
        d_dim = d_dimn
        offset = offsetn
        d_dimn = dm(id)%d_dim
        offsetn= dm(id)%offset
        call asigma_product_subtract_local_transpose(vect((offset+1):(offset+d_dim)),&
               vect((offsetn+1):(offsetn+d_dimn)), id+1, id, sigma)
        if (grd(id)%mattype.eq.'FULL') then
          call ZGETRS('C',d_dimn,1,asigma(id)%mat,d_dimn, &
                      asigma(id)%ipiv,                    &
                      vect((offsetn+1):(offsetn+d_dimn)), &
                      d_dimn,info_lapack)
        else
          ku = dm(id)%ku
          kl = dm(id)%kl
          lda = 2*kl + ku + 1
          call ZGBTRS('C',d_dimn,kl,ku,1,asigma(id)%mat,lda, &
                      asigma(id)%ipiv,                       &
                      vect((offsetn+1):(offsetn+d_dimn)),    &
                      d_dimn,info_lapack)
        endif
      enddo

      end subroutine solve_multi_domain_transpose
!--------------------------------------------------------------
      end module
