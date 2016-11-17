module abstract_model_mod
    implicit none

!------------------------------------------------------------------------
    type, abstract :: abstract_model
        private
        character(len=64) :: model_name
        contains
            procedure(init_interface), deferred :: init
            procedure(get_field_interface), deferred :: get_field
            procedure(get_grid_interface), deferred :: get_grid
            procedure(get_grid_size_interface), deferred :: get_grid_size
    end type abstract_model
!------------------------------------------------------------------------

    class(abstract_model), pointer :: model_ptr

    abstract interface
!------------------------------------------------------------------------
    subroutine init_interface(this, filename, ierr)
        import abstract_model
        class(abstract_model), target :: this
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ierr
    end subroutine init_interface
!------------------------------------------------------------------------
    subroutine get_field_interface(this, fname, field)
        import abstract_model
        class(abstract_model) :: this
        character(len=*), intent(in) :: fname
        real(kind=8), allocatable, intent(out) :: field(:, :)
    end subroutine get_field_interface
!------------------------------------------------------------------------
    subroutine get_grid_interface(this, r, th, nr, nt)
        import abstract_model
        class(abstract_model) :: this
        real(kind=8), intent(out) :: r(nr, nt), th(nt)
        integer, intent(in) :: nr, nt
    end subroutine get_grid_interface
!------------------------------------------------------------------------
    subroutine get_grid_size_interface(this, n_r, n_t)
        import abstract_model
        class(abstract_model) :: this
        integer, intent(out) :: n_r, n_t
    end subroutine get_grid_size_interface
!------------------------------------------------------------------------
    end interface

end module abstract_model_mod
