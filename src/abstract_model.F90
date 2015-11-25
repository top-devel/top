module abstract_model_mod
    implicit none

!------------------------------------------------------------------------
    type, abstract :: abstract_model
        private
        character(len=64) :: model_name
        contains
            procedure(init_interface), deferred :: init
            procedure(get_field_interface), deferred :: get_field
    end type abstract_model
!------------------------------------------------------------------------

    class(abstract_model), pointer :: model_ptr

    abstract interface
!------------------------------------------------------------------------
    subroutine init_interface(this, filename)
        import abstract_model
        class(abstract_model), target :: this
        character(len=*), intent(in) :: filename
    end subroutine init_interface
!------------------------------------------------------------------------
    subroutine get_field_interface(this, fname, field)
        import abstract_model
        class(abstract_model) :: this
        character(len=*), intent(in) :: fname
        real(kind=8), allocatable, intent(out) :: field(:, :)
    end subroutine get_field_interface
!------------------------------------------------------------------------
    end interface

end module abstract_model_mod
