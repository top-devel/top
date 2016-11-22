module cfg

    implicit none

    logical, save :: dump_asigma
    logical, save :: stop_after_dump
    logical, save :: dump_terms
    character(256), save :: tag

contains

    subroutine set_tag(thetag)
        character(len=*), intent(in) :: thetag

        tag = thetag
    end subroutine set_tag

    subroutine init_config()
        dump_asigma = .false.
        dump_terms = .false.
        stop_after_dump = .false.
        tag = ""
    end subroutine init_config

end module cfg
