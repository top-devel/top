#include "config.h"

module string
    use iso_c_binding

    interface cpy_str
        module procedure c_f_string, f_c_string
    end interface

contains
    function c_f_string(c_str) result(f_str)
        implicit none
        character(kind=c_char) :: c_str(:)
        character(kind=c_char, len=size(c_str)) :: f_str

        integer i

        do i=1, size(c_str)
            f_str(i:i) = c_str(i)
        end do
    end function

    function f_c_string(f_str) result(c_str)
        implicit none
        character(kind=c_char, len=*) :: f_str
        character(kind=c_char) :: c_str(len(f_str))

        integer i

        do i=1, len(f_str)
            c_str(i) = f_str(i:i)
        end do
    end function
end module string
