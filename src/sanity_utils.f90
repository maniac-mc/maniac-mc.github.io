module sanity_utils

    use output_utils
    
    use iso_fortran_env, only: real64

    implicit none

contains

    !------------------------------------------------------------
    ! Checks that a position is inside given bounds
    !------------------------------------------------------------
    subroutine check_inside_bounds(pos, lower, upper, message)
        
        ! Input parameters
        real(real64), intent(in) :: pos(:)
        real(real64), intent(in) :: lower(:)
        real(real64), intent(in) :: upper(:)
        character(*), intent(in) :: message

        ! Local variable
        integer :: dim

        do dim = 1, size(pos)
            if (pos(dim) < lower(dim) .or. pos(dim) > upper(dim)) then
                call WarnUser(message)
            end if
        end do

    end subroutine check_inside_bounds

end module sanity_utils
