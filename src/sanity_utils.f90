module sanity_utils

    use output_utils

    use iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

contains

    !------------------------------------------------------------
    ! Checks that all elements of a vector are finite (no NaN/Inf)
    !------------------------------------------------------------
    subroutine check_finite_vector(vec, message)

        ! Input parameters
        real(real64), intent(in) :: vec(:)
        character(*), intent(in) :: message

        if (.not. all(ieee_is_finite(vec))) then
            call abort_run(message)
        end if

    end subroutine check_finite_vector

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

    !------------------------------------------------------------
    ! Check if the center of mass is unusually far from all atoms
    ! using coords(3,n_atoms).
    !------------------------------------------------------------
    subroutine check_com_distance(coords, n_atoms, com, threshold, message)
        use iso_fortran_env, only: real64
        implicit none

        ! Inputs
        real(real64), intent(in) :: coords(3, n_atoms)   ! atomic coordinates
        integer,     intent(in) :: n_atoms
        real(real64), intent(in) :: com(3)               ! center of mass
        real(real64), intent(in) :: threshold
        character(*), intent(in) :: message

        ! Local
        real(real64) :: dx, dy, dz, dist
        real(real64) :: min_dist
        integer :: i

        min_dist = huge(1.0_real64)

        do i = 1, n_atoms
            dx = coords(1,i) - com(1)
            dy = coords(2,i) - com(2)
            dz = coords(3,i) - com(3)
            dist = sqrt(dx*dx + dy*dy + dz*dz)
            if (dist < min_dist) min_dist = dist
        end do

        if (min_dist > threshold) then
            call WarnUser(message)
        end if

    end subroutine check_com_distance

end module sanity_utils
