module tabulated_utils

    use simulation_state
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

contains

    !---------------------------------------------------------------------------
    ! Initialize all precomputed tables for the simulation
    !---------------------------------------------------------------------------
    subroutine precompute_table()

        ! Initialise table for erfc_r_table, r6_table, and r12_table
        call initialize_tabulated_erfc(erfc_r_table, ewald%param%alpha, mc_input%real_space_cutoff)
        call initialize_tabulated_rpower(r6_table, mc_input%real_space_cutoff, 6)
        call initialize_tabulated_rpower(r12_table, mc_input%real_space_cutoff, 12)
    
    end subroutine precompute_table

    !---------------------------------------------------------------------------
    ! Precomputes a tabulated version of the function erfc(alpha*r)/r for faster
    ! evaluation during simulations. Useful for Ewald direct-space Coulomb energy.
    !---------------------------------------------------------------------------
    subroutine initialize_tabulated_erfc(table, alpha, r_cut)

        ! Input variables
        type(tabulated), intent(inout) :: table  ! Tabulated function structure to store x and f arrays
        real(real64), intent(in) :: alpha        ! Ewald screening parameter for erfc(alpha*r)/r
        real(real64), intent(in) :: r_cut        ! Maximum distance (real-space cutoff) for the table

        ! Local variables
        integer :: i                              ! Loop index over table points
        real(real64) :: r                         ! Current distance corresponding to table grid point

        ! Allocate arrays for grid points (x) and function values (f)
        allocate(table%x(0:TABULATED_POINTS))
        allocate(table%f(0:TABULATED_POINTS))
        table%npoint = TABULATED_POINTS
        table%dx = r_cut / real(TABULATED_POINTS, real64)

        ! Fill the table
        do i = 0, TABULATED_POINTS
            r = i * table%dx
            table%x(i) = r
            if (r < error) then
                table%f(i) = two * alpha / sqrt(pi)
            else
                table%f(i) = erfc(alpha*r) / r
            end if
        end do

        ! Mark table as initialized
        table%initialized = .true.
    end subroutine initialize_tabulated_erfc

    !---------------------------------------------------------------------------
    ! Precompute a tabulated version of r^power for faster evaluation
    !---------------------------------------------------------------------------
    subroutine initialize_tabulated_rpower(table, r_cut, power)

        ! Input/Output
        type(tabulated), intent(inout) :: table ! Table structure to store x and f arrays
        real(real64), intent(in) :: r_cut    ! Maximum distance to tabulate
        integer, intent(in) :: power         ! Power to raise r

        ! Local variables
        integer :: i
        real(real64) :: r

        ! Allocate arrays
        allocate(table%x(0:TABULATED_POINTS))
        allocate(table%f(0:TABULATED_POINTS))
        table%npoint = TABULATED_POINTS
        table%dx = r_cut / real(TABULATED_POINTS, real64)

        ! Fill the table
        do i = 0, TABULATED_POINTS
            r = i * table%dx
            table%x(i) = r
            if (r < error) then
                ! Avoid 0^power singularities
                table%f(i) = zero
            else
                table%f(i) = r**power
            end if
        end do

        ! Mark table as initialized
        table%initialized = .true.

    end subroutine initialize_tabulated_rpower

    ! Linear interpolation in a tabulated function
    pure function lookup_tabulated(table, r) result(f_r)

        ! Input variables
        type(tabulated), intent(in) :: table    ! Tabulated function structure to store x and f arrays
        real(real64), intent(in) :: r           ! Query distance

        ! Local variables
        real(real64) :: f_r                     ! Interpolated function value
        real(real64) :: t, f1, f2               ! Interpolation weight and function values
        integer :: i                            ! Index of lower grid point

        if (r <= zero) then
            f_r = table%f(0)
            return
        end if

        if (r >= table%x(table%npoint)) then
            f_r = zero
            return
        end if

        i = int(r / table%dx)
        f1 = table%f(i)
        f2 = table%f(i+1)
        t = (r - table%x(i)) / table%dx
        f_r = (one - t) * f1 + t * f2
        
    end function lookup_tabulated

end module tabulated_utils
