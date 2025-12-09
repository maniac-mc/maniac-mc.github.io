program test_scan_files

    use, intrinsic :: iso_fortran_env, only: real64
    use simulation_state
    use data_parser
    use prescan_files
    use input_parser

    implicit none

    character(len=LENPATH) :: base_main     ! Path for main file
    character(len=LENPATH) :: base_res      ! Path for reservoir file

    !---------------------------------------------------------------------------
    ! Test 1 : ZIF8-H2O without reservoir
    !---------------------------------------------------------------------------
    base_main = "../../../mc-topology/testcase-energy/ZIF8-H2O/"
    base_res = ""
    call run_test_zif8_h2o(base_main, base_res, .false.)

contains

    subroutine run_test_zif8_h2o(base_main, base_res, reservoir)

        ! Inputs
        character(len=*), intent(in) :: base_main
        character(len=*), intent(in) :: base_res
        logical, intent(in) :: reservoir
        integer :: max_atoms_per_residue

        ! Local variables
        logical :: pass, atoms_ok, connectivity_ok
        integer :: i, j, k
        character(len=200) :: msg
        real(real64), dimension(3,2) :: expected_bounds
        integer :: expected_atoms, expected_bonds, expected_angles
        integer :: expected_dihedrals, expected_impropers

        ! Set expected values
        expected_bounds(1,1) = -17.01162d0  ! xlo
        expected_bounds(1,2) =  17.01162d0  ! xhi
        expected_bounds(2,1) = -17.01162d0  ! ylo
        expected_bounds(2,2) =  17.01162d0  ! yhi
        expected_bounds(3,1) = -17.01162d0  ! zlo
        expected_bounds(3,2) =  17.01162d0  ! zhi
        expected_atoms = 2220
        expected_bonds = 9
        expected_angles = 9
        expected_dihedrals = 0
        expected_impropers = 0

        !-------------------------------------------------------------------
        ! Setup paths
        !-------------------------------------------------------------------
        path%input = trim(base_main) // "input.maniac"
        path%topology = trim(base_main) // "topology.data"
        path%parameters = trim(base_main) // "parameters.inc"
        status%reservoir_provided = reservoir
        if (reservoir) then
            path%reservoir = trim(base_res) // "topology.data"
        else
            path%reservoir = trim(base_res)
        end if

        !-------------------------------------------------------------------
        ! Read the system
        !-------------------------------------------------------------------
        call prescan_inputs()
        call read_input_file()
        call read_system_data()

        ! Test Expected values from the LAMMPS data file
        call assert_real_matrix_equal(primary%cell%bounds, expected_bounds, 1.0d-6)
        call assert_int_equal(primary%num%atoms, expected_atoms, "Number of atoms")
        call assert_int_equal(primary%num%bonds, expected_bonds, "Number of bonds")
        call assert_int_equal(primary%num%angles, expected_angles, "Number of angles")
        call assert_int_equal(primary%num%dihedrals, expected_dihedrals, "Number of dihedrals")
        call assert_int_equal(primary%num%impropers, expected_impropers, "Number of impropers")

    end subroutine run_test_zif8_h2o

    subroutine assert_real_matrix_equal(a, b, tol)

        real(real64), intent(in) :: a(:,:), b(:,:)
        real(real64), intent(in) :: tol
        integer :: i, j

        do i = 1, size(a,1)
            do j = 1, size(a,2)
                if (abs(a(i,j) - b(i,j)) > tol) then
                    print *, "ASSERTION FAILED at (", i, ",", j, ")"
                    print *, "  got     =", a(i,j)
                    print *, "  expected=", b(i,j)
                    call exit(1)
                end if
            end do
        end do
    end subroutine assert_real_matrix_equal

    subroutine assert_int_equal(actual, expected, msg)

        integer, intent(in) :: actual, expected
        character(len=*), intent(in) :: msg

        if (actual /= expected) then
            print *, "ASSERTION FAILED:", trim(msg)
            print *, "  got     =", actual
            print *, "  expected=", expected
            call exit(1)
        end if
    end subroutine assert_int_equal

end program test_scan_files
