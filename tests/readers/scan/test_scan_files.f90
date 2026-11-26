program test_scan_files
    use, intrinsic :: iso_fortran_env, only: real64
    use simulation_state
    use data_parser
    use prescan_files
    use input_parser
    implicit none

    character(len=LENPATH), parameter :: case_dir = "../../../mc-topology/testcase-energy/ZIF8-H2O/"

    ! Setup paths and read the system (no reservoir)
    path%input      = trim(case_dir) // "input.maniac"
    path%topology   = trim(case_dir) // "topology.data"
    path%parameters = trim(case_dir) // "parameters.inc"
    path%reservoir  = ""
    status%reservoir_provided = .false.

    call prescan_inputs()
    call read_input_file()
    call read_system_data()

    ! Expected values from the LAMMPS data file
    call assert_bounds(primary%cell%bounds)
    call assert_int("atoms",     primary%num%atoms,     2220)
    call assert_int("bonds",     primary%num%bonds,     9)
    call assert_int("angles",    primary%num%angles,    9)
    call assert_int("dihedrals", primary%num%dihedrals, 0)
    call assert_int("impropers", primary%num%impropers, 0)

    print *, "✅ all scan_files tests passed"

contains

    subroutine assert_int(name, actual, expected)
        character(len=*), intent(in) :: name
        integer, intent(in) :: actual, expected
        if (actual /= expected) then
            print *, "FAIL ", name, ": got ", actual, ", expected ", expected
            error stop 1
        end if
    end subroutine assert_int

    subroutine assert_bounds(bounds)
        real(real64), intent(in) :: bounds(3,2)
        real(real64), parameter :: tol = 1.0d-6, L = 17.01162_real64
        if (any(abs(bounds(:,1) + L) > tol) .or. any(abs(bounds(:,2) - L) > tol)) then
            print *, "FAIL bounds: got ", bounds
            error stop 1
        end if
    end subroutine assert_bounds

end program test_scan_files