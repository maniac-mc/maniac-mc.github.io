program test_prescan_files

    use, intrinsic :: iso_fortran_env, only: real64

    use prescan_files

    implicit none

    character(len=256) :: topology_path, input_file, data_file, inc_file
    integer :: expected_active, expected_inactive
    logical :: pass
    character(len=LENPATH) :: base

    ! Setup paths
    base = "../../../mc-topology/testcase-energy/ZIF8-H2O/"
    path%input = trim(base) // "input.maniac"
    path%topology = trim(base) // "topology.data"
    path%parameters = trim(base) // "parameters.inc"

    ! Expectation in ZIF8-H2O system
    expected_active = 5000
    expected_inactive = 1

    ! Run prescan
    call prescan_inputs()

    ! Validate results
    pass = (nmax%active_residues == expected_active) .and. &
           (nmax%inactive_residues == expected_inactive)

    if (.not. pass) then
        print *, "❌ Test FAILED!"
        if (nmax%active_residues /= expected_active) then
            print *, "   active_residues expected=", expected_active, " got=", nmax%active_residues
        end if
        if (nmax%inactive_residues /= expected_inactive) then
            print *, "   inactive_residues expected=", expected_inactive, " got=", nmax%inactive_residues
        end if
        stop 1
    end if

end program test_prescan_files
