program test_prescan_files

    use, intrinsic :: iso_fortran_env, only: real64

    use prescan_files

    implicit none

    character(len=256) :: topology_path, input_file, data_file, inc_file
    integer :: expected_active, expected_inactive
    logical :: pass
    character(len=LENPATH) :: base

    !---------------------------------------------------------------------------
    ! ZIF8-H2O system without reservoir
    !---------------------------------------------------------------------------

    ! Setup paths
    base = "../../../mc-topology/testcase-energy/ZIF8-H2O/"
    path%input = trim(base) // "input.maniac"
    path%topology = trim(base) // "topology.data"
    path%parameters = trim(base) // "parameters.inc"

    expected_active = NB_MAX_MOLECULE ! In absence of reservoir, use NB_MAX_MOLECULE
    expected_inactive = 1 ! One host is expected

    ! Run prescan
    call prescan_inputs()

    ! Validate results
    pass = (nmax%active_residues == expected_active) .and. &
           (nmax%inactive_residues == expected_inactive)

    if (.not. pass) then
        if (nmax%active_residues /= expected_active) then
            print *, "   active_residues expected=", expected_active, " got=", nmax%active_residues
        end if
        if (nmax%inactive_residues /= expected_inactive) then
            print *, "   inactive_residues expected=", expected_inactive, " got=", nmax%inactive_residues
        end if
        stop 1
    end if

    !---------------------------------------------------------------------------
    ! ZIF8-CH4O system with reservoir
    !---------------------------------------------------------------------------

    ! Setup paths
    base = "../../../mc-topology/testcase-adsorption/ZIF8-CH4O/"
    path%input = trim(base) // "input.maniac"
    path%topology = trim(base) // "topology.data"
    path%parameters = trim(base) // "parameters.inc"
    base = "../../../mc-topology//molecule-reservoir/CH4O-H2O/"
    path%reservoir = trim(base) // "topology.data"
    status%reservoir_provided = .true.

    expected_active = 5+1500 ! 5 molecules in main, 1500 in reservoir
    expected_inactive = 1 ! One host is expected

    ! Run prescan
    call prescan_inputs()

    ! Validate results
    pass = (nmax%active_residues == expected_active) .and. &
           (nmax%inactive_residues == expected_inactive)

    if (.not. pass) then
        if (nmax%active_residues /= expected_active) then
            print *, "   active_residues expected=", expected_active, " got=", nmax%active_residues
        end if
        if (nmax%inactive_residues /= expected_inactive) then
            print *, "   inactive_residues expected=", expected_inactive, " got=", nmax%inactive_residues
        end if
        stop 1
    end if


end program test_prescan_files
