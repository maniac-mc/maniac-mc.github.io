program test_prescan_files

    use, intrinsic :: iso_fortran_env, only: real64
    use data_parser

    implicit none

    character(len=LENPATH) :: base_main     ! Path for main file
    character(len=LENPATH) :: base_res      ! Path for reservoir file

    !---------------------------------------------------------------------------
    ! Test 1 : ZIF8-H2O without reservoir
    !---------------------------------------------------------------------------

    base_main = "../../../mc-topology/testcase-energy/ZIF8-H2O/"
    base_res = ""

    ! Expected : in absence of reservoir, use NB_MAX_MOLECULE
    call run_test(base_main, base_res, .false., NB_MAX_MOLECULE, 1)

    !---------------------------------------------------------------------------
    ! Test 2 : ZIF8-CH4O system with reservoir
    !---------------------------------------------------------------------------

    base_main = "../../../mc-topology/testcase-adsorption/ZIF8-CH4O/"
    base_res  = "../../../mc-topology/molecule-reservoir/CH4O-H2O/"

    ! Expected : 5 molecules in main, 1500 in reservoir
    call run_test(base_main, base_res, .true., 5+1500, 1)

    !---------------------------------------------------------------------------
    ! Test 3 : Bulk CO2 system without reservoir
    !---------------------------------------------------------------------------

    base_main = "../../../mc-topology/molecule-reservoir/CO2/"
    base_res = ""

    ! Expected : in absence of reservoir, use NB_MAX_MOLECULE
    call run_test(base_main, base_res, .false., NB_MAX_MOLECULE, 0)

    !---------------------------------------------------------------------------
    ! Test 4 : Bulk CO2 system with reservoir
    !---------------------------------------------------------------------------

    base_main = "../../../mc-topology/molecule-reservoir/CO2/"

    ! Expected : 50 molecules in main, 50 in reservoir
    call run_test(base_main, base_main, .true., 100, 0)

contains

    subroutine run_test(base_main, base_res, reservoir, exp_active, exp_inactive)

        ! Input parameters
        character(len=*), intent(in) :: base_main   ! Path for main file
        character(len=*), intent(in) :: base_res    ! Path for reservoir file (if reservoir is present)
        logical, intent(in) :: reservoir            ! Is reservoir present
        integer, intent(in) :: exp_active           ! Expected guest residue
        integer, intent(in) :: exp_inactive         ! Expected host residue

        ! Local variable
        logical :: pass                             ! Did the test pass

        ! Set paths
        path%input = trim(base_main) // "input.maniac"
        path%topology = trim(base_main) // "topology.data"
        path%parameters = trim(base_main) // "parameters.inc"

        status%reservoir_provided = reservoir

        if (reservoir) path%reservoir = trim(base_res) // "topology.data"

        ! Call data reader
        call read_system_data()

        ! Check results
        pass = (nmax%active_residues == exp_active) .and. &
            (nmax%inactive_residues == exp_inactive)

        if (.not. pass) then
            print *, "Test FAILED"
            print *, " expected active =", exp_active, " got=", nmax%active_residues
            print *, " expected inactive =", exp_inactive, " got=", nmax%inactive_residues
            stop 1
        end if

    end subroutine run_test

end program test_prescan_files
