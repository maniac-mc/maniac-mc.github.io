program test_scan_files

    use, intrinsic :: iso_fortran_env, only: real64
    use data_parser
    use simulation_state

    implicit none

    character(len=LENPATH) :: base_main     ! Path for main file
    character(len=LENPATH) :: base_res      ! Path for reservoir file

    !---------------------------------------------------------------------------
    ! Test 1 : ZIF8-H2O without reservoir
    !---------------------------------------------------------------------------
    base_main = "../../../mc-topology/testcase-energy/ZIF8-H2O/"
    base_res = ""
    call run_test(base_main, base_res, .false., NB_MAX_MOLECULE, 1)

    !---------------------------------------------------------------------------
    ! Test 2 : ZIF8-CH4O system with reservoir
    !---------------------------------------------------------------------------
    base_main = "../../../mc-topology/testcase-adsorption/ZIF8-CH4O/"
    base_res  = "../../../mc-topology/molecule-reservoir/CH4O-H2O/"
    call run_test(base_main, base_res, .true., 5+1500, 1)

    !---------------------------------------------------------------------------
    ! Test 3 : Bulk CO2 system without reservoir
    !---------------------------------------------------------------------------
    base_main = "../../../mc-topology/molecule-reservoir/CO2/"
    base_res = ""
    call run_test(base_main, base_res, .false., NB_MAX_MOLECULE, 0)

    !---------------------------------------------------------------------------
    ! Test 4 : Bulk CO2 system with reservoir
    !---------------------------------------------------------------------------
    base_main = "../../../mc-topology/molecule-reservoir/CO2/"
    call run_test(base_main, base_main, .true., 100, 0)

contains

    subroutine run_test(base_main, base_res, reservoir, exp_active, exp_inactive)
        ! Inputs
        character(len=*), intent(in) :: base_main
        character(len=*), intent(in) :: base_res
        logical, intent(in) :: reservoir
        integer, intent(in) :: exp_active, exp_inactive
        integer :: max_atoms_per_residue

        ! Local variables
        logical :: pass, atoms_ok, connectivity_ok
        integer :: i, j, k
        character(len=200) :: msg

        !-------------------------------------------------------------------
        ! Setup paths
        !-------------------------------------------------------------------
        path%input = trim(base_main) // "input.maniac"
        path%topology = trim(base_main) // "topology.data"
        path%parameters = trim(base_main) // "parameters.inc"
        status%reservoir_provided = reservoir
        if (reservoir) path%reservoir = trim(base_res) // "topology.data"

        !-------------------------------------------------------------------
        ! Allocate arrays inside box
        !-------------------------------------------------------------------
        max_atoms_per_residue = 1000
        allocate(primary%atoms%types(res%number, max_atoms_per_residue))
        allocate(primary%atoms%ids(res%number, max_atoms_per_residue))
        allocate(primary%atoms%names(res%number, max_atoms_per_residue))
        allocate(primary%atoms%charges(res%number, max_atoms_per_residue))

        ! Initialize
        primary%atoms%types = 0
        primary%atoms%ids = 0
        primary%atoms%names = ""
        primary%atoms%charges = 0.0_real64

        !-------------------------------------------------------------------
        ! Read the system
        !-------------------------------------------------------------------
        call read_system_data()

        !-------------------------------------------------------------------
        ! Check residue counts
        !-------------------------------------------------------------------
        pass = (nmax%active_residues == exp_active) .and. &
               (nmax%inactive_residues == exp_inactive)
        if (.not. pass) then
            write(*,*) "Test FAILED: residue count mismatch"
            write(*,*) " Expected active =", exp_active, " got=", nmax%active_residues
            write(*,*) " Expected inactive =", exp_inactive, " got=", nmax%inactive_residues
        end if

        !-------------------------------------------------------------------
        ! Check total atom count matches sum over residues
        !-------------------------------------------------------------------
        atoms_ok = (primary%num%atoms == sum(res%atom))
        if (.not. atoms_ok) then
            write(*,*) "Test FAILED: total atom count mismatch"
            write(*,*) " primary%num%atoms =", primary%num%atoms
            write(*,*) " sum(res%atom) =", sum(res%atom)
        end if

        !-------------------------------------------------------------------
        ! Check connectivity: all bonded atoms exist within their residue
        !-------------------------------------------------------------------
        connectivity_ok = .true.
        do i = 1, res%number
            ! Bonds
            do j = 1, res%bonds(i)
                do k = 1, 2
                    if (.not. isInResidue(primary, i, connect%bonds(i,j,k))) then
                        connectivity_ok = .false.
                        write(msg, '(A,I0,A,I0)') "Bond atom ID ", connect%bonds(i,j,k), &
                            " not found in residue ", i
                        write(*,*) trim(msg)
                    end if
                end do
            end do
            ! Angles
            do j = 1, res%angles(i)
                do k = 1, 3
                    if (.not. isInResidue(primary, i, connect%angles(i,j,k))) then
                        connectivity_ok = .false.
                        write(msg, '(A,I0,A,I0)') "Angle atom ID ", connect%angles(i,j,k), &
                            " not found in residue ", i
                        write(*,*) trim(msg)
                    end if
                end do
            end do
            ! Dihedrals
            do j = 1, res%dihedrals(i)
                do k = 1, 4
                    if (.not. isInResidue(primary, i, connect%dihedrals(i,j,k))) then
                        connectivity_ok = .false.
                        write(msg, '(A,I0,A,I0)') "Dihedral atom ID ", connect%dihedrals(i,j,k), &
                            " not found in residue ", i
                        write(*,*) trim(msg)
                    end if
                end do
            end do
            ! Impropers
            do j = 1, res%impropers(i)
                do k = 1, 4
                    if (.not. isInResidue(primary, i, connect%impropers(i,j,k))) then
                        connectivity_ok = .false.
                        write(msg, '(A,I0,A,I0)') "Improper atom ID ", connect%impropers(i,j,k), &
                            " not found in residue ", i
                        write(*,*) trim(msg)
                    end if
                end do
            end do
        end do

        if (.not. connectivity_ok) then
            write(*,*) "Test FAILED: connectivity errors detected"
        end if

        !-------------------------------------------------------------------
        ! Final verdict
        !-------------------------------------------------------------------
        if (.not. pass .or. .not. atoms_ok .or. .not. connectivity_ok) then
            stop 1
        else
            write(*,*) "Test PASSED"
        end if

    end subroutine run_test

end program test_scan_files
