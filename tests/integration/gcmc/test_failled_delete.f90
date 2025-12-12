program test_delete_and_create
    use, intrinsic :: iso_fortran_env, only: real64

    ! Use the same modules your codebase uses for initialization & energy
    use prescan_files
    use input_parser
    use data_parser
    use parameters_parser
    use prepare_utils
    use energy_utils          ! <-- provides the module-level `energy` variable
    use simulation_state      ! primary, thermo, status, etc.

    use molecule_creation
    use molecule_deletion
    use monte_carlo_utils
    use random_utils
    use output_utils

    implicit none

    character(len=LENPATH) :: base                  ! Base path to the test input folder containing input files
    real(real64) :: e_total, e_after_delete, e_after_create  ! Total system energies at different stages: initial, after deletion, and after creation
    integer :: res_type, mol_index                  ! Residue type and molecule index used for deletion and creation moves
    logical :: passed                               ! Flag indicating whether the test passed
    real(real64) :: e_coul                          ! Coulombic contribution to the system energy
    real(real64) :: e_long                          ! Long-range contribution to the system energy (Ewald reciprocal + self)
    real(real64), parameter :: ref_total = -0.012882149     ! Reference total energy from LAMMPS for validation
    real(real64), parameter :: ref_ecoul = 28.911538        ! Reference Coulomb energy from LAMMPS for validation
    real(real64), parameter :: ref_elong = -28.924421      ! Reference long-range (Ewald) energy from LAMMPS for validation
    real(real64), parameter :: tol = 0.1                    ! Tolerance used for comparing energy values in the test
    integer :: i

    ! Path to folder that contains input.maniac, topology.data, parameters.inc
    path%input = "input.maniac"
    path%topology = "topology.data"
    path%parameters = "parameters.inc"
    path%reservoir = ""

    status%reservoir_provided = .false.

    ! ---------- Load & initialize system ----------
    call prescan_inputs()
    call read_input_file()
    call read_system_data()
    call read_parameters()
    call setup_simulation_parameters()

    ! ---------- Compute initial energy (module-level `energy`) ----------
    call compute_system_energy(primary)

    e_total = energy%total
    e_coul = energy%coulomb + energy%intra_coulomb
    e_long = energy%recip_coulomb + energy%ewald_self

    if (abs(e_coul - ref_ecoul) > tol) then
        print *, "FAIL: Coulomb energy mismatch"
        stop 1
    end if

    if (abs(e_long - ref_elong) > tol) then
        print *, "FAIL: Long-range energy mismatch"
        stop 1
    end if

    if (abs(e_total - ref_total) > tol) then
        print *, "FAIL: total energy mismatch"
        stop 1
    end if

    ! ---------- Identify the molecule to delete ----------
    res_type = pick_random_residue_type(thermo%is_active)
    mol_index = 1

    if (primary%num%residues(res_type) /= 1) then
        print *, "ERROR: topology must contain exactly 1 molecule"
        stop 1
    end if

    ! ---------- Attempt deleting the molecule ----------
    thermo%chemical_potential = -0.1 ! Probability of deleting the molecule is 0.0

    do i = 1, 10
        call attempt_deletion_move(res_type, mol_index)
    end do

    if (primary%num%residues(res_type) /= 1) then
        print *, "ERROR: Molecule was deleted"
        stop 1
    end if

    ! Energy after deletion
    e_after_delete = energy%total

    if (abs(e_after_delete - ref_total) > tol) then
        print *, "FAIL: Coulomb energy mismatch"
        stop 1
    end if

    ! Full recomputed energy
    call compute_system_energy(primary)
    e_after_delete = energy%total

    if (abs(e_after_delete - ref_total) > tol) then
        print *, "FAIL: Coulomb energy mismatch"
        stop 1
    end if

end program test_delete_and_create
