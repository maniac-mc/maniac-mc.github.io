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

    character(len=LENPATH) :: base
    real(real64) :: e_total, e_after_delete, e_after_create
    integer :: res_type, mol_index
    logical :: passed
    real(real64) :: e_coul     ! Coulomb contribution
    real(real64) :: e_long     ! Long range contribution
    real(real64), parameter :: ref_total = -0.012882149     ! LAMMPS reference value
    real(real64), parameter :: ref_ecoul =  28.911538     ! LAMMPS reference value
    real(real64), parameter :: ref_elong = -28.924421     ! LAMMPS reference value
    real(real64), parameter :: tol = 0.1

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

    ! ---------- Delete the molecule ----------
    thermo%chemical_potential = -10 ! Probability of deleting the molecule is 1.0
    call attempt_deletion_move(res_type, mol_index)

    ! Recompute energy
    call compute_system_energy(primary)
    e_after_delete = energy%total

    if (abs(e_after_delete - 0.0) > tol) then
        print *, "FAIL: total energy mismatch"
        stop 1
    end if

    ! ---------- Recreate the molecule ----------
    mol_index = primary%num%residues(res_type) + 1   ! if deletion removed it, create at new index
    thermo%chemical_potential = -0.1 ! Probability of creating the molecule is 1.0
    call attempt_creation_move(res_type, mol_index)

    ! Recompute energy after creation
    call compute_system_energy(primary)
    e_after_create = energy%total
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

    if (abs(e_after_create - ref_total) > tol) then
        print *, "FAIL: total energy mismatch"
        stop 1
    end if

end program test_delete_and_create
