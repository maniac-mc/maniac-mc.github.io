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
    use molecule_translation
    use molecule_rotation
    use monte_carlo_utils
    use random_utils
    use output_utils
    use write_utils

    implicit none

    character(len=LENPATH) :: base                  ! Base path to the test input folder containing input files
    real(real64) :: e_total, e_recip_after_move     ! Total system energies at different stages: initial, after deletion, and after creation
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
    path%outputs = "outputs/"

    ! ---------- Load & initialize system ----------
    call prescan_inputs()
    call read_input_file()
    call read_system_data()
    call read_parameters()
    call setup_simulation_parameters()

    ! ---------- Compute initial energy (module-level `energy`) ----------
    call update_system_energy(primary)

    ! ---------- Identify the molecule to delete ----------
    res_type = pick_random_residue_type(thermo%is_active)

    call update_output_files(.false.)

    ! ---------- ENERGY CONSERVATION DURING TRANSLATION ----------

    ! Attempt 40 molecule creation
    thermo%chemical_potential = -0.1 ! Probability of creating the molecule is 1.0
    do i = 1, 40
        status%desired_block = i
        mol_index = primary%num%residues(res_type) + 1
        call attempt_creation_move(res_type, mol_index)
        call update_output_files(.true.)
    end do

    write (*,*) "Number of residue after creation :", primary%num%residues(res_type)

    ! Attempt 10 molecule deletion
    thermo%chemical_potential = -10.0 ! Probability of deleting the molecule is 1.0
    do i = 41, 50
        status%desired_block = i
        mol_index = pick_random_molecule_index(primary%num%residues(res_type))
        call attempt_deletion_move(res_type, mol_index)
        call update_output_files(.true.)
    end do

    write (*,*) "Number of residue after deletion :", primary%num%residues(res_type)

    ! Attemps 50 translation move
    do i = 51, 100
        status%desired_block = i
        mol_index = pick_random_molecule_index(primary%num%residues(res_type))
        call attempt_translation_move(res_type, mol_index)
        call update_output_files(.true.)
    end do

    ! Attemps 50 rotation move
    do i = 101, 150
        status%desired_block = i
        mol_index = pick_random_molecule_index(primary%num%residues(res_type))
        call attempt_rotation_move(res_type, mol_index)
        call update_output_files(.true.)
    end do

    e_recip_after_move = energy%recip_coulomb

    ! Full recalculation of energy
    call update_system_energy(primary)

    if (abs(energy%recip_coulomb - e_recip_after_move) > tol) then
        write (*,*) "FAIL: energy drift during monte carlo move"
        write (*,*) "energy%recip_coulomb", energy%recip_coulomb
        write (*,*) "e_recip_after_move", e_recip_after_move
        stop 1
    end if

end program test_delete_and_create
