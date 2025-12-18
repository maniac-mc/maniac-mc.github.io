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
    path%outputs = "outputs/"

    ! ---------- Load & initialize system ----------
    call prescan_inputs()
    call read_input_file()
    call read_system_data()
    call read_parameters()
    call setup_simulation_parameters()

    ! ---------- Compute initial energy (module-level `energy`) ----------
    call update_system_energy(primary)
    ! e_total = energy%total
    ! e_coul = energy%coulomb + energy%intra_coulomb
    ! e_long = energy%recip_coulomb + energy%ewald_self

    ! ---------- Identify the molecule to delete ----------
    res_type = pick_random_residue_type(thermo%is_active)

    call update_output_files(.false.)

    ! ---------- ENERGY CONSERVATION DURING TRANSLATION ----------
    ! Create 20 molecules
    thermo%chemical_potential = -0.1 ! Probability of creating the molecule is 1.0
    do i = 1, 50
        mol_index = primary%num%residues(res_type) + 1   ! if deletion removed it, create at new index
        call attempt_creation_move(res_type, mol_index)
        call update_output_files(.true.)
    end do

    write (*,*) "Inserted", primary%num%residues(res_type), "molecules"

    ! Compute energy
    call update_system_energy(primary)

    write (*,*)
    write (*,*) "energy before a few translation move", energy%total
    write (*,*)

    status%desired_block = 0


    ! Attemps translation move
    do i = 1, 5

        status%desired_block = i

        mol_index = pick_random_molecule_index(primary%num%residues(res_type))
        call attempt_translation_move(res_type, mol_index)
        call update_output_files(.true.)
        
        e_coul = energy%coulomb + energy%intra_coulomb
        e_long = energy%recip_coulomb + energy%ewald_self

        write (*,*) "energy%recip_coulomb after a few translation move", "recip_coulomb",   energy%recip_coulomb

    end do

    ! e_coul = energy%coulomb + energy%intra_coulomb
    ! e_long = energy%recip_coulomb + energy%ewald_self
    ! write (*,*) "energy after a few translation move", energy%total, e_coul, e_long

    call update_system_energy(primary)

    write (*,*) "recip_coulomb", energy%recip_coulomb

    ! write (*,*) "energy after recalcualtion", energy%total

    ! if (primary%num%residues(res_type) /= 1) then
    !     print *, "ERROR: topology must contain exactly 1 molecule"
    !     stop 1
    ! end if

    ! ! ---------- Delete the molecule ----------
    ! thermo%chemical_potential = -10 ! Probability of deleting the molecule is 1.0
    ! call attempt_deletion_move(res_type, mol_index)

    ! ! Recompute energy
    ! call compute_system_energy(primary)
    ! e_after_delete = energy%total

    ! if (abs(e_after_delete - 0.0) > tol) then
    !     print *, "FAIL: total energy mismatch"
    !     stop 1
    ! end if

!         mol_index = primary%num%residues(res_type) + 1   ! if deletion removed it, create at new index
   ! thermo%chemical_potential = -0.1 ! Probability of creating the molecule is 1.0
   ! thermo%chemical_potential = -10 ! Probability of deleting the molecule is 1.0

        !call attempt_creation_move(res_type, mol_index)
        ! call attempt_deletion_move(res_type, mol_index)

!         call attempt_rotation_move(res_type, mol_index)

    ! Recompute energy after creation
    ! call compute_system_energy(primary)
    ! e_after_create = energy%total
    ! e_coul = energy%coulomb + energy%intra_coulomb
    ! e_long = energy%recip_coulomb + energy%ewald_self

    ! if (abs(e_coul - ref_ecoul) > tol) then
    !     print *, "FAIL: Coulomb energy mismatch"
    !     stop 1
    ! end if

    ! if (abs(e_long - ref_elong) > tol) then
    !     print *, "FAIL: Long-range energy mismatch"
    !     stop 1
    ! end if

    ! if (abs(e_after_create - ref_total) > tol) then
    !     print *, "FAIL: total energy mismatch"
    !     stop 1
    ! end if

end program test_delete_and_create
