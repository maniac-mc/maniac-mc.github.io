module monte_carlo_utils

    use simulation_state
    use random_utils
    use constants
    use ewald_kvectors
    use ewald_phase
    use ewald_energy
    use output_utils
    use energy_utils
    use helper_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !========================================================
    ! Subroutine: ApplyRandomRotation
    !
    ! Rotates all atoms of a residue around a randomly chosen
    ! axis (X, Y, or Z) by a random angle. Angle can be a
    ! small perturbation or a full 0-2π rotation.
    !
    ! Inputs:
    !   res_type - integer, index of the residue type
    !   mol_index - integer, index of the molecule
    !   full_rotation - optional logical, if .true. use full rotation
    !========================================================
    subroutine ApplyRandomRotation(res_type, mol_index, full_rotation)

        implicit none

        ! Input arguments
        integer, intent(in) :: res_type           ! Index of the residue type
        integer, intent(in) :: mol_index          ! Index of the molecule
        logical, intent(in), optional :: full_rotation ! Flag for full or small-step rotation

        ! Local variables
        integer :: rotation_axis                  ! Chosen rotation axis (1=X, 2=Y, 3=Z)
        integer :: n_atoms                        ! Number of atoms in residue
        logical :: use_full_rotation              ! Actual value of full_rotation (defaults to .false.)
        real(real64) :: rotation_matrix(3,3)      ! 3x3 rotation matrix for the rotation
        real(real64) :: theta                     ! Rotation angle in radians

        ! Handle optional argument: default = .false.
        use_full_rotation = .false.
        if (present(full_rotation)) use_full_rotation = full_rotation

        ! Exit if single-atom residue (nothing to rotate)
        n_atoms = nb%atom_in_residue(res_type)
        if (n_atoms == 1) return

        ! Choose rotation angle
        theta = ChooseRotationAngle(use_full_rotation)

        ! Choose random axis (1=X, 2=Y, 3=Z)
        rotation_axis = int(rand_uniform() * three) + 1 ! Random integer in [1,3]

        ! Set rotation matrix based on axis
        rotation_matrix = RotationMatrix(rotation_axis, theta)

        ! Apply rotation to all atoms in the residue
        primary%site_offset(:, res_type, mol_index, 1:n_atoms) = &
            matmul(rotation_matrix, primary%site_offset(:, res_type, mol_index, 1:n_atoms))

    end subroutine ApplyRandomRotation

    ! Returns a random rotation angle (radians); small-step if use_full_rotation=.false., full [0,2π] if .true.
    function ChooseRotationAngle(use_full_rotation) result(theta)

        implicit none

        ! Input variables
        logical, intent(in) :: use_full_rotation

        ! Output variables
        real(real64) :: theta

        if (.not. use_full_rotation) then

            ! For small-step mode, make sure rotation_step_angle is reasonable
            if (input%rotation_step_angle <= zero .or. input%rotation_step_angle > TWOPI) then
                call AbortRun('Invalid rotation_step_angle in ChooseRotationAngle')
            end if

            ! Use small rotation
            theta = (rand_uniform() - half) * input%rotation_step_angle
        else
            ! Use large rotation
            theta = rand_uniform() * TWOPI
        end if

    end function ChooseRotationAngle

    !----------------------------------------------------------------------
    ! Adjusts the Monte Carlo translation and rotation step sizes dynamically
    ! This subroutine recalibrates the translational and rotational move steps
    ! based on the observed acceptance ratios of recent MC trials. If the
    ! acceptance is too high, the step size is increased (up to a defined maximum).
    ! If the acceptance is too low, the step size is decreased (down to a defined minimum).
    !----------------------------------------------------------------------
    subroutine AdjustMoveStepSizes()

        implicit none

        real(real64) :: acc_trans, acc_rot

        if (.not. input%recalibrate_moves) return

        ! Adjust translation step
        if (counter%trial_translations > MIN_TRIALS_FOR_RECALIBRATION) then

            acc_trans = real(counter%translations, real64) / &
                real(counter%trial_translations, real64)

            if (acc_trans - TARGET_ACCEPTANCE > TOL_ACCEPTANCE) then

                ! acceptance too high → increase step
                input%translation_step = min(input%translation_step * 1.05d0, MAX_TRANSLATION_STEP)

            else if (acc_trans - TARGET_ACCEPTANCE < TOL_ACCEPTANCE) then

                ! acceptance too low → decrease step
                input%translation_step = max(input%translation_step * 0.95d0, MIN_TRANSLATION_STEP)

            end if

        end if

        ! Adjust rotational step
        if (counter%trial_rotations > MIN_TRIALS_FOR_RECALIBRATION) then

            acc_rot = real(counter%rotations,real64) / &
                    real(counter%trial_rotations,real64)

            if (acc_rot - TARGET_ACCEPTANCE > TOL_ACCEPTANCE) then

                ! acceptance too high → increase step
                input%rotation_step_angle = min(input%rotation_step_angle * 1.05d0, MAX_ROTATION_ANGLE)

            else if (acc_rot - TARGET_ACCEPTANCE < TOL_ACCEPTANCE) then

                ! acceptance too low → decrease step
                input%rotation_step_angle = max(input%rotation_step_angle * 0.95d0, MIN_ROTATION_ANGLE)
            
            end if

        end if

    end subroutine AdjustMoveStepSizes

    !----------------------------------------------------------------------
    ! PickRandomResidueType: randomly selects an active residue type from the
    ! available types in [1, active_residue_count].
    !----------------------------------------------------------------------
    function PickRandomResidueType(is_active) result(residue_type)
        integer, dimension(:), intent(in) :: is_active
        integer :: residue_type
        integer :: i, n_active
        integer, allocatable :: active_indices(:)

        ! Count active residues
        n_active = count(is_active == 1)

        if (n_active == 0) then
            residue_type = 0  ! or some error code
            return
        end if

        ! Collect indices of active residues
        allocate(active_indices(n_active))
        n_active = 0
        do i = 1, size(is_active)
            if (is_active(i) == 1) then
                n_active = n_active + 1
                active_indices(n_active) = i
            end if
        end do

        ! Pick a random index among active ones
        residue_type = active_indices(INT(rand_uniform() * n_active) + 1)

        deallocate(active_indices)
    end function PickRandomResidueType

    !----------------------------------------------------------------------
    ! PickRandomMoleculeIndex: randomly selects a molecule index within a given
    ! residue type. Returns 0 if there are no molecules of that type.
    !----------------------------------------------------------------------
    function PickRandomMoleculeIndex(residue_count_for_type) result(molecule_index)

        integer, intent(in) :: residue_count_for_type
        integer :: molecule_index

        if (residue_count_for_type == 0) then
            molecule_index = 0
        else
            molecule_index = INT(rand_uniform() * residue_count_for_type) + 1
            if (molecule_index > residue_count_for_type) molecule_index = residue_count_for_type
        end if

    end function PickRandomMoleculeIndex

    !----------------------------------------------------------------------
    ! Compute the Metropolis acceptance probability for a Monte Carlo move.
    ! 
    ! This routine evaluates the acceptance probability for creation,
    ! deletion, translation, and rotation moves in a grand-canonical or
    ! canonical Monte Carlo simulation.  The expression uses:
    !
    !   - ΔE  = new%total – old%total          (energy difference)
    !   - β   = 1 / (kB * T)                   (inverse temperature)
    !   - φ   = fugacity of the residue type
    !   - V   = simulation box volume
    !   - N   = current number of residues of that type
    !
    ! Move-dependent acceptance rules:
    !   * Creation:   min(1, (φ V / (N+1)) * exp(-β ΔE))
    !   * Deletion:   min(1, ((N+1) / (φ V)) * exp(-β ΔE))
    !   * Translation/Rotation:  min(1, exp(-β ΔE))
    !----------------------------------------------------------------------
    function mc_acceptance_probability(old, new, residue_type, move_type) result(probability)

        implicit none

        ! Input arguments
        type(energy_state), intent(in) :: old   ! Old energy states
        type(energy_state), intent(in) :: new   ! New energy states
        integer, intent(in) :: move_type        ! MC move type (TYPE_CREATION, TYPE_DELETION, TYPE_TRANSLATION, or TYPE_ROTATION)
        integer, intent(in) :: residue_type     ! Index of the residue type
        
        ! Local variables
        real(real64) :: N                       ! Number of residues of this type (local copy)
        real(real64) :: V                       ! Simulation box volume (local copy)
        real(real64) :: phi                     ! Fugacity of the residue type (local copy)
        real(real64) :: T                       ! Temperature in units of kB*T (optional local copy for clarity)
        real(real64) :: delta_e                 ! Energy difference ΔE between trial and current state
        real(real64) :: beta                    ! 1/kB T

        ! Return value
        real(real64) :: probability             ! Acceptance probability (0 <= P <= 1)

        N   = real(primary%num_residues(residue_type), real64)
        V   = primary%volume                    ! Angstrom^3
        phi = input%fugacity(residue_type)      ! Fugacity in Angstrom^-3
        T = input%temp_K                        ! Kelvin
        beta = 1/(KB_kcalmol*T)                 ! 1/(kcal/mol)
        delta_e = new%total - old%total         ! kcal/mol

        ! Compute factor based on move type
        select case (move_type)
            case (TYPE_CREATION)

                probability = min(one, (phi * V / (N + one)) * exp(-beta * delta_e)) ! Note: N+1 instead of N to avoid division by zero

            case (TYPE_DELETION)

                probability = min(one, ((N + one) / (phi * V)) * exp(-beta * delta_e))

            case (TYPE_TRANSLATION, TYPE_ROTATION)

                probability = min(one, exp(-beta * delta_e))

            case default

                call AbortRun("Unknown move_type in mc_acceptance_probability!", 1)

        end select

    end function mc_acceptance_probability

    !----------------------------------------------------------------------
    ! Compute the Metropolis acceptance probability for a residue-type swap.
    !
    ! This routine evaluates the acceptance probability for replacing a
    ! molecule of type_old with one of type_new in a grand-canonical or
    ! semi-grand Monte Carlo simulation.  The expression uses:
    !
    !   - ΔE  = new%total - old%total          (energy difference)
    !   - β   = 1 / (kB * T)                   (inverse temperature)
    !   - φ   = fugacity of each residue type
    !   - N   = current number of residues of each type
    !
    ! Swap move acceptance rule:
    !
    !   P = min(1,
    !           (φ_new / φ_old)
    !         * (N_old / (N_new + 1))
    !         * exp(-β ΔE) )
    !
    !---------------------------------------------------------------------- 
    function swap_acceptance_probability(old, new, type_old, type_new) result(probability)

        implicit none

        ! Arguments
        type(energy_state), intent(in) :: old   ! Energy of system with old molecule
        type(energy_state), intent(in) :: new   ! Energy of system with new molecule
        integer, intent(in) :: type_old         ! Residue type being removed
        integer, intent(in) :: type_new         ! Residue type being inserted

        ! Local valriables
        real(real64) :: delta_e                 ! Energy difference
        real(real64) :: phi_old, phi_new        ! Fugacities of species
        real(real64) :: T                       ! Temperature in units of kB*T
        real(real64) :: beta                    ! 1/kB T
        real(real64) :: N_old, N_new            ! Number of molecule per types

        ! Return value
        real(real64) :: probability             ! Acceptance probability (0 <= P <= 1)

        N_new = real(primary%num_residues(type_new), real64)
        N_old = real(primary%num_residues(type_old), real64)
        phi_old = input%fugacity(type_old)      ! Fugacity in Angstrom^-3
        phi_new = input%fugacity(type_new)      ! Fugacity in Angstrom^-3
        T = input%temp_K                        ! Kelvin
        beta = 1/(KB_kcalmol*T)                 ! 1/(kcal/mol)
        delta_e = new%total - old%total         ! kcal/mol

        ! Swap acceptance probability
        probability = min(one, (phi_new / phi_old) * (N_old / (N_new + one)) * exp(-beta * delta_e))

    end function swap_acceptance_probability

    !---------------------------------------------------------------------------
    ! Purpose:
    !   Compute the updated energy of a single molecule after a trial move
    !   for use in the Monte Carlo acceptance test.
    !---------------------------------------------------------------------------
    subroutine ComputeNewEnergy(residue_type, molecule_index, new, is_creation, is_deletion)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type         ! Residue type to be moved
        integer, intent(in) :: molecule_index       ! Index of the molecule to move
        type(energy_state), intent(out) :: new      ! New energy states
        logical, intent(in), optional :: is_creation
        logical, intent(in), optional :: is_deletion

        ! Local variables
        logical :: creation_flag
        logical :: deletion_flag

        creation_flag = present_or_false(is_creation)
        deletion_flag = present_or_false(is_deletion)

        if (creation_flag) then
        
            ! Note: In creation scenario, compute all energy components

            call SingleMolFourierTerms(residue_type, molecule_index)
            call ComputeRecipEnergySingleMol(residue_type, molecule_index, new%recip_coulomb, is_creation = creation_flag)
            call ComputePairInteractionEnergy_singlemol(primary, residue_type, molecule_index, new%non_coulomb, new%coulomb)
            call ComputeEwaldSelfInteractionSingleMol(residue_type, new%ewald_self)
            call ComputeIntraResidueRealCoulombEnergySingleMol(residue_type, molecule_index, new%intra_coulomb)
        
            ! Recalculate total energy
            new%total = new%non_coulomb + new%coulomb + new%recip_coulomb + new%ewald_self + new%intra_coulomb
        
        else if (deletion_flag) then

            ! Note: Most energy terms in the absence of a molecule are 0
            ! --> One must only recalculate energy%recip_coulomb

            new%non_coulomb = 0
            new%coulomb = 0
            new%ewald_self = 0
            new%intra_coulomb = 0
            call ComputeRecipEnergySingleMol(residue_type, molecule_index, new%recip_coulomb, is_creation = deletion_flag)

            ! Recalculate total energy
            new%total = new%non_coulomb + new%coulomb + new%recip_coulomb + new%ewald_self + new%intra_coulomb
        
        else

            ! Note, for simple move (translation or rotation), one only needs to
            ! recompute reciprocal and pairwise interactions

            call SingleMolFourierTerms(residue_type, molecule_index)
            call ComputeRecipEnergySingleMol(residue_type, molecule_index, new%recip_coulomb)
            call ComputePairInteractionEnergy_singlemol(primary, residue_type, molecule_index, new%non_coulomb, new%coulomb)

            ! Recalculate total energy
            new%total = new%non_coulomb + new%coulomb + new%recip_coulomb

        end if

    end subroutine ComputeNewEnergy

    !---------------------------------------------------------------------------
    ! Purpose:
    !   Compute the previous energy of a single molecule before a trial move
    !   for use in the Monte Carlo acceptance test.
    !---------------------------------------------------------------------------
    subroutine ComputeOldEnergy(residue_type, molecule_index, old, is_creation, is_deletion)

        implicit none

        ! Input variables
        integer, intent(in) :: residue_type         ! Residue type to be moved
        integer, intent(in) :: molecule_index       ! Index of the molecule to move
        type(energy_state), intent(out) :: old      ! Old energy states
        logical, intent(in), optional :: is_creation
        logical, intent(in), optional :: is_deletion

        ! Local variables
        logical :: creation_flag
        logical :: deletion_flag

        creation_flag = present_or_false(is_creation)
        deletion_flag = present_or_false(is_deletion)

        if (creation_flag) then

            ! Note: Most energy terms in the absence of a molecule are 0
            ! --> One must only recalculate energy%recip_coulomb

            old%non_coulomb = zero
            old%coulomb = zero
            old%ewald_self = zero
            old%intra_coulomb = zero
            old%recip_coulomb = energy%recip_coulomb ! global system reciprocal energy

            ! Recalculate total energy
            old%total = old%non_coulomb + old%coulomb + old%recip_coulomb + old%ewald_self + old%intra_coulomb

        else if (deletion_flag) then

            ! Note: In deletion scenario, compute all energy components

            call ComputeEwaldSelfInteractionSingleMol(residue_type, old%ewald_self)
            call ComputeIntraResidueRealCoulombEnergySingleMol(residue_type, molecule_index, old%intra_coulomb)
            call ComputePairInteractionEnergy_singlemol(primary, residue_type, molecule_index, old%non_coulomb, old%coulomb)
            old%recip_coulomb = energy%recip_coulomb

            ! Recalculate total energy
            old%total = old%non_coulomb + old%coulomb + old%recip_coulomb + old%ewald_self + old%intra_coulomb

        else

            ! Note, for simple move (translation or rotation), one only needs to
            ! recompute reciprocal and pairwise interactions

            call ComputeRecipEnergySingleMol(residue_type, molecule_index, old%recip_coulomb)
            call ComputePairInteractionEnergy_singlemol(primary, residue_type, molecule_index, old%non_coulomb, old%coulomb)

            ! Recalculate total energy
            old%total = old%non_coulomb + old%coulomb + old%recip_coulomb

        end if

    end subroutine ComputeOldEnergy

    !---------------------------------------------------------------------------
    ! Purpose:
    !   Update the global system energy and Monte Carlo counters after
    !   accepting a trial move (translation or rotation) of a molecule.
    !
    ! Inputs:
    !   old         - Type(energy_state): Energy of the molecule before the move
    !   new         - Type(energy_state): Energy of the molecule after the move
    !
    ! Input/Output:
    !   counter_var - Integer: Monte Carlo counter for successful moves
    !                 (e.g., translations or rotations), incremented if move accepted
    !---------------------------------------------------------------------------
    subroutine AcceptMove(old, new, counter_var)

        ! Input arguments
        type(energy_state), intent(in) :: old, new  ! Old and new energy states
        integer, intent(inout) :: counter_var       ! Counter for succesfull move

        energy%recip_coulomb = energy%recip_coulomb + new%recip_coulomb - old%recip_coulomb
        energy%non_coulomb = energy%non_coulomb + new%non_coulomb - old%non_coulomb
        energy%coulomb = energy%coulomb + new%coulomb - old%coulomb
        energy%total = energy%total + new%total - old%total
        counter_var = counter_var + 1

    end subroutine

    !---------------------------------------------------------------------------
    ! Subroutine: SaveMoleculeState
    !
    ! Purpose:
    !   Save the current state of a molecule for a Monte Carlo move.
    !   - For translation: saves center-of-mass (COM)
    !   - For rotation: saves site offsets
    !   - Also saves Fourier terms for later restoration if the move is rejected
    !
    ! Inputs:
    !   residue_type   - Integer: Residue type of the molecule
    !   molecule_index - Integer: Index of the molecule
    !
    ! Outputs:
    !   com_old         - Real(3) array: center-of-mass (for translation)
    !   site_offset_old - Real(3, natoms) array: site offsets (for rotation)
    ! Notes:
    !   - Only one of com_old or site_offset_old is allocated by the caller
    !   - Fourier terms are always saved
    !---------------------------------------------------------------------------
    subroutine SaveMoleculeState(res_type, mol_index, com_old, offset_old)

        implicit none

        ! Input arguments
        integer, intent(in) :: res_type         ! Residue type of the molecule
        integer, intent(in) :: mol_index        ! Index of the molecule
        real(real64), intent(out), optional :: com_old(3) ! Center-of-mass before move (for translation)
        real(real64), intent(out), optional :: offset_old(:, :) ! Site offsets before move (for rotation)

        ! Local variable
        integer :: natoms

        ! Save Fourier terms
        call SaveSingleMolFourierTerms(res_type, mol_index)

        ! Save center-of-mass if requested (translation)
        if (present(com_old)) then
            com_old(:) = primary%mol_com(:, res_type, mol_index)
        end if

        ! Save site offsets if requested (rotation)
        if (present(offset_old)) then
            natoms = nb%atom_in_residue(res_type)
            offset_old(:, 1:natoms) = primary%site_offset(:, res_type, mol_index, 1:natoms)
        end if

    end subroutine SaveMoleculeState

    !---------------------------------------------------------------------------
    ! Subroutine: RejectMoleculeMove
    !
    ! Purpose:
    !   Restore a molecule's previous state (COM or site offsets) and Fourier terms
    !   if a Monte Carlo move is rejected.
    !
    ! Inputs:
    !   res_type - Residue type of the molecule
    !   mol_index - Index of the molecule
    !   com_old - Previous COM (for translation, optional)
    !   site_offset_old - Previous site offsets (for rotation, optional)
    !---------------------------------------------------------------------------
    subroutine RejectMoleculeMove(res_type, mol_index, com_old, site_offset_old)
        implicit none

        integer, intent(in) :: res_type       ! Residue type of the molecule
        integer, intent(in) :: mol_index      ! Index of the molecule
        real(real64), intent(in), optional :: com_old(3) ! Previous COM (translation)
        real(real64), intent(in), optional :: site_offset_old(:, :) ! Previous site offsets (rotation)

        ! Local variable
        integer :: natoms

        ! Restore COM if present (translation)
        if (present(com_old)) then
            primary%mol_com(:, res_type, mol_index) = com_old(:)
        end if

        ! Restore site offsets if present (rotation)
        if (present(site_offset_old)) then
            natoms = nb%atom_in_residue(res_type)
            primary%site_offset(:, res_type, mol_index, 1:natoms) = &
                site_offset_old(:, 1:natoms)
        end if

        ! Restore Fourier states
        call RestoreSingleMolFourier(res_type, mol_index)

    end subroutine RejectMoleculeMove

    !---------------------------------------------------------------------------
    ! Subroutine: InsertAndOrientMolecule
    !
    ! Purpose:
    !   Generates a random position for the new molecule and copies/orients
    !   its atomic geometry. Can take geometry from a reservoir or apply
    !   a random rotation if no reservoir exists.
    !
    !---------------------------------------------------------------------------
    subroutine InsertAndOrientMolecule(residue_type, molecule_index, rand_mol_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type     ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID
        integer, intent(out) :: rand_mol_index  ! Randomly selected molecule index from the reservoir

        ! Local variables
        logical :: full_rotation                ! Flag indicating whether a full 360° random rotation should be applied
        real(real64) :: random_nmb              ! Uniform random number in [0,1), used for random index selection
        real(real64) :: trial_pos(3)            ! Random numbers for initial molecule position in the box

        ! Generate a random position in the simulation box
        call random_number(trial_pos) ! Random numbers in [0,1)
        primary%mol_com(:, residue_type, molecule_index) = primary%bounds(:,1) &
            + matmul(primary%matrix, trial_pos)

        ! Copy geometry from reservoir or rotate if no reservoir
        if (has_reservoir) then

            ! Pick a random (and existing) molecule in the reservoir
            call random_number(random_nmb) ! generates random_nmb in [0,1)
            rand_mol_index = int(random_nmb * reservoir%num_residues(residue_type)) + 1 ! random integer in [1, N]

            ! Copy site offsets from the chosen molecule
            primary%site_offset(:, residue_type, molecule_index, 1:nb%atom_in_residue(residue_type)) = &
                reservoir%site_offset(:, residue_type, rand_mol_index, 1:nb%atom_in_residue(residue_type))

        else

            ! Copy site offsets from the first molecule
            primary%site_offset(:, residue_type, molecule_index, 1:nb%atom_in_residue(residue_type)) = &
                primary%site_offset(:, residue_type, 1, 1:nb%atom_in_residue(residue_type))

            ! Rotate the new molecule randomly (using full 360° rotation)
            full_rotation = .true.
            call ApplyRandomRotation(residue_type, molecule_index, full_rotation)

        end if

    end subroutine InsertAndOrientMolecule

    !---------------------------------------------------------------------------
    ! Subroutine: RejectCreationMove
    !
    ! Purpose:
    !   Restores molecule and atom counts, and resets Fourier states if a
    !   creation move is rejected.
    !
    !---------------------------------------------------------------------------
    subroutine RejectCreationMove(residue_type, molecule_index)

        implicit none

        integer, intent(in) :: residue_type     ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID

        ! Restore previous residue/atom numbers
        primary%num_atoms = primary%num_atoms - nb%atom_in_residue(residue_type)
        primary%num_residues(residue_type) = primary%num_residues(residue_type) - 1

        ! Restore Fourier states (ik_alloc and dk_alloc, all zeros)
        call RestoreSingleMolFourier(residue_type, molecule_index)

    end subroutine RejectCreationMove

    !------------------------------------------------------------------------------
    ! Subroutine: CalculateExcessMu
    !
    ! Purpose:
    !   Compute the excess chemical potential (μ_ex) for each residue type
    !   using Widom particle insertion method, and optionally the ideal chemical
    !   potential (μ_ideal) for reporting purposes.
    !
    ! Formulas:
    !   1. Excess chemical potential (Widom):
    !        μ_ex = - k_B * T * ln(<exp(-β ΔU)>)
    !          - <exp(-β ΔU)> = average Boltzmann weight from Widom sampling
    !          - k_B: Boltzmann constant in kcal/mol/K
    !          - T: temperature in Kelvin
    !
    !   2. Ideal chemical potential:
    !        μ_ideal = k_B * T * ln(ρ * Λ^3)
    !          - ρ = N / V = number density (molecules/m^3)
    !          - Λ = hbar / sqrt(2 * π * m * k_B * T) = thermal de Broglie wavelength (m)
    !          - m: mass per molecule in kg
    !------------------------------------------------------------------------------
    subroutine CalculateExcessMu()

        implicit none

        ! Local variables
        integer :: type_residue     ! Residue type index
        integer :: N                ! Number of molecules of current residue
        real(real64) :: avg_weight  ! Average Boltzmann weight for Widom sampling
        real(real64) :: mu_ideal    ! Ideal chemical potential (kcal/mol)
        real(real64) :: m           ! Mass per molecule (kg)
        real(real64) :: Lambda      ! Thermal de Broglie wavelength (m)
        real(real64) :: T           ! Temperature (K)
        real(real64) :: V           ! Simulation box volume (m^3)
        real(real64) :: rho         ! Number density (molecules/m^3)

        ! Loop over all residue types
        do type_residue = 1, nb%type_residue
            if (widom_stat%sample(type_residue) > 0) then

                ! ----------------------------------------------------------------
                ! 1. Compute average Boltzmann factor from Widom sampling
                !    <exp(-β ΔU)> = sum_weights / N_samples
                ! ----------------------------------------------------------------
                avg_weight = widom_stat%weight(type_residue) / real(widom_stat%sample(type_residue), kind=real64)

                ! ----------------------------------------------------------------
                ! 2. Compute excess chemical potential (kcal/mol)
                !    μ_ex = - k_B * T * ln(<exp(-β ΔU)>)
                ! ----------------------------------------------------------------
                widom_stat%mu_ex(type_residue) = - KB_kcalmol * input%temp_K * log(avg_weight) ! kcal/mol

                ! ----------------------------------------------------------------
                ! 3. Compute ideal gas chemical potential (kcal/mol)
                !    μ_ideal = k_B * T * ln(ρ * Λ^3)
                ! ----------------------------------------------------------------
                T = input%temp_K                                  ! Temperature (K)
                m = res%mass_residue(type_residue) * G_TO_KG / NA ! Mass per molecule (kg)
                Lambda = Hplank / sqrt(TWOPI * m * KB_JK * T)     ! Thermal de Broglie wavelength (m)

                N = primary%num_residues(type_residue)           ! Number of molecules
                V = primary%volume * A3_TO_M3                    ! Box volume (m^3)
                rho = real(N, kind=real64) / V                   ! Number density (molecules/m^3)

                mu_ideal = KB_kcalmol * T * log(rho * Lambda**3) ! Ideal chemical potential (kcal/mol)

                widom_stat%mu_tot(type_residue) = mu_ideal + widom_stat%mu_ex(type_residue)

            end if
        end do

    end subroutine CalculateExcessMu


end module monte_carlo_utils
