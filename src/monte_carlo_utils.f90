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
        rotation_axis = int(rand_uniform() * three) + 1 ! Random integer in range [1,3]

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
    ! Adjust the Monte Carlo translation and rotation step sizes using a
    ! Robbins–Monro adaptive scheme.
    !
    ! This routine updates the translational and rotational move amplitudes
    ! based on the acceptance ratios observed over recent trial moves.
    ! The adjustment follows the Robbins–Monro stochastic approximation:
    !
    !     step_new = step_old * exp( γ * (acc - TARGET) )
    !
    ! where:
    !     acc   = observed acceptance ratio for the move type
    !     γ     = learning rate (small positive constant)
    !     TARGET = desired acceptance probability
    !----------------------------------------------------------------------
    subroutine AdjustMoveStepSizes()

        implicit none

        ! Local variables
        real(real64) :: acc_trans, acc_rot
        real(real64), parameter :: gamma = 0.10d0  ! learning rate

        if (.not. input%recalibrate_moves) return

        ! Adjust translation step
        if (counter%trial_translations > MIN_TRIALS_FOR_RECALIBRATION) then

            acc_trans = real(counter%translations, real64) / &
                real(counter%trial_translations, real64)

            input%translation_step = input%translation_step * &
                                    exp(gamma * (acc_trans - TARGET_ACCEPTANCE))

            input%translation_step = max(MIN_TRANSLATION_STEP, &
                                    min(input%translation_step, MAX_TRANSLATION_STEP))

        end if

        ! Adjust rotational step
        if (counter%trial_rotations > MIN_TRIALS_FOR_RECALIBRATION) then

            acc_rot = real(counter%rotations,real64) / &
                    real(counter%trial_rotations,real64)

            input%rotation_step_angle = input%rotation_step_angle * &
                                        exp(gamma * (acc_rot - TARGET_ACCEPTANCE))

            input%rotation_step_angle = max(MIN_ROTATION_ANGLE, &
                                        min(input%rotation_step_angle, MAX_ROTATION_ANGLE))

        end if

    end subroutine AdjustMoveStepSizes

    !----------------------------------------------------------------------
    ! PickRandomResidueType: randomly selects an active residue type from the
    ! available types in [1, active_residue_count].
    !----------------------------------------------------------------------
    function PickRandomResidueType(is_active) result(residue_type)

        implicit none

        integer, dimension(:), intent(in) :: is_active
        integer :: residue_type
        integer :: i, n_active
        integer, allocatable :: active_indices(:)

        ! Count active residues
        n_active = count(is_active == 1)

        if (n_active == 0) then ! No active residue to pick from
            residue_type = 0
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
    ! PickRandomMoleculeIndex: randomly select a molecule index within a given
    ! residue type. Returns 0 if there are no molecules of that type.
    !----------------------------------------------------------------------
    function PickRandomMoleculeIndex(residue_count_for_type) result(molecule_index)

        integer, intent(in) :: residue_count_for_type
        integer :: molecule_index

        if (residue_count_for_type == 0) then

            molecule_index = 0
        
        else
        
            molecule_index = int(rand_uniform() * residue_count_for_type) + 1

            if (molecule_index > residue_count_for_type) then
            
                molecule_index = residue_count_for_type
            
            end if

        end if

    end function PickRandomMoleculeIndex

    !----------------------------------------------------------------------
    ! Compute the Metropolis acceptance probability for a Monte Carlo move.
    ! 
    ! This routine evaluates the acceptance probability for creation,
    ! deletion, translation, and rotation moves in a grand-canonical or
    ! canonical Monte Carlo simulation.
    !----------------------------------------------------------------------
    function compute_acceptance_probability(old, new, residue_type, move_type) result(probability)

        implicit none

        ! Input arguments
        type(energy_state), intent(in) :: old   ! Old energy states
        type(energy_state), intent(in) :: new   ! New energy states
        integer, intent(in) :: move_type        ! MC move type (TYPE_CREATION, TYPE_DELETION, TYPE_TRANSLATION, or TYPE_ROTATION)
        integer, intent(in) :: residue_type     ! Index of the residue type

        ! Local variables
        real(real64) :: N, Nplus1               ! Number of residues of this type, and Number + 1
        real(real64) :: mu                      ! Chemical potential
        real(real64) :: deltaU                  ! Energy difference ΔE between trial and current state
        real(real64) :: prefactor               ! Prefactor for probability calculation
        real(real64) :: lambda                  ! de Broglie wavelength in m

        ! Return value
        real(real64) :: probability             ! Acceptance probability (0 <= P <= 1)

        N = real(primary%num_residues(residue_type), real64)
        Nplus1 = N + 1.0_real64
        deltaU = new%total - old%total                  ! kcal/mol
        mu = input%chemical_potential(residue_type)     ! kcal/mol

        ! Compute factor based on move type
        select case (move_type)
            case (TYPE_CREATION)

                ! P_acc(N -> N+1) = min[1, (V / ((N+1) λ³)) * exp(-β * (ΔU - μ))]
                ! V in Å³, λ in Å, ΔU and μ in kcal/mol, β = 1/(kB T)
                ! Note: N+1 instead of N to avoid division by zero
                lambda = res%lambda(residue_type) ! Thermal de Broglie wavelength (A)
                prefactor = primary%volume / Nplus1 / lambda**3
                probability = min(1.0_real64, prefactor * exp(-beta * (deltaU - mu)))

                ! write (*,*) "creation"
                ! write (*,*) "beta", beta
                ! write (*,*) "lambda", lambda
                ! write (*,*) "deltaU", deltaU
                ! write (*,*) "mu", mu
                ! write (*,*) "prefactor", prefactor
                ! write (*,*) "probability", probability
                ! write (*,*)

            case (TYPE_DELETION)

                ! P_acc(N -> N-1) = min[1, (N λ³ / V) * exp(-β * (ΔU + μ))]
                ! λ in Å, V in Å³, ΔU and μ in kcal/mol, β = 1/(kB T)
                lambda = res%lambda(residue_type) ! Thermal de Broglie wavelength (A)
                prefactor = Nplus1 * lambda**3 / (primary%volume)
                probability = min(1.0_real64, prefactor * exp(-beta * (deltaU + mu)))

                ! write (*,*) N
                ! write (*,*) lambda, prefactor, probability
            
                ! write (*,*) "deletion"
                ! write (*,*) "beta", beta
                ! write (*,*) "lambda", lambda
                ! write (*,*) "deltaU", deltaU
                ! write (*,*) "mu", mu
                ! write (*,*) "prefactor", prefactor
                ! write (*,*) "probability", probability
                ! write (*,*)

            case (TYPE_TRANSLATION, TYPE_ROTATION)

                ! P_acc = min[1, exp(-β * ΔU)]
                probability = min(1.0_real64, exp(-beta * deltaU))

                ! write (*,*) "move"
                ! write (*,*) "beta", beta
                ! write (*,*) "deltaU", deltaU
                ! write (*,*) "probability", probability
                ! write (*,*)


            case default
                call AbortRun("Unknown move_type in compute_acceptance_probability!", 1)
        end select

    end function compute_acceptance_probability

    !----------------------------------------------------------------------
    ! Compute the Metropolis acceptance probability for a residue-type swap.
    !
    ! This routine evaluates the acceptance probability for replacing a
    ! molecule of type_old with one of type_new in a grand-canonical or
    ! semi-grand Monte Carlo simulation.
    !---------------------------------------------------------------------- 
    function swap_acceptance_probability(old, new, type_old, type_new) result(probability)

        implicit none

        ! Arguments
        type(energy_state), intent(in) :: old   ! Energy of system with old molecule
        type(energy_state), intent(in) :: new   ! Energy of system with new molecule
        integer, intent(in) :: type_old         ! Residue type being removed
        integer, intent(in) :: type_new         ! Residue type being inserted

        ! Local variables
        real(real64) :: deltaU                  ! Energy difference
        real(real64) :: mu_old, mu_new          ! Chemical potentials (kcal/mol)
        real(real64) :: N_old, N_new            ! Number of molecules per type
        real(real64) :: Nplus1

        ! Return value
        real(real64) :: probability             ! Acceptance probability (0 <= P <= 1)

        N_new = real(primary%num_residues(type_new), real64)
        N_old = real(primary%num_residues(type_old), real64)
        Nplus1 = N_new + 1.0_real64

        ! Chemical potentials
        mu_old = input%chemical_potential(type_old) ! kcal/mol
        mu_new = input%chemical_potential(type_new) ! kcal/mol
        deltaU = new%total - old%total         ! kcal/mol

        ! Swap acceptance probability:
        ! P_acc = min[1, (N_old / (N_new + 1)) * exp(-β (ΔE + μ_new - μ_old))]
        probability = min(1.0_real64, (N_old / Nplus1) * exp(-beta * (deltaU + mu_new - mu_old)))

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

            ! #tocheck
            !call ComputeRecipEnergySingleMol(residue_type, molecule_index, old%recip_coulomb)
            old%recip_coulomb = energy%recip_coulomb

            ! Recalculate total energy
            old%total = old%non_coulomb + old%coulomb + old%recip_coulomb + old%ewald_self + old%intra_coulomb

        else if (deletion_flag) then

            ! Note: In deletion scenario, compute all energy components
            call ComputeEwaldSelfInteractionSingleMol(residue_type, old%ewald_self)
            call ComputeIntraResidueRealCoulombEnergySingleMol(residue_type, molecule_index, old%intra_coulomb)
            call ComputePairInteractionEnergy_singlemol(primary, residue_type, molecule_index, old%non_coulomb, old%coulomb)
            call ComputeRecipEnergySingleMol(residue_type, molecule_index, old%recip_coulomb)

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
    ! Save the current state of a molecule for a Monte Carlo move. For
    ! translation: saves center-of-mass (COM). For rotation: saves site offsets.
    ! Also saves Fourier terms for later restoration if the move is rejected
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
    ! Restore a molecule's previous state (COM or site offsets) and Fourier terms
    ! if a Monte Carlo move is rejected.
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
    ! Generates a random position for the new molecule and copies/orients
    ! its atomic geometry. Can take geometry from a reservoir or apply
    ! a random rotation if no reservoir exists.
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
    ! Restores molecule and atom counts, and resets Fourier states if a
    ! creation move is rejected.
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
    ! Compute the excess chemical potential (μ_ex) for each residue type
    ! using Widom particle insertion method, and optionally the ideal chemical
    ! potential (μ_ideal) for reporting purposes.
    !------------------------------------------------------------------------------
    subroutine CalculateExcessMu()

        implicit none

        ! Local variables
        integer :: type_residue     ! Residue type index
        integer :: N                ! Number of molecules of current residue
        real(real64) :: avg_weight  ! Average Boltzmann weight for Widom sampling
        real(real64) :: mu_ideal    ! Ideal chemical potential (kcal/mol)
        real(real64) :: lambda      ! Thermal de Broglie wavelength (m)
        real(real64) :: temperature ! Temperature (K)
        real(real64) :: volume      ! Simulation box volume (m^3)
        real(real64) :: rho         ! Number density (molecules/m^3)

        ! Loop over all residue types
        do type_residue = 1, nb%type_residue
            if (widom_stat%sample(type_residue) > 0) then

                ! Compute average Boltzmann factor from Widom sampling
                ! <exp(-β ΔU)> = sum_weights / N_samples
                avg_weight = widom_stat%weight(type_residue) / real(widom_stat%sample(type_residue), kind=real64)

                ! Compute excess chemical potential (kcal/mol)
                ! μ_ex = - k_B * T * ln(<exp(-β ΔU)>)
                temperature = input%temperature
                widom_stat%mu_ex(type_residue) = - KB_kcalmol * temperature * log(avg_weight) ! kcal/mol

                ! Compute ideal gas chemical potential (kcal/mol)
                ! μ_ideal = k_B * T * ln(ρ * Λ^3)
                lambda = res%lambda(type_residue)                   ! Thermal de Broglie wavelength (A)
                N = primary%num_residues(type_residue)              ! Number of molecules
                volume = primary%volume                             ! Box volume (A^3)
                rho = real(N, kind=real64) / volume                 ! Number density (molecules/A^3)
                mu_ideal = KB_kcalmol * temperature * log(rho * lambda**3) ! Ideal chemical potential (kcal/mol)

                widom_stat%mu_tot(type_residue) = mu_ideal + widom_stat%mu_ex(type_residue)

            end if
        end do

    end subroutine CalculateExcessMu

end module monte_carlo_utils
