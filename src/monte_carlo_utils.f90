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

    !---------------------------------------------------------------------------
    ! Rotates all atoms of a residue around a randomly chosen
    ! axis (X, Y, or Z) by a random angle. Angle can be a
    ! small perturbation or a full 0-2π rotation.
    !---------------------------------------------------------------------------
    subroutine apply_random_rotation(res_type, mol_index, full_rotation)

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
        type(type_coordinate), pointer :: coord   ! Pointer for host or guest coordinate

        ! Handle optional argument
        use_full_rotation = .false.
        if (present(full_rotation)) use_full_rotation = full_rotation

        ! Exit if single-atom residue (nothing to rotate)
        n_atoms = res%atom(res_type)
        if (n_atoms == 1) return

        ! Choose rotation angle
        theta = choose_rotation_angle(use_full_rotation)

        ! Choose random axis (1=X, 2=Y, 3=Z)
        rotation_axis = int(rand_uniform() * three) + 1 ! Random integer in range [1,3]

        ! Set rotation matrix based on axis
        rotation_matrix = return_rotation_matrix(rotation_axis, theta)

        ! Return the correct pointer (host or guest)
        coord => get_coord(res_type)

        ! Apply rotation to all atoms in the residue
        coord%offset(:, res_type, mol_index, 1:n_atoms) = &
            matmul(rotation_matrix, coord%offset(:, res_type, mol_index, 1:n_atoms))

    end subroutine apply_random_rotation

    !---------------------------------------------------------------------------
    ! Returns a random rotation angle (radians), 
    ! small-step if use_full_rotation=.false., full [0,2π] if use_full_rotation =.true.
    !---------------------------------------------------------------------------
    function choose_rotation_angle(use_full_rotation) result(theta)

        ! Input variables
        logical, intent(in) :: use_full_rotation

        ! Output variables
        real(real64) :: theta

        if (.not. use_full_rotation) then

            ! For small-step mode, make sure rotation_step_angle is reasonable
            if (mc_input%rotation_step_angle <= zero .or. mc_input%rotation_step_angle > TWOPI) then
                call abort_run('Invalid rotation_step_angle in choose_rotation_angle')
            end if

            ! Use small rotation
            theta = (rand_uniform() - half) * mc_input%rotation_step_angle
        else
            ! Use large rotation
            theta = rand_uniform() * TWOPI
        end if

    end function choose_rotation_angle

    !----------------------------------------------------------------------
    ! Adjust the Monte Carlo translation and rotation step sizes using a
    ! Robbins–Monro adaptive scheme.
    !----------------------------------------------------------------------
    subroutine adjust_move_step_sizes()

        ! Local variables
        real(real64) :: acc_trans, acc_rot
        real(real64), parameter :: gamma = 0.10d0  ! learning rate

        if (.not. mc_input%recalibrate_moves) return

        ! Adjust translation step
        if (counter%translations(1) > MIN_TRIALS_FOR_RECALIBRATION) then

            acc_trans = real(counter%translations(2), real64) / &
                real(counter%translations(1), real64)

            mc_input%translation_step = mc_input%translation_step * &
                                    exp(gamma * (acc_trans - TARGET_ACCEPTANCE))

            mc_input%translation_step = max(MIN_TRANSLATION_STEP, &
                                    min(mc_input%translation_step, MAX_TRANSLATION_STEP))

        end if

        ! Adjust rotational step
        if (counter%rotations(1) > MIN_TRIALS_FOR_RECALIBRATION) then

            acc_rot = real(counter%rotations(2),real64) / &
                    real(counter%rotations(1),real64)

            mc_input%rotation_step_angle = mc_input%rotation_step_angle * &
                                        exp(gamma * (acc_rot - TARGET_ACCEPTANCE))

            mc_input%rotation_step_angle = max(MIN_ROTATION_ANGLE, &
                                        min(mc_input%rotation_step_angle, MAX_ROTATION_ANGLE))

        end if

    end subroutine adjust_move_step_sizes

    !----------------------------------------------------------------------
    ! Randomly selects an active residue type from the
    ! available types in [1, active_residue_count].
    !----------------------------------------------------------------------
    function pick_random_residue_type(is_active) result(res_type)

        ! Input parameter
        logical, dimension(:), intent(in) :: is_active
        
        ! Output parameter
        integer :: res_type

        ! Local variables
        integer :: i, n_active
        integer, allocatable :: active_indices(:)

        ! Count active residues
        n_active = count(is_active)

        if (n_active == 0) then ! No active residue to pick from
            res_type = 0
            return
        end if

        ! Collect indices of active residues
        allocate(active_indices(n_active))
        n_active = 0
        do i = 1, size(is_active)
            if (is_active(i)) then
                n_active = n_active + 1
                active_indices(n_active) = i
            end if
        end do

        ! Pick a random index among active ones
        res_type = active_indices(INT(rand_uniform() * n_active) + 1)

        deallocate(active_indices)

    end function pick_random_residue_type

    !----------------------------------------------------------------------
    ! Randomly select a molecule index within a given
    ! residue type. Returns 0 if there are no molecules of that type.
    !----------------------------------------------------------------------
    function pick_random_molecule_index(residue_count_for_type) result(mol_index)

        ! Input parameter
        integer, intent(in) :: residue_count_for_type
        
        ! Output parameter
        integer :: mol_index

        if (residue_count_for_type == 0) then
            mol_index = 0
        else
        
            mol_index = int(rand_uniform() * residue_count_for_type) + 1

            if (mol_index > residue_count_for_type) mol_index = residue_count_for_type

        end if

    end function pick_random_molecule_index

    !----------------------------------------------------------------------
    ! Compute the Metropolis acceptance probability for a Monte Carlo move.
    !----------------------------------------------------------------------
    function compute_acceptance_probability(old, new, res_type, move_type) result(probability)

        ! Input arguments
        type(energy_type), intent(in) :: old   ! Old energy states
        type(energy_type), intent(in) :: new   ! New energy states
        integer, intent(in) :: move_type        ! MC move type (TYPE_CREATION, TYPE_DELETION, TYPE_TRANSLATION, or TYPE_ROTATION)
        integer, intent(in) :: res_type     ! Index of the residue type

        ! Local variables
        real(real64) :: N, Nplus1               ! Number of residues of this type, and Number + 1
        real(real64) :: mu                      ! Chemical potential
        real(real64) :: deltaU                  ! Energy difference ΔE between trial and current state
        real(real64) :: prefactor               ! Prefactor for probability calculation
        real(real64) :: lambda                  ! de Broglie wavelength in m

        ! Return value
        real(real64) :: probability             ! Acceptance probability (0 <= P <= 1)

        N = real(primary%num%residues(res_type), real64)
        Nplus1 = N + 1.0_real64
        deltaU = new%total - old%total                  ! kcal/mol
        mu = thermo%chemical_potential(res_type)     ! kcal/mol

        ! Compute factor based on move type
        select case (move_type)
            case (TYPE_CREATION)

                ! P_acc(N -> N+1) = min[1, (V / ((N+1) λ³)) * exp(-β * (ΔU - μ))]
                ! V in Å³, λ in Å, ΔU and μ in kcal/mol, β = 1/(kB T)
                ! Note: N+1 instead of N to avoid division by zero
                lambda = res%lambda(res_type) ! Thermal de Broglie wavelength (A)
                prefactor = primary%cell%volume / N / lambda**3 ! note: N must be used because the residue was already added
                probability = min(1.0_real64, prefactor * exp(-beta * (deltaU - mu)))

            case (TYPE_DELETION)

                ! P_acc(N -> N-1) = min[1, (N λ³ / V) * exp(-β * (ΔU + μ))]
                ! λ in Å, V in Å³, ΔU and μ in kcal/mol, β = 1/(kB T)
                lambda = res%lambda(res_type) ! Thermal de Broglie wavelength (A)
                prefactor = Nplus1 * lambda**3 / (primary%cell%volume) ! note: Nplus1 must be used because the residue was already removed
                probability = min(1.0_real64, prefactor * exp(-beta * (deltaU + mu)))

            case (TYPE_TRANSLATION, TYPE_ROTATION)

                ! P_acc = min[1, exp(-β * ΔU)]
                probability = min(1.0_real64, exp(-beta * deltaU))

            case default
                call abort_run("Unknown move_type in compute_acceptance_probability!", 1)
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

        ! Arguments
        type(energy_type), intent(in) :: old   ! Energy of system with old molecule
        type(energy_type), intent(in) :: new   ! Energy of system with new molecule
        integer, intent(in) :: type_old         ! Residue type being removed
        integer, intent(in) :: type_new         ! Residue type being inserted

        ! Local variables
        real(real64) :: deltaU                  ! Energy difference
        real(real64) :: mu_old, mu_new          ! Chemical potentials (kcal/mol)
        real(real64) :: N_old, N_new            ! Number of molecules per type
        real(real64) :: Nplus1

        ! Return value
        real(real64) :: probability             ! Acceptance probability (0 <= P <= 1)

        N_new = real(primary%num%residues(type_new), real64)
        N_old = real(primary%num%residues(type_old), real64)
        Nplus1 = N_new + 1.0_real64

        ! Chemical potentials
        mu_old = thermo%chemical_potential(type_old) ! kcal/mol
        mu_new = thermo%chemical_potential(type_new) ! kcal/mol
        deltaU = new%total - old%total         ! kcal/mol

        ! Swap acceptance probability:
        ! P_acc = min[1, (N_old / (N_new + 1)) * exp(-β (ΔE + μ_new - μ_old))]
        probability = min(1.0_real64, (N_old / Nplus1) * exp(-beta * (deltaU + mu_new - mu_old)))

    end function swap_acceptance_probability

    !---------------------------------------------------------------------------
    ! Compute the updated energy of a single molecule after a trial move
    ! for use in the Monte Carlo acceptance test.
    !---------------------------------------------------------------------------
    subroutine compute_new_energy(res_type, mol_index, is_creation, is_deletion)

        ! Input arguments
        integer, intent(in) :: res_type         ! Residue type to be moved
        integer, intent(in) :: mol_index       ! Index of the molecule to move
        logical, intent(in), optional :: is_creation
        logical, intent(in), optional :: is_deletion

        ! Local variables
        logical :: creation_flag
        logical :: deletion_flag

        creation_flag = present_or_false(is_creation)
        deletion_flag = present_or_false(is_deletion)

        if (creation_flag) then
        
            ! Note: In creation scenario, compute all energy components

            call compute_ewald_phase_factors(res_type, mol_index)
            call update_reciprocal_amplitude_single_mol(res_type, mol_index, new%recip_coulomb, is_creation = creation_flag)
            call pairwise_energy_for_molecule(primary, res_type, mol_index, new%non_coulomb, new%coulomb)
            call compute_ewald_self_interaction_single_mol(res_type, new%ewald_self)
            new%intra_coulomb = intra_res_real_coulomb_energy(res_type, mol_index)
        
            ! Recalculate total energy
            new%total = new%non_coulomb + new%coulomb + new%recip_coulomb + new%ewald_self + new%intra_coulomb
        
        else if (deletion_flag) then

            ! Note: Most energy terms in the absence of a molecule are 0
            ! --> One must only recalculate energy%recip_coulomb

            new%non_coulomb = 0
            new%coulomb = 0
            new%ewald_self = 0
            new%intra_coulomb = 0
            call update_reciprocal_amplitude_single_mol(res_type, mol_index, new%recip_coulomb, is_deletion = deletion_flag)

            ! Recalculate total energy
            new%total = new%non_coulomb + new%coulomb + new%recip_coulomb + new%ewald_self + new%intra_coulomb
        
        else

            ! Note, for simple move (translation or rotation), one only needs to
            ! recompute reciprocal and pairwise interactions

            call compute_ewald_phase_factors(res_type, mol_index)
            call update_reciprocal_amplitude_single_mol(res_type, mol_index, new%recip_coulomb)
            call pairwise_energy_for_molecule(primary, res_type, mol_index, new%non_coulomb, new%coulomb)

            ! Recalculate total energy
            new%total = new%non_coulomb + new%coulomb + new%recip_coulomb

        end if

    end subroutine compute_new_energy

    !---------------------------------------------------------------------------
    ! Compute the previous energy of a single molecule before a trial move
    ! for use in the Monte Carlo acceptance test.
    !---------------------------------------------------------------------------
    subroutine compute_old_energy(res_type, mol_index, is_creation, is_deletion)

        ! Input variables
        integer, intent(in) :: res_type         ! Residue type to be moved
        integer, intent(in) :: mol_index       ! Index of the molecule to move
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

            ! #TODO CHECK
            ! call update_reciprocal_amplitude_single_mol(res_type, mol_index, old%recip_coulomb)
            old%recip_coulomb = energy%recip_coulomb

            ! Recalculate total energy
            old%total = old%non_coulomb + old%coulomb + old%recip_coulomb + old%ewald_self + old%intra_coulomb

        else if (deletion_flag) then

            ! Note: In deletion scenario, compute all energy components
            call compute_ewald_self_interaction_single_mol(res_type, old%ewald_self)
            old%intra_coulomb = intra_res_real_coulomb_energy(res_type, mol_index)
            call pairwise_energy_for_molecule(primary, res_type, mol_index, old%non_coulomb, old%coulomb)
            
            ! #TODO CHECK
            ! call update_reciprocal_amplitude_single_mol(res_type, mol_index, old%recip_coulomb)
            old%recip_coulomb = energy%recip_coulomb

            ! Recalculate total energy
            old%total = old%non_coulomb + old%coulomb + old%recip_coulomb + old%ewald_self + old%intra_coulomb

        else

            ! Note, for simple move (translation or rotation), one only needs to
            ! recompute pairwise interactions. Current value for the reciprocal
            ! energy can be used.
            old%recip_coulomb = energy%recip_coulomb
            call pairwise_energy_for_molecule(primary, res_type, mol_index, old%non_coulomb, old%coulomb)

            ! Recalculate total energy
            old%total = old%non_coulomb + old%coulomb + old%recip_coulomb

        end if

    end subroutine compute_old_energy

    !---------------------------------------------------------------------------
    ! Update the global system energy and Monte Carlo counters after
    ! accepting a trial move (translation or rotation) of a molecule.
    !---------------------------------------------------------------------------
    subroutine accept_molecule_move(old, new, counter_var)

        ! Input arguments
        type(energy_type), intent(in) :: old, new  ! Old and new energy states
        integer, intent(inout) :: counter_var(2)    ! Counter for succesfull move

        energy%recip_coulomb = new%recip_coulomb
        energy%non_coulomb = energy%non_coulomb + new%non_coulomb - old%non_coulomb
        energy%coulomb = energy%coulomb + new%coulomb - old%coulomb

        energy%total = energy%total + new%total - old%total
        counter_var(2) = counter_var(2) + 1

    end subroutine accept_molecule_move

    !---------------------------------------------------------------------------
    ! Save the current state of a molecule for a Monte Carlo move. For
    ! translation: saves center-of-mass (COM). For rotation: saves site offsets.
    ! Also saves Fourier terms for later restoration if the move is rejected
    !---------------------------------------------------------------------------
    subroutine save_molecule_state(res_type, mol_index, com_old, offset_old)

        ! Input arguments
        integer, intent(in) :: res_type         ! Residue type of the molecule
        integer, intent(in) :: mol_index        ! Index of the molecule
        real(real64), intent(out), optional :: com_old(3) ! Center-of-mass before move (for translation)
        real(real64), intent(out), optional :: offset_old(:, :) ! Site offsets before move (for rotation)

        ! Local variable
        integer :: natoms

        ! Save Fourier terms
        call save_single_mol_fourier_terms(res_type, mol_index)

        ! Save center-of-mass if requested (translation)
        if (present(com_old)) then
            com_old(:) = guest%com(:, res_type, mol_index)
        end if

        ! Save site offsets if requested (rotation)
        if (present(offset_old)) then
            natoms = res%atom(res_type)
            offset_old(:, 1:natoms) = guest%offset(:, res_type, mol_index, 1:natoms)
        end if

    end subroutine save_molecule_state

    !---------------------------------------------------------------------------
    ! Restore a molecule's previous state (COM or site offsets) and Fourier terms
    ! if a Monte Carlo move is rejected.
    !---------------------------------------------------------------------------
    subroutine reject_molecule_move(res_type, mol_index, com_old, site_offset_old)

        ! Input parameters
        integer, intent(in) :: res_type     ! Residue type of the molecule
        integer, intent(in) :: mol_index    ! Index of the molecule
        real(real64), intent(in), optional :: com_old(3) ! Previous COM (translation)
        real(real64), intent(in), optional :: site_offset_old(:, :) ! Previous site offsets (rotation)

        ! Local variable
        integer :: natoms                   ! Number of atoms in the residue

        ! Restore COM if present (translation)
        if (present(com_old)) then
            guest%com(:, res_type, mol_index) = com_old(:)
        end if

        ! Restore site offsets if present (rotation)
        if (present(site_offset_old)) then
            natoms = res%atom(res_type)
            guest%offset(:, res_type, mol_index, 1:natoms) = site_offset_old(:, 1:natoms)
        end if

        ! Restore Fourier states
        call restore_single_mol_fourier(res_type, mol_index)

    end subroutine reject_molecule_move

    !---------------------------------------------------------------------------
    ! Generates a random position for the new molecule and copies/orients
    ! its atomic geometry. Can take geometry from a reservoir or apply
    ! a random rotation if no reservoir exists.
    !---------------------------------------------------------------------------
    subroutine insert_and_orient_molecule(res_type, molecule_index, reservoir_index, place_random_com)

        ! Input arguments
        integer, intent(in) :: res_type     ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID
        integer, intent(out) :: reservoir_index ! Randomly selected molecule index from the reservoir
        logical, intent(in), optional :: place_random_com ! To control if the COM must be picked

        ! Local variables
        logical :: do_place_com                 ! Internal flag to control COM generation
        logical :: full_rotation                ! Flag indicating whether a full 360° random rotation should be applied
        real(real64) :: random_nmb              ! Uniform random number in [0,1), used for random index selection
        real(real64) :: trial_pos(3)            ! Random numbers for initial molecule position in the box

        ! Decide whether to place a random COM
        if (present(place_random_com)) then
            do_place_com = place_random_com
        else
            do_place_com = .true. ! Default behavior
        end if

        ! Generate a random position in the simulation box (if enabled)
        if (do_place_com) then
            call random_number(trial_pos) ! Random numbers in [0,1)
            guest%com(:, res_type, molecule_index) = primary%cell%bounds(:,1) &
                + matmul(primary%cell%matrix, trial_pos)
        end if

        ! Copy geometry from reservoir or rotate if no reservoir
        if (status%reservoir_provided) then

            ! Pick a random (and existing) molecule in the reservoir
            call random_number(random_nmb)

            if (present(place_random_com)) then
                reservoir_index = int(random_nmb * reservoir%num%residues(res_type)) + 1  
                ! Copy site offsets from the chosen molecule
                guest%offset(:, res_type, molecule_index, 1:res%atom(res_type)) = &
                    gas%offset(:, res_type, reservoir_index, 1:res%atom(res_type))
            else
                reservoir_index = int(random_nmb * reservoir%num%residues(res_type)) + 1  
                ! Copy site offsets from the chosen molecule
                guest%offset(:, res_type, molecule_index, 1:res%atom(res_type)) = &
                    gas%offset(:, res_type, reservoir_index, 1:res%atom(res_type))
            end if  

        else

            ! No reservoir provided
            reservoir_index = 0

            ! Copy site offsets from the first molecule
            guest%offset(:, res_type, molecule_index, 1:res%atom(res_type)) = &
                guest%offset(:, res_type, 1, 1:res%atom(res_type))

            ! Rotate the new molecule randomly (using full 360° rotation)
            full_rotation = .true.
            call apply_random_rotation(res_type, molecule_index, full_rotation)

        end if

    end subroutine insert_and_orient_molecule

    !---------------------------------------------------------------------------
    ! Restores molecule and atom counts, and resets Fourier states if a
    ! creation move is rejected.
    !---------------------------------------------------------------------------
    subroutine reject_creation_move(res_type, mol_index)

        ! Input parameters
        integer, intent(in) :: res_type     ! Residue type to be moved
        integer, intent(in) :: mol_index   ! Molecule ID

        ! Restore previous residue/atom numbers
        call update_counts(primary, res_type, -1)

        ! Restore Fourier states (ik_alloc and dk_alloc, all zeros)
        call restore_single_mol_fourier(res_type, mol_index)

    end subroutine reject_creation_move

    !------------------------------------------------------------------------------
    ! Compute the excess chemical potential (μ_ex) for each residue type
    ! using Widom particle insertion method, and optionally the ideal chemical
    ! potential (μ_ideal) for reporting purposes.
    !------------------------------------------------------------------------------
    subroutine calculate_excess_mu()

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
        do type_residue = 1, res%number
            if (statistic%sample(type_residue) > 0) then

                ! Compute average Boltzmann factor from Widom sampling
                ! <exp(-β ΔU)> = sum_weights / N_samples
                avg_weight = statistic%weight(type_residue) / real(statistic%sample(type_residue), kind=real64)

                ! Compute excess chemical potential (kcal/mol)
                ! μ_ex = - k_B * T * ln(<exp(-β ΔU)>)
                temperature = thermo%temperature
                statistic%mu_ex(type_residue) = - KB_kcalmol * temperature * log(avg_weight) ! kcal/mol

                ! Compute ideal gas chemical potential (kcal/mol)
                ! μ_ideal = k_B * T * ln(ρ * Λ^3)
                lambda = res%lambda(type_residue)                   ! Thermal de Broglie wavelength (A)
                N = primary%num%residues(type_residue)              ! Number of molecules
                volume = primary%cell%volume                        ! Box volume (A^3)
                rho = real(N, kind=real64) / volume                 ! Number density (molecules/A^3)
                mu_ideal = KB_kcalmol * temperature * log(rho * lambda**3) ! Ideal chemical potential (kcal/mol)

                statistic%mu_tot(type_residue) = mu_ideal + statistic%mu_ex(type_residue)

            end if
        end do

    end subroutine calculate_excess_mu

    !---------------------------------------------------------------------------
    ! Physically removes a molecule from the simulation box by replacing it
    ! with the last molecule in the array and updating Fourier terms.
    !---------------------------------------------------------------------------
    subroutine remove_molecule(res_type, molecule_index, last_molecule_index)

        integer, intent(in) :: res_type      ! Residue type to remove
        integer, intent(in) :: molecule_index    ! Molecule index to remove
        integer, intent(in):: last_molecule_index ! Index of the last molecule in the primary box

        ! Replace with the last molecule
        guest%com(:, res_type, molecule_index) = &
            guest%com(:, res_type, last_molecule_index)
        guest%offset(:, res_type, molecule_index, 1:res%atom(res_type)) = &
            guest%offset(:, res_type, last_molecule_index, 1:res%atom(res_type))

        ! Replace Fourier terms
        call replace_fourier_terms_single_mol(res_type, molecule_index, last_molecule_index)

    end subroutine remove_molecule

    !---------------------------------------------------------------------------
    ! Updates residue and atom counters with sign = +1 or -1
    !---------------------------------------------------------------------------
    subroutine update_counts(box, res_type, sign)

        ! Input parameters
        type(type_box), intent(inout) :: box
        integer, intent(in) :: res_type
        integer, intent(in) :: sign

        box%num%residues(res_type) = box%num%residues(res_type) + sign
        box%num%atoms = box%num%atoms + sign * res%atom(res_type)

    end subroutine update_counts

end module monte_carlo_utils
