module molecule_swap

    !===========================================================================
    ! Module: molecule_swap
    !
    ! Purpose:
    !   Implements Monte Carlo identity–swap (transmutation) moves between 
    !   different molecular residue types. A selected molecule of type A is 
    !   converted into a molecule of type B while preserving its center-of-mass 
    !   position.
    !
    !======

    use monte_carlo_utils
    use ewald_kvectors
    use ewald_phase
    use ewald_energy
    use simulation_state
    use random_utils
    use energy_utils
    use geometry_utils
    use constants
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    private ! Hide everything by default
    public :: SwapMolecules  ! Only expose the main entry point

contains

    !---------------------------------------------------------------------------
    ! Subroutine: SwapMolecules
    !
    !---------------------------------------------------------------------------
    subroutine SwapMolecules(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type     ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID

        ! Local variables
        integer :: residue_type_bis     ! Residue type to be swapped
        integer :: molecule_index_bis   ! Index of molecule to be swapped
        integer :: rand_mol_index       ! Randomly selected molecule index from the reservoir for copying geometry
        real(real64) :: probability     ! Acceptance probability of creation move
        integer :: last_molecule_index  ! Index of the last molecule in the primary box
        logical :: is_deletion          ! Flag indicating creation
        logical :: is_creation          ! Flag indicating creation

        ! Pick randomly a second residue type
        residue_type_bis = PickDifferentResidueType(residue_type)

        ! If no valid different residue type was found, skip this move
        if (residue_type_bis == -1) return

        ! Check that there is at least one molecule of the selected type to swap
        if (primary%num_residues(residue_type_bis)==0) return ! #todo : Is this really necessary ?

        ! Count trial move (success + fail)
        counter%trial_swap = counter%trial_swap + 1

        ! Pick a molecule ID for the second type
        molecule_index_bis = primary%num_residues(residue_type_bis) + 1
        
        ! STEP 1 - Delete a molecule

        ! Energy of the previous configuration
        is_deletion = .true.
        call ComputeOldEnergy(residue_type, molecule_index, old, is_deletion = is_deletion)
        call SaveMoleculeState(residue_type, molecule_index, com_old = res%mol_com_old, offset_old = res%site_offset_old)

        ! Record the index of the last molecule of type "residue_type"
        last_molecule_index = primary%num_residues(residue_type)

        ! Delete molecule
        call RemoveMolecule(residue_type, molecule_index, last_molecule_index)

        ! Update molecule and atom counts
        primary%num_residues(residue_type) = primary%num_residues(residue_type) - 1
        primary%num_atoms = primary%num_atoms - nb%atom_in_residue(residue_type)

        ! STEP 2 - Place a new molecule at the same location

        ! Update molecule and atom counts (second time)
        primary%num_residues(residue_type_bis) = primary%num_residues(residue_type_bis) + 1
        primary%num_atoms = primary%num_atoms + nb%atom_in_residue(residue_type_bis)

        ! Use the CoM of the deleted molecule
        primary%mol_com(:, residue_type_bis, molecule_index_bis) = res%mol_com_old

        ! Generate or pick orientation for the new molecule
        call OrientMolecule(residue_type_bis, molecule_index_bis, rand_mol_index)

        ! Compute new energy
        is_creation = .true.
        call ComputeNewEnergy(residue_type_bis, molecule_index_bis, new, is_creation = is_creation)

        ! STEP 3 - Accept or reject move

        ! Compute acceptance probability for the move
        probability = swap_acceptance_probability(old, new, residue_type, residue_type_bis)

        ! Accept or reject
        if (rand_uniform() <= probability) then ! Accept move
            call AcceptSwapMove(old, new)
        else ! Reject move
            call RejectSwapMove(residue_type, molecule_index, residue_type_bis, res%mol_com_old, res%site_offset_old)
        end if

    end subroutine SwapMolecules

    subroutine RejectSwapMove(residue_type, molecule_index, residue_type_bis, mol_com_old, site_offset_old)

        implicit none

        integer, intent(in) :: residue_type, residue_type_bis       ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID
        real(real64), dimension(3) :: mol_com_old ! For storing old molecule center-of-mass
        real(real64), dimension(:, :) :: site_offset_old

        ! Restore previous residue/atom numbers
        primary%num_atoms = primary%num_atoms + nb%atom_in_residue(residue_type)
        primary%num_residues(residue_type) = primary%num_residues(residue_type) + 1
        primary%num_atoms = primary%num_atoms - nb%atom_in_residue(residue_type_bis)
        primary%num_residues(residue_type_bis) = primary%num_residues(residue_type_bis) - 1

        ! Restore previous positions and orientation
        primary%mol_com(:, residue_type, molecule_index) = mol_com_old(:)
        primary%site_offset(:, residue_type, molecule_index, 1:nb%atom_in_residue(residue_type)) = &
            site_offset_old(:, 1:nb%atom_in_residue(residue_type))

        ! Restore Fourier states (ik_alloc and dk_alloc, all zeros)
        call RestoreSingleMolFourier(residue_type, molecule_index)

    end subroutine RejectSwapMove

    subroutine AcceptSwapMove(old, new)

        implicit none

        type(energy_state), intent(in) :: old   ! Previous energy states
        type(energy_state), intent(in) :: new   ! New energy states

        energy%recip_coulomb    = new%recip_coulomb
        energy%non_coulomb      = energy%non_coulomb    + new%non_coulomb   - old%non_coulomb
        energy%coulomb          = energy%coulomb        + new%coulomb       - old%coulomb
        energy%ewald_self       = energy%ewald_self     + new%ewald_self    - old%ewald_self
        energy%intra_coulomb    = energy%intra_coulomb  + new%intra_coulomb - old%intra_coulomb
        energy%total            = energy%total          + new%total         - old%total 

        ! Count successful move
        counter%swap = counter%swap + 1

    end subroutine AcceptSwapMove

    !---------------------------------------------------------------------------
    ! Function: PickDifferentResidueType
    !
    ! Purpose:
    !   Attempts to select a residue type different from 'current_type'. The
    !   selection is made among the residue types marked as active in 
    !   input%is_active. A maximum number of attempts is performed to avoid 
    !   infinite loops.
    !
    ! Returns:
    !   new_type  - A residue type different from current_type, or -1 on failure.
    !---------------------------------------------------------------------------
    function PickDifferentResidueType(current_type, max_attempts) result(new_type)

        implicit none
    
        ! Input arguments
        integer, intent(in) :: current_type        ! Original residue type to avoid
        integer, intent(in), optional :: max_attempts ! Max attempts to pick a different type
        
        ! Local variables
        integer :: new_type                         ! Selected residue type (or -1 on failure)
        integer :: attempt                          ! Current attempt counter
        integer :: n_attempts                       ! Maximum number of attempts

        ! Set maximum attempts (default = 10)
        n_attempts = 10
        if (present(max_attempts)) n_attempts = max_attempts

        ! Initialize
        new_type = -1
        attempt  = 0

        ! Try to pick a different type
        do while (attempt < n_attempts)
            new_type = PickRandomResidueType(input%is_active)

            if (new_type /= current_type) then
                return
            end if

            attempt = attempt + 1
        end do

        ! If all attempts failed, return failure code
        new_type = -1

    end function PickDifferentResidueType

    subroutine RemoveMolecule(residue_type, molecule_index, last_molecule_index)

        implicit none

        integer, intent(in) :: residue_type      ! Residue type to remove
        integer, intent(in) :: molecule_index    ! Molecule index to remove
        integer, intent(in):: last_molecule_index ! Index of the last molecule in the primary box

        ! Replace with the last molecule
        primary%mol_com(:, residue_type, molecule_index) = &
            primary%mol_com(:, residue_type, last_molecule_index)
        primary%site_offset(:, residue_type, molecule_index, 1:nb%atom_in_residue(residue_type)) = &
            primary%site_offset(:, residue_type, last_molecule_index, 1:nb%atom_in_residue(residue_type))

        ! Replace Fourier terms
        call ReplaceFourierTermsSingleMol(residue_type, molecule_index, last_molecule_index)

    end subroutine RemoveMolecule

    subroutine OrientMolecule(residue_type, molecule_index, rand_mol_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type     ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID
        integer :: rand_mol_index               ! Randomly selected molecule index from the reservoir for copying geometry

        ! Local variables
        logical :: full_rotation                ! Flag indicating whether a full 360° random rotation should be applied
        real(real64) :: random_nmb              ! Uniform random number in [0,1), used for random index selection

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
            full_rotation = .True.
            call ApplyRandomRotation(residue_type, molecule_index, full_rotation)

        end if

    end subroutine OrientMolecule

end module molecule_swap
