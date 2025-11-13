module molecule_swap

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
        real(real64), dimension(:, :), allocatable :: site_offset_old ! For storing old site offsets
        real(real64), dimension(3) :: mol_com_old ! Mol com old
        integer :: last_molecule_index  ! Index of the last molecule in the primary box

        ! Pick randomly a second residue type
        residue_type_bis = PickDifferentResidueType(residue_type)

        ! If we couldn't find a different residue type, skip this move
        if (residue_type_bis == residue_type) return

        ! Lack of molecule for swapping
        if (primary%num_residues(residue_type_bis)==0) return

        ! Count trial move (success + fail)
        counter%trial_swap = counter%trial_swap + 1

        ! Pick a molecule ID for the second type
        molecule_index_bis = primary%num_residues(residue_type_bis) + 1
        
        allocate(site_offset_old(3, nb%max_atom_in_residue))

        ! STEP 1 - Delete a molecule

        ! Compute old energies
        call ComputeOldEnergy_swap(residue_type, molecule_index, old)

        ! Save molecule for aborted move
        call SaveMoleculeState_swap(residue_type, molecule_index, mol_com_old, site_offset_old)

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
        primary%mol_com(:, residue_type_bis, molecule_index_bis) = mol_com_old

        ! Generate random orientation for the new molecule
        call OrientMolecule(residue_type_bis, molecule_index_bis, rand_mol_index)

        ! Compute new energy
        call ComputeNewEnergy_swap(residue_type_bis, molecule_index_bis, new)

        probability = mc_acceptance_probability_swap(old, new, residue_type, residue_type_bis)

        ! Accept or reject
        if (rand_uniform() <= probability) then ! Accept move
            call AcceptSwapMove(old, new)
        else ! Reject move
            call RejectSwapMove(residue_type, molecule_index, residue_type_bis, mol_com_old, site_offset_old)
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
        call RestoreFourierState_singlemol(residue_type, molecule_index)

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
    ! Returns a residue type different from 'current_type'
    !---------------------------------------------------------------------------
    function PickDifferentResidueType(current_type, max_attempts) result(new_type)

        implicit none
    
        integer, intent(in) :: current_type
        integer, intent(in), optional :: max_attempts
        integer :: new_type
        integer :: attempt, n_attempts

        n_attempts = 10
        if (present(max_attempts)) n_attempts = max_attempts

        new_type = current_type
        attempt = 0
        do while (new_type == current_type .and. attempt < n_attempts)
            new_type = PickRandomResidueType(input%is_active)
            attempt = attempt + 1
        end do
    end function PickDifferentResidueType

    subroutine ComputeOldEnergy_swap(residue_type, molecule_index, old)

        implicit none

        type(energy_state), intent(inout) :: old    ! old energy states
        integer, intent(in) :: residue_type         ! Residue type to be moved
        integer, intent(in) :: molecule_index       ! Molecule ID

        call ComputeEwaldSelfInteraction_singlemol(residue_type, old%ewald_self)
        call ComputeIntraResidueRealCoulombEnergy_singlemol(residue_type, molecule_index, old%intra_coulomb)
        call ComputePairInteractionEnergy_singlemol(primary, residue_type, molecule_index, old%non_coulomb, old%coulomb)
        old%recip_coulomb = energy%recip_coulomb
        old%total = old%non_coulomb + old%coulomb + old%recip_coulomb + old%ewald_self + old%intra_coulomb

    end subroutine ComputeOldEnergy_swap

    subroutine SaveMoleculeState_swap(residue_type, molecule_index, mol_com_old, site_offset_old)
    
        implicit none

        integer, intent(in) :: residue_type      ! Residue type to remove
        integer, intent(in) :: molecule_index    ! Molecule index to remove
        real(real64), intent(out) :: mol_com_old(3) ! For storing old molecule center-of-mass
        real(real64), intent(out), dimension(:, :) :: site_offset_old

        ! Store positions and site offsets
        mol_com_old(:) = primary%mol_com(:, residue_type, molecule_index)
        site_offset_old(:, 1:nb%atom_in_residue(residue_type)) = &
            primary%site_offset(:, residue_type, molecule_index, 1:nb%atom_in_residue(residue_type))

        ! Save Fourier terms
        call SaveFourierTerms_singlemol(residue_type, molecule_index)

    end subroutine SaveMoleculeState_swap

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
        call ReplaceFourierTerms_singlemol(residue_type, molecule_index, last_molecule_index)

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

    subroutine ComputeNewEnergy_swap(residue_type, molecule_index, new)

        implicit none

        type(energy_state), intent(inout) :: new    ! New energy states
        integer, intent(in) :: residue_type         ! Residue type to be moved
        integer, intent(in) :: molecule_index       ! Molecule ID

        ! Compute Fourier terms for new molecules
        call ComputeFourierTerms_singlemol(residue_type, molecule_index)

        ! Compute new energy
        call UpdateReciprocalEnergy_creation(residue_type, molecule_index, new%recip_coulomb)
        call ComputePairInteractionEnergy_singlemol(primary, residue_type, molecule_index, new%non_coulomb, new%coulomb)
        call ComputeEwaldSelfInteraction_singlemol(residue_type, new%ewald_self)
        call ComputeIntraResidueRealCoulombEnergy_singlemol(residue_type, molecule_index, new%intra_coulomb)
        new%total = new%non_coulomb + new%coulomb + new%recip_coulomb + new%ewald_self + new%intra_coulomb

    end subroutine ComputeNewEnergy_swap

end module molecule_swap
