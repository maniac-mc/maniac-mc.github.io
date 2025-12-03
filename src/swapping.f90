module molecule_swap

    !===========================================================================
    ! Implements Monte Carlo identity–swap (transmutation) moves between 
    ! different molecular residue types. A selected molecule of type A is 
    ! converted into a molecule of type B while preserving its center-of-mass 
    ! position.
    !===========================================================================

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
    
contains

    !---------------------------------------------------------------------------
    ! Attempts an identity-swap (transmutation) move between two molecular
    ! residue types. A molecule of type 'residue_type' is deleted while a new
    ! molecule of a different type is inserted at the same center-of-mass
    ! position. The routine evaluates the old and new energies, constructs the
    ! trial configuration, computes the Metropolis acceptance probability, and
    ! finally accepts or rejects the swap while restoring or updating all state
    ! variables accordingly.
    !---------------------------------------------------------------------------
    subroutine attempt_swap_move(residue_type, molecule_index)

        ! Input arguments
        integer, intent(in) :: residue_type     ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID

        ! Local variables
        integer :: residue_type_bis     ! Residue type to be swapped
        integer :: molecule_index_bis   ! Index of molecule to be swapped
        integer :: rand_mol_index       ! Randomly selected molecule index from the reservoir for copying geometry
        real(real64) :: probability     ! Acceptance probability of creation move
        integer :: last_molecule_index  ! Index of the last molecule in the primary box

        ! Pick randomly a second residue type
        residue_type_bis = pick_different_residue_type(residue_type)

        ! If no valid different residue type was found, skip this move
        if (residue_type_bis == -1) return

        ! Check that there is at least one molecule of the selected type to swap
        if (primary%num%residues(residue_type_bis)==0) return ! #todo : Is this really necessary ?

        ! Count trial move
        counter%swaps(1) = counter%swaps(1) + 1

        ! Pick a molecule ID for the second type
        molecule_index_bis = primary%num%residues(residue_type_bis) + 1
        
        ! STEP 1 - Delete a molecule

        ! Energy of the previous configuration
        call compute_old_energy(residue_type, molecule_index, is_deletion = .true.)
        call save_molecule_state(residue_type, molecule_index, com_old = saved%com, offset_old = saved%offset)

        ! Record the index of the last molecule of type "residue_type"
        last_molecule_index = primary%num%residues(residue_type)

        ! Delete molecule
        call remove_molecule(residue_type, molecule_index, last_molecule_index)

        ! Update molecule and atom counts
        primary%num%residues(residue_type) = primary%num%residues(residue_type) - 1
        primary%num%atoms = primary%num%atoms - res%atom(residue_type)

        ! STEP 2 - Place a new molecule at the same location

        ! Update molecule and atom counts (second time)
        primary%num%residues(residue_type_bis) = primary%num%residues(residue_type_bis) + 1
        primary%num%atoms = primary%num%atoms + res%atom(residue_type_bis)

        ! Use the CoM of the deleted molecule
        guest%com(:, residue_type_bis, molecule_index_bis) = saved%com

        ! Generate or pick orientation for the new molecule
        call insert_and_orient_molecule(residue_type_bis, molecule_index_bis, rand_mol_index, place_random_com = .false.)

        ! Compute new energy
        call compute_new_energy(residue_type_bis, molecule_index_bis, is_creation = .true.)

        ! STEP 3 - Accept or reject move

        ! Compute acceptance probability for the move
        probability = swap_acceptance_probability(old, new, residue_type, residue_type_bis)

        ! Accept or reject
        if (rand_uniform() <= probability) then ! Accept move
            call accept_swap_move()
        else ! Reject move
            call reject_swap_move(residue_type, molecule_index, residue_type_bis, saved%com, saved%offset)
        end if

    end subroutine attempt_swap_move

    !---------------------------------------------------------------------------
    ! Restores the system to its state prior to a rejected swap move. The
    ! subroutine reinstates the original residue and atom counts, restores the
    ! molecule’s center-of-mass and internal geometry, and rebuilds the associated
    ! Fourier (reciprocal-space) contributions.
    !---------------------------------------------------------------------------
    subroutine reject_swap_move(residue_type, molecule_index, residue_type_bis, mol_com_old, site_offset_old)

        ! Input parameters
        integer, intent(in) :: residue_type, residue_type_bis       ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID
        real(real64), dimension(3) :: mol_com_old ! For storing old molecule center-of-mass
        real(real64), dimension(:, :) :: site_offset_old

        ! Restore previous residue/atom numbers
        primary%num%atoms = primary%num%atoms + res%atom(residue_type)
        primary%num%residues(residue_type) = primary%num%residues(residue_type) + 1

        primary%num%atoms = primary%num%atoms - res%atom(residue_type_bis)
        primary%num%residues(residue_type_bis) = primary%num%residues(residue_type_bis) - 1

        ! Restore previous positions and orientation
        guest%com(:, residue_type, molecule_index) = mol_com_old(:)
        guest%offset(:, residue_type, molecule_index, 1:res%atom(residue_type)) = &
            site_offset_old(:, 1:res%atom(residue_type))

        ! Restore Fourier states (ik_alloc and dk_alloc, all zeros)
        call restore_single_mol_fourier(residue_type, molecule_index)

    end subroutine reject_swap_move

    !---------------------------------------------------------------------------
    ! Commits an accepted swap move by updating global energy components with the
    ! differences between the new and old configurations, and increments the
    ! counter tracking successful identity–swap moves.
    !---------------------------------------------------------------------------
    subroutine accept_swap_move()

        energy%recip_coulomb    = new%recip_coulomb
        energy%non_coulomb      = energy%non_coulomb    + new%non_coulomb   - old%non_coulomb
        energy%coulomb          = energy%coulomb        + new%coulomb       - old%coulomb
        energy%ewald_self       = energy%ewald_self     + new%ewald_self    - old%ewald_self
        energy%intra_coulomb    = energy%intra_coulomb  + new%intra_coulomb - old%intra_coulomb
        energy%total            = energy%total          + new%total         - old%total 

        ! Count successful move
        counter%swaps = counter%swaps + 1

    end subroutine accept_swap_move

    !---------------------------------------------------------------------------
    ! Attempts to select a residue type different from 'current_type'. The
    ! selection is made among the residue types marked as active in 
    ! thermo%is_active. A maximum number of attempts is performed to avoid 
    ! infinite loops.
    !---------------------------------------------------------------------------
    function pick_different_residue_type(current_type, max_attempts) result(new_type)
    
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
        
            new_type = pick_random_residue_type(thermo%is_active)

            if (new_type /= current_type) then
                return
            end if

            attempt = attempt + 1
        end do

        ! If all attempts failed, return failure code
        new_type = -1

    end function pick_different_residue_type

end module molecule_swap
