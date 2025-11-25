module molecule_translation

    use monte_carlo_utils
    use simulation_state
    use random_utils
    use energy_utils
    use ewald_kvectors
    use ewald_phase
    use ewald_energy
    use constants
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !---------------------------------------------------------------------------
    ! Perform a trial translation of a molecule by a random displacement
    ! vector.
    !---------------------------------------------------------------------------
    subroutine attempt_translation_move(residue_type, molecule_index)

        ! Input arguments
        integer, intent(in) :: residue_type         ! Residue type to be moved
        integer, intent(in) :: molecule_index       ! Molecule ID

        ! Local variables
        real(real64), dimension(3) :: com_old       ! Previous X, Y, Z coordinate of molecule
        real(real64) :: probability                 ! Acceptance probability of move

        ! Exit early if molecule index is zero
        if (molecule_index == 0) return

        ! Increment trial move counter
        counter%translations(1) = counter%translations(1) + 1

        ! Save the current state of the molecule
        call save_molecule_state(residue_type, molecule_index, com_old = com_old)

        ! Compute old energy of the molecule/system
        call compute_old_energy(residue_type, molecule_index)

        ! Propose a random translation move
        call propose_translation_move(residue_type, molecule_index)

        ! Compute new energy after the proposed move
        call compute_new_energy(residue_type, molecule_index)

        ! Compute Metropolis acceptance probability
        ! probability = min(1, exp(-(new - old)/kT))
        probability = compute_acceptance_probability(old, new, residue_type, TYPE_TRANSLATION)

        ! Accept or reject
        if (rand_uniform() <= probability) then
            ! Accept move: update system state
            call accept_molecule_move(old, new, counter%translations)
        else
            ! Reject move: restore previous position
            call reject_molecule_move(residue_type, molecule_index, com_old = com_old)
        end if

    end subroutine attempt_translation_move

    !-----------------------------------------------------------------------
    ! Generate a trial translation move for a selected molecule of a given
    ! residue type. A random displacement vector is drawn uniformly from
    ! [-translation_step/2, translation_step/2) along each coordinate axis,
    ! added to the molecule's previous center of mass, and wrapped into
    ! the simulation box using the minimum image convention.
    !-----------------------------------------------------------------------
    subroutine propose_translation_move(res_type, mol_index)

        ! Input variables
        integer, intent(in) :: res_type             ! Residue type to be moved
        integer, intent(in) :: mol_index            ! Molecule ID
        
        ! Local variables
        real(real64) :: trial_pos(3)
        type(type_coordinate), pointer :: coord     ! Pointer for host or guest coordinate

        ! Generate random move of max size "translation_step/2"
        trial_pos = rand_symmetric(3) * input%translation_step

        ! Return the correct pointer (host or guest)
        coord => get_coord(res_type)

        ! Apply translation to previous COM position
        coord%com(:, res_type, mol_index) = coord%com(:, res_type, mol_index) + trial_pos(:)

        ! Apply minimum image convension
        call apply_PBC(coord%com(:, res_type, mol_index), primary)

    end subroutine propose_translation_move

end module molecule_translation
