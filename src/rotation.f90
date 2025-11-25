module molecule_rotation

    use monte_carlo_utils
    use simulation_state
    use ewald_kvectors
    use ewald_phase
    use ewald_energy
    use random_utils
    use energy_utils
    use constants
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !---------------------------------------------------------------------------    !
    ! Performs a trial rotation of a molecule around a random axis by a random
    ! angle. The move is accepted or rejected based on the Metropolis
    ! criterion using energy differences.
    !---------------------------------------------------------------------------
    subroutine attempt_rotation_move(residue_type, molecule_index)

        ! Input arguments
        integer, intent(in) :: residue_type             ! Residue type to be moved
        integer, intent(in) :: molecule_index           ! Molecule ID

        ! Local variables
        real(real64) :: probability                     ! Acceptance probability of rotation move

        ! Exit early if molecule cannot rotate
        if ((nb%atom_in_residue(residue_type) == 1) .or. (molecule_index == 0)) return

        ! Count trial move
        counter%rotations(1) = counter%rotations(1) + 1

        call save_molecule_state(residue_type, molecule_index, offset_old = res%site_offset_old)

        ! Compute previous energy
        call compute_old_energy(residue_type, molecule_index)

        ! Rotate the molecule randomly
        call apply_random_rotation(residue_type, molecule_index)

        ! Recompute energies
        call compute_new_energy(residue_type, molecule_index)

        ! Compute acceptance probability for the move
        probability = compute_acceptance_probability(old, new, residue_type, TYPE_ROTATION)

        ! Accept or reject
        if (rand_uniform() <= probability) then
            ! Accept move: update system state
            call accept_molecule_move(old, new, counter%rotations)
        else ! Reject move
            ! Reject move: restore previous orientation
            call reject_molecule_move(residue_type, molecule_index, site_offset_old = res%site_offset_old)
        end if

    end subroutine attempt_rotation_move

end module molecule_rotation
