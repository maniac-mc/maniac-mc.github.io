module widom

    !===========================================================================
    ! Implements Widom particle insertion (trial molecule creation) moves for
    ! Monte Carlo simulations. The module generates random positions and 
    ! orientations for test particles, computes their energy contributions, 
    ! and accumulates statistics for chemical potential calculations.
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
    ! Attempts a Widom trial insertion of a molecule of the given residue type.
    ! Computes energy contributions and updates Widom statistics. The move is 
    ! systematically rejected to avoid modifying the system state.
    !---------------------------------------------------------------------------
    subroutine widom_trial(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type     ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID

        ! Local variables
        real(real64) :: deltaU                  ! Energy difference
        real(real64) :: weight                  ! Bolzmann weight

        call check_molecule_index(molecule_index)

        ! Count trial move (success + fail)
        counter%trial_widom = counter%trial_widom + 1

        ! Compute old energy
        call compute_old_energy(residue_type, molecule_index, is_creation = .true.)

        ! Increase the residue and atom counts
        primary%num_residues(residue_type) = primary%num_residues(residue_type) + 1
        primary%num_atoms = primary%num_atoms + nb%atom_in_residue(residue_type)

        ! Save current Fourier terms (should be all zeros here)
        call SaveSingleMolFourierTerms(residue_type, molecule_index)

        ! Generate random insertion position within the simulation box
        call insert_and_orient_molecule(residue_type, molecule_index)

        ! Compute new energy
        call ComputeNewEnergy(residue_type, molecule_index, new, is_creation = .true.)

        ! Reject systematically (the system state is preserved)
        call reject_creation_move(residue_type, molecule_index)

        ! Compute Boltzmann weight for Widom insertion
        ! Formula: = exp(-ΔU / (k_B * T)), beta = 1/(kB T)
        ! where ΔU = U_new - U_old, kB is the Boltzmann constant, and T is the temperature
        deltaU = new%total - old%total ! In K units
        weight = exp(-deltaU * beta)

        widom_stat%weight(residue_type) = widom_stat%weight(residue_type) + weight
        widom_stat%sample(residue_type) = widom_stat%sample(residue_type) + 1

    end subroutine widom_trial

end module widom
