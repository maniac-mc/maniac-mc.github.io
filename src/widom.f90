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

        ! Input arguments
        integer, intent(in) :: residue_type     ! Residue type to be inserted
        integer, intent(in) :: molecule_index   ! Molecule ID

        ! Local variables
        real(real64) :: deltaU                  ! Energy difference

        call check_molecule_index(molecule_index)

        ! Count trial move
        counter%widom(1) = counter%widom(1) + 1

        ! Compute old energy
        call compute_old_energy(residue_type, molecule_index, is_creation = .true.)

        ! Increase the residue and atom counts
        primary%num_residues(residue_type) = primary%num_residues(residue_type) + 1
        primary%num_atoms = primary%num_atoms + nb%atom_in_residue(residue_type)

        ! Save current Fourier terms (should be all zeros here)
        call save_single_mol_fourier_terms(residue_type, molecule_index)

        ! Generate random insertion position within the simulation box
        call insert_and_orient_molecule(residue_type, molecule_index)

        ! Compute new energy
        call compute_new_energy(residue_type, molecule_index, is_creation = .true.)

        ! Reject systematically (the system state is preserved)
        call reject_creation_move(residue_type, molecule_index)

        ! Compute Boltzmann weight for Widom insertion
        deltaU = new%total - old%total
        call accumulate_widom_weight(residue_type, deltaU)

    end subroutine widom_trial

    !---------------------------------------------------------------------------
    ! Computes the Widom Boltzmann weight w = exp(-ΔU / (kB*T)) and updates
    ! accumulated statistics for the given residue type.
    !---------------------------------------------------------------------------
    subroutine accumulate_widom_weight(residue_type, deltaU)
        
        ! Input parameters
        integer, intent(in) :: residue_type     ! Residue type to be inserted
        real(real64), intent(in) :: deltaU      ! Energy difference

        ! Local variables
        real(real64) :: weight                  ! Boltzmann weight

        weight = exp(-deltaU * beta)

        if (weight > error) then ! Correspond to a success
            counter%widom(2) = counter%widom(2) + 1
        end if

        statistic%weight(residue_type) = statistic%weight(residue_type) + weight
        statistic%sample(residue_type) = statistic%sample(residue_type) + 1

    end subroutine accumulate_widom_weight

end module widom
