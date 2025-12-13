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
    subroutine widom_trial(res_type, mol_index)

        ! Input arguments
        integer, intent(in) :: res_type     ! Residue type to be inserted
        integer, intent(in) :: mol_index   ! Molecule ID

        ! Local variables
        real(real64) :: deltaU                  ! Energy difference
        integer :: rand_mol_index               ! Randomly selected molecule index from the reservoir for copying geometry

        call check_molecule_index(mol_index)

        ! Count trial move
        counter%widom(1) = counter%widom(1) + 1

        ! Compute old energy
        call compute_old_energy(res_type, mol_index, is_creation = .true.)

        ! Increase the residue and atom counts
        call update_counts(primary, res_type, +1)

        ! Save current Fourier terms
        call save_single_mol_fourier_terms(res_type, mol_index)

        ! Generate random insertion position within the simulation box
        call insert_and_orient_molecule(res_type, mol_index, rand_mol_index)

        ! Compute new energy
        call compute_new_energy(res_type, mol_index, is_creation = .true.)

        ! Reject systematically (the system state is preserved during Widom trial)
        call reject_creation_move(res_type, mol_index)

        ! Compute Boltzmann weight for Widom insertion
        deltaU = new%total - old%total

        call accumulate_widom_weight(res_type, deltaU)

    end subroutine widom_trial

    !---------------------------------------------------------------------------
    ! Computes the Widom Boltzmann weight w = exp(-ΔU / (kB*T)) and updates
    ! accumulated statistics for the given residue type.
    !---------------------------------------------------------------------------
    subroutine accumulate_widom_weight(res_type, deltaU)
        
        ! Input parameters
        integer, intent(in) :: res_type     ! Residue type to be inserted
        real(real64), intent(in) :: deltaU      ! Energy difference

        ! Local variables
        real(real64) :: weight                  ! Boltzmann weight

        weight = exp(-deltaU * beta)

        if (weight > error) then ! Correspond to a success
            counter%widom(2) = counter%widom(2) + 1
            statistic%weight(res_type) = statistic%weight(res_type) + weight
        end if

        statistic%sample(res_type) = statistic%sample(res_type) + 1

    end subroutine accumulate_widom_weight

end module widom
