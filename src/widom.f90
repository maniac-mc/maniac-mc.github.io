module widom

    !===========================================================================
    ! Module: widom
    !
    ! Purpose:
    !   Implements Widom particle insertion (trial molecule creation) moves for
    !   Monte Carlo simulations. The module generates random positions and 
    !   orientations for test particles, computes their energy contributions, 
    !   and accumulates statistics for chemical potential calculations.
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

contains

    !---------------------------------------------------------------------------
    ! Subroutine: WidomTest
    !
    ! Purpose:
    !   Attempts a Widom trial insertion of a molecule of the given residue type.
    !   Computes energy contributions and updates Widom statistics. The move is 
    !   systematically rejected to avoid modifying the system state.
    !
    !---------------------------------------------------------------------------
    subroutine WidomTest(residue_type, molecule_index)

        implicit none

        ! Input arguments
        integer, intent(in) :: residue_type     ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID

        ! Local variables
        integer :: rand_mol_index               ! Randomly selected molecule index from the reservoir for copying geometry
        logical :: is_creation                  ! Flag indicating creation

        real(real64) :: deltaU, weight, T

        call CheckMoleculeIndex(molecule_index, NB_MAX_MOLECULE)

        ! Count trial move (success + fail)
        counter%trial_widom = counter%trial_widom + 1

        ! Compute old energy
        is_creation = .true.
        call ComputeOldEnergy(residue_type, molecule_index, old, is_creation = is_creation)

        ! Increase the residue and atom counts
        primary%num_residues(residue_type) = primary%num_residues(residue_type) + 1
        primary%num_atoms = primary%num_atoms + nb%atom_in_residue(residue_type)

        ! Save current Fourier terms (should be all zeros here)
        call SaveSingleMolFourierTerms(residue_type, molecule_index)

        ! Generate random insertion position within the simulation box
        call InsertAndOrientMolecule(residue_type, molecule_index, rand_mol_index)

        ! Compute new energy
        call ComputeNewEnergy(residue_type, molecule_index, new, is_creation = is_creation)

        ! Compute acceptance probability for the move
        write(*,*) "probability", mc_acceptance_probability(old, new, residue_type, TYPE_CREATION)

        ! Reject systematically
        call RejectCreationMove(residue_type, molecule_index)

        ! Bolzmann weight 
        deltaU = new%total - old%total ! In K units
        T = input%temp_K ! In K units

        if (new%total < old%total) then
            write(*,*) "WIDOM", deltaU, new%total, old%total, T
            stop 99
        end if

        weight = exp(-deltaU / T)

        widom_stat%weight(residue_type) = widom_stat%weight(residue_type) + weight
        widom_stat%sample(residue_type) = widom_stat%sample(residue_type) + 1

    end subroutine WidomTest

end module widom
