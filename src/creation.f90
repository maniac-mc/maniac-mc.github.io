module molecule_creation

    !===========================================================================
    ! Implements Monte Carlo creation (insertion) moves for molecules in a
    ! simulation box. Handles inserting molecules, updating energies,
    ! Fourier terms, and optionally removing molecules from a reservoir.
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
    ! Attempts to insert a new molecule of a given residue type at a random
    ! position and orientation inside the simulation box. Computes energy
    ! changes, applies the Metropolis criterion, and either accepts or
    ! rejects the creation.
    !---------------------------------------------------------------------------
    subroutine attempt_creation_move(res_type, mol_index)

        ! Input arguments
        integer, intent(in) :: res_type     ! Residue type to be moved
        integer, intent(in) :: mol_index   ! Molecule ID

        ! Local variables
        integer :: rand_mol_index               ! Randomly selected molecule index from the reservoir for copying geometry
        real(real64) :: probability             ! Acceptance probability of creation move

        call check_molecule_index(mol_index)

        ! Count trial move
        counter%creations(1) = counter%creations(1) + 1

        ! Compute old energy
        call compute_old_energy(res_type, mol_index, is_creation = .true.)

        ! Increase the residue and atom counts
        call update_counts(primary, res_type, +1)

        ! Save current Fourier terms (should be all zeros here)
        call save_single_mol_fourier_terms(res_type, mol_index)

        ! Generate random insertion position within the simulation box
        call insert_and_orient_molecule(res_type, mol_index, rand_mol_index)

        ! Compute new energy
        call compute_new_energy(res_type, mol_index, is_creation = .true.)

        ! Compute acceptance probability for the move
        probability = compute_acceptance_probability(old, new, res_type, TYPE_CREATION)

        ! Accept or reject
        if (rand_uniform() <= probability) then
        
            ! Accept move
            call accept_creation_move(res_type, rand_mol_index)

        else
        
            ! Reject move
            call reject_creation_move(res_type, mol_index)

        end if

    end subroutine attempt_creation_move

    !---------------------------------------------------------------------------
    ! Updates system energies and counters when a creation move is accepted.
    ! Optionally removes the inserted molecule from the reservoir.
    !---------------------------------------------------------------------------
    subroutine accept_creation_move(res_type, rand_mol_index)

        ! Input parameters
        integer, intent(in) :: res_type     ! Residue type to be moved
        integer, intent(in) :: rand_mol_index               ! Randomly selected molecule index from the reservoir for copying geometry

        ! Local variable
        integer :: last_mol_index          ! Index of the last molecule in the reservoir (used when removing a molecule)
        
        ! Update total energies
        energy%recip_coulomb = new%recip_coulomb
        energy%non_coulomb = energy%non_coulomb + new%non_coulomb - old%non_coulomb
        energy%coulomb = energy%coulomb + new%coulomb - old%coulomb
        energy%ewald_self = energy%ewald_self + new%ewald_self - old%ewald_self
        energy%intra_coulomb = energy%intra_coulomb + new%intra_coulomb - old%intra_coulomb
        energy%total = energy%total + new%total - old%total

        ! Count successful move
        counter%creations = counter%creations + 1

        ! Remove molecule from reservoir if present
        if (status%reservoir_provided) then

            ! Replace mol_index with the last molecule in the list to maintain continuity
            last_mol_index = reservoir%num%residues(res_type)
            gas%com(:, res_type, rand_mol_index) = gas%com(:, res_type, last_mol_index)
            gas%offset(:, res_type, rand_mol_index, 1:res%atom(res_type)) = &
                gas%offset(:, res_type, last_mol_index, 1:res%atom(res_type))

            ! Increase the residue and atom counts
            call update_counts(reservoir, res_type, -1)

        end if

    end subroutine accept_creation_move

end module molecule_creation
