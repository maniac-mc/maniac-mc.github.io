module molecule_deletion

    !===========================================================================
    ! Implements Monte Carlo deletion moves for molecules in a simulation
    ! box. Handles removal of molecules, updating of energies, Fourier terms,
    ! and optional placement into a reservoir.
    !===========================================================================

    use monte_carlo_utils
    use simulation_state
    use geometry_utils
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
    ! Attempts to remove a molecule of a given residue type from the simulation
    ! box. Computes energy changes, applies the Metropolis criterion, and
    ! either accepts or rejects the deletion.
    !---------------------------------------------------------------------------
    subroutine attempt_deletion_move(res_type, molecule_index)

        ! Input arguments
        integer, intent(in) :: res_type     ! Residue type to be moved
        integer, intent(in) :: molecule_index   ! Molecule ID

        ! Local variables
        real(real64) :: probability             ! Acceptance probability of deletion move
        integer :: last_molecule_index          ! Index of the last molecule in the primary box

        ! Return immediately if no molecules of type res_type exist
        if (primary%num%residues(res_type)==0) return

        ! Count trial move
        counter%deletions(1) = counter%deletions(1) + 1

        ! Energy of the previous configuration
        call compute_old_energy(res_type, molecule_index, is_deletion = .true.)
        call save_molecule_state(res_type, molecule_index, com_old = saved%com, offset_old = saved%offset)

        ! Record the index of the last molecule
        last_molecule_index = primary%num%residues(res_type)

        ! Delete molecule
        call remove_molecule(res_type, molecule_index, last_molecule_index)

        ! Update molecule and atom counts
        call update_counts(primary, res_type, -1)

        ! Calculate new energy
        call compute_new_energy(res_type, molecule_index, is_deletion = .true.)

        ! Compute acceptance probability for the move
        probability = compute_acceptance_probability(old, new, res_type, TYPE_DELETION)

        ! Accept or reject
        if (rand_uniform() <= probability) then
        
            ! Accept move
            call accept_deletion_move(res_type, last_molecule_index)
            
        else
        
            ! Reject move
            call reject_deletion_move(res_type, molecule_index, saved%com, saved%offset)

        end if

    end subroutine attempt_deletion_move

    !---------------------------------------------------------------------------
    ! Updates system energy and counters when a deletion move is accepted.
    ! Optionally adds the molecule to a reservoir if one exists.
    !---------------------------------------------------------------------------
    subroutine accept_deletion_move(res_type, last_molecule_index)

        ! Input parameters
        integer, intent(in) :: res_type         ! Residue type to be moved
        integer, intent(in) :: last_molecule_index  ! Index of the last molecule in the reservoir

        ! Local variable
        real(real64) :: trial_pos(3)                ! Random numbers for initial molecule position in the box

        ! Update total energies
        energy%recip_coulomb    = new%recip_coulomb
        energy%non_coulomb      = energy%non_coulomb    + new%non_coulomb   - old%non_coulomb
        energy%coulomb          = energy%coulomb        + new%coulomb       - old%coulomb
        energy%ewald_self       = energy%ewald_self     + new%ewald_self    - old%ewald_self
        energy%intra_coulomb    = energy%intra_coulomb  + new%intra_coulomb - old%intra_coulomb
        energy%total            = energy%total          + new%total         - old%total

        ! Count succesful move
        counter%deletions = counter%deletions + 1

        ! Add the molecule to the reservoir
        if (status%reservoir_provided) then

            ! Generate three random numbers in [0,1) and shift to [-0.5,0.5)
            call random_number(trial_pos)
            trial_pos = trial_pos - half

            ! Place the deleted molecule randomly in the reservoir
            gas%com(:, res_type, reservoir%num%residues(res_type)+1) = &
                trial_pos(1)*reservoir%cell%matrix(:, 1) + &
                trial_pos(2)*reservoir%cell%matrix(:, 2) + &
                trial_pos(3)*reservoir%cell%matrix(:, 3)
            gas%offset(:, res_type, reservoir%num%residues(res_type)+1, &
                1:res%atom(res_type)) = &
                gas%offset(:, res_type, last_molecule_index, 1:res%atom(res_type))

            call update_counts(reservoir, res_type, +1)

        end if

    end subroutine accept_deletion_move

    !---------------------------------------------------------------------------
    ! Restores molecule positions, orientations, counts, and Fourier states
    ! when a deletion move is rejected.
    !---------------------------------------------------------------------------
    subroutine reject_deletion_move(res_type, molecule_index, mol_com_old, site_offset_old)

        ! Input parameters
        integer, intent(in) :: res_type         ! Residue type to be moved
        integer, intent(in) :: molecule_index       ! Molecule ID
        real(real64), dimension(3) :: mol_com_old   ! For storing old molecule center-of-mass
        real(real64), dimension(:, :) :: site_offset_old ! For storing old molecule offset

        ! Restore previous residue/atom numbers
        call update_counts(primary, res_type, +1)


        ! Restore previous positions and orientation
        guest%com(:, res_type, molecule_index) = mol_com_old(:)
        guest%offset(:, res_type, molecule_index, 1:res%atom(res_type)) = &
            site_offset_old(:, 1:res%atom(res_type))

        ! Restore Fourier states (ik_alloc and dk_alloc)
        call restore_single_mol_fourier(res_type, molecule_index)

    end subroutine reject_deletion_move

end module molecule_deletion
