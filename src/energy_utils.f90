module energy_utils

    use simulation_state
    use tabulated_utils
    use geometry_utils
    use constants
    use ewald_kvectors
    use ewald_phase
    use ewald_energy
    use pairwise_energy_utils

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !------------------------------------------------------------------------------
    ! Compute the total energy of the system
    !------------------------------------------------------------------------------
    subroutine update_system_energy(box)

        ! Input arguments
        type(type_box), intent(inout) :: box

        ! Compute each energy components
        call evaluate_pairwise_energy(box)
        call compute_ewald_self()
        call compute_ewald_recip()
        call evaluate_intra_real_coulomb_energy(box)

        ! Calculate total energy
        energy%total = energy%recip_coulomb + energy%non_coulomb + energy%coulomb + &
            energy%ewald_self + energy%intra_coulomb ! In kcal/mol

        call energy_report(final=.false.) ! Print energy in log

    end subroutine update_system_energy

    !------------------------------------------------------------------------------
    ! Loops over all active residue types and their molecules, computes the
    ! intra-residue Coulomb energy for each molecule, and accumulates it into
    ! the total intra-residue energy in the simulation.
    !------------------------------------------------------------------------------
    subroutine evaluate_intra_real_coulomb_energy(box)

        ! Input arguments
        type(type_box), intent(inout) :: box    ! Box (reservoir or primary)

        ! Local variables
        integer :: res_type                     ! Type of residue
        integer :: mol_index                    ! Index of the molecule

        ! Initialize total intra-residue Coulomb energy
        energy%intra_coulomb = zero ! In kcal/mol

        ! Loop over all residue types
        do res_type = 1, res%number
            
            ! Skip inactive residues
            if (.not. thermo%is_active(res_type)) cycle

            ! Loop over all molecules of this residue type
            do mol_index = 1, box%num%residues(res_type)

                ! Compute intra-residue Coulomb energy for this molecule
                ! and accumulate it into total intra-residue energy
                energy%intra_coulomb = energy%intra_coulomb + &
                    intra_res_real_coulomb_energy(res_type, mol_index) ! In kcal/mol

            end do
        end do

    end subroutine evaluate_intra_real_coulomb_energy

    subroutine evaluate_pairwise_energy(box)

        ! Input arguments
        type(type_box), intent(inout) :: box    ! Box (reservoir or primary)

        ! Local variables
        integer :: res_type
        integer :: mol_index
        real(real64) :: e_non_coulomb
        real(real64) :: e_coulomb

        ! Initialize energies to zero
        energy%non_coulomb = zero
        energy%coulomb = zero

        ! Loop over all residue types
        do res_type = 1, res%number

            ! Loop over all molecule of type "res_type"
            do mol_index = 1, box%num%residues(res_type)

                ! Compute the energy for residue_1, molecule
                call pairwise_energy_for_molecule(box, res_type, &
                    mol_index, e_non_coulomb, e_coulomb)

                ! Add residue_1, molecule_1 energy to the total pairwise energy
                energy%non_coulomb = energy%non_coulomb + e_non_coulomb ! In kcal/mol
                energy%coulomb = energy%coulomb + e_coulomb ! In kcal/mol

            end do
        end do

    end subroutine evaluate_pairwise_energy

    !------------------------------------------------------------------------------
    ! Computes the reciprocal-space contribution to the Ewald electrostatic energy.
    ! This part handles the long-range component of Coulomb interactions by summing
    ! over reciprocal lattice vectors. The algorithm proceeds in three steps:
    !   1. Initialize weighting coefficients for each reciprocal lattice vector.
    !   2. Compute Fourier structure factors (∑ q_j exp(i·k·r_j)) for each molecule.
    !   3. Accumulate the reciprocal-space energy using precomputed factors.
    !------------------------------------------------------------------------------
    subroutine compute_ewald_recip()

        ! Step 1: Precompute weighting coefficients that depend only on |k|-vectors.
        ! These account for the Gaussian charge screening used in the Ewald method.
        call compute_reciprocal_weights()

        ! Step 2: Build Fourier terms e^(i·k·r) for every atom inthe system.
        ! This avoids recomputing expensive exponentials during the k-sum.
        call compute_all_ewald_phase_factors()

        ! Step 3: Compute reciprocal-space electrostatic energy using the structure
        ! factors and the precomputed reciprocal weighting coefficients.
        call compute_total_recip_energy(energy%recip_coulomb)

    end subroutine compute_ewald_recip

    !------------------------------------------------------------------------------
    ! Computes the Ewald self-interaction correction.
    ! In the Ewald summation, each point charge is artificially spread out by a
    ! Gaussian distribution. This leads to an unphysical interaction of each charge
    ! with its own Gaussian "image". The self-energy term removes this contribution
    ! to avoid overcounting.
    !------------------------------------------------------------------------------
    subroutine compute_ewald_self()

        ! Local variables
        integer :: res_type
        real(real64) :: e_ewald_self

        energy%ewald_self = zero ! Initialise ewald_self

        ! Loop over all residue types
        do res_type = 1, res%number

            ! Compute self-energy for a single molecule of this residue type
            e_ewald_self = zero
            call single_mol_ewald_self(res_type, e_ewald_self)

            ! Multiply by the number of molecules of this residue type
            e_ewald_self = e_ewald_self * primary%num%residues(res_type)

            ! Accumulate into total self-energy
            energy%ewald_self = energy%ewald_self + e_ewald_self    ! In kcal/mol

        end do

    end subroutine compute_ewald_self

    !------------------------------------------------------------------------------
    ! Computes the Ewald self-energy correction for a single molecule.
    !
    ! In the Ewald summation, each point charge interacts with an artificial
    ! Gaussian charge distribution representing itself. This leads to an unphysical
    ! self-interaction energy that must be subtracted to obtain the correct total
    ! electrostatic energy.
    !------------------------------------------------------------------------------
    subroutine single_mol_ewald_self(res_type, self_energy_1)

        ! Input arguments
        integer, intent(in) :: res_type           ! Residue type for the molecule
        real(real64), intent(out) :: self_energy_1    ! Computed self-energy for this molecule

        ! Local variables
        integer :: atom_index
        real(real64) :: charge_1

        ! Initialize self-energy accumulator
        self_energy_1 = zero

        ! Loop over all atoms in the residue
        do atom_index = 1, res%atom(res_type)

            charge_1 = primary%atoms%charges(res_type, atom_index)

            ! Skip atoms with negligible charge
            if (abs(charge_1) < error) cycle

            ! Add the self-energy contribution of this atom
            self_energy_1 = self_energy_1 - ewald%param%alpha / SQRTPI * charge_1**2  ! In units of e^2/Å

        end do

        ! Convert to kcal/mol at the end
        self_energy_1 = self_energy_1 * EPS0_INV_real  ! In kcal/mol

    end subroutine single_mol_ewald_self

end module energy_utils
