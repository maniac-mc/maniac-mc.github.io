module energy_utils

    use simulation_state
    use tabulated_utils
    use geometry_utils
    use constants
    use ewald_kvectors
    use ewald_phase
    use ewald_energy
    use pairwise_energy_utils
    use self_energy_utils

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
        call evaluate_ewald_self_energy()
        call evaluate_ewald_recip_energy()
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

                ! Compute the energy for the molecule
                call pairwise_energy_for_molecule(box, res_type, mol_index, e_non_coulomb, e_coulomb)

                ! Add the molecule energy to the total pairwise energy
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
    subroutine evaluate_ewald_recip_energy()

        ! Precompute weighting coefficients that depend only on |k|-vectors.
        ! These account for the Gaussian charge screening used in the Ewald method.
        call compute_reciprocal_weights()

        ! Build Fourier terms e^(i·k·r) for every atom inthe system.
        ! This avoids recomputing expensive exponentials during the k-sum.
        call compute_all_ewald_phase_factors()

        ! Compute reciprocal-space electrostatic energy using the structure
        ! factors and the precomputed reciprocal weighting coefficients.
        call compute_total_reciprocal_energy()

    end subroutine evaluate_ewald_recip_energy

end module energy_utils
