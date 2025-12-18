module self_energy_utils

    use simulation_state
    use tabulated_utils
    use geometry_utils
    use constants
    use ewald_kvectors
    use ewald_phase
    use ewald_energy

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !------------------------------------------------------------------------------
    ! Computes the Ewald self-interaction correction.
    ! In the Ewald summation, each point charge is artificially spread out by a
    ! Gaussian distribution. This leads to an unphysical interaction of each charge
    ! with its own Gaussian "image". The self-energy term removes this contribution
    ! to avoid overcounting.
    !------------------------------------------------------------------------------
    subroutine evaluate_ewald_self_energy()

        ! Local variables
        integer :: res_i                    ! Index of residue type
        real(real64) :: self_energy_mol     ! Self-energy of a single molecule (kcal/mol)

        energy%ewald_self = zero            ! Initialise ewald_self

        ! Loop over all residue types
        do res_i = 1, res%number

            ! Compute self-energy for a single molecule of this residue type
            self_energy_mol = single_mol_ewald_self(res_i)

            ! Multiply by the number of molecules of this residue type
            ! Total self-energy contribution from this residue type
            self_energy_mol = self_energy_mol * primary%num%residues(res_i)

            ! Accumulate into total self-energy
            energy%ewald_self = energy%ewald_self + self_energy_mol ! In kcal/mol

        end do

    end subroutine evaluate_ewald_self_energy

    !------------------------------------------------------------------------------
    ! Computes the Ewald self-energy correction for a single molecule.
    ! Description : In the Ewald summation, each point charge interacts with an artificial
    ! Gaussian charge distribution representing itself. This leads to an unphysical
    ! self-interaction energy that must be subtracted to obtain the correct total
    ! electrostatic energy.
    !------------------------------------------------------------------------------
    pure function single_mol_ewald_self(res_i) result(self_energy)

        ! Input arguments
        integer, intent(in) :: res_i                ! Residue type for the molecule

        ! Output value
        real(real64) :: self_energy                 ! Computed self-energy (kcal/mol)

        ! Local variables
        integer :: atom_i                           ! Index of current atom
        real(real64) :: q_i                         ! Charge of atom_i (in elementary charge units, e)

        ! Initialize self-energy accumulator
        self_energy = zero

        ! Loop over all atoms in the residue
        do atom_i = 1, res%atom(res_i)

            q_i = primary%atoms%charges(res_i, atom_i)

            ! Skip atoms with negligible charge
            if (abs(q_i) < error) cycle

            ! Add the self-energy contribution of this atom
            self_energy = self_energy - ewald%param%alpha / SQRTPI * q_i**2  ! In units of e^2/Å

        end do

        ! Convert from e^2/Å to kcal/mol
        self_energy = self_energy * EPS0_INV_real  ! In kcal/mol

    end function single_mol_ewald_self

end module self_energy_utils
