module energy_utils

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
    ! Compute the total energy of the system
    !------------------------------------------------------------------------------
    subroutine compute_system_energy(box)

        ! Input arguments
        type(type_box), intent(inout) :: box

        ! Compute each energy components
        call compute_pairwise_energy(box)
        call compute_ewald_self()
        call compute_ewald_recip()
        call compute_total_intra_residue_coulomb_energy()

        ! Calculate total energy
        energy%total = energy%recip_coulomb + energy%non_coulomb + energy%coulomb + &
            energy%ewald_self + energy%intra_coulomb ! In kcal/mol

        call energy_report(final=.false.) ! Print energy in log

    end subroutine compute_system_energy

    !------------------------------------------------------------------------------
    ! Loops over all active residue types and their molecules, computes the
    ! intra-residue Coulomb energy for each molecule, and accumulates it into
    ! the total intra-residue energy in the simulation.
    !------------------------------------------------------------------------------
    subroutine compute_total_intra_residue_coulomb_energy()

        ! Local variables
        integer :: res_type_1
        integer :: mol_index_1
        real(real64) :: e_intra_coulomb

        ! Initialize total intra-residue Coulomb energy
        energy%intra_coulomb = zero ! In kcal/mol

        ! Loop over all residue types
        do res_type_1 = 1, res%number
            
            ! Skip inactive residues
            if (thermo%is_active(res_type_1)) then

                ! Loop over all molecules of this residue type
                do mol_index_1 = 1, primary%num%residues(res_type_1)

                    ! Compute intra-residue Coulomb energy for this molecule
                    call compute_intra_real_energy_residue(res_type_1, mol_index_1, e_intra_coulomb)

                    ! Accumulate into total intra-residue energy
                    energy%intra_coulomb = energy%intra_coulomb + e_intra_coulomb ! In kcal/mol

                end do

            end if

        end do

    end subroutine compute_total_intra_residue_coulomb_energy

    subroutine compute_pairwise_energy(box)

        ! Input arguments
        type(type_box), intent(inout) :: box    ! Box (reservoir or primary)

        ! Local variables
        integer :: res_type_1
        integer :: mol_index_1
        real(real64) :: e_non_coulomb
        real(real64) :: e_coulomb

        ! Initialize energies to zero
        energy%non_coulomb = zero
        energy%coulomb = zero

        ! Loop over all residue types
        do res_type_1 = 1, res%number
            ! Loop over all molecule of type "res_type_1"
            do mol_index_1 = 1, box%num%residues(res_type_1)

                ! Compute the energy for residue_1, molecule_1
                call single_mol_pairwise_energy(box, res_type_1, &
                    mol_index_1, e_non_coulomb, e_coulomb)

                ! Add residue_1, molecule_1 energy to the total pairwise energy
                energy%non_coulomb = energy%non_coulomb + e_non_coulomb ! In kcal/mol
                energy%coulomb = energy%coulomb + e_coulomb ! In kcal/mol

            end do
        end do

    end subroutine compute_pairwise_energy

    !------------------------------------------------------------------------------
    ! Calculates the non-Coulombian and Coulomb (direct space)
    !------------------------------------------------------------------------------
    subroutine single_mol_pairwise_energy(box, res_type_1, mol_index_1, e_non_coulomb, e_coulomb)

        ! Input arguments
        type(type_box), intent(inout) :: box
        integer, intent(in) :: res_type_1        ! Residue type to be moved
        integer, intent(in) :: mol_index_1      ! Molecule ID
        real(real64), intent(out) :: e_non_coulomb
        real(real64), intent(out) :: e_coulomb

        ! Local variables
        integer :: atom_index_1, atom_index_2
        integer :: mol_index_2
        integer :: res_type_2
        real(real64) :: distance
        real(real64) :: sigma, epsilon              ! Epsilon and sigma LJ potential calculations
        real(real64) :: charge_1, charge_2          ! Charge for Coulomb interactions

        e_non_coulomb = zero
        e_coulomb = zero

        ! Loop over sites in molecule res_type
        do atom_index_1 = 1, res%atom(res_type_1)

            ! Loop over all molecule types 2
            do res_type_2 = 1, res%number

                ! Loop over all molecule index 2
                do mol_index_2 = 1, box%num%residues(res_type_2)

                    ! Remove intra molecular contribution
                    if ((mol_index_1 == mol_index_2) .and. &
                        (res_type_1 == res_type_2)) cycle

                    ! Enforce ordering to avoid double-counting
                    if ((res_type_2 < res_type_1) .or. &
                        ((res_type_2 == res_type_1) .and. &
                        (mol_index_2 <= mol_index_1))) cycle

                    ! Loop over all side of the selected molecule 2
                    do atom_index_2 = 1, res%atom(res_type_2)

                        ! Read pair parameters
                        sigma = coeff%sigma(res_type_1, res_type_2, atom_index_1, atom_index_2) ! In Angstrom
                        epsilon = coeff%epsilon(res_type_1, res_type_2, atom_index_1, atom_index_2) ! In kcal/mol
                        charge_1 = primary%atoms%charges(res_type_1, atom_index_1)                   ! In units of e
                        charge_2 = primary%atoms%charges(res_type_2, atom_index_2)                   ! In units of e

                        ! Calculate the distance, accouting for periodic boundary conditions
                        distance = minimum_image_distance(box, res_type_1, mol_index_1, atom_index_1, &
                                   res_type_2, mol_index_2, atom_index_2)                      ! In Angstrom

                        ! Update non-Coulomb energy
                        e_non_coulomb = e_non_coulomb + lennard_jones_energy(distance, sigma, epsilon)    ! In kcal/mol

                        ! Update Coulomb energy
                        e_coulomb = e_coulomb + coulomb_energy(distance, charge_1, charge_2)             ! In e^2/Å

                    end do
                end do
            end do
        end do

        ! Re-scale energy
        e_coulomb = e_coulomb * EPS0_INV_real   ! In kcal/mol

    end subroutine single_mol_pairwise_energy

    !------------------------------------------------------------------------------
    ! Function to compute Lennard-Jones interaction energy
    !------------------------------------------------------------------------------
    pure function lennard_jones_energy(r, sigma, epsilon) result(energy)

        ! Input variables
        real(real64), intent(in) :: r           ! Distance between the two atoms (in Å)
        real(real64), intent(in) :: sigma       ! Lennard-Jones sigma parameter (in Å)
        real(real64), intent(in) :: epsilon     ! Lennard-Jones epsilon parameter (in kJ/mol)

        ! Local variables
        real(real64) :: r6                      ! (sigma / r)^6 term of the LJ potential
        real(real64) :: r12                     ! (sigma / r)^12 term of the LJ potential
        real(real64) :: energy                  ! Lennard-Jones energy contribution (kcal/mol)

        if (r >= mc_input%real_space_cutoff) then

            ! Return 0 if distance larger than cutoff
            energy = zero                                           ! kcal/mol

        else

            if (use_table .and. r6_table%initialized .and. r12_table%initialized) then

                ! Use tabulated r^6 and r^12 if available and requested
                r6 = sigma**6 / lookup_tabulated(r6_table, r)        ! No units
                r12 = sigma**12 / lookup_tabulated(r12_table, r)     ! No units

            else

                ! Calculate r^6 and r^12
                r6 = (sigma / r)**6                                 ! No units
                r12 = r6 * r6                                       ! No units

            end if

            energy = four * epsilon * (r12 - r6)                    ! kcal/mol

        end if

    end function lennard_jones_energy

    !------------------------------------------------------------------------------
    ! Function to compute Coulomb interaction energy (Ewald direct-space term)
    !------------------------------------------------------------------------------
    pure function coulomb_energy(r, q1, q2) result(energy)

        ! Input variables
        real(real64), intent(in) :: r       ! Distance between the two atoms (in Å)
        real(real64), intent(in) :: q1      ! Atomic partial charge of atom 1 (in e)
        real(real64), intent(in) :: q2      ! Atomic partial charge of atom 2 (in e)

        ! Local variable
        real(real64) :: energy   ! Computed Coulomb energy contribution in units of e^2/Å

        ! Initialize energy
        energy = zero           ! In units of e^2/Å

        ! Skip negligible charges (and therefore returns an energy of 0)
        if ((abs(q1) < error) .or. (abs(q2) < error)) return

        ! Avoid division by zero (and therefore returns an energy of 0)
        if (r < error) return

        ! Compute Coulomb energy (tabulated or direct)
        if (use_table .and. erfc_r_table%initialized) then
            ! energy = q1*q2 * f(r)  , f(r) is erfc(r) / r from lookup table
            energy = q1 * q2 * lookup_tabulated(erfc_r_table, r)     ! In units of e^2/Å
        else
            ! Direct-space Coulomb potential with Ewald damping
            ! V(r) = (q1*q2) * erfc(alpha * r) / r
            energy = q1 * q2 * erfc(ewald%param%alpha * r) / r            ! In units of e^2/Å
        end if

    end function coulomb_energy

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
        integer :: res_type_1
        real(real64) :: e_ewald_self

        energy%ewald_self = zero ! Initialise ewald_self

        ! Loop over all residue types
        do res_type_1 = 1, res%number

            ! Compute self-energy for a single molecule of this residue type
            e_ewald_self = zero
            call single_mol_ewald_self(res_type_1, e_ewald_self)

            ! Multiply by the number of molecules of this residue type
            e_ewald_self = e_ewald_self * primary%num%residues(res_type_1)

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
        integer :: atom_index_1
        real(real64) :: charge_1

        ! Initialize self-energy accumulator
        self_energy_1 = zero

        ! Loop over all atoms in the residue
        do atom_index_1 = 1, res%atom(res_type)

            charge_1 = primary%atoms%charges(res_type, atom_index_1)

            ! Skip atoms with negligible charge
            if (abs(charge_1) < error) cycle

            ! Add the self-energy contribution of this atom
            self_energy_1 = self_energy_1 - ewald%param%alpha / SQRTPI * charge_1**2  ! In units of e^2/Å

        end do

        ! Convert to kcal/mol at the end
        self_energy_1 = self_energy_1 * EPS0_INV_real  ! In kcal/mol

    end subroutine single_mol_ewald_self

    !------------------------------------------------------------------------------
    ! Calculates the non-Coulombian and Coulomb (direct space)
    !------------------------------------------------------------------------------
    subroutine compute_pair_interaction_energy_singlemol(box, res_type_1, mol_index_1, e_non_coulomb, e_coulomb)

        ! Input arguments
        type(type_box), intent(inout) :: box
        integer, intent(in) :: res_type_1        ! Residue type 1
        integer, intent(in) :: mol_index_1      ! Molecule ID 1
        real(real64), intent(out) :: e_non_coulomb
        real(real64), intent(out) :: e_coulomb

        ! Local variables
        integer :: atom_index_1, atom_index_2       ! Atom index 1 et 2  
        integer :: mol_index_2                 ! Molecule ID 2
        integer :: res_type_2                   ! Residue type 1
        real(real64) :: distance                    ! Distance in Angstrom
        real(real64) :: r6, r12                     ! r^n for LJ potential calculations
        real(real64) :: sigma, epsilon              ! Epsilon and sigma LJ potential calculations
        real(real64) :: charge_1, charge_2          ! Charge for Coulomb interactions

        e_non_coulomb = zero
        e_coulomb = zero

        ! Loop over sites in molecule res_type
        do atom_index_1 = 1, res%atom(res_type_1)

            ! Loop over all molecule types 2
            do res_type_2 = 1, res%number

                ! Loop over all molecule index 2
                do mol_index_2 = 1, box%num%residues(res_type_2)

                    ! Remove intra molecular contribution
                    if ((mol_index_1 == mol_index_2) .and. &
                        (res_type_1 == res_type_2)) cycle

                    ! Loop over all side of the selected molecule 2
                    do atom_index_2 = 1, res%atom(res_type_2)

                        distance = minimum_image_distance(box, res_type_1, mol_index_1, atom_index_1, &
                                   res_type_2, mol_index_2, atom_index_2)

                        if (distance < mc_input%real_space_cutoff) then

                            ! LJ potential
                            sigma = coeff%sigma(res_type_1, res_type_2, atom_index_1, atom_index_2) ! In Å
                            epsilon = coeff%epsilon(res_type_1, res_type_2, atom_index_1, atom_index_2) ! In kcal/mol
                            r6 = (sigma / distance)**6                                  ! No units
                            r12 = r6 * r6                                               ! No units
                            e_non_coulomb = e_non_coulomb + four * epsilon * (r12 - r6) ! In kcal/mol

                        end if

                        ! Use Coulomb potential
                        charge_1 = primary%atoms%charges(res_type_1, atom_index_1)
                        charge_2 = primary%atoms%charges(res_type_2, atom_index_2)

                        ! Skip calculations if one charge is too small
                        if ((abs(charge_1) < error) .or. (abs(charge_2) < error)) cycle

                        e_coulomb = e_coulomb + charge_1 * charge_2 * (erfc(ewald%param%alpha * distance)) / distance ! e**2/Å

                    end do
                end do
            end do
        end do

        ! Convert to kcal/mol at the end
        e_coulomb = e_coulomb * EPS0_INV_real                          ! In kcal/mol

    end subroutine compute_pair_interaction_energy_singlemol

end module energy_utils
