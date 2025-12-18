module pairwise_energy_utils

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
    ! Calculates the non-Coulombian and Coulomb (direct space) energies
    ! for a given molecule
    !------------------------------------------------------------------------------
    subroutine pairwise_energy_for_molecule(box, res_i, mol_i, e_non_coulomb, e_coulomb)

        ! Input arguments
        type(type_box), intent(inout) :: box        ! Simulation box (primary or reservoir)
        integer, intent(in) :: res_i                ! Residue type index of the current molecule
        integer, intent(in) :: mol_i                ! Molecule index of the current residue

        ! Output arguments
        real(real64), intent(out) :: e_non_coulomb  ! Accumulated non-Coulomb (LJ) energy for this molecule
        real(real64), intent(out) :: e_coulomb      ! Accumulated Coulomb (electrostatic) energy for this molecule

        ! Local variables
        integer :: atom_i                            ! Index of current atom in the first molecule
        integer :: atom_j                            ! Index of current atom in the second molecule
        integer :: mol_j                             ! Molecule index in the second residue type
        integer :: res_j                             ! Residue type index for interacting molecule
        real(real64) :: r                            ! Interatomic distance (Å) with periodic boundary correction
        real(real64) :: sigma                        ! Lennard-Jones sigma parameter for atom pair (Å)
        real(real64) :: epsilon                      ! Lennard-Jones epsilon parameter for atom pair (kcal/mol)
        real(real64) :: q_i                          ! Charge of atom_i (in elementary charge units, e)
        real(real64) :: q_j                          ! Charge of atom_j (in elementary charge units, e)

        e_non_coulomb = zero
        e_coulomb = zero

        ! Loop over sites in molecule res_i
        do atom_i = 1, res%atom(res_i)

            ! Loop over all molecule types j
            do res_j = 1, res%number

                ! Loop over all molecule index j
                do mol_j = 1, box%num%residues(res_j)

                    ! Remove intra molecular contribution
                    if ((mol_i == mol_j) .and. (res_i == res_j)) cycle

                    ! Enforce ordering to avoid double-counting
                    if ((res_j < res_i) .or. ((res_j == res_i) .and. (mol_j <= mol_i))) cycle

                    ! Loop over all side of the selected molecule j
                    do atom_j = 1, res%atom(res_j)

                        ! Read pair parameters
                        sigma = coeff%sigma(res_i, res_j, atom_i, atom_j)           ! In Angstrom
                        epsilon = coeff%epsilon(res_i, res_j, atom_i, atom_j)       ! In kcal/mol
                        q_i = primary%atoms%charges(res_i, atom_i)                  ! In units of e
                        q_j = primary%atoms%charges(res_j, atom_j)                  ! In units of e

                        ! Calculate the distance, accouting for periodic boundary conditions
                        r = minimum_image_distance(box, res_i, mol_i, atom_i, res_j, mol_j, atom_j)                      ! In Angstrom

                        ! Update non-Coulomb energy
                        e_non_coulomb = e_non_coulomb + pairwise_lj_energy(r, sigma, epsilon)    ! In kcal/mol

                        ! Update Coulomb energy
                        e_coulomb = e_coulomb + pairwise_coulomb_energy(r, q_i, q_j)             ! In e^2/Å

                    end do
                end do
            end do
        end do

        ! Re-scale energy from e^2/Å to kcal/mol
        e_coulomb = e_coulomb * EPS0_INV_real   ! In kcal/mol

    end subroutine pairwise_energy_for_molecule

    !------------------------------------------------------------------------------
    ! Function to compute Lennard-Jones interaction energy
    !------------------------------------------------------------------------------
    pure function pairwise_lj_energy(r, sigma, epsilon) result(energy)

        ! Input variables
        real(real64), intent(in) :: r           ! Distance between the two atoms (in Å)
        real(real64), intent(in) :: sigma       ! Lennard-Jones sigma parameter (in Å)
        real(real64), intent(in) :: epsilon     ! Lennard-Jones epsilon parameter (in kJ/mol)

        ! Local variables
        real(real64) :: r6                      ! (sigma / r)^6 term of the LJ potential
        real(real64) :: r12                     ! (sigma / r)^12 term of the LJ potential

        ! Output value
        real(real64) :: energy                  ! Lennard-Jones energy contribution (kcal/mol)

        ! Initialize energy
        energy = zero   

        ! Return 0 if distance is larger than cutoff and avoid division by zero
        if (r >= mc_input%real_space_cutoff) return

        ! Return large value if overlap
        if (r < error) then
            energy = 1.0e20 ! Very large repulsive energy
            return
        end if

        ! Evaluate r^6 and r^12 either from table or directly
        if (use_table .and. r6_table%initialized .and. r12_table%initialized) then

            ! Use tabulated r^6 and r^12 if available and requested
            r6 = sigma**6 / lookup_tabulated(r6_table, r)           ! No units
            r12 = sigma**12 / lookup_tabulated(r12_table, r)        ! No units

        else

            ! Calculate r^6 and r^12
            r6 = (sigma / r)**6                                     ! No units
            r12 = r6 * r6                                           ! No units

        end if

        energy = four * epsilon * (r12 - r6)                        ! In kcal/mol

    end function pairwise_lj_energy

    !------------------------------------------------------------------------------
    ! Function to compute Coulomb interaction energy (Ewald direct-space term)
    !------------------------------------------------------------------------------
    pure function pairwise_coulomb_energy(r, q_i, q_j) result(energy)

        ! Input variables
        real(real64), intent(in) :: r           ! Distance between the two atoms (in Å)
        real(real64), intent(in) :: q_i         ! Atomic partial charge of atom 1 (in e)
        real(real64), intent(in) :: q_j         ! Atomic partial charge of atom 2 (in e)

        ! Output value
        real(real64) :: energy                  ! Computed Coulomb energy contribution in units of e^2/Å

        ! Initialize energy
        energy = zero

        ! Skip negligible charges and avoid division by zero
        if (abs(q_i) < error .or. abs(q_j) < error) return

        ! Return large value if overlap
        if (r < error) then
            energy = 1.0e20 ! Very large repulsive energy
            return
        end if

        ! Compute Coulomb energy (tabulated or direct)
        if (use_table .and. erfc_r_table%initialized) then

            ! energy = q_i*q_j * f(r)  , f(r) is erfc(r) / r from lookup table
            energy = q_i * q_j * lookup_tabulated(erfc_r_table, r)          ! In units of e^2/Å
        
        else
        
            ! Direct-space Coulomb potential with Ewald damping
            ! V(r) = (q_i*q_j) * erfc(alpha * r) / r
            energy = q_i * q_j * erfc(ewald%param%alpha * r) / r            ! In units of e^2/Å
        
        end if

    end function pairwise_coulomb_energy

end module pairwise_energy_utils
