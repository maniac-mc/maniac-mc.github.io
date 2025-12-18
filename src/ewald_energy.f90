module ewald_energy

    use simulation_state
    use output_utils
    use constants
    use geometry_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !--------------------------------------------------------------------
    ! Computes the total reciprocal-space Coulomb energy for the system.
    ! Loops over all precomputed k-vectors and evaluates:
    ! E(k) = form_factor * W(k) * |A(k)|^2
    ! where A(k) is the structure factor amplitude, W(k) is the reciprocal-space weight
    ! and the form_factor accounts for k vs -k symmetry
    !--------------------------------------------------------------------
    subroutine compute_total_reciprocal_energy()

        ! Internal variables
        real(real64) :: Ek              ! Ewald recip energy Ek in e^2 x Å^2
        integer :: idx                  ! Index over precomputed reciprocal vectors
        real(real64) :: form_factor     ! Factor to account for symmetry (k vs -k)
        real(real64) :: Wk              ! Precomputed Ewald reciprocal-space weight
        real(real64) :: Ak_square       ! Squared modulus of the structure factor amplitude

        ! Initialize
        Ek = zero

        ! Loop over all precomputed reciprocal lattice vectors
        do idx = 1, ewald%param%nkvec

            ! Use pre-computed form factor
            form_factor = ewald%form_factor(idx)

            ! Compute the structure factor (complex amplitude) for this k-vector
            ewald%Ak(idx) = compute_all_recip_amplitude(ewald%kvectors(idx)%kx, &
                ewald%kvectors(idx)%ky, ewald%kvectors(idx)%kz)

            ! Compute squared modulus of the structure factor amplitude
            Ak_square = amplitude_squared(ewald%Ak(idx)) ! e^2

            ! Retrieve the precomputed reciprocal-space weight
            Wk = ewald%kweights(idx)! Å^2

            ! Accumulate reciprocal-space energy:
            ! form_factor * reciprocal constant * |amplitude|^2
            ! E_k = form_factor * Wk * |Ak|^2
            Ek = Ek + form_factor * Wk * Ak_square ! In e^2 x Å^2

        end do

        ! Convert accumulated energy to correct units (kcal/mol)
        energy%recip_coulomb = Ek * EPS0_INV_real * TWOPI / primary%cell%volume ! In kcal/mol

    end subroutine compute_total_reciprocal_energy

    !--------------------------------------------------------------------
    ! Computes the reciprocal-space Coulomb energy contribution from a
    ! single molecule or residue using the Ewald summation method.
    !--------------------------------------------------------------------
    subroutine update_reciprocal_amplitude_single_mol(res_type, mol_index, is_creation, is_deletion)

        ! Input arguments
        integer, intent(in) :: res_type         ! Index of the residue type
        integer, intent(in) :: mol_index        ! Index of the molecule in the system
        logical, intent(in), optional :: is_creation
        logical, intent(in), optional :: is_deletion

        ! Local variables
        integer :: kx_idx, ky_idx, kz_idx      ! Components of current reciprocal lattice vector
        integer :: idx                         ! Loop index over precomputed k-vectors
        integer :: natoms                      ! Number of atoms in this residue type
        real(real64) :: form_factor            ! Symmetry factor: 1 if kx=0, 2 otherwise
        logical :: creation_flag
        logical :: deletion_flag

        ! Ensure that the routine deals with guest, not host
        if (res%role(res_type) == TYPE_HOST) call abort_run("Inconsistence in fourier routine")

        creation_flag = present_or_false(is_creation)
        deletion_flag = present_or_false(is_deletion)

        ! Atom charges in this residue
        natoms = res%atom(res_type)
        ewald%q_buffer(1:natoms) = primary%atoms%charges(res_type, 1:natoms)

        ! Loop over all precomputed reciprocal lattice vectors
        do idx = 1, ewald%param%nkvec

            ! Use pre-computed form factor
            form_factor = ewald%form_factor(idx)

            ! Current k-vector indices
            kx_idx = ewald%kvectors(idx)%kx
            ky_idx = ewald%kvectors(idx)%ky
            kz_idx = ewald%kvectors(idx)%kz

            ! Compute phase factors for the current residue
            ewald%phase%new(1:natoms) = ewald%phase%factor_guest(1, res_type, mol_index, 1:natoms, kx_idx) * &
                ewald%phase%factor_guest(2, res_type, mol_index, 1:natoms, ky_idx) * &
                ewald%phase%factor_guest(3, res_type, mol_index, 1:natoms, kz_idx)

            ewald%phase%old(1:natoms) = ewald%phase%factor_old(1, 1:natoms, kx_idx) * &
                ewald%phase%factor_old(2, 1:natoms, ky_idx) * &
                ewald%phase%factor_old(3, 1:natoms, kz_idx)

            ! Update Fourier coefficient A(k)
            if (creation_flag) then

                ! Molecule creation
                ! A(k) ← A(k) + Σ q_i [ e^(i k·r_i,new) ]
                ewald%Ak(idx) = ewald%Ak(idx) + sum(ewald%q_buffer(1:natoms) * ewald%phase%new(1:natoms))
            
            else if (deletion_flag) then
            
                ! Molecule deletion
                ! A(k) ← A(k) + Σ q_i [ - e^(i k·r_i,old) ]
                ewald%Ak(idx) = ewald%Ak(idx) - sum(ewald%q_buffer(1:natoms) * ewald%phase%old(1:natoms))
            
            else
            
                ! Standard move (translation, rotation)
                ! A(k) ← A(k) + Σ q_i [ e^(i k·r_i,new) - e^(i k·r_i,old) ]
                ewald%Ak(idx) = ewald%Ak(idx) + sum(ewald%q_buffer(1:natoms) * &
                    (ewald%phase%new(1:natoms) - ewald%phase%old(1:natoms)))
            
            end if

        end do

    end subroutine update_reciprocal_amplitude_single_mol

    !--------------------------------------------------------------------
    ! Compute reciprocal energy from current A(k)
    !--------------------------------------------------------------------
    pure function reciprocal_ewald_energy() result(energy)

        real(real64) :: energy
        integer :: idx
        real(real64) :: amplitude_sq

        energy = zero

        do idx = 1, ewald%param%nkvec

            ! Compute squared modulus of the structure factor amplitude
            amplitude_sq = amplitude_squared(ewald%Ak(idx))

            !----------------------------------------------
            ! Accumulate reciprocal-space energy:
            ! E_k = form_factor * W(k) * |A(k)|^2
            ! where W(k) is the precomputed reciprocal constant
            !----------------------------------------------
            energy = energy + ewald%form_factor(idx) * &
                    ewald%kweights(idx) * amplitude_sq
        end do

        ! Convert accumulated energy to physical units
        energy = energy * EPS0_INV_real * TWOPI / primary%cell%volume ! In kcal/mol

    end function reciprocal_ewald_energy


    !--------------------------------------------------------------------
    ! Computes the self-interaction correction for the reciprocal-space
    ! Coulomb energy of a single molecule/residue in the Ewald summation.
    ! In Ewald summation, each charge interacts with all periodic images
    ! of the system. This includes a spurious interaction of a charge
    ! with itself, which must be subtracted to obtain the correct energy.
    ! The self-energy for a single charge q_i is :
    !    E_self(i) = - (α / √π) * q_i^2
    ! where α is the Ewald screening parameter. The total self-energy
    ! for a molecule/residue is the sum over all atoms in that molecule.
    !--------------------------------------------------------------------
    pure function ewald_self_energy_single_mol(res_type) result(self_energy)

        ! Input arguments
        integer, intent(in) :: res_type        ! Residue type index

        ! Result
        real(real64) :: self_energy            ! Self-energy (kcal/mol)

        ! Local variables
        integer :: atom_index
        real(real64) :: charge

        ! Initialize self-energy accumulator
        self_energy = zero

        ! Loop over all atoms in the molecule/residue
        do atom_index = 1, res%atom(res_type)
            charge = primary%atoms%charges(res_type, atom_index)

            ! Skip atoms with negligible charge
            if (abs(charge) < error) cycle

            ! Add the individual atomic self-energy contribution
            self_energy = self_energy - ewald%param%alpha / SQRTPI * charge**2
        end do

        self_energy = self_energy * EPS0_INV_real ! In kcal/mol

    end function ewald_self_energy_single_mol

    !------------------------------------------------------------------------------
    ! Computes the real-space intramolecular Coulomb energy for a single
    ! molecule or residue using the Ewald summation method:
    !    E_ij_real = q_i * q_j * (erfc(α * r_ij) / r_ij)
    !------------------------------------------------------------------------------
    function intra_res_real_coulomb_energy(res_type, mol_index) result(u_intraCoulomb)
        
        ! Return value
        real(real64) :: u_intraCoulomb
        
        ! Input arguments
        integer, intent(in) :: res_type
        integer, intent(in) :: mol_index

        ! Local variables
        integer :: atom_index_1, atom_index_2
        real(real64) :: distance
        real(real64) :: charge_1, charge_2

        ! Initialize energy accumulator
        u_intraCoulomb = 0

        ! Loop over all unique atom pairs in the molecule
        do atom_index_1 = 1, res%atom(res_type)-1
            charge_1 = primary%atoms%charges(res_type, atom_index_1)

            do atom_index_2 = atom_index_1+1, res%atom(res_type)
                charge_2 = primary%atoms%charges(res_type, atom_index_2)

                ! Compute interatomic distance
                distance = minimum_image_distance(primary, res_type, mol_index, atom_index_1, &
                            res_type, mol_index, atom_index_2)

                ! Skip extremely small distances to avoid singularity
                if (distance < error) cycle

                ! Add the real-space intramolecular contribution (erfc(α * r) - 1) / r
                u_intraCoulomb = u_intraCoulomb + charge_1 * charge_2 * &
                    (erfc(ewald%param%alpha * distance) - one) / distance

            end do
        end do

        u_intraCoulomb = u_intraCoulomb * EPS0_INV_real ! In kcal/mol

    end function intra_res_real_coulomb_energy

    !--------------------------------------------------------------------
    ! Computes the reciprocal-space structure factor amplitude A(k) for
    ! a specified reciprocal lattice vector (kx, ky, kz):
    !     A(k) = Σ q_i * exp(i k · r_i)
    ! where q_i is the atomic charge and the exponential is represented
    ! by precomputed phase factors in x, y, and z directions.
    !--------------------------------------------------------------------
    pure function compute_all_recip_amplitude(kx_idx, ky_idx, kz_idx) result(Ak)

        ! Input arguments
        integer, intent(in) :: kx_idx, ky_idx, kz_idx       ! Reciprocal lattice vector indices

        ! Output
        complex(real64) :: Ak                               ! Accumulated amplitude A(k)

        ! Internal variables
        complex(real64) :: product                          ! Phase factor product for a single atom
        real(real64) :: charges                             ! Partial charge of the current atom
        integer :: res_type                                 ! Index of the current residue type
        integer :: mol_index                                ! Index of the current molecule
        integer :: atom_index                               ! Index of the current atom

        ! Initialize amplitude to zero (complex)
        Ak = (zero, zero)

        ! Loop over all residue types
        do res_type = 1, res%number

            ! For host
            if (res%role(res_type) == TYPE_HOST) then

                ! Loop over all molecules of this residue type
                do mol_index = 1, primary%num%residues(res_type)
                    ! Loop over sites in molecule
                    do atom_index = 1, res%atom(res_type)

                        ! Extract charge and phase product
                        charges = primary%atoms%charges(res_type, atom_index)
                        product = ewald%phase%factor_host(1, res_type, mol_index, atom_index, kx_idx) * &
                                ewald%phase%factor_host(2, res_type, mol_index, atom_index, ky_idx) * &
                                ewald%phase%factor_host(3, res_type, mol_index, atom_index, kz_idx)

                        ! Accumulate contribution from this atom:
                        ! charge * exp(i k · r) factor in x, y, z directions
                        Ak = Ak + charges * product

                    end do
                end do 

            ! For guest
            else if (res%role(res_type) == TYPE_GUEST) then

                ! Loop over all molecules of this residue type
                do mol_index = 1, primary%num%residues(res_type)
                    ! Loop over sites in molecule
                    do atom_index = 1, res%atom(res_type)

                        ! Extract charge and phase product
                        charges = primary%atoms%charges(res_type, atom_index)
                        product = ewald%phase%factor_guest(1, res_type, mol_index, atom_index, kx_idx) * &
                                ewald%phase%factor_guest(2, res_type, mol_index, atom_index, ky_idx) * &
                                ewald%phase%factor_guest(3, res_type, mol_index, atom_index, kz_idx)

                        ! Accumulate contribution from this atom:
                        ! charge * exp(i k · r) factor in x, y, z directions
                        Ak = Ak + charges * product

                    end do
                end do 

            end if

        end do

    end function compute_all_recip_amplitude

end module ewald_energy
