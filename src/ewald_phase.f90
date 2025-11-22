module ewald_phase

    use simulation_state
    use output_utils
    use constants
    use geometry_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none
    
contains

    !--------------------------------------------------------------------
    ! Saves the current Fourier-space phase factors and reciprocal amplitudes
    ! for a specific molecule or residue. This enables rollback after a
    ! rejected Monte Carlo or molecular dynamics move.
    !--------------------------------------------------------------------
    subroutine save_single_mol_fourier_terms(residue_type, molecule_index, symmetrize_x)

        ! Input arguments
        integer, intent(in) :: residue_type     ! Residue type identifier
        integer, intent(in) :: molecule_index   ! Index of the molecule to save
        logical, intent(in), optional :: symmetrize_x  ! If true, save negative kx terms

        ! Local variables
        integer :: atom_index                   ! Atom index within the residue
        integer :: kx_idx, ky_idx, kz_idx       ! Reciprocal vector indices
        integer :: idx                          ! Index over precomputed k-vectors
        logical :: do_sym                       ! Local flag

        ! Determine if we should symmetrize X
        do_sym = .false.
        if (present(symmetrize_x)) do_sym = symmetrize_x

        do atom_index = 1, nb%atom_in_residue(residue_type)

            ! Save IKX terms (kx_idx from 0 to kmax(1))
            do kx_idx = 0, ewald%kmax(1)
                ewald%phase_factor_old(1, atom_index, kx_idx) = &
                    ewald%phase_factor(1, residue_type, molecule_index, atom_index, kx_idx)
                    ! Optionally save negative kx for symmetry
                    if (do_sym .and. kx_idx /= 0) then
                        ewald%phase_factor_old(1, atom_index, -kx_idx) = &
                            ewald%phase_factor(1, residue_type, molecule_index, atom_index, -kx_idx)
                    end if
            end do

            ! Save IKY terms (KY from 0 to kmax(2)), include negative KY only if non-zero
            do ky_idx = 0, ewald%kmax(2)
                ewald%phase_factor_old(2, atom_index, ky_idx) = &
                    ewald%phase_factor(2, residue_type, molecule_index, atom_index, ky_idx)
                if (ky_idx /= 0) then
                    ewald%phase_factor_old(2, atom_index, -ky_idx) = &
                        ewald%phase_factor(2, residue_type, molecule_index, atom_index, -ky_idx)
                end if
            end do

            ! Save IKZ terms (KZ from 0 to kmax(3)), include negative KZ only if non-zero
            do kz_idx = 0, ewald%kmax(3)
                ewald%phase_factor_old(3, atom_index, kz_idx) = &
                    ewald%phase_factor(3, residue_type, molecule_index, atom_index, kz_idx)
                if (kz_idx /= 0) then
                    ewald%phase_factor_old(3, atom_index, -kz_idx) = &
                        ewald%phase_factor(3, residue_type, molecule_index, atom_index, -kz_idx)
                end if
            end do
        end do

        !------------------------------------------------------
        ! Save reciprocal amplitudes A(k)
        ! Loop directly over precomputed valid k-vectors
        !------------------------------------------------------
        do idx = 1, ewald%num_kvectors
            ewald%recip_amplitude_old(idx) = ewald%recip_amplitude(idx)
        end do

    end subroutine save_single_mol_fourier_terms

    !--------------------------------------------------------------------
    ! Restores previously saved Fourier-space phase factors and reciprocal
    ! amplitudes for a molecule or residue after a rejected move.
    !--------------------------------------------------------------------
    subroutine restore_single_mol_fourier(residue_type, molecule_index, symmetrize_x)

        ! Input arguments
        integer, intent(in) :: residue_type    ! Residue type identifier
        integer, intent(in) :: molecule_index  ! Index of the molecule to restore
        logical, intent(in), optional :: symmetrize_x  ! If true, restore negative kx terms

        ! Local variables
        integer :: atom_index_1                ! Atom index within the residue
        integer :: kx_idx, ky_idx, kz_idx      ! Reciprocal vector indices
        integer :: idx                         ! Index over precomputed k-vectors
        logical :: do_sym                       ! Local flag

        ! Determine if we should symmetrize X
        do_sym = .false.
        if (present(symmetrize_x)) do_sym = symmetrize_x

        ! Restore per-atom phase factors (IKX, IKY, IKZ)
        do atom_index_1 = 1, nb%atom_in_residue(residue_type)

            ! Restore IKX
            do kx_idx = 0, ewald%kmax(1)
                ewald%phase_factor(1, residue_type, molecule_index, atom_index_1, kx_idx) = &
                    ewald%phase_factor_old(1, atom_index_1, kx_idx)
                    ! Optionally restore negative kx for symmetry
                    if (do_sym .and. kx_idx /= 0) then
                        ewald%phase_factor(1, residue_type, molecule_index, atom_index_1, -kx_idx) = &
                            ewald%phase_factor_old(1, atom_index_1, -kx_idx)
                    end if
            end do

            ! Restore IKY (include negative KY only if non-zero)
            do ky_idx = 0, ewald%kmax(2)
                ewald%phase_factor(2, residue_type, molecule_index, atom_index_1, ky_idx) = &
                    ewald%phase_factor_old(2, atom_index_1, ky_idx)
                if (ky_idx /= 0) then
                    ewald%phase_factor(2, residue_type, molecule_index, atom_index_1, -ky_idx) = &
                        ewald%phase_factor_old(2, atom_index_1, -ky_idx)
                end if
            end do

            ! Restore IKZ (include negative KZ only if non-zero)
            do kz_idx = 0, ewald%kmax(3)
                ewald%phase_factor(3, residue_type, molecule_index, atom_index_1, kz_idx) = &
                    ewald%phase_factor_old(3, atom_index_1, kz_idx)
                if (kz_idx /= 0) then
                    ewald%phase_factor(3, residue_type, molecule_index, atom_index_1, -kz_idx) = &
                        ewald%phase_factor_old(3, atom_index_1, -kz_idx)
                end if
            end do
        end do

        !------------------------------------------------------
        ! Restore reciprocal amplitudes A(k) from backup
        ! Loop directly over precomputed valid k-vectors
        !------------------------------------------------------
        do idx = 1, ewald%num_kvectors
            ewald%recip_amplitude(idx) = ewald%recip_amplitude_old(idx)
        end do

    end subroutine restore_single_mol_fourier

    !--------------------------------------------------------------------
    ! Replaces the Fourier-space phase factors of one molecule (`index_1`)
    ! with those from another molecule (`index_2`) of the same residue type.
    ! Used when molecule identities or positions are swapped (e.g., in
    ! exchange Monte Carlo or symmetry operations).
    !--------------------------------------------------------------------
    subroutine replace_fourier_terms_single_mol(residue_type, index_1, index_2, symmetrize_x)

        ! Input arguments
        integer, intent(in) :: residue_type
        integer, intent(in) :: index_1
        integer, intent(in) :: index_2
        logical, intent(in), optional :: symmetrize_x ! If true, replace negative kx terms

        ! Local variables
        integer :: kx_idx, ky_idx, kz_idx
        integer :: atom_index_1
        logical :: do_sym                       ! Local flag

        ! Determine if we should symmetrize X
        do_sym = .false.
        if (present(symmetrize_x)) do_sym = symmetrize_x

        ! Restore IKX, IKY, IKZ from backups
        do atom_index_1 = 1, nb%atom_in_residue(residue_type)

            ! Restore IKX
            do kx_idx = 0, ewald%kmax(1)
                ewald%phase_factor(1, residue_type, index_1, atom_index_1, kx_idx) = &
                    ewald%phase_factor(1, residue_type, index_2, atom_index_1, kx_idx)
                   ! Optionally replace negative kx for symmetry
                    if (do_sym .and. kx_idx /= 0) then
                        ewald%phase_factor(1, residue_type, index_1, atom_index_1, -kx_idx) = &
                            ewald%phase_factor(1, residue_type, index_2, atom_index_1, -kx_idx)
                    end if
            end do

            ! Restore IKY (include negative KY only if non-zero)
            do ky_idx = 0, ewald%kmax(2)
                ewald%phase_factor(2, residue_type, index_1, atom_index_1, ky_idx) = &
                    ewald%phase_factor(2, residue_type, index_2, atom_index_1, ky_idx)
                if (ky_idx /= 0) then
                    ewald%phase_factor(2, residue_type, index_1, atom_index_1, -ky_idx) = &
                        ewald%phase_factor(2, residue_type, index_2, atom_index_1, -ky_idx)
                end if
            end do

            ! Restore IKZ (include negative KZ only if non-zero)
            do kz_idx = 0, ewald%kmax(3)
                ewald%phase_factor(3, residue_type, index_1, atom_index_1, kz_idx) = &
                    ewald%phase_factor(3, residue_type, index_2, atom_index_1, kz_idx)
                if (kz_idx /= 0) then
                    ewald%phase_factor(3, residue_type, index_1, atom_index_1, -kz_idx) = &
                        ewald%phase_factor(3, residue_type, index_2, atom_index_1, -kz_idx)
                end if
            end do
        end do

    end subroutine replace_fourier_terms_single_mol

    !--------------------------------------------------------------------
    ! Computes and stores the Fourier structure factors exp(i * k · r)
    ! for all atoms in all molecules across all residue types.
    !--------------------------------------------------------------------
    subroutine compute_all_fourier_terms()

        ! Local variables
        integer :: residue_type_1
        integer :: molecule_index_1

        ! Loop over all residue types
        do residue_type_1 = 1, nb%type_residue

            ! Loop over all molecules of this residue type
            do molecule_index_1 = 1, primary%num_residues(residue_type_1)

                ! Compute Fourier terms for a single molecule
                call single_mol_fourier_terms(residue_type_1, molecule_index_1)

            end do
        end do

    end subroutine compute_all_fourier_terms

    !--------------------------------------------------------------------
    ! Computes the per-atom Fourier-space phase factors exp(i * k · r)
    ! for all k-indices along each Cartesian direction for a given molecule.
    !--------------------------------------------------------------------
    subroutine single_mol_fourier_terms(res_type, mol_index)

        ! Input arguments
        integer, intent(in) :: res_type
        integer, intent(in) :: mol_index
        
        ! Local variables
        integer :: atom_index_1                 ! Atom index
        real(real64), dimension(3) :: atom      ! Atom coordinates in real space
        real(real64), dimension(3) :: phase     ! Phase factors for Fourier terms
        integer :: idim  ! component index: 1=X, 2=Y, 3=Z

        do atom_index_1 = 1, nb%atom_in_residue(res_type)

            atom = primary%mol_com(:, res_type, mol_index) + &
                primary%site_offset(:, res_type, mol_index, atom_index_1)

            ! Compute the phase vector components as the dot product of the atom position
            ! with each reciprocal lattice vector (columns of reciprocal_box), scaled by 2π.
            ! Compute the phase vector (2π * reciprocal_boxᵀ · atom)
            phase = compute_atom_phase(atom, primary%reciprocal)

            ! Precompute the complex exponential (phase) factors for this atom
            ! along each Cartesian direction. These factors will be used repeatedly
            ! in the reciprocal-space sum for the Ewald energy.
            do idim = 1, 3
                call compute_phase_factor(ewald%temp(idim,:), phase(idim), ewald%kmax(idim))
                ewald%phase_factor(idim, res_type, mol_index, atom_index_1, :) = ewald%temp(idim,:)
            end do

        end do

    end subroutine single_mol_fourier_terms

    !--------------------------------------------------------------------
    ! Computes the phase vector (k · r) for a single atom in reciprocal space.
    !--------------------------------------------------------------------
    pure function compute_atom_phase(atom_pos, reciprocal_box) result(phase)

        ! Input arguments
        real(real64), intent(in)  :: atom_pos(3)         ! Cartesian position of the atom in real space (x, y, z)
        real(real64), intent(in)  :: reciprocal_box(3,3) ! Reciprocal lattice vectors as columns of a 3×3 matrix

        ! Local variables
        integer :: i, j                                 ! Loop indices over spatial dimensions

        ! Function result
        real(real64) :: phase(3)                        ! Returned phase components (2π * k·r)

        ! For each reciprocal direction i ∈ {1,2,3}:
        ! phase(i) = 2π * Σ_j [ reciprocal_box(j,i) * atom_pos(j) ]
        do i = 1, 3
            phase(i) = zero
            do j = 1, 3
                ! Accumulate dot product of atom position (r_j)
                ! with i-th reciprocal lattice vector (k_i = column i of reciprocal_box)
                phase(i) = phase(i) + reciprocal_box(j,i) * atom_pos(j)
            end do
            ! Multiply by 2π to convert to the proper phase angle (k·r in radians)
            phase(i) = TWOPI * phase(i)
        end do

    end function compute_atom_phase

    !--------------------------------------------------------------------
    ! Computes 1D complex exponential phase factors exp(i * k * θ)
    ! along a single reciprocal axis for one atom.
    !--------------------------------------------------------------------
    pure subroutine compute_phase_factor(phase_factor_axis, phase_component, kmax)

        ! Input arguments
        complex(real64), intent(inout) :: phase_factor_axis(-kmax:kmax) ! Array of complex phase factors for all k indices along one axis
        real(real64), intent(in) :: phase_component                  ! Phase angle for this atom along the current axis
        integer, intent(in) :: kmax                                  ! Maximum k-index in the positive direction

        ! Local arguments
        integer :: k                                                 ! Loop index over k-values
        complex(real64) :: phase_exp                                 ! Temporary complex exponential for current k

        !   For k ∈ [-kmax, kmax]:
        !       phase_factor_axis(k) = exp(i * k * phase_component)
        do k = 0, kmax

            ! Compute complex exponential for positive k
            phase_exp = cmplx(dcos(k*phase_component), dsin(k*phase_component), kind=real64)
            phase_factor_axis(k) = phase_exp

            ! Exploit conjugate symmetry for negative k values
            if (k /= 0) then
                phase_factor_axis(-k) = conjg(phase_exp)
            end if

        end do

    end subroutine compute_phase_factor

end module ewald_phase
