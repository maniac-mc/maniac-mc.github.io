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
    ! for a specific molecule or residue, enabling rollback after a rejected move.
    !--------------------------------------------------------------------
    subroutine save_single_mol_fourier_terms(residue_type, molecule_index, symmetrize_x)

        ! Input arguments
        integer, intent(in) :: residue_type       ! Residue type identifier
        integer, intent(in) :: molecule_index     ! Molecule index to save
        logical, intent(in), optional :: symmetrize_x   ! Whether to save negative kx

        ! Local variables
        integer :: atom_index                     ! Atom index within the residue
        integer :: dim                            ! Cartesian dimension (1=x,2=y,3=z)
        integer :: k_idx                          ! Reciprocal index
        integer :: idx                            ! Index over precomputed k-vectors
        logical :: do_sym                         ! Local symmetry flag

        ! Determine if saving negative kx is requested
        do_sym = .false.
        if (present(symmetrize_x)) do_sym = symmetrize_x

        !----------------------------------------------------------
        ! Save per-atom Fourier phase factors (IKX, IKY, IKZ)
        !----------------------------------------------------------
        do atom_index = 1, nb%atom_in_residue(residue_type)

            do dim = 1, 3

                do k_idx = 0, ewald%kmax(dim)

                    ! Positive k
                    ewald%phase_factor_old(dim, atom_index, k_idx) = &
                        ewald%phase_factor(dim, residue_type, molecule_index, atom_index, k_idx)

                    ! Negative k (only if non-zero)
                    if (k_idx /= 0) then
                        if (dim /= 1 .or. do_sym) then
                            ewald%phase_factor_old(dim, atom_index, -k_idx) = &
                                ewald%phase_factor(dim, residue_type, molecule_index, atom_index, -k_idx)
                        end if
                    end if

                end do  ! k_idx

            end do  ! dim

        end do  ! atom_index

        !------------------------------------------------------
        ! Save reciprocal amplitudes A(k)
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
        integer, intent(in) :: residue_type       ! Residue type index
        integer, intent(in) :: molecule_index     ! Molecule index to restore
        logical, intent(in), optional :: symmetrize_x ! Restore negative kx if true

        ! Local variables
        integer :: atom_index_1                   ! Atom index
        integer :: dim                            ! Cartesian dimension (1=x,2=y,3=z)
        integer :: k_idx                          ! Reciprocal index
        integer :: idx                            ! Index over k-vector list
        logical :: do_sym                         ! Whether to symmetrize X

        ! Determine if we symmetrize kx only
        do_sym = .false.
        if (present(symmetrize_x)) do_sym = symmetrize_x

        !----------------------------------------------------------
        ! Restore per-atom Fourier phase factors (IKX, IKY, IKZ)
        !----------------------------------------------------------
        do atom_index_1 = 1, nb%atom_in_residue(residue_type)

            do dim = 1, 3

                do k_idx = 0, ewald%kmax(dim)

                    ! Positive k
                    ewald%phase_factor(dim, residue_type, molecule_index, atom_index_1, k_idx) = &
                        ewald%phase_factor_old(dim, atom_index_1, k_idx)

                    ! Negative k
                    if (k_idx /= 0) then
                        if (dim /= 1 .or. do_sym) then
                            ewald%phase_factor(dim, residue_type, molecule_index, atom_index_1, -k_idx) = &
                                ewald%phase_factor_old(dim, atom_index_1, -k_idx)
                        end if
                    end if

                end do
            end do
        end do

        !------------------------------------------------------
        ! Restore reciprocal amplitudes A(k)
        !------------------------------------------------------
        do idx = 1, ewald%num_kvectors
            ewald%recip_amplitude(idx) = ewald%recip_amplitude_old(idx)
        end do

    end subroutine restore_single_mol_fourier

    !--------------------------------------------------------------------
    ! Replaces the Fourier-space phase factors of one molecule (`index_1`)
    ! with those from another molecule (`index_2`) of the same residue type.
    !--------------------------------------------------------------------
    subroutine replace_fourier_terms_single_mol(residue_type, index_1, index_2, symmetrize_x)

        ! Input arguments
        integer, intent(in) :: residue_type     ! Residue type index
        integer, intent(in) :: index_1          ! Destination molecule index
        integer, intent(in) :: index_2          ! Source molecule index
        logical, intent(in), optional :: symmetrize_x   ! Copy negative kx if true

        ! Local variables
        integer :: k_idx                        ! Wave index (k)
        integer :: atom_index_1                 ! Atom index within residue
        logical :: do_sym                       ! Whether to symmetrize X
        integer :: dim                           ! Spatial dimension (1..3)

        ! Determine if we should symmetrize x
        do_sym = .false.
        if (present(symmetrize_x)) do_sym = symmetrize_x

        ! Restore IKX, IKY, IKZ from backups
        do atom_index_1 = 1, nb%atom_in_residue(residue_type)

            do dim = 1, 3

                do k_idx = 0, ewald%kmax(dim)

                    ! Always copy positive (and zero) k
                    ewald%phase_factor(dim, residue_type, index_1, atom_index_1, k_idx) = &
                        ewald%phase_factor(dim, residue_type, index_2, atom_index_1, k_idx)

                    ! Copy negative k only when allowed (when d=1, only if symmetrize_x is true, or when d=2,3)
                    if (k_idx /= 0) then
                        if (dim /= 1 .or. do_sym) then
                            ewald%phase_factor(dim, residue_type, index_1, atom_index_1, -k_idx) = &
                                ewald%phase_factor(dim, residue_type, index_2, atom_index_1, -k_idx)
                        end if
                    end if

                end do
            end do
        end do

    end subroutine replace_fourier_terms_single_mol

    !--------------------------------------------------------------------
    ! Computes and stores the Fourier structure factors exp(i * k · r)
    ! for all atoms in all molecules across all residue types.
    !--------------------------------------------------------------------
    subroutine compute_all_fourier_terms()

        ! Local variables
        integer :: residue_type
        integer :: molecule_index

        ! Loop over all residue types
        do residue_type = 1, nb%type_residue

            ! Loop over all molecules of this residue type
            do molecule_index = 1, primary%num%residues(residue_type)

                ! Compute Fourier terms for a single molecule
                call single_mol_fourier_terms(residue_type, molecule_index)

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
        integer :: idim                         ! component index: 1=X, 2=Y, 3=Z
        type(type_coordinate), pointer :: coord ! Pointer for host or guest coordinate

        ! Return the correct pointer (host, guest, or gas)
        coord => get_coord(res_type)

        do atom_index_1 = 1, nb%atom_in_residue(res_type)

            ! Atom coordinate
            atom = coord%com(:, res_type, mol_index) + &
                coord%offset(:, res_type, mol_index, atom_index_1)

            ! Compute the phase vector components as the dot product of the atom position
            ! with each reciprocal lattice vector (columns of reciprocal_box), scaled by 2π.
            ! Compute the phase vector (2π * reciprocal_boxᵀ · atom)
            phase = compute_atom_phase(atom, primary%cell%reciprocal)

            ! Precompute the complex exponential (phase) factors for this atom
            ! along each Cartesian direction. These factors will be used repeatedly
            ! in the reciprocal-space sum for the Ewald energy.
            do idim = 1, 3
                ewald%temp_1d(:) = ewald%temp(idim, :)
                call compute_phase_factor(ewald%temp_1d(:), phase(idim), ewald%kmax(idim))
                ewald%phase_factor(idim, res_type, mol_index, atom_index_1, -ewald%kmax(idim):ewald%kmax(idim)) = ewald%temp_1d(:)
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
