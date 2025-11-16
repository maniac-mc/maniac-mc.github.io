program test_ApplyRandomRotation
    use iso_fortran_env, only: real64
    use simulation_state
    use monte_carlo_utils
    implicit none

    integer :: res_type, mol_index, n_atoms
    real(real64), allocatable :: original(:,:), rotated(:,:)
    logical :: pass1, pass2
    real(real64) :: d, tol

    ! Initialize input parameters for rotation
    input%rotation_step_angle = 0.1_real64  ! must be > 0 and <= TWOPI
    input%temp_K = 300.0_real64            ! Kelvin
    input%fugacity = [0.01_real64]         ! allocate and assign if needed

    !--------------------------------------------------
    ! Tolerance for floating-point comparison
    !--------------------------------------------------
    tol = 1.0e-12_real64

    !--------------------------------------------------
    ! Minimal setup: 2 atoms residue
    !--------------------------------------------------
    res_type = 1
    mol_index = 1
    n_atoms = 2

    ! Allocate nb%atom_in_residue and set number of atoms
    if (allocated(nb%atom_in_residue)) deallocate(nb%atom_in_residue)
    allocate(nb%atom_in_residue(1))
    nb%atom_in_residue(1) = n_atoms

    ! Allocate primary%site_offset (3 x n_residues x n_molecules x n_atoms)
    if (allocated(primary%site_offset)) deallocate(primary%site_offset)
    allocate(primary%site_offset(3, 1, 1, n_atoms))
    primary%site_offset(:,:,:,:) = 0.0_real64

    ! Set positions for 2 atoms
    primary%site_offset(:, res_type, mol_index, 1) = [0.0_real64, 0.0_real64, 0.0_real64]
    primary%site_offset(:, res_type, mol_index, 2) = [1.0_real64, 0.0_real64, 0.0_real64]

    ! Copy original coordinates
    allocate(original(3,n_atoms))
    original = primary%site_offset(:, res_type, mol_index, 1:n_atoms)

    !--------------------------------------------------
    ! Apply rotation
    !--------------------------------------------------
    call ApplyRandomRotation(res_type, mol_index)

    ! Copy rotated coordinates
    allocate(rotated(3,n_atoms))
    rotated = primary%site_offset(:, res_type, mol_index, 1:n_atoms)

    !--------------------------------------------------
    ! Test 1: distances preserved
    !--------------------------------------------------
    d = sqrt(sum((rotated(:,1) - rotated(:,2))**2))
    if (abs(d - 1.0_real64) > tol) then
        print *, "FAILED: Inter-atomic distance changed by rotation."
        stop 1
    end if

    !--------------------------------------------------
    ! Test 2: coordinates changed (rotation actually occurred)
    !--------------------------------------------------
    pass1 = any(abs(rotated - original) > tol)
    if (.not. pass1) then
        print *, "FAILED: Atom positions did not change after rotation."
        stop 1
    end if

    !--------------------------------------------------
    ! Test 3: single-atom residue unchanged
    !--------------------------------------------------
    n_atoms = 1
    nb%atom_in_residue(1) = n_atoms
    if (allocated(primary%site_offset)) deallocate(primary%site_offset)
    allocate(primary%site_offset(3, 1, 1, n_atoms))
    primary%site_offset(:,:,:,:) = 42.0_real64

    call ApplyRandomRotation(res_type, mol_index)

    pass2 = all(abs(primary%site_offset(:, res_type, mol_index, 1) - 42.0_real64) < tol)
    if (.not. pass2) then
        print *, "FAILED: Single-atom residue should remain unchanged."
        stop 1
    end if

end program test_ApplyRandomRotation
