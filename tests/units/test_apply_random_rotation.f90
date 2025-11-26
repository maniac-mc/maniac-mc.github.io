program test_apply_random_rotation
    use iso_fortran_env, only: real64
    use simulation_state
    use monte_carlo_utils
    implicit none

    integer :: res_type, mol_index, n_atoms
    real(real64), allocatable :: original(:,:), rotated(:,:)
    logical :: pass2
    real(real64) :: d, local_error

    !--------------------------------------------------
    ! Tolerance for floating point comparisons
    !--------------------------------------------------
    local_error = 1.0d-12

    !--------------------------------------------------
    ! Allocate res%role
    !--------------------------------------------------
    if (allocated(res%role)) deallocate(res%role)
    allocate(res%role(1))
    res%role = TYPE_GUEST

    !--------------------------------------------------
    ! Initialize input parameters for rotation
    !--------------------------------------------------
    mc_input%rotation_step_angle = 0.1_real64  ! must be > 0 and <= TWOPI
    thermo%temperature = 300.0_real64       ! Kelvin
    thermo%fugacity = [0.01_real64]         ! allocate and assign if needed

    !--------------------------------------------------
    ! Minimal setup: 2 atoms residue
    !--------------------------------------------------
    res_type = 1
    mol_index = 1
    n_atoms = 2

    ! Allocate res%atom and set number of atoms
    if (allocated(res%atom)) deallocate(res%atom)
    allocate(res%atom(1))
    res%atom(1) = n_atoms

    ! Allocate guest%offset (3 x n_residues x n_molecules x n_atoms)
    if (allocated(guest%offset)) deallocate(guest%offset)
    allocate(guest%offset(3, 1, 1, n_atoms))
    guest%offset(:,:,:,:) = 0.0_real64

    ! Set positions for 2 atoms
    guest%offset(:, res_type, mol_index, 1) = [0.0_real64, 0.0_real64, 0.0_real64]
    guest%offset(:, res_type, mol_index, 2) = [1.0_real64, 0.0_real64, 0.0_real64]

    ! Copy original coordinates
    allocate(original(3,n_atoms))
    original = guest%offset(:, res_type, mol_index, 1:n_atoms)

    !--------------------------------------------------
    ! Apply rotation
    !--------------------------------------------------
    call apply_random_rotation(res_type, mol_index)

    ! Copy rotated coordinates
    allocate(rotated(3,n_atoms))
    rotated = guest%offset(:, res_type, mol_index, 1:n_atoms)

    !--------------------------------------------------
    ! Test 1: distances preserved
    !--------------------------------------------------
    d = sqrt(sum((rotated(:,1) - rotated(:,2))**2))
    if (abs(d - 1.0_real64) > local_error) then
        print *, "FAILED: Inter-atomic distance changed by rotation."
        stop 1
    end if

    !--------------------------------------------------
    ! Test 2: single-atom residue unchanged
    !--------------------------------------------------
    n_atoms = 1
    res%atom(1) = n_atoms
    if (allocated(guest%offset)) deallocate(guest%offset)
    allocate(guest%offset(3, 1, 1, n_atoms))
    guest%offset(:,:,:,:) = 42.0_real64

    call apply_random_rotation(res_type, mol_index)

    pass2 = all(abs(guest%offset(:, res_type, mol_index, 1) - 42.0_real64) < local_error)
    if (.not. pass2) then
        print *, "FAILED: Single-atom residue should remain unchanged."
        stop 1
    end if

    !--------------------------------------------------
    ! Clean up
    !--------------------------------------------------
    if (allocated(res%role)) deallocate(res%role)
    if (allocated(res%atom)) deallocate(res%atom)
    if (allocated(guest%offset)) deallocate(guest%offset)

end program test_apply_random_rotation
