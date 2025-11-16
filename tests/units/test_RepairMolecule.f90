program test_RepairMolecule
    use, intrinsic :: iso_fortran_env, only: real64
    use readers_utils
    implicit none

    type(type_box) :: box
    integer :: nb_atoms
    real(real64), allocatable :: atom_x(:), atom_y(:), atom_z(:)
    real(real64) :: tol
    logical :: pass1, pass2, pass3
    real(real64), dimension(3) :: dr
    real(real64), dimension(3,3) :: recip
    real(real64), dimension(:,:), allocatable :: expected

    tol = 1.0d-12

    !======================================================
    ! Test 1: Cubic box
    !======================================================
    box%type = 1
    box%matrix = 0.0_real64
    box%matrix(1,1) = 1.0_real64
    box%matrix(2,2) = 1.0_real64
    box%matrix(3,3) = 1.0_real64

    nb_atoms = 2
    allocate(atom_x(nb_atoms), atom_y(nb_atoms), atom_z(nb_atoms))
    allocate(expected(3,nb_atoms))

    ! Atom 1 inside box, atom 2 just across boundary
    atom_x = [0.9_real64, -0.9_real64]
    atom_y = [0.0_real64,  0.0_real64]
    atom_z = [0.0_real64,  0.0_real64]

    call RepairMolecule(atom_x, atom_y, atom_z, nb_atoms, box)

    ! Expected: atom 2 should be at 1.1 (contiguous with atom 1)
    expected(:,1) = [0.9_real64, 0.0_real64, 0.0_real64]
    expected(:,2) = [1.1_real64, 0.0_real64, 0.0_real64]

    pass1 = all(abs([atom_x(1),atom_y(1),atom_z(1)]-expected(:,1))<tol) .and. &
           all(abs([atom_x(2),atom_y(2),atom_z(2)]-expected(:,2))<tol)

    if (.not. pass1) then
        print *, "Cubic RepairMolecule test FAILED:", atom_x, atom_y, atom_z
        stop 1
    end if

    !======================================================
    ! Test 2: Orthorhombic box
    !======================================================
    box%type = 2
    box%matrix = 0.0_real64
    box%matrix(1,1) = 2.0_real64
    box%matrix(2,2) = 3.0_real64
    box%matrix(3,3) = 4.0_real64

    atom_x = [1.9_real64, -1.9_real64]
    atom_y = [1.0_real64,  1.0_real64]
    atom_z = [2.0_real64,  2.0_real64]

    call RepairMolecule(atom_x, atom_y, atom_z, nb_atoms, box)

    expected(:,1) = [1.9_real64, 1.0_real64, 2.0_real64]
    expected(:,2) = [2.1_real64, 1.0_real64, 2.0_real64]

    pass2 = all(abs([atom_x(1),atom_y(1),atom_z(1)]-expected(:,1))<tol) .and. &
           all(abs([atom_x(2),atom_y(2),atom_z(2)]-expected(:,2))<tol)

    if (.not. pass2) then
        print *, "Orthorhombic RepairMolecule test FAILED:", atom_x, atom_y, atom_z
        stop 1
    end if

    !======================================================
    ! Test 3: Triclinic box
    !======================================================
    box%type = 3
    box%matrix(:,1) = [1.0_real64, 0.0_real64, 0.0_real64]
    box%matrix(:,2) = [0.5_real64, 1.0_real64, 0.0_real64]
    box%matrix(:,3) = [0.0_real64, 0.5_real64, 1.0_real64]

    ! Compute reciprocal
    recip = inverse(box%matrix)
    box%reciprocal = recip

    atom_x = [0.9_real64, -0.8_real64]
    atom_y = [0.1_real64,  0.2_real64]
    atom_z = [0.0_real64,  0.1_real64]

    call RepairMolecule(atom_x, atom_y, atom_z, nb_atoms, box)

    ! Expected: second atom should be shifted close to first
    ! (we just check continuity, not absolute values)
    dr = [atom_x(2)-atom_x(1), atom_y(2)-atom_y(1), atom_z(2)-atom_z(1)]
    pass3 = maxval(abs(dr)) < 1.0_real64  ! bond length is small

    if (.not. pass3) then
        print *, "Triclinic RepairMolecule test FAILED:", atom_x, atom_y, atom_z
        stop 1
    end if

    !-------------------------------------------------------------
    ! Final summary
    !-------------------------------------------------------------
    if (pass1 .and. pass2 .and. pass3) then
        print *, 'RepairMolecule test PASSED'
    end if

contains
    function inverse(mat) result(inv)
        real(real64), intent(in) :: mat(3,3)
        real(real64) :: inv(3,3)
        real(real64) :: det

        det = mat(1,1)*(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)) &
            - mat(1,2)*(mat(2,1)*mat(3,3)-mat(2,3)*mat(3,1)) &
            + mat(1,3)*(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))

        if (abs(det) < 1.0d-12) stop "Singular matrix in inverse"

        inv(1,1) =  (mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2))/det
        inv(1,2) = -(mat(1,2)*mat(3,3)-mat(1,3)*mat(3,2))/det
        inv(1,3) =  (mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2))/det

        inv(2,1) = -(mat(2,1)*mat(3,3)-mat(2,3)*mat(3,1))/det
        inv(2,2) =  (mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1))/det
        inv(2,3) = -(mat(1,1)*mat(2,3)-mat(1,3)*mat(2,1))/det

        inv(3,1) =  (mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))/det
        inv(3,2) = -(mat(1,1)*mat(3,2)-mat(1,2)*mat(3,1))/det
        inv(3,3) =  (mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1))/det
    end function inverse

end program test_RepairMolecule
