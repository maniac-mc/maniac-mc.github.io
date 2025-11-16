program test_WrapIntoBox
    use, intrinsic :: iso_fortran_env, only: real64
    use geometry_utils
    implicit none

    type(type_box) :: box
    real(real64), dimension(3) :: pos, expected
    logical :: pass1, pass2, pass3
    real(real64) :: tol

    tol = 1.0d-12

    !======================================================
    ! Test 1: Cubic box
    !======================================================
    box%type = 1
    box%matrix = 0.0_real64
    box%matrix(1,1) = 2.0_real64
    box%matrix(2,2) = 2.0_real64
    box%matrix(3,3) = 2.0_real64

    pos = [3.0_real64, -2.5_real64, 4.1_real64]
    call WrapIntoBox(pos, box)

    expected = [-1.0_real64, -0.5_real64, 0.1_real64]
    pass1 = all(abs(pos - expected) < tol)
    if (.not. pass1) then
        print *, "Cubic box WrapIntoBox test FAILED: pos = ", pos, " expected = ", expected
        stop 1
    end if

    !======================================================
    ! Test 2: Orthorhombic box
    !======================================================
    box%type = 2
    box%matrix = 0.0_real64
    box%matrix(1,1) = 1.0_real64
    box%matrix(2,2) = 2.0_real64
    box%matrix(3,3) = 3.0_real64

    pos = [-1.2_real64, 2.5_real64, -3.7_real64]
    call WrapIntoBox(pos, box)

    expected = [-0.2_real64, 0.5_real64, -0.7_real64]  ! wrapped into [-L/2, L/2]
    pass2 = all(abs(pos - expected) < tol)
    if (.not. pass2) then
        print *, "Orthorhombic box WrapIntoBox test FAILED: pos = ", pos, " expected = ", expected
        stop 1
    end if

    !======================================================
    ! Test 3: Triclinic box
    !======================================================
    box%type = 3
    ! Example triclinic box
    box%matrix(:,1) = [2.0_real64, 0.1_real64, 0.05_real64]
    box%matrix(:,2) = [0.1_real64, 3.0_real64, 0.1_real64]
    box%matrix(:,3) = [0.05_real64, 0.1_real64, 4.0_real64]

    ! Compute reciprocal matrix
    call ComputeInverse(box)

    pos = [3.0_real64, -4.0_real64, 5.0_real64]
    call WrapIntoBox(pos, box)

    ! Convert to fractional coordinates
    pos = matmul(box%reciprocal, pos)
    expected = mod(pos + 0.5_real64, 1.0_real64) - 0.5_real64
    pass3 = all(abs(pos - expected) < tol)
    if (.not. pass3) then
        print *, "Triclinic box WrapIntoBox test FAILED: pos = ", pos, " expected = ", expected
        stop 1
    end if

    !-------------------------------------------------------------
    ! Final summary
    !-------------------------------------------------------------
    if (pass1 .and. pass2 .and. pass3) then
        print *, 'WrapIntoBox test PASSED'
    end if

end program test_WrapIntoBox
