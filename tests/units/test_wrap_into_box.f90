program test_wrap_into_box
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
    box%cell%shape = 1
    box%cell%matrix = 0.0_real64
    box%cell%matrix(1,1) = 2.0_real64
    box%cell%matrix(2,2) = 2.0_real64
    box%cell%matrix(3,3) = 2.0_real64

    pos = [3.0_real64, -2.5_real64, 4.1_real64]
    call wrap_into_box(pos, box)

    expected = [-1.0_real64, -0.5_real64, 0.1_real64]
    pass1 = all(abs(pos - expected) < tol)
    if (.not. pass1) then
        print *, "Cubic box wrap_into_box test FAILED: pos = ", pos, " expected = ", expected
        stop 1
    end if

    !======================================================
    ! Test 2: Orthorhombic box
    !======================================================
    box%cell%shape = 1
    box%cell%matrix = 0.0_real64
    box%cell%matrix(1,1) = 1.0_real64
    box%cell%matrix(2,2) = 2.0_real64
    box%cell%matrix(3,3) = 3.0_real64

    pos = [-1.2_real64, 2.5_real64, -3.7_real64]
    call wrap_into_box(pos, box)

    expected = [-0.2_real64, 0.5_real64, -0.7_real64]
    pass2 = all(abs(pos - expected) < tol)
    if (.not. pass2) then
        print *, "Orthorhombic box wrap_into_box test FAILED: pos = ", pos, " expected = ", expected
        stop 1
    end if

    !======================================================
    ! Test 3: Triclinic box
    !======================================================
    box%cell%shape = 2
    box%cell%matrix(:,1) = [2.0_real64, 0.1_real64, 0.05_real64]
    box%cell%matrix(:,2) = [0.1_real64, 3.0_real64, 0.1_real64]
    box%cell%matrix(:,3) = [0.05_real64, 0.1_real64, 4.0_real64]

    call compute_box_determinant_and_inverse(box)

    pos = [3.0_real64, -4.0_real64, 5.0_real64]
    call wrap_into_box(pos, box)

    ! Convert to fractional coordinates
    pos = matmul(box%cell%reciprocal, pos)
    expected = mod(pos + 0.5_real64, 1.0_real64) - 0.5_real64

    pass3 = all(abs(pos - expected) < tol)
    if (.not. pass3) then
        print *, "Triclinic box wrap_into_box test FAILED: pos = ", pos, " expected = ", expected
        stop 1
    end if

end program test_wrap_into_box
