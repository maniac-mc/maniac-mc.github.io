program test_ComputeCellProperties
    use, intrinsic :: iso_fortran_env, only: real64
    use geometry_utils
    implicit none

    type(type_box) :: box
    logical :: pass1, pass2, pass3
    real(real64) :: expected, tol

    tol = 1.0d-10   ! Numerical tolerance

    !======================================================
    ! Test 1: Cubic box (1x1x1 cube)
    !======================================================
    box%matrix = 0.0_real64
    box%matrix(1,1) = 1.0_real64
    box%matrix(2,2) = 1.0_real64
    box%matrix(3,3) = 1.0_real64

    call ComputeCellProperties(box)

    pass1 = abs(box%volume - 1.0_real64) < tol
    if (.not. pass1) then
        print *, "Cubic box test FAILED (volume = ", box%volume, ")"
        stop 1
    end if

    !======================================================
    ! Test 2: Orthorhombic box (2x3x4)
    !======================================================
    box%matrix = 0.0_real64
    box%matrix(1,1) = 2.0_real64
    box%matrix(2,2) = 3.0_real64
    box%matrix(3,3) = 4.0_real64

    call ComputeCellProperties(box)

    expected = 2.0_real64 * 3.0_real64 * 4.0_real64
    pass2 = abs(box%volume - expected) < tol
    if (.not. pass2) then
        print *, "Orthorhombic box test FAILED (volume = ", box%volume, ")"
        stop 1
    end if

    ! Check cosines are zero (orthogonal)
    if (.not. all(abs(box%metrics(4:6)) < tol)) then
        print *, "Orthorhombic cosines FAILED:", box%metrics(4:6)
        stop 1
    end if

    !======================================================
    ! Test 3: Triclinic box
    !======================================================
    ! Skewed box: a=(1,0,0), b=(0.5,1,0), c=(0,0.5,1)
    box%matrix(:,1) = [1.0_real64, 0.0_real64, 0.0_real64]
    box%matrix(:,2) = [0.5_real64, 1.0_real64, 0.0_real64]
    box%matrix(:,3) = [0.0_real64, 0.5_real64, 1.0_real64]

    call ComputeCellProperties(box)

    expected = box%matrix(1,1) * (box%matrix(2,2)*box%matrix(3,3) - box%matrix(2,3)*box%matrix(3,2)) &
             - box%matrix(1,2) * (box%matrix(2,1)*box%matrix(3,3) - box%matrix(2,3)*box%matrix(3,1)) &
             + box%matrix(1,3) * (box%matrix(2,1)*box%matrix(3,2) - box%matrix(2,2)*box%matrix(3,1))
    expected = abs(expected)

    pass3 = abs(box%volume - expected) < 1.0d-8
    if (.not. pass3) then
        print *, "Triclinic box test FAILED (volume = ", box%volume, ", expected=", expected, ")"
        stop 1
    end if

end program test_ComputeCellProperties
