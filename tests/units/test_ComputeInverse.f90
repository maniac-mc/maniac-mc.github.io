program test_ComputeInverse

    use, intrinsic :: iso_fortran_env, only: real64
    use geometry_utils
    implicit none

    type(type_box) :: test_box
    real(real64), dimension(3,3) :: identity_check
    logical :: pass1, pass2, pass3
    integer :: i, j
    real(real64), parameter :: tol = 1.0D-10  ! tolerance for numerical checks

    !-------------------------------------------------------------
    ! Test 1: Cubic box (diagonal with equal sides)
    !-------------------------------------------------------------
    test_box%matrix = 0.0_real64
    test_box%matrix(1,1) = 2.0_real64
    test_box%matrix(2,2) = 2.0_real64
    test_box%matrix(3,3) = 2.0_real64

    call ComputeInverse(test_box)

    ! Check determinant
    pass1 = abs(test_box%determinant - 8.0_real64) < tol

    ! Check inverse correctness: A * A_inv ≈ I
    identity_check = matmul(test_box%matrix, test_box%reciprocal)
    do i = 1,3
        do j = 1,3
            if (i == j) then
                pass1 = pass1 .and. abs(identity_check(i,j) - 1.0_real64) < tol
            else
                pass1 = pass1 .and. abs(identity_check(i,j)) < tol
            end if
        end do
    end do

    if (.not. pass1) then
        print *, 'ComputeInverse : Cubic box test FAILED'
        stop 1
    end if

    !-------------------------------------------------------------
    ! Test 2: Orthorhombic box (diagonal with unequal sides)
    !-------------------------------------------------------------
    test_box%matrix = 0.0_real64
    test_box%matrix(1,1) = 1.0_real64
    test_box%matrix(2,2) = 2.0_real64
    test_box%matrix(3,3) = 3.0_real64

    call ComputeInverse(test_box)

    pass2 = abs(test_box%determinant - 6.0_real64) < tol
    identity_check = matmul(test_box%matrix, test_box%reciprocal)
    do i = 1,3
        do j = 1,3
            if (i == j) then
                pass2 = pass2 .and. abs(identity_check(i,j) - 1.0_real64) < tol
            else
                pass2 = pass2 .and. abs(identity_check(i,j)) < tol
            end if
        end do
    end do

    if (.not. pass2) then
        print *, 'ComputeInverse : Orthorhombic box test FAILED'
        stop 1
    end if

    !-------------------------------------------------------------
    ! Test 3: General (triclinic) box with off-diagonal elements
    !-------------------------------------------------------------
    test_box%matrix(:,1) = [2.0_real64, 0.1_real64, 0.05_real64]
    test_box%matrix(:,2) = [0.1_real64, 3.0_real64, 0.1_real64]
    test_box%matrix(:,3) = [0.05_real64, 0.1_real64, 4.0_real64]

    call ComputeInverse(test_box)

    pass3 = abs(test_box%determinant) > tol
    identity_check = matmul(test_box%matrix, test_box%reciprocal)
    do i = 1,3
        do j = 1,3
            if (i == j) then
                pass3 = pass3 .and. abs(identity_check(i,j) - 1.0_real64) < tol
            else
                pass3 = pass3 .and. abs(identity_check(i,j)) < tol
            end if
        end do
    end do

    if (.not. pass3) then
        print *, 'ComputeInverse : Triclinic box test FAILED'
        stop 1
    end if

    !-------------------------------------------------------------
    ! Final summary
    !-------------------------------------------------------------
    if (pass1 .and. pass2 .and. pass3) then
        print *, 'ComputeInverse test PASSED'
    end if

end program test_ComputeInverse
