program test_compute_cross_product

    use geometry_utils
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    real(real64) :: a(3), b(3), c(3), expected(3)
    real(real64), parameter :: tol = 1.0d-12
    logical :: pass1, pass2, pass3

    !---------------------------------------
    ! Test 1: Standard basis vectors
    ! a = x̂, b = ŷ → a × b = ẑ
    !---------------------------------------
    a = [1.0_real64, 0.0_real64, 0.0_real64]
    b = [0.0_real64, 1.0_real64, 0.0_real64]
    expected = [0.0_real64, 0.0_real64, 1.0_real64]

    c = compute_cross_product(a, b)
    pass1 = all(abs(c - expected) < tol)

    !---------------------------------------
    ! Test 2: Parallel vectors → cross product = 0
    ! a = [2,3,4], b = 2*a → a × b = 0
    !---------------------------------------
    a = [2.0_real64, 3.0_real64, 4.0_real64]
    b = 2.0_real64 * a
    expected = [0.0_real64, 0.0_real64, 0.0_real64]

    c = compute_cross_product(a, b)
    pass2 = all(abs(c - expected) < tol)

    !---------------------------------------
    ! Test 3: Arbitrary vectors
    ! a = [1,2,3], b = [4,5,6] → a×b = [-3,6,-3]
    !---------------------------------------
    a = [1.0_real64, 2.0_real64, 3.0_real64]
    b = [4.0_real64, 5.0_real64, 6.0_real64]
    expected = [-3.0_real64, 6.0_real64, -3.0_real64]

    c = compute_cross_product(a, b)
    pass3 = all(abs(c - expected) < tol)

    !---------------------------------------
    ! Final result
    !---------------------------------------
    if (.not. pass1) then
        print *, ' Test 1 failed: ', c
        stop 1
    end if
    if (.not. pass2) then
        print *, ' Test 2 failed: ', c
        stop 1
    end if
    if (.not. pass3) then
        print *, ' Test 3 failed: ', c
        stop 1
    end if

end program test_compute_cross_product
