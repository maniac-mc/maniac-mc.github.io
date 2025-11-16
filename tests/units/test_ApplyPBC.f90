program test_ApplyPBC
    use, intrinsic :: iso_fortran_env, only: real64
    use geometry_utils
    implicit none

    type(type_box) :: box
    real(real64), dimension(3) :: pos, expected
    logical :: pass1, pass2, pass3
    real(real64) :: tol
    real(real64), dimension(3,3) :: recip
    integer :: i

    tol = 1.0d-12

    !======================================================
    ! Test 1: Cubic box
    !======================================================
    box%type = 1
    box%is_triclinic = .false.
    box%bounds(:,1) = [0.0_real64, 0.0_real64, 0.0_real64]
    box%bounds(:,2) = [1.0_real64, 1.0_real64, 1.0_real64]
    box%matrix = 0.0_real64
    box%matrix(1,1) = 1.0_real64
    box%matrix(2,2) = 1.0_real64
    box%matrix(3,3) = 1.0_real64

    ! Position outside [0,1)
    pos = [1.2_real64, -0.3_real64, 2.7_real64]
    call ApplyPBC(pos, box)

    expected = modulo(pos, 1.0_real64)
    pass1 = all(abs(pos - expected) < tol)
    if (.not. pass1) then
        print *, "Cubic box PBC test FAILED: pos = ", pos, " expected = ", expected
        stop 1
    end if

    !======================================================
    ! Test 2: Orthorhombic box
    !======================================================
    box%type = 2
    box%is_triclinic = .false.
    box%bounds(:,1) = [0.0_real64, -1.5_real64, 0.0_real64]
    box%bounds(:,2) = [2.0_real64, 1.5_real64, 4.0_real64]
    box%matrix = 0.0_real64
    box%matrix(1,1) = 2.0_real64
    box%matrix(2,2) = 3.0_real64
    box%matrix(3,3) = 4.0_real64

    pos = [2.5_real64, -1.6_real64, 9.0_real64]
    call ApplyPBC(pos, box)

    do i = 1,3
        expected(i) = box%bounds(i,1) + modulo(pos(i) - box%bounds(i,1), box%matrix(i,i))
    end do

    pass2 = all(abs(pos - expected) < tol)
    if (.not. pass2) then
        print *, "Orthorhombic box PBC test FAILED: pos = ", pos, " expected = ", expected
        stop 1
    end if

    !======================================================
    ! Test 3: Triclinic box
    !======================================================
    box%type = 3
    box%is_triclinic = .true.
    ! Define lower bounds
    box%bounds(:,1) = [0.0_real64, 0.0_real64, 0.0_real64]
    ! Box vectors: a=(1,0,0), b=(0.5,1,0), c=(0,0.5,1)
    box%matrix(:,1) = [1.0_real64, 0.0_real64, 0.0_real64]
    box%matrix(:,2) = [0.5_real64, 1.0_real64, 0.0_real64]
    box%matrix(:,3) = [0.0_real64, 0.5_real64, 1.0_real64]

    ! Compute reciprocal matrix
    recip = inverse(box%matrix)
    box%reciprocal = recip

    pos = [1.2_real64, 0.3_real64, -0.6_real64]
    call ApplyPBC(pos, box)

    ! Compute expected using fractional coordinates
    expected = box%bounds(:,1) + matmul(box%matrix, modulo(matmul(box%reciprocal, pos - box%bounds(:,1)), 1.0_real64))

    pass3 = all(abs(pos - expected) < tol)
    if (.not. pass3) then
        print *, "Triclinic box PBC test FAILED: pos = ", pos, " expected = ", expected
        stop 1
    end if

    !-------------------------------------------------------------
    ! Final summary
    !-------------------------------------------------------------
    if (pass1 .and. pass2 .and. pass3) then
        print *, 'ApplyPBC test PASSED'
    end if

contains
    !-----------------------------------------------------------------------
    ! Simple 3x3 matrix inverse (for triclinic test)
    !-----------------------------------------------------------------------
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

end program test_ApplyPBC
