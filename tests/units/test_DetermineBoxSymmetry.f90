program test_DetermineBoxSymmetry

    use, intrinsic :: iso_fortran_env, only: real64
    use geometry_utils
    implicit none

    type(type_box) :: test_box
    logical :: pass1, pass2, pass3

    ! Test Cubic box
    test_box%matrix = 0.0_real64
    test_box%matrix(1,1) = 1.0_real64
    test_box%matrix(2,2) = 1.0_real64
    test_box%matrix(3,3) = 1.0_real64

    call DetermineBoxSymmetry(test_box)
    pass1 = (test_box%type == 1)
    if (.not. pass1) then
        print *, 'DetermineBoxSymmetry : Cubic box test FAILED: type =', test_box%type
        stop 1
    end if

    ! Test Orthorhombic box
    test_box%matrix = 0.0_real64
    test_box%matrix(1,1) = 1.0_real64
    test_box%matrix(2,2) = 2.0_real64
    test_box%matrix(3,3) = 3.0_real64

    call DetermineBoxSymmetry(test_box)
    pass2 = (test_box%type == 2)
    if (.not. pass2) then
        print *, 'DetermineBoxSymmetry : Orthorhombic box test FAILED: type =', test_box%type
        stop 1
    end if

    ! Test Triclinic box
    test_box%matrix = 0.0_real64
    test_box%matrix(1,1) = 1.0_real64
    test_box%matrix(2,2) = 1.0_real64
    test_box%matrix(3,3) = 1.0_real64
    test_box%matrix(1,2) = 0.1_real64 ! Nonzero off-diagonal

    call DetermineBoxSymmetry(test_box)
    pass3 = (test_box%type == 3)
    if (.not. pass3) then
        print *, 'DetermineBoxSymmetry : Triclinic box test FAILED: type =', test_box%type
        stop 1
    end if

    if ((pass3) .and. (pass2) .and. (pass3)) then
        print *, 'DetermineBoxSymmetry test PASSED'
    end if

end program test_DetermineBoxSymmetry
