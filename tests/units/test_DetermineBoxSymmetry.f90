program test_determine_box_symmetry

    use, intrinsic :: iso_fortran_env, only: real64
    use geometry_utils
    implicit none

    type(type_box) :: test_box
    logical :: pass1, pass2, pass3

    ! Test Cubic box
    test_box%cell%matrix = 0.0_real64
    test_box%cell%matrix(1,1) = 1.0_real64
    test_box%cell%matrix(2,2) = 1.0_real64
    test_box%cell%matrix(3,3) = 1.0_real64

    call determine_box_symmetry(test_box)
    pass1 = (test_box%cell%shape == 1)
    if (.not. pass1) then
        print *, 'determine_box_symmetry : Cubic box test FAILED: type =', test_box%cell%shape
        stop 1
    end if

    ! Test Orthorhombic box
    test_box%cell%matrix = 0.0_real64
    test_box%cell%matrix(1,1) = 1.0_real64
    test_box%cell%matrix(2,2) = 2.0_real64
    test_box%cell%matrix(3,3) = 3.0_real64

    call determine_box_symmetry(test_box)
    pass2 = (test_box%cell%shape == 1)
    if (.not. pass2) then
        print *, 'determine_box_symmetry : Orthorhombic box test FAILED: type =', test_box%cell%shape
        stop 1
    end if

    ! Test Triclinic box
    test_box%cell%matrix = 0.0_real64
    test_box%cell%matrix(1,1) = 1.0_real64
    test_box%cell%matrix(2,2) = 1.0_real64
    test_box%cell%matrix(3,3) = 1.0_real64
    test_box%cell%matrix(1,2) = 0.1_real64 ! Nonzero off-diagonal

    call determine_box_symmetry(test_box)
    pass3 = (test_box%cell%shape == 2)
    if (.not. pass3) then
        print *, 'determine_box_symmetry : Triclinic box test FAILED: type =', test_box%cell%shape
        stop 1
    end if

end program test_determine_box_symmetry
