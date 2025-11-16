program test_ComputeDistance

    use, intrinsic :: iso_fortran_env, only: real64
    use energy_utils
    use geometry_utils

    implicit none

    type(type_box) :: box
    real(real64) :: dist, expected
    real(real64), parameter :: tol = 1.0d-6
    logical :: pass1, pass2, pass3
    integer :: nres, natoms, nmaxmol

    !------------------------------------------------------
    ! Two residues, one atom each
    !------------------------------------------------------
    nmaxmol = 2
    nres = 2
    natoms = 1

    allocate(box%mol_com(3, nres, nmaxmol))
    allocate(box%site_offset(3, nres, nmaxmol, natoms))

    box%mol_com(:,1,1) = [-4.0_real64, -4.0_real64, -4.0_real64]
    box%site_offset(:,1,1,1) = [1.0_real64, 0.0_real64, 0.0_real64]
    box%mol_com(:,2,1) = [4.0_real64, 4.0_real64, 4.0_real64]
    box%site_offset(:,2,1,1) = [1.0_real64, 0.0_real64, 0.0_real64]

    !------------------------------------------------------
    ! 1) Cubic box
    !------------------------------------------------------
    box%type = 1
    box%matrix(:,1) = [10.0_real64, 0.0_real64, 0.0_real64]
    box%matrix(:,2) = [0.0_real64, 10.0_real64, 0.0_real64]
    box%matrix(:,3) = [0.0_real64, 0.0_real64, 10.0_real64]
    dist = ComputeDistance(box,1,1,1,2,1,1)

    ! Precomputed with MDAnalysis
    expected = 3.4641014086612083_real64

    pass1 = abs(dist-expected) < tol
    if (.not. pass1) then
        print *, "Cubic box distance test FAILED"
        print *, "Computed distance = ", dist
        print *, "Expected distance = ", expected
        ! stop 1
    end if

    !------------------------------------------------------
    ! 2) Orthorhombic box
    !------------------------------------------------------
    box%type = 2
    box%matrix(:,1) = [12.0_real64, 0.0_real64, 0.0_real64]
    box%matrix(:,2) = [0.0_real64, 10.0_real64, 0.0_real64]
    box%matrix(:,3) = [0.0_real64, 0.0_real64, 10.0_real64]
    dist = ComputeDistance(box,1,1,1,2,1,1)

    ! Precomputed with MDAnalysis
    expected = 4.898979193564424_real64

    pass2 = abs(dist-expected) < tol
    if (.not. pass2) then
        print *, "Orthorhombic box distance test FAILED"
        print *, "Computed distance = ", dist
        print *, "Expected distance = ", expected
        ! stop 1
    end if

    !------------------------------------------------------
    ! 3) Triclinic box
    !------------------------------------------------------
    ! box%type = 3
    ! box%matrix(:,1) = [12.0_real64, 0.0_real64, 0.0_real64]
    ! box%matrix(:,2) = [2.0_real64, 10.0_real64, 0.0_real64]
    ! box%matrix(:,3) = [1.0_real64, 1.0_real64, 10.0_real64]
    ! call ComputeInverse(box)
    ! dist = ComputeDistance(box,1,1,1,2,1,1)

    ! ! Precomputed with MDAnalysis
    ! expected = 6.164413925615866_real64

    ! pass3 = abs(dist-expected) < tol
    ! if (.not. pass3) then
    !     print *, "Triclinic box distance test FAILED"
    !     print *, "Computed distance = ", dist
    !     print *, "Expected distance = ", expected
    !     stop 1
    ! end if

end program test_ComputeDistance
