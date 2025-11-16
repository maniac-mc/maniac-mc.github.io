program test_PickRandomMoleculeIndex

    use, intrinsic :: iso_fortran_env, only: real64
    use monte_carlo_utils
    use random_utils

    implicit none

    integer, parameter :: n_trials     = 2000
    integer, parameter :: n_molecules  = 7

    integer :: i, pick
    logical :: pass1, pass2, pass3
    integer, allocatable :: seen(:)

    !---------------------------------------
    ! Test 1: Valid molecule count (>0)
    !---------------------------------------
    pass1 = .true.
    do i = 1, n_trials
        pick = PickRandomMoleculeIndex(n_molecules)
        if (pick < 1 .or. pick > n_molecules) then
            pass1 = .false.
            exit
        end if
    end do

    !---------------------------------------
    ! Test 2: Zero molecules
    !---------------------------------------
    pick = PickRandomMoleculeIndex(0)
    pass2 = (pick == 0)

    !---------------------------------------
    ! Test 3: Coverage of all indices
    !---------------------------------------
    allocate(seen(n_molecules))
    seen = 0

    do i = 1, n_trials
        pick = PickRandomMoleculeIndex(n_molecules)
        if (pick >= 1 .and. pick <= n_molecules) then
            seen(pick) = 1
        end if
    end do

    pass3 = all(seen == 1)

    !---------------------------------------
    ! Final output
    !---------------------------------------
    if (pass1 .and. pass2 .and. pass3) then
        print *, 'PickRandomMoleculeIndex test PASSED'
    else
        print *, 'PickRandomMoleculeIndex test FAILED'
        if (.not. pass1) print *, ' Error: returned index out of range'
        if (.not. pass2) print *, ' Error: did not return 0 for zero molecules'
        if (.not. pass3) print *, ' Error: random sampling missed some indices'
        stop 1
    end if

end program test_PickRandomMoleculeIndex
