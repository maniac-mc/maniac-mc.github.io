program test_PickRandomResidueType

    use, intrinsic :: iso_fortran_env, only: real64
    use monte_carlo_utils
    use random_utils

    implicit none

    integer, allocatable :: is_active(:)
    integer :: i, pick
    logical :: pass1, pass2
    integer, parameter :: n_trials = 1000

    !---------------------------------------
    ! Test 1: Some active residues
    !---------------------------------------
    allocate(is_active(5))
    is_active = [1, 0, 1, 0, 1]  ! residues 1,3,5 are active
    pass1 = .true.

    do i = 1, n_trials
        pick = PickRandomResidueType(is_active)
        if (.not. any(pick == [1,3,5])) then
            pass1 = .false.
            exit
        end if
    end do

    !---------------------------------------
    ! Test 2: No active residues
    !---------------------------------------
    is_active = [0,0,0,0,0]
    pick = PickRandomResidueType(is_active)
    pass2 = (pick == 0)

    !---------------------------------------
    ! Final result
    !---------------------------------------
    if (pass1 .and. pass2) then
        print *, 'PickRandomResidueType test PASSED'
    else
        print *, 'PickRandomResidueType test FAILED'
        if (.not. pass1) print *, ' Error: function returned inactive residue'
        if (.not. pass2) print *, ' Error: function did not return 0 for no active residues'
        stop 1
    end if

end program test_PickRandomResidueType
