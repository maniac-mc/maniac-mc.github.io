program test_pick_random_residue_type

    use, intrinsic :: iso_fortran_env, only: real64
    use monte_carlo_utils
    use random_utils

    implicit none

    integer :: i, pick
    logical, allocatable :: is_active(:)
    logical :: pass1, pass2
    integer, parameter :: n_trials = 1000

    !---------------------------------------
    ! Test 1: Some active residues
    !---------------------------------------
    allocate(is_active(5))
    is_active = [.true., .false., .true., .false., .true.]  ! residues 1,3,5 are active
    pass1 = .true.

    do i = 1, n_trials
        pick = pick_random_residue_type(is_active)
        if (.not. any(pick == [1,3,5])) then
            pass1 = .false.
            exit
        end if
    end do

    !---------------------------------------
    ! Test 2: No active residues
    !---------------------------------------
    is_active = [.false., .false., .false., .false., .false.]
    pick = pick_random_residue_type(is_active)
    pass2 = (pick == 0)

    !---------------------------------------
    ! Final result
    !---------------------------------------
    if (.not. pass1) then
        print *, 'pick_random_residue_type test FAILED'
        print *, ' Error: function returned inactive residue'
        stop 1
    end if

    if (.not. pass2) then
        print *, 'pick_random_residue_type test FAILED'
        print *, ' Error: function did not return 0 for no active residues'
        stop 1
    end if

end program test_pick_random_residue_type
