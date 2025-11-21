module random_utils

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    ! Preallocate put array for RNG seeding
    integer, allocatable :: seed_array(:)

contains

    !------------------------------------------------------------------------------
    ! Returns double precision random number uniformly distributed in [0,1)
    !------------------------------------------------------------------------------
    function rand_uniform() result(val)

        ! Output variable
        real(real64) :: val                         ! Result in [-0.5, 0.5)

        call random_number(val)

    end function rand_uniform

    !------------------------------------------------------------------------------
    ! Returns an array of double precision random numbers in [-0.5, 0.5)
    !------------------------------------------------------------------------------
    function rand_symmetric(n) result(vals)

        ! Input variable
        integer, intent(in) :: n                    ! Number of values to generate

        ! Output variable
        real(real64) :: vals(n)                     ! Result array in [-0.5, 0.5)

        call random_number(vals)                    ! Fills entire array with [0,1)
        vals = vals - 0.5_real64                    ! Shift to [-0.5,0.5)

    end function rand_symmetric

    !------------------------------------------------------------------------------
    ! Seeds the RNG: if seed=0, generate a random seed from system_clock
    !------------------------------------------------------------------------------
    subroutine seed_rng(seed)

        integer, intent(in) :: seed                 ! User-provided seed (0 = auto)
        integer :: seed_len                         ! Size of RNG seed array
        integer :: seed_idx                         ! Index for filling seed array
        integer :: useed                            ! Effective seed value
        integer :: count                            ! System clock counter
        integer :: count_rate                       ! Clock ticks per second
        integer :: count_max                        ! Maximum clock count before rollover

        if (seed == 0) then
            ! Get system clock count as a seed
            call system_clock(count, count_rate, count_max)
            useed = abs(mod(count, 2147483647))     ! keep in 32-bit integer range
        else
            useed = seed
        end if

        call random_seed(size=seed_len)

        ! Allocate seed_array only once
        if (.not. allocated(seed_array)) allocate(seed_array(seed_len))

        seed_array = useed + 104729 * [(seed_idx-1, seed_idx=1, seed_len)]
        call random_seed(put=seed_array)

    end subroutine seed_rng

end module random_utils
