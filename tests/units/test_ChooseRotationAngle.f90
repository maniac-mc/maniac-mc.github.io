program test_ChooseRotationAngle

    use iso_fortran_env, only: real64
    use simulation_state
    use monte_carlo_utils

    implicit none

    real(real64) :: theta

    !--------------------------------------------
    ! Setup input parameters
    !--------------------------------------------
    input%rotation_step_angle = 0.5_real64  ! small-step max rotation (radians)

    !--------------------------------------------
    ! Test 1: small-step rotation
    !--------------------------------------------
    theta = ChooseRotationAngle(.false.)
    if (theta < -input%rotation_step_angle - error .or. theta > input%rotation_step_angle + error) then
        print *, "FAILED: small-step rotation out of range."
        stop 1
    end if

    !--------------------------------------------
    ! Test 2: full rotation [0, 2π]
    !--------------------------------------------
    theta = ChooseRotationAngle(.true.)
    if (theta < 0.0_real64 - error .or. theta > 2.0_real64*acos(-1.0_real64) + error) then
        print *, "FAILED: full rotation out of range."
        stop 1
    end if

end program test_ChooseRotationAngle
