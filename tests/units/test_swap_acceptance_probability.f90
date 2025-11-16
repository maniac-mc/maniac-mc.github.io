program test_swap_acceptance_probability
    
    use iso_fortran_env, only: real64
    use simulation_state
    use constants
    use parameters
    use monte_carlo_utils

    implicit none

    type(energy_state) :: old_state, new_state
    real(real64) :: prob, expected
    real(real64) :: beta, delta_e
    real(real64) :: N_old, N_new, phi_old, phi_new, T
    real(real64), parameter :: tol = 1.0e-12_real64

    integer :: type_old, type_new

    !--------------------------------------------------------------
    ! Construct a controlled scenario with 2 species
    !--------------------------------------------------------------
    allocate(primary%num_residues(2))
    allocate(input%fugacity(2))

    type_old = 1
    type_new = 2

    primary%num_residues = [10, 5]           ! N_old = 10, N_new = 5
    input%fugacity       = [0.02_real64, 0.05_real64]  ! φ_old < φ_new
    input%temp_K = 300.0_real64

    old_state%total = 6.0_real64
    new_state%total = 8.0_real64       ! ΔE = +2

    N_old = 10.0_real64
    N_new = 5.0_real64
    phi_old = 0.02_real64
    phi_new = 0.05_real64
    T = input%temp_K

    delta_e = new_state%total - old_state%total
    beta = 1.0_real64 / (KB_kcalmol * T)

    !==============================================================
    ! 1. BASELINE TEST (positive ΔE, φ_new > φ_old)
    !==============================================================
    expected = min(one, &
        (phi_new / phi_old) * (N_old / (N_new + one)) * exp(-beta * delta_e) )

    prob = swap_acceptance_probability(old_state, new_state, type_old, type_new)

    if (abs(prob - expected) > tol) then
        print *, "FAILED: Basic swap acceptance probability incorrect."
        print *, "Expected =", expected
        print *, "Got      =", prob
        stop 1
    end if


    !==============================================================
    ! 2. NEGATIVE ΔE TEST — MUST ACCEPT
    !==============================================================
    new_state%total = 5.0_real64      ! ΔE = -1

    prob = swap_acceptance_probability(old_state, new_state, type_old, type_new)

    if (prob < 0.999999_real64) then
        print *, "FAILED: Negative ΔE should yield P = 1."
        print *, "Got =", prob
        stop 1
    end if


    !==============================================================
    ! 3. EDGE CASE: N_new = 0 (rare but must work)
    !==============================================================
    primary%num_residues(type_new) = 0
    N_new = 0.0_real64

    new_state%total = 8.0_real64      ! restore ΔE = +2

    expected = min(one, &
        (phi_new / phi_old) * (N_old / (N_new + one)) * exp(-beta * delta_e))

    prob = swap_acceptance_probability(old_state, new_state, type_old, type_new)

    if (abs(prob - expected) > tol) then
        print *, "FAILED: N_new = 0 case incorrect."
        print *, "Expected =", expected
        print *, "Got      =", prob
        stop 1
    end if


    !==============================================================
    ! 4. EDGE CASE: N_old = 1 (last molecule of its type)
    !==============================================================
    primary%num_residues(type_old) = 1
    N_old = 1.0_real64

    expected = min(one, &
        (phi_new / phi_old) * (N_old / (N_new + one)) * exp(-beta * delta_e))

    prob = swap_acceptance_probability(old_state, new_state, type_old, type_new)

    if (abs(prob - expected) > tol) then
        print *, "FAILED: N_old = 1 case incorrect."
        print *, "Expected =", expected
        print *, "Got      =", prob
        stop 1
    end if

end program test_swap_acceptance_probability
