program test_compute_acceptance_probability

    use iso_fortran_env, only: real64
    use simulation_state
    use constants
    use parameters
    use monte_carlo_utils

    implicit none

    type(energy_state) :: old_state, new_state
    real(real64) :: prob, expected, beta, delta_e, N, V, phi, T
    integer :: residue_type

    ! ========== Prepare minimal input ==========
    residue_type = 1
    allocate(primary%num_residues(1))
    allocate(input%fugacity(1))

    primary%num_residues(1) = 10
    primary%volume = 1000.0_real64
    input%fugacity(1) = 0.01_real64
    input%temp_K = 300.0_real64

    old_state%total = 5.0_real64
    new_state%total = 7.0_real64

    ! Local copies to compute expected values
    N = real(primary%num_residues(1), real64)
    V = primary%volume
    phi = input%fugacity(1)
    T = input%temp_K

    delta_e = new_state%total - old_state%total   ! = +2
    beta = 1.0_real64/(KB_kcalmol * T)

    ! ======================================================
    ! 1. TRANSLATION / ROTATION TEST
    ! Expected: exp(-beta * ΔE)
    ! ======================================================
    expected = exp(-beta * delta_e)

    prob = compute_acceptance_probability(old_state, new_state, residue_type, TYPE_TRANSLATION)

    if (abs(prob - expected) > error) then
        print *, "FAILED: Translation probability incorrect."
        print *, "Expected =", expected
        print *, "Got      =", prob
        stop 1
    end if

    ! ======================================================
    ! 2. CREATION TEST
    ! Expected: min(1, phi*V/(N+1) * exp(-beta * ΔE))
    ! ======================================================
    expected = min(one, (phi*V/(N+one)) * exp(-beta * delta_e))

    prob = compute_acceptance_probability(old_state, new_state, residue_type, TYPE_CREATION)

    if (abs(prob - expected) > error) then
        print *, "FAILED: Creation probability incorrect."
        print *, "Expected =", expected
        print *, "Got      =", prob
        stop 1
    end if

    ! ======================================================
    ! 3. NEGATIVE ΔE TEST (should always accept)
    ! ======================================================
    new_state%total = 4.0_real64          ! ΔE = -1

    prob = compute_acceptance_probability(old_state, new_state, residue_type, TYPE_TRANSLATION)

    if (prob < 0.999999_real64) then
        print *, "FAILED: Negative ΔE should produce P=1."
        print *, "Got =", prob
        stop 1
    end if

end program test_compute_acceptance_probability
