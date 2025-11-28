module ewald_kvectors

    !----------------------------------------------------------------------
    ! Provides routines to precompute all reciprocal lattice vectors
    ! and associated quantities for the reciprocal-space part of the
    ! Ewald summation method.
    !----------------------------------------------------------------------

    use simulation_state
    use output_utils
    use constants
    use geometry_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !-------------------------------------------------------------------
    ! Generates all valid reciprocal lattice vectors for the system
    ! and stores them in ewald%kvectors. Also precomputes the normalized
    ! squared magnitude (k^2) and Cartesian squared magnitude.
    !-------------------------------------------------------------------
    subroutine precompute_valid_reciprocal_vectors()

        integer :: kx_idx, ky_idx, kz_idx   ! Loop indices for each k-component
        integer :: count                    ! Counter for valid k-vectors
        real(real64) :: k_squared           ! Normalized squared magnitude of k-vector
        real(real64) :: k_squared_mag       ! Squared magnitude of the k-vector
        real(real64), dimension(3,3) :: kvec_matrix ! Columns are reciprocal lattice vectors b1, b2, b3

        ! Store reciprocal lattice vectors as columns of a 3x3 matrix
        kvec_matrix = TWOPI * reshape(primary%cell%reciprocal, shape(kvec_matrix))

        ! Fill the array with the actual k-vectors
        count = 0
        do kx_idx = 0, ewald%param%kmax(1)
            do ky_idx = -ewald%param%kmax(2), ewald%param%kmax(2)
                do kz_idx = -ewald%param%kmax(3), ewald%param%kmax(3)

                    if (kx_idx == 0 .and. ky_idx == 0 .and. kz_idx == 0) cycle

                    ! Compute normalized k^2 again
                    k_squared = normalized_K_squared(kx_idx, ky_idx, kz_idx, ewald%param%kmax)

                    ! Skip invalid k-vectors
                    if (.not. check_valid_reciprocal_vector(k_squared)) cycle

                    k_squared_mag = compute_cartesian_k_squared(kx_idx, ky_idx, kz_idx, kvec_matrix)

                    ! Increment counter and store k-vector components
                    count = count + 1
                    ewald%kvectors(count)%kx = kx_idx
                    ewald%kvectors(count)%ky = ky_idx
                    ewald%kvectors(count)%kz = kz_idx
                    ewald%kvectors(count)%k_squared = k_squared
                    ewald%kvectors(count)%k_squared_mag = k_squared_mag

                    ! Precompute form factor for current kx index: 1 for zero, 2 otherwise
                    ewald%form_factor(count) = compute_symmetry_form_factor(kx_idx)
                end do
            end do
        end do

    end subroutine precompute_valid_reciprocal_vectors

    !-------------------------------------------------------------------
    ! Computes the squared magnitude of a k-vector in normalized index space:
    ! k^2_normalized = (kx/kmax_x)^2 + (ky/kmax_y)^2 + (kz/kmax_z)^2
    !-------------------------------------------------------------------
    pure function normalized_K_squared(kx, ky, kz, kmax) result(k_squared)

        ! Input arguments
        integer, intent(in) :: kx, ky, kz
        integer, intent(in) :: kmax(3)
        
        ! Output argument
        real(real64) :: k_squared

        ! Compute normalized squared magnitude in index space
        k_squared = (dble(kx)/dble(kmax(1)))**2 + &
                    (dble(ky)/dble(kmax(2)))**2 + &
                    (dble(kz)/dble(kmax(3)))**2

    end function normalized_K_squared

    !-------------------------------------------------------------------
    ! Returns the multiplicative symmetry factor for a k-index along one axis.
    ! 1 for zero index, 2 for non-zero.
    !-------------------------------------------------------------------
    pure function compute_symmetry_form_factor(idx) result(factor)

        ! Input arguments
        integer, intent(in) :: idx

        ! Output argument
        real(real64) :: factor

        if (idx == 0) then
            factor = one
        else
            factor = two
        end if

    end function compute_symmetry_form_factor

    !-------------------------------------------------------------------
    ! Determines whether a normalized k-vector is valid for inclusion
    ! in the reciprocal-space Ewald summation.
    !-------------------------------------------------------------------
    pure function check_valid_reciprocal_vector(k_squared) result(valid)

        ! Input arguments
        real(real64), intent(in) :: k_squared

        ! Output argument
        logical :: valid

        ! Reject near-zero k-vectors (avoid singularity at k=0)
        ! and any vectors outside the normalized unit sphere.
        valid = (abs(k_squared) >= error) .and. (k_squared <= one)

    end function check_valid_reciprocal_vector

    !-------------------------------------------------------------------
    ! Precomputes reciprocal-space weighting factors W(k) for the
    ! Ewald sum, using the standard formula:
    !   W(k) = exp(-|k|^2 / (4 * alpha^2)) / |k|^2
    !-------------------------------------------------------------------
    subroutine compute_reciprocal_weights()

        ! Local variables
        integer :: idx                  ! Loop index over precomputed k-vectors
        real(real64) :: k_squared_mag   ! Squared magnitude of the k-vector
        real(real64) :: alpha_squared   ! Precompute alpha^2 for efficiency (alpha = screening parameter)

        ! Calculate the square of the screening parameter
        alpha_squared = ewald%param%alpha**2

        ! Loop over all precomputed reciprocal lattice vectors
        do idx = 1, ewald%param%nkvec

            ! Compute the reciprocal-space weighting factor for this k-vector
            k_squared_mag = ewald%kvectors(idx)%k_squared_mag

            if (k_squared_mag <= 0.0_real64) then
                write(*,*) "FATAL: compute_reciprocal_weights: k_squared_mag <= 0 for idx=", idx, &
                            "  kx,ky,kz=", ewald%kvectors(idx)%kx, &
                            ewald%kvectors(idx)%ky, ewald%kvectors(idx)%kz
                call flush(0)
                stop 1
            end if

            ! Guard against extremely small k^2 that could cause overflow in 1/k^2
            if (k_squared_mag < 1.0e-300_real64) then
                ! Very small: set weight to 0 or a safe finite number
                ewald%kweights(idx) = 0.0_real64
            else
                ewald%kweights(idx) = exp(-k_squared_mag / (four * alpha_squared)) / k_squared_mag
            end if

        end do

    end subroutine compute_reciprocal_weights

    !-------------------------------------------------------------------
    ! Computes the squared magnitude of a reciprocal lattice vector
    ! in Cartesian space using the reciprocal lattice vectors.
    !-------------------------------------------------------------------
    pure function compute_cartesian_k_squared(kx, ky, kz, kvec_matrix) result(k2_mag)

        ! Input arguments
        integer, intent(in) :: kx, ky, kz
        real(real64), intent(in) :: kvec_matrix(3,3)

        ! Output argument
        real(real64) :: k2_mag

        ! Local variable
        real(real64) :: kvec(3)

        ! Build the 3D k-vector
        kvec = dble(kx) * kvec_matrix(:,1) + &
               dble(ky) * kvec_matrix(:,2) + &
               dble(kz) * kvec_matrix(:,3)

        ! Return squared magnitude
        k2_mag = dot_product(kvec, kvec)

    end function compute_cartesian_k_squared

end module ewald_kvectors
