module prepare_utils

    use simulation_state
    use output_utils
    use ewald_kvectors
    use ewald_phase
    use ewald_energy
    use tabulated_utils

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !-----------------------------------------------------------
    ! Convert fugacities to dimensionless activities and
    ! initialize Ewald summation parameters for the simulation.
    !-----------------------------------------------------------
    subroutine setup_simulation_parameters()

        ! Initialize Ewald summation parameters
        ! (cutoff, precision, reciprocal space, etc.)
        call setup_ewald()

        ! Allocate memory for arrays needed by the Ewald
        ! method (reciprocal vectors, coefficients, etc.)
        call allocate_array()

        ! Precompute k vectors
        call precompute_valid_reciprocal_vectors()

        ! Some precalculations for Monte Carlo
        call prepare_monte_carlo()

        ! Log the Ewald parameters and settings for
        ! reproducibility and debugging
        call log_ewald_parameters()

        ! Precompute tables for faster calculation
        call precompute_table()

    end subroutine setup_simulation_parameters

    !-------------------------------------------------------------------
    ! Allocate arrays used for Fourier components and related data
    !-------------------------------------------------------------------
    subroutine allocate_array()

        ! Local variable
        integer :: kmax_max

        allocate(saved%offset(3, nmax%atoms_per_residue))
        allocate(saved%com(3))

        ! Allocate real arrays for coefficients (dimension ewald%param%nkvec)
        allocate(ewald%kweights(1:ewald%param%nkvec))
        allocate(ewald%Ak(1:ewald%param%nkvec))
        allocate(ewald%Ak_old(1:ewald%param%nkvec))
        allocate(ewald%form_factor(ewald%param%nkvec)) ! Note, there is no need for such as large vector

        ! Allocate complex arrays for wave vector components
        kmax_max = maxval(ewald%param%kmax)
        allocate(ewald%phase%factor(3, res%number, 0:NB_MAX_MOLECULE, 1:nmax%atoms_per_residue, -kmax_max:kmax_max)) ! TOFIX, do not use NB_MAX_MOLECULE systematical ? 
        allocate(ewald%phase%factor_old(3, 1:nmax%atoms_per_residue, -kmax_max:kmax_max))

        ! Allocate temporary arrays once
        allocate(ewald%phase_axis(-kmax_max:kmax_max))
        allocate(ewald%phase%new(nmax%atoms_per_residue))
        allocate(ewald%phase%old(nmax%atoms_per_residue))
        allocate(ewald%q_buffer(nmax%atoms_per_residue))

        ! Allocate kvectors
        allocate(ewald%kvectors(ewald%param%nkvec))

        ! Allocate widom
        if (proba%widom > 0) then
            allocate(statistic%weight(res%number))
            allocate(statistic%sample(res%number))
            allocate(statistic%mu_ex(res%number))
            allocate(statistic%mu_tot(res%number))
            statistic%weight(:) = 0
            statistic%sample(:) = 0
        end if

    end subroutine allocate_array

    !--------------------------------------------------------------------
    ! Ensure that the user-specified real-space cutoff fits
    ! within the periodic simulation box. If it is too large,
    ! reduce it to half the minimum box dimension.
    !--------------------------------------------------------------------
    subroutine adjust_real_space_cutoff(do_log)

        ! Input parameter
        logical, intent(in) :: do_log

        ! Local variables
        character(200) :: msg           ! Buffer for logging
        real(real64) :: safe_cutoff     ! Adjusted real-space cutoff length to fit inside the simulation box safely

        if (mc_input%real_space_cutoff > primary%cell%metrics(1) &
                .or. mc_input%real_space_cutoff > primary%cell%metrics(2) &
                .or. mc_input%real_space_cutoff > primary%cell%metrics(3)) then

            if (do_log) then

                write(msg, '(A)') 'WARNING: real_space_cutoff too large for box. Reducing to safe value.'
                call log_message(msg)
            
            end if
            
            safe_cutoff = min(primary%cell%metrics(1), primary%cell%metrics(2), primary%cell%metrics(3)) / two
            mc_input%real_space_cutoff = safe_cutoff
        
        end if

    end subroutine adjust_real_space_cutoff

    !--------------------------------------------------------------------
    ! Compute the maximum Fourier indices (kmax) in the X, Y, Z
    ! directions, and the total number of reciprocal lattice vectors.
    !--------------------------------------------------------------------
    subroutine compute_fourier_indices()

        integer :: kx_idx, ky_idx, kz_idx   ! Loop indices for each k-component
        integer :: count                    ! Counter for valid k-vectors
        real(real64) :: k_squared           ! Normalized squared magnitude of k-vector

        ! Compute maximum Fourier indices in X, Y, Z directions
        ewald%param%kmax = nint(quarter + primary%cell%metrics(1:3) * ewald%param%alpha * ewald%param%fprecision / PI)

        ! Count the total number of valid reciprocal vectors
        count = 0
        do kx_idx = 0, ewald%param%kmax(1)
            do ky_idx = -ewald%param%kmax(2), ewald%param%kmax(2)
                do kz_idx = -ewald%param%kmax(3), ewald%param%kmax(3)

                    if (kx_idx == 0 .and. ky_idx == 0 .and. kz_idx == 0) cycle
                
                    k_squared = normalized_K_squared(kx_idx, ky_idx, kz_idx, ewald%param%kmax)
                    if (check_valid_reciprocal_vector(k_squared)) then
                        count = count + 1
                    end if
                
                end do
            end do
        end do

        ewald%param%nkvec = count

    end subroutine compute_fourier_indices

    !--------------------------------------------------------------------
    ! Initialize parameters for the Ewald summation method.
    !--------------------------------------------------------------------
    subroutine setup_ewald(verbose)

        ! Input parameter
        logical, optional, intent(in) :: verbose   ! optional flag for logging

        ! Local variable
        logical :: do_log                          ! Internal flag: whether to log messages in this routine

        ! Default: log unless explicitly disabled
        do_log = .true.; if (present(verbose)) do_log = verbose

        ! Adjust cutoff if too large
        call adjust_real_space_cutoff(do_log)

        ! Clamp tolerance
        call clamp_tolerance()

        ! Compute Ewald parameters
        call compute_ewald_parameters()

        ! Compute Fourier indices
        call compute_fourier_indices()

    end subroutine setup_ewald

    !--------------------------------------------------------------------
    ! Limit the Ewald accuracy tolerance to a maximum value of 0.5.
    !--------------------------------------------------------------------
    subroutine clamp_tolerance()

        ! Clamp accuracy tolerance to max 0.5
        ewald%param%tolerance = min(abs(ewald%param%tolerance), half)
    
    end subroutine clamp_tolerance

    !--------------------------------------------------------------------
    ! Compute the main Ewald summation parameters; the screening factor, the damping parameter (alpha)
    ! and the fourier-space precision
    !--------------------------------------------------------------------
    subroutine compute_ewald_parameters()
    
        ! Intermediate tolerance factor for screening width
        ewald%param%screen = sqrt(abs(log(ewald%param%tolerance * mc_input%real_space_cutoff)))

        ! Compute Ewald damping parameter
        ewald%param%alpha = sqrt(abs(log(ewald%param%tolerance * mc_input%real_space_cutoff * ewald%param%screen))) / &
                    mc_input%real_space_cutoff

        ! Estimate needed Fourier-space precision
        ewald%param%fprecision = sqrt(-log(ewald%param%tolerance * mc_input%real_space_cutoff * &
                                (two * ewald%param%screen * ewald%param%alpha)**2))

    end subroutine compute_ewald_parameters

    !-------------------------------------------------------------------
    ! Compute beta and de Broglie
    !-------------------------------------------------------------------
    subroutine prepare_monte_carlo()

        ! Local variables
        real(real64) :: mass                   ! Mass of residue
        integer :: val_int                     ! Integer value read from input

        ! Compute inverse thermal energy β = 1/(k_B T)
        beta = 1/(KB_kcalmol*thermo%temperature) ! 1/(kB T) in 1/(kcal/mol)

        do val_int = 1, res%number

            if (.not. thermo%is_active(val_int)) cycle

            if (thermo%fugacity(val_int) >= zero) then
                thermo%chemical_potential(val_int) = log(thermo%fugacity(val_int)) / beta
            end if

            ! Compute thermal de Broglie wavelength λ for each active residue:
            ! λ = h / sqrt(2 π m k_B T)
            ! H_PLANCK in J s, mass in kg, KB in J/K, T in K
            mass = res%mass(val_int) * G_TO_KG / NA                                         ! Mass per residue (kg)
            res%lambda(val_int) = H_PLANCK / sqrt(TWOPI * mass * KB * thermo%temperature)    ! Thermal de Broglie wavelength
            
            ! Convert lambda to Å
            res%lambda(val_int) = res%lambda(val_int) * M_TO_A

        end do

    end subroutine prepare_monte_carlo

end module prepare_utils
