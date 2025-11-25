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

        allocate(res%site_offset_old(3, nb%max_atom_in_any_residue))

        ! Allocate real arrays for coefficients (dimension ewald%num_kvectors)
        allocate(ewald%recip_constants(1:ewald%num_kvectors))
        allocate(ewald%recip_amplitude(1:ewald%num_kvectors))
        allocate(ewald%recip_amplitude_old(1:ewald%num_kvectors))
        allocate(ewald%form_factor(ewald%num_kvectors)) ! Note, there is no need for such as large vector

        ! Allocate complex arrays for wave vector components
        kmax_max = maxval(ewald%kmax)
        allocate(ewald%phase_factor(3, nb%type_residue, 0:NB_MAX_MOLECULE, 1:nb%max_atom_in_any_residue, -kmax_max:kmax_max))
        allocate(ewald%phase_factor_old(3, 1:nb%max_atom_in_any_residue, -kmax_max:kmax_max))

        allocate(ewald%temp_1d(-kmax_max:kmax_max))

        ! Allocate temporary arrays once
        allocate(ewald%temp(3, -kmax_max:kmax_max))
        allocate(ewald%phase_new(nb%max_atom_in_any_residue))
        allocate(ewald%phase_old(nb%max_atom_in_any_residue))
        allocate(ewald%charges(nb%max_atom_in_any_residue))

        ! Allocate kvectors
        allocate(ewald%kvectors(ewald%num_kvectors))

        ! Allocate widom
        if (proba%widom > 0) then
            allocate(statistic%weight(nb%type_residue))
            allocate(statistic%sample(nb%type_residue))
            allocate(statistic%mu_ex(nb%type_residue))
            allocate(statistic%mu_tot(nb%type_residue))
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

        if (input%real_space_cutoff > primary%metrics(1) &
                .or. input%real_space_cutoff > primary%metrics(2) &
                .or. input%real_space_cutoff > primary%metrics(3)) then

            if (do_log) then

                write(msg, '(A)') 'WARNING: real_space_cutoff too large for box. Reducing to safe value.'
                call log_message(msg)
            
            end if
            
            safe_cutoff = min(primary%metrics(1), primary%metrics(2), primary%metrics(3)) / two
            input%real_space_cutoff = safe_cutoff
        
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
        ewald%kmax = nint(quarter + primary%metrics(1:3) * ewald%alpha * ewald%fourier_precision / PI)

        ! Count the total number of valid reciprocal vectors
        count = 0
        do kx_idx = 0, ewald%kmax(1)
            do ky_idx = -ewald%kmax(2), ewald%kmax(2)
                do kz_idx = -ewald%kmax(3), ewald%kmax(3)

                    if (kx_idx == 0 .and. ky_idx == 0 .and. kz_idx == 0) cycle
                
                    k_squared = normalized_K_squared(kx_idx, ky_idx, kz_idx, ewald%kmax)
                    if (check_valid_reciprocal_vector(k_squared)) then
                        count = count + 1
                    end if
                
                end do
            end do
        end do

        ewald%num_kvectors = count

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
        ewald%tolerance = min(abs(ewald%tolerance), half)
    
    end subroutine clamp_tolerance

    !--------------------------------------------------------------------
    ! Compute the main Ewald summation parameters; the screening factor, the damping parameter (alpha)
    ! and the fourier-space precision
    !--------------------------------------------------------------------
    subroutine compute_ewald_parameters()
    
        ! Intermediate tolerance factor for screening width
        ewald%screening_factor = sqrt(abs(log(ewald%tolerance * input%real_space_cutoff)))

        ! Compute Ewald damping parameter
        ewald%alpha = sqrt(abs(log(ewald%tolerance * input%real_space_cutoff * ewald%screening_factor))) / &
                    input%real_space_cutoff

        ! Estimate needed Fourier-space precision
        ewald%fourier_precision = sqrt(-log(ewald%tolerance * input%real_space_cutoff * &
                                (two * ewald%screening_factor * ewald%alpha)**2))

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

        do val_int = 1, nb%type_residue

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
