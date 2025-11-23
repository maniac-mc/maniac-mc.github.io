module input_parser

    use parameters
    use simulation_state
    use output_utils
    use check_utils
    use readers_utils
    use random_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !-----------------------------------------------------------------------------
    ! Reads the main MANIAC input file in two passes (prescan + actual read).
    !-----------------------------------------------------------------------------
    subroutine read_input_file()

        ! Pre-scan input file
        call prescan_input_file(path%input)

        ! Allocate required arrays
        call allocate_atom_arrays()
        
        ! Parse input file and populate data
        call ReadFullInputFile(path%input)
        
        ! Validate and rescale move probabilities
        call ValidateAndRescaleMoveProbabilities()
        
        ! Print summary to log
        call PrintInputSummary()

    end subroutine read_input_file

    !---------------------------------------------------------------------------
    ! Opens the full MANIAC input file, reads and parses its contents by calling
    ! ParseInputFile, then closes the file.
    !---------------------------------------------------------------------------
    subroutine ReadFullInputFile(filename)

        ! Input parameter
        character(len=*), intent(in) :: filename

        ! Local variables
        integer :: unit     ! Fortran unit number used to open the input file
        integer :: ios      ! I/O status returned by open/read operations

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
            call check_IO_status(filename, ios)
            call ParseInputFile(unit)
        close(unit)

    end subroutine ReadFullInputFile

    !---------------------------------------------------------------------------
    ! Open the main input file, pre-detect sizes and counts for allocation,
    ! then close the file.
    !---------------------------------------------------------------------------
    subroutine prescan_input_file(filename)

        ! Input parameter
        character(len=*), intent(in) :: filename

        ! Local variables
        integer :: unit                             ! Fortran unit number used to open the input file
        integer :: ios                              ! I/O status returned by open/read operations

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)

            call check_IO_status(filename, ios)
            call predetect_number_info(unit)
        
        close(unit)

    end subroutine prescan_input_file

    !---------------------------------------------------------------------------
    ! Ensures the sum of move probabilities equals 1.0.
    ! If necessary, rescales all probabilities and warns the user.
    ! Aborts execution if probabilities are invalid.
    !---------------------------------------------------------------------------
    subroutine ValidateAndRescaleMoveProbabilities()

        real(real64) :: proba_total         ! Sum of all move probabilities (translation, rotation, etc.)
        real(real64) :: scale_factor        ! Factor used to rescale probabilities so they sum to 1.0

        ! Compute total probability
        proba_total = proba%translation + proba%rotation + proba%insertion_deletion + proba%swap + proba%widom

        ! Rescale if not exactly 1.0
        if (abs(proba_total - one) > error) then

            scale_factor = one / proba_total
            proba%translation = proba%translation * scale_factor
            proba%rotation = proba%rotation * scale_factor
            proba%insertion_deletion = proba%insertion_deletion * scale_factor
            proba%swap = proba%swap * scale_factor
            proba%widom = proba%widom * scale_factor
        
            call WarnUser("Move probabilities rescaled to sum to 1.0")
        
        end if

        ! Recompute total to validate
        proba_total = proba%translation + proba%rotation + proba%insertion_deletion + proba%swap + proba%widom

        ! Abort if proba invalid
        if (abs(proba_total - one) > error) then
            call AbortRun("Invalid move probabilities: must sum to 1.0", 1)
        end if

    end subroutine ValidateAndRescaleMoveProbabilities

    !-----------------------------------------------------------------------------
    ! Pre-scans a residue input file to determine the maximum sizes needed
    ! for array allocations
    !-----------------------------------------------------------------------------
    subroutine predetect_number_info(INFILE)

        ! Input parameters
        integer, intent(in) :: INFILE

        ! Local variables
        character(len=256) :: line, token, rest_line
        integer :: ios, ios_val
        integer :: pos, len_rest, val, pos_space
        integer :: nb_type_per_residue
        integer :: max_atom_in_residue
        logical :: in_residue_block
        logical :: is_active_residue

        ! Initialize counters
        nb%type_residue = 0
        nb%max_type_per_residue = 0
        nb%max_atom_in_residue = 0
        in_residue_block = .false.

        ! Loop over lines in input file
        do
            read(INFILE, '(A)', iostat=ios) line
            if (ios /= 0) exit  ! End of file or read error

            line = adjustl(trim(line))

            ! Detect start/end of residue block
            select case (trim(line))
            case ('begin_residue')
                in_residue_block = .true.
                cycle
            case ('end_residue')
                in_residue_block = .false.
                nb%type_residue = nb%type_residue + 1
                cycle
            end select

            if (.not. in_residue_block) cycle

            ! Tokenize line: first word and remainder
            pos_space = index(line, ' ')
            if (pos_space > 0) then
                token = line(1:pos_space-1)
                rest_line = adjustl(line(pos_space+1:))
            else
                token = line
                rest_line = ''
            end if

            ! Count atom types per residue
            if (trim(token) == 'types') then
                len_rest = len_trim(rest_line)
                pos = 1
                nb_type_per_residue = 0

                do
                    read(rest_line(pos:), *, iostat=ios_val) val
                    if (ios_val /= 0) exit
                    nb_type_per_residue = nb_type_per_residue + 1

                    ! Advance to next word
                    pos = pos + index(rest_line(pos:), ' ')
                    if (pos == 0 .or. pos > len_rest) exit
                end do

                ! Update global max
                if (nb_type_per_residue > nb%max_type_per_residue) then
                    nb%max_type_per_residue = nb_type_per_residue
                end if
            end if

            ! Determine max number of atoms per residue
            pos = index(line, ' ')
            if (pos > 0) then
                token = line(1:pos-1)
                rest_line = adjustl(line(pos+1:))
            else
                token = trim(line)
                rest_line = ''
            end if

            if (trim(token) == 'state') then
                if (trim(rest_line) == 'actif') then
                    is_active_residue = .true.
                else
                    is_active_residue = .false.
                end if
            end if

            if (trim(token) == 'nb-atoms') then
                read(rest_line, *, iostat=ios) max_atom_in_residue

                ! Update global maximum
                if (max_atom_in_residue > nb%max_atom_in_residue) then
                    nb%max_atom_in_residue = max_atom_in_residue
                end if

                ! Update active / inactive maxima
                if (is_active_residue) then
                    if (max_atom_in_residue > nb%max_atom_in_residue_active) then
                        nb%max_atom_in_residue_active = max_atom_in_residue
                    end if
                else
                    if (max_atom_in_residue > nb%max_atom_in_residue_inactive) then
                        nb%max_atom_in_residue_inactive = max_atom_in_residue
                    end if
                end if
            end if
        end do

    end subroutine predetect_number_info

    !-----------------------------------------------------------------------------
    ! Allocates all arrays for atoms, residues, molecules, and interaction parameters
    ! for both primary system and reservoir, based on system sizes and residue types.
    !-----------------------------------------------------------------------------
    subroutine allocate_atom_arrays()

        ! Allocate atom arrays for primary and reservoir

        call allocate_coordinate(host, .true.)
        call allocate_coordinate(guest, .false.)
        call allocate_coordinate(gas, .false.)
        allocate(resid_location(nb%type_residue))

        call allocate_atom_block(primary)
        call allocate_atom_block(reservoir)

        ! Allocate residue-level arrays
        allocate(res%mass(nb%type_residue))
        allocate(res%lambda(nb%type_residue))
        allocate(res%types_2d(nb%type_residue, nb%max_type_per_residue))
        allocate(res%names_2d(nb%type_residue, nb%max_type_per_residue))
        allocate(res%names_1d(nb%type_residue))
        allocate(res%bond_type_2d(nb%type_residue, NB_MAX_BOND, 3))
        allocate(res%angle_type_2d(nb%type_residue, NB_MAX_ANGLE, 4))
        allocate(res%dihedral_type_2d(nb%type_residue, NB_MAX_DIHEDRAL, 5))
        allocate(res%improper_type_2d(nb%type_residue, NB_MAX_IMPROPER, 5))

        ! Allocate interaction parameters
        allocate(coeff%sigma(nb%type_residue, nb%type_residue, &
                            nb%max_atom_in_residue, nb%max_atom_in_residue))
        allocate(coeff%epsilon(nb%type_residue, nb%type_residue, &
                            nb%max_atom_in_residue, nb%max_atom_in_residue))

        ! Allocate input arrays
        allocate(input%fugacity(nb%type_residue))
        allocate(input%chemical_potential(nb%type_residue))
        allocate(input%is_active(nb%type_residue))

        ! Allocate system bookkeeping arrays
        allocate(nb%types_per_residue(nb%type_residue))
        allocate(nb%atom_in_residue(nb%type_residue))
        allocate(nb%bonds_per_residue(nb%type_residue))
        allocate(nb%angles_per_residue(nb%type_residue))
        allocate(nb%dihedrals_per_residue(nb%type_residue))
        allocate(nb%impropers_per_residue(nb%type_residue))
        allocate(nb%types_pattern(nb%type_residue, nb%max_atom_in_residue))

    end subroutine allocate_atom_arrays

    !-----------------------------------------------------------------------------
    ! Allocates all atom-level arrays for a given system, including charges, masses,
    ! types, IDs, molecular coordinates, and residue counters.
    ! This subroutine is used for both primary system and reservoir.
    !-----------------------------------------------------------------------------
    subroutine allocate_atom_block(system)

        type(type_box), intent(inout) :: system

        ! Allocate basic atom properties
        allocate(system%atom_charges(nb%type_residue, nb%max_atom_in_residue))
        allocate(system%atom_masses(nb%type_residue, nb%max_atom_in_residue))
        allocate(system%atom_names(nb%type_residue, nb%max_atom_in_residue))
        allocate(system%atom_types(nb%type_residue, nb%max_atom_in_residue))
        allocate(system%atom_ids(nb%type_residue, nb%max_atom_in_residue))

        ! Allocate molecular coordinates
        allocate(system%mol_com(3, nb%type_residue, NB_MAX_MOLECULE))
        allocate(system%site_offset(3, nb%type_residue, NB_MAX_MOLECULE, nb%max_atom_in_residue))

        ! Allocate residue counters
        allocate(system%num_residues(nb%type_residue))

    end subroutine allocate_atom_block

    !-----------------------------------------------------------------------------
    ! Allocate coordinate array for host, guest, and gas
    !-----------------------------------------------------------------------------
    subroutine allocate_coordinate(type, is_host)

        ! Input parameters
        type(type_coordinate), intent(inout) :: type    ! host, guest, gas
        logical :: is_host                              ! To differentiate host

        ! Local variables
        integer :: dim = 3                              ! Number of physical dimensions
        integer :: max_atom                             ! Maximum number of atom in the residue
        integer :: max_molecule                         ! Maximum number of molecule

        ! Detect max number of atom (generaly much larger for host)
        if (is_host) then
            max_atom = nb%max_atom_in_residue_inactive
            max_molecule = 1 ! Only one host (may be too rigid?)
        else
            max_atom = nb%max_atom_in_residue_active
            max_molecule = NB_MAX_MOLECULE
        end if

        ! Allocate molecular coordinates
        allocate(type%residue_exists(nb%type_residue))
        allocate(type%com(dim, nb%type_residue, max_molecule))


        write (*,*) "dim, nb%type_residue, max_molecule, max_atom", dim, nb%type_residue, max_molecule, max_atom

        stop 778

        allocate(type%offset(dim, nb%type_residue, max_molecule, max_atom))
        host%max_nb_atom = max_atom
        host%max_nb_molecule = max_molecule
        type%residue_exists = .false.

    end subroutine allocate_coordinate

    !-----------------------------------------------------------------------------
    ! Reads the MANIAC input file, parses global parameters and residue blocks,
    ! validates required values, and initializes simulation probabilities and states.
    !-----------------------------------------------------------------------------
    subroutine ParseInputFile(INFILE)

        ! Input parameter
        integer, intent(in) :: INFILE

        ! Local variables
        character(len=256) :: line           ! Current line read from input file
        character(len=256) :: rest_line      ! Portion of line after the keyword
        character(len=256) :: keyword        ! Keyword parsed from input line
        character(len=256) :: token          ! Token parsed from line (e.g., parameter name)
        character(len=100) :: val_cha        ! String value read from input (longer strings)
        character(len=10)  :: val_str        ! Short string value from input (e.g., site names)
        integer :: ios                        ! I/O status for read operations
        integer :: val_int                     ! Integer value read from input
        integer :: pos                         ! Position index within a string
        integer :: pos_val                     ! Start position of value in line after keyword
        integer :: n_ids                       ! Counter for number of IDs parsed (e.g., types/names)
        integer :: len_rest                    ! Length of rest_line
        real(real64) :: val_real               ! Real value read from input
        real(real64) :: proba_total            ! Sum of move probabilities for normalization
        logical :: in_residue_block            ! Flag indicating parsing inside a residue block
        logical :: val_bool                     ! Logical value read from input
        integer :: val                         ! Temporary integer for parsing site/type IDs
        logical :: has_nb_block, has_nb_step, has_temp, has_seed
        logical :: has_cutoff, has_tolerance
        logical :: has_translation_step, has_rotation_step, has_recalibrate
        logical :: has_translation_proba, has_rotation_proba
        logical :: has_insertdel_proba, has_swap_proba, has_widom_proba
        logical :: has_fugacity, has_chemical_potential

        has_nb_block = .false.
        has_nb_step = .false.
        has_temp = .false.
        has_seed = .false.
        has_cutoff = .false.
        has_tolerance = .false.
        has_translation_step = .false.
        has_rotation_step = .false.
        has_recalibrate = .false.
        has_translation_proba = .false.
        has_rotation_proba = .false.
        has_insertdel_proba = .false.
        has_swap_proba = .false.
        has_widom_proba = .false.

        ! Initialize state
        in_residue_block = .false.
        nb%type_residue = 0

        ! Initialize all fugacities to -1.0 (indicating unset)
        input%fugacity(:) = -one
        ! Initialize all chemical_potential to 0.0 (indicating unset)
        input%chemical_potential(:) = zero

        do
            read(INFILE, '(A)', IOSTAT=ios) line
            if (ios /= 0) exit
            if (line(1:1) == '#' .or. len_trim(line) == 0) cycle

            read(line, *) keyword
            pos_val = index(line, keyword) + len_trim(keyword) + 1
            if (pos_val > len(line)) then
                rest_line = ''
            else
                rest_line = adjustl(line(pos_val:))
            end if

            select case (trim(adjustl(keyword)))

            case ("nb_block")
                read(rest_line, *, iostat=ios) val_int
                if (ios /= 0) error stop "Error reading nb_block"
                status%desired_block = val_int
                has_nb_block = .true.

            case ("nb_step")
                read(rest_line, *, iostat=ios) val_int
                if (ios /= 0) error stop "Error reading nb_step"
                status%desired_step = val_int
                has_nb_step = .true.

            case ("temperature")
                read(rest_line, *, iostat=ios) val_real
                if (ios /= 0) error stop "Error reading temperature"
                if (val_real <= zero) error stop "Invalid temperature: must be > 0"
                input%temperature = val_real
                has_temp = .true.

            case ("seed")
                read(rest_line, *, iostat=ios) val_int
                if (ios /= 0) error stop "Error reading seed"
                input%seed = val_int
                has_seed = .true.

            case ("ewald_tolerance")
                read(rest_line, *, iostat=ios) val_real
                if (ios /= 0) error stop "Error reading ewald_tolerance"
                if (val_real <= zero) error stop "Invalid ewald_tolerance: must be > 0"
                ewald%tolerance = val_real
                has_tolerance = .true.

            case ("real_space_cutoff")
                read(rest_line, *, iostat=ios) val_real
                if (ios /= 0) error stop "Error reading real_space_cutoff"
                if (val_real <= zero) error stop "Invalid real_space_cutoff: must be > 0"
                input%real_space_cutoff = val_real
                has_cutoff = .true.

            case ("translation_step")
                read(rest_line, *, iostat=ios) val_real
                if (ios /= 0) error stop "Error reading translation_step"
                if (val_real <= zero) error stop "Invalid translation_step: must be > 0"
                input%translation_step = val_real
                has_translation_step = .true.

            case ("rotation_step_angle")
                read(rest_line, *, iostat=ios) val_real
                if (ios /= 0) error stop "Error reading rotation_step_angle"
                if (val_real <= zero) error stop "Invalid rotation_step_angle: must be > 0"
                input%rotation_step_angle = val_real
                has_rotation_step = .true.

            case ("recalibrate_moves")
                read(rest_line, *, iostat=ios) val_bool
                if (ios /= 0) error stop "Error reading recalibrate_moves"
                input%recalibrate_moves = val_bool
                has_recalibrate = .true.

            case ("translation_proba")
                read(rest_line, *, iostat=ios) val_real
                if (ios /= 0) error stop "Error reading translation_proba"
                if (val_real < zero .or. val_real > one) error stop &
                    "Invalid translation_proba: must be in [0,1]"
                proba%translation = val_real
                has_translation_proba = .true.

            case ("rotation_proba")
                read(rest_line, *, iostat=ios) val_real
                if (ios /= 0) error stop "Error reading rotation_proba"
                if (val_real < zero .or. val_real > one) error stop &
                    "Invalid rotation_proba: must be in [0,1]"
                proba%rotation = val_real
                has_rotation_proba = .true.

            case ("insertion_deletion_proba")
                read(rest_line, *, iostat=ios) val_real
                if (ios /= 0) error stop "Error reading insertion_deletion_proba"
                if (val_real < zero .or. val_real > one) error stop &
                    "Invalid insertion_deletion_proba: must be in [0,1]"
                proba%insertion_deletion = val_real
                has_insertdel_proba = .true.

            case ("swap_proba")
                read(rest_line, *, iostat=ios) val_real
                if (ios /= 0) error stop "Error reading swap_proba"
                if (val_real < zero .or. val_real > one) error stop &
                    "Invalid swap_proba: must be in [0,1]"
                proba%swap = val_real
                has_swap_proba = .true.

            case ("widom_proba")
                read(rest_line, *, iostat=ios) val_real
                if (ios /= 0) error stop "Error reading widom_proba"
                if (val_real < zero .or. val_real > one) error stop &
                    "Invalid widom_proba: must be in [0,1]"
                proba%widom = val_real
                has_widom_proba = .true.

            case ("begin_residue")
                in_residue_block = .true.
                cycle

            case ("end_residue")
                in_residue_block = .false.
                nb%type_residue = nb%type_residue + 1
                cycle
            end select

            ! Residue parsing (types, names, fugacity, chemical potential nb-atoms)
            if (in_residue_block) then

                ! First, try to read the first two words from the line:
                ! 'token' will hold the keyword (e.g., 'name', 'state')
                ! 'val_cha' will hold the string value after the keyword
                read(line, *) token, val_cha

                ! Parse residue name or state
                if (trim(token) == "name") then

                    ! Store the residue name in the 1D array of residue names
                    res%names_1d(nb%type_residue + 1) = trim(val_cha)

                ! Store the activity state of the residue: active (1) or inactive (0)
                else if (trim(token) == "state") then
                    select case (trim(val_cha))
                        case ("actif")
                            input%is_active(nb%type_residue + 1) = 1
                        case ("inactif")
                            input%is_active(nb%type_residue + 1) = 0
                        case default
                            call WarnUser("Unknown state: " // trim(val_cha))
                            call AbortRun("Unknown residue state")
                    end select
                end if

                ! Parse fugacity and number of atoms
                ! Adjust line for further parsing
                line = adjustl(trim(line))

                ! Split line into token and the rest of the line
                pos = index(line, ' ')
                if (pos > 0) then
                    token = line(1:pos-1)
                    rest_line = adjustl(line(pos+1:))
                else
                    token = trim(line)
                    rest_line = ''
                end if

                ! Check if the line specifies the fugacity
                if (trim(token) == 'fugacity') then
                    read(rest_line, *, iostat=ios) val_real
                    input%fugacity(nb%type_residue + 1) = val_real

                ! Check if the line specifies the chemical potential
                else if (trim(token) == 'chemical_potential') then
                    read(rest_line, *, iostat=ios) val_real
                    input%chemical_potential(nb%type_residue + 1) = val_real

                ! Check if the line specifies the number of atoms
                else if (trim(token) == 'nb-atoms') then
                    read(rest_line, *, iostat=ios) val_int
                    nb%atom_in_residue(nb%type_residue + 1) = val_int
                end if
            end if

            ! Parse site types
            line = adjustl(trim(line))
            pos = index(line, ' ')
            if (pos > 0) then
                token = line(1:pos-1)
                rest_line = adjustl(line(pos+1:))
            else
                token = line
                rest_line = ''
            end if

            ! Parse the types of atoms for this residue
            if (trim(token) == 'types') then
                len_rest = len_trim(rest_line)
                pos = 1
                n_ids = 0
                do while (pos <= len_rest)
                    read(rest_line(pos:), *, iostat=ios) val
                    if (ios /= 0) exit
                    n_ids = n_ids + 1
                    res%types_2d(nb%type_residue + 1, n_ids) = val
                    pos = pos + index(rest_line(pos:), ' ')
                end do
                nb%types_per_residue(nb%type_residue + 1) = n_ids

            ! Parse the names of atom types for this residue
            else if (trim(token) == 'names') then
                len_rest = len_trim(rest_line)
                pos = 1
                n_ids = 0
                do while (pos <= len_rest)
                    read(rest_line(pos:), *, iostat=ios) val_str
                    if (ios /= 0) exit
                    n_ids = n_ids + 1
                    res%names_2d(nb%type_residue + 1, n_ids) = val_str
                    pos = pos + index(rest_line(pos:), ' ')
                end do
            end if
        end do

        ! Validate fugacity for active residues
        do val_int = 1, nb%type_residue

            if (input%is_active(val_int) == 0) cycle

            has_fugacity = (input%fugacity(val_int) >= zero)
            has_chemical_potential = (input%chemical_potential(val_int) < zero)

            ! Rule 1: at least one must be provided
            if (.not.(has_fugacity .or. has_chemical_potential)) then
                call AbortRun("Neither fugacity nor chemical potential provided for active residue: " // &
                            trim(res%names_1d(val_int)))
            end if

            ! Rule 2: cannot both be provided
            if (has_fugacity .and. has_chemical_potential) then
                call AbortRun("Both fugacity and chemical potential were specified for active residue: " // &
                            trim(res%names_1d(val_int)))
            end if

        end do

        ! === Validation of required parameters ===
        if (.not. has_nb_block)         error stop "Missing required parameter: nb_block"
        if (.not. has_nb_step)          error stop "Missing required parameter: nb_step"
        if (.not. has_temp)             error stop "Missing required parameter: temperature"
        if (.not. has_cutoff)           error stop "Missing required parameter: real_space_cutoff"
        if (.not. has_tolerance)        error stop "Missing required parameter: ewald_tolerance"
        if (.not. has_translation_step) error stop "Missing required parameter: translation_step"
        if (.not. has_rotation_step)    error stop "Missing required parameter: rotation_step_angle"

        ! Non mandatoray parameters
        if (.not. has_recalibrate) input%recalibrate_moves = .false.
        ! === Ensure only enabled moves are counted ===
        if (.not. has_insertdel_proba) proba%insertion_deletion = zero
        if (.not. has_rotation_proba) proba%rotation = zero
        if (.not. has_translation_proba) proba%translation = zero
        if (.not. has_swap_proba) proba%swap = zero
        if (.not. has_widom_proba) proba%widom = zero

        ! === Sum of enabled probabilities ===
        proba_total = proba%translation + &
            proba%rotation + &
            proba%insertion_deletion + &
            proba%swap + &
            proba%widom

        ! === Check for impossible case ===
        if (proba_total < error) then
            call AbortRun("Invalid move probabilities: all enabled moves have zero probability", 1)
        end if

        if (proba%widom > 0 .and. proba%insertion_deletion > 0) then
            call AbortRun("Cannot enable both Widom insertions and physical insertion/deletion moves", 1)
        end if

        ! Use provided seed, or generate a new one
        if (.not. has_seed) call seed_rng(input%seed)

        call SortResidues()

    end subroutine ParseInputFile

    subroutine SortResidues()

        integer :: i, j, k, n
        integer, allocatable :: keys(:), order(:)
        character(len=10), allocatable :: tmp_names_1d(:)
        integer, allocatable :: tmp_is_active(:), tmp_atom_in_residue(:), tmp_types_per_residue(:)
        real(real64), allocatable :: tmp_fugacity(:), tmp_chemical_potential(:)
        integer, allocatable :: tmp_types_2d(:,:)
        character(len=10), allocatable :: tmp_names_2d(:,:)

        ! Number of residues
        n = nb%type_residue

        ! Allocate temporary arrays
        allocate(keys(n), order(n))
        allocate(tmp_names_1d(n))
        allocate(tmp_is_active(n))
        allocate(tmp_fugacity(n))
        allocate(tmp_chemical_potential(n))
        allocate(tmp_atom_in_residue(n))
        allocate(tmp_types_per_residue(n))
        allocate(tmp_types_2d(size(res%types_2d,1), size(res%types_2d,2)))
        allocate(tmp_names_2d(size(res%names_2d,1), size(res%names_2d,2)))

        ! Compute sorting key for each residue: minimum atom type ID
        do i = 1, n
            keys(i) = minval(res%types_2d(i, 1:nb%types_per_residue(i)))
        end do

        ! Initialize the order array to [1,2,3,...,n]
        order = [(i, i=1,n)]

        ! Simple insertion sort based on keys (min atom type)
        do i = 2, n
            k = order(i)
            j = i - 1
            do while (j >= 1 .and. keys(order(j)) > keys(k))
                order(j+1) = order(j)
                j = j - 1
            end do
            order(j+1) = k
        end do

        ! Reorder all residue arrays based on the computed order
        do i = 1, n
            k = order(i)
            tmp_names_1d(i)          = res%names_1d(k)
            tmp_is_active(i)         = input%is_active(k)
            tmp_fugacity(i)          = input%fugacity(k)
            tmp_chemical_potential(i) = input%chemical_potential(k)
            tmp_atom_in_residue(i)   = nb%atom_in_residue(k)
            tmp_types_per_residue(i) = nb%types_per_residue(k)
            tmp_types_2d(i,:)        = res%types_2d(k,:)
            tmp_names_2d(i,:)        = res%names_2d(k,:)
            ! tmp_types_2d(i, 1:nb%types_per_residue(k)) = res%types_2d(k, 1:nb%types_per_residue(k))
            ! tmp_names_2d(i, 1:nb%types_per_residue(k)) = res%names_2d(k, 1:nb%types_per_residue(k))
        end do

        ! Copy temporary arrays back into original arrays
        res%names_1d          = tmp_names_1d
        input%is_active       = tmp_is_active
        input%fugacity        = tmp_fugacity
        input%chemical_potential = tmp_chemical_potential
        nb%atom_in_residue    = tmp_atom_in_residue
        nb%types_per_residue  = tmp_types_per_residue
        res%types_2d          = tmp_types_2d
        res%names_2d          = tmp_names_2d

        ! Deallocate temporary arrays
        deallocate(keys, order, tmp_names_1d, tmp_is_active, tmp_fugacity, tmp_chemical_potential, &
                tmp_atom_in_residue, tmp_types_per_residue, tmp_types_2d, tmp_names_2d)

    end subroutine SortResidues

    !-----------------------------------------------------------------------------    !
    ! Checks the I/O status returned by Fortran read/open/write operations
    ! and aborts the program with a message if an error occurred.
    !-----------------------------------------------------------------------------
    subroutine check_IO_status(filename, ios)

        ! Input parameters
        character(len=*), intent(in) :: filename
        integer, intent(in) :: ios

        if (ios /= 0) then
            call AbortRun("I/O error on file: "//trim(filename), ios)
        end if

    end subroutine check_IO_status

end module input_parser
