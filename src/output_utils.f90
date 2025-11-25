!-----------------------------------------------------------------------------
! This module handles all output-related operations for the Monte Carlo
! simulation.
!-----------------------------------------------------------------------------

module output_utils

    use simulation_state
    use helper_utils
    use, intrinsic :: iso_fortran_env, only: real64, output_unit

    implicit none

    integer :: out_unit = 10                        ! Default log file unit

contains

    !---------------------------------------------------------------------------
    ! Close the main output file.
    !---------------------------------------------------------------------------
    subroutine CloseOutput()

        call PrintTerminationMessage()

        close(out_unit)

    end subroutine CloseOutput

    !---------------------------------------------------------------------------
    ! Write a custom message to the specified output unit.
    !---------------------------------------------------------------------------
    subroutine log_message(msg)

        ! Input parameter
        character(*), intent(in) :: msg

        write(*,*) trim(msg)
        write(out_unit,*) trim(msg)
        flush(out_unit) ! Forces the data to be written to disk

    end subroutine log_message

    !---------------------------------------------------------------------------
    !  Log the start of the Monte Carlo loop with a modern ASCII style
    !---------------------------------------------------------------------------
    subroutine LogStartMC()

        ! Blank line before message
        call log_message("")

        ! Top border
        call log_message("+" // repeat_char("-", BOX_WIDTH-2) // "+")

        ! Message line centered
        call BoxLine("Started Monte Carlo Loop", BOX_WIDTH)

        ! Bottom border
        call log_message("+" // repeat_char("-", BOX_WIDTH-2) // "+")

        ! Blank line after message
        call log_message("")

    end subroutine LogStartMC

    !-----------------------------------------------------------------------
    ! Helper: BoxLine
    ! Writes a line inside the ASCII box, padded to box_width
    !-----------------------------------------------------------------------
    subroutine BoxLine(text, width)

        ! Input variables
        character(len=*), intent(in) :: text
        integer, intent(in) :: width

        ! Local variables
        character(len=width-4) :: padded

        padded = adjustl(text)
        if (len_trim(padded) < width-4) then
            padded(len_trim(padded)+1:) = ' '  ! pad with spaces
        end if

        call log_message("| " // padded(1:width-4) // " |")

    end subroutine BoxLine

    !---------------------------------------------------------------------------
    ! Print a formatted termination footer with Monte Carlo summary
    !---------------------------------------------------------------------------
    subroutine PrintTerminationMessage()

        integer :: type_residue
        character(len=256) :: line

        ! Blank line before footer
        call log_message("")

        ! Top border
        call log_message("+" // repeat_char("-", BOX_WIDTH-2) // "+")

        ! Title
        call BoxLine("MANIAC-MC simulation completed", BOX_WIDTH)
        call BoxLine("", BOX_WIDTH)  ! blank line inside box

        ! Output path information
        call BoxLine("All output files have been written to:", BOX_WIDTH)
        call BoxLine(trim(path%outputs), BOX_WIDTH)

        ! Summary statistics (Trial / Accepted moves)
        if (proba%translation > 0) then
            write(line,'(A,I8,A,I8)') "  Translations (Trial/Accepted): ", &
                counter%translations(1), " / ", counter%translations(2)
            call BoxLine(trim(line), BOX_WIDTH)
        end if

        if (proba%rotation > 0) then
            write(line,'(A,I8,A,I8)') "  Rotations    (Trial/Accepted): ", &
                counter%rotations(1), " / ", counter%rotations(2)
            call BoxLine(trim(line), BOX_WIDTH)
        end if

        if (proba%insertion_deletion > 0) then
            write(line,'(A,I8,A,I8)') "  Creations    (Trial/Accepted): ", &
                counter%creations(1), " / ", counter%creations(2)
            call BoxLine(trim(line), BOX_WIDTH)
            write(line,'(A,I8,A,I8)') "  Deletions    (Trial/Accepted): ", &
                counter%deletions(1), " / ", counter%deletions(2)
            call BoxLine(trim(line), BOX_WIDTH)
        end if

        if (proba%swap > 0) then
            write(line,'(A,I8,A,I8)') "  Swap         (Trial/Accepted): ", &
                counter%swaps(1), " / ", counter%swaps(2)
            call BoxLine(trim(line), BOX_WIDTH)
        end if

        if (proba%widom > 0) then
            do type_residue = 1, nb%type_residue
                if (statistic%sample(type_residue) > 0) then
                    write(line,'(A,I3,A,I8)') "  Widom trials for residue ", type_residue, ": ", statistic%sample(type_residue)
                    call BoxLine(trim(line), BOX_WIDTH)
                    write(line,'(A,F12.5)') "    Excess chemical potential (kcal/mol): ", statistic%mu_ex(type_residue)
                    call BoxLine(trim(line), BOX_WIDTH)
                    write(line,'(A,F12.5)') "    Total chemical potential (kcal/mol): ", statistic%mu_tot(type_residue)
                    call BoxLine(trim(line), BOX_WIDTH)
                end if
            end do
        end if

        call BoxLine("", BOX_WIDTH)  ! blank line inside box

        ! Bottom border
        call log_message("+" // repeat_char("-", BOX_WIDTH-2) // "+")

        ! Blank line after footer
        call log_message("")

    end subroutine PrintTerminationMessage

    subroutine PrintStatus()

        ! Local variables
        integer :: nb_type_residue
        character(len=64)  :: tmp
        character(len=1024) :: header_msg
        character(len=1024) :: move_msg
        character(len=256) :: numeric_msg

        real(real64) :: e_tot, e_coul, e_long

        ! Blank line before status
        call log_message("")

        ! -----------------------
        ! Active molecule summary
        ! -----------------------
        header_msg = "  Energy report | Active molecules: "
        do nb_type_residue = 1, nb%type_residue
            if (thermo%is_active(nb_type_residue)) then
                write(tmp,'(A,"=",I0)') trim(res%names_1d(nb_type_residue)), primary%num%residues(nb_type_residue)
                if (len_trim(header_msg) > 0) then
                    header_msg = trim(header_msg)//" "//trim(tmp)
                else
                    header_msg = trim(tmp)
                end if
            end if
        end do
        call log_message(header_msg)

        ! Composite energies
        call compute_composite_energies(e_tot, e_coul, e_long)

        ! -----------------------
        ! Print header for energy and MC moves
        ! -----------------------
        write(header_msg,'(A10,1X,A14,1X,A14,1X,A14,1X,A14,2X,A10,2X,A10,2X,A20)') &
            'Step','TotEng','E_vdwl','E_coul','E_long','TransStep','RotAngle','MC (acc/trial)'
        call log_message(header_msg)

        ! -----------------------
        ! Build dynamic MC move statistics into move_msg
        ! -----------------------
        move_msg = ""

        if (proba%translation > 0) then
            write(tmp,'("T(",I0,"/",I0,")")') counter%translations(2), counter%translations(1)
            move_msg = trim(move_msg)//" "//trim(tmp)
        end if

        if (proba%rotation > 0) then
            write(tmp,'("R(",I0,"/",I0,")")') counter%rotations(2), counter%rotations(1)
            move_msg = trim(move_msg)//" "//trim(tmp)
        end if

        if (proba%insertion_deletion > 0) then
            write(tmp,'("C(",I0,"/",I0,")")') counter%creations(2), counter%creations(1)
            move_msg = trim(move_msg)//" "//trim(tmp)
            write(tmp,'("D(",I0,"/",I0,")")') counter%deletions(2), counter%deletions(1)
            move_msg = trim(move_msg)//" "//trim(tmp)
        end if

        if (proba%swap > 0) then
            write(tmp,'("S(",I0,"/",I0,")")') counter%swaps(2), counter%swaps(1)
            move_msg = trim(move_msg)//" "//trim(tmp)
        end if

        if (proba%widom > 0) then
            write(tmp,'("S(",I0,"/",I0,")")') counter%widom(2), counter%widom(1)
            move_msg = trim(move_msg)//" "//trim(tmp)
        end if

        ! -----------------------
        ! Build the complete single-line status (numbers + MC moves)
        ! -----------------------
        write(numeric_msg,'(I10,1X,F14.4,1X,F14.4,1X,F14.4,1X,F14.4,2X,F10.4,2X,F10.4)') &
            status%block, e_tot, energy%non_coulomb, e_coul, e_long, &
            input%translation_step, input%rotation_step_angle

        ! Append MC move string to the same line
        if (len_trim(move_msg) > 0) then
            numeric_msg = trim(numeric_msg) // "  " // trim(move_msg)
        endif

        ! Print final combined line
        call log_message(numeric_msg)

    end subroutine PrintStatus

    !----------------------------------------------------------------------
    ! Computes composite energies from the energy components in 'energy'.
    !----------------------------------------------------------------------
    subroutine compute_composite_energies(e_tot, e_coul, e_long)

        ! Output parameter
        real(real64), intent(out) :: e_tot, e_coul, e_long

        ! Compute composite energies
        e_tot  = energy%non_coulomb + energy%recip_coulomb + energy%coulomb + &
                energy%ewald_self + energy%intra_coulomb
        e_coul = energy%coulomb + energy%intra_coulomb
        e_long = energy%recip_coulomb + energy%ewald_self

    end subroutine compute_composite_energies

    !------------------------------------------------------------------------------
    ! Subroutine: final_report
    !------------------------------------------------------------------------------
    subroutine final_report()

        ! Energy components in kcal/mol
        real(real64) :: e_tot          ! Total energy for reporting (computed)
        real(real64) :: e_coul         ! Coulombic energy including intra-molecular interactions
        real(real64) :: e_long         ! Long-range Coulombic energy (reciprocal + self)
        character(LEN=1024) :: formatted_msg   ! Formatted message for logging

        ! Compute combined components
        call compute_composite_energies(e_tot, e_coul, e_long)

        ! Blank line before box
        call log_message("")

        ! Top border
        call log_message("+" // repeat_char("-", BOX_WIDTH-2) // "+")

        ! Box title
        call BoxLine("Final Energy Report", BOX_WIDTH)
        call BoxLine("", BOX_WIDTH)

        ! Column headers
        call BoxLine("  Step        TotEng        E_vdwl        E_coul        E_long", BOX_WIDTH)

        ! Energies line
        write(formatted_msg,'(I10,1X,F15.6,1X,F15.6,1X,F15.6,1X,F15.6)') &
            status%block, e_tot, energy%non_coulomb, e_coul, e_long
        call BoxLine(trim(formatted_msg), BOX_WIDTH)

        call BoxLine("", BOX_WIDTH)  ! blank line inside box

        ! Bottom border
        call log_message("+" // repeat_char("-", BOX_WIDTH-2) // "+")

        ! Blank line after box
        call log_message("")

        call CloseOutput()                  ! Close files and finalize

    end subroutine final_report

    subroutine log_parameters(input_file_name, n_pairs, pair1, pair2, epsilons, sigmas)

        character(len=*), intent(in) :: input_file_name         ! Name of the parameter input file being processed
        integer, intent(out) :: n_pairs                         ! Counter for the number of unique atom type pairs
        integer, intent(out) :: pair1(:), pair2(:)              ! Arrays to store the first and second atom type indices of unique pairs
        real(real64), intent(out) :: epsilons(:), sigmas(:)     ! Arrays to store epsilon (kcal/mol) and sigma (Å) for unique pairs

        integer :: i, j, k, l, m, n                             ! Loop indices for iterating over residues and atom types
        integer :: type1, type2, type_m, type_n                 ! Atom type indices for site and atom pairs
        integer :: it1, it2                                     ! Sorted atom type indices (min and max) for unique pair identification
        logical :: found_pair                                   ! Flag to indicate if a pair was found (unused but set for clarity)
        character(len=200) :: formatted_msg                     ! String to hold formatted log messages for output
        logical, allocatable :: printed(:,:)                    ! 2D array to track which atom type pairs have been recorded
        integer :: max_atom_type                                ! Maximum atom type index from atom_types_2d array

        if (status%reservoir_provided) then
            max_atom_type = max(maxval(primary%atoms%types(:,:)), maxval(reservoir%atoms%types(:,:)) )
        else
            max_atom_type = maxval(primary%atoms%types(:,:))
        end if
        allocate(printed(max_atom_type, max_atom_type))
        printed = .false.

        ! Log the start of parameter file import
        call log_message("")
        call log_message("====== Import parameter file ======")
        call log_message("")
        write(formatted_msg, '("Reading file ", A)') trim(input_file_name) ! Format message with input file name
        call log_message(formatted_msg)                                     ! Log the input file name

        n_pairs = 0 ! Initialize the pair counter
        do i = 1, nb%type_residue
            do j = 1, nb%types_per_residue(i)
                do k = 1, nb%type_residue
                    do l = 1, nb%types_per_residue(k)
                        type1 = res%types_2d(i, j) ! Get atom type for site j in residue i
                        type2 = res%types_2d(k, l) ! Get atom type for site l in residue k
                        found_pair = .false.        ! Reset found_pair flag (unused in logic)
                        if (type1 <= type2) then    ! Process pairs where type1 <= type2 to avoid duplicates
                            do m = 1, nb%atom_in_residue(i)
                                do n = 1, nb%atom_in_residue(k)
                                    type_m = primary%atoms%types(i, m) ! Get atom type for atom m in residue i
                                    type_n = primary%atoms%types(k, n) ! Get atom type for atom n in residue k
                                    it1 = min(type_m, type_n)    ! Get minimum atom type for pair
                                    it2 = max(type_m, type_n)    ! Get maximum atom type for pair
                                    if (it1 > 0 .and. it2 > 0) then
                                        if (.not. printed(it1, it2)) then ! Check if pair hasn’t been recorded
                                            n_pairs = n_pairs + 1         ! Increment pair counter
                                            pair1(n_pairs) = it1          ! Store first atom type
                                            pair2(n_pairs) = it2          ! Store second atom type
                                            epsilons(n_pairs) = coeff%epsilon(i,k,m,n)              ! Store epsilon in kcal/mol
                                            sigmas(n_pairs) = coeff%sigma(i,k,m,n)                  ! Store sigma in Å
                                            printed(it1, it2) = .true.    ! Mark pair as recorded
                                            printed(it2, it1) = .true.    ! Mark symmetric pair as recorded
                                        end if
                                    end if
                                end do
                            end do
                        end if
                    end do
                end do
            end do
        end do
    end subroutine log_parameters

    subroutine LogData(data_file_name, box, is_primary)

        ! Input parameters
        type(type_box), intent(inout) :: box
        logical, intent(in) :: is_primary

        ! Local variables
        character(len=*), intent(in) :: data_file_name
        character(len=200) :: formatted_msg
        integer :: i, atom_id
        integer :: active_molecule_count ! Number of species that will be moved/inserted/deleted

        ! === Step 9: Logging ===
        call log_message("")
        call log_message("====== Import data file ======")
        write(formatted_msg, '("Reading file ", A)') trim(data_file_name)
        call log_message(formatted_msg)
        call log_message("")
        write(formatted_msg, '("Number of atoms: ", I0)') box%num%atoms
        call log_message(formatted_msg)
        write(formatted_msg, '("Number of type of residues: ", I0)') nb%type_residue
        call log_message(formatted_msg)
        write(formatted_msg, '("Number of type of atoms: ", I0)') box%num%atomtypes
        call log_message(formatted_msg)

        do i = 1, nb%type_residue
            if ((box%num%residues(i) /= 0) .and. (thermo%is_active(i))) then
                ! Active residue present in data file
                active_molecule_count = active_molecule_count + box%num%residues(i)
                write(formatted_msg, '("Active residue ", A, " found in the data file: ", I0)') &
                    trim(res%names_1d(i)), box%num%residues(i)
            else if ((box%num%residues(i) /= 0) .and. (thermo%is_active(i))) then
                ! Inactive residue present in data file
                active_molecule_count = active_molecule_count + box%num%residues(i)
                write(formatted_msg, '("Inactive residue ", A, " found in the data file: ", I0)') &
                    trim(res%names_1d(i)), box%num%residues(i)
            else if ((box%num%residues(i) == 0) .and. (thermo%is_active(i)) .and. (is_primary)) then
                ! Inactive residue defined in input but not present in data file
                call abort_run("Inactive residue '" // trim(res%names_1d(i)) // "' (ID=" // &
                            trim(adjustl(to_string(i))) // ") defined in input file but not present in data file.", 1)
            end if
            call log_message(formatted_msg)
        end do

        call log_message("")
        call log_message("Simulation box (rows):")
        write(formatted_msg, '(3F12.6)') box%cell%matrix(1,1), box%cell%matrix(1,2), box%cell%matrix(1,3)
        call log_message(formatted_msg)
        write(formatted_msg, '(3F12.6)') box%cell%matrix(2,1), box%cell%matrix(2,2), box%cell%matrix(2,3)
        call log_message(formatted_msg)
        write(formatted_msg, '(3F12.6)') box%cell%matrix(3,1), box%cell%matrix(3,2), box%cell%matrix(3,3)
        call log_message(formatted_msg)

        call log_message("")
        call log_message("Atoms masses (g/mol):")
        do atom_id = 1, primary%num%atomtypes
            write(formatted_msg, '(I5, 2X, F12.6)') atom_id, box%atoms%masses_vec(atom_id)
            call log_message(formatted_msg)
        end do

    end subroutine LogData

    subroutine LogConnectivity(box)

        ! Input
        type(type_box), intent(inout) :: box

        ! Local
        character(len=200) :: formatted_msg
        integer :: i, j
        integer, parameter :: MAX_PRINT = 6

        if ((box%num%bonds > 0) .or. (box%num%angles > 0)) then

            call log_message("")
            call log_message("===== Connectivity summary =====")

            ! --- Bonds ---
            call log_message("")
            do i = 1, nb%type_residue
                if (box%num%residues(i) > 0) then
                    write(formatted_msg, '("Residue ", A, ": ", I0, " bonds")') &
                        trim(res%names_1d(i)), nb%bonds_per_residue(i)
                    call log_message(formatted_msg)

                    ! Print up to MAX_PRINT bonds
                    do j = 1, min(nb%bonds_per_residue(i), MAX_PRINT)
                        write(formatted_msg, '("   bond type ", I0, ": atoms [", I0, ",", I0, "]")') &
                            res%bond_type_2d(i,j,1), res%bond_type_2d(i,j,2), res%bond_type_2d(i,j,3)
                        call log_message(formatted_msg)
                    end do

                    ! If more bonds exist, indicate truncation
                    if (nb%bonds_per_residue(i) > MAX_PRINT) then
                        write(formatted_msg, '("   ... ", I0, " more bonds not shown")') &
                            nb%bonds_per_residue(i) - MAX_PRINT
                        call log_message(formatted_msg)
                    end if
                end if
            end do

            ! --- Angles ---
            call log_message("")
            do i = 1, nb%type_residue

                if (box%num%residues(i) > 0) then
                    write(formatted_msg, '("Residue ", A, ": ", I0, " angles")') &
                        trim(res%names_1d(i)), nb%angles_per_residue(i)
                    call log_message(formatted_msg)

                    ! Print up to MAX_PRINT angles
                    do j = 1, min(nb%angles_per_residue(i), MAX_PRINT)
                        write(formatted_msg, '("   angle type ", I0, ": atoms [", I0, ",", I0, ",", I0, "]")') &
                            res%angle_type_2d(i,j,1), res%angle_type_2d(i,j,2), &
                            res%angle_type_2d(i,j,3), res%angle_type_2d(i,j,4)
                        call log_message(formatted_msg)
                    end do

                    ! If more angles exist, indicate truncation
                    if (nb%angles_per_residue(i) > MAX_PRINT) then
                        write(formatted_msg, '("   ... ", I0, " more angles not shown")') &
                            nb%angles_per_residue(i) - MAX_PRINT
                        call log_message(formatted_msg)
                    end if
                end if
            end do

            ! --- Dihedrals ---
            call log_message("")
            do i = 1, nb%type_residue

                if (box%num%residues(i) > 0) then
                    write(formatted_msg, '("Residue ", A, ": ", I0, " dihedrals")') &
                        trim(res%names_1d(i)), nb%dihedrals_per_residue(i)
                    call log_message(formatted_msg)

                    ! Print up to MAX_PRINT dihedrals
                    do j = 1, min(nb%dihedrals_per_residue(i), MAX_PRINT)
                        write(formatted_msg, &
                        '("   dihedral type ", I0, ": atoms [", I0, ",", I0, ",", I0, ",", I0, "]")') &
                            res%dihedral_type_2d(i,j,1), res%dihedral_type_2d(i,j,2), &
                            res%dihedral_type_2d(i,j,3), res%dihedral_type_2d(i,j,4), &
                            res%dihedral_type_2d(i,j,5)
                        call log_message(formatted_msg)
                    end do

                    ! If more dihedrals exist, indicate truncation
                    if (nb%dihedrals_per_residue(i) > MAX_PRINT) then
                        write(formatted_msg, '("   ... ", I0, " more dihedrals not shown")') &
                            nb%dihedrals_per_residue(i) - MAX_PRINT
                        call log_message(formatted_msg)
                    end if
                end if
            end do

            ! --- Impropers ---
            call log_message("")
            do i = 1, nb%type_residue

                if (box%num%residues(i) > 0) then
                    write(formatted_msg, '("Residue ", A, ": ", I0, " impropers")') &
                        trim(res%names_1d(i)), nb%impropers_per_residue(i)
                    call log_message(formatted_msg)

                    ! Print up to MAX_PRINT impropers
                    do j = 1, min(nb%impropers_per_residue(i), MAX_PRINT)
                        write(formatted_msg, &
                        '("   improper type ", I0, ": atoms [", I0, ",", I0, ",", I0, ",", I0, "]")') &
                            res%improper_type_2d(i,j,1), res%improper_type_2d(i,j,2), &
                            res%improper_type_2d(i,j,3), res%improper_type_2d(i,j,4), &
                            res%improper_type_2d(i,j,5)
                        call log_message(formatted_msg)
                    end do

                    ! If more impropers exist, indicate truncation
                    if (nb%impropers_per_residue(i) > MAX_PRINT) then
                        write(formatted_msg, '("   ... ", I0, " more impropers not shown")') &
                            nb%impropers_per_residue(i) - MAX_PRINT
                        call log_message(formatted_msg)
                    end if
                end if
            end do

        end if

    end subroutine LogConnectivity

    subroutine abort_run(error_msg, exit_code)

        ! Input parameters
        character(len=*), intent(in) :: error_msg   ! Error description
        integer, intent(in), optional :: exit_code  ! Exit code (default = 1)

        ! Local variable
        integer :: code

        if (present(exit_code)) then
            code = exit_code
        else
            code = 1
        end if

        call log_message("--------------------------------------------------")
        call log_message("FATAL ERROR:")
        call log_message("" // trim(error_msg))
        call log_message("Simulation will now terminate.")
        call log_message("--------------------------------------------------")

        ! Ensure buffers/files are flushed if log_message writes to files
        call flush(output_unit)

        stop code

    end subroutine abort_run

    !-----------------------------------------------------------------------------
    ! Print a standardized warning message without terminating
    !-----------------------------------------------------------------------------
    subroutine warn_user(warn_msg)

        character(len=*), intent(in) :: warn_msg   ! Warning description

        call log_message("--------------------------------------------------")
        call log_message("WARNING:")
        call log_message("" // trim(warn_msg))
        call log_message("Execution will continue.")
        call log_message("--------------------------------------------------")

        ! Ensure message is flushed immediately
        call flush(output_unit)

    end subroutine warn_user

    !-----------------------------------------------------------------------------
    ! Print a standardized informational message
    !-----------------------------------------------------------------------------
    subroutine InfoMessage(info_msg)

        character(len=*), intent(in) :: info_msg   ! Information description

        call log_message("INFO: " // trim(info_msg))
        call flush(output_unit)

    end subroutine InfoMessage

    !-----------------------------------------------------------------------------
    ! Validates that the given molecule index does not exceed the maximum allowed
    ! number of molecules (NB_MAX_MOLECULE). Aborts the run with an informative
    ! message if the index is out of bounds.
    !-----------------------------------------------------------------------------
    subroutine check_molecule_index(molecule_index)

        ! Input argument
        integer, intent(in) :: molecule_index       ! Index of molecule being created

        ! Local variable
        character(len=200) :: formatted_msg         ! Buffer for constructing formatted log or error messages

        if (molecule_index > NB_MAX_MOLECULE) then
            write(formatted_msg, '(A, I0)') &
                "Trying to insert a molecule with index = " , molecule_index
            call abort_run(trim(formatted_msg)//new_line('a')// &
                        "This exceeds the maximum allowed number of molecules = "// &
                        trim(adjustl(to_string(NB_MAX_MOLECULE)))//new_line('a')// &
                        "Increase 'NB_MAX_MOLECULE' in src/parameters.f90 and recompile.")
        end if

    end subroutine check_molecule_index

    !-----------------------------------------------------------------------------
    ! Converts an integer input to a trimmed character string.
    ! Returns an allocatable string of minimal length.
    !-----------------------------------------------------------------------------
    pure function to_string(i) result(str)

        integer, intent(in) :: i
        character(len=:), allocatable :: str
        character(len=32) :: tmp                ! big enough for most integers

        write(tmp, '(I0)') i
        str = trim(tmp)

    end function to_string

    !-----------------------------------------------------------------------------
    ! Logs a summary of the parsed input, including global simulation parameters,
    ! Monte Carlo settings, and detailed residue information.
    !-----------------------------------------------------------------------------
    subroutine print_input_summary()

        ! Locals
        integer :: i, j
        character(len=200) :: msg, temp

        ! ===== HEADER =====
        call log_message("====== Import input file ======")
        call log_message("")
        write(msg, '("Reading file ", A)') trim(path%input)
        call log_message(msg)
        call log_message("")

        ! ===== GENERIC PARAMETERS =====
        call log_message("=== Generic parameters")
        write(msg, '("Number of blocks: ", I0)') status%desired_block
        call log_message(msg)
        write(msg, '("Number of steps: ", I0)') status%desired_step
        call log_message(msg)
        write(msg, '("Temperature (K): ", F10.2)') thermo%temperature
        call log_message(msg)
        call log_message("")

        ! ===== ELECTROSTATIC =====
        call log_message("=== Electrostatic interactions")
        write(msg, '("Ewald tolerance: ", F15.8)') ewald%tolerance
        call log_message(msg)
        write(msg, '("Cutoff (Å): ", F10.2)') input%real_space_cutoff
        call log_message(msg)
        call log_message("")

        ! ===== MONTE CARLO =====
        call log_message("=== Monte carlo move")
        write(msg, '("Translation step (Å): ", F10.2)') input%translation_step
        call log_message(msg)
        write(msg, '("Rotation step angle (radian): ", F10.2)') input%rotation_step_angle
        call log_message(msg)

        write(msg, '("Translation proba: ", F10.2)') proba%translation
        call log_message(msg)
        write(msg, '("Rotation proba: ", F10.2)') proba%rotation
        call log_message(msg)
        write(msg, '("Insertion deletion proba: ", F10.2)') proba%insertion_deletion
        call log_message(msg)
        write(msg, '("Swap proba: ", F10.2)') proba%swap
        call log_message(msg)
        call log_message("")

        ! ===== RESIDUES =====
        call log_message("=== Residue information")
        call log_message("")
        write(msg, '("Number of type of residue found: ", I0)') nb%type_residue
        call log_message(msg)
        call log_message("")

        do i = 1, nb%type_residue

            write(msg, '("  Residue ", A)') trim(res%names_1d(i))
            call log_message(msg)

            write(msg, '("  Is active: ", A)') merge("yes", "no ", thermo%is_active(i))
            call log_message(msg)

            if (thermo%is_active(i)) then
                if (thermo%fugacity(i) > 0) then
                    write(msg, '("  Fugacity (atm): ", F10.2)') thermo%fugacity(i)
                    call log_message(msg)
                end if
            end if

            if (thermo%is_active(i)) then
                write(msg, '("  Chemical potential (kcal/mol): ", F10.2)') thermo%chemical_potential(i)
                call log_message(msg)
            end if

            write(msg, '("  Number of atoms in residue: ", I0)') nb%atom_in_residue(i)
            call log_message(msg)

            write(msg, '("  Number of atom types in residue: ", I0)') nb%types_per_residue(i)
            call log_message(msg)

            ! Types
            msg = "  Types: "
            do j = 1, nb%types_per_residue(i)
                write(temp, '(I0)') res%types_2d(i, j)
                msg = trim(msg) // " " // trim(temp)
            end do
            call log_message(msg)

            ! Names
            msg = "  Names: "
            do j = 1, nb%types_per_residue(i)
                temp = res%names_2d(i, j)
                msg = trim(msg) // " " // trim(temp)
            end do
            call log_message(msg)

            call log_message("")

        end do

    end subroutine print_input_summary

    !-----------------------------------------------------------
    ! Print Ewald information
    !-----------------------------------------------------------
    subroutine log_ewald_parameters()

        character(200) :: formatted_msg

        write(formatted_msg, '(A, F10.4)') 'Real-space cutoff (Å): ', input%real_space_cutoff
        call log_message(formatted_msg)
        write(formatted_msg, '(A, ES12.5)') 'Ewald accuracy tolerance: ', ewald%tolerance
        call log_message(formatted_msg)
        write(formatted_msg, '(A, F10.4)') 'Screening factor (dimensionless): ', ewald%screening_factor
        call log_message(formatted_msg)
        write(formatted_msg, '(A, F10.4)') 'Ewald damping parameter alpha (1/Å): ', ewald%alpha
        call log_message(formatted_msg)
        write(formatted_msg, '(A, F10.4)') 'Fourier-space precision parameter: ', ewald%fourier_precision
        call log_message(formatted_msg)
        write(formatted_msg, '(A, I5, A, I5, A, I5)') 'Max Fourier index (kmax(1), kmax(2), kmax(3)): ', &
            ewald%kmax(1), ', ', ewald%kmax(2), ', ', ewald%kmax(3)
        call log_message(formatted_msg)
        write(formatted_msg, '(A, I10)') 'Total reciprocal lattice vectors: ', ewald%num_kvectors
        call log_message(formatted_msg)

    end subroutine log_ewald_parameters

end module output_utils
