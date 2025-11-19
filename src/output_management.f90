module output_management

    use simulation_state
    use output_utils
    use helper_utils
    use version_module

    implicit none

contains

    !---------------------------------------------------------------------------
    ! Open the main output file for writing and create the output
    ! directory if necessary.
    !---------------------------------------------------------------------------
    subroutine setup_output_files()

        implicit none

        ! Ensure output directory exists
        call ensure_directory_exists(output_path)

        ! Open log file for writing, replacing any existing file
        call open_log_file(output_path, 'log.maniac')

        ! Write header to log
        call write_header()

    end subroutine setup_output_files

    !-------------------------------------------------------
    ! Ensure directory exists (create if missing)
    !-------------------------------------------------------
    subroutine ensure_directory_exists(path)

        implicit none

        ! Local variables
        character(len=*), intent(in) :: path
        character(len=200) :: command
        integer :: exit_status
        logical :: exists
        
        inquire(file=trim(path), exist=exists)

        if (.not. exists) then

            command = 'mkdir -p ' // trim(path)
        
            call execute_command_line(command, exitstat=exit_status)
        
            if (exit_status /= 0) call AbortRun("Failed to create output directory: "//trim(path), exit_status)
        end if

    end subroutine ensure_directory_exists

    !-------------------------------------------------------
    ! Open log file
    !-------------------------------------------------------
    subroutine open_log_file(path, filename, unit_number)

        implicit none

        ! Input parameters
        character(len=*), intent(in) :: path               ! Directory path for log file
        character(len=*), intent(in), optional :: filename ! Optional log file name
        integer, intent(in), optional :: unit_number       ! Optional Fortran unit number for log

        ! Local variables
        character(len=200) :: logname  ! Full log file name
        integer :: ios                 ! I/O status from OPEN
        integer :: unit_to_use         ! Actual unit used for opening the log

        if (present(unit_number)) then
            unit_to_use = unit_number
        else
            unit_to_use = out_unit
        end if

        if (present(filename)) then
            logname = filename
        else
            logname = 'log.maniac'
        end if

        open(unit=unit_to_use, file=join_path(path, logname), status='replace', iostat=ios)

        if (ios /= 0) then
            call AbortRun("Failed to open log file: "//trim(path)//trim(logname), ios)
        end if
    
    end subroutine open_log_file

    !---------------------------------------------------------------------------
    ! Write a modern ASCII header with project info (dynamic version)
    !---------------------------------------------------------------------------
    subroutine write_header()

        implicit none

        integer, parameter :: box_width = 78 ! width of the box

        ! Blank line before header
        call LogMessage("")

        ! Top border
        call LogMessage("+" // repeat_char("-", box_width-2) // "+")

        ! Version info lines
        call BoxLine("MANIAC-MC - Version " // version, box_width)

        ! Credits (optional, you can also keep them in version_module if desired)
        call BoxLine("Code written and maintained by Simon Gravelle, LIPhy, CNRS", box_width)

        ! Bottom border
        call LogMessage("+" // repeat_char("-", box_width-2) // "+")

        ! Blank line after header
        call LogMessage("")

    end subroutine write_header

end module output_management
