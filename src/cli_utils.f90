module cli_utils

    use output_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !-----------------------------------------------------------------------------
    ! Reads and processes the command-line arguments
    !-----------------------------------------------------------------------------
    subroutine parse_command_line_arguments()

        ! Set default CLI values
        call set_default_CLI_values()

        ! Read command-line arguments into global variables
        call read_command_line_args()

        ! Validate CLI arguments and check file existence
        call validate_CLI_arguments()

        ! Ensure output path ends with a trailing slash
        call normalize_output_path()

    end subroutine parse_command_line_arguments

    !-----------------------------------------------------------------
    ! Set default CLI values
    !-----------------------------------------------------------------
    subroutine set_default_CLI_values()

        path%input = ''
        path%topology = ''
        path%parameters = ''
        path%reservoir = ''
        path%outputs = 'outputs/'

    end subroutine set_default_CLI_values

    !-----------------------------------------------------------------
    ! Read command-line arguments into global variables
    !-----------------------------------------------------------------
    subroutine read_command_line_args()
    
        ! Local variables
        character(len=256) :: arg
        integer :: i, nargs
        logical :: seen_maniac_file = .false.
        logical :: seen_data_file   = .false.
        logical :: seen_inc_file    = .false.
        logical :: seen_res_file    = .false.
        logical :: seen_output_path = .false.

        nargs = command_argument_count()

        do i = 1, nargs

            call get_command_argument(i, arg)
            
            select case(trim(arg))

            case ('-i')
                if (seen_maniac_file) then
                    call abort_run("Duplicate option: -i", 1)
                end if
                call expect_value(i, nargs, "-i", path%input)
                seen_maniac_file = .true.
                cycle

            case ('-d')
                if (seen_data_file) then
                    call abort_run("Duplicate option: -d", 1)
                end if
                call expect_value(i, nargs, "-d", path%topology)
                seen_data_file = .true.
                cycle

            case ('-p')
                if (seen_inc_file) then
                    call abort_run("Duplicate option: -p", 1)
                end if
                call expect_value(i, nargs, "-p", path%parameters)
                seen_inc_file = .true.
                cycle

            case ('-r')
                if (seen_res_file) then
                    call abort_run("Duplicate option: -r", 1)
                end if
                call expect_value(i, nargs, "-r", path%reservoir)
                seen_res_file = .true.
                cycle

            case ('-o')
                if (seen_output_path) then
                    call abort_run("Duplicate option: -o", 1)
                end if
                call expect_value(i, nargs, "-o", path%outputs)
                seen_output_path = .true.
                cycle
                
            case default ! If argument starts with '-', it's an unknown flag
                if (arg(1:1) == '-') then
                    call abort_run("Unknown option: "//trim(arg), 1)
                end if

            end select
        end do

        if (seen_res_file) status%reservoir_provided = .true.

    end subroutine read_command_line_args

    !-----------------------------------------------------------------
    ! Validate CLI arguments and check file existence
    !-----------------------------------------------------------------
    subroutine validate_CLI_arguments()

        if (trim(path%input) == '' .or. trim(path%topology) == '' .or. trim(path%parameters) == '') then
            call abort_run("Missing mandatory input arguments: -i, -d, -p required.", 1)
        end if

        if (.not. file_exists(path%input)) then
            call abort_run("Input file not found: "//trim(path%input), 1)
        end if
        
        if (.not. file_exists(path%topology)) then
            call abort_run("Data file not found: "//trim(path%topology), 1)
        end if 

        if (.not. file_exists(path%parameters)) then
            call abort_run("Parameter file not found: "//trim(path%parameters), 1)
        end if

        if (trim(path%reservoir) /= '' .and. .not. file_exists(path%reservoir)) then
            call abort_run("Reservoir file not found: "//trim(path%reservoir), 1)
        end if

    end subroutine validate_CLI_arguments

    !-----------------------------------------------------------------
    ! Ensure output path ends with a trailing slash
    !-----------------------------------------------------------------
    subroutine normalize_output_path()

        if (len_trim(path%outputs) > 0) then

            if (path%outputs(len_trim(path%outputs):len_trim(path%outputs)) /= '/') then
                path%outputs = trim(path%outputs) // '/'
            end if

        end if
    end subroutine normalize_output_path

    !-----------------------------------------------------------------
    ! Utility function to check if a file exists
    !-----------------------------------------------------------------
    logical function file_exists(filename)

        ! Input parameters
        character(len=*), intent(in) :: filename
    
        ! Local variables
        logical :: exists

        inquire(file=trim(filename), exist=exists)
        file_exists = exists

    end function file_exists

    !-----------------------------------------------------------------
    ! Reads the value following a CLI flag and validates it.
    !-----------------------------------------------------------------
    subroutine expect_value(idx, nargs, flag, value_out)
    
        ! Input parameters
        integer, intent(in)  :: idx, nargs
        character(len=*), intent(in)  :: flag
        character(len=*), intent(out) :: value_out
        
        ! Local variables
        character(len=256) :: tmp

        ! Ensure value exists after the flag
        if (idx >= nargs) then
            call abort_run("Missing value after option "//trim(flag), 1)
        end if

        ! Read next argument
        call get_command_argument(idx+1, tmp)

        ! Check empty value (e.g. "-i  ")
        if (len_trim(tmp) == 0) then
            call abort_run("Empty value after option "//trim(flag), 1)
        end if

        ! Check empty value (e.g. "-i  ")
        if (len_trim(tmp) == 0) then
            call abort_run("Empty value after option "//trim(flag), 1)
        end if

        value_out = trim(tmp)

    end subroutine expect_value

end module cli_utils
