module prescan_files

    use parameters
    use simulation_state
    use output_utils
    use check_utils
    use readers_utils
    use random_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !---------------------------------------------------------------------------
    ! Prescan the input file and the data file
    !---------------------------------------------------------------------------
    subroutine prescan_inputs()

        ! Pre-scan input files
        call prescan_input_file(path%input)
        call prescan_topology(path%topology, is_reservoir = .false.)
        call prescan_topology(path%reservoir, is_reservoir = .true.)

    end subroutine prescan_inputs

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
            call predetect_type(unit)
        
        close(unit)

    end subroutine prescan_input_file

    !-----------------------------------------------------------------------------
    ! Pre-scans a residue input file to determine the maximum sizes needed
    ! for array allocations
    !-----------------------------------------------------------------------------
    subroutine predetect_number_info(infile)

        ! Input parameters
        integer, intent(in) :: infile                 ! File unit to read residue definitions

        ! Local variables
        character(len=256) :: line                    ! Current line read from input file
        character(len=256) :: token                   ! First word of the line
        character(len=256) :: rest_line               ! Remaining part of the line after token
        integer :: ios                                ! I/O status of read operation
        integer :: ios_val                            ! Temporary I/O status for internal reads
        integer :: pos                                ! Position of first space in line
        integer :: len_rest                            ! Length of rest_line
        integer :: val                                ! Temporary value parsed from line
        integer :: pos_space                           ! Position of space for tokenization
        integer :: nb_type_per_residue                ! Number of types in current residue
        integer :: max_atom_in_residue                ! Max number of atoms in current residue
        logical :: in_residue_block                   ! True if inside a residue block
        logical :: is_active_residue                  ! True if residue is active

        ! Initialize counters
        nb%type_residue = 0
        nb%max_type_per_residue = 0
        nb%max_atom_in_any_residue = 0
        in_residue_block = .false.

        ! Loop over lines in input file
        do
            read(infile, '(A)', iostat=ios) line
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
                if (max_atom_in_residue > nb%max_atom_in_any_residue) then
                    nb%max_atom_in_any_residue = max_atom_in_residue
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
    ! Extract the atom-type list for each residue definition in the input file
    !-----------------------------------------------------------------------------
    subroutine predetect_type(infile)

        ! Input argument
        integer, intent(in) :: infile           ! Input file unit number

        ! Local variables
        character(len=256) :: line              ! Current line read from file
        character(len=256) :: token             ! First word of the line
        character(len=256) :: rest_line         ! Remaining part of the line
        integer :: ios                          ! I/O status of the read statement
        integer :: ios_val                      ! I/O status for numeric reads
        integer :: pos_space                    ! Position of first space in line
        integer :: pos                          ! Current parsing position in rest_line
        integer :: len_rest                     ! Length of rest_line
        integer :: residue_id                   ! Current residue index
        integer :: val                          ! Temporary value read from line
        logical :: in_residue_block             ! Flag: inside a residue block

        ! Safety: allocate res_infos
        if (allocated(res_infos)) deallocate(res_infos)
        allocate(res_infos(nb%type_residue))

        do residue_id = 1, nb%type_residue
            res_infos(residue_id)%nb_types = 0
            res_infos(residue_id)%nb_atoms = 0
            if (allocated(res_infos(residue_id)%types)) &
                deallocate(res_infos(residue_id)%types)
        end do

        rewind(infile)
        in_residue_block = .false.
        residue_id = 0

        ! Parse again to extract "types" entries
        do
            read(infile, '(A)', iostat=ios) line
            if (ios /= 0) exit

            line = adjustl(trim(line))

            select case (trim(line))
            case ('begin_residue')
                in_residue_block = .true.
                residue_id = residue_id + 1
                cycle

            case ('end_residue')
                in_residue_block = .false.
                cycle
            end select

            if (.not. in_residue_block) cycle

            ! -------------------------------------------------------
            ! Tokenize line
            ! -------------------------------------------------------
            pos_space = index(line, ' ')
            if (pos_space > 0) then
                token = line(1:pos_space-1)
                rest_line = adjustl(line(pos_space+1:))
            else
                token = trim(line)
                rest_line = ''
            end if

            ! -------------------------------------------------------
            ! Extract list of atom types
            ! -------------------------------------------------------
            if (trim(token) == 'types') then

                ! First count how many values
                pos = 1
                len_rest = len_trim(rest_line)
                res_infos(residue_id)%nb_types = 0

                do
                    read(rest_line(pos:), *, iostat=ios_val) val
                    if (ios_val /= 0) exit
                    res_infos(residue_id)%nb_types = res_infos(residue_id)%nb_types + 1
                    pos = pos + index(rest_line(pos:), ' ')
                    if (pos == 0 .or. pos > len_rest) exit
                end do

                ! Allocate the type array
                allocate(res_infos(residue_id)%types(res_infos(residue_id)%nb_types))

                ! Second pass: fill values
                pos = 1
                do val = 1, res_infos(residue_id)%nb_types
                    read(rest_line(pos:), *, iostat=ios_val) res_infos(residue_id)%types(val)
                    if (ios_val /= 0) exit
                    pos = pos + index(rest_line(pos:), ' ')
                    if (pos == 0) exit
                end do

            end if

            if (trim(token) == 'nb-atoms') then

                read(rest_line,*) res_infos(residue_id)%nb_atoms

            end if

        end do

    end subroutine predetect_type

    subroutine prescan_topology(filename, is_reservoir)

        ! Input parameters
        character(len=*), intent(in) :: filename      ! Name of the topology file
        logical, intent(in) :: is_reservoir           ! Flag: true if reservoir file

        ! Local variables
        integer :: unit                                ! Fortran unit number for file
        integer :: ios                                 ! I/O status
        character(len=256) :: line                     ! Current line read from file
        character(len=256) :: token                    ! First word of the line
        character(len=256) :: rest                     ! Remaining part of the line
        integer :: pos                                 ! Position of first space in line
        integer :: atom_id                             ! Atom index from file
        integer :: atom_type                           ! Atom type from file
        integer :: mol_id                              ! Molecule index from file
        integer :: residue_id                          ! Current residue index
        integer :: type_id                             ! Loop counter

        ! Temporary dynamic table for residue atom count
        integer, allocatable :: residue_atom_count(:)
        allocate(residue_atom_count(nb%type_residue))
        residue_atom_count = 0

        ! ---------------------------------------------------------
        ! Open topology file
        ! ---------------------------------------------------------
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        call check_IO_status(filename, ios)

        ! ---------------------------------------------------------
        ! Scan file
        ! ---------------------------------------------------------
        do
            read(unit,'(A)', iostat=ios) line
            if (ios /= 0) exit

            line = adjustl(trim(line))
            if (len_trim(line) == 0) cycle

            ! Extract token + remainder
            pos = index(line,' ')
            if (pos > 0) then
                token = trim(line(1:pos-1))
                rest  = adjustl(line(pos+1:))
            else
                cycle
            end if

            ! Note, the atome type is expected to be in second position
            ! Read atom-ID, molecule-ID, atom-type
            read(line,*,iostat=ios) atom_id, mol_id, atom_type
            if (ios /= 0) cycle

            ! ---------------------------------------------------------
            ! Classify this atom according to which residue type it belongs to
            ! ---------------------------------------------------------
            do residue_id = 1, nb%type_residue
                do type_id = 1, res_infos(residue_id)%nb_types
                    if (atom_type == res_infos(residue_id)%types(type_id)) then
                        residue_atom_count(residue_id) = residue_atom_count(residue_id) + 1
                    end if
                end do
            end do
        end do

        close(unit)

        ! ---------------------------------------------------------
        ! Compute number of residues of each type
        ! ---------------------------------------------------------
        if (is_reservoir) then
            do residue_id = 1, nb%type_residue
                if (res_infos(residue_id)%nb_atoms > 0) then
                    res_infos(residue_id)%nb_res(2) = residue_atom_count(residue_id) / res_infos(residue_id)%nb_atoms
                else
                    res_infos(residue_id)%nb_res(2) = 0
                end if
            end do
        else
            do residue_id = 1, nb%type_residue
                if (res_infos(residue_id)%nb_atoms > 0) then
                    res_infos(residue_id)%nb_res(1) = residue_atom_count(residue_id) / res_infos(residue_id)%nb_atoms
                else
                    res_infos(residue_id)%nb_res(1) = 0
                end if
            end do
        end if

    end subroutine prescan_topology

end module prescan_files
