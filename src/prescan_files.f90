module prescan_files

    use parameters
    use simulation_state
    use output_utils
    use check_utils
    use readers_utils
    use random_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    !---------------------------------------------------------------------------
    ! Information extracted during the prescan for each residue type.
    ! Stores the list of atom types, the number of types, and the expected atom
    ! count for that residue as defined in the Maniac input file.
    !---------------------------------------------------------------------------
    type residue
        integer, allocatable :: types(:)    ! List of atom types belonging to this residue
        integer :: nb_types                 ! Number of atom types in this residue
        integer :: nb_atoms                 ! Total number of atoms in this residue
        integer :: nb_res(2)                ! Total number of residue in primary and reservoir 
        logical :: is_active                ! The residue is active
    end type residue
    type(residue), allocatable :: res_infos(:)

contains

    !---------------------------------------------------------------------------
    ! Prescan the input file and the data file
    !---------------------------------------------------------------------------
    subroutine prescan_inputs()

        ! Initialize max number of residue
        nmax%active_residues = 0
        nmax%inactive_residues = 0

        ! Pre-scan input files
        call prescan_input_file(path%input) ! Detect residue definitions and sizes
        call prescan_topology(path%topology, is_reservoir = .false.) ! Count residues in primary topology
        if (status%reservoir_provided) then
            call prescan_topology(path%reservoir, is_reservoir = .true.) ! Count residues in reservoir topology
        else
            nmax%active_residues = NB_MAX_MOLECULE
        end if

    end subroutine prescan_inputs

    !---------------------------------------------------------------------------
    ! Open the main input file, pre-detect sizes and counts for allocation,
    ! then close the file.
    !---------------------------------------------------------------------------
    subroutine prescan_input_file(filename)

        ! Input parameter
        character(len=*), intent(in) :: filename

        ! Local variables
        integer :: unit                         ! Fortran unit number used to open the input file
        integer :: ios                          ! I/O status returned by open/read operations

        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)

        call check_IO_status(filename, ios)     ! Verify file opened correctly
        call predetect_number_info(unit)        ! Scan residue blocks for size limits
        call predetect_type(unit)               ! Extract per-residue atom-type lists

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
        res%number = 0
        nmax%types_per_residue = 0
        nmax%atoms_per_residue = 0
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
                res%number = res%number + 1
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
                if (nb_type_per_residue > nmax%types_per_residue) then
                    nmax%types_per_residue = nb_type_per_residue
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
                if (max_atom_in_residue > nmax%atoms_per_residue) then
                    nmax%atoms_per_residue = max_atom_in_residue
                end if

                ! Update active / inactive maxima
                if (is_active_residue) then
                    if (max_atom_in_residue > nmax%atoms_active_residue) then
                        nmax%atoms_active_residue = max_atom_in_residue
                    end if
                else
                    if (max_atom_in_residue > nmax%atoms_inactive_residue) then
                        nmax%atoms_inactive_residue = max_atom_in_residue
                    end if
                end if
            end if
        end do

    end subroutine predetect_number_info

    !-----------------------------------------------------------------------------
    ! Scans the residue definition blocks again to extract the explicit list of
    ! atom types, activation state, and declared atom counts for each residue.
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
        integer :: res_id                   ! Current residue index
        integer :: val                          ! Temporary value read from line
        logical :: in_residue_block             ! Flag: inside a residue block

        if (allocated(res_infos)) deallocate(res_infos)
        allocate(res_infos(res%number))

        do res_id = 1, res%number
            res_infos(res_id)%nb_res = 0
            res_infos(res_id)%nb_types = 0
            res_infos(res_id)%nb_atoms = 0
            if (allocated(res_infos(res_id)%types)) deallocate(res_infos(res_id)%types)
        end do

        rewind(infile)
        in_residue_block = .false.
        res_id = 0

        ! Parse again to extract "types" entries
        do
            read(infile, '(A)', iostat=ios) line
            if (ios /= 0) exit

            line = adjustl(trim(line))

            select case (trim(line))
            case ('begin_residue')
                in_residue_block = .true.
                res_id = res_id + 1
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
                res_infos(res_id)%nb_types = 0

                do
                    read(rest_line(pos:), *, iostat=ios_val) val
                    if (ios_val /= 0) exit
                    res_infos(res_id)%nb_types = res_infos(res_id)%nb_types + 1
                    pos = pos + index(rest_line(pos:), ' ')
                    if (pos == 0 .or. pos > len_rest) exit
                end do

                ! Allocate the type array
                allocate(res_infos(res_id)%types(res_infos(res_id)%nb_types))

                ! Second pass: fill values
                pos = 1
                do val = 1, res_infos(res_id)%nb_types
                    read(rest_line(pos:), *, iostat=ios_val) res_infos(res_id)%types(val)
                    if (ios_val /= 0) exit
                    pos = pos + index(rest_line(pos:), ' ')
                    if (pos == 0) exit
                end do

            end if

            if (trim(token) == 'state') then
                if (trim(rest_line) == 'actif') then
                    res_infos(res_id)%is_active = .true.
                else
                    res_infos(res_id)%is_active = .false.
                end if
            end if

            if (trim(token) == 'nb-atoms') then
                read(rest_line,*) res_infos(res_id)%nb_atoms
            end if

        end do

    end subroutine predetect_type

    !-----------------------------------------------------------------------------
    ! Parses the topology file to count how many atoms of each residue type
    ! actually appear in the system.
    !-----------------------------------------------------------------------------
    subroutine prescan_topology(filename, is_reservoir)

        ! Input parameters
        character(len=*), intent(in) :: filename        ! Name of the topology file
        logical, intent(in) :: is_reservoir             ! Flag: true if reservoir file

        ! Local variables
        integer :: unit                                 ! Fortran unit number for file
        integer :: ios                                  ! I/O status
        character(len=256) :: line                      ! Current line read from file
        character(len=256) :: token                     ! First word of the line
        character(len=256) :: rest                      ! Remaining part of the line
        integer :: pos                                  ! Position of first space in line
        integer :: atom_id                              ! Atom index from file
        integer :: atom_type                            ! Atom type from file
        integer :: mol_id                               ! Molecule index from file
        real(real64) :: q, x, y, z                      ! Charge and coordinate
        integer :: res_id                               ! Current residue index
        integer :: type_id                              ! Loop counter
        logical :: found                                ! Logical indicator for atom counting
        integer, allocatable :: residue_atom_count(:)   ! Temporary dynamic table for residue atom count

        allocate(residue_atom_count(res%number))
        residue_atom_count = 0

        ! Open topology file
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        call check_IO_status(filename, ios)

        ! Scan file
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

            ! Read line (the atom type is expected to be in second position)
            read(line,*,iostat=ios) atom_id, mol_id, atom_type, q, x, y, z
            if (ios /= 0) cycle

            ! Classify this atom according to which residue type it belongs to
            found = .false.
            do res_id = 1, res%number
                do type_id = 1, res_infos(res_id)%nb_types
                    if (atom_type == res_infos(res_id)%types(type_id)) then
                        residue_atom_count(res_id) = residue_atom_count(res_id) + 1
                        found = .true.
                        exit
                    end if
                end do
                if (found) exit
            end do
        end do

        close(unit)

        call compute_max_residue(residue_atom_count, is_reservoir)

    end subroutine prescan_topology

    !-----------------------------------------------------------------------------
    ! Compute max numbers of residues used for memory allocation.
    !-----------------------------------------------------------------------------
    subroutine compute_max_residue(residue_atom_count, is_reservoir)

        ! Input parameters
        integer, intent(in) :: residue_atom_count(:)    ! Temporary dynamic table for residue atom count
        logical, intent(in) :: is_reservoir             ! Flag: true if reservoir file

        ! Local variable
        integer :: res_id                               ! Current residue index
        integer :: nb_residue                           ! Local variable for number of residue

        ! Compute number of residues
        do res_id = 1, res%number

            if (res_infos(res_id)%nb_atoms <= 0) cycle  ! Skip invalid residue

            ! Compute main/reservoir counts as real then round
            if (is_reservoir) then
                res_infos(res_id)%nb_res(2) = nint( real(residue_atom_count(res_id), kind=real64) &
                                                    / real(res_infos(res_id)%nb_atoms, kind=real64) )
            else
                res_infos(res_id)%nb_res(1) = nint( real(residue_atom_count(res_id), kind=real64) &
                                                    / real(res_infos(res_id)%nb_atoms, kind=real64) )
            end if

            ! Check for anomalous values
            if (res_infos(res_id)%nb_res(1) < 0 .or. res_infos(res_id)%nb_res(2) < 0) then
                call abort_run("Error: Determinant fell into denormal/underflow range")
            end if

        end do

        ! Determine maximum active and inactive residues
        do res_id = 1, res%number
            nb_residue = res_infos(res_id)%nb_res(1) + res_infos(res_id)%nb_res(2)            
            if (res_infos(res_id)%is_active) then
                if (nb_residue > nmax%active_residues) then
                    nmax%active_residues = nb_residue
                end if
            else
                if (nb_residue > nmax%inactive_residues) then
                    nmax%inactive_residues = nb_residue
                end if
            end if
        end do

    end subroutine compute_max_residue

end module prescan_files
