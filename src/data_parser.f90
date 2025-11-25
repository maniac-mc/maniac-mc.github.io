module data_parser

    use parameters
    use simulation_state
    use output_utils
    use geometry_utils
    use check_utils
    use readers_utils
    use sanity_utils

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    ! 1D atom information
    real(real64), allocatable :: atom_xyz(:,:)                              ! Array coordinate of atoms
    real(real64), allocatable :: atom_charges_1d(:)                         ! Partial charges on sites

    real(real64), allocatable :: atom_masses_1d(:)                          ! Masses of atoms
    integer, allocatable :: atom_types_1d(:)                                ! Array of atom types
    integer, allocatable :: atom_ids_1d(:)                                  ! Atom_id
    integer, allocatable :: atom_original_1d(:)                             ! Atom_id (from data file)
    character(10), allocatable :: atom_names_1d(:)                          ! Array of atom names

    integer, allocatable :: bond_ids_1d(:)
    integer, allocatable :: bond_types_1d(:)
    integer, allocatable :: bond_atoms_1d(:, :)

    integer, allocatable :: angle_ids_1d(:)
    integer, allocatable :: angle_types_1d(:)
    integer, allocatable :: angle_atoms_1d(:, :)

    integer, allocatable :: dihedral_ids_1d(:)
    integer, allocatable :: dihedral_types_1d(:)
    integer, allocatable :: dihedral_atoms_1d(:, :)

    integer, allocatable :: improper_ids_1d(:)
    integer, allocatable :: improper_types_1d(:)
    integer, allocatable :: improper_atoms_1d(:, :)

contains

    !-----------------------------------------------------------------------------
    ! Reads the main system (primary) data file and optionally a reservoir file.
    ! Ensures consistency between primary and reservoir if applicable.
    !-----------------------------------------------------------------------------
    subroutine read_system_data()

        ! Read primary system
        call read_lmp_data(path%topology, primary, is_primary = .true.)

        ! Read reservoir if provided
        if (trim(path%reservoir) == '') then
            status%reservoir_provided = .false.
        else
            status%reservoir_provided = .true.
            call read_lmp_data(path%reservoir, reservoir, is_primary = .false.)
            call AssertMassConsistency()                ! Ensure masses in reservoir are consistent with primary
        end if

    end subroutine read_system_data

    subroutine read_lmp_data(data_file_name, box, is_primary)

        ! Input parameters
        character(len=*), intent(in) :: data_file_name  ! Path to the LAMMPS data file to read
        type(type_box), intent(inout) :: box            ! Type box (primary or reservoir)
        logical, intent(in) :: is_primary               ! Flag for primary/reservoir reading

        ! Local variables
        integer :: infile = 100                         ! Fortran unit number used to open the input file
        integer :: ios                                  ! I/O status flag (nonzero if file open/read fails)

        ! === Step 1: Open file ===
        open(UNIT=infile, FILE=data_file_name, STATUS='OLD', ACTION='read', IOSTAT=ios)

        if (ios /= 0) then
            call abort_run("Error opening file: " // trim(data_file_name))
        end if

        ! === Step 2: Read header ===
        call read_lmp_header_info(infile, box)

        ! === Step 3: Allocate arrays ===
        allocate(atom_ids_1d(box%num_atoms))
        allocate(atom_original_1d(box%num_atoms))
        allocate(atom_types_1d(box%num_atoms))
        allocate(atom_charges_1d(box%num_atoms))
        allocate(atom_xyz(3, box%num_atoms))
        allocate(atom_names_1d(box%num_atoms))
        allocate(bond_ids_1d(box%num_bonds))
        allocate(bond_types_1d(box%num_bonds))
        allocate(bond_atoms_1d(box%num_bonds,2))
        allocate(angle_ids_1d(box%num_angles))
        allocate(angle_types_1d(box%num_angles))
        allocate(angle_atoms_1d(box%num_angles,3))
        allocate(dihedral_ids_1d(box%num_dihedrals))
        allocate(dihedral_types_1d(box%num_dihedrals))
        allocate(dihedral_atoms_1d(box%num_dihedrals,4))
        allocate(improper_ids_1d(box%num_impropers))
        allocate(improper_types_1d(box%num_impropers))
        allocate(improper_atoms_1d(box%num_impropers,4))

        ! === Step 4: Read box dimensions ===
        rewind(infile)
        call parse_lammps_box(infile, box)

        ! Detect box type and other informations
        call prepare_simulation_box(box)

        ! === Step 5: Read masses ===
        rewind(infile)
        call read_lammps_masses(infile, box)

        ! === Step 6: Read atoms ===
        rewind(infile)
        call read_lammps_atoms(infile, data_file_name, box)
        call assign_atom_names(box)

        ! === Step 6: Read bonds ===
        rewind(infile)
        call read_lammps_bonds(infile, data_file_name, box)

        ! === Step 6: Read angles ===
        rewind(infile)
        call read_lammps_angles(infile, data_file_name, box)

        ! === Step 7: Read dihedrals ===
        rewind(infile)
        call read_lammps_dihedrals(infile, data_file_name, box)

        ! === Step 7: Read improper ===
        rewind(infile)
        call read_lammps_impropers(infile, data_file_name, box)

        close(infile)

        ! === Step 7: Post-processing ===
        call sort_atoms_by_original_ID(box)
        call detect_molecules(box)
        call repair_active_molecules(box)

        ! #tofix : simplify by only using "is_reservoir" everywhere
        if (is_primary) then
            call transform_to_COM_frame(box, .false.)
            call validate_molecule_geometry(box, .false.)
        else
            call transform_to_COM_frame(box, .true.)
            call validate_molecule_geometry(box, .true.)
        end if 

        ! === Step 7: Detect bond, angle, dihedral, improper per residue ===
        call DetectBondPerResidue(box)
        call DetectAnglePerResidue(box)
        call DetectDihedralPerResidue(box)
        call DetectImproperPerResidue(box)

        ! === Step 8: Deallocate arrays ===
        if (allocated(atom_ids_1d)) deallocate(atom_ids_1d)
        if (allocated(atom_original_1d)) deallocate(atom_original_1d)
        if (allocated(atom_types_1d)) deallocate(atom_types_1d)
        if (allocated(atom_charges_1d)) deallocate(atom_charges_1d)
        if (allocated(atom_xyz)) deallocate(atom_xyz)
        if (allocated(atom_names_1d)) deallocate(atom_names_1d)
        if (allocated(bond_ids_1d)) deallocate(bond_ids_1d)
        if (allocated(bond_types_1d)) deallocate(bond_types_1d)
        if (allocated(bond_atoms_1d)) deallocate(bond_atoms_1d)
        if (allocated(angle_ids_1d)) deallocate(angle_ids_1d)
        if (allocated(angle_types_1d)) deallocate(angle_types_1d)
        if (allocated(angle_atoms_1d)) deallocate(angle_atoms_1d)
        if (allocated(dihedral_ids_1d)) deallocate(dihedral_ids_1d)
        if (allocated(dihedral_types_1d)) deallocate(dihedral_types_1d)
        if (allocated(dihedral_atoms_1d)) deallocate(dihedral_atoms_1d)
        if (allocated(improper_ids_1d)) deallocate(improper_ids_1d)
        if (allocated(improper_types_1d)) deallocate(improper_types_1d)
        if (allocated(improper_atoms_1d)) deallocate(improper_atoms_1d)

        ! === Setep 9: Log data ===
        call LogData(data_file_name, box, is_primary)
        call LogConnectivity(box)

    end subroutine read_lmp_data

    !-----------------------------------------------------------------------------
    ! Reads atomic masses from the "Masses" section of a LAMMPS data file and maps them
    ! into a 2D mass array for each atom type and residue type.
    !-----------------------------------------------------------------------------
    subroutine read_lammps_masses(infile, box)

        ! Input parameters
        type(type_box), intent(inout) :: box         ! Simulation box containing atom/mass info
        integer, intent(in) :: infile                ! LAMMPS data file unit to read from

        ! Local variables
        character(len=256) :: line, trimmed_line     ! Raw and trimmed lines from file
        integer :: ios                               ! I/O status of read operations
        integer :: i, j, k                           ! Loop indices for atoms, residues, types
        integer :: tmp_int                           ! Temporary integer for atom type ID
        integer :: mass_found                        ! Counter for number of masses successfully read
        real(real64) :: tmp_flt                      ! Temporary variable for atomic mass

        ! Initialize
        mass_found = 0

        allocate(box%site_masses_vector(box%num_atomtypes))
        box%site_masses_vector = 0.0d0

        ! Rewind file to start (optional, remove if file pointer is already positioned)
        rewind(infile)

        ! Look for "Masses" section
        do
            read(infile, '(A)', IOSTAT=ios) line
            if (ios /= 0) then
                exit
            end if

            ! Trim leading and trailing spaces
            trimmed_line = adjustl(trim(line))

            ! Skip empty lines or comments
            if (len_trim(trimmed_line) == 0 .or. trimmed_line(1:1) == '!') then
                cycle
            end if

            ! Check for "Masses" keyword
            if (trim(trimmed_line) == "Masses") then
                ! Skip next line (usually blank or a comment)
                read(infile, '(A)', IOSTAT=ios) line
                if (ios /= 0) then
                    exit
                end if
                trimmed_line = adjustl(trim(line))
                if (len_trim(trimmed_line) > 0 .and. trimmed_line(1:1) /= '!') then
                end if

                ! Read num_atoms(SYSTEM) lines for masses
                do i = 1, box%num_atoms
                    read(infile, '(A)', IOSTAT=ios) line
                    if (ios /= 0) then
                        exit
                    end if

                    trimmed_line = adjustl(trim(line))
                    ! Skip empty or comment lines
                    if (len_trim(trimmed_line) == 0 .or. trimmed_line(1:1) == '!') then
                        cycle
                    end if

                    ! Parse species ID and mass
                    read(trimmed_line, *, IOSTAT=ios) tmp_int, tmp_flt
                    if (ios /= 0) then
                        exit
                    end if

                    if (tmp_int < 1 .or. tmp_int > box%num_atoms) then
                        exit
                    end if

                    box%site_masses_vector(tmp_int) = tmp_flt
                    mass_found = mass_found + 1
                end do

                ! Map site masses into 2D mass array (assuming nb%type_residue and nb%types_per_residue are defined)
                i = 1
                do j = 1, nb%type_residue
                    do k = 1, nb%types_per_residue(j)
                        if ((i <= box%num_atoms) .and. (i <= box%num_atomtypes)) then
                            box%atom_masses(j, k) = box%site_masses_vector(i)
                            i = i + 1
                        else
                            exit
                        end if
                    end do
                end do

                exit  ! Exit after processing Masses section
            end if
        end do

        ! Trigger error if no masses found or mismatch with declared atom types
        if (mass_found == 0) then
            call abort_run("No masses found in data file", 12)
        else if (mass_found /= box%num_atomtypes) then
            call abort_run("Number of masses found in data file differs from declared atom types", 13)
        end if

    end subroutine read_lammps_masses

    !-----------------------------------------------------------------------------
    ! Assigns atom names based on their types by matching each atom's type with
    ! the residue type definitions. If no match is found, assigns a default name
    ! "Unknown".
    !-----------------------------------------------------------------------------
    subroutine assign_atom_names(box)

        ! Input/Output
        type(type_box), intent(inout) :: box         ! Simulation box containing atom types and names

        ! Local variables
        integer :: i, j, k                            ! Loop indices for residues, types, and atoms

        ! Loop over all atoms
        do k = 1, box%num_atoms
            atom_names_1d(k) = '' ! Initialize to empty string
            do i = 1, nb%type_residue
                do j = 1, nb%types_per_residue(i)
                    if (res%types_2d(i, j) == atom_types_1d(k)) then
                        atom_names_1d(k) = res%names_2d(i, j)
                        exit  ! Exit inner loop once match is found
                    end if
                end do
                if (len_trim(atom_names_1d(k)) > 0) exit  ! Exit outer loop if found
            end do
            if (len_trim(atom_names_1d(k)) == 0) then
                atom_names_1d(k) = 'Unknown'  ! Optional default value
            end if
        end do

    end subroutine assign_atom_names

    subroutine DetectBondPerResidue(box)

        type(type_box), intent(inout) :: box

        integer :: i, j
        integer :: id1, id2
        integer :: k, l
        integer :: nb_bond
        character(len=100) :: formatted_msg

        do i = 1, nb%type_residue
            nb_bond = 0

            ! Loop directly over the bonds
            do j = 1, box%num_bonds

                ! Ids in bonds
                id1 = bond_atoms_1d(j,1)
                id2 = bond_atoms_1d(j,2)

                ! Check if both atoms belong to residue i
                if (isInResidue(box, i, id1) .and. isInResidue(box, i, id2)) then

                    ! Get atom local index in the residue
                    k = atomIndexInResidue(box, i, id1)
                    l = atomIndexInResidue(box, i, id2)
                    ! Store bond info
                    nb_bond = nb_bond + 1

                    if (nb_bond > NB_MAX_BOND) then
                        write(formatted_msg, '(A, I0)') &
                            "The number of bonds exceeds the maximum allowed = ", NB_MAX_BOND
                        call abort_run(trim(formatted_msg) // &
                            " Please increase 'NB_MAX_BOND' in src/parameters.f90 and recompile.", 11)
                    end if

                    res%bond_type_2d(i, nb_bond, 1) = bond_types_1d(j)

                    if (l < k) then
                        res%bond_type_2d(i, nb_bond, 2) = l
                        res%bond_type_2d(i, nb_bond, 3) = k
                    else
                        res%bond_type_2d(i, nb_bond, 2) = k
                        res%bond_type_2d(i, nb_bond, 3) = l
                    end if
                end if

            end do

            nb%bonds_per_residue(i) = nb_bond

        end do

    end subroutine DetectBondPerResidue

    subroutine DetectAnglePerResidue(box)

        type(type_box), intent(inout) :: box

        integer :: i, j, k, l, m
        integer :: id1, id2, id3
        integer :: nb_angle
        character(len=100) :: formatted_msg

        do i = 1, nb%type_residue
            nb_angle = 0

            ! Loop directly over the angles
            do j = 1, box%num_angles
                id1 = angle_atoms_1d(j,1)
                id2 = angle_atoms_1d(j,2)
                id3 = angle_atoms_1d(j,3)

                ! Check if all 3 atoms belong to residue i
                if (isInResidue(box, i, id1) .and. isInResidue(box, i, id2) .and. isInResidue(box, i, id3)) then

                    ! Get local indices inside the residue
                    k = atomIndexInResidue(box, i, id1)
                    l = atomIndexInResidue(box, i, id2)
                    m = atomIndexInResidue(box, i, id3)

                    ! Store angle info
                    nb_angle = nb_angle + 1

                    if (nb_angle > NB_MAX_ANGLE) then
                        write(formatted_msg, '(A, I0)') &
                            "The number of angles exceeds the maximum allowed = ", NB_MAX_ANGLE
                        call abort_run(trim(formatted_msg) // &
                            " Please increase 'NB_MAX_ANGLE' in src/parameters.f90 and recompile.", 11)
                    end if

                    res%angle_type_2d(i, nb_angle, 1) = angle_types_1d(j)
                    if (k < m) then
                        res%angle_type_2d(i, nb_angle, 2) = k
                        res%angle_type_2d(i, nb_angle, 3) = l
                        res%angle_type_2d(i, nb_angle, 4) = m
                    else
                        res%angle_type_2d(i, nb_angle, 2) = m
                        res%angle_type_2d(i, nb_angle, 3) = l
                        res%angle_type_2d(i, nb_angle, 4) = k
                    end if
                end if
            end do
            nb%angles_per_residue(i) = nb_angle
        end do

    end subroutine DetectAnglePerResidue

    subroutine DetectDihedralPerResidue(box)

        type(type_box), intent(inout) :: box

        integer :: i, j, k, l, m, n
        integer :: id1, id2, id3, id4
        integer :: nb_dihedral
        character(len=100) :: formatted_msg

        do i = 1, nb%type_residue
            nb_dihedral = 0

            ! Loop directly over the dihedrals
            do j = 1, box%num_dihedrals
                id1 = dihedral_atoms_1d(j,1)
                id2 = dihedral_atoms_1d(j,2)
                id3 = dihedral_atoms_1d(j,3)
                id4 = dihedral_atoms_1d(j,4)

                ! Check if all 4 atoms belong to residue i
                if (isInResidue(box, i, id1) .and. isInResidue(box, i, id2) .and. &
                    isInResidue(box, i, id3) .and. isInResidue(box, i, id4)) then

                    ! Get local indices inside the residue
                    k = atomIndexInResidue(box, i, id1)
                    l = atomIndexInResidue(box, i, id2)
                    m = atomIndexInResidue(box, i, id3)
                    n = atomIndexInResidue(box, i, id4)

                    ! Store dihedral info
                    nb_dihedral = nb_dihedral + 1

                    if (nb_dihedral > NB_MAX_DIHEDRAL) then
                        write(formatted_msg, '(A, I0)') &
                            "The number of dihedrals exceeds the maximum allowed = ", NB_MAX_DIHEDRAL
                        call abort_run(trim(formatted_msg) // &
                            " Please increase 'NB_MAX_DIHEDRAL' in src/parameters.f90 and recompile.", 11)
                    end if

                    res%dihedral_type_2d(i, nb_dihedral, 1) = dihedral_types_1d(j)
                    if (k < n) then
                        res%dihedral_type_2d(i, nb_dihedral, 2) = k
                        res%dihedral_type_2d(i, nb_dihedral, 3) = l
                        res%dihedral_type_2d(i, nb_dihedral, 4) = m
                        res%dihedral_type_2d(i, nb_dihedral, 5) = n
                    else
                        res%dihedral_type_2d(i, nb_dihedral, 2) = n
                        res%dihedral_type_2d(i, nb_dihedral, 3) = m
                        res%dihedral_type_2d(i, nb_dihedral, 4) = l
                        res%dihedral_type_2d(i, nb_dihedral, 5) = k
                    end if
                end if
            end do
            nb%dihedrals_per_residue(i) = nb_dihedral
        end do

    end subroutine DetectDihedralPerResidue

    subroutine DetectImproperPerResidue(box)

        type(type_box), intent(inout) :: box

        integer :: i, j, k, l, m, n
        integer :: id1, id2, id3, id4
        integer :: nb_improper
        character(len=100) :: formatted_msg

        do i = 1, nb%type_residue
            nb_improper = 0

            ! Loop directly over the impropers
            do j = 1, box%num_impropers
                id1 = improper_atoms_1d(j,1)
                id2 = improper_atoms_1d(j,2)
                id3 = improper_atoms_1d(j,3)
                id4 = improper_atoms_1d(j,4)

                ! Check if all 4 atoms belong to residue i
                if (isInResidue(box, i, id1) .and. isInResidue(box, i, id2) .and. &
                    isInResidue(box, i, id3) .and. isInResidue(box, i, id4)) then

                    ! Get local indices inside the residue
                    k = atomIndexInResidue(box, i, id1)
                    l = atomIndexInResidue(box, i, id2)
                    m = atomIndexInResidue(box, i, id3)
                    n = atomIndexInResidue(box, i, id4)

                    ! Store improper info
                    nb_improper = nb_improper + 1

                    if (nb_improper > NB_MAX_IMPROPER) then
                        write(formatted_msg, '(A, I0)') &
                            "The number of impropers exceeds the maximum allowed = ", NB_MAX_IMPROPER
                        call abort_run(trim(formatted_msg) // &
                            " Please increase 'NB_MAX_IMPROPER' in src/parameters.f90 and recompile.", 11)
                    end if

                    res%improper_type_2d(i, nb_improper, 1) = improper_types_1d(j)
                    if (k < n) then
                        res%improper_type_2d(i, nb_improper, 2) = k
                        res%improper_type_2d(i, nb_improper, 3) = l
                        res%improper_type_2d(i, nb_improper, 4) = m
                        res%improper_type_2d(i, nb_improper, 5) = n
                    else
                        res%improper_type_2d(i, nb_improper, 2) = n
                        res%improper_type_2d(i, nb_improper, 3) = m
                        res%improper_type_2d(i, nb_improper, 4) = l
                        res%improper_type_2d(i, nb_improper, 5) = k
                    end if
                end if
            end do
            nb%impropers_per_residue(i) = nb_improper
        end do

    end subroutine DetectImproperPerResidue

    !-----------------------------------------------------------------------------
    ! Reads the "Atoms" section from a LAMMPS data file and fills the simulation box
    ! with atom IDs, types, charges, and positions. Performs extensive error checks
    ! for missing sections, invalid atom types, parsing errors, and array bounds.
    !-----------------------------------------------------------------------------
    subroutine read_lammps_atoms(infile, data_file_name, box)

        ! Input parameters
        integer, intent(in) :: infile                   ! LAMMPS data file unit to read atom data from
        character(len=*), intent(in) :: data_file_name  ! Name of the LAMMPS data file
        type(type_box), intent(inout) :: box           ! Simulation box to store atom properties

        ! Local variables
        character(len=200) :: line, formatted_msg       ! Raw line read and formatted error messages
        integer :: ios, k                               ! I/O status and loop index for atoms
        integer :: atom_found                            ! Counter for successfully read atoms
        integer :: tmp1_int, tmp2_int, tmp3_int         ! Temporary integers: atom ID, mol-ID, type
        real(real64) :: tmp1_flt, tmp2_flt, tmp3_flt, tmp4_flt  ! Temporary floats: charge, x, y, z
        integer, parameter :: ERROR_NO_ATOMS = 13       ! Error code: no atoms found
        integer, parameter :: ERROR_READ = 14           ! Error code: read failure
        integer, parameter :: ERROR_PARSE = 15          ! Error code: parse failure
        integer, parameter :: ERROR_ATOM_TYPE = 16      ! Error code: invalid atom type
        integer, parameter :: ERROR_BOUNDS = 17         ! Error code: array index out of bounds

        atom_found = 0

        ! Read until "Atoms" section or end of file
        do
            read(infile, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, A)') "No atoms found in data file: ", trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_NO_ATOMS)
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, A)') "I/O error reading data file: ", trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            ! Check for "Atoms" or "Atoms # full" (case-insensitive)
            if (index(adjustl(line), 'Atoms') == 1) then
                exit
            end if
        end do

        ! Read and discard a blank line, if present
        read(infile, '(A)', IOSTAT=ios) line
        if (ios < 0) then
            write(formatted_msg, '(A, A)') "No atoms found in data file: ", trim(data_file_name)
            call abort_run(trim(formatted_msg), ERROR_READ)
        end if

        if (TRIM(adjustl(line)) /= "") then
            ! If not blank, rewind to re-read as atom data
            backspace(infile, IOSTAT=ios)
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to rewind file after non-blank line in: ", trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if
        end if

        ! Loop over all atoms
        do k = 1, box%num_atoms
            read(infile, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, I0, A)') &
                    "Unexpected end of file at atom line ", k, " in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, I0, A)') &
                    "I/O error reading atom line ", k, " in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            ! Parse atom info: id, molecule-ID, type, q, x, y, z
            read(line, *, IOSTAT=ios) tmp1_int, tmp2_int, tmp3_int, &
                                    tmp1_flt, tmp2_flt, tmp3_flt, tmp4_flt

            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to parse atom line: '" &
                    // TRIM(line) // "' in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_PARSE)
            end if

            ! Validate atom type
            if (tmp3_int < 1 .or. tmp3_int > box%num_atomtypes) then
                write(formatted_msg, '(A, I0, A, I0, A)') "Invalid atom type ", tmp3_int, &
                                                        " (max allowed: ", box%num_atomtypes, &
                                                        ") in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_ATOM_TYPE)
            end if

            ! Validate array bounds
            if (k < 1 .or. k > size(atom_ids_1d)) then
                write(formatted_msg, '(A, I0, A, I0, A)') "Invalid array index k=", k, &
                                                        " (array size: ", size(atom_ids_1d), &
                                                        ") in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_BOUNDS)
            end if

            ! Store atom properties
            atom_ids_1d(k) = k
            atom_original_1d(k) = tmp1_int
            atom_types_1d(k) = tmp3_int
            atom_charges_1d(k) = tmp1_flt
            atom_xyz(1, k) = tmp2_flt
            atom_xyz(2, k) = tmp3_flt
            atom_xyz(3, k) = tmp4_flt
            atom_found = atom_found + 1
        end do

        ! Check if the number of atoms found matches the expected number
        if (atom_found /= box%num_atoms) then
            write(formatted_msg, '(A, I0, A, I0, A)') "Found ", atom_found, &
                                                    " atoms, expected ", box%num_atoms, &
                                                    " in: " // trim(data_file_name)
            call abort_run(trim(formatted_msg), ERROR_NO_ATOMS)
        end if

    end subroutine read_lammps_atoms

    !-----------------------------------------------------------------------------
    ! Reads bond definitions from the "Bonds" section of a LAMMPS data file and
    ! stores bond IDs, types, and connected atom indices in the corresponding arrays.
    ! If no bonds are expected, the routine returns immediately.
    !-----------------------------------------------------------------------------
    subroutine read_lammps_bonds(infile, data_file_name, box)

        ! Input
        integer, intent(in) :: infile                    ! Fortran unit number of the LAMMPS data file
        character(len=*), intent(in) :: data_file_name   ! Path to the LAMMPS data file
        type(type_box), intent(inout) :: box             ! Box structure containing bond counts and arrays

        ! Local variables
        character(len=200) :: line, formatted_msg        ! Current line and formatted error messages
        integer :: ios, k                                 ! I/O status and loop index
        integer :: bond_found                              ! Counter for bonds read
        integer :: tmp1_int, tmp2_int, tmp3_int, tmp4_int ! Temporary variables for bond id, type, atom1, atom2
        integer, parameter :: ERROR_READ = 24             ! Error code for I/O errors
        integer, parameter :: ERROR_PARSE = 25            ! Error code for parsing errors

        bond_found = 0

        ! If no bonds are expected, skip the whole routine
        if (box%num_bonds == 0) then
            call InfoMessage("No bonds expected in data file: " // trim(data_file_name))
            return
        end if

        ! Read until "Bonds" section or end of file
        do
            read(infile, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, A)') "No bonds found in data file: ", trim(data_file_name)
                call InfoMessage(trim(formatted_msg))
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, A)') "I/O error reading data file: ", trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            ! Check for "Bonds" or "Bonds # full" (case-insensitive)
            if (index(adjustl(line), 'Bonds') == 1) then
                exit
            end if

        end do

        ! Read and discard a blank line, if present
        read(infile, '(A)', IOSTAT=ios) line

        if (ios < 0) then
            write(formatted_msg, '(A, A)') "No bonds found in data file: ", trim(data_file_name)
            call abort_run(trim(formatted_msg), ERROR_READ)
        end if

        if (TRIM(adjustl(line)) /= "") then
            ! If not blank, rewind to re-read as bond data
            backspace(infile, IOSTAT=ios)
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to rewind file after non-blank line in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if
        end if

        ! Loop over all bonds
        do k = 1, box%num_bonds
            read(infile, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, I0, A)') "Unexpected end of file at bond line ", k, " in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, I0, A)') "I/O error reading bond line ", k, " in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            ! Parse bond info: id, type, atom1, atom2
            read(line, *, IOSTAT=ios) tmp1_int, tmp2_int, tmp3_int, tmp4_int
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to parse bond line: '" // TRIM(line) // "' in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_PARSE)
            end if

            ! Store bond properties
            bond_ids_1d(k) = tmp1_int
            bond_types_1d(k) = tmp2_int
            bond_atoms_1d(k,1) = tmp3_int
            bond_atoms_1d(k,2) = tmp4_int
            bond_found = bond_found + 1
        end do

    end subroutine read_lammps_bonds

    !-----------------------------------------------------------------------------
    ! Reads angle definitions from the "Angles" section of a LAMMPS data file
    ! and stores angle IDs, types, and connected atom indices in the corresponding arrays.
    ! If no angles are expected, the routine returns immediately.
    !-----------------------------------------------------------------------------
    subroutine read_lammps_angles(infile, data_file_name, box)

        ! Input
        integer, intent(in) :: infile                    ! Fortran unit number of the LAMMPS data file
        character(len=*), intent(in) :: data_file_name   ! Path to the LAMMPS data file
        type(type_box), intent(inout) :: box             ! Box structure containing angle counts and arrays

        ! Local variables
        character(len=200) :: line, formatted_msg        ! Current line and formatted error messages
        integer :: ios, k                                 ! I/O status and loop index
        integer :: angle_found                             ! Counter for angles read
        integer :: tmp1_int, tmp2_int, tmp3_int, tmp4_int, tmp5_int ! Temporary variables for angle id, type, atoms
        integer, parameter :: ERROR_READ = 24             ! Error code for I/O errors
        integer, parameter :: ERROR_PARSE = 25            ! Error code for parsing errors

        angle_found = 0

        ! If no angles are expected, skip the whole routine
        if (box%num_angles == 0) then
            write(formatted_msg, '(A, A)') "No angles expected in data file: ", trim(data_file_name)
            call InfoMessage(trim(formatted_msg))
            return
        end if

        ! Read until "Angles" section or end of file
        do
            read(infile, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, A)') "No angles found in data file: ", trim(data_file_name)
                call InfoMessage(trim(formatted_msg))
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, A)') "I/O error reading data file: ", trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            ! Check for "Angles" or "Angles # full" (case-insensitive)
            if (index(adjustl(line), 'Angles') == 1) then
                exit
            end if
        end do

        ! Read and discard a blank line, if present
        read(infile, '(A)', IOSTAT=ios) line

        if (ios < 0) then
            write(formatted_msg, '(A, A)') "No angles found in data file: ", trim(data_file_name)
            call abort_run(trim(formatted_msg), ERROR_READ)
        end if

        if (TRIM(adjustl(line)) /= "") then
            ! If not blank, rewind to re-read as angle data
            backspace(infile, IOSTAT=ios)
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to rewind file after non-blank line in: ", trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if
        end if

        ! Loop over all angles
        do k = 1, box%num_angles
            read(infile, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, I0, A)') "Unexpected end of file at angle line ", k, " in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, I0, A)') "I/O error reading angle line ", k, " in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            ! Parse angle info: id, type, atom1, atom2
            read(line, *, IOSTAT=ios) tmp1_int, tmp2_int, tmp3_int, tmp4_int, tmp5_int
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to parse angle line: '" // TRIM(line) // "' in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_PARSE)
            end if

            ! Store bond properties
            angle_ids_1d(k) = tmp1_int
            angle_types_1d(k) = tmp2_int
            angle_atoms_1d(k,1) = tmp3_int
            angle_atoms_1d(k,2) = tmp4_int
            angle_atoms_1d(k,3) = tmp5_int
            angle_found = angle_found + 1
        end do

    end subroutine read_lammps_angles

    !-----------------------------------------------------------------------------
    ! Reads dihedral definitions from the "Dihedrals" section of a LAMMPS data file
    ! and stores dihedral IDs, types, and connected atom indices in the corresponding arrays.
    ! If no dihedrals are expected, the routine returns immediately.
    !-----------------------------------------------------------------------------
    subroutine read_lammps_dihedrals(infile, data_file_name, box)

        ! Input
        integer, intent(in) :: infile                    ! Fortran unit number of the LAMMPS data file
        character(len=*), intent(in) :: data_file_name   ! Path to the LAMMPS data file
        type(type_box), intent(inout) :: box             ! Box structure containing dihedral counts and arrays

        ! Local variables
        character(len=200) :: line, formatted_msg        ! Current line and formatted error messages
        integer :: ios, k                                 ! I/O status and loop index
        integer :: dihedral_found                          ! Counter for dihedrals read
        integer :: tmp1_int, tmp2_int, tmp3_int, tmp4_int, tmp5_int, tmp6_int ! Temporary variables for dihedral id, type, atoms
        integer, parameter :: ERROR_READ = 24             ! Error code for I/O errors
        integer, parameter :: ERROR_PARSE = 25            ! Error code for parsing errors

        dihedral_found = 0

        ! If no dihedral are expected, skip the whole routine
        if (box%num_dihedrals == 0) then
            write(formatted_msg, '(A, A)') "No dihedrals expected in data file: ", trim(data_file_name)
            call InfoMessage(trim(formatted_msg))
            return
        end if

        ! Read until "Dihedrals" section or end of file
        do
            read(infile, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, A)') "No dihedrals found in data file: ", trim(data_file_name)
                call InfoMessage(trim(formatted_msg))
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, A)') "I/O error reading data file: ", trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            ! Check for "Dihedrals" or "Dihedrals # full" (case-insensitive)
            if (index(adjustl(line), 'Dihedrals') == 1) then
                exit
            end if
        end do

        ! Read and discard a blank line, if present
        read(infile, '(A)', IOSTAT=ios) line

        if (ios < 0) then
            write(formatted_msg, '(A, A)') "No dihedrals found in data file: ", trim(data_file_name)
            call abort_run(trim(formatted_msg), ERROR_READ)
        end if

        if (TRIM(adjustl(line)) /= "") then
            ! If not blank, rewind to re-read as bond data
            backspace(infile, IOSTAT=ios)
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to rewind file after non-blank line in: ", trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if
        end if

        ! Loop over all dihedrals
        do k = 1, box%num_dihedrals
            read(infile, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, I0, A)') "Unexpected end of file at dihedrals line ", &
                    k, " in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, I0, A)') "I/O error reading dihedrals line ", &
                    k, " in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            ! Parse dihedrals info: id, type, atom1, atom2
            read(line, *, IOSTAT=ios) tmp1_int, tmp2_int, tmp3_int, tmp4_int, tmp5_int, tmp6_int
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to parse dihedrals line: '" // &
                    TRIM(line) // "' in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_PARSE)
            end if

            ! Store bond properties
            dihedral_ids_1d(k) = tmp1_int
            dihedral_types_1d(k) = tmp2_int
            dihedral_atoms_1d(k,1) = tmp3_int
            dihedral_atoms_1d(k,2) = tmp4_int
            dihedral_atoms_1d(k,3) = tmp5_int
            dihedral_atoms_1d(k,4) = tmp6_int
            dihedral_found = dihedral_found + 1
        end do

    end subroutine read_lammps_dihedrals

    !-----------------------------------------------------------------------------
    ! Reads improper dihedral definitions from the "Impropers" section of a LAMMPS data file
    ! and stores improper IDs, types, and connected atom indices in the corresponding arrays.
    ! If no impropers are expected, the routine returns immediately.
    !-----------------------------------------------------------------------------
    subroutine read_lammps_impropers(infile, data_file_name, box)

        ! Input
        integer, intent(in) :: infile                    ! Fortran unit number of the LAMMPS data file
        character(len=*), intent(in) :: data_file_name   ! Path to the LAMMPS data file
        type(type_box), intent(inout) :: box             ! Box structure containing improper counts and arrays

        ! Local variables
        character(len=200) :: line, formatted_msg        ! Current line and formatted error messages
        integer :: ios, k                                 ! I/O status and loop index
        integer :: improper_found                          ! Counter for impropers read
        integer :: tmp1_int, tmp2_int, tmp3_int, tmp4_int, tmp5_int, tmp6_int ! Temporary variables for improper id, type, atoms
        integer, parameter :: ERROR_READ = 24             ! Error code for I/O errors
        integer, parameter :: ERROR_PARSE = 25            ! Error code for parsing errors

        improper_found = 0

        ! If no improper are expected, skip the whole routine
        if (box%num_impropers == 0) then
            write(formatted_msg, '(A, A)') "No impropers expected in data file: ", trim(data_file_name)
            call InfoMessage(trim(formatted_msg))
            return
        end if

        ! Read until "Impropers" section or end of file
        do
            read(infile, '(A)', IOSTAT=ios) line
            if (ios < 0) then
                write(formatted_msg, '(A, A)') "No impropers found in data file: ", trim(data_file_name)
                call InfoMessage(trim(formatted_msg))
            end if
            if (ios > 0) then
                write(formatted_msg, '(A, A)') "I/O error reading data file: ", trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            ! Check for "Impropers" or "Impropers # full" (case-insensitive)
            if (index(adjustl(line), 'Impropers') == 1) then
                exit
            end if
        end do

        ! Read and discard a blank line, if present
        read(infile, '(A)', IOSTAT=ios) line

        if (ios < 0) then
            write(formatted_msg, '(A, A)') "No impropers found in data file: ", trim(data_file_name)
            call abort_run(trim(formatted_msg), ERROR_READ)
        end if

        if (TRIM(adjustl(line)) /= "") then
            ! If not blank, rewind to re-read as bond data
            backspace(infile, IOSTAT=ios)
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to rewind file after non-blank line in: ", trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if
        end if

        ! Loop over all impropers
        do k = 1, box%num_impropers
            read(infile, '(A)', IOSTAT=ios) line

            if (ios < 0) then
                write(formatted_msg, '(A, I0, A)') &
                    "Unexpected end of file at impropers line ", k, " in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            if (ios > 0) then
                write(formatted_msg, '(A, I0, A)') "I/O error reading impropers line ", &
                    k, " in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_READ)
            end if

            ! Parse impropers info: id, type, atom1, atom2
            read(line, *, IOSTAT=ios) tmp1_int, tmp2_int, tmp3_int, tmp4_int, tmp5_int, tmp6_int
            if (ios /= 0) then
                write(formatted_msg, '(A, A)') "Failed to parse impropers line: '" &
                    // TRIM(line) // "' in: " // trim(data_file_name)
                call abort_run(trim(formatted_msg), ERROR_PARSE)
            end if

            ! Store bond properties
            improper_ids_1d(k) = tmp1_int
            improper_types_1d(k) = tmp2_int
            improper_atoms_1d(k,1) = tmp3_int
            improper_atoms_1d(k,2) = tmp4_int
            improper_atoms_1d(k,3) = tmp5_int
            improper_atoms_1d(k,4) = tmp6_int
            improper_found = improper_found + 1
        end do

    end subroutine read_lammps_impropers

    !-----------------------------------------------------------------------------
    ! Sorts atoms in the box structure by their original IDs (atom_original_1d)
    ! using an insertion sort algorithm. Updates all related arrays including
    ! atom_names_1d, atom_charges_1d, atom_xyz, atom_types_1d, atom_ids_1d.
    ! Calls detect_residue_pattern() after sorting to update residue mapping.
    !-----------------------------------------------------------------------------
    subroutine sort_atoms_by_original_ID(box)

        ! Input parameters
        type(type_box), intent(inout) :: box           ! Box structure containing atom data arrays

        ! Local variables
        integer :: i, j, temp_index, tmp_orig            ! Loop indices and temporary value
        integer, allocatable :: sort_index(:)           ! Sorting index array
        character(10), allocatable :: tmp_atom_names_1d(:) ! Temporary atom names array
        real(real64), allocatable :: tmp_atom_charges_1d(:) ! Temporary atom charges array
        real(real64), allocatable :: tmp_atom_xyz(:, :)   ! Temporary atom coordinates array
        integer, allocatable :: tmp_atom_types_1d(:)    ! Temporary atom types array
        integer, allocatable :: tmp_atom_original_1d(:) ! Temporary original IDs array

        ! Allocate sorting index array
        allocate(sort_index(box%num_atoms))
        do i = 1, box%num_atoms
            sort_index(i) = i
        end do

        ! Insertion sort based on atom_original_1d
        do i = 2, box%num_atoms
            temp_index = sort_index(i)
            tmp_orig = atom_original_1d(temp_index)
            j = i - 1

            do while (j >= 1)
                if (atom_original_1d(sort_index(j)) > tmp_orig) then
                    sort_index(j + 1) = sort_index(j)
                    j = j - 1
                else
                    exit
                end if
            end do

            sort_index(j + 1) = temp_index
        end do

        ! Allocate temporary arrays
        allocate(tmp_atom_names_1d(box%num_atoms))
        allocate(tmp_atom_charges_1d(box%num_atoms))
        allocate(tmp_atom_xyz(3, box%num_atoms))
        allocate(tmp_atom_types_1d(box%num_atoms))
        allocate(tmp_atom_original_1d(box%num_atoms))

        ! Apply sorted indices
        do i = 1, box%num_atoms
            tmp_atom_names_1d(i)   = atom_names_1d(sort_index(i))
            tmp_atom_charges_1d(i) = atom_charges_1d(sort_index(i))
            tmp_atom_xyz(:, i) = atom_xyz(:, sort_index(i))
            tmp_atom_types_1d(i)   = atom_types_1d(sort_index(i))
            tmp_atom_original_1d(i) = atom_original_1d(sort_index(i))
        end do

        ! Update global arrays with sorted data
        do i = 1, box%num_atoms
            atom_names_1d(i) = tmp_atom_names_1d(i)
            atom_charges_1d(i) = tmp_atom_charges_1d(i)
            atom_xyz(:, i) = tmp_atom_xyz(:, i)
            atom_types_1d(i) = tmp_atom_types_1d(i)
            atom_ids_1d(i) = i
            atom_original_1d(i) = tmp_atom_original_1d(i)
        end do

        call detect_residue_pattern(box, tmp_atom_types_1d)

        if (allocated(tmp_atom_names_1d)) deallocate(tmp_atom_names_1d)
        if (allocated(tmp_atom_charges_1d)) deallocate(tmp_atom_charges_1d)
        if (allocated(tmp_atom_xyz)) deallocate(tmp_atom_xyz)
        if (allocated(tmp_atom_types_1d)) deallocate(tmp_atom_types_1d)
        if (allocated(tmp_atom_original_1d)) deallocate(tmp_atom_original_1d)

    end subroutine sort_atoms_by_original_ID

    !-----------------------------------------------------------------------------
    ! Maps atoms to their residue types based on atom types and updates the
    ! types_pattern array in the nb structure. Each residue type maintains a
    ! counter to track placement of atoms within the residue template.
    !-----------------------------------------------------------------------------
    subroutine detect_residue_pattern(box, tmp_atom_types_1d)

        type(type_box), intent(inout) :: box                  ! Box structure containing atom data
        integer, allocatable, intent(in) :: tmp_atom_types_1d(:) ! Array of atom types sorted by original ID

        integer :: idx, i, j, residue                         ! Loop indices and residue index
        logical :: found                                      ! Flag to indicate atom matched a residue
        integer, allocatable :: cpt_per_residue(:)            ! Counters for each residue type

        ! one counter per residue
        allocate(cpt_per_residue(nb%type_residue))
        cpt_per_residue = 1

        do idx = 1, box%num_atoms
            found = .false.
            residue = -1

            ! find which residue this atom belongs to
            do i = 1, nb%type_residue
                do j = 1, nb%types_per_residue(i)
                    if (res%types_2d(i,j) == tmp_atom_types_1d(idx)) then
                        found = .true.
                        residue = i
                        exit  ! break from j-loop
                    end if
                end do
                if (found) exit  ! break from i-loop
            end do

            if (found) then
                nb%types_pattern(residue, cpt_per_residue(residue)) = tmp_atom_types_1d(idx)
                cpt_per_residue(residue) = cpt_per_residue(residue) + 1
                if (cpt_per_residue(residue) > nb%atom_in_residue(residue)) then
                    cpt_per_residue(residue) = 1
                end if
            end if
        end do

    end subroutine detect_residue_pattern

    !-----------------------------------------------------------------------------
    ! Checks if a given atom name exists in a specific row of a 2D array of expected names.
    ! Returns .true. if found, .false. otherwise.
    !-----------------------------------------------------------------------------
    logical function name_in_expected(atom_name, n_atoms, names2d, row)

        character(len=*), intent(in) :: atom_name               ! Atom name to search for
        integer, intent(in) :: n_atoms, row                     ! Number of atoms in row, row index
        character(len=*), dimension(:,:), intent(in) :: names2d ! 2D array of expected names
        integer :: k                                            ! Loop index

        name_in_expected = .false.
        do k = 1, n_atoms
            if (trim(atom_name) == trim(names2d(row,k))) then

                name_in_expected = .true.
                return
            
            end if
        end do

    end function name_in_expected

    !-----------------------------------------------------------------------
    ! Scan through the list of atoms in the simulation box and group them
    ! into molecules (residues) according to predefined residue templates.
    !-----------------------------------------------------------------------
    subroutine detect_molecules(box)

        ! Input parameters
        type(type_box), intent(inout) :: box    ! Box structure containing atom data

        ! Local variables
        integer :: i, j, k, cpt                 ! loop counters
        integer :: detected_nb_max_molecule     ! track max molecules found
        character(len=200) :: formatted_msg     ! buffer for error/info messages

        detected_nb_max_molecule = 0

        ! Initialize arrays inside box
        box%atom_types = 0
        box%atom_ids = 0
        box%atom_names = ""
        box%atom_charges = 0.0D0

        ! Loop over residue types
        do i = 1, nb%type_residue

            box%num_residues(i) = 0
            k = 1

            ! Scan atoms until we reach the end of the box
            do
                if (k > box%num_atoms) exit ! stop when no atoms left

                ! Reset atom counter for matching residue pattern
                cpt = 1

                ! If atom at position k matches the first atom of residue i
                if (atom_types_1d(k) == nb%types_pattern(i, cpt)) then

                    ! Ensure there are enough atoms left to build a full residue
                    if (k + nb%atom_in_residue(i) - 1 > box%num_atoms) then
                        call abort_run("Not enough atoms left in box to complete residue type ")
                    end if

                    ! Copy all atoms belonging to this residue
                    do j = 1, nb%atom_in_residue(i)
                        box%atom_types(i, j) = atom_types_1d(k)
                        box%atom_ids(i, j) = atom_original_1d(k)
                        box%atom_names(i, j) = atom_names_1d(k)
                        box%atom_charges(i, j) = atom_charges_1d(k)

                        ! Check atom order against residue pattern
                        if (thermo%is_active(i)) then
                            if (atom_types_1d(k) /= nb%types_pattern(i, cpt)) then
                                call abort_run("Issue with atom order in data file")
                            end if
                        end if

                        ! Move to next atom and increment residue pattern index
                        k = k + 1
                        cpt = cpt + 1

                    end do

                    ! Increment residue counter for type i
                    box%num_residues(i) = box%num_residues(i) + 1

                else
                    k = k + 1
                end if
            end do

            ! Track maximum number of molecules detected among all types
            if (box%num_residues(i) > detected_nb_max_molecule) detected_nb_max_molecule = box%num_residues(i)

        end do

        ! Global check: maximum allowed molecules
        if (detected_nb_max_molecule > NB_MAX_MOLECULE) then
            write(formatted_msg, '(A, I0)') "The number of molecules exceeds the maximum allowed = ", NB_MAX_MOLECULE
            call InfoMessage(trim(formatted_msg))
            write(formatted_msg, '(A, I0)') "Number of molecules found = ", detected_nb_max_molecule
            call abort_run(trim(formatted_msg) // " Please increase 'NB_MAX_MOLECULE' in src/parameters.f90 and recompile.", 11)
        end if

    end subroutine detect_molecules

    !-----------------------------------------------------------------------
    ! Loop through all molecules in the box and, for those marked as active,
    ! extract their coordinates, apply a repair procedure, and write the
    ! updated coordinates back into the global arrays.
    !-----------------------------------------------------------------------
    subroutine repair_active_molecules(box)

        ! Input parameters
        type(type_box), intent(inout) :: box

        ! Local variables
        integer :: natoms, nmolecules             ! atoms per residue, molecules per type
        integer :: i, j, k, l, m, k_temp          ! loop counters and atom index tracking
        real(real64), allocatable :: tmp_xyz(:, :) ! temporary xyz coordinates
        real(real64) :: dist                      ! Intramoleclar distance for sanity check

        ! Allocate temporary arrays for one molecule's atom positions
        allocate(tmp_xyz(3, maxval(nb%atom_in_residue)))

        k = 1  ! Global atom index

        ! Loop over all residue types
        do i = 1, nb%type_residue

            natoms = nb%atom_in_residue(i)      ! number of atoms per molecule of type i
            nmolecules = box%num_residues(i)    ! number of molecules of this type

            ! Loop over each molecule of this residue type
            do j = 1, nmolecules

                if (thermo%is_active(i)) then ! ACTIVE molecule → repair coordinates

                    k_temp = k ! Save starting index for this molecule

                    ! Copy coordinates into temporary arrays
                    tmp_xyz(:, 1:natoms) = atom_xyz(:, k:k+natoms-1)
                    k = k + natoms

                    ! Repair the molecule if its marked as active
                    call repair_molecule(tmp_xyz, natoms, box)

                    k = k_temp  ! Reset k to overwrite the same atoms

                    ! Sanity check : make sure distance between atoms of a same
                    ! molecules are not too large, and are not overlapping
                    do l = 1, natoms-1
                        do m = l+1, natoms
                            
                            dist = sqrt(sum((tmp_xyz(:, l) - tmp_xyz(:, m))**2))
                            ! do to : do not hardcode these value

                            if (dist > 10.0_real64) then
                                call warn_user("Unusually large distance (> 1 nm) detected in active residue")
                            else if (dist < 1.0e-5_real64) then
                                call warn_user("Overlapping atoms detected in molecule")
                            end if

                        end do
                    end do

                    ! Update coordinates with repaired values
                    atom_xyz(:, k:k+natoms-1) = tmp_xyz(:, 1:natoms)
                    k = k + natoms

                ! For inactve residue, update counter
                else
                    k = k + natoms ! INACTIVE molecule → skip without modification
                end if
            end do
        end do

        ! Free temporary arrays
        deallocate(tmp_xyz)

    end subroutine repair_active_molecules

    !---------------------------------------------------------------
    ! Transform absolute atom coordinates into relative coordinates
    ! with respect to the molecule's center of mass (CoM).
    !---------------------------------------------------------------
    subroutine transform_to_COM_frame(box, is_reservoir)

        ! Input parameters
        type(type_box), intent(inout) :: box
        logical, intent(in) :: is_reservoir                             ! To indicate if reservoir

        ! Local variables
        integer :: i, j, k, l, m, cpt
        real(real64), dimension(3) :: com, original_com
        integer :: nb_res                                               ! Number of molecules of this type
        real(real64) :: total_mass
        type(type_coordinate), pointer :: coord                         ! Pointer for host or guest coordinate

        ! Temporary arrays (assumed declared module-wide or can be moved here)
        real(real64), allocatable :: tmp_atom_masses_1d(:)
        real(real64), allocatable :: tmp_atom_xyz(:,:)

        ! Allocate temporaries (size depends on max atoms per residue)
        allocate(tmp_atom_masses_1d(box%num_atoms))

        if (.not. allocated(res%masses_1d)) then
            allocate(res%masses_1d(primary%num_atomtypes))
        end if

        allocate(tmp_atom_xyz(3, box%num_atoms))

        cpt = 1
        do i = 1, nb%type_residue

            ! #tofix, should have been done before
            ! Store CoM position
            if (thermo%is_active(i)) then
                resid_location(i) = TYPE_GUEST
            else
                resid_location(i) = TYPE_HOST
            end if

            if (is_reservoir) then
                coord => gas
            else
                coord => get_coord(i)
            end if

            nb_res = 0

            ! Build mass vector for the residue type i
            do j = 1, nb%atom_in_residue(i)
                do l = 1, nb%types_per_residue(i)
                    if (res%types_2d(i, l) == primary%atom_types(i, j)) then
                        res%masses_1d(primary%atom_types(i, j)) = primary%atom_masses(i, l)
                        tmp_atom_masses_1d = primary%atom_masses(i, l)
                        cpt = cpt + 1
                    end if
                end do
            end do

            ! Sanity check
            if (any(tmp_atom_masses_1d(1:nb%atom_in_residue(i)) <= 0.0_real64)) then
                call warn_user("Zero or negative atomic mass detected in residue type")
            else if (sum(tmp_atom_masses_1d(1:nb%atom_in_residue(i))) < 1.0e-6_real64) then
                call warn_user("Total molecular mass nearly zero in residue")
            end if

            ! Loop over all atoms and group them into molecules
            k = 1
            l = 1
            do while (k <= box%num_atoms)
                if (atom_types_1d(k) == box%atom_types(i, 1)) then

                    nb_res = nb_res + 1

                    ! Extract atom coordinates for this molecule
                    tmp_atom_xyz(:, 1:nb%atom_in_residue(i)) = atom_xyz(:, k:k+nb%atom_in_residue(i)-1)
                    k = k + nb%atom_in_residue(i)

                    ! Compute Center of Mass
                    call compute_COM(tmp_atom_xyz, nb%atom_in_residue(i), tmp_atom_masses_1d, com)

                    total_mass = compute_mass(nb%atom_in_residue(i), tmp_atom_masses_1d)
                    res%mass(i) = total_mass 

                    ! Save original CoM before applying periodic boundary conditions
                    original_com = com ! Save the CoM

                    ! Ensure CoM is inside simulation box
                    call apply_PBC(com, box)

                    ! Perform sanity checks
                    call check_finite_vector(com, "Invalid (NaN/Inf) CoM detected in residue")
                    call check_inside_bounds(com, box%bounds(:,1), box%bounds(:,2), &
                        "Molecule COM outside simulation box")
                    if (thermo%is_active(i)) then
                        call check_com_distance(tmp_atom_xyz, nb%atom_in_residue(i), &
                            original_com, 10.0_real64, "CoM unusually far from all atoms in residue type")
                    end if

                    coord%residue_exists(i) = .true.
                    coord%com(:, i, l) = com(:)
                    do m = 1, nb%atom_in_residue(i)
                        coord%offset(:, i, l, m) = tmp_atom_xyz(:, m) - original_com
                    end do

                    l = l + 1 ! Move to next molecule slot
                else
                    k = k + 1 ! Not start of this residue, keep scanning
                end if
            end do
        end do

        if (allocated(tmp_atom_masses_1d)) deallocate(tmp_atom_masses_1d)

    end subroutine transform_to_COM_frame

    !---------------------------------------------------------------
    ! Check if a given atom ID belongs to a specific residue.
    !---------------------------------------------------------------
    logical function isInResidue(box, res, atom_id)

        integer, intent(in) :: res, atom_id
        type(type_box), intent(inout) :: box

        integer :: n

        isInResidue = .false.

        do n = 1, nb%atom_in_residue(res)

            if (box%atom_ids(res, n) == atom_id) then

                isInResidue = .true.
                return

            end if
        end do

    end function isInResidue

    !---------------------------------------------------------------
    !Get the local index of an atom within a residue.
    !---------------------------------------------------------------
    integer function atomIndexInResidue(box, res, atom_id)

        integer, intent(in) :: res, atom_id
        type(type_box), intent(inout) :: box

        integer :: n

        atomIndexInResidue = -1
        do n = 1, nb%atom_in_residue(res)
            if (box%atom_ids(res, n) == atom_id) then

                atomIndexInResidue = n
                return
            
            end if
        end do

    end function atomIndexInResidue

end module data_parser
