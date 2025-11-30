module write_utils

    use parameters
    use geometry_utils
    use simulation_state
    use monte_carlo_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !------------------------------------------------------------------------------
    ! Updates all relevant output files for the current Monte Carlo step.
    !------------------------------------------------------------------------------
    subroutine update_output_files(later_step)

        ! Input coefficient
        logical, intent(in) :: later_step       ! Flag to distinguish first vs. later steps

        ! Write LAMMPS trajectory (main and reservoir)
        call write_dump_lammpstrj(primary, "trajectory.lammpstrj", later_step, is_reservoir = .false.)
        if (status%reservoir_provided) then
            call write_dump_lammpstrj(reservoir, "reservoir.lammpstrj", later_step, is_reservoir = .true.)
        end if

        ! Write energies and move counts to data file
        call write_dat_info()

        ! Write LAMMPS data file
        call write_topology_data(primary, is_reservoir = .false.)

    end subroutine update_output_files

    !------------------------------------------------------------------------------
    ! Write trajectory in lammpstrj format (see https://docs.lammps.org for details)
    !------------------------------------------------------------------------------
    subroutine write_dump_lammpstrj(box, lammpstrj_filename, append_mode, is_reservoir)

        ! Input parameters
        type(type_box), intent(inout) :: box                            ! Type for box
        character(len=*), intent(in), optional :: lammpstrj_filename    ! Output trajectory in lammpstrj format
        logical, intent(in), optional :: append_mode                    ! Flag to append to existing files
        logical, intent(in) :: is_reservoir                             ! To indicate if reservoir

        ! Local variables
        integer :: res_type, mol_index, atom_index                      ! Loop indices over residue types, residues, and atoms
        integer :: atom_id, atom_type                                   ! Atom type and id
        integer :: UNIT_LMP = 18
        logical :: do_append                                            ! Internal flag controlling append/overwrite
        character(len=20) :: file_position                              ! 'APPEND' or 'ASIS' for file positioning
        real(real64), dimension(3) :: pos, com                          ! Center of mass and atomic position vectors
        type(type_coordinate), pointer :: coord                         ! Pointer for host or guest coordinate

        ! Determine whether to append or overwrite
        do_append = .false.
        if (PRESENT(append_mode)) do_append = append_mode
        if (do_append) then
            file_position = 'APPEND'
        else
            file_position = 'ASIS'
        end if

        ! Open file for writing
        open(UNIT=UNIT_LMP, FILE=trim(path%outputs) // lammpstrj_filename, &
            STATUS='UNKNOWN', ACTION='write', POSITION=file_position)

        ! Write LAMMPS-style header
        write(UNIT_LMP, '(A)') "ITEM: TIMESTEP"
        write(UNIT_LMP, '(I10)') status%desired_block
        write(UNIT_LMP, '(A)') "ITEM: NUMBER OF ATOMS"
        write(UNIT_LMP, '(I10)') box%num%atoms
        write(UNIT_LMP, '(A)') "ITEM: BOX BOUNDS pp pp pp"
        write(UNIT_LMP, '(F15.8,1X,F15.8)') -box%cell%matrix(1, 1)/2, box%cell%matrix(1, 1)/2
        write(UNIT_LMP, '(F15.8,1X,F15.8)') -box%cell%matrix(2, 2)/2, box%cell%matrix(2, 2)/2
        write(UNIT_LMP, '(F15.8,1X,F15.8)') -box%cell%matrix(3, 3)/2, box%cell%matrix(3, 3)/2

        ! Atom data header
        write(UNIT_LMP, '(A)') "ITEM: ATOMS id type x y z"

        atom_id = 0
        do res_type = 1, res%number

            if (is_reservoir) then
                coord => gas
            else
                coord => get_coord(res_type)
            end if

            do mol_index = 1, box%num%residues(res_type)

                ! Extract CoM
                com(:) = coord%com(:, res_type, mol_index)

                ! Wrap CoM into box for active molecules
                if (thermo%is_active(res_type)) then
                    call wrap_into_box(com, box)
                end if

                do atom_index = 1, res%atom(res_type)

                    atom_id = atom_id + 1
                    atom_type = box%atoms%types(res_type, atom_index)

                    pos(:) = com(:) + coord%offset(:, res_type, mol_index, atom_index)
                    
                    ! Wrap position for inactive structure
                    if (.not. thermo%is_active(res_type)) then
                        call wrap_into_box(pos, box)
                    end if

                    write(UNIT_LMP, '(I6,1X,I4,3(1X,F12.7))') atom_id, atom_type, pos(1), pos(2), pos(3)
                end do
            end do
        end do

        close(UNIT_LMP)

    end subroutine write_dump_lammpstrj

    !------------------------------------------------------------------------------
    ! Write all simulation output files to .dat file for the current block.
    !------------------------------------------------------------------------------
    subroutine write_dat_info()

        ! Local variables
        character(len=8) :: file_status
        integer :: type_residue

        ! Decide whether to create new files or append
        ! For nb_block == 0, we want to recreate the file from scratch
        ! For later blocks, we append to existing files
        if (status%block == 0) then
            file_status = 'REPLACE'
        else
            file_status = 'OLD'
        end if

        call write_dat_energy(status%block, file_status)

        call write_dat_number(status%block, file_status)

        call write_dat_mcmove(status%block, file_status)

        call write_dat_widom(status%block, file_status)

    end subroutine write_dat_info

    !------------------------------------------------------------------------------
    ! Write energy.dat. Outputs total energy, Coulomb, non-Coulomb,
    ! intramolecular, and Ewald contributions.
    !------------------------------------------------------------------------------
    subroutine write_dat_energy(nb_block, file_status)

        ! Input arguments
        integer, intent(in) :: nb_block
        character(len=*), intent(in) :: file_status

        ! Local variables
        character(len=512) :: line
        character(len=256) :: header
        integer :: UNIT_ENERGY = 18

        ! Open energy.dat in append mode
        open(unit=UNIT_ENERGY, file=trim(path%outputs) // 'energy.dat', &
            status=file_status, action='write', position='append')

        ! Write header only at first block
        if (nb_block == 0) then
            header = '#    block        total        recipCoulomb' // &
                    '     non-coulomb      coulomb     ewald_self    intramolecular-coulomb'
            write(UNIT_ENERGY, '(A)') trim(header)
        end if

        ! Build data line with proper formatting
        write(line,'(I10,1X,F16.6,1X,F16.6,1X,F16.6,1X,F16.6,1X,F16.6,1X,F16.6)') &
            nb_block, energy%total, energy%recip_coulomb, energy%non_coulomb, &
            energy%coulomb, energy%ewald_self, energy%intra_coulomb
        write(UNIT_ENERGY,'(A)') trim(line)

        ! Close file
        close(UNIT_ENERGY)

    end subroutine write_dat_energy

    !------------------------------------------------------------------------------
    ! Write widom_RESNAME.dat for active residues if Widom sampling is enabled.
    ! Outputs block number, excess chemical potential, total chemical potential,
    ! and number of Widom samples. Header written only for first block.
    !------------------------------------------------------------------------------
    subroutine write_dat_widom(nb_block, file_status)

        ! Input arguments
        integer, intent(in) :: nb_block
        character(len=*), intent(in) :: file_status

        ! Local variables
        integer :: type_residue
        character(len=256) :: filename
        character(len=64) :: line
        integer :: UNIT_WIDOM = 21

        ! Check if Widom calculation is active
        if (proba%widom > 0) then

            ! Compute excess and total chemical potentials
            call calculate_excess_mu()

            ! Loop over residue types
            do type_residue = 1, res%number

                if (thermo%is_active(type_residue)) then

                    ! Construct the file name
                    filename = trim(path%outputs) // 'widom_' // trim(res%names(type_residue)) // '.dat'

                    ! Open file in append mode
                    open(unit=UNIT_WIDOM, file=filename, status=file_status, action='write', position='append')

                    ! Write header only at the first block
                    if (nb_block == 0) then
                        write(line,'(A10,1X,A16,1X,A16,1X,A16)') '# Block', 'Excess_Mu_kcalmol', 'Total_Mu_kcalmol', 'Widom_Samples'
                        write(UNIT_WIDOM,'(A)') trim(line)
                    end if

                    ! Write the data line: block, excess mu, total mu, number of samples
                    write(line,'(I10,1X,F16.6,1X,F16.6,1X,I12)') nb_block, &
                        statistic%mu_ex(type_residue), statistic%mu_tot(type_residue), &
                        statistic%sample(type_residue)
                    write(UNIT_WIDOM,'(A)') trim(line)

                    ! Close file
                    close(UNIT_WIDOM)

                end if

            end do

        end if

    end subroutine write_dat_widom

    !------------------------------------------------------------------------------
    ! Write number_RESNAME.dat for each active residue. Tracks block number
    ! and number of active molecules. Header is written only for the first block.
    !------------------------------------------------------------------------------
    subroutine write_dat_number(nb_block, file_status)

        ! Input arguments
        integer, intent(in) :: nb_block
        character(len=*), intent(in) :: file_status

        ! Local variables
        integer :: resi
        character(len=256) :: filename
        character(len=32) :: line
        integer :: UNIT_COUNT  = 19

        ! Loop over residues
        do resi = 1, res%number

            if (thermo%is_active(resi)) then
                ! Construct the filename for this residue
                filename = trim(path%outputs) // 'number_' // trim(res%names(resi)) // '.dat'

                ! Open file for append (create if not exists)
                open(unit=UNIT_COUNT, file=filename, status=file_status, action='write', position='append')

                ! Write header only once (when nb_block == 0)
                if (nb_block == 0) then
                    write(line,'(A10,1X,A10)') '# Block', 'Active_Molecules'
                    write(UNIT_COUNT,'(A)') trim(line)
                end if

                ! Write the block number and count for this residue
                write(line,'(I10,1X,I10)') nb_block, primary%num%residues(resi)
                write(UNIT_COUNT,'(A)') trim(line)

                close(UNIT_COUNT)
            end if

        end do

    end subroutine write_dat_number

    !------------------------------------------------------------------------------
    ! Write moves_new.dat. Tracks move accuracies and trials for translation,
    ! rotation, creation/deletion, swap, and Widom moves. Header generated
    ! dynamically based on which move types are enabled.
    !------------------------------------------------------------------------------
    subroutine write_dat_mcmove(nb_block, file_status)

        ! Input arguments
        integer, intent(in) :: nb_block
        character(len=*), intent(in) :: file_status

        ! Local variables
        character(len=512) :: header, line, tmp
        logical :: first_block
        integer :: UNIT_MOVES  = 20

        ! Determine if this is the first block written
        first_block = (nb_block == 0)

        ! Open file
        open(unit=UNIT_MOVES, file=trim(path%outputs)//'moves.dat', &
            status=file_status, action='write', position='append')

        ! -----------------------
        ! Build header dynamically
        ! -----------------------
        if (first_block) then
            header = ""
            ! Block column matches I12 numeric
            write(tmp,'(A12)') "Block"
            header = trim(tmp)

            ! All other columns: 1 space + 12 characters to match I12 numeric
            if (proba%translation > 0) then
                write(tmp,'(1X,A12)') "Trans_Acc"
                header = trim(header)//tmp
                write(tmp,'(1X,A12)') "Trans_Trial"
                header = trim(header)//tmp
            end if

            if (proba%rotation > 0) then
                write(tmp,'(1X,A12)') "Rot_Acc"
                header = trim(header)//tmp
                write(tmp,'(1X,A12)') "Rot_Trial"
                header = trim(header)//tmp
            end if

            if (proba%insertion_deletion > 0) then
                write(tmp,'(1X,A12)') "Create_Acc"
                header = trim(header)//tmp
                write(tmp,'(1X,A12)') "Create_Trial"
                header = trim(header)//tmp
                write(tmp,'(1X,A12)') "Delete_Acc"
                header = trim(header)//tmp
                write(tmp,'(1X,A12)') "Delete_Trial"
                header = trim(header)//tmp
            end if

            if (proba%swap > 0) then
                write(tmp,'(1X,A12)') "Swap_Acc"
                header = trim(header)//tmp
                write(tmp,'(1X,A12)') "Swap_Trial"
                header = trim(header)//tmp
            end if

            if (proba%widom > 0) then
                write(tmp,'(1X,A12)') "Widom_Trial"
                header = trim(header)//tmp
            end if

            write(UNIT_MOVES,'(A)') trim(header)
        end if

        ! -----------------------
        ! Build numeric line
        ! -----------------------

        ! Block
        write(line,'(I12)') nb_block

        ! Translation
        if (proba%translation > 0) then
            write(tmp,'(1X,I12,1X,I12)') counter%translations(2), counter%translations(1)
            line = trim(line)//tmp
        end if

        ! Rotation
        if (proba%rotation > 0) then
            write(tmp,'(1X,I12,1X,I12)') counter%rotations(2), counter%rotations(1)
            line = trim(line)//tmp
        end if

        ! Insertion / Deletion
        if (proba%insertion_deletion > 0) then
            write(tmp,'(1X,I12,1X,I12)') counter%creations(2), counter%creations(1)
            line = trim(line)//tmp
            write(tmp,'(1X,I12,1X,I12)') counter%deletions(2), counter%deletions(1)
            line = trim(line)//tmp
        end if

        ! Swap
        if (proba%swap > 0) then
            write(tmp,'(1X,I12,1X,I12)') counter%swaps(2), counter%swaps(1)
            line = trim(line)//tmp
        end if

        ! Widom
        if (proba%widom > 0) then
            write(tmp,'(1X,I12,1X,I12)') counter%widom(2), counter%widom(1)
            line = trim(line)//tmp
        end if

        ! Write the completed line
        write(UNIT_MOVES,'(A)') trim(line)

        close(UNIT_MOVES)

    end subroutine write_dat_mcmove

    !------------------------------------------------------------------------------
    ! Write topology in LAMMPS data format (see https://docs.lammps.org for details)
    !------------------------------------------------------------------------------
    subroutine write_topology_data(box, is_reservoir)

        ! Input parameters
        type(type_box), intent(inout) :: box
        logical, intent(in) :: is_reservoir                             ! To indicate if reservoir

        ! Local variables
        integer :: i, j, k, atom_id, mol_id, atom_type
        integer :: cpt_bond, cpt_atom, cpt_angle, cpt_dihedral, cpt_improper
        integer :: unit_data = 19
        real(real64) :: charge
        real(real64), dimension(3) :: pos
        integer :: dim ! Integer for looping over dimensions
        type(type_coordinate), pointer :: coord                         ! Pointer for host or guest coordinate

        ! Update bond number count
        cpt_bond = 0
        do i = 1, res%number
            do j = 1, box%num%residues(i)
                do k = 1, res%bonds(i)
                    cpt_bond = cpt_bond + 1
                end do
            end do
        end do
        box%num%bonds = cpt_bond

        ! Update angle number count
        cpt_angle = 0
        do i = 1, res%number
            do j = 1, box%num%residues(i)
                do k = 1, res%angles(i)
                    cpt_angle = cpt_angle + 1
                end do
            end do
        end do
        box%num%angles = cpt_angle

        ! Update dihedral number count
        cpt_dihedral = 0
        do i = 1, res%number
            do j = 1, box%num%residues(i)
                do k = 1, res%dihedrals(i)
                    cpt_dihedral = cpt_dihedral + 1
                end do
            end do
        end do
        box%num%dihedrals = cpt_dihedral

        ! Update improper number count
        cpt_improper = 0
        do i = 1, res%number
            do j = 1, box%num%residues(i)
                do k = 1, res%impropers(i)
                    cpt_improper = cpt_improper + 1
                end do
            end do
        end do
        box%num%impropers = cpt_improper

        ! Open file
        open(UNIT=unit_data, FILE=trim(path%outputs) // data_filename, STATUS='REPLACE', ACTION='write')

        write(unit_data, *) "! LAMMPS data file (atom_style full)"
        write(unit_data, *) box%num%atoms, " atoms"
        write(unit_data, *) box%num%atomtypes, " atom types"
        write(unit_data, *) box%num%bonds, " bonds"
        write(unit_data, *) box%num%bondtypes, " bond types"
        write(unit_data, *) box%num%angles, " angles"
        write(unit_data, *) box%num%angletypes, " angle types"
        write(unit_data, *) box%num%dihedrals, " dihedrals"
        write(unit_data, *) box%num%dihedraltypes, " dihedral types"
        write(unit_data, *) box%num%impropers, " impropers"
        write(unit_data, *) box%num%impropertypes, " improper types"
        write(unit_data, *)

        ! X bounds
        write(unit_data, '(2(F15.8,1X))', ADVANCE='NO') box%cell%bounds(1,1), box%cell%bounds(1,2)
        write(unit_data, '(A)') "xlo xhi"

        ! Y bounds
        write(unit_data, '(2(F15.8,1X))', ADVANCE='NO') box%cell%bounds(2,1), box%cell%bounds(2,2)
        write(unit_data, '(A)') "ylo yhi"

        ! Z bounds
        write(unit_data, '(2(F15.8,1X))', ADVANCE='NO') box%cell%bounds(3,1), box%cell%bounds(3,2)
        write(unit_data, '(A)') "zlo zhi"

        if (box%cell%shape == TRICLINIC) then
            write(unit_data, '(3(F15.8,1X))') box%cell%tilt(1), box%cell%tilt(2), box%cell%tilt(3)
            write(unit_data, '(A)') "xy xz yz"
        end if

        write(unit_data, *)

        ! Masses section (assumes atomic mass array `atom_masses`)
        write(unit_data, *) "Masses"
        write(unit_data, *)
        do atom_id = 1, primary%num%atomtypes
            write(unit_data, '(I5, 1X, F12.6)') atom_id, box%atoms%masses_vec(atom_id)
        end do

        write(unit_data, *)

        ! Atoms section
        write(unit_data, *) "Atoms"
        write(unit_data, *)

        atom_id = 0
        mol_id = 0
        do i = 1, res%number

            if (is_reservoir) then
                coord => gas
            else
                coord => get_coord(i)
            end if

            do j = 1, box%num%residues(i)
                mol_id = mol_id + 1
                do k = 1, res%atom(i)

                    atom_id = atom_id + 1
                    atom_type = box%atoms%types(i,k)
                    charge = box%atoms%charges(i,k)

                    do dim = 1, 3
                        pos(dim) = coord%com(dim, i, j) + coord%offset(dim, i, j, k)
                    end do

                    ! Only wrap atom of inative species (i.e. leave active molecules continuous at the pbc)
                    if (.not. thermo%is_active(i)) then
                        call wrap_into_box(pos, box)
                    end if

                    write(unit_data, '(I6,1X,I6,1X,I4,1X,F12.8,3(1X,F12.7))') &
                        atom_id, mol_id, atom_type, charge, pos(1), pos(2), pos(3)

                end do
            end do
        end do

        if (primary%num%bonds > 0 .or. reservoir%num%bonds > 0) then
            ! Atoms section
            write(unit_data, *)
            write(unit_data, *) "Bonds"
            write(unit_data, *)
            cpt_bond = 1
            cpt_atom = 0
            do i = 1, res%number
                do j = 1, box%num%residues(i)
                    do k = 1, res%bonds(i)
                        write(unit_data, *) cpt_bond, connect%bonds(i, k, 1), &
                            cpt_atom + connect%bonds(i, k, 2), &
                            cpt_atom + connect%bonds(i, k, 3)
                        cpt_bond = cpt_bond + 1
                    end do
                    cpt_atom = cpt_atom + res%atom(i)
                end do
            end do
        end if

        if (primary%num%angles > 0 .or. reservoir%num%angles > 0) then
            ! Atoms section
            write(unit_data, *)
            write(unit_data, *) "Angles"
            write(unit_data, *)
            cpt_angle = 1
            cpt_atom = 0
            do i = 1, res%number
                do j = 1, box%num%residues(i)
                    do k = 1, res%angles(i)
                        write(unit_data, *) cpt_angle, connect%angles(i, k, 1), &
                            cpt_atom + connect%angles(i, k, 2), &
                            cpt_atom + connect%angles(i, k, 3), &
                            cpt_atom + connect%angles(i, k, 4)
                        cpt_angle = cpt_angle + 1
                    end do
                    cpt_atom = cpt_atom + res%atom(i)
                end do
            end do
        end if

        if (primary%num%dihedrals > 0 .or. reservoir%num%dihedrals > 0) then
            ! Atoms section
            write(unit_data, *)
            write(unit_data, *) "Dihedrals"
            write(unit_data, *)
            cpt_dihedral = 1
            cpt_atom = 0
            do i = 1, res%number
                do j = 1, box%num%residues(i)
                    do k = 1, res%dihedrals(i)
                        write(unit_data, *) cpt_dihedral, connect%dihedrals(i, k, 1), &
                            cpt_atom + connect%dihedrals(i, k, 2), &
                            cpt_atom + connect%dihedrals(i, k, 3), &
                            cpt_atom + connect%dihedrals(i, k, 4), &
                            cpt_atom + connect%dihedrals(i, k, 5)
                        cpt_dihedral= cpt_dihedral + 1
                    end do
                    cpt_atom = cpt_atom + res%atom(i)
                end do
            end do
        end if

        if (primary%num%impropers > 0 .or. reservoir%num%impropers > 0) then
            ! Atoms section
            write(unit_data, *)
            write(unit_data, *) "Impropers"
            write(unit_data, *)
            cpt_improper = 1
            cpt_atom = 0
            do i = 1, res%number
                do j = 1, box%num%residues(i)
                    do k = 1, res%impropers(i)
                        write(unit_data, *) cpt_improper, connect%impropers(i, k, 1), &
                            cpt_atom + connect%impropers(i, k, 2), &
                            cpt_atom + connect%impropers(i, k, 3), &
                            cpt_atom + connect%impropers(i, k, 4), &
                            cpt_atom + connect%impropers(i, k, 5)
                        cpt_improper= cpt_improper + 1
                    end do
                    cpt_atom = cpt_atom + res%atom(i)
                end do
            end do
        end if

        close(unit_data)

    end subroutine write_topology_data

end module write_utils
