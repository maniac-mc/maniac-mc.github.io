module write_utils

    use parameters
    use geometry_utils
    use simulation_state
    use monte_carlo_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    private
    public :: update_output_files

contains

    !------------------------------------------------------------------------------
    ! Updates all relevant output files for the current Monte Carlo step.
    !------------------------------------------------------------------------------
    subroutine update_output_files(later_step)

        logical, intent(in) :: later_step       ! Flag to distinguish first vs. later steps

        ! Write LAMMPS trajectory (main and reservoir)
        call write_dump_lammpstrj(primary, "trajectory.lammpstrj", later_step)
        if (has_reservoir) call write_dump_lammpstrj(reservoir, "reservoir.lammpstrj", later_step)

        ! Write energies and move counts to data file
        call write_dat_info()

        ! Write LAMMPS data file
        call write_topology_data(primary)

    end subroutine update_output_files

    !------------------------------------------------------------------------------
    ! Write trajectory in lammpstrj format (see https://docs.lammps.org for details)
    !------------------------------------------------------------------------------
    subroutine write_dump_lammpstrj(box, lammpstrj_filename, append_mode)

        ! Input parameters
        type(type_box), intent(inout) :: box                            ! Type for box
        character(len=*), intent(in), optional :: lammpstrj_filename    ! Output trajectory in lammpstrj format
        logical, optional :: append_mode                                ! Flag to append to existing files

        ! Local variables
        integer :: i, j, k                      ! Loop indices over residue types, residues, and atoms
        integer :: atom_id, atom_type
        integer :: UNIT_LMP = 18
        logical :: do_append                    ! Internal flag controlling append/overwrite
        integer :: dim                          ! Integer for looping over dimensions
        character(len=20) :: file_position      ! 'APPEND' or 'ASIS' for file positioning
        real(real64), dimension(3) :: pos, com  ! Center of mass and atomic position vectors

        ! Determine whether to append or overwrite
        do_append = .false.
        if (PRESENT(append_mode)) do_append = append_mode
        if (do_append) then
            file_position = 'APPEND'
        else
            file_position = 'ASIS'
        end if

        ! Open file for writing
        open(UNIT=UNIT_LMP, FILE=trim(output_path) // lammpstrj_filename, &
            STATUS='UNKNOWN', ACTION='write', POSITION=file_position)

        ! Write LAMMPS-style header
        write(UNIT_LMP, '(A)') "ITEM: TIMESTEP"
        write(UNIT_LMP, '(I10)') input%nb_block
        write(UNIT_LMP, '(A)') "ITEM: NUMBER OF ATOMS"
        write(UNIT_LMP, '(I10)') box%num_atoms
        write(UNIT_LMP, '(A)') "ITEM: BOX BOUNDS pp pp pp"
        write(UNIT_LMP, '(F15.8,1X,F15.8)') -box%matrix(1, 1)/2, box%matrix(1, 1)/2
        write(UNIT_LMP, '(F15.8,1X,F15.8)') -box%matrix(2, 2)/2, box%matrix(2, 2)/2
        write(UNIT_LMP, '(F15.8,1X,F15.8)') -box%matrix(3, 3)/2, box%matrix(3, 3)/2

        ! Atom data header
        write(UNIT_LMP, '(A)') "ITEM: ATOMS id type x y z"

        atom_id = 0
        do i = 1, nb%type_residue
            do j = 1, box%num_residues(i)

                ! Extract CoM
                do dim = 1, 3
                    com(dim) = box%mol_com(dim, i, j)
                end do

                ! Wrap CoM into box for active molecules
                if (input%is_active(i) == 1) then
                    call wrap_into_box(com, box)
                end if

                do k = 1, nb%atom_in_residue(i)

                    atom_id = atom_id + 1
                    atom_type = box%atom_types(i, k)

                    do dim = 1, 3
                        pos(dim) = com(dim) + box%site_offset(dim, i, j, k)
                    end do

                    ! Wrap position for inactive structure
                    if (input%is_active(i) == 0) then
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
        ! For current_block == 0, we want to recreate the file from scratch
        ! For later blocks, we append to existing files
        if (current_block == 0) then
            file_status = 'REPLACE'
        else
            file_status = 'OLD'
        end if

        call write_dat_energy(current_block, file_status)

        call write_dat_number(current_block, file_status)

        call write_dat_mcmove(current_block, file_status)

        call write_dat_widom(current_block, file_status)

    end subroutine write_dat_info

    !------------------------------------------------------------------------------
    ! Write energy.dat. Outputs total energy, Coulomb, non-Coulomb,
    ! intramolecular, and Ewald contributions.
    !------------------------------------------------------------------------------
    subroutine write_dat_energy(current_block, file_status)

        ! Input arguments
        integer, intent(in) :: current_block
        character(len=*), intent(in) :: file_status

        ! Local variables
        character(len=512) :: line
        character(len=256) :: header
        integer :: UNIT_ENERGY = 18

        ! Open energy.dat in append mode
        open(unit=UNIT_ENERGY, file=trim(output_path) // 'energy.dat', &
            status=file_status, action='write', position='append')

        ! Write header only at first block
        if (current_block == 0) then
            header = '#    block        total        recipCoulomb' // &
                    '     non-coulomb      coulomb     ewald_self    intramolecular-coulomb'
            write(UNIT_ENERGY, '(A)') trim(header)
        end if

        ! Build data line with proper formatting
        write(line,'(I10,1X,F16.6,1X,F16.6,1X,F16.6,1X,F16.6,1X,F16.6,1X,F16.6)') &
            current_block, energy%total, energy%recip_coulomb, energy%non_coulomb, &
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
    subroutine write_dat_widom(current_block, file_status)

        ! Input arguments
        integer, intent(in) :: current_block
        character(len=*), intent(in) :: file_status

        ! Local variables
        integer :: type_residue
        character(len=256) :: filename
        character(len=64) :: line
        integer :: UNIT_WIDOM = 21

        ! Check if Widom calculation is active
        if (proba%widom > 0) then

            ! Compute excess and total chemical potentials
            call CalculateExcessMu()

            ! Loop over residue types
            do type_residue = 1, nb%type_residue

                if (input%is_active(type_residue) > 0) then

                    ! Construct the file name
                    filename = trim(output_path) // 'widom_' // trim(res%names_1d(type_residue)) // '.dat'

                    ! Open file in append mode
                    open(unit=UNIT_WIDOM, file=filename, status=file_status, action='write', position='append')

                    ! Write header only at the first block
                    if (current_block == 0) then
                        write(line,'(A10,1X,A16,1X,A16,1X,A16)') '# Block', 'Excess_Mu_kcalmol', 'Total_Mu_kcalmol', 'Widom_Samples'
                        write(UNIT_WIDOM,'(A)') trim(line)
                    end if

                    ! Write the data line: block, excess mu, total mu, number of samples
                    write(line,'(I10,1X,F16.6,1X,F16.6,1X,I12)') current_block, &
                        widom_stat%mu_ex(type_residue), widom_stat%mu_tot(type_residue), &
                        widom_stat%sample(type_residue)
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
    subroutine write_dat_number(current_block, file_status)

        ! Input arguments
        integer, intent(in) :: current_block
        character(len=*), intent(in) :: file_status

        ! Local variables
        integer :: resi
        character(len=256) :: filename
        character(len=32) :: line
        integer :: UNIT_COUNT  = 19

        ! Loop over residues
        do resi = 1, nb%type_residue

            if (input%is_active(resi) == 1) then
                ! Construct the filename for this residue
                filename = trim(output_path) // 'number_' // trim(res%names_1d(resi)) // '.dat'

                ! Open file for append (create if not exists)
                open(unit=UNIT_COUNT, file=filename, status=file_status, action='write', position='append')

                ! Write header only once (when current_block == 0)
                if (current_block == 0) then
                    write(line,'(A10,1X,A10)') '# Block', 'Active_Molecules'
                    write(UNIT_COUNT,'(A)') trim(line)
                end if

                ! Write the block number and count for this residue
                write(line,'(I10,1X,I10)') current_block, primary%num_residues(resi)
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
    subroutine write_dat_mcmove(current_block, file_status)

        ! Input arguments
        integer, intent(in) :: current_block
        character(len=*), intent(in) :: file_status

        ! Local variables
        character(len=512) :: header, line, tmp
        logical :: first_block
        integer :: UNIT_MOVES  = 20

        ! Determine if this is the first block written
        first_block = (current_block == 0)

        ! Open file
        open(unit=UNIT_MOVES, file=trim(output_path)//'moves.dat', &
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
        write(line,'(I12)') current_block

        ! Translation
        if (proba%translation > 0) then
            write(tmp,'(1X,I12,1X,I12)') counter%translations, counter%trial_translations
            line = trim(line)//tmp
        end if

        ! Rotation
        if (proba%rotation > 0) then
            write(tmp,'(1X,I12,1X,I12)') counter%rotations, counter%trial_rotations
            line = trim(line)//tmp
        end if

        ! Insertion / Deletion
        if (proba%insertion_deletion > 0) then
            write(tmp,'(1X,I12,1X,I12)') counter%creations, counter%trial_creations
            line = trim(line)//tmp
            write(tmp,'(1X,I12,1X,I12)') counter%deletions, counter%trial_deletions
            line = trim(line)//tmp
        end if

        ! Swap
        if (proba%swap > 0) then
            write(tmp,'(1X,I12,1X,I12)') counter%swaps, counter%trial_swaps
            line = trim(line)//tmp
        end if

        ! Widom
        if (proba%widom > 0) then
            write(tmp,'(1X,I12)') counter%trial_widom
            line = trim(line)//tmp
        end if

        ! Write the completed line
        write(UNIT_MOVES,'(A)') trim(line)

        close(UNIT_MOVES)

    end subroutine write_dat_mcmove

    !------------------------------------------------------------------------------
    ! Write topology in LAMMPS data format (see https://docs.lammps.org for details)
    !------------------------------------------------------------------------------
    subroutine write_topology_data(box)

        implicit none

        ! Input parameters
        type(type_box), intent(inout) :: box

        ! Local variables
        integer :: i, j, k, atom_id, mol_id, atom_type
        integer :: cpt_bond, cpt_atom, cpt_angle, cpt_dihedral, cpt_improper
        integer :: UNIT_DATA = 19
        real(real64) :: charge
        real(real64), dimension(3) :: pos
        integer :: dim ! Integer for looping over dimensions

        ! Update bond number count
        cpt_bond = 0
        do i = 1, nb%type_residue
            do j = 1, box%num_residues(i)
                do k = 1, nb%bonds_per_residue(i)
                    cpt_bond = cpt_bond + 1
                end do
            end do
        end do
        box%num_bonds = cpt_bond

        ! Update angle number count
        cpt_angle = 0
        do i = 1, nb%type_residue
            do j = 1, box%num_residues(i)
                do k = 1, nb%angles_per_residue(i)
                    cpt_angle = cpt_angle + 1
                end do
            end do
        end do
        box%num_angles = cpt_angle

        ! Update dihedral number count
        cpt_dihedral = 0
        do i = 1, nb%type_residue
            do j = 1, box%num_residues(i)
                do k = 1, nb%dihedrals_per_residue(i)
                    cpt_dihedral = cpt_dihedral + 1
                end do
            end do
        end do
        box%num_dihedrals = cpt_dihedral

        ! Update improper number count
        cpt_improper = 0
        do i = 1, nb%type_residue
            do j = 1, box%num_residues(i)
                do k = 1, nb%impropers_per_residue(i)
                    cpt_improper = cpt_improper + 1
                end do
            end do
        end do
        box%num_impropers = cpt_improper

        ! Open file
        open(UNIT=UNIT_DATA, FILE=trim(output_path) // data_filename, STATUS='REPLACE', ACTION='write')

        write(UNIT_DATA, *) "! LAMMPS data file (atom_style full)"
        write(UNIT_DATA, *) box%num_atoms, " atoms"
        write(UNIT_DATA, *) box%num_atomtypes, " atom types"
        write(UNIT_DATA, *) box%num_bonds, " bonds"
        write(UNIT_DATA, *) box%num_bondtypes, " bond types"
        write(UNIT_DATA, *) box%num_angles, " angles"
        write(UNIT_DATA, *) box%num_angletypes, " angle types"
        write(UNIT_DATA, *) box%num_dihedrals, " dihedrals"
        write(UNIT_DATA, *) box%num_dihedraltypes, " dihedral types"
        write(UNIT_DATA, *) box%num_impropers, " impropers"
        write(UNIT_DATA, *) box%num_impropertypes, " improper types"
        write(UNIT_DATA, *)

        ! X bounds
        write(UNIT_DATA, '(2(F15.8,1X))', ADVANCE='NO') box%bounds(1,1), box%bounds(1,2)
        write(UNIT_DATA, '(A)') "xlo xhi"

        ! Y bounds
        write(UNIT_DATA, '(2(F15.8,1X))', ADVANCE='NO') box%bounds(2,1), box%bounds(2,2)
        write(UNIT_DATA, '(A)') "ylo yhi"

        ! Z bounds
        write(UNIT_DATA, '(2(F15.8,1X))', ADVANCE='NO') box%bounds(3,1), box%bounds(3,2)
        write(UNIT_DATA, '(A)') "zlo zhi"

        if (box%is_triclinic) then
            write(UNIT_DATA, '(3(F15.8,1X))') box%tilt(1), box%tilt(2), box%tilt(3)
            write(UNIT_DATA, '(A)') "xy xz yz"
        end if

        write(UNIT_DATA, *)

        ! Masses section (assumes atomic mass array `atom_masses`)
        write(UNIT_DATA, *) "Masses"
        write(UNIT_DATA, *)
        do atom_id = 1, primary%num_atomtypes
            write(UNIT_DATA, '(I5, 1X, F12.6)') atom_id, box%site_masses_vector(atom_id)
        end do

        write(UNIT_DATA, *)

        ! Atoms section
        write(UNIT_DATA, *) "Atoms"
        write(UNIT_DATA, *)

        atom_id = 0
        mol_id = 0
        do i = 1, nb%type_residue
            do j = 1, box%num_residues(i)
                mol_id = mol_id + 1
                do k = 1, nb%atom_in_residue(i)

                    atom_id = atom_id + 1
                    atom_type = box%atom_types(i,k)
                    charge = box%atom_charges(i,k)

                    do dim = 1, 3
                        pos(dim) = box%mol_com(dim, i, j) + box%site_offset(dim, i, j, k)
                    end do

                    ! Only wrap atom of inative species (i.e. leave active molecules continuous at the pbc)
                    if (input%is_active(i) == 0) then
                        call wrap_into_box(pos, box)
                    end if

                    write(UNIT_DATA, '(I6,1X,I6,1X,I4,1X,F12.8,3(1X,F12.7))') &
                        atom_id, mol_id, atom_type, charge, pos(1), pos(2), pos(3)

                end do
            end do
        end do

        if (primary%num_bonds > 0 .or. reservoir%num_bonds > 0) then
            ! Atoms section
            write(UNIT_DATA, *)
            write(UNIT_DATA, *) "Bonds"
            write(UNIT_DATA, *)
            cpt_bond = 1
            cpt_atom = 0
            do i = 1, nb%type_residue
                do j = 1, box%num_residues(i)
                    do k = 1, nb%bonds_per_residue(i)
                        write(UNIT_DATA, *) cpt_bond, res%bond_type_2d(i, k, 1), &
                            cpt_atom + res%bond_type_2d(i, k, 2), &
                            cpt_atom + res%bond_type_2d(i, k, 3)
                        cpt_bond = cpt_bond + 1
                    end do
                    cpt_atom = cpt_atom + nb%atom_in_residue(i)
                end do
            end do
        end if

        if (primary%num_angles > 0 .or. reservoir%num_angles > 0) then
            ! Atoms section
            write(UNIT_DATA, *)
            write(UNIT_DATA, *) "Angles"
            write(UNIT_DATA, *)
            cpt_angle = 1
            cpt_atom = 0
            do i = 1, nb%type_residue
                do j = 1, box%num_residues(i)
                    do k = 1, nb%angles_per_residue(i)
                        write(UNIT_DATA, *) cpt_angle, res%angle_type_2d(i, k, 1), &
                            cpt_atom + res%angle_type_2d(i, k, 2), &
                            cpt_atom + res%angle_type_2d(i, k, 3), &
                            cpt_atom + res%angle_type_2d(i, k, 4)
                        cpt_angle = cpt_angle + 1
                    end do
                    cpt_atom = cpt_atom + nb%atom_in_residue(i)
                end do
            end do
        end if

        if (primary%num_dihedrals > 0 .or. reservoir%num_dihedrals > 0) then
            ! Atoms section
            write(UNIT_DATA, *)
            write(UNIT_DATA, *) "Dihedrals"
            write(UNIT_DATA, *)
            cpt_dihedral = 1
            cpt_atom = 0
            do i = 1, nb%type_residue
                do j = 1, box%num_residues(i)
                    do k = 1, nb%dihedrals_per_residue(i)
                        write(UNIT_DATA, *) cpt_dihedral, res%dihedral_type_2d(i, k, 1), &
                            cpt_atom + res%dihedral_type_2d(i, k, 2), &
                            cpt_atom + res%dihedral_type_2d(i, k, 3), &
                            cpt_atom + res%dihedral_type_2d(i, k, 4), &
                            cpt_atom + res%dihedral_type_2d(i, k, 5)
                        cpt_dihedral= cpt_dihedral + 1
                    end do
                    cpt_atom = cpt_atom + nb%atom_in_residue(i)
                end do
            end do
        end if

        if (primary%num_impropers > 0 .or. reservoir%num_impropers > 0) then
            ! Atoms section
            write(UNIT_DATA, *)
            write(UNIT_DATA, *) "Impropers"
            write(UNIT_DATA, *)
            cpt_improper = 1
            cpt_atom = 0
            do i = 1, nb%type_residue
                do j = 1, box%num_residues(i)
                    do k = 1, nb%impropers_per_residue(i)
                        write(UNIT_DATA, *) cpt_improper, res%improper_type_2d(i, k, 1), &
                            cpt_atom + res%improper_type_2d(i, k, 2), &
                            cpt_atom + res%improper_type_2d(i, k, 3), &
                            cpt_atom + res%improper_type_2d(i, k, 4), &
                            cpt_atom + res%improper_type_2d(i, k, 5)
                        cpt_improper= cpt_improper + 1
                    end do
                    cpt_atom = cpt_atom + nb%atom_in_residue(i)
                end do
            end do
        end if

        close(UNIT_DATA)

    end subroutine write_topology_data

end module write_utils
