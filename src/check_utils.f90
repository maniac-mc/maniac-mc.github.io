module check_utils

    use output_utils
    use simulation_state
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

contains

    !------------------------------------------------------------------------------
    ! Verify the integrity of molecules in the simulation box.
    ! Performs two main checks:
    ! 1. The center-of-mass (COM) of each molecule is inside
    !    the simulation box boundaries.
    ! 2. Atom offsets from the COM do not exceed a physically
    !    reasonable threshold (10 Å) for active molecules.
    !------------------------------------------------------------------------------
    subroutine validate_molecule_geometry(box, is_reservoir)

        ! Input parameters
        type(type_box), intent(inout) :: box
        logical, intent(in) :: is_reservoir                             ! To indicate if reservoir

        ! Local variables
        integer :: res_type, mol_index, atom_index, dim
        real(real64), dimension(3) :: offset, com
        type(type_coordinate), pointer :: coord                         ! Pointer for host or guest coordinate
        logical :: warned_large_offset
        warned_large_offset = .false.

        do res_type = 1, nb%type_residue

            if (is_reservoir) then
                coord => gas
            else
                coord => get_coord(res_type)
            end if

            do mol_index = 1, box%num_residues(res_type)

                ! Read molecule center of mass
                com = coord%com(:, res_type, mol_index)

                ! Assert COM is inside the box
                do dim = 1,3
                    if (com(dim) < box%bounds(dim,1) .or. com(dim) > box%bounds(dim,2)) then
                        call warn_user("Error: Molecule COM outside simulation box!")
                    end if
                end do

                ! Check atom offsets
                do atom_index = 1, nb%atom_in_residue(res_type)
                    
                    ! Read molecule offsets
                    offset = coord%offset(:, res_type, mol_index, atom_index)
                    
                    if (any(abs(offset) > ten) .and. (thermo%is_active(res_type))) then
                        if (.not. warned_large_offset) then
                            call warn_user("One of the active molecules has an offset larger than 1 nanometer.")
                            warned_large_offset = .true.
                        end if
                    end if
                end do
            end do
        end do

    end subroutine validate_molecule_geometry

    !------------------------------------------------------------------------------
    ! Ensure consistency of atomic masses between the primary
    ! simulation system and an external reservoir.
    !------------------------------------------------------------------------------
    subroutine AssertMassConsistency()

        integer :: j, k
        character(len=128) :: msg

        if (status%reservoir_provided) then
            ! Check that system and reservoir atom masses are consistent
            do j = 1, nb%type_residue
                do k = 1, nb%types_per_residue(j)
                    if (abs(primary%atom_masses(j, k) - reservoir%atom_masses(j, k)) > error) then

                        ! Generic header warning
                        call warn_user("Reservoir and system mass don't match.")

                        ! Detailed mismatch info
                        write(msg, '(A,I3,A,I3,A,F10.5,A,F10.5)') &
                            "Mismatch at residue ", j, ", site ", k, &
                            ": system mass = ", primary%atom_masses(j, k), &
                            ", reservoir mass = ", reservoir%atom_masses(j, k)
                        call warn_user(trim(msg))

                        ! Extra guidance
                        call warn_user("Make sure atom types are consistent between the two files.")
                    end if
                end do
            end do
        end if

    end subroutine AssertMassConsistency

end module check_utils
