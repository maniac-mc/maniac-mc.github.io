module geometry_utils

    use simulation_state
    use output_utils
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

contains

    !-----------------------------------------------------------
    ! Initializes a simulation box for molecular simulations:
    !  1. Determine box symmetry (cubic, orthorhombic, triclinic)
    !  2. Compute geometric properties: cell lengths, angles, volume
    !  3. Calculate the inverse box matrix for reciprocal-space operations
    !-----------------------------------------------------------
    subroutine prepare_simulation_box(box)

        ! Input argument
        type(type_box), intent(inout) :: box

        ! Step 1: Determine box symmetry
        call determine_box_symmetry(box)

        ! Step 2: Compute geometric properties
        call compute_cell_properties(box)

        ! Step 3: Compute inverse matrix for reciprocal-space operations
        call compute_box_determinant_and_inverse(box)

        ! Log the information
        call log_geometry_parameters(box)

    end subroutine prepare_simulation_box

    !-----------------------------------------------------------    !
    ! Wraps a Cartesian position into the simulation box using
    ! periodic boundary conditions (PBC). Supports orthogonal
    ! and triclinic boxes.
    !
    ! Orthogonal: x' = lower_bound + mod(x - lower_bound, L)
    ! Triclinic : fractional coordinates s = H^{-1}·(r - r0),
    !             wrap s into [0,1), then r' = r0 + H·s
    !-----------------------------------------------------------
    subroutine apply_PBC(pos, box)

        ! Input parameters
        type(type_box), intent(in) :: box
        real(real64), dimension(3), intent(inout) :: pos

        ! Local variables
        real(real64), dimension(3) :: frac_coords    ! Fractional coordinates inside the box
        real(real64) :: box_length                   ! Box length
        real(real64) :: lower_bound                  ! Low bound
        integer :: dim                               ! Loop index

        ! Orthogonal / Rectangular Box
        if (box%type <= 3) then

            do dim = 1,3

                box_length   = box%matrix(dim, dim)     ! Box length along this dimension
                lower_bound  = box%bounds(dim, 1)      ! Lower bound (minimum coordinate) along this dimension

                ! Wrap coordinate into [lower_bound, lower_bound + box_length):
                ! x' = lower_bound + mod(x - lower_bound, box_length)
                pos(dim) = lower_bound + modulo(pos(dim) - lower_bound, box_length)

            end do

        ! Triclinic Box
        else

            ! Convert Cartesian -> fractional coordinates
            ! Using the box matrix H (columns = box vectors):
            !   s = H^{-1} * (r - r0)
            ! where:
            !   r  = Cartesian position
            !   r0 = lower corner of the box (box%bounds(:,1))
            !   s  = fractional coordinates in [0,1)^3
            frac_coords = matmul(box%reciprocal, pos - box%bounds(:,1))

            ! Wrap fractional coords into [0,1)
            ! Equation: s' = mod(s, 1)
            ! Ensures each fractional component lies inside the unit cube
            do dim = 1,3
                frac_coords(dim) = modulo(frac_coords(dim), one)
            end do

            ! Convert back to Cartesian
            ! Equation: r' = r0 + H * s'
            ! Resulting position r' is wrapped inside the triclinic box
            pos = box%bounds(:,1) + matmul(box%matrix, frac_coords)

        end if

    end subroutine apply_PBC

    !-----------------------------------------------------------
    ! Ensures a position is inside the simulation box using the
    ! following conventions:
    !  - Cubic/Orthorhombic: [-L/2, L/2]
    !  - Triclinic (fractional space): [-0.5, 0.5)
    !-----------------------------------------------------------
    subroutine wrap_into_box(pos, box)

        ! Input parameters
        real(real64), dimension(3), intent(inout) :: pos   ! Cartesian position to wrap into box
        type(type_box), intent(inout) :: box               ! Simulation box (geometry + PBC matrices)

        ! Local variables
        real(real64), dimension(3) :: fractional_pos    ! Position in fractional coordinates (triclinic)
        integer :: dim                                  ! Loop index over Cartesian dimensions

        ! Cubic or Orthorhombic box
        if (box%type == 1 .or. box%type == 2) then

            do dim = 1, 3
                ! Shift position into [-L/2, L/2]
                pos(dim) = pos(dim) - box%matrix(dim, dim) * &
                        nint(pos(dim) / box%matrix(dim, dim))
            end do

        ! Triclinic box
        else if (box%type == 3) then

            ! Convert Cartesian to fractional coordinates
            fractional_pos = matmul(box%reciprocal, pos)

            ! Wrap fractional coordinates into [-0.5, 0.5)
            do dim = 1, 3
                fractional_pos(dim) = fractional_pos(dim) - nint(fractional_pos(dim))
            end do

            ! Convert back to Cartesian
            pos = matmul(box%matrix, fractional_pos)

        end if

    end subroutine wrap_into_box

    !----------------------------------------------------------------------------
    ! Computes the inverse and determinant of a 3x3 box matrix.
    ! Uses the adjugate (cofactor transpose) method:
    !   H^{-1} = adj(H) / det(H)
    ! Also checks for near-zero determinant to avoid numerical issues.
    !----------------------------------------------------------------------------
    subroutine compute_box_determinant_and_inverse(box)

        ! Input parameters
        type(type_box), intent(inout) :: box        ! Box structure (matrix, determinant, inverse) to be read and updated

        ! Local variables
        real(real64), dimension(3,3) :: adjugate    ! Adjugate (cofactor transpose) matrix of box_mat
        real(real64) :: reciprocal                  ! Reciprocal of the determinant of box_mat
        real(real64), dimension(3) :: vec1, vec2    ! Temporary vectors for cross products
        integer :: i, j                             ! Loop indices for matrix element iteration

        ! Compute the adjugate (cofactor transpose) using cross products
        ! First column of adjugate
        vec1 = box%matrix(:, 2)
        vec2 = box%matrix(:, 3)
        adjugate(:,1) = compute_cross_product(vec1, vec2)

        ! Second column
        vec1 = box%matrix(:, 3)
        vec2 = box%matrix(:, 1)
        adjugate(:,2) = compute_cross_product(vec1, vec2)

        ! Third column
        vec1 = box%matrix(:, 1)
        vec2 = box%matrix(:, 2)
        adjugate(:,3) = compute_cross_product(vec1, vec2)

        ! Compute determinant from first row
        box%determinant = dot_product(box%matrix(:,1), adjugate(:,1))

        ! Check for denormal/underflow values in determinant
        if (abs(box%determinant) < one) then
            call abort_run("Error: Determinant fell into denormal/underflow range")
        end if

        ! Calculate reciprocal of determinant if non-zero
        reciprocal = zero
        if (abs(box%determinant) > error) then  ! small tolerance to avoid near-zero
            reciprocal = one / box%determinant
        else
            ! Handle singular matrix case if needed, e.g. print error or return
            call WarnUser("Determinant is zero or near zero, inverse not computed.")
            return
        end if

        ! Store inverse matrix in reciprocal
        do i = 1, 3
            do j = 1, 3
                box%reciprocal(i, j) = reciprocal * adjugate(i, j)
            end do
        end do

    end subroutine compute_box_determinant_and_inverse

    !---------------------------------------------------------------------
    ! Compute the minimum-image distance between two atoms in a periodic box
    !
    ! For cubic and orthorhombic boxes, the minimum-image convention is
    ! applied per Cartesian component. For triclinic boxes, all 27
    ! neighboring periodic images are checked to ensure the true
    ! shortest distance is returned.
    !---------------------------------------------------------------------
    function minimum_image_distance(box, residue_type_1, molecule_index_1, atom_index_1, &
                            residue_type_2, molecule_index_2, atom_index_2) result(distance)

        ! Input parameters
        type(type_box), intent(in) :: box                 ! simulation box (geometry + PBC matrices)
        integer, intent(in) :: residue_type_1             ! residue type ID for atom 1
        integer, intent(in) :: molecule_index_1           ! molecule index within that residue type
        integer, intent(in) :: atom_index_1               ! atom index within the molecule (atom 1)
        integer, intent(in) :: residue_type_2             ! residue type ID for atom 2
        integer, intent(in) :: molecule_index_2           ! molecule index within that residue type
        integer, intent(in) :: atom_index_2               ! atom index within the molecule (atom 2)

        ! Local variables
        integer :: shift_x, shift_y, shift_z   ! periodic image shifts in triclinic search
        integer :: dim                         ! Cartesian dimension index (1..3)
        real(real64) :: delta(3)               ! displacement vector between atoms (with PBC)
        real(real64) :: trial_delta(3)         ! displacement to a specific periodic image
        real(real64) :: distance               ! final minimum-image distance
        real(real64) :: min_dist2              ! smallest squared distance found
        real(real64) :: trial_dist2            ! squared distance for current image
        real(real64) :: pos1(3), pos2(3)       ! absolute Cartesian positions of the two atoms

        pos1 = box%mol_com(:,residue_type_1,molecule_index_1) + &
            box%site_offset(:,residue_type_1,molecule_index_1,atom_index_1)

        pos2 = box%mol_com(:,residue_type_2,molecule_index_2) + &
            box%site_offset(:,residue_type_2,molecule_index_2,atom_index_2)

        ! Compute Cartesian difference between the two atoms
        delta = pos2 - pos1

        ! Cubic or Orthorhombic box
        if (box%type == 1 .or. box%type == 2) then

            ! Apply PBC to delta vector directly
            do dim = 1,3
                delta(dim) = modulo(delta(dim) + half*box%matrix(dim,dim), &
                    box%matrix(dim,dim)) - half*box%matrix(dim,dim)
            end do

            ! Compute distance
            distance = vector_norm(delta)

        ! Triclinic box
        else if (box%type == 3) then

            min_dist2 = huge(one)

            do shift_x = -1,1
                do shift_y = -1,1
                    do shift_z = -1,1

                        trial_delta = delta + shift_x*box%matrix(:,1) + &
                                            shift_y*box%matrix(:,2) + &
                                            shift_z*box%matrix(:,3)

                        trial_dist2 = sum(trial_delta * trial_delta)

                        if (trial_dist2 < min_dist2) then
                            min_dist2 = trial_dist2
                        end if

                    end do
                end do
            end do

            distance = sqrt(min_dist2)

        end if

    end function minimum_image_distance

    !-----------------------------------------------------------
    ! Computes geometric properties of the simulation box:
    !  - Lengths of cell vectors
    !  - Cosines of angles between vectors
    !  - Cell volume
    !  - Perpendicular widths along each axis
    !-----------------------------------------------------------
    subroutine compute_cell_properties(box)

        ! Input parameters
        type(type_box), intent(inout) :: box

        ! Local variables
        real(real64), dimension(3) :: vec_lengths ! Lengths of the three cell vectors
        real(real64), dimension(3) :: axb ! Cross product of vector a and b
        real(real64), dimension(3) :: bxc ! Cross product of vector b and c
        real(real64), dimension(3) :: cxa ! Cross product of vector c and a
        real(real64) :: a(3), b(3), c(3) ! Local vectors for cross products

        ! Calculate lengths of cell vectors
        box%metrics(1:3) = sqrt(sum(box%matrix(:, 1:3)**2, dim=1))

        ! Copy cell vectors from the box matrix
        a = box%matrix(:, 1)
        b = box%matrix(:, 2)
        c = box%matrix(:, 3)

        ! Calculate vector lengths
        vec_lengths(1) = vector_norm(a)
        vec_lengths(2) = vector_norm(b)
        vec_lengths(3) = vector_norm(c)

        ! Calculate cosines of angles using dot products divided by product of lengths
        box%metrics(4) = dot_product(a, b) / (vec_lengths(1) * vec_lengths(2))
        box%metrics(5) = dot_product(a, c) / (vec_lengths(1) * vec_lengths(3))
        box%metrics(6) = dot_product(b, c) / (vec_lengths(2) * vec_lengths(3))

        ! Calculate cross products
        axb = compute_cross_product(a, b)
        bxc = compute_cross_product(b, c)
        cxa = compute_cross_product(c, a)

        ! Calculate volume of cell
        box%volume = abs(dot_product(box%matrix(:,1), bxc))

        ! Calculate cell perpendicular widths
        box%metrics(7) = box%volume / vector_norm(bxc)
        box%metrics(8) = box%volume / vector_norm(cxa)
        box%metrics(9) = box%volume / vector_norm(axb)

    end subroutine compute_cell_properties

    !-----------------------------------------------------------    !
    ! Determines the symmetry type of a simulation
    !-----------------------------------------------------------
    subroutine determine_box_symmetry(box)

        ! Input parameters
        type(type_box), intent(inout) :: box

        ! Local variables
        real(real64), dimension(6) :: off_diag

        ! Collect off-diagonal elements
        off_diag = [box%matrix(1,2), box%matrix(1,3), &
                    box%matrix(2,1), box%matrix(2,3), &
                    box%matrix(3,1), box%matrix(3,2)]

        ! Check for triclinic (any off-diagonal > tolerance)
        if (maxval(abs(off_diag)) > error) then
            box%type = 3
        ! Check for orthorhombic (unequal diagonals)
        else if (abs(box%matrix(1,1) - box%matrix(2,2)) > error .or. &
                 abs(box%matrix(1,1) - box%matrix(3,3)) > error) then
            box%type = 2
        ! Otherwise, cubic
        else
            box%type = 1
        end if

    end subroutine determine_box_symmetry

    !-----------------------------------------------------------
    ! Print geometry information
    !-----------------------------------------------------------
    subroutine log_geometry_parameters(box)

        ! Input argument
        type(type_box), intent(inout) :: box

        ! Local variable
        character(200) :: formatted_msg ! Buffer for formatted output messages

        call LogMessage("====== Simulation preparation ======")
        call LogMessage("")

        select case (box%type)
            case (1)

                call LogMessage("Box symmetry type: Cubic")
            
            case (2)
            
                call LogMessage("Box symmetry type: Orthorhombic")
            
            case (3)
            
                call LogMessage("Box symmetry type: Triclinic")
            
            case default
            
                write(formatted_msg, '(A, I0)') 'Box symmetry type determined: ', box%type
                call LogMessage(formatted_msg)
        
        end select

        write(formatted_msg, '(A, F20.4)') 'Cell volume (Å^3): ', box%volume
        call LogMessage(formatted_msg)

    end subroutine log_geometry_parameters

end module geometry_utils
