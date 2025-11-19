module helper_utils

    use simulation_state
    use constants
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

contains

    !-------------------------------------------------------
    ! Function: repeat_char
    ! Creates a string of repeated characters
    !-------------------------------------------------------
    pure function repeat_char(ch, n) result(res)

        character(len=*), intent(in) :: ch
        integer, intent(in) :: n
        character(len=n) :: res
        integer :: i

        do i = 1, n
            res(i:i) = ch
        end do
    end function repeat_char

    !========================================================
    ! Function: RotationMatrix
    !
    ! Returns a 3x3 rotation matrix for a given axis (X=1, Y=2, Z=3)
    ! and rotation angle theta (radians).
    !
    ! Inputs:
    !   axis  - integer, rotation axis (1=X, 2=Y, 3=Z)
    !   theta - real(real64), rotation angle in radians
    !
    ! Output:
    !   rotation_matrix - real(real64) 3x3 rotation matrix
    !========================================================
    function RotationMatrix(axis, theta) result(rotation_matrix)

        implicit none

        integer, intent(in) :: axis               ! Rotation axis (1=X, 2=Y, 3=Z)
        real(real64), intent(in) :: theta         ! Rotation angle in radians
        real(real64) :: rotation_matrix(3,3)      ! 3x3 rotation matrix to be returned
        real(real64) :: cos_theta, sin_theta      ! Cosine and sine of theta

        ! Compute trigonometric values
        cos_theta = cos(theta)
        sin_theta = sin(theta)

        ! Initialize as identity
        rotation_matrix = zero
        rotation_matrix(1,1) = one
        rotation_matrix(2,2) = one
        rotation_matrix(3,3) = one

        ! Fill rotation matrix based on axis
        select case(axis)

        case(1) ! X-axis

            rotation_matrix(2,2) = cos_theta
            rotation_matrix(2,3) = -sin_theta
            rotation_matrix(3,2) = sin_theta
            rotation_matrix(3,3) = cos_theta

        case(2) ! Y-axis

            rotation_matrix(1,1) = cos_theta
            rotation_matrix(1,3) = sin_theta
            rotation_matrix(3,1) = -sin_theta
            rotation_matrix(3,3) = cos_theta

        case(3) ! Z-axis

            rotation_matrix(1,1) = cos_theta
            rotation_matrix(1,2) = -sin_theta
            rotation_matrix(2,1) = sin_theta
            rotation_matrix(2,2) = cos_theta
            
        end select

    end function RotationMatrix

    !--------------------------------------------------------------------
    ! present_or_false
    !
    ! Utility function to safely handle optional logical arguments.
    !
    ! If the optional argument `opt_flag` is present, its value is returned.
    ! If it is not present, the function returns `.false.` by default.
    !--------------------------------------------------------------------
    pure logical function present_or_false(opt_flag)

        logical, intent(in), optional :: opt_flag

        if (present(opt_flag)) then
            present_or_false = opt_flag
        else
            present_or_false = .false.
        end if

    end function present_or_false

    !--------------------------------------------------------------------
    ! amplitude_squared
    !
    ! Purpose:
    !   Computes the **squared modulus** (magnitude squared) of a complex number.
    !
    ! Description:
    !   For a complex number z = x + i*y, the squared modulus is defined as:
    !
    !       |z|^2 = x^2 + y^2 = z * conjg(z)
    !
    !   This is commonly used in Fourier/Ewald calculations to compute
    !   |A(k)|^2 for structure factor amplitudes.
    !--------------------------------------------------------------------
    pure function amplitude_squared(z) result(val)

        ! Input argument
        complex(real64), intent(in) :: z
        ! Output rgument
        real(real64) :: val

        val = real(z*conjg(z), kind=real64)

    end function amplitude_squared

    !----------------------------------------------------------------------------
    ! Function: CrossProduct
    !
    ! Computes the cross product of two 3D vectors a and b.
    ! Given vectors a = (a1, a2, a3) and b = (b1, b2, b3), the cross product
    ! c = a × b is defined as:
    !
    !   c1 = a2*b3 - a3*b2
    !   c2 = a3*b1 - a1*b3
    !   c3 = a1*b2 - a2*b1
    !
    ! The resulting vector c is perpendicular to both a and b, and its magnitude
    ! equals |a||b|sin(θ), where θ is the angle between a and b.
    !----------------------------------------------------------------------------
    function CrossProduct(a, b) result(c)

        ! Input argument
        real(real64), intent(in) :: a(3), b(3)  ! Input vectors
        real(real64) :: c(3)                    ! Resulting cross product vector

        ! Compute cross product components using standard determinant formula
        c(1) = a(2)*b(3) - a(3)*b(2)            ! x-component
        c(2) = a(3)*b(1) - a(1)*b(3)            ! y-component
        c(3) = a(1)*b(2) - a(2)*b(1)            ! z-component

    end function CrossProduct

    !--------------------------------------------------------------------
    ! Function: vector_norm
    !
    ! Computes the Euclidean (L2) norm of a 3-component real vector.
    !
    ! For a vector v = (v1, v2, v3), the norm is defined as:
    !
    !       |v| = sqrt( v1^2 + v2^2 + v3^2 )
    !
    ! This is used frequently in geometry, box operations, and
    ! minimum-image distance calculations.
    !--------------------------------------------------------------------
    pure real(real64) function vector_norm(v)

        real(real64), intent(in) :: v(3)

        vector_norm = sqrt(sum(v * v))

    end function vector_norm

    !--------------------------------------------------------------------
    ! Pads the given column title to 12 characters and adds one space (A12,1X),
    ! ensuring consistent alignment with numeric columns printed using I12.
    ! Appends the formatted field to the growing header line.
    !--------------------------------------------------------------------
    subroutine add_col(line, title)

        ! Input argument
        character(len=*), intent(inout) :: line
        character(len=*), intent(in)    :: title

        ! Local variable
        character(len=16) :: tmp

        ! EXACT match to "(1X,I12)" used by data columns
        write(tmp,'(1X, A12)') title

        ! DO NOT trim, append raw padded field
        line = line // tmp

    end subroutine add_col

    subroutine add_first_col(line, title)

        ! Input arguments
        character(len=*), intent(inout) :: line
        character(len=*), intent(in)    :: title

        ! Local variable
        character(len=16) :: tmp

        ! FIRST column: matches (I12)
        write(tmp,'(A12)') title
        line = line // tmp

    end subroutine add_first_col


end module helper_utils
