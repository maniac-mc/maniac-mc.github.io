module parameters

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    integer, parameter :: TYPE_HOST = 1, TYPE_GUEST = 2, TYPE_GAS = 3

    ! Integer parameters defining maximum sizes and counts
    integer, parameter :: NB_MAX_MOLECULE = 5000    ! Maximum number of molecules of a given type
    integer, parameter :: NB_MAX_BOND = 20000       ! Maximum number of bond per residue
    integer, parameter :: NB_MAX_ANGLE = 20000      ! Maximum number of angle per residue
    integer, parameter :: NB_MAX_DIHEDRAL = 20000   ! Maximum number of dihedral per residue
    integer, parameter :: NB_MAX_IMPROPER = 20000   ! Maximum number of improper per residue

    ! Parameters controlling the MC adjustement
    real(real64), parameter :: TARGET_ACCEPTANCE = 0.40d0    ! Target acceptance ratio for adaptive MC moves
    real(real64), parameter :: TOL_ACCEPTANCE    = 0.05d0    ! Tolerance window around target acceptance
    real(real64), parameter :: MIN_TRANSLATION_STEP = 1.0d-3 ! Minimum allowed translation displacement (Å)
    real(real64), parameter :: MAX_TRANSLATION_STEP = 3.0d0  ! Maximum allowed translation displacement (Å)
    real(real64), parameter :: MIN_ROTATION_ANGLE  = 1.0d-3  ! Minimum allowed rotation angle (rad)
    real(real64), parameter :: MAX_ROTATION_ANGLE  = 0.78d0  ! Maximum allowed rotation angle (~45° in rad)
    real(real64), parameter :: PROB_CREATE_DELETE = 0.5d0    ! Probability for insertion/deletion
    integer, parameter :: MIN_TRIALS_FOR_RECALIBRATION = 500 ! Minimum trials before recalibrating step sizes

    ! Parameters controlling character length
    integer, parameter :: BOX_WIDTH = 78                    ! Box in outputed log
    integer, parameter :: LENPATH = 200                     ! Path and names length

    ! Parameters for unit conversion
    real(real64), parameter :: A_TO_M = 1.0d-10         ! Conversion factor: 1 Å → 1.0E-10 m
    real(real64), parameter :: M_TO_A = 1.0d10          ! Conversion factor: 1 m → 1.0E10 Å
    real(real64), parameter :: A3_TO_M3 = 1.0d-30       ! Conversion factor: 1 Å^3 → 1.0E-30 m^3
    real(real64), parameter :: ATM_TO_PA = 1.01325d5    ! Conversion factor: 1 atm → 1.01325 × 10^5 Pa
    real(real64), parameter :: G_TO_KG = 1.0d-3         ! Conversion factor: 1 g → 1.0E-3 km

    real(real64), parameter :: J_to_kcal  = 0.000239005736  ! 1 J = 0.239005736 cal
    real(real64), parameter :: cal_to_kcal = 1.0e-3

    ! File names
    character(len=*), parameter :: data_filename = 'topology.data' ! Topology data file

    ! Define constants for move types
    integer, parameter :: TYPE_CREATION = 1             ! Move type: molecule creation
    integer, parameter :: TYPE_DELETION = 2             ! Move type: molecule deletion
    integer, parameter :: TYPE_TRANSLATION = 3          ! Move type: molecular translation
    integer, parameter :: TYPE_ROTATION = 4             ! Move type: molecular rotation
    
    ! Define constants for box types
    integer, parameter :: ORTHORHOMBIC = 1
    integer, parameter :: TRICLINIC = 2

    ! In module parameters
    integer, parameter :: TABULATED_POINTS = 5000 ! Default number of points in tabulated functions
    logical, parameter :: use_table = .false. ! #TODO: make it an input parameter for user

    real(real64), parameter :: error = 1.0D-10              ! Small number for error calculation

end module parameters
