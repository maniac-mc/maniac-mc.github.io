module constants

    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    real(real64), parameter :: PI = 3.14159265358979323846_real64 ! Constant Pi
    real(real64), parameter :: TWOPI = 2.0_real64 * PI ! Two times Pi
    real(real64), parameter :: SQRTPI = SQRT(PI) ! Square root of Pi
    real(real64), parameter :: KB_JK = 1.380658D-23 ! Boltzmann constant (J/K)
    real(real64), parameter :: KB_kcalmol = 0.0019872041_real64 ! Boltzmann constant (kcal/mol)
    real(real64), parameter :: EPS0_INV_eVA = 14.40198_real64 ! e^2 / 4 pi epsilon_0 (eV Å / e^2)
    real(real64), parameter :: EPS0_INV_kcalA = 332.0637_real64 ! e^2 / 4 pi epsilon_0 (kcal Å/(mol e^2))
    real(real64), parameter :: KB_eVK = 8.6173852D-5    ! Boltzmann constant (eV/K)
    real(real64), parameter :: zero = 0.0_real64        ! Zero in real64
    real(real64), parameter :: quarter = 0.25_real64    ! 1/4 in real64
    real(real64), parameter :: half = 0.5_real64        ! 1/2 in real64
    real(real64), parameter :: one = 1.0_real64         ! One in real64
    real(real64), parameter :: two = 2.0_real64         ! Two in real64
    real(real64), parameter :: three = 3.0_real64       ! Three in real64
    real(real64), parameter :: four = 4.0_real64        ! Four in real64
    real(real64), parameter :: ten = 10.0_real64        ! Ten in real64
    real(real64), parameter :: error = 1.0D-10          ! Small number for error calculation
    real(real64), parameter :: Hplank = 6.62607015D-34  ! Planck constant (J s)
    real(real64), parameter :: HBAR = 1.05457182D-34    ! hbar (m^2 kg / s)
    real(real64), parameter :: NA = 6.022e23            ! Na (mol-1)
end module constants
