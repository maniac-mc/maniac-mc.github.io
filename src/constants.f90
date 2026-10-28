!===========================================================================
! Physical and mathematical constants used throughout MANIAC
!===========================================================================
module constants

    use, intrinsic :: iso_fortran_env, only: real64
    use parameters, only: J_to_kcal, m_to_A

    implicit none

    ! --- Mathematical constants ---
    real(real64), parameter :: PI = 3.141592653589793238462643383279502884_real64 ! Constant Pi
    real(real64), parameter :: TWOPI = 2.0_real64 * PI      ! Two times Pi
    real(real64), parameter :: SQRTPI = SQRT(PI)            ! Square root of Pi

    ! --- Physical constants (SI) ---
    real(real64), parameter :: H_PLANCK = 6.62607015e-34_real64 ! Planck constant (J s)
    real(real64), parameter :: EPS0 = 8.8541878128e-12_real64 ! F/m = C^2/(N m^2)
    real(real64), parameter :: NA = 6.02214076e23_real64    ! Na (mol-1)
    real(real64), parameter :: KB = 1.380649e-23_real64     ! Boltzmann constant (J/K)
    real(real64), parameter :: E_CHARGE  = 1.602176634e-19_real64 ! Elementary charge (C)

    ! --- Derived constants ---
    ! Coulomb prefactor e^2 / (4 pi eps0) in SI units (J m or N m^2)
    real(real64), parameter :: EPS0_INV = E_CHARGE**2 / (4.0_real64 * PI * EPS0)
    ! Coulomb prefactor in simulation units (kcal Angstrom / (mol e^2))
    real(real64), parameter :: EPS0_INV_real = EPS0_INV * J_to_kcal * m_to_A * NA
    ! Boltzmann constant in kcal/(mol K)
    real(real64), parameter :: KB_kcalmol = KB * NA * J_to_kcal

    ! --- Runtime variable (NOT a constant, #tocheck : move somewhere else ?) ---
    real(real64) :: beta ! 1 / (kB T)

    ! --- Numeric literals in real64 (to avoid writing 2.0_real64 everywhere) ---
    real(real64), parameter :: zero = 0.0_real64
    real(real64), parameter :: quarter = 0.25_real64
    real(real64), parameter :: half = 0.5_real64
    real(real64), parameter :: one = 1.0_real64
    real(real64), parameter :: two = 2.0_real64
    real(real64), parameter :: three = 3.0_real64
    real(real64), parameter :: four = 4.0_real64
    real(real64), parameter :: ten = 10.0_real64

end module constants
