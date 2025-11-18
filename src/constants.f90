module constants

    use, intrinsic :: iso_fortran_env, only: real64
    use parameters

    implicit none

    real(real64), parameter :: PI = 3.1415926536_real64     ! Constant Pi
    real(real64), parameter :: TWOPI = 2.0_real64 * PI      ! Two times Pi
    real(real64), parameter :: SQRTPI = SQRT(PI)            ! Square root of Pi
    real(real64), parameter :: H_PLANCK = 6.62607015d-34    ! Planck constant (J s)
    real(real64), parameter :: EPS0 = 8.854187817d-12       ! F/m = C^2/(N m^2)
    real(real64), parameter :: NA = 6.02214076e23           ! Na (mol-1)
    real(real64), parameter :: KB = 1.380658d-23            ! Boltzmann constant (J/K)
    real(real64), parameter :: E_CHARGE  = 1.602176634d-19  ! Elementary charge (C)

    real(real64), parameter :: EPS0_INV = E_CHARGE**2 / (4.0_real64 * PI * EPS0) ! Coulomb prefactor e^2 / 4 pi epsilon_0, SI units, N m^2
    real(real64), parameter :: EPS0_INV_real = EPS0_INV * J_to_kcal * m_to_A * NA ! Coulomb prefactor in real units, 

    ! TO REMOVE -----------------
    real(real64), parameter :: KB_JK = 1.380658D-23 ! Boltzmann constant (J/K)
    real(real64), parameter :: KB_kcalmol = 0.0019872041_real64 ! Boltzmann constant (kcal/mol)
    ! real(real64), parameter :: EPS0_INV_kcalA = 332.0637_real64 ! e^2 / 4 pi epsilon_0 (kcal Å/(mol e^2))
    real(real64), parameter :: HBAR = 1.05457182D-34    ! hbar (m^2 kg / s)
    ! TO REMOVE -----------------

    real(real64), parameter :: zero = 0.0_real64        ! Zero in real64
    real(real64), parameter :: quarter = 0.25_real64    ! 1/4 in real64
    real(real64), parameter :: half = 0.5_real64        ! 1/2 in real64
    real(real64), parameter :: one = 1.0_real64         ! One in real64
    real(real64), parameter :: two = 2.0_real64         ! Two in real64
    real(real64), parameter :: three = 3.0_real64       ! Three in real64
    real(real64), parameter :: four = 4.0_real64        ! Four in real64
    real(real64), parameter :: ten = 10.0_real64        ! Ten in real64

end module constants
