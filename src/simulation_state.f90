module simulation_state

    use constants
    use parameters
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    !---------------------------------------------------------------------------
    ! Path and file names
    !---------------------------------------------------------------------------
    type path_type
        character(len=LENPATH) :: input                 ! Main input file
        character(len=LENPATH) :: outputs               ! Folder for saving outputs
        character(len=LENPATH) :: topology              ! Topology/data file
        character(len=LENPATH) :: parameters            ! Parameters include file
        character(len=LENPATH) :: reservoir             ! Optional reservoir file
    end type path_type
    type(path_type) :: path

    !---------------------------------------------------------------------------
    ! Generic simulation status information
    !---------------------------------------------------------------------------
    type status_type
        integer :: step                             ! Monte Carlo step within the block
        integer :: block                            ! Monte Carlo block number
        integer :: desired_step                     ! Desired Monte Carlo step
        integer :: desired_block                    ! Desired Monte Carlo block number
        logical :: reservoir_provided               ! Is reservoir provided by the user
    end type status_type
    type(status_type) :: status

    !---------------------------------------------------------------------------
    ! Counters for Monte Carlo move (trial, success)
    !---------------------------------------------------------------------------
    type counter_type
        integer :: translations(2) = 0              ! Translational Monte Carlo moves
        integer :: rotations(2) = 0                 ! Rotational Monte Carlo moves
        integer :: creations(2) = 0                 ! Creation moves
        integer :: deletions(2) = 0                 ! Deletion moves
        integer :: swaps(2) = 0                     ! Swap moves
        integer :: widom(2) = 0                     ! Widom moves
    end type counter_type
    type(counter_type) :: counter

    !---------------------------------------------------------------------------
    ! Desired Monte Carlo move probability
    !---------------------------------------------------------------------------
    type proba_type
        real(real64) :: insertion_deletion          ! Probability of attempting an insertion or a deletion
        real(real64) :: translation                 ! Probability of attempting a translation move
        real(real64) :: rotation                    ! Probability of attempting a rotation move
        real(real64) :: widom                       ! Probability of attempting a widom move
        real(real64) :: swap                        ! Probability of attempting a swap move
    end type proba_type
    type(proba_type) :: proba

    !---------------------------------------------------------------------------
    ! Energy terms
    !---------------------------------------------------------------------------
    type energy_type
        real(real64) :: non_coulomb                 ! Neutral-charged interaction energy
        real(real64) :: coulomb                     ! Charged-electrostatic interaction energy
        real(real64) :: recip_coulomb               ! Reciprocal-space Coulomb contribution
        real(real64) :: ewald_self                  ! Ewald self-interaction energy
        real(real64) :: intra_coulomb               ! Intramolecular Coulomb energy
        real(real64) :: total                       ! Total energies
    end type energy_type
    type(energy_type) :: energy, old, new           ! Note: old and new are used during Monte Carlo move

    !---------------------------------------------------------------------------
    ! Thermodynamic input parameters
    !---------------------------------------------------------------------------
    type input_thermo_type
        real(real64) :: temperature                     ! System temperature (K)
        real(real64), allocatable :: fugacity(:)        ! Species fugacity for GCMC reservoir (dimensionless)
        real(real64), allocatable :: chemical_potential(:) ! Species chemical potential of reservoir (kcal/mol)
        logical, allocatable :: is_active(:)            ! Activity flag: 1 = species active, 0 = inactive
    end type input_thermo_type
    type(input_thermo_type) :: thermo

    !---------------------------------------------------------------------------
    ! Coordinate for host, guest, and gas residue
    !---------------------------------------------------------------------------
    type type_coordinate
        real(real64), dimension(:, :, :), allocatable :: com            ! X Y Z coordinate of molecule centers or atoms
        real(real64), dimension(:, :, :, :), allocatable :: offset      ! Local site X Y Z displacements from molecule center
    end type type_coordinate
    type(type_coordinate), target :: host, guest, gas

    !---------------------------------------------------------------------------
    ! Old coordinate used in Monte Carlo move
    !---------------------------------------------------------------------------
    type saved_coordinate
        real(real64), dimension(:), allocatable :: com                  ! X Y Z coordinate of molecule centers or atoms
        real(real64), dimension(:, :), allocatable :: offset            ! Local site X Y Z displacements from molecule center
    end type saved_coordinate
    type(saved_coordinate), target :: saved

    !---------------------------------------------------------------------------
    ! Simulation cell geometrical properties
    !---------------------------------------------------------------------------
    type type_cell
        real(real64) :: reciprocal(3,3)             ! Reciprocal box matrix
        real(real64) :: matrix(3,3)                 ! Simulation box matrix
        real(real64) :: bounds(3,2)                 ! Box dimensions (lo, hi)
        real(real64) :: metrics(9)                  ! Misc. box properties (e.g., lengths, cosines of angles)
        real(real64) :: tilt(3)                     ! Tilt factors (xy, xz, yz)
        real(real64) :: determinant                 ! Volume scaling factor of a linear transformation
        real(real64) :: volume                      ! Box volume
        integer :: shape                            ! To differentiace between triclinic and orthorhombic boc
    end type type_cell

    !---------------------------------------------------------------------------
    ! Simulation number topology
    !---------------------------------------------------------------------------
    type type_topology
        integer :: atoms                                ! Number of atoms
        integer :: bonds                                ! Number of bonds
        integer :: angles                               ! Number of angles
        integer :: dihedrals                            ! Number of dihedrals
        integer :: impropers                            ! Number of impropers
        integer :: atomtypes                            ! Number of atom types
        integer :: bondtypes                            ! Number of bond types
        integer :: angletypes                           ! Number of angle types
        integer :: dihedraltypes                        ! Number of dihedral types
        integer :: impropertypes                        ! Number of improper types
        integer, allocatable :: residues(:)             ! Number of residues of each type
    end type type_topology

    !---------------------------------------------------------------------------
    ! Atom-specific data
    !---------------------------------------------------------------------------
    type type_atomdata
        real(real64), allocatable :: masses_vec(:)      ! Mass vector for all atom types
        real(real64), allocatable :: charges(:,:)       ! Partial charges on sites (system/reservoir)
        real(real64), allocatable :: masses(:,:)        ! Masses of atoms (system/reservoir)
        character(len=10), allocatable :: names(:,:)    ! Atom names for each residue (system/reservoir)
        integer, allocatable :: types(:,:)              ! Atom types for each residue (system/reservoir)
        integer, allocatable :: ids(:,:)                ! Atom ids for each residue (system/reservoir)
    end type type_atomdata

    !---------------------------------------------------------------------------
    ! Simulation box definition
    !---------------------------------------------------------------------------
    type type_box
        type(type_cell) :: cell                         ! Geometrical properties of the simulation cell
        type(type_topology) :: num                      ! Number of atoms, bonds, angles, etc.
        type(type_atomdata) :: atoms                    ! Atom-specific data (masses, charges, types, ids)
    end type type_box
    type(type_box) :: primary, reservoir                ! Primary simulation box and optional reservoir box

    !---------------------------------------------------------------------------
    ! Maximum counts for residue and atom-type characteristics
    !---------------------------------------------------------------------------
    type type_max
        integer :: atoms_per_residue          ! Max number of atoms in any residue
        integer :: types_per_residue          ! Max number of atom types in any residue
        integer :: active_residues            ! Max number of active residues
        integer :: inactive_residues          ! Max number of inactive residues
        integer :: atoms_active_residue       ! Max number of atoms in active residues
        integer :: atoms_inactive_residue     ! Max number of atoms in inactive residues
    end type type_max
    type (type_max) :: nmax

    !---------------------------------------------------------------------------
    ! Statistic accumulator
    !---------------------------------------------------------------------------
    type mc_stat_type
        integer, dimension(:), allocatable :: sample        ! Indices or count of Widom trial samples
        real(real64), dimension(:), allocatable :: mu_ex    ! Excess chemical potential
        real(real64), dimension(:), allocatable :: mu_tot   ! Total chemical potential
        real(real64), dimension(:), allocatable :: weight   ! Accumulated Boltzmann weight sum for chemical potential
    end type mc_stat_type
    type(mc_stat_type) :: statistic

    !---------------------------------------------------------------------------
    ! Interaction arrays
    !---------------------------------------------------------------------------
    type type_coeff
        real(real64), dimension(:, :, :, :), allocatable :: sigma ! Lennard-Jones coefficients
        real(real64), dimension(:, :, :, :), allocatable :: epsilon ! Lennard-Jones oefficients
    end type type_coeff
    type(type_coeff) :: coeff

   !---------------------------------------------------------------------------
    ! For tabulating potential
    !---------------------------------------------------------------------------
    type tabulated
        ! Table properties
        real(real64), allocatable :: x(:)           ! Grid points (r, k^2, etc.)
        real(real64), allocatable :: f(:)           ! Function values
        real(real64) :: dx                          ! Grid spacing
        integer :: npoint                           ! Number of points
        logical :: initialized = .false.            ! Flag to indicate table is ready
    end type tabulated
    type(tabulated) :: erfc_r_table                 ! For precomputed erfc(r) / r
    type(tabulated) :: r6_table                     ! For precomputed r**6
    type(tabulated) :: r12_table                    ! For precomputed r**12

    !---------------------------------------------------------------------------
    ! Intramolecular connectivity of a residue
    !---------------------------------------------------------------------------
    type type_connections
        integer, dimension(:, :, :), allocatable :: bonds       ! Site bonds for each residue
        integer, dimension(:, :, :), allocatable :: angles      ! Site angles for each residue
        integer, dimension(:, :, :), allocatable :: dihedrals   ! Site dihedrals for each residue
        integer, dimension(:, :, :), allocatable :: impropers   ! Site impropers for each residue
    end type type_connections
    type(type_connections) :: connect                           ! Array of residue definitions and intramolecular connectivity

    !---------------------------------------------------------------------------
    ! Per residues information
    !---------------------------------------------------------------------------
    type type_residue
        integer :: number                                       ! Total number of residues
        integer, allocatable :: role(:)                         ! Guest or host
        integer, dimension(:), allocatable :: atom              ! Number of atoms in the residue
        real(real64), dimension(:), allocatable :: mass         ! Array of mass of residue
        real(real64), dimension(:), allocatable :: lambda       ! De Broglie length
        integer, dimension(:, :), allocatable :: site_types(:, :) ! Site types for each residue
        integer, dimension(:, :), allocatable :: pattern_types(:, :) ! Type pattern in residue (eg, for TIP4P water 1 2 3 3)
        integer, dimension(:), allocatable :: types(:)          ! Number of atom types in the residue
        integer, allocatable :: bonds(:)                        ! Number of bonds in the residue
        integer, allocatable :: angles(:)                       ! Number of angles in the residue
        integer, allocatable :: dihedrals(:)                    ! Number of dihedrals in the residue
        integer, allocatable :: impropers(:)                    ! Number of impropers in the residue
        character(len=10), allocatable :: names(:)              ! Array of residue names
        character(len=10), dimension(:, :), allocatable :: site_names ! Site names for each residue
    end type type_residue
    type(type_residue) :: res

    !---------------------------------------------------------------------------
    ! Monte Carlo parameters provided in the input file
    !---------------------------------------------------------------------------
    type input_type
        real(real64) :: translation_step            ! Maximum displacement for MC moves
        real(real64) :: rotation_step_angle         ! Maximum rotation for MC moves
        real(real64) :: real_space_cutoff           ! Cutoff radius - maximum interaction distance in real space
        logical :: recalibrate_moves                ! Enable automatic recalibration of move steps (true/false)
        integer :: seed                             ! Initial seed for the random number generator
    end type input_type
    type(input_type) :: mc_input

    !---------------------------------------------------------------------------
    ! Type to store precomputed reciprocal vectors (used in ewald)
    !---------------------------------------------------------------------------
    type kvector_type
        integer :: kx
        integer :: ky
        integer :: kz
        real(real64) :: k_squared                   ! Normalized k^2 (for validity)
        real(real64) :: k_squared_mag               ! Cartesian squared magnitude (for W(k))
    end type kvector_type

    !---------------------------------------------------------------------------
    ! Allocatable array for phase calculation
    !---------------------------------------------------------------------------
    type phase_ewald
        complex(real64), allocatable :: factor(:,:,:,:,:)       ! Complex exponentials for reciprocal space
        complex(real64), allocatable :: factor_old(:,:,:)       ! Old complex exponential terms
        complex(real64), allocatable :: new(:)                  ! Temporary array for new configuration phases
        complex(real64), allocatable :: old(:)                  ! Temporary array for old configuration phases
    end type phase_ewald

    !---------------------------------------------------------------------------
    ! Input parameters for ewald calculations
    !---------------------------------------------------------------------------
    type ewald_parameters
        integer :: nkvec                            ! Number of precomputed vectors
        integer :: kmax(3)                          ! Maximum index for reciprocal lattice vector x y z component
        real(real64) :: alpha                       ! Ewald summation alpha parameter (screening parameter)
        real(real64) :: tolerance                   ! Numerical accuracy for Ewald summation
        real(real64) :: fprecision                  ! Estimated precision for reciprocal space summation
        real(real64) :: screen                      ! Real-space screening factor
    end type ewald_parameters

    ! !---------------------------------------------------------------------------
    ! ! Runtime arrays for Ewald calculations
    ! !---------------------------------------------------------------------------
    ! type ewald_runtime
    !     real(real64), allocatable :: recip_constants(:)      ! Constants for reciprocal space summations
    !     complex(real64), allocatable :: recip_amplitude(:)   ! Fourier coefficients of charge density or potential
    !     complex(real64), allocatable :: recip_amplitude_old(:) ! Old Fourier coefficients
    !     real(real64), allocatable :: form_factor(:)          ! Symmetry factor for k vs -k
    !     complex(real64), allocatable :: fft2d(:, :)         ! 2D FFT temporary array
    !     complex(real64), allocatable :: fft1d(:)            ! 1D FFT temporary array
    !     real(real64), allocatable :: charges(:)             ! Temporary array for atom charges
    ! end type ewald_runtime

    !---------------------------------------------------------------------------
    ! For Ewald calculation
    !---------------------------------------------------------------------------
    type type_ewald
        real(real64), allocatable :: recip_constants(:) ! Constants for reciprocal space summations
        complex(real64), allocatable :: recip_amplitude(:) ! Fourier coefficients of charge density or potential
        complex(real64), allocatable :: recip_amplitude_old(:) ! Old fourier coefficients of charge density or potential
        real(real64), allocatable :: form_factor(:) ! Factor to account for symmetry (k vs -k)
        complex(real64), allocatable :: temp(:, :) ! Temporary Fourier array
        complex(real64), allocatable :: temp_1d(:) ! Temporary Fourier array
        real(real64), allocatable :: charges(:) ! Temporary array for atom charges    
        type(ewald_parameters) :: param
        type(phase_ewald) :: phase
        type(kvector_type), allocatable :: kvectors(:) ! Precomputed reciprocal vectors
    end type type_ewald
    type(type_ewald) :: ewald
 
end module simulation_state
