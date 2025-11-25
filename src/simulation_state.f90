module simulation_state

    use constants
    use parameters
    use, intrinsic :: iso_fortran_env, only: real64

    implicit none

    !---------------------------------------------------------------------------
    ! Path and file names
    !---------------------------------------------------------------------------
    type path_type
        character(len=200) :: outputs               ! Folder for saving outputs
        character(len=200) :: input                 ! Main input file
        character(len=200) :: topology              ! Topology/data file
        character(len=200) :: parameters            ! Parameters include file
        character(len=200) :: reservoir             ! Optional reservoir file
    end type path_type
    type(path_type) :: path

    !---------------------------------------------------------------------------
    ! Generic status information
    !---------------------------------------------------------------------------
    type status_type
        integer :: desired_block                    ! Desired Monte Carlo block number
        integer :: desired_step                     ! Desired Monte Carlo step
        integer :: block                            ! Monte Carlo block number
        integer :: step                             ! Monte Carlo step within the block
        logical :: reservoir_provided               ! Is reservoir provided by the user
    end type status_type
    type(status_type) :: status

    !---------------------------------------------------------------------------
    ! Counters for Monte Carlo move (trial, success)
    !---------------------------------------------------------------------------
    type :: counter_type
        integer :: creations(2) = 0                 ! Counter for creation moves
        integer :: deletions(2) = 0                 ! Counter for deletion moves
        integer :: translations(2) = 0              ! Counter for translational Monte Carlo moves
        integer :: rotations(2) = 0                 ! Counter for rotational Monte Carlo moves
        integer :: swaps(2) = 0                     ! Counter for swap moves
        integer :: widom(2) = 0                     ! Counter for widom moves
    end type counter_type
    type(counter_type) :: counter

    !---------------------------------------------------------------------------
    ! Desired Monte Carlo move probability
    !---------------------------------------------------------------------------
    type :: proba_type
        real(real64) :: insertion_deletion          ! Probability of attempting an insertion or a deletion
        real(real64) :: translation                 ! Probability of attempting a translation move
        real(real64) :: rotation                    ! Probability of attempting a rotation move
        real(real64) :: swap                        ! Probability of attempting a swap move
        real(real64) :: widom                       ! Probability of attempting a widom move
    end type proba_type
    type(proba_type) :: proba

    !---------------------------------------------------------------------------
    ! Energy terms
    !---------------------------------------------------------------------------
    type :: energy_type
        ! Main terms
        real(real64) :: non_coulomb                 ! Neutral-charged interaction energy
        real(real64) :: coulomb                     ! Charged-electrostatic interaction energy
        real(real64) :: recip_coulomb               ! Reciprocal-space Coulomb contribution
        real(real64) :: ewald_self                  ! Ewald self-interaction energy
        real(real64) :: intra_coulomb               ! Intramolecular Coulomb energy (alternative)
        real(real64) :: total                       ! Total energies
        ! Additional terms
        real(real64) :: self_interaction            ! Site-site short-range energy
        real(real64) :: ke_reciprocal               ! Reciprocal-space (k-space) energy
        real(real64) :: total_coulomb               ! Total Coulomb energy
        real(real64) :: total_non_coulomb           ! Total non-Coulomb energy
    end type energy_type
    type(energy_type) :: energy, old, new           ! Note: old and new are used during Monte Carlo move



    ! Parameters provided in the input file
    type :: input_type
        real(real64), dimension(:), allocatable :: fugacity ! Fugacity of the GCMC reservoir, unitless (for each species)
        real(real64), dimension(:), allocatable :: chemical_potential ! Chemical potential of the GCMC reservoir, kcal/mol (for each species)
        integer, dimension(:), allocatable :: is_active ! Activity flags or counts for each molecule type
        real(real64) :: temperature                 ! Temperature in Kelvin
        real(real64) :: translation_step            ! Maximum displacement for MC moves
        real(real64) :: rotation_step_angle         ! Maximum rotation for MC moves
        real(real64) :: real_space_cutoff           ! Cutoff radius - maximum interaction distance in real space
        integer :: seed                             ! Initial seed for the random number generator
        logical :: recalibrate_moves                ! Enable automatic recalibration of move steps (true/false)
    end type input_type
    type(input_type) :: input

    ! Simulation box definition
    type :: type_box
        ! Generic box geometry parameters
        real(real64) :: reciprocal(3,3)             ! Reciprocal box matrix
        real(real64) :: matrix(3,3)                 ! Simulation box matrix
        real(real64) :: bounds(3,2)                 ! Box dimensions (lo, hi)
        real(real64) :: metrics(9)                  ! Misc. box properties (e.g., lengths, cosines of angles)
        real(real64) :: tilt(3)                     ! Tilt factors (xy, xz, yz)
        real(real64) :: determinant                 ! Volume scaling factor of a linear transformation
        real(real64) :: volume                      ! Box volume
        integer :: type = 0                         ! 0 = unset, 1 = cubic, 2 = orthorhombic, 3 = triclinic
        ! About box content
        integer :: num_atoms                        ! Number of atoms in the box
        integer :: num_atomtypes                    ! Number of atom types
        integer :: num_bonds                        ! Number of bonds
        integer :: num_bondtypes                    ! Number of bond types
        integer :: num_angles                       ! Number of angles
        integer :: num_angletypes                   ! Number of angle types
        integer :: num_dihedrals                    ! Number of dihedrals
        integer :: num_dihedraltypes                ! Number of dihedral types
        integer :: num_impropers                    ! Number of impropers
        integer :: num_impropertypes                ! Number of improper types
        integer, allocatable :: num_residues(:)     ! Number of residue of each type
        ! Atom information
        real(real64), dimension(:), allocatable :: site_masses_vector ! Mass vector for all atom types in initial inputs
        real(real64), dimension(:, :), allocatable :: atom_charges ! Partial charges on sites (1,:,:)=system, (2,:,:)=reservoir
        real(real64), dimension(:, :), allocatable :: atom_masses ! Masses of atoms (1,:,:)=system, (2,:,:)=reservoir
        character(len=10), dimension(:, :), allocatable :: atom_names ! Atom names for each residue (1,:,:)=system, (2,:,:)=reservoir
        integer, dimension(:, :), allocatable :: atom_types ! Atom types for each residue (1,:,:)=system, (2,:,:)=reservoir
        integer, dimension(:, :), allocatable :: atom_ids ! Atom ids for each residue (1,:,:)=system, (2,:,:)=reservoir
    end type type_box
    type(type_box) :: primary, reservoir

    !---------------------------------------------------------------------------
    ! Separate coordinate for host, guest, and gas residue
    !---------------------------------------------------------------------------
    type :: type_coordinate
        integer :: max_nb_atom
        integer :: max_nb_molecule
        logical, dimension(:), allocatable :: residue_exists            ! #tofix : remove?
        real(real64), dimension(:, :, :), allocatable :: com            ! X Y Z coordinate of molecule centers or atoms
        real(real64), dimension(:, :, :, :), allocatable :: offset      ! Local site X Y Z displacements from molecule center
    end type type_coordinate
    ! type(type_coordinate) :: host, guest, gas                           ! Host and gest from the main, gas from reservoir
    type(type_coordinate), target :: host, guest, gas

    integer, allocatable :: resid_location(:) ! size = nb%type_residue

    ! Simulation box definition
    type :: type_number
        integer :: type_residue                     ! Total number of residues
        integer :: max_atom_in_any_residue          ! Max number of atoms in the largest residue
        integer :: max_type_per_residue             ! Max number of type per residue
        integer, dimension(:), allocatable :: atom_in_residue ! Number of atoms in the residue
        integer, dimension(:), allocatable :: types_per_residue ! Number of atom types in the residue
        integer, dimension(:), allocatable :: bonds_per_residue ! Number of bonds in the residue
        integer, dimension(:), allocatable :: angles_per_residue ! Number of angles in the residue
        integer, dimension(:), allocatable :: dihedrals_per_residue ! Number of dihedrals in the residue
        integer, dimension(:), allocatable :: impropers_per_residue ! Number of impropers in the residue
        integer, dimension(:, :), allocatable :: types_pattern ! Type pattern in residue (eg, for TIP4P water 1 2 3 3)
        integer :: max_atom_in_residue_active
        integer :: max_atom_in_residue_inactive
        integer :: max_active_residue
        integer :: max_inactive_residue
        ! integer, allocatable :: residue_count(:)
    end type type_number
    type(type_number) :: nb

    !---------------------------------------------------------------------------
    ! Information extracted during the prescan for each residue type.
    ! Stores the list of atom types, the number of types, and the expected atom
    ! count for that residue as defined in the Maniac input file.
    !---------------------------------------------------------------------------
    type residue
        integer, allocatable :: types(:)    ! List of atom types belonging to this residue
        integer :: nb_types                 ! Number of atom types in this residue
        integer :: nb_atoms                 ! Total number of atoms in this residue
        integer :: nb_res(2)                ! Total number of residue in primary and reservoir 
        logical :: is_active                ! The residue is active
    end type residue
    type(residue), allocatable :: res_infos(:)

    ! Residues information
    type :: type_residue
        real(real64), dimension(:), allocatable :: mass ! de Broglie length
        real(real64), dimension(:), allocatable :: lambda ! Array of mass of residue
        real(real64), dimension(:), allocatable :: masses_1d ! Array of atoms masses
        character(len=10), dimension(:), allocatable :: names_1d ! Array of residue names
        character(len=10), dimension(:, :), allocatable :: names_2d ! Site names for each residue
        integer, dimension(:, :), allocatable :: types_2d ! Site types for each residue
        integer, dimension(:, :, :), allocatable :: bond_type_2d ! Site bonds for each residue
        integer, dimension(:, :, :), allocatable :: angle_type_2d ! Site angles for each residue
        integer, dimension(:, :, :), allocatable :: dihedral_type_2d ! Site dihedrals for each residue
        integer, dimension(:, :, :), allocatable :: improper_type_2d ! Site impropers for each residue
        real(real64), dimension(:, :), allocatable :: site_offset_old ! Local site X Y Z displacements from molecule center
        real(real64), dimension(3) :: mol_com_old   ! For storing old molecule center-of-mass
    end type type_residue
    type(type_residue) :: res

    ! Widom statistic
    type :: type_widom
        real(real64), dimension(:), allocatable :: weight ! Accumulated Boltzmann weight sum for chemical potential
        real(real64), dimension(:), allocatable :: mu_ex ! Excess chemical potential
        real(real64), dimension(:), allocatable :: mu_tot ! Total chemical potential
        integer, dimension(:), allocatable :: sample ! Indices or count of Widom trial samples
    end type type_widom
    type(type_widom) :: widom_stat

    ! Interaction arrays
    type :: type_coeff
        real(real64), dimension(:, :, :, :), allocatable :: sigma ! Lennard-Jones coefficients
        real(real64), dimension(:, :, :, :), allocatable :: epsilon ! Lennard-Jones oefficients
    end type type_coeff
    type(type_coeff) :: coeff

    ! Type to store precomputed reciprocal vectors
    type :: kvector_type
        integer :: kx
        integer :: ky
        integer :: kz
        real(real64) :: k_squared                   ! Normalized k^2 (for validity)
        real(real64) :: k_squared_mag               ! Cartesian squared magnitude (for W(k))
    end type kvector_type

    ! For Ewald calculation
    type :: type_ewald
        integer :: num_kvectors                     ! Number of precomputed vectors
        integer :: kmax(3)                          ! Maximum index for reciprocal lattice vector x y z component
        real(real64) :: alpha                       ! Ewald summation alpha parameter (screening parameter)
        real(real64) :: fourier_precision           ! Estimated precision required for reciprocal (Fourier) space summation
        real(real64) :: screening_factor            ! Intermediate tolerance factor for real-space screening width calculation
        real(real64), dimension(:), allocatable :: recip_constants ! Constants for reciprocal space summations
        complex(real64), dimension(:, :, :, :, :), allocatable :: phase_factor   ! Complex exponentials for reciprocal space
        complex(real64), dimension(:, :, :), allocatable :: phase_factor_old ! Old complex exponential terms
        complex(real64), dimension(:), allocatable :: recip_amplitude ! Fourier coefficients of charge density or potential
        complex(real64), dimension(:), allocatable :: recip_amplitude_old ! Old fourier coefficients of charge density or potential
        real(real64), dimension(:), allocatable :: form_factor ! Factor to account for symmetry (k vs -k)
        type(kvector_type), allocatable :: kvectors(:) ! Precomputed reciprocal vectors
        complex(real64), dimension(:, :), allocatable :: temp ! Temporary Fourier array
        complex(real64), dimension(:), allocatable :: temp_1d ! Temporary Fourier array
        complex(real64), dimension(:), allocatable :: phase_new  ! Temporary array for new configuration phases
        complex(real64), dimension(:), allocatable :: phase_old  ! Temporary array for old configuration phases
        real(real64), dimension(:), allocatable :: charges ! Temporary array for atom charges    
        real(real64) :: tolerance                   ! Numerical accuracy for Ewald summation,
    end type type_ewald
    type(type_ewald) :: ewald

    ! For tabulating potential
    type :: tabulated
        ! Table properties
        real(real64), allocatable :: x(:)           ! Grid points (r, k^2, etc.)
        real(real64), allocatable :: f(:)           ! Function values
        real(real64) :: dx                          ! Grid spacing
        integer :: n                                ! Number of points
        logical :: initialized = .false.            ! Flag to indicate table is ready
    end type tabulated
    type(tabulated) :: erfc_r_table                 ! For precomputed erfc(r) / r
    type(tabulated) :: r6_table                     ! For precomputed r**6
    type(tabulated) :: r12_table                    ! For precomputed r**12

end module simulation_state
