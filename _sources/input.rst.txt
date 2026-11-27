Input Files
===========

This page describes the structure and parameters of the **input file**
passed to **MANIAC** with the ``-i`` command-line flag (see :doc:`build` for
the full command line). Each line specifies a parameter and its value.
Lines starting with ``#`` are comments, blank lines are ignored, and any
unknown keyword outside a residue block makes the program abort with an
error.

All quantities use the following units: lengths in Å, energies in kcal/mol,
temperatures in K, angles in radians.

Minimal Example
---------------

Below is a minimal example for water adsorption in a rigid ZIF-8 framework
(the active species is water, the host framework is rigid).

.. code-block:: text

    # Input for maniac program

    # Generic parameters
    nb_block 200                     # Number of Monte Carlo blocks
    nb_step 1000                     # Number of Monte Carlo steps per block
    temperature 300                  # Simulation temperature (K)
    seed 8743                        # RNG seed (omit for a random seed)

    # Electrostatic interactions
    ewald_tolerance 1e-5             # Ewald convergence tolerance
    real_space_cutoff 17             # Real-space cutoff for electrostatics (Å)

    # Monte Carlo moves
    translation_step 1               # Maximum translational displacement (Å)
    rotation_step_angle 0.685        # Maximum rotation angle (radians)
    recalibrate_moves true           # Automatically tune step sizes for target acceptance

    translation_proba 0.4            # Relative weight of translation moves
    rotation_proba 0.4               # Relative weight of rotation moves
    insertion_deletion_proba 0.2     # Relative weight of insertion/deletion moves

    # Define residue composition and activity
    begin_residue
      name wat                       # Residue name: water
      state actif                    # "actif" = subject to MC moves
      fugacity 0.5                   # Reservoir fugacity (activity) of water
      types 1 2 3                    # Atom type IDs present in the residue
      names Ow Mw Hw                 # Atom names: O, virtual site, H
      nb-atoms 4                     # Atoms per molecule (TIP4P: Ow Mw Hw Hw)
    end_residue

    begin_residue
      name zif                       # Residue name: ZIF-8 framework
      state inactif                  # "inactif" = rigid, no MC moves applied
      types 4 5 6 7 8 9 10 11        # Atom type IDs of the framework
      names Zn C1 C2 H1 C3 H2 H3 N
      nb-atoms 2208                  # Total number of atoms in the framework
    end_residue

This example corresponds to water adsorption in a rigid ZIF-8:

- The first residue block defines a TIP4P water molecule (``wat``).
- It is marked as **actif**, meaning Monte Carlo moves (translation,
  rotation, insertion, deletion) will be attempted on it.
- ``fugacity`` is the reservoir fugacity, used as a dimensionless activity:
  internally ``mu = ln(fugacity) / beta`` (kcal/mol). No unit conversion is
  applied, so do not read the value as a pressure in atm or bar.
- ``types 1 2 3`` and ``names Ow Mw Hw`` list the *distinct* site types of the
  molecule (oxygen Ow, virtual site Mw, hydrogen Hw), while ``nb-atoms 4``
  counts all atoms including the repeated Hw site.
- The second residue block defines the **ZIF-8 framework** (``zif``). It is
  marked as **inactif**, so it remains rigid and is not subject to Monte
  Carlo moves. Its types/names describe the crystallographic composition,
  and ``nb-atoms 2208`` is the total number of framework atoms.

Simulation Control
------------------

**nb_block** (integer, *required*)
  Number of Monte Carlo **blocks** to run. A block is a group of MC steps
  used for statistical averaging (one line per block in ``energy.dat``,
  ``moves.dat``, ``number_<name>.dat``).

**nb_step** (integer, *required*)
  Number of Monte Carlo **steps** per block. One move is attempted per step.

**temperature** (float, *required*, K, must be > 0)
  Target simulation temperature.

**seed** (integer, *optional*, default: random)
  Initial seed of the random number generator. The same seed reproduces the
  same trajectory; if omitted, a random seed is drawn.

Electrostatics
--------------

**ewald_tolerance** (float, *required*, must be > 0)
  Convergence tolerance of the Ewald summation for electrostatic
  interactions. Smaller values are more accurate but more expensive.

**real_space_cutoff** (float, *required*, Å, must be > 0)
  Real-space cutoff distance for electrostatics (and for Lennard-Jones
  interactions).

Monte Carlo Moves
-----------------

**translation_step** (float, *required*, Å)
  Maximum translational displacement of a trial move.

**rotation_step_angle** (float, *required*, radians)
  Maximum rotation angle of a trial move.

**recalibrate_moves** (boolean, *optional*, default: ``false``)
  If ``true``, step sizes are automatically adjusted during the run to keep
  the acceptance ratio near its target.

**histogram_every** (integer, *optional*, default: ``nb_step``, must be >= 1)
  Sampling interval, in MC steps, of the particle-number histogram
  (``histogram_<name>.dat``). The default samples once per block.

Move probabilities
##################

The five move types are selected with relative weights that do **not** need
to sum to 1: when all probabilities are provided, **MANIAC** automatically
rescales them (with a warning) so that their total equals 1. This ensures
that at each Monte Carlo step, one move is chosen according to the relative
weights of the provided probabilities.

- **translation_proba** (float, *optional*, default: 0.0) — translation moves.
- **rotation_proba** (float, *optional*, default: 0.0) — rotation moves.
- **insertion_deletion_proba** (float, *optional*, default: 0.0) — insertion
  or deletion of a whole molecule (grand-canonical move).
- **swap_proba** (float, *optional*, default: 0.0) — swap of a molecule with
  another species of the reservoir (semi-grand move).
- **widom_proba** (float, *optional*, default: 0.0) — Widom test-particle
  insertions (used to measure the chemical potential).

Rules enforced at parse time:

- At least one probability must be non-zero, otherwise the run aborts.
- ``widom_proba`` and ``insertion_deletion_proba`` cannot both be non-zero.

Residue Definition
------------------

A residue block defines one molecule type (or one rigid framework). Each
block starts with ``begin_residue`` and ends with ``end_residue``. After
parsing, residues are sorted internally by their smallest atom type ID.

Residue parameters:

- **name** (string) — residue name, used in output file names
  (``number_<name>.dat``, ``histogram_<name>.dat``).
- **state** (string) — ``actif`` (subject to MC moves) or ``inactif``
  (rigid; only wrapped at the periodic boundaries). Any other value aborts.
- **types** (list of integers) — atom type IDs present in the residue. They
  must match the type IDs used for this residue in the LAMMPS data file
  (``-d``). Note this lists the *distinct* site types: for TIP4P water,
  ``types 1 2 3`` covers Ow, Mw and Hw even though the molecule has 4 atoms.
- **names** (list of strings) — atom names corresponding to ``types``
  (e.g. ``Ow Mw Hw``).
- **nb-atoms** (integer) — total number of atoms in one molecule, including
  repeated sites (so ``nb-atoms 4`` for TIP4P water with 3 listed types).
- **fugacity** (float, must be > 0) — reservoir fugacity (activity) of the
  residue. Internally the chemical potential is derived as
  ``mu = ln(fugacity) / beta`` (kcal/mol).
- **chemical_potential** (float, kcal/mol) — alternative to ``fugacity``:
  specify the chemical potential directly.

Rules for ``fugacity`` / ``chemical_potential``:

- They are mutually exclusive for a given residue (specifying both aborts).
- For an ``actif`` residue, one of them is *required* whenever
  ``insertion_deletion_proba > 0`` or ``swap_proba > 0`` (the moves need a
  reservoir chemical potential). Pure NVT runs (translation/rotation only)
  do not need either.

Residue Definition
------------------

A residue block defines one molecule type (or one rigid framework). Each
block starts with ``begin_residue`` and ends with ``end_residue``. After
parsing, residues are sorted internally by their smallest atom type ID.

Residue parameters:

- **name** (string) — residue name, used in output file names
  (``number_<name>.dat``, ``histogram_<name>.dat``).
- **state** (string) — ``actif`` (subject to MC moves) or ``inactif``
  (rigid; only wrapped at the periodic boundaries). Any other value aborts.
- **types** (list of integers) — atom type IDs present in the residue. They
  must match the type IDs used for this residue in the LAMMPS data file
  (``-d``). Note this lists the *distinct* site types: for TIP4P water,
  ``types 1 2 3`` covers Ow, Mw and Hw even though the molecule has 4 atoms.
- **names** (list of strings) — atom names corresponding to ``types``
  (e.g. ``Ow Mw Hw``).
- **nb-atoms** (integer) — total number of atoms in one molecule, including
  repeated sites (so ``nb-atoms 4`` for TIP4P water with 3 listed types).
- **fugacity** (float, must be > 0) — reservoir fugacity (activity) of the
  residue. Internally the chemical potential is derived as
  ``mu = ln(fugacity) / beta`` (kcal/mol).
- **chemical_potential** (float, kcal/mol) — alternative to ``fugacity``:
  specify the chemical potential directly.

Rules for ``fugacity`` / ``chemical_potential``:

- They are mutually exclusive for a given residue (specifying both aborts).
- For an ``actif`` residue, one of them is *required* whenever
  ``insertion_deletion_proba > 0`` or ``swap_proba > 0`` (the moves need a
  reservoir chemical potential). Pure NVT runs (translation/rotation only)
  do not need either.

Files and outputs
-----------------

The input file itself only holds the parameters above; the simulation files
are given on the command line (see :doc:`build`):

- ``-i`` : this input file (``.maniac``)
- ``-d`` : LAMMPS data file of the simulation box (``topology.data``)
- ``-p`` : LAMMPS parameter file (``parameters.inc``: masses, pair/bond/angle
  coefficients)
- ``-o`` : output folder (created if missing); receives ``log.maniac``,
  ``outputs/topology.data``, ``trajectory.lammpstrj``, ``energy.dat``,
  ``moves.dat``, ``number_<name>.dat``, ``histogram_<name>.dat``, ...
- ``-r`` : optional LAMMPS data file of the molecule reservoir (required for
  swap moves)

Notes
-----
- Parameter order is not strict, but each residue must be fully enclosed
  between ``begin_residue`` and ``end_residue``.
- Values must be positive where physically required (``temperature``,
  ``ewald_tolerance``, ``real_space_cutoff``, ``fugacity``); violations abort
  with an explicit error message.
