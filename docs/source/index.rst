MANIAC documentation
====================

**MANIAC** is a lightweight Monte Carlo simulation code written in modern Fortran.
It is designed for grand canonical Monte Carlo (GCMC) simulations used for constant
chemical potential calculations, particularly useful in adsorption studies. The main
difference between MANIAC and other Monte Carlo codes is that it reads basic
LAMMPS-style topology files.

**MANIAC** supports the following Monte Carlo moves:

- **Translation -** A molecule is randomly selected, and a translational move of
  random distance (shorter than a specified maximum value) and random direction is
  attempted. The move is accepted according to the Metropolis criterion:

  .. math::
    p_{\rm accept} = \min \left( 1, \mathrm{e}^{-\beta \Delta U} \right)

  where :math:`\Delta U` is the change in potential energy due to the move and
  :math:`\beta = 1 / (k_\mathrm{B} T)`.

- **Rotation -** A molecule is randomly chosen, and a random rotational move around
  a randomly selected axis is attempted, within a defined maximum rotation angle.
  The move is accepted using the same Metropolis criterion:

  .. math::
    p_{\rm accept} = \min \left( 1, \mathrm{e}^{-\beta \Delta U} \right)

- **Insertion -** A new molecule is randomly inserted into the simulation box at a
  random position and orientation. The move is accepted according to the grand
  canonical acceptance criterion:

  .. math::
    p_{\rm accept} = \min \left( 1, \frac{V}{\Lambda^3 (N+1)}
    \mathrm{e}^{\beta \mu} \mathrm{e}^{-\beta \Delta U} \right)

  where :math:`V` is the box volume, :math:`\Lambda` is the thermal wavelength,
  :math:`N` is the number of molecules before insertion, :math:`\mu` is the
  chemical potential of the reservoir, and :math:`\Delta U` is the energy change.

- **Deletion -** An existing molecule is randomly selected and removed from the
  system. The acceptance probability is:

  .. math::
    p_{\rm accept} = \min \left( 1, \frac{\Lambda^3 N}{V}
    \mathrm{e}^{-\beta \mu} \mathrm{e}^{-\beta \Delta U} \right).

- **Swap -** A molecule of type A is randomly selected and replaced with a molecule
  of type B at the same position with a random orientation. In the semi-grand canonical
  ensemble, the acceptance criterion for swapping type A to B is:

  .. math::
    p_{\rm accept} = \min \left( 1, \frac{N_A}{N_B + 1}
    \mathrm{e}^{-\beta (\Delta U + \mu_A - \mu_B)} \right)

  where :math:`N_A` is the number of molecules of type A before the swap,
  :math:`N_B` is the number of molecules of type B, :math:`\Delta U` is the energy
  change, and :math:`\beta = 1 / (k_\mathrm{B} T)`, and where :math:`\mu_A` and
  :math:`\mu_B` are the chemical potentials of types A and B.

Why the name **MANIAC**?
------------------------

The name **MANIAC**, or **MANIAC-MC**, is an homage to the original MANIAC computer (Mathematical
Analyzer, Numerical Integrator, and Computer), built in the early 1950s at Los
Alamos National Laboratory. This pioneering machine was among the first to perform
Monte Carlo simulations for statistical physics and nuclear research. The
**MC** in **MANIAC** reflects its core use of Monte Carlo methods.

LAMMPS compatibility
--------------------

**MANIAC** uses the same `.data` file format as LAMMPS for molecular
topology. MANIAC assumes that the real units system is used, and that
a pair style of the family `lj/cut/coul/long` is used.

.. toctree::
    :maxdepth: 2
    :caption: Contents:
    :hidden:

    input
    build
    unit
    ewald
    credit
