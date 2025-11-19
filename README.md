# MANIAC

![Last Commit](https://img.shields.io/github/last-commit/maniac-mc/maniac-mc.github.io?color=%2328a745)
![License: MIT](https://img.shields.io/badge/License-MIT-%2328a745)
![GitHub release](https://img.shields.io/github/v/release/maniac-mc/maniac-mc.github.io?color=%2328a745)
[![Build](https://github.com/maniac-mc/maniac-mc.github.io/actions/workflows/tests.yml/badge.svg)](https://github.com/maniac-mc/maniac-mc.github.io/actions/workflows/tests.yml)
[![Docs](https://github.com/maniac-mc/maniac-mc.github.io/actions/workflows/docs.yml/badge.svg)](https://github.com/maniac-mc/maniac-mc.github.io/actions/workflows/docs.yml)
[![Lint](https://github.com/maniac-mc/maniac-mc.github.io/actions/workflows/lint.yml/badge.svg)](https://github.com/maniac-mc/maniac-mc.github.io/actions/workflows/lint.yml)

<img
    src="https://raw.githubusercontent.com/maniac-mc/mc-visuals/refs/heads/main/gallery/ZIF8-H2O/system.png"
    width="30%" align="right"/>
</a>

MANIAC is a lightweight Monte Carlo simulation code written in Fortran,
designed for GCMC and adsorption studies. It reads basic LAMMPS-style topology
files and supports the following Monte Carlo moves:

- Translation  
- Rotation  
- Insertion  
- Deletion  
- Swap  

## Why the name MANIAC?

The original MANIAC computer (for Mathematical Analyzer, Numerical Integrator, and
Computer) was built in the early 1950s at Los Alamos National Laboratory. It
was one of the first machines used to perform Monte Carlo simulations in
statistical physics and nuclear research.

## Compile

To build MANIAC from source, you will need the following tools installed:

```bash
sudo apt update
sudo apt install gfortran make makedepf90
```

Then compile the project by running:

```bash
./build.sh
```

This will:

- Generate the Fortran version module from the template.
- Clean previous build artifacts.
- Rebuild the dependency file using makedepf90.
- Compile all Fortran source files with gfortran.

The resulting executable will be located at:

```bash
build/maniac
```

## LAMMPS compatibility

MANIAC uses the same `.data` file format as LAMMPS for molecular
topology. MANIAC assumes that the real units system is used, and that
a pair style of the family `lj/cut/coul/long` is used.

## Build the documentation

To build the documentation locally, install doxygen:

```bash
sudo apt install doxygen
sudo apt install python3-sphinx
pip install furo
```

Then, navigate to the `docs` directory and run:

```bash
cd docs
./build.sh
```

This will create the complete documentation in the `build/html` folder, that can
be open using:

```bash
firefox build/html/index.html
```

The documentation is also visible on [maniac-mc.github.io](https://maniac-mc.github.io).

## Examples and tests

Several example systems are provided in a separate repository:
[mc-topology]([topology-gallery/](https://github.com/maniac-mc/mc-topology)).
The code has been validated against LAMMPS and RASPA for several example cases
located in [mc-topology]([topology-gallery/](https://github.com/maniac-mc/mc-topology)).
Basic adsorption and energy tests are available in the [tests/](tests/) folder.

## Credit

This code was written by Simon Gravelle, who currently maintains the
code, documentation, and the associated [GitHub organization](https://github.com/maniac-mc).
