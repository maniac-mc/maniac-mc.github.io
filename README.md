# MANIAC

![Last Commit](https://img.shields.io/github/last-commit/maniac-mc/maniac-mc.github.io?color=%2328a745)
![License: MIT](https://img.shields.io/badge/License-MIT-%2328a745)
![GitHub release](https://img.shields.io/github/v/release/maniac-mc/maniac-mc.github.io?color=%2328a745)
[![Build](https://github.com/maniac-mc/maniac-mc.github.io/actions/workflows/tests.yml/badge.svg)](https://github.com/maniac-mc/maniac-mc.github.io/actions/workflows/tests.yml)
[![Docs](https://github.com/maniac-mc/maniac-mc.github.io/actions/workflows/docs.yml/badge.svg)](https://github.com/maniac-mc/maniac-mc.github.io/actions/workflows/docs.yml)
[![Lint](https://github.com/maniac-mc/maniac-mc.github.io/actions/workflows/lint.yml/badge.svg)](https://github.com/maniac-mc/maniac-mc.github.io/actions/workflows/lint.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24052849.svg)](https://doi.org/10.5281/zenodo.24052849)

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

The resulting executable is located at:

```bash
build/maniac
```

## Quickstart

Clone the repository together with the `mc-topology` submodule that contains
the example systems, build, and run the ZIF8-CH4O grand canonical example:

```bash
git clone --recurse-submodules https://github.com/maniac-mc/maniac-mc.github.io.git
cd maniac-mc.github.io
./build.sh
./run.sh
```

`./run.sh` executes the ZIF8-CH4O case and writes the results to `outputs/`.
To run another system, call the executable directly:

```bash
build/maniac -i input.maniac -d topology.data -p parameters.inc -o outputs/ -r reservoir.data
```

| Option | Meaning | Required |
| ------ | ------- | -------- |
| `-i` | input file (`.maniac`) | yes |
| `-d` | topology data file (LAMMPS `.data`) | yes |
| `-p` | parameter include file (`parameters.inc`) | yes |
| `-r` | reservoir topology (for insertion/swap moves) | no |
| `-o` | output directory (default: `outputs/`) | no |

The full list of input-file keywords is described in the
[input file documentation](https://maniac-mc.github.io/input.html).

## Output files

All results are written to the output directory (`outputs/` by default):

- `log.maniac`: run log
- `topology.data`: final configuration, in LAMMPS data format
- `trajectory.lammpstrj` (and `reservoir.lammpstrj`): trajectories in LAMMPS dump format
- `energy.dat`: energy as a function of the Monte Carlo step
- `moves.dat`: attempted/accepted move counters
- `number_<residue>.dat`, `histogram_<residue>.dat`, `widom_<residue>.dat`: per-residue observables

## LAMMPS compatibility

MANIAC uses the same `.data` file format as LAMMPS for molecular
topology. MANIAC assumes that the real units system is used, and that
a pair style of the family `lj/cut/coul/long` is used.

## Build the documentation

The documentation is built with Sphinx. Install the requirements:

```bash
pip install sphinx furo
```

Then, navigate to the `docs` directory and run:

```bash
cd docs
./build.sh
```

This reads the version from `version.txt` and creates the complete
documentation in the `docs/build/html` folder, which can be opened with:

```bash
firefox build/html/index.html
```

The documentation is also visible online at
[maniac-mc.github.io](https://maniac-mc.github.io).

## Examples and tests

Several example systems are provided in a separate repository:
[mc-topology](https://github.com/maniac-mc/mc-topology).
The code has been validated against LAMMPS and RASPA for the example cases
located in [mc-topology](https://github.com/maniac-mc/mc-topology).

The automated test suite (unit, reader, writer, and integration tests) can be
run with a single command:

```bash
bash tests/run-all.sh
```

See [tests/README.md](tests/README.md) for details.

## How to cite

If you use MANIAC in your research, please cite it via its Zenodo DOI:
[10.5281/zenodo.24052849](https://doi.org/10.5281/zenodo.24052849).

## Credit

This code was written by Simon Gravelle, who currently maintains the
code, documentation, and the associated [GitHub organization](https://github.com/maniac-mc).
