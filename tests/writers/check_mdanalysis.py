#!/usr/bin/env python3
"""Independent LAMMPS-compatibility gate: parse the file with MDAnalysis.

Usage:
    python3 check_mdanalysis.py outputs/topology.data

MDAnalysis ships a strict LAMMPS data-file parser used by the whole MD
ecosystem, so if it reads the file, real LAMMPS-style tools will too.
We additionally cross-check the parsed counts against the file header.

Requires MDAnalysis (CI installs it in a venv). Exit 0 = OK, 1 = failure.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_output import parse_data_file


def main():
    if len(sys.argv) != 2:
        print("usage: check_mdanalysis.py FILE.data", file=sys.stderr)
        sys.exit(1)
    path = sys.argv[1]

    d = parse_data_file(path)
    h = d['header']

    try:
        import MDAnalysis as mda
    except ImportError:
        print(f"FAIL {path}: MDAnalysis is not installed", file=sys.stderr)
        sys.exit(1)

    try:
        u = mda.Universe(path)
    except Exception as exc:
        print(f"FAIL {path}: MDAnalysis could not parse the file: {exc}", file=sys.stderr)
        sys.exit(1)

    errs = []
    got = {
        'atoms': len(u.atoms),
        'bonds': len(u.bonds),
        'angles': len(u.angles),
        'dihedrals': len(u.dihedrals),
        'impropers': len(u.impropers),
        'atom types': len(set(u.atoms.types)),
    }
    for key, want in got.items():
        declared = h.get(key, 0)
        if want != declared:
            errs.append(f"MDAnalysis parsed {want} {key} but header declares {declared}")
    if len(u.atoms.masses) != len(u.atoms):
        errs.append("MDAnalysis: not every atom got a mass")

    if errs:
        for e in errs:
            print(f"FAIL {path}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
