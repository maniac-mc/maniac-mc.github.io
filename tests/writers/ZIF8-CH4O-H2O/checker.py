#!/usr/bin/env python
# coding: utf-8

from checker_utilities import read_lammps_data, check_bonds, check_angles, \
    check_dihedrals, check_impropers

atoms, bonds, angles, dihedrals, impropers = read_lammps_data("outputs/topology.data")
from parameters import allowed_bonds, allowed_angles, allowed_dihedrals, allowed_impropers

check_bonds(atoms, bonds, allowed_bonds)
check_angles(atoms, angles, allowed_angles)
check_dihedrals(atoms, dihedrals, allowed_dihedrals)
check_impropers(atoms, impropers, allowed_impropers)

