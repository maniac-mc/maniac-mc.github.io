#!/bin/bash
set -euo pipefail

topology_path="../../../mc-topology/testcase-adsorption/H2O/"
build_path="../../../build/maniac"

inc=$topology_path"parameters.inc"
input=$topology_path"input.maniac"
outputs="outputs/"
data=$topology_path"topology.data"

$build_path -i "$input" -d "$data" -p "$inc" -o "$outputs" > log.maniac 2>&1

# Allow GitHub Actions to override LAMMPS executable
lmp=${LMP:-lmp}

${lmp} -in input.lmp > log.lammps 2>&1
