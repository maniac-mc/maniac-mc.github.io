#!/bin/bash
set -euo pipefail

topology_path="../../../mc-topology/molecule-reservoir/CH4O-H2O/"
build_path="../../../build/maniac"

inc=$topology_path"parameters.inc"
input="input.maniac"
outputs="outputs/"
data=$topology_path"topology.data"

$build_path -i "$input" -d "$data" -p "$inc" -o "$outputs" > log.maniac 2>&1

lmp=${LMP:-lmp}

${lmp} -in input.lmp > log.lammps 2>&1
