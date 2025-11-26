#!/bin/bash
set -e  # Exit on error

case="WIDOM"

base_energy="mc-topology/testcase-energy"
base_adsorption="mc-topology/testcase-adsorption"
base_reservoir="mc-topology/molecule-reservoir"
base_widom="mc-topology/testcase-widom"

case "$case" in
  "ZIF8-CH4O")
    folder="$base_adsorption/ZIF8-CH4O"
    ;;

  "WIDOM")
    folder="$base_widom/ZIF8-MET"
    ;;

  "SLIT")
    folder="$base_adsorption/SLIT"
    ;;

  *)
    echo "Unknown case: $case"
    exit 1
    ;;
esac

# --- Shared filenames (automatic for all cases) ---
input="$folder/input.maniac"
data="$folder/topology.data"
inc="$folder/parameters.inc"

outputs="outputs/"

echo "Running maniac with system=$case ..."
./build/maniac -i "$input" -d "$data" -p "$inc" -o "$outputs" ${res:+-r "$res"}
