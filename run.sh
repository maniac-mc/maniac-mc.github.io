#!/bin/bash
set -e  # Exit on error

case="SLIT"

base_energy="mc-topology/testcase-energy"
base_adsorption="mc-topology/testcase-adsorption"
base_reservoir="mc-topology/molecule-reservoir"
base_widom="mc-topology/testcase-widom"

case "$case" in
  "ZIF8-CH4O")
    folder="$base_adsorption/ZIF8-CH4O"
    reservoir="$base_reservoir/CH4O-H2O"
    ;;

  "WIDOM")
    folder="$base_widom/ZIF8-MET"
    reservoir="$base_reservoir/CH4O"
    ;;

  "CO2")
    folder="$base_reservoir/CO2"
    reservoir="$base_reservoir/CO2"
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

# If reservoir was set, define res
if [[ -n "$reservoir" ]]; then
    res="$reservoir/topology.data"
fi

outputs="outputs/"

echo "Running maniac with system=$case ..."
./build/maniac -i "$input" -d "$data" -p "$inc" -o "$outputs" ${res:+-r "$res"}
