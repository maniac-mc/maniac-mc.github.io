#!/bin/bash
set -euo pipefail

# Path to LAMMPS executable (default: "lmp" in PATH)
LMP=${LMP:-lmp}


# Loop over everything that is a directory (H2O, CO2, ... later)
for testdir in */ ; do
    (
        cd "$testdir"

        # Every test directory must contain its own inner run-test.sh
        if [[ ! -x ./run-test.sh ]]; then
            echo "❌ ERROR: $testdir/run-test.sh is missing or not executable"
            exit 1
        fi

        # Export LMP so test scripts can use it
        LMP="$LMP" ./run-test.sh
    )
    echo "✅ Test passed: $testdir"
done
