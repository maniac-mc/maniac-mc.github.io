#!/bin/bash
# run-all.sh - run every MANIAC test suite and report a global result.
#
# Usage:  tests/run-all.sh        (from anywhere; requires ./build.sh to have run)
#         LMP=/path/to/lmp tests/run-all.sh   to point at a LAMMPS binary
#
# Exits non-zero if any suite fails.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if [[ ! -x build/maniac ]]; then
    echo "ERROR: build/maniac not found - run ./build.sh first" >&2
    exit 1
fi

failures=0

run_suite() {
    local dir="$1"
    echo "================ $dir ================"
    if ( cd "$dir" && ./run-test.sh ); then
        echo "✅ suite OK: $dir"
    else
        echo "❌ suite FAILED: $dir"
        failures=$((failures + 1))
    fi
}

# Compiled Fortran unit tests
run_suite tests/units

# Reader tests (parsing of .maniac / .data files)
for d in tests/readers/*/; do
    run_suite "${d%/}"
done

# Writer tests (require LAMMPS)
if command -v "${LMP:-lmp}" >/dev/null 2>&1; then
    run_suite tests/writers
else
    echo "⚠️  LAMMPS ('${LMP:-lmp}') not found - skipping tests/writers"
fi

# Integration tests: end-to-end runs of the maniac binary
for d in tests/integration/energy/*/ tests/integration/energy-analytical/*/; do
    run_suite "${d%/}"
done
run_suite tests/integration/gcmc
run_suite tests/integration/energy-drift

# CLI tests
run_suite tests/cli

echo
if [[ $failures -gt 0 ]]; then
    echo "❌ $failures test suite(s) FAILED"
    exit 1
fi
echo "✅ All test suites passed"
