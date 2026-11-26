#!/bin/bash
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # tests/writers
REPO_ROOT="$(cd "$HERE/../.." && pwd)"
BUILD="$REPO_ROOT/build/maniac"
TOPO="$REPO_ROOT/mc-topology"
CHECK_OUT="$HERE/check_output.py"
CHECK_MDA="$HERE/check_mdanalysis.py"

if [[ ! -x "$BUILD" ]]; then
    echo "ERROR: $BUILD not found - run ./build.sh first" >&2
    exit 1
fi

failures=0

# run_case NAME TOPOLOGY_DIR INPUT_MANIAC [CHECKER_DIR]
run_case() {
    local name="$1" topo_dir="$2" input="$3" checker_dir="${4:-}"
    local work
    work="$(mktemp -d)"
    if ! ( cd "$work" && "$BUILD" -i "$input" -d "$topo_dir/topology.data" \
              -p "$topo_dir/parameters.inc" -o outputs > log.maniac 2>&1 ); then
        echo "❌ $name: maniac run failed (artifacts kept at $work)"
        failures=$((failures + 1)); return
    fi
    if ! python3 "$CHECK_OUT" "$work/outputs/topology.data"; then
        echo "❌ $name: structural check failed (artifacts kept at $work)"
        failures=$((failures + 1)); return
    fi
    if python3 -c "import MDAnalysis" >/dev/null 2>&1; then
        if ! python3 "$CHECK_MDA" "$work/outputs/topology.data"; then
            echo "❌ $name: MDAnalysis check failed (artifacts kept at $work)"
            failures=$((failures + 1)); return
        fi
    else
        echo "SKIP MDAnalysis check: MDAnalysis not installed"
    fi
    if [[ -n "$checker_dir" ]]; then
        if ! ( cd "$work" && PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$checker_dir" \
                  python3 "$checker_dir/checker.py" ); then
            echo "❌ $name: connectivity pattern check failed (artifacts kept at $work)"
            failures=$((failures + 1)); return
        fi
    fi
    rm -rf "$work"
    echo "✅ $name"
}

run_case H2O            "$TOPO/testcase-adsorption/H2O"           "$TOPO/testcase-adsorption/H2O/input.maniac"
run_case LJ-gas         "$TOPO/testcase-energy/LJ-gas"             "$TOPO/testcase-energy/LJ-gas/input.maniac"
run_case ZIF8-CH4O-H2O  "$TOPO/testcase-adsorption/ZIF8-CH4O-H2O"  "$HERE/ZIF8-CH4O-H2O/input.maniac" "$HERE/ZIF8-CH4O-H2O"

if [[ $failures -gt 0 ]]; then
    echo "❌ $failures writer test(s) FAILED"
    exit 1
fi
echo "✅ All writer tests passed"
