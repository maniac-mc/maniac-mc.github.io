#!/bin/bash
# lib.sh - shared helpers for MANIAC tests.
#
# This file is a LIBRARY: source it from any tests/**/run-test.sh, e.g.:
#   source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"        (from tests/units)
#   source "$(dirname "${BASH_SOURCE[0]}")/../../lib.sh"     (from tests/readers/<name>)
#   source "$(dirname "${BASH_SOURCE[0]}")/../../../lib.sh"  (from tests/integration/<name>)
#
# To try it manually from the directory that contains test_*.f90 files:
#   cd tests/units && source ../lib.sh && compile_and_run_all_tests

# Guard: refuse to run when executed directly instead of sourced
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "lib.sh is a library, not a script." >&2
    echo "Source it from a run-test.sh, or manually:" >&2
    echo "  cd tests/units && source ../lib.sh && compile_and_run_all_tests" >&2
    exit 1
fi

# Find the repository root (the directory containing Makefile and src/)
# by walking up from the sourcing script.
_maniac_find_repo_root() {
    local d
    d="$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)"
    while [[ "$d" != "/" ]]; do
        if [[ -f "$d/Makefile" && -d "$d/src" ]]; then
            printf '%s\n' "$d"
            return 0
        fi
        d="$(dirname "$d")"
    done
    return 1
}

MANIAC_REPO_ROOT="$(_maniac_find_repo_root)" || {
    echo "ERROR: could not locate repo root (Makefile + src/) above ${BASH_SOURCE[1]}" >&2
    exit 1
}
MANIAC_BUILD_DIR="$MANIAC_REPO_ROOT/build"
MANIAC_INCLUDE_DIR="$MANIAC_REPO_ROOT/include"

# compile_and_run_all_tests [glob]
#
# Compiles every file matching the glob (default: test_*.f90) found in the
# CURRENT directory, links it against all maniac objects except main.o,
# runs the resulting binary, and returns non-zero on the first failure
# (compilation, linking or runtime).
compile_and_run_all_tests() {
    local pattern="${1:-test_*.f90}"
    local bin_dir="./build"
    mkdir -p "$bin_dir"

    local objs=() obj
    for obj in "$MANIAC_BUILD_DIR"/*.o; do
        [[ -e "$obj" ]] || continue
        [[ "$(basename "$obj")" == "main.o" ]] && continue
        objs+=("$obj")
    done
    if [[ ${#objs[@]} -eq 0 ]]; then
        echo "ERROR: no object files in $MANIAC_BUILD_DIR - run ./build.sh first" >&2
        return 1
    fi

    local test_src test_name out_exe found=false
    for test_src in $pattern; do
        [[ -f "$test_src" ]] || continue
        found=true
        test_name="$(basename "$test_src" .f90)"
        out_exe="$bin_dir/$test_name"

        if ! gfortran -O2 -I"$MANIAC_INCLUDE_DIR" -I"$MANIAC_BUILD_DIR" -c "$test_src" -o "$test_name.o"; then
            echo "❌ [FAIL] $test_name (compilation error)"
            return 1
        fi
        if ! gfortran -O2 -o "$out_exe" "$test_name.o" "${objs[@]}"; then
            echo "❌ [FAIL] $test_name (linking error)"
            rm -f "$test_name.o"
            return 1
        fi
        if "$out_exe"; then
            echo "✅ [PASS] $test_name"
        else
            echo "❌ [FAIL] $test_name (runtime error)"
            rm -f "$test_name.o"
            return 1
        fi
        rm -f "$test_name.o"
    done

    if ! $found; then
        echo "ERROR: no files matched '$pattern' in $(pwd)" >&2
        return 1
    fi
    rm -f fort.10 2>/dev/null || true
    return 0
}
