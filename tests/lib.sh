#!/bin/bash
# lib.sh - shared helpers for MANIAC tests.
#
# This file is a LIBRARY: source it from any tests/**/run-test.sh, e.g.:
#   source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"     (from tests/units)
#   source "$(dirname "${BASH_SOURCE[0]}")/../../lib.sh"  (from tests/readers/<name>, tests/integration/<name>)
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
#
# All build artifacts (objects, executables, logs) go to a temporary
# scratch directory that is removed when the function returns, so test
# directories stay clean. Test output is captured to a log inside that
# scratch dir and shown only on failure.
compile_and_run_all_tests() {
    local pattern="${1:-test_*.f90}"
    local scratch rc
    scratch="$(mktemp -d "${TMPDIR:-/tmp}/maniac-test.XXXXXX")"
    _maniac_run_tests "$pattern" "$scratch"
    rc=$?
    rm -rf "$scratch"
    return "$rc"
}

_maniac_run_tests() {
    local pattern="$1" bin_dir="$2"

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

    local test_src test_name out_exe log_file found=false
    for test_src in $pattern; do
        [[ -f "$test_src" ]] || continue
        found=true
        test_name="$(basename "$test_src" .f90)"
        out_exe="$bin_dir/$test_name"
        log_file="$bin_dir/$test_name.log"

        if ! gfortran -O2 -I"$MANIAC_INCLUDE_DIR" -I"$MANIAC_BUILD_DIR" -c "$test_src" -o "$bin_dir/$test_name.o"; then
            echo "❌ [FAIL] $test_name (compilation error)"
            return 1
        fi
        if ! gfortran -O2 -o "$out_exe" "$bin_dir/$test_name.o" "${objs[@]}"; then
            echo "❌ [FAIL] $test_name (linking error)"
            return 1
        fi
        if "$out_exe" > "$log_file" 2>&1; then
            echo "✅ [PASS] $test_name"
        else
            echo "❌ [FAIL] $test_name (runtime error) - output:"
            cat "$log_file"
            return 1
        fi
    done

    if ! $found; then
        echo "ERROR: no files matched '$pattern' in $(pwd)" >&2
        return 1
    fi
    rm -f fort.10 2>/dev/null || true   # fort.10 lands in the current directory
    return 0
}

# assert_toteng_within <maniac_log> <ref_energy> <tolerance>
#
# Extracts the last TotEng value from a maniac log and returns non-zero if
# it differs from <ref_energy> by more than <tolerance>. The comparison is
# done in awk on purpose: bc cannot parse scientific notation (e.g. 4.8e-05)
# and would silently treat any comparison as false.
assert_toteng_within() {
    local log_file="$1" ref="$2" tol="$3"
    local last
    last="$(awk '/TotEng/ {getline; gsub(/^[| ]+/,""); print $2; exit}' "$log_file")"
    if [[ -z "$last" ]]; then
        echo "❌ TotEng not found in $log_file" >&2
        return 1
    fi
    if awk -v a="$last" -v b="$ref" -v t="$tol" 'BEGIN { d = (a > b ? a - b : b - a); exit !(d <= t) }'; then
        awk -v a="$last" -v b="$ref" -v t="$tol" 'BEGIN { d = (a > b ? a - b : b - a); printf "✅ TotEng=%s (diff=%g, tolerance %g)\n", a, d, t }'
        return 0
    else
        awk -v a="$last" -v b="$ref" -v t="$tol" 'BEGIN { d = (a > b ? a - b : b - a); printf "❌ TotEng=%s differs from ref=%s by %g (tolerance %g)\n", a, b, d, t }' >&2
        return 1
    fi
}
