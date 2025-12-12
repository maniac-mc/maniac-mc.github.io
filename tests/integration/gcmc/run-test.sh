#!/bin/bash
set -euo pipefail

# ------------------------------------------------------------
#  Paths
# ------------------------------------------------------------
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../" && pwd)"

BUILD_DIR="$ROOT_DIR/build"
MODULE_DIR="$ROOT_DIR/include"
BIN_DIR="./build"

mkdir -p "$BIN_DIR"

# ------------------------------------------------------------
#  Collect Maniac object files (excluding main.o)
# ------------------------------------------------------------
shopt -s nullglob
MANIAC_OBJS=()
for obj in "$BUILD_DIR"/*.o; do
    [[ $(basename "$obj") == "main.o" ]] && continue
    MANIAC_OBJS+=("$obj")
done
shopt -u nullglob

if [[ ${#MANIAC_OBJS[@]} -eq 0 ]]; then
    echo "ERROR: No Maniac object files found in $BUILD_DIR"
    echo "Run:  (from project root)"
    echo "   mkdir -p build && cd build && cmake .. && make -j"
    exit 1
fi

# ------------------------------------------------------------
#  Compile & run each test_*.f90 file
# ------------------------------------------------------------
found_tests=false
for TEST_SRC in test_*.f90; do

    [[ ! -f "$TEST_SRC" ]] && continue
    found_tests=true

    TEST_NAME="${TEST_SRC%.f90}"
    OUT_EXE="$BIN_DIR/$TEST_NAME"

    # Compile test file
    gfortran -O2 -I"$MODULE_DIR" -I"$BUILD_DIR" -c "$TEST_SRC" -o "$TEST_NAME.o"

    # Link against Maniac core object files
    gfortran -O2 -o "$OUT_EXE" "$TEST_NAME.o" "${MANIAC_OBJS[@]}"

    # Run test
    if "$OUT_EXE" &>/dev/null; then
        echo "✅ [PASS] $TEST_NAME"
    else
        echo "❌ [FAIL] $TEST_NAME (runtime error)"
        exit 1
    fi

    rm -f "$TEST_NAME.o"
done

if ! $found_tests; then
    echo "No tests found (no test_*.f90 files)."
fi

rm -f fort.10 2>/dev/null || true