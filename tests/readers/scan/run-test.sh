#!/bin/bash
set -euo pipefail

# Path to maniac files
BUILD_DIR="../../../build"
MODULE_DIR="../../../include"

# Path to test executables
BIN_DIR="./build"
mkdir -p "$BIN_DIR"

# Collect Maniac object files (excluding main.o)
shopt -s nullglob
MANIAC_OBJS=()
for obj in "$BUILD_DIR"/*.o; do
    [[ $(basename "$obj") == "main.o" ]] && continue
    MANIAC_OBJS+=("$obj")
done
shopt -u nullglob

# COMPILE & RUN TESTS
for TEST_SRC in test_*.f90; do
    [[ ! -f "$TEST_SRC" ]] && {
        echo "No test_*.f90 files found — nothing to run."
        exit 0
    }

    TEST_NAME="${TEST_SRC%.f90}"
    OUT_EXE="$BIN_DIR/$TEST_NAME"

    gfortran -O2 -I"$MODULE_DIR" -I"$BUILD_DIR" -c "$TEST_SRC" -o "$TEST_NAME.o"
    gfortran -O2 -o "$OUT_EXE" "$TEST_NAME.o" "${MANIAC_OBJS[@]}"

    if "$OUT_EXE"; then
        echo "✅ [PASS] $TEST_NAME"
    else
        echo "❌ [FAIL] $TEST_NAME (runtime error)"
        exit 1
    fi

    rm -f "$TEST_NAME.o"
done

rm -f fort.10 2>/dev/null || true
