#!/bin/bash
set -euo pipefail

BUILD_DIR="../../build"
INCLUDE_DIR="../../include"
BIN_DIR="./build"
mkdir -p "$BIN_DIR"

# Loop over all test_*.f90 files in this directory
for TEST_SRC in test_*.f90; do
    TEST_NAME=$(basename "$TEST_SRC" .f90)
    OUT_EXE="$BIN_DIR/$TEST_NAME"

    # Compile
    gfortran -O2 -I"$INCLUDE_DIR" -c "$TEST_SRC"

    # Gather all object files except main.o
    OBJS=$(ls "$BUILD_DIR"/*.o | grep -v "main.o" || true)

    # Link
    gfortran -O2 -o "$OUT_EXE" "$TEST_NAME.o" $OBJS

    # Run and report
    if "$OUT_EXE"; then
        echo "✅ [PASS] $TEST_NAME ran successfully"
    else
        echo "❌ [FAIL] $TEST_NAME crashed or returned non-zero"
        exit 1
    fi

    # Clean up test-specific object files
    rm -f "$TEST_NAME.o"
done

# Cleanup remaining build artifacts
rm -f fort.10
