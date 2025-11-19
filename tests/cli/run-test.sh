#!/bin/bash

topology_path="../../mc-topology/testcase-energy/H2O-gas/"
reservoir_path="../../mc-topology/molecule-reservoir/tip4p-water/"
maniac="../../build/maniac"

input=$topology_path"input.maniac"
data=$topology_path"topology.data"
inc=$topology_path"parameters.inc"
reservoir=$reservoir_path"water.data"
outputs="outputs/"

# Color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Helper function to print pass/fail
pass() { echo -e "${GREEN}✅ [PASS] $1${NC}"; }
fail() { echo -e "${RED}❌ [FAIL] $1${NC}"; }

# ---------------------------
# Test 1: Basic valid CLI call
# ---------------------------
$maniac -i $input -d $data -p $inc -o $outputs > log.maniac 2>&1
if [ $? -eq 0 ]; then pass "Basic CLI"; else fail "Basic CLI"; fi

# --------------------------------
# Test 2: Missing mandatory flag -i
# --------------------------------
if $maniac -d $data -p $inc -o $outputs > log.maniac 2>&1 ; then
    fail "Missing -i not detected"
else
    grep -q "Missing mandatory input arguments" log.maniac && pass "Missing -i detected"
fi

# ---------------------------
# Test 3: Unknown flag
# ---------------------------
if $maniac -i $input -d $data -p $inc -x unknown > log.maniac 2>&1; then
    fail "Unknown flag not detected"
else
    grep -q "Unknown option" log.maniac && pass "Unknown flag detected"
fi

# ---------------------------
# Test 4: Missing value after flag
# ---------------------------
if $maniac -i -d $data -p $inc > log.maniac 2>&1; then
    fail "Missing value after -i not detected"
else
    grep -q "Missing value after option -i" log.maniac && pass "Missing value after -i detected"
fi

# ---------------------------
# Test 5: Duplicate flag
# ---------------------------
if $maniac -i $input -i $input -d $data -p $inc > log.maniac 2>&1; then
    fail "Duplicate -i flag not detected"
else
    grep -q "Duplicate option: -i" log.maniac && pass "Duplicate -i flag detected"
fi

# ---------------------------
# Test 6: Optional -r flag
# ---------------------------
if $maniac -i $input -d $data -p $inc -r $reservoir -o $outputs > log.maniac 2>&1; then
    pass "Optional -r flag"
else
    fail "Optional -r flag"
fi

# ---------------------------
# Test 7: Output path normalization
# ---------------------------
$maniac -i $input -d $data -p $inc -o "myoutput" > log.maniac 2>&1
grep -q "myoutput/" log.maniac && pass "Output path normalized" || fail "Output path normalization"

# ---------------------------
# Test 8: Help message
# ---------------------------
if $maniac --help > log.maniac 2>&1; then
    grep -q "Usage" log.maniac && pass "Help message displayed" || fail "Help message not displayed"
fi

# Delete extra folder 
rm -rf myoutput
rm -rf outputs