#!/bin/bash
set -e  # Exit on error

build_path="../../../../build/maniac"

input="input.maniac"
data="topology.data"
inc="parameters.inc"
outputs="outputs/" # optional

maniac_log="log.maniac"
$build_path -i "$input" -d "$data" -p "$inc" -o "$outputs" > "$maniac_log" 2>&1

# Theoretical energy
ref_energy=16.167421

tolerance=0.001 # allowed deviation

# Get the last TotEng value from the log
last_toteng=$(awk '/TotEng/ {getline; gsub(/^[| ]+/,""); print $2; exit}' "$maniac_log")

if [ -z "$last_toteng" ]; then
    echo "❌ Test failed: 'TotEng' not found in $maniac_log"
    exit 1
fi

# Compute absolute difference
diff=$(awk -v a="$last_toteng" -v b="$ref_energy" 'BEGIN { print (a > b ? a - b : b - a) }')

# Check tolerance
if (( $(echo "$diff > $tolerance" | bc -l) )); then
    echo "❌ Test failed: TotEng=$last_toteng differs from ref=$ref_energy by $diff"
    exit 1
else
    echo "✅ Test passed: TotEng=$last_toteng (diff=$diff within tolerance $tolerance)"
fi
