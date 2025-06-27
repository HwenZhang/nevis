#!/bin/bash

# original_file="nevis_regional_blister.m"
original_file="nevis_1d_expanding_blister_const_influx.m"
original_file="nevis_2d_example_hewitt_geometry.m"

# Debug: Print lines that mention oo.casename
echo "Looking for oo.casename in $original_file..."
grep "oo.casename" "$original_file"

# Improved sed pattern: handles arbitrary spaces/tabs
casename=$(sed -nE "s/^[[:space:]]*oo\.casename[[:space:]]*=[[:space:]]*['\"]([^'\"]+)['\"].*/\1/p" "$original_file")

# Check and copy
if [ -z "$casename" ]; then
    echo "Error: Could not find oo.casename in $original_file"
    exit 1
fi

new_file="${casename}.m"
cp "$original_file" "$new_file"
echo "Copied $original_file to $new_file"