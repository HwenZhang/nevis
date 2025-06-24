#!/bin/bash

# Create logs and parameter_sweep directories
mkdir -p logs parameter_sweep

# Automatically find all MATLAB scripts starting with "n2d"
scripts=()
for file in n2d*.m; do
    if [ -e "$file" ]; then
        scripts+=("${file%.m}")
    fi
done

# Check if any scripts were found
if [ ${#scripts[@]} -eq 0 ]; then
    echo "No MATLAB scripts starting with 'n2d' found."
    exit 1
else
    # Output all found script names before execution
    echo "Found MATLAB scripts to run:"
    for script in "${scripts[@]}"; do
        echo "  - ${script}.m"
    done
fi

# Loop through and run each MATLAB script
for script in "${scripts[@]}"
do
    echo "Launching ${script}"
    matlab -batch "${script}" > "./logs/temp_${script}.log" 2>&1

    # After the script finishes, move its .m file to the parameter_sweep directory
    if [ -f "${script}.m" ]; then
        mv "${script}.m" parameter_sweep/
    else
        echo "Warning: ${script}.m file not found, skipping move."
    fi
done

echo "All MATLAB scripts completed."