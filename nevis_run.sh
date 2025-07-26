#!/bin/bash

# Create logs and parameter_sweep directories
mkdir -p logs parameter_sweep

# Automatically find all MATLAB scripts starting with "n2d"
scripts=()
for file in n2d*.m n1d*.m; do
    if [ -e "$file" ]; then
        scripts+=("${file%.m}")
    fi
done

# Check if any scripts were found
if [ ${#scripts[@]} -eq 0 ]; then
    echo "No MATLAB scripts starting with 'n2d' found."
    exit 1
fi

# Sort the scripts array to make sure spinup scripts are run first
for script in "${scripts[@]}"; do
    if [[ "$script" == *"spinup"* ]]; then
        sorted_scripts+=("$script")
    fi
done
for script in "${scripts[@]}"; do
    if [[ "$script" != *"spinup"* ]]; then
        sorted_scripts+=("$script")
    fi
done
# print sorted scripts
echo "Sorted MATLAB scripts to run:"
for script in "${sorted_scripts[@]}"; do
    echo "  - ${script}.m"
done
# Loop through and run each MATLAB script
for script in "${sorted_scripts[@]}"
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