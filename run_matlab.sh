#!/bin/bash

# Create logs and parameter_sweep directories
mkdir -p logs parameter_sweep

# List all MATLAB scripts to run
scripts=("nreg_0mm_cg0_00_a0_1_kh0_ks1_mu5e0_c1_V0e8")

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