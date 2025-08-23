#!/bin/bash

mv ./generated_scripts/spinup/* ./generated_scripts/output
mv ./generated_scripts/drainage/* ./generated_scripts/output
mv ./generated_scripts/output/* ./

# count the number of spinup cases
num_spinup=$(ls n2d*spinup*.m 2>/dev/null | wc -l)
echo "Number of spinup cases: $num_spinup"
# set the number of parallel jobs to the number of spinup cases, if less than 16
parallel_jobs=${parallel_jobs:-$num_spinup}
parallel_jobs=$(( parallel_jobs < 12 ? parallel_jobs : 12 ))
    
cleanup() {
    echo -e "\nInterrupt received. Killing all processes in this group (PGID: $$)..."
    # send SIGTERM to every process in our process group
    kill 0
    echo "Cleanup complete."
    exit 1
}

trap cleanup INT TERM

echo "Running up to $parallel_jobs cases in parallel ..."
# Create logs and parameter_sweep directories
mkdir -p logs parameter_sweep

# Find all MATLAB scripts starting with "n2d" or "n1d"
scripts=()
for file in n2d*.m n1d*.m; do
    [ -e "$file" ] && scripts+=("${file%.m}")
done
[ ${#scripts[@]} -ne 0 ] || { echo "No MATLAB scripts found."; exit 1; }

# Sort spinup first, then drainage
sorted_scripts=()
for s in "${scripts[@]}"; do [[ "$s" == *spinup* ]] && sorted_scripts+=("$s"); done
for s in "${scripts[@]}"; do [[ "$s" != *spinup* ]] && sorted_scripts+=("$s"); done

echo "Running up to $parallel_jobs cases in parallel ..."
echo "Scripts to run:"
printf '  - %s.m\n' "${sorted_scripts[@]}"

# Loop and launch jobs in background, throttled by $parallel_jobs
for script in "${sorted_scripts[@]}"; do
    echo "Launching ${script}.m"
    (
        start_time=$(date +%s)
        matlab -batch "${script}" > "logs/temp_${script}.log" 2>&1
        end_time=$(date +%s)
        elapsed=$(( end_time - start_time ))
        echo "Case ${script} completed in ${elapsed} seconds"
        echo "Case ${script} completed in ${elapsed} seconds" >> "logs/temp_${script}.log"
        # move .m file to parameter_sweep
        # [ -f "${script}.m" ] && mv "${script}.m" parameter_sweep/
    ) &

    # if we’ve reached the parallel_jobs limit, wait for at least one to finish
    while [ "$(jobs -rp | wc -l)" -ge "$parallel_jobs" ]; do
        sleep 1
    done
done

# wait for all background jobs to finish
wait
echo "All MATLAB scripts completed."