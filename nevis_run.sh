#!/bin/bash

# Move generated scripts to the root directory
# mv ./generated_scripts/spinup/* ./generated_scripts/output
# mv ./generated_scripts/drainage/* ./generated_scripts/output
# mv ./generated_scripts/output/* ./
# mv ./generated_scripts/convergence_tests/* ./
mv ./generated_scripts/case_studys/* ./

# --- Configuration ---
# Set the maximum number of parallel jobs. Default is 12.
MAX_PARALLEL_JOBS=${MAX_PARALLEL_JOBS:-9}

# --- Global Variables ---
# This variable will hold the PIDs of the running subshells and must be accessible by cleanup
running_pids=""

# --- Functions ---

# Cleanup function to kill all child processes on exit
cleanup() {
    echo -e "\nInterrupt received. Killing all background jobs..."
    
    # Terminate the specific MATLAB processes spawned by this script's subshells.
    # This is a more targeted and robust approach.
    if [ -n "$(echo $running_pids | xargs)" ]; then
        echo "Terminating running MATLAB jobs..."
        for pid_job in $running_pids; do
            subshell_pid=$(echo $pid_job | cut -d: -f1)
            job_name=$(echo $pid_job | cut -d: -f2)
            
            # Find the matlab process that is a child of our subshell
            matlab_pid=$(pgrep -P $subshell_pid)
            
            if [ -n "$matlab_pid" ]; then
                echo "Killing MATLAB process with PID $matlab_pid (from job $job_name)"
                # Kill the entire process group starting with the MATLAB process
                kill -9 -$matlab_pid 2>/dev/null || kill -9 $matlab_pid 2>/dev/null
            fi
            # Also kill the subshell itself
            kill -9 $subshell_pid 2>/dev/null
        done
    fi

    # As a final safeguard, a general pkill for any orphaned processes.
    echo "Final cleanup check for any remaining MATLAB processes..."
    pkill matlab

    echo "Cleanup complete."
    exit 1
}

# Function to run a single MATLAB script
run_matlab_script() {
    local script=$1
    echo "Launching ${script}.m"
    
    local start_time=$(date +%s)
    # By using 'exec', the shell process is replaced by the matlab process.
    # This makes PID tracking and cleanup much more reliable.
    exec matlab -batch "${script}" > "logs/${script}.log" 2>&1
    # The code from here on will not be executed because of 'exec'
}

# --- Main Script ---

# Set trap for interruption signals
trap cleanup INT TERM

echo "Starting NEVIS job runner..."

# Create logs directory
mkdir -p logs

# --- Task Queue Preparation ---
all_spinup_scripts=$(ls n[12]d*spinup.m 2>/dev/null | sed 's/\.m$//')
all_drainage_scripts=$(ls n[12]d*drainage.m 2>/dev/null | sed 's/\.m$//')
all_standalone_scripts=$(ls n[12]d*.m 2>/dev/null | grep -v -e "_spinup.m" -e "_drainage.m" | sed 's/\.m$//')

initial_queue="$all_standalone_scripts $all_spinup_scripts"
total_jobs=$(echo $initial_queue $all_drainage_scripts | wc -w)

if [ $total_jobs -eq 0 ]; then
    echo "No MATLAB scripts found. Exiting."
    exit 1
fi

# --- Job Execution Engine ---
pending_queue="$initial_queue"
completed_jobs=0

echo "Found a total of $total_jobs jobs to run."
echo "Running up to $MAX_PARALLEL_JOBS cases in parallel..."

while [ $completed_jobs -lt $total_jobs ]; do
    # --- 1. Launch new jobs if there are free slots ---
    num_running=$(echo $running_pids | wc -w)
    while [ $num_running -lt $MAX_PARALLEL_JOBS ] && [ -n "$(echo $pending_queue | xargs)" ]; do
        next_job=$(echo $pending_queue | awk '{print $1}')
        pending_queue=$(echo $pending_queue | awk '{$1=""; print $0}' | xargs)

        # Launch the job in a subshell that will 'exec' matlab
        (run_matlab_script "$next_job") &
        pid=$!
        running_pids="$running_pids $pid:$next_job"
        echo "Launched job: $next_job (PID: $pid)"
        num_running=$(echo $running_pids | wc -w)
    done

    # --- 2. Wait and check for completed jobs ---
    if [ -z "$(echo $running_pids | xargs)" ]; then
        if [ $completed_jobs -lt $total_jobs ]; then
             echo "Warning: Job execution stalled. Check for pending drainage cases whose spinups may have failed."
        fi
        break
    fi
    
    sleep 5

    # --- 3. Process finished jobs and update queues ---
    new_running_pids=""
    for pid_job in $running_pids; do
        pid=$(echo $pid_job | cut -d: -f1)
        job_name=$(echo $pid_job | cut -d: -f2)

        if kill -0 $pid 2>/dev/null; then
            new_running_pids="$new_running_pids $pid_job"
        else
            wait $pid
            exit_code=$?
            ((completed_jobs++))
            echo "Job finished: $job_name (Completed: $completed_jobs/$total_jobs)"
            # move the finished .m scripts to ./parameter_sweep
            # mv "n[12]d*${job_name}.m" "parameter_sweep/"

            if [[ "$job_name" == *spinup ]] && [ $exit_code -eq 0 ]; then
                base_name=${job_name%_spinup}
                matching_drainage_cases=$(echo "$all_drainage_scripts" | grep "^${base_name}" || true)
                
                if [ -n "$matching_drainage_cases" ]; then
                    echo "Spinup for $base_name succeeded. Queueing corresponding drainage case(s):"
                    echo "$matching_drainage_cases"
                    pending_queue="$matching_drainage_cases $pending_queue"
                fi
            elif [ $exit_code -ne 0 ]; then
                echo "Job $job_name failed with exit code $exit_code. Corresponding drainage case(s) will be skipped."
            fi
        fi
    done
    running_pids=$(echo $new_running_pids | xargs)
done

echo "All MATLAB scripts completed."