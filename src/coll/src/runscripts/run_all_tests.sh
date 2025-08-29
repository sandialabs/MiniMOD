#!/bin/bash

# Set the project root directory
PROJECT_ROOT=$(dirname $(dirname $(realpath $0)))

# Change to the project root directory
pushd $PROJECT_ROOT

# Remove existing module and mpicc path files
rm -f runscripts/loaded_modules.txt runscripts/mpicc_path.txt

# Cancel all current jobs (not neccicary but keeping everything clean for now)
scancel -u $USER

# Save the loaded modules
module list 2>&1 | tee runscripts/output/loaded_modules.txt

# Save the mpicc being used
which mpicc | tee runscripts/output/mpicc_path.txt

# Build the tests
make

# Create the output directory if it doesn't exist
mkdir -p runscripts/output

# Remove existing output files
rm -f runscripts/output/*.out runscripts/output/*.err

# Return to the original directory
popd

# Get the name of this script
this_script=$(basename "$0")

# Submit all .sh files in the runscripts directory to SLURM and capture job IDs
job_ids=()
for script in *.sh; do
    if [ "$(basename "$script")" != "$this_script" ]; then
        job_id=$(sbatch $script | awk '{print $4}')
        job_ids+=($job_id)
    fi
done

# Function to print the status of jobs
print_status() {
    echo "Current job status for user mdosanj:"
    squeue -u mdosanj
}

# Wait for all jobs to finish
for job_id in "${job_ids[@]}"; do
    while squeue -j $job_id &> /dev/null; do
        print_status
        sleep 10 # Check status every 10 seconds
    done
done

echo "All batch jobs have completed."

# Check the output files for test results
passed=0
failed=0

for output_file in runscripts/output/*.out; do
    if grep -q "verification passed" "$output_file" && ! grep -q "verification failed" "$output_file"; then
        ((passed++))
    else
        ((failed++))
    fi
done

# Print the summary
echo "Summary of test results:"
echo "Passed: $passed"
echo "Failed: $failed"

if [ $failed -gt 0 ]; then
    echo "Some tests failed. Please check the output files for details."
fi

cat output/test_allgather*
