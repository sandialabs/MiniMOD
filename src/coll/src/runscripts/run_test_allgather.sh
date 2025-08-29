#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=test_allgather
#SBATCH --output=output/test_allgather_%j.out
#SBATCH --error=output/test_allgather_%j.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:01:00
#SBATCH --account=FY140001
#SBATCH --licenses=pscratch,tscratch

# Load necessary modules (if any)

#export OMPI_MCA_osc_base_verbose=100

export OMPI_MCA_osc=rdma,sm,ucx

# Executable name
EXEC=../test_allgather

# Run the All-to-All unit tests
echo "Running Direct..."
#srun $EXEC direct
echo "Running Recoursive Doubling..."
#srun $EXEC recursive_doubling
echo "Running Bruck..."
#srun $EXEC bruck
echo "Running Direct Datacopy..."
#srun $EXEC datacopy
echo "Running Recursive Doubling Direct Copy with Rounds..."
srun $EXEC rd_dc_rnds
echo "Running Bruck Direct Copy with Rounds..."
#srun $EXEC b_dc_rnds
