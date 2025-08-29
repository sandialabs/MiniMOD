#!/bin/bash
#SBATCH --job-name=earlycoll_benchmark
#SBATCH --output=results/earlycoll_benchmark_%j.out
#SBATCH --error=results/earlycoll_benchmark_%j.err
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --account=FY140001

# Load necessary modules (if any)
# module load mpi

# Executable name
EXEC=earlycoll

# Operation to benchmark
OPERATION=alltoall

# Chunk size range
MIN_CHUNK_SIZE=4
MAX_CHUNK_SIZE=16

# Number of iterations
ITERATIONS=10

# Maximum number of nodes
MAX_NODES=$SLURM_NNODES

# Run the benchmark for an increasing number of nodes
for (( nodes=1; nodes<=MAX_NODES; nodes*=2 )); do
    echo "Running benchmark with $nodes nodes"
    srun --nodes=$nodes --ntasks-per-node=4 --cpus-per-task=1 $EXEC $OPERATION $MIN_CHUNK_SIZE $MAX_CHUNK_SIZE $ITERATIONS 1 0
done
