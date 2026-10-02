#!/bin/bash
#SBATCH --job-name=retry_501_1000
#SBATCH --partition=k2-hipri
#SBATCH --array=1-121
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/retry_501_1000/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/retry_501_1000/slurm_%A_%a.err

# Retry pass for the 501-1000 taxa tier's "likely to finish" families: 242
# families (never-started + killed-with-<=1h-remaining + Domain1-already-
# done), excluding the 22 killed with >1h still remaining.
#
# Chunk size 2/chunk (smaller than 201-500's 3, per this tier's higher
# cost): real timing data from the 15 families that completed last round
# gave mean 4489s, max 9113s. 2x mean = ~2.5h.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
CHUNK_FILE=$(printf "%s/manifest/chunks_501_1000_retry/chunk_%04d.txt" "$DATA_DIR" "$SLURM_ARRAY_TASK_ID")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family_chunk.sh "$CHUNK_FILE" "$MANIFEST" "$OUT_DIR" config/pipeline.conf
