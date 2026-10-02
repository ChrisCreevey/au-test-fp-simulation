#!/bin/bash
#SBATCH --job-name=retry_le200
#SBATCH --partition=k2-hipri
#SBATCH --array=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/retry_le200/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/retry_le200/slurm_%A_%a.err

# Trivial final mop-up for the le200 tranche: 8 families left, in one
# chunk -- this size class has a strong track record (435s/family
# average), no special handling needed.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
CHUNK_FILE="${DATA_DIR}/manifest/chunks_le200_retry/chunk_0001.txt"

cd "$PROJECT_DIR"
echo "Host: $(hostname), allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family_chunk.sh "$CHUNK_FILE" "$MANIFEST" "$OUT_DIR" config/pipeline.conf
