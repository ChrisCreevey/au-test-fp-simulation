#!/bin/bash
#SBATCH --job-name=retry_201_500
#SBATCH --partition=k2-hipri
#SBATCH --array=1-175
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/retry_201_500/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/retry_201_500/slurm_%A_%a.err

# Retry pass for the 201-500 taxa tier's "likely to finish" families: 523
# families (never-started + killed-with-<=1h-remaining + Domain1-already-
# done), explicitly EXCLUDING the 15 killed with >1h still remaining
# (real evidence they're slow) -- see project history for the full
# evidence-based categorisation.
#
# Chunk size 3/chunk (vs the first attempt's 15, which had 0/42 chunks
# finish), grounded in real timing data from the 89 families that DID
# complete last round: mean 2803s, max 8693s. 3x mean = ~2.3h, leaving
# margin even if one chunk member lands near the observed max.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
CHUNK_FILE=$(printf "%s/manifest/chunks_201_500_retry/chunk_%04d.txt" "$DATA_DIR" "$SLURM_ARRAY_TASK_ID")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family_chunk.sh "$CHUNK_FILE" "$MANIFEST" "$OUT_DIR" config/pipeline.conf
