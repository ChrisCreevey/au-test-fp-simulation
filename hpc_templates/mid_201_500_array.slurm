#!/bin/bash
#SBATCH --job-name=mid_201_500
#SBATCH --partition=k2-hipri
#SBATCH --array=1-42
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/mid_201_500/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/mid_201_500/slurm_%A_%a.err

# Second bulk run: the 627 families where the largest domain is 201-500
# taxa (the tier we're keeping ourselves; >500 taxa goes to Peter -- see
# docs/ADAPTING_TO_YOUR_CLUSTER.md and the family list already sent).
#
# Chunk size (15 families/chunk, 42 chunks) was sized conservatively off
# the <=200-taxa run's OWN measured shortfall, not a fresh calibration on
# this tier (none exists yet): that run's 36-family chunks only achieved
# 62.8% real completion in 3h (actual throughput ~23 families/chunk at
# <=200-taxa speeds, not the 36 assumed). 15/chunk applies a further
# safety margin on top of that correction, since 201-500-taxa families are
# themselves slower per-family than the <=200-taxa data this is grounded
# in, and that slowdown isn't quantified by any real timing data yet.
# Expect this run to still need a mop-up pass -- check_completion.py
# against chunks_201_500 the same way as the <=200 tranche.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
CHUNK_FILE=$(printf "%s/manifest/chunks_201_500/chunk_%04d.txt" "$DATA_DIR" "$SLURM_ARRAY_TASK_ID")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family_chunk.sh "$CHUNK_FILE" "$MANIFEST" "$OUT_DIR" config/pipeline.conf
