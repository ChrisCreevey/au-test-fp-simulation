#!/bin/bash
#SBATCH --job-name=mid_501_1000
#SBATCH --partition=k2-hipri
#SBATCH --array=1-35
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/mid_501_1000/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/mid_501_1000/slurm_%A_%a.err

# Third bulk run: the 279 families where the largest domain is 501-1000
# taxa -- the band we'd flagged as borderline (see docs/METHODS.md /
# project history) and are now taking on ourselves rather than handing to
# Peter, since queue conditions on k2-hipri are currently favourable
# (fast starts, cheap mop-up iteration) -- NOT because the per-family
# compute cost has changed. >1000 taxa stays with Peter: that tier risks
# a single replicate's search alone exceeding the 3h cap, which no amount
# of fast requeuing fixes (IQTree doesn't checkpoint mid-search).
#
# Chunk size (8 families/chunk, 35 chunks) is a conservative guess, NOT
# calibrated against real 501-1000-taxa timing data (none exists) --
# explicitly skipped calibration this round in favour of a smaller chunk
# size than the 201-500 tier's 15/chunk, scaled down for this band's
# roughly 2x higher taxon ceiling. Expect a mop-up pass may still be
# needed; check_completion.py against chunks_501_1000 the same way as the
# other tranches.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
CHUNK_FILE=$(printf "%s/manifest/chunks_501_1000/chunk_%04d.txt" "$DATA_DIR" "$SLURM_ARRAY_TASK_ID")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family_chunk.sh "$CHUNK_FILE" "$MANIFEST" "$OUT_DIR" config/pipeline.conf
