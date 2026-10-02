#!/bin/bash
#SBATCH --job-name=mopup5_medpri
#SBATCH --partition=k2-medpri
#SBATCH --array=1-36
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/mopup5_medpri/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/mopup5_medpri/slurm_%A_%a.err

# Fifth mop-up: 36 families that got through Steps 1-3 and started Step 4
# (tree reconstruction) in a prior k2-hipri attempt, but didn't finish
# within the 3h cap. All are 339-1000 taxa on Domain1 (Domain2 is small,
# 8-221 taxa) -- not an extreme size, but several have high-order FreeRate
# models (LG+R10, WAG+F+R9, LG+R8 etc.) that cost much more per search
# iteration regardless of taxon count.
#
# IQTree's own in-progress time estimates (pulled from the killed attempts'
# logs) put most of these at roughly 2-4.5h total for the dominant Step 4
# cost -- comfortably within k2-medpri's 24h limit. --time=12h here is a
# generous margin above that, not a tight estimate -- idle allocated time
# costs nothing, under-provisioning does.
#
# One family per task, same resource profile as the hipri attempts
# (REPLICATE_JOBS=10, single-threaded replicates) -- this is a partition
# change, not a parallelism strategy change.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/manifest/mopup5_medpri.txt")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family.sh -f "$FAMILY" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf
