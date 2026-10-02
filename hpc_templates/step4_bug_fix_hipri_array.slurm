#!/bin/bash
#SBATCH --job-name=step4_bugfix_hi
#SBATCH --partition=k2-hipri
#SBATCH --array=1-137
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/step4_bug_fix_hipri/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/step4_bug_fix_hipri/slurm_%A_%a.err

# Step 4 bug remediation, <=1000-taxa half (137 of the 160 affected
# families -- see step4_bug_fix_rerun_array.slurm's original comment for
# the full bug description). This size range mostly completes well within
# 3h per our own real timing data (201-500 mean 2803s/max 8693s; 501-1000
# mean 4489s/max 9113s), so k2-hipri is the right partition here -- only
# the 23 genuinely-large (>1000 taxa) families went to k2-medpri instead.
#
# One family per task, force-rerun (-x) -- see the medpri template for why
# -x is needed (Step 5 can't detect that its Step 4 input changed).

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/manifest/step4_bug_affected_le1000.txt")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family.sh -f "$FAMILY" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf -x
