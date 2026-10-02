#!/bin/bash
#SBATCH --job-name=step4_bugfix_med
#SBATCH --partition=k2-medpri
#SBATCH --array=1-23
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/step4_bug_fix_medpri/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/step4_bug_fix_medpri/slurm_%A_%a.err

# Step 4 bug remediation, >1000-taxa half (23 of the 160 affected
# families -- see step4_bug_fix_rerun_array.slurm's original comment for
# the full bug description). Only these are routed to k2-medpri with a
# generous 12h budget; the other 137 (<=1000 taxa) went to k2-hipri
# instead since that size range mostly finishes well within 3h.
#
# One family per task, force-rerun (-x) -- Step 5's own skip-check has no
# way to detect that its Step 4 input changed, so a plain rerun would
# leave stale AU-test results in place even after Step 4 redoes the tree.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/manifest/step4_bug_affected_gt1000.txt")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family.sh -f "$FAMILY" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf -x
