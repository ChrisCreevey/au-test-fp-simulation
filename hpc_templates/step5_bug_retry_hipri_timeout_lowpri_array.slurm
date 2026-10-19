#!/bin/bash
#SBATCH --job-name=step5_bugfix_to
#SBATCH --partition=k2-lowpri
#SBATCH --array=1-15
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=5-00:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/step5_bug_retry_hipri_timeout_lowpri/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/step5_bug_retry_hipri_timeout_lowpri/slurm_%A_%a.err

# Step 5 bug remediation -- re-route for the 15 families that TIMED OUT on
# k2-hipri's 3h cap (job 10149856; see step5_bug_retry_hipri_array.slurm).
# Confirmed directly before this submission: every replicate that already
# finished under the timed-out run has a genuine, complete .iqtree report
# (ends with IQTree's own "Total wall-clock time" line) -- the timeout only
# killed whichever replicate(s) were still mid-search, nothing already
# written was corrupted. Step 4's skip-check (.iqtree-based, fixed
# 2026-10-02) means this resubmission -- same run_family.sh, no -x --
# will skip every already-complete replicate and redo ONLY the ones still
# missing. No repeated work, no re-run of the finished reps, and the
# family's other (unaffected) domain is untouched as before.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/manifest/step5_bug_retry_hipri_timeout_lowpri.txt")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family.sh -f "$FAMILY" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf
