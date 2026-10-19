#!/bin/bash
#SBATCH --job-name=step5_bugfix_lo
#SBATCH --partition=k2-lowpri
#SBATCH --array=1-104
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=5-00:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/step5_bug_retry_lowpri/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/step5_bug_retry_lowpri/slurm_%A_%a.err

# Step 5 bug remediation (lowpri half) -- see step5_bug_retry_hipri_array.slurm
# for the full bug description. This half (104 of 250 families): either a
# sibling replicate in this family-domain genuinely took >2h (own timing
# evidence, not a taxa-count guess), or every replicate in the domain was
# affected so there's no sibling to check at all (6 families) -- too few
# to bother guessing on, defaulted here to be safe. One of these (F121365
# Domain1) only needs an AU-test-only redo (no missing Step 4 tree), so it
# finishes in seconds regardless of partition -- included here anyway for
# simplicity rather than a one-off third job.
#
# k2-lowpri, no effective time cap -- step5/ and step6/ for the affected
# domain(s) were deleted ahead of this submission. Run WITHOUT -x: Step 4's
# already-fixed skip-check auto-detects and redoes ONLY the specific
# replicates missing their .iqtree report; Step 5 and Step 6 then redo for
# real since their directories are gone.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/manifest/step5_bug_retry_lowpri.txt")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family.sh -f "$FAMILY" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf
