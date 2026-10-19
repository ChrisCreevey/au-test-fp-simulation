#!/bin/bash
#SBATCH --job-name=step5_bugfix_hi
#SBATCH --partition=k2-hipri
#SBATCH --array=1-146
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/step5_bug_retry_hipri/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/step5_bug_retry_hipri/slurm_%A_%a.err

# Step 5 bug remediation (hipri half): 05_run_au_test.sh used to check only
# masked_unconstrained.treefile existence, not masked_unconstrained.iqtree
# (the genuine-completion marker) -- the same bug class fixed in Step 4's
# skip-check on 2026-10-02 (9cdb1ee), but never patched here too. Result:
# Step 5 silently ran the AU test against truncated/killed mid-search trees
# for 1,117 replicates across 250 family-domains, left over from before
# that fix. Now fixed in 05_run_au_test.sh itself (this submission already
# includes that fix).
#
# This half (146 of 250 families): every family-domain here has at least
# one OWN successfully-completed sibling replicate that finished in <=2h
# (own real timing evidence, not a taxa-count guess -- taxa count alone
# doesn't reliably predict duration here, confirmed again on this batch).
# The other 104 (genuine >2h sibling evidence, or no sibling to check at
# all) went to the lowpri array instead.
#
# step5/ and step6/ for the affected domain(s) were deleted ahead of this
# submission (not here). Run WITHOUT -x: Step 4's already-fixed skip-check
# auto-detects and redoes ONLY the specific replicates missing their
# .iqtree report; Step 5 and Step 6 then redo for real since their
# directories are gone. Unaffected domains/replicates are never touched.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/manifest/step5_bug_retry_hipri.txt")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family.sh -f "$FAMILY" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf
