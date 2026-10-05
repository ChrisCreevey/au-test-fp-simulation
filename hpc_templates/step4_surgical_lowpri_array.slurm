#!/bin/bash
#SBATCH --job-name=step4_surgical
#SBATCH --partition=k2-lowpri
#SBATCH --array=1-72
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=5-00:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/step4_surgical_lowpri/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/step4_surgical_lowpri/slurm_%A_%a.err

# Step 4 bug remediation, round 2 -- surgical approach. The previous round
# (-x, full-family force-rerun) proved far more expensive than expected
# (many timeouts even at 3h/12h) because it redid BOTH domains from
# scratch, when only one domain is actually affected in every one of
# these 72 families (confirmed: zero families have both domains
# affected). This round instead: the affected domain's step5/, step6/,
# and the family's step6_family/ were deleted ahead of this submission
# (not here -- already done), leaving Step 1-3 and the OTHER domain
# completely untouched. Run WITHOUT -x: the fixed Step 4 skip-check will
# automatically redo only the specific replicates missing their .iqtree
# report and skip the rest, and the unaffected domain is never touched at
# all since its outputs still exist.
#
# k2-lowpri (no time cap) rather than another hipri/medpri time-limit
# guess -- round 1 showed even "safe" <=1000-taxa cases can need >12h
# with complex (high-order FreeRate) models, independent of taxa count.
# --mem=32G since 2 of these 72 (F679, F1243) previously hit OOM at 8G.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/manifest/step4_surgical_retry.txt")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family.sh -f "$FAMILY" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf
