#!/bin/bash
#SBATCH --job-name=mopup3_per_fam
#SBATCH --partition=k2-hipri
#SBATCH --array=1-349
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/mopup3_per_family/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/mopup3_per_family/slurm_%A_%a.err

# Third mop-up pass: the 349 families still not done from our own
# responsibility (everything in the manifest minus the 334 given to
# Peter) -- the leftovers from the 201-500/501-1000/>1000-Domain1-done
# retry batches, plus the 15 201-500-taxa families originally excluded as
# "killed with >1h remaining" (that evidence was chunk-based, not a solo
# attempt, so worth one real per-family try here rather than handing to
# Peter -- that tier was never his to begin with).
#
# One family per task, no batching -- per-family dispatch has
# consistently outperformed chunking across every round so far.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/manifest/mopup3_remaining.txt")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family.sh -f "$FAMILY" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf
