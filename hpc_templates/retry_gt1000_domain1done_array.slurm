#!/bin/bash
#SBATCH --job-name=retry_gt1000_d1
#SBATCH --partition=k2-hipri
#SBATCH --array=1-39
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/retry_gt1000_d1/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/retry_gt1000_d1/slurm_%A_%a.err

# Retry the 39 >1000-taxa families whose Domain1 already finished last
# round (Domain2 is all that's left) -- one family per task, no batching,
# matching what already worked best for this tier. The 312 families where
# a full solo 3h didn't even finish Domain1 are deliberately NOT included
# here -- real evidence a 3h cap isn't enough for them.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/manifest/families_gt1000_d1done_retry.txt")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family.sh -f "$FAMILY" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf
