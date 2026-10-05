#!/bin/bash
#SBATCH --job-name=oom_retry
#SBATCH --partition=k2-hipri
#SBATCH --array=1-2
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/oom_retry_hipri/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/oom_retry_hipri/slurm_%A_%a.err

# Retry for F2402 and F5610, which hit OUT_OF_MEMORY at --mem=8G despite
# not being especially large (181-405 taxa) -- not part of the Step 4
# bug-affected set, just an apparent occasional memory spike from 10
# concurrent single-threaded replicates. Bumped to 32G, kept on k2-hipri
# (3h) since there's no evidence these need more time, just more memory.
# No -x: the Step 4 fix already handles any partial state from the OOM
# correctly (missing .iqtree -> redone; complete ones skipped).

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/manifest/oom_only_retry.txt")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family.sh -f "$FAMILY" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf
