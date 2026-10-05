#!/bin/bash
#SBATCH --job-name=slow_lowpri
#SBATCH --partition=k2-lowpri
#SBATCH --array=1-4
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=5-00:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/slow_lowpri/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/slow_lowpri/slurm_%A_%a.err

# F116, F43163, F57821, F121365 -- genuinely slow Domain1 searches (not
# Step4-bug related), driven by high-order FreeRate models (e.g.
# F121365's JTT+R10) more than taxa count alone. Didn't finish at 3h
# (hipri) or 12h (medpri) -- rather than guess a third time limit,
# k2-lowpri (no time cap, same node pool, just lower scheduling priority)
# lets these finish in one shot.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/manifest/slow_lowpri_retry.txt")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family.sh -f "$FAMILY" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf
