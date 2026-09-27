#!/bin/bash
#SBATCH --job-name=bulk_le200
#SBATCH --partition=k2-hipri
#SBATCH --array=1-113
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/bulk_le200/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/bulk_le200/slurm_%A_%a.err

# Bulk run: every family with both domains <=200 taxa (4,062 families,
# 113 chunks of ~36), ONE array submission on k2-hipri.
#
# Chunk size (36 families) was deliberately sized against a 3.5h NOMINAL
# budget even though the real --time cap here is 3h -- the calibration
# sample showed most families finish in well under a minute to a few
# minutes (avg ~348s), so packing chunks tightly against the real 3h limit
# would waste most of it; sizing against an optimistic 3.5h assumption
# means typical chunks finish comfortably while still using most of the
# real budget, and a minority of unlucky chunks (more large/slow families
# than typical) may not fully finish -- that's expected and fine, since
# every step already skips completed work: scripts/check_completion.py
# finds whichever families didn't finish across ALL chunks and compiles
# them into a mop-up chunk set to resubmit, rather than over-provisioning
# every chunk "just in case".

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
CHUNK_FILE=$(printf "%s/manifest/chunks_le200/chunk_%04d.txt" "$DATA_DIR" "$SLURM_ARRAY_TASK_ID")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family_chunk.sh "$CHUNK_FILE" "$MANIFEST" "$OUT_DIR" config/pipeline.conf
