#!/bin/bash
#SBATCH --job-name=mopup_le200
#SBATCH --partition=k2-hipri
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/mopup_le200/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/mopup_le200/slurm_%A_%a.err

# Mop-up pass for whichever families the first bulk_le200_array.slurm pass
# didn't finish. --array is set at submit time from check_completion.py's
# reported chunk count, e.g.:
#   sbatch --array=1-<n_chunks> hpc_templates/mopup_le200_array.slurm
#
# Chunks here are smaller (default 12 families, see check_completion.py
# --mopup-chunk-size) than the first-pass 36, since this set is enriched for
# the larger/slower families that didn't make it through in time -- and
# every step still skips work already completed on the first attempt, so
# this resumes rather than redoing anything.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
CHUNK_FILE=$(printf "%s/manifest/chunks_le200_mopup/chunk_%04d.txt" "$DATA_DIR" "$SLURM_ARRAY_TASK_ID")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family_chunk.sh "$CHUNK_FILE" "$MANIFEST" "$OUT_DIR" config/pipeline.conf
