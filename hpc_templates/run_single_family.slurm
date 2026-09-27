#!/bin/bash
#SBATCH --job-name=au_test_sim
#SBATCH --partition=k2-sandbox
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=30:00
#SBATCH --output=slurm_%j.log
#SBATCH --error=slurm_%j.err

# Run the full pipeline (Steps 1-6, both/all domains, family-level join) for
# ONE family. Useful for testing a new family, rerunning after a failure, or
# a real-data comparison pass (-r on real_data_delta_logl.sh, off by default
# for bulk runs -- see run_family.sh).
#
# For many families at once, see bulk_le200_array.slurm (batches many
# families per array task, sized for k2-hipri's 3h cap) and
# mopup_le200_array.slurm (resubmits whatever didn't finish).
#
# Edit FAMILY_ID below, or override at submit time:
#   sbatch --export=FAMILY_ID=F100058 hpc_templates/run_single_family.slurm

set -euo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs"
FAMILY_ID="${FAMILY_ID:-F100058}"

cd "$PROJECT_DIR"
echo "Host: $(hostname), family: ${FAMILY_ID}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

bash scripts/run_family.sh -f "$FAMILY_ID" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf
