#!/bin/bash
#SBATCH --job-name=au_test_sim
#SBATCH --array=1-100  # Change this to number of RHGs to process
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/slurm_%A_%a.log
#SBATCH --error=logs/slurm_%A_%a.err

# SLURM submission script for AU test false positive simulation pipeline
# Submits one job per RHG (parallelizes across RHGs using job arrays)
#
# Usage:
#   # Modify --array range to match number of RHGs
#   sbatch hpc_templates/submit_pipeline.slurm

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_FILE="${PROJECT_DIR}/config/pipeline.conf"
INPUT_DIR="${PROJECT_DIR}/inputs"
SCRIPT_DIR="${PROJECT_DIR}/scripts"

# Load config
source "$CONFIG_FILE"

# Create logs directory
mkdir -p logs

# ============================================================================
# Get RHG to Process (from job array index)
# ============================================================================

# Get list of RHG directories
RHG_DIRS=($(ls -d "$INPUT_DIR"/RHG_* 2>/dev/null | sort))
TOTAL_RHGS=${#RHG_DIRS[@]}

if [[ $TOTAL_RHGS -eq 0 ]]; then
    echo "ERROR: No RHG directories found in $INPUT_DIR"
    exit 1
fi

# Get current RHG (SLURM_ARRAY_TASK_ID is 1-indexed)
if [[ $SLURM_ARRAY_TASK_ID -gt $TOTAL_RHGS ]]; then
    echo "ERROR: Array task ID $SLURM_ARRAY_TASK_ID exceeds number of RHGs ($TOTAL_RHGS)"
    exit 1
fi

CURRENT_RHG="${RHG_DIRS[$((SLURM_ARRAY_TASK_ID - 1))]}"
RHG_NAME=$(basename "$CURRENT_RHG")

echo "========================================================================"
echo "Starting pipeline for RHG: $RHG_NAME"
echo "Job ID: $SLURM_JOB_ID, Array Task: $SLURM_ARRAY_TASK_ID"
echo "CPUs: $SLURM_CPUS_PER_TASK, Memory: ${SLURM_MEM_PER_NODE:-$HPC_MEMORY_PER_JOB}G"
echo "========================================================================"

# ============================================================================
# Run Pipeline
# ============================================================================

cd "$PROJECT_DIR"

bash "$SCRIPT_DIR/master_pipeline.sh" \
    -i "$CURRENT_RHG" \
    -c "$CONFIG_FILE" \
    || {
        echo "ERROR: Pipeline failed for $RHG_NAME"
        exit 1
    }

echo "========================================================================"
echo "Completed pipeline for RHG: $RHG_NAME"
echo "========================================================================"
