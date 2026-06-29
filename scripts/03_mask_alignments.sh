#!/bin/bash
set -euo pipefail

# Step 3: Apply single-origin fusion pattern (masking step)
# For each replicate:
#   - Fusion clade taxa: concatenate block1 + block2 sequences
#   - Non-fusion taxa: keep block1 and block2 separate
# Input: simulated_block1.phy, simulated_block2.phy from Step 2
# Output: rep_001_block1.aln, rep_001_block2.aln, ..., rep_100_block1.aln, rep_100_block2.aln

# ============================================================================
# Setup
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -i <input_dir> -o <output_dir> -c <config_file> [-f]

Required:
  -i <input_dir>        Input RHG directory (contains fusion_taxa.txt)
  -o <output_dir>       Output directory for step 3 results
  -c <config_file>      Config file

Optional:
  -f                    Force re-run
  -h                    Show help

EOF
    exit 1
}

# ============================================================================
# Parse Arguments
# ============================================================================

INPUT_DIR=""
OUTPUT_DIR=""
CONFIG_FILE=""
FORCE_RERUN=false
STEP2_OUTPUT_DIR=""

while getopts "i:o:c:s:fh" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) CONFIG_FILE="$OPTARG" ;;
        s) STEP2_OUTPUT_DIR="$OPTARG" ;;
        f) FORCE_RERUN=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$CONFIG_FILE" ]]; then
    usage
fi

# ============================================================================
# Load Configuration
# ============================================================================

source "$CONFIG_FILE"
mkdir -p "$OUTPUT_DIR"

LOG="${OUTPUT_DIR}/step3.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"
}

log "Starting Step 3: Apply Single-Origin Fusion Masking"
log "Output directory: $OUTPUT_DIR"

# Find Step 2 output directory
if [[ -z "$STEP2_OUTPUT_DIR" ]]; then
    PARENT_OUTPUT_DIR=$(dirname "$OUTPUT_DIR")
    STEP2_OUTPUT_DIR=$(ls -d "$PARENT_OUTPUT_DIR"/step2_* 2>/dev/null | head -1)

    if [[ -z "$STEP2_OUTPUT_DIR" ]]; then
        log "ERROR: Could not find Step 2 output directory"
        exit 1
    fi
fi

log "Using Step 2 output from: $STEP2_OUTPUT_DIR"

# ============================================================================
# Check Inputs
# ============================================================================

PA_FILE="${INPUT_DIR}/presence_absence.txt"
SIM1="${STEP2_OUTPUT_DIR}/simulated_block1.phy"
SIM2="${STEP2_OUTPUT_DIR}/simulated_block2.phy"

for file in "$PA_FILE" "$SIM1" "$SIM2"; do
    if [[ ! -f "$file" ]]; then
        log "ERROR: Required file not found: $file"
        exit 1
    fi
done

log "Presence/absence file: $PA_FILE"

# ============================================================================
# Call Python Masking Script
# ============================================================================

log "Running Python masking script..."

"$PYTHON_BIN" "$SRC_DIR/apply_fusion_mask.py" \
    --sim-block1 "$SIM1" \
    --sim-block2 "$SIM2" \
    --presence-absence-file "$PA_FILE" \
    --output-dir "$OUTPUT_DIR" \
    --n-replicates "$N_REPLICATES" \
    --log "$LOG"

if [[ $? -ne 0 ]]; then
    log "ERROR: Masking script failed"
    exit 1
fi

# ============================================================================
# Verify Output
# ============================================================================

log "Verifying output..."

expected_files=$N_REPLICATES  # 1 mixed alignment per replicate
actual_files=$(ls -1 "$OUTPUT_DIR"/rep_*_mixed.aln 2>/dev/null | wc -l)

log "Expected files: $expected_files (1 mixed alignment per replicate)"
log "Actual files: $actual_files"

if [[ $actual_files -lt $expected_files ]]; then
    log "WARN: Not all expected files were created (may be intentional if resuming)"
fi

if [[ $actual_files -gt 0 ]]; then
    log "SUCCESS: Step 3 completed"
    log "Mixed alignments (one per replicate) ready for tree reconstruction in Step 4"
else
    log "ERROR: No output files created"
    exit 1
fi
