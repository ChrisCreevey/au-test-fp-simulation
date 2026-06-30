#!/bin/bash
set -euo pipefail

# Step 3: Apply gap pattern from original alignment to simulated sequences
# Input: simulated_combined.phy (ungapped), original combined.aln (with gaps)
# Output: rep_001_masked.aln, rep_002_masked.aln, ... (with gap pattern applied)

# ============================================================================
# Setup
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -i <input_dir> -o <output_dir> -c <config_file> [-f]

Required:
  -i <input_dir>        Input RHG directory
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

while getopts "i:o:c:fh" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) CONFIG_FILE="$OPTARG" ;;
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

log "Starting Step 3: Apply Gap Pattern from Original Alignment"
log "Output directory: $OUTPUT_DIR"

# Find Step 2 output directory
PARENT_OUTPUT_DIR=$(dirname "$OUTPUT_DIR")
STEP2_OUTPUT_DIR=$(ls -d "$PARENT_OUTPUT_DIR"/step2_* 2>/dev/null | head -1)

if [[ -z "$STEP2_OUTPUT_DIR" ]]; then
    log "ERROR: Could not find Step 2 output directory"
    exit 1
fi

log "Using Step 2 output from: $STEP2_OUTPUT_DIR"

# ============================================================================
# Check Inputs
# ============================================================================

ORIGINAL_ALN="${INPUT_DIR}/combined.aln"
SIM="${STEP2_OUTPUT_DIR}/simulated_combined.phy"

for file in "$ORIGINAL_ALN" "$SIM"; do
    if [[ ! -f "$file" ]]; then
        log "ERROR: Required file not found: $file"
        exit 1
    fi
done

log "Original alignment: $ORIGINAL_ALN"
log "Simulated sequences: $SIM"

# ============================================================================
# Call Python Masking Script
# ============================================================================

log "Running Python masking script..."

"$PYTHON_BIN" "$SRC_DIR/apply_gap_pattern.py" \
    --original-aln "$ORIGINAL_ALN" \
    --simulated-seqs "$SIM" \
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

expected_files=$N_REPLICATES
actual_files=$(ls -1 "$OUTPUT_DIR"/rep_*_masked.aln 2>/dev/null | wc -l)

log "Expected files: $expected_files"
log "Actual files: $actual_files"

if [[ $actual_files -lt $expected_files ]]; then
    log "WARN: Not all expected files were created (may be intentional if resuming)"
fi

if [[ $actual_files -gt 0 ]]; then
    log "SUCCESS: Step 3 completed"
    log "Masked alignments ready for tree reconstruction in Step 4"
else
    log "ERROR: No output files created"
    exit 1
fi
