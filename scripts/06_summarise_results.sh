#!/bin/bash
set -euo pipefail

# Step 6: Summarize results across all replicates
# Aggregate false positive rates per RHG and generate null distribution
# Input: outputs from Steps 1-5
# Output: results.tsv (per RHG), SUMMARY.txt (null distribution across families)

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
  -o <output_dir>       Output directory for step 6 results
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
STEP5_OUTPUT_DIR=""

while getopts "i:o:c:s:fh" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) CONFIG_FILE="$OPTARG" ;;
        s) STEP5_OUTPUT_DIR="$OPTARG" ;;
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

LOG="${OUTPUT_DIR}/step6.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"
}

log "Starting Step 6: Summarize Results"
log "Output directory: $OUTPUT_DIR"

# Find Step 5 output directory
if [[ -z "$STEP5_OUTPUT_DIR" ]]; then
    PARENT_OUTPUT_DIR=$(dirname "$OUTPUT_DIR")
    STEP5_OUTPUT_DIR=$(ls -d "$PARENT_OUTPUT_DIR"/step5_* 2>/dev/null | head -1)

    if [[ -z "$STEP5_OUTPUT_DIR" ]]; then
        log "ERROR: Could not find Step 5 output directory"
        exit 1
    fi
fi

log "Using Step 5 output from: $STEP5_OUTPUT_DIR"

# ============================================================================
# Call Python Summarization Script
# ============================================================================

log "Running Python summarization script..."

RHG_NAME=$(basename "$INPUT_DIR")

"$PYTHON_BIN" "$SRC_DIR/summarise_results.py" \
    --step5-output "$STEP5_OUTPUT_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --rhg-id "$RHG_NAME" \
    --log "$LOG"

if [[ $? -ne 0 ]]; then
    log "ERROR: Summarization script failed"
    exit 1
fi

# ============================================================================
# Verify and Report Results
# ============================================================================

RESULTS_FILE="${OUTPUT_DIR}/results.tsv"

if [[ -f "$RESULTS_FILE" ]]; then
    log "Results file created: $RESULTS_FILE"
    log ""
    log "===== FINAL RESULTS FOR $RHG_NAME ====="

    # Display summary
    cat "$RESULTS_FILE" | tee -a "$LOG"

    log ""
    log "SUCCESS: Step 6 completed"
    log "Per-RHG results saved to: $RESULTS_FILE"
else
    log "ERROR: Results file was not created"
    exit 1
fi
