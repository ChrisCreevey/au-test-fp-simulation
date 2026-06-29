#!/bin/bash
set -euo pipefail

# Step 5: Run AU test on each replicate
# For each replicate, run AU test comparing constrained vs unconstrained trees
# Input: rep_*/mixed_unconstrained.treefile, rep_*/mixed_constrained.treefile
# Output: rep_*/au_test_result.txt (AU test results and p-values)

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
  -o <output_dir>       Output directory for step 5 results
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
STEP3_OUTPUT_DIR=""
STEP4_OUTPUT_DIR=""

while getopts "i:o:c:s:fh" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) CONFIG_FILE="$OPTARG" ;;
        s) STEP3_OUTPUT_DIR="$OPTARG" ;;
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

LOG="${OUTPUT_DIR}/step5.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"
}

log "Starting Step 5: Run AU Test"
log "Output directory: $OUTPUT_DIR"

# Find Step 3 and Step 4 output directories
if [[ -z "$STEP3_OUTPUT_DIR" ]]; then
    PARENT_OUTPUT_DIR=$(dirname "$OUTPUT_DIR")
    STEP3_OUTPUT_DIR=$(ls -d "$PARENT_OUTPUT_DIR"/step3_* 2>/dev/null | head -1)
fi

if [[ -z "$STEP4_OUTPUT_DIR" ]]; then
    PARENT_OUTPUT_DIR=$(dirname "$OUTPUT_DIR")
    STEP4_OUTPUT_DIR=$(ls -d "$PARENT_OUTPUT_DIR"/step4_* 2>/dev/null | head -1)
fi

if [[ -z "$STEP3_OUTPUT_DIR" || -z "$STEP4_OUTPUT_DIR" ]]; then
    log "ERROR: Could not find Step 3/4 output directories"
    exit 1
fi

log "Using Step 3 output from: $STEP3_OUTPUT_DIR"
log "Using Step 4 output from: $STEP4_OUTPUT_DIR"

# ============================================================================
# Check Inputs
# ============================================================================

MODEL=$(cat "${INPUT_DIR}/model_block1.txt" | tr -d '[:space:]')
log "Using model: $MODEL"

# ============================================================================
# Process Each Replicate
# ============================================================================

rep_count=0
au_reject_count=0

for rep_dir in "$STEP4_OUTPUT_DIR"/rep_*; do
    rep_name=$(basename "$rep_dir")

    if [[ ! -d "$rep_dir" ]]; then
        continue
    fi

    ((rep_count++))
    if [[ $((rep_count % 10)) -eq 0 ]]; then
        log "Processing AU test for replicate $rep_count"
    fi

    REP_OUTPUT_DIR="${OUTPUT_DIR}/${rep_name}"
    mkdir -p "$REP_OUTPUT_DIR"

    # ====================================================================
    # Get alignment and tree files
    # ====================================================================

    MIXED_ALN=$(ls "$STEP3_OUTPUT_DIR"/${rep_name}_mixed.aln 2>/dev/null || echo "")

    if [[ ! -f "$MIXED_ALN" ]]; then
        log "WARN: Mixed alignment not found for $rep_name, skipping"
        continue
    fi

    UNCONSTRAINED="${rep_dir}/mixed_unconstrained.treefile"
    CONSTRAINED="${rep_dir}/mixed_constrained.treefile"

    if [[ ! -f "$UNCONSTRAINED" ]] || [[ ! -f "$CONSTRAINED" ]]; then
        log "WARN: Missing tree files for $rep_name, skipping"
        continue
    fi

    # ====================================================================
    # Run AU Test
    # ====================================================================

    AU_OUT="${REP_OUTPUT_DIR}/au_test_result"

    if [[ ! -f "${AU_OUT}.txt" || "$FORCE_RERUN" == "true" ]]; then
        $IQTREE_BIN \
            -s "$MIXED_ALN" \
            -m "$MODEL" \
            -z "$CONSTRAINED" \
            -te "$UNCONSTRAINED" \
            -au \
            -redo \
            -prefix "$AU_OUT" \
            >> "$LOG" 2>&1
    fi

    # Extract AU test p-value
    if [[ -f "${AU_OUT}.txt" ]]; then
        au_pval=$(grep -oP "AU\s+\K[0-9.]+" "${AU_OUT}.txt" | head -1 || echo "NA")

        if [[ "$au_pval" != "NA" ]]; then
            # Check if p-value < 0.05 (reject monophyly)
            if (( $(echo "$au_pval < 0.05" | bc -l) )); then
                ((au_reject_count++))
            fi
        fi
    fi
done

log "Processed $rep_count replicates for AU testing"
log "AU test rejections (p < 0.05): $au_reject_count / $rep_count"

# ============================================================================
# Calculate and Save Summary
# ============================================================================

SUMMARY_FILE="${OUTPUT_DIR}/au_test_summary.txt"

{
    echo "AU Test Summary (Mixed Alignment)"
    echo "=================================="
    echo "Total replicates: $rep_count"
    echo "AU test rejections (p < 0.05): $au_reject_count"

    if [[ $rep_count -gt 0 ]]; then
        fp_rate=$(echo "scale=4; 100 * $au_reject_count / $rep_count" | bc)
        echo "False positive rate: ${fp_rate}%"
    fi

    echo ""
    echo "Interpretation:"
    echo "  If FP rate ~5%: AU test is correctly calibrated (alpha=0.05)"
    echo "  If FP rate << 5%: AU test is conservative"
    echo "  If FP rate >> 5%: AU test may be liberal"
} > "$SUMMARY_FILE"

log "Summary written to: $SUMMARY_FILE"
log "SUCCESS: Step 5 completed"
log "Results ready for summarization in Step 6"
