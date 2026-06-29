#!/bin/bash
set -euo pipefail

# Step 5: Run AU test on each replicate
# For each replicate, run AU test comparing constrained vs unconstrained trees
# Input: rep_*/block1_unconstrained.treefile, rep_*/block1_constrained.treefile, etc.
# Output: rep_*/au_test_block1.txt, rep_*/au_test_block2.txt (AU test results)

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

MODEL1=$(cat "${INPUT_DIR}/model_block1.txt" | tr -d '[:space:]')
MODEL2=$(cat "${INPUT_DIR}/model_block2.txt" | tr -d '[:space:]')

log "Models: Block 1: $MODEL1, Block 2: $MODEL2"

# ============================================================================
# Process Each Replicate
# ============================================================================

rep_count=0
au_reject_count1=0
au_reject_count2=0

for rep_dir in "$STEP4_OUTPUT_DIR"/rep_*; do
    rep_name=$(basename "$rep_dir")
    rep_num=$(echo "$rep_name" | grep -oP '(?<=rep_)\d+')

    if [[ -z "$rep_num" ]]; then
        log "WARN: Could not parse replicate number from $rep_name, skipping"
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

    ALN1=$(ls "$STEP3_OUTPUT_DIR"/*_block1.aln 2>/dev/null | sed -n "${rep_num}p")
    ALN2=$(ls "$STEP3_OUTPUT_DIR"/*_block2.aln 2>/dev/null | sed -n "${rep_num}p")

    UNCONSTRAINED1="${rep_dir}/block1_unconstrained.treefile"
    UNCONSTRAINED2="${rep_dir}/block2_unconstrained.treefile"
    CONSTRAINED1="${rep_dir}/block1_constrained.treefile"
    CONSTRAINED2="${rep_dir}/block2_constrained.treefile"

    if [[ ! -f "$ALN1" ]] || [[ ! -f "$UNCONSTRAINED1" ]] || [[ ! -f "$CONSTRAINED1" ]]; then
        log "WARN: Missing files for replicate $rep_num, skipping"
        continue
    fi

    # ====================================================================
    # Block 1 AU Test
    # ====================================================================

    AU_OUT1="${REP_OUTPUT_DIR}/au_test_block1"

    if [[ ! -f "${AU_OUT1}.txt" || "$FORCE_RERUN" == "true" ]]; then
        $IQTREE_BIN \
            -s "$ALN1" \
            -m "$MODEL1" \
            -z "$CONSTRAINED1" \
            -te "$UNCONSTRAINED1" \
            -au \
            -redo \
            -prefix "$AU_OUT1" \
            >> "$LOG" 2>&1
    fi

    # Extract AU test p-value for block 1
    if [[ -f "${AU_OUT1}.txt" ]]; then
        au_pval1=$(grep -oP "AU\s+\K[0-9.]+" "${AU_OUT1}.txt" | head -1 || echo "NA")

        if [[ "$au_pval1" != "NA" && $(echo "$au_pval1 < 0.05" | bc) -eq 1 ]]; then
            ((au_reject_count1++))
        fi
    fi

    # ====================================================================
    # Block 2 AU Test
    # ====================================================================

    if [[ -f "$ALN2" ]]; then
        AU_OUT2="${REP_OUTPUT_DIR}/au_test_block2"

        if [[ ! -f "${AU_OUT2}.txt" || "$FORCE_RERUN" == "true" ]]; then
            $IQTREE_BIN \
                -s "$ALN2" \
                -m "$MODEL2" \
                -z "$CONSTRAINED2" \
                -te "$UNCONSTRAINED2" \
                -au \
                -redo \
                -prefix "$AU_OUT2" \
                >> "$LOG" 2>&1
        fi

        # Extract AU test p-value for block 2
        if [[ -f "${AU_OUT2}.txt" ]]; then
            au_pval2=$(grep -oP "AU\s+\K[0-9.]+" "${AU_OUT2}.txt" | head -1 || echo "NA")

            if [[ "$au_pval2" != "NA" && $(echo "$au_pval2 < 0.05" | bc) -eq 1 ]]; then
                ((au_reject_count2++))
            fi
        fi
    fi
done

log "Processed $rep_count replicates for AU testing"
log "Block 1 AU rejections (p < 0.05): $au_reject_count1 / $rep_count"
log "Block 2 AU rejections (p < 0.05): $au_reject_count2 / $rep_count"

# ============================================================================
# Calculate and Save Summary
# ============================================================================

SUMMARY_FILE="${OUTPUT_DIR}/au_test_summary.txt"

{
    echo "AU Test Summary"
    echo "==============="
    echo "Total replicates: $rep_count"
    echo ""
    echo "Block 1:"
    echo "  AU test rejections (p < 0.05): $au_reject_count1"
    if [[ $rep_count -gt 0 ]]; then
        fp_rate1=$(echo "scale=4; 100 * $au_reject_count1 / $rep_count" | bc)
        echo "  False positive rate: ${fp_rate1}%"
    fi
    echo ""
    echo "Block 2:"
    echo "  AU test rejections (p < 0.05): $au_reject_count2"
    if [[ $rep_count -gt 0 ]]; then
        fp_rate2=$(echo "scale=4; 100 * $au_reject_count2 / $rep_count" | bc)
        echo "  False positive rate: ${fp_rate2}%"
    fi
} > "$SUMMARY_FILE"

log "Summary written to: $SUMMARY_FILE"
log "SUCCESS: Step 5 completed"
log "Results ready for summarization in Step 6"
