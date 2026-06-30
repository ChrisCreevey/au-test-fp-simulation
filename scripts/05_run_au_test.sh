#!/bin/bash
set -euo pipefail

# Step 5: Run AU test on each replicate
# For each replicate, run AU test comparing constrained vs unconstrained trees
# Input: rep_*/masked_unconstrained.treefile, rep_*/masked_constrained.treefile
# Output: rep_*/au_test_result.txt (AU test results)

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

LOG="${OUTPUT_DIR}/step5.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"
}

log "Starting Step 5: Run AU Test"
log "Output directory: $OUTPUT_DIR"

# Find Step 3 and Step 4 output directories
PARENT_OUTPUT_DIR=$(dirname "$OUTPUT_DIR")
STEP3_OUTPUT_DIR=$(ls -d "$PARENT_OUTPUT_DIR"/step3_* 2>/dev/null | head -1)
STEP4_OUTPUT_DIR=$(ls -d "$PARENT_OUTPUT_DIR"/step4_* 2>/dev/null | head -1)

if [[ -z "$STEP3_OUTPUT_DIR" || -z "$STEP4_OUTPUT_DIR" ]]; then
    log "ERROR: Could not find Step 3/4 output directories"
    exit 1
fi

log "Using Step 3 output from: $STEP3_OUTPUT_DIR"
log "Using Step 4 output from: $STEP4_OUTPUT_DIR"

# ============================================================================
# Check Inputs
# ============================================================================

MODEL=$(cat "${INPUT_DIR}/model.txt" | tr -d '[:space:]')
log "Using model: $MODEL"

# ============================================================================
# Read Presence/Absence to Get Fused Taxa
# ============================================================================

PA_FILE="${INPUT_DIR}/presence_absence.txt"

if [[ ! -f "$PA_FILE" ]]; then
    log "ERROR: presence_absence.txt not found"
    exit 1
fi

# Extract fused taxa (skip header)
FUSED_TAXA=()
while IFS=$'\t' read -r species gene_type; do
    if [[ "$gene_type" == "fused" ]]; then
        FUSED_TAXA+=("$species")
    fi
done < <(tail -n +2 "$PA_FILE")

log "Fused taxa (count: ${#FUSED_TAXA[@]}): ${FUSED_TAXA[@]}"

if [[ ${#FUSED_TAXA[@]} -eq 0 ]]; then
    log "ERROR: No fused taxa found in presence_absence.txt"
    exit 1
fi

# ============================================================================
# Process Each Replicate
# ============================================================================

rep_count=0
au_reject_count=0
au_reject_no_monophyly=0
au_reject_other_reasons=0

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

    MASKED_ALN=$(ls "$STEP3_OUTPUT_DIR"/${rep_name}_masked.aln 2>/dev/null || echo "")

    if [[ ! -f "$MASKED_ALN" ]]; then
        log "WARN: Masked alignment not found for $rep_name, skipping"
        continue
    fi

    UNCONSTRAINED="${rep_dir}/masked_unconstrained.treefile"
    CONSTRAINED="${rep_dir}/masked_constrained.treefile"

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
            -s "$MASKED_ALN" \
            -m "$MODEL" \
            -z "$CONSTRAINED" \
            -te "$UNCONSTRAINED" \
            -au \
            -redo \
            -prefix "$AU_OUT" \
            >> "$LOG" 2>&1
    fi

    # Extract AU test p-value
    au_rejected=0
    if [[ -f "${AU_OUT}.txt" ]]; then
        au_pval=$(grep -oP "AU\s+\K[0-9.]+" "${AU_OUT}.txt" | head -1 || echo "NA")

        if [[ "$au_pval" != "NA" ]]; then
            # Check if p-value < 0.05 (reject constraint)
            if (( $(echo "$au_pval < 0.05" | bc -l) )); then
                ((au_reject_count++))
                au_rejected=1
            fi
        fi
    fi

    # ====================================================================
    # Check Monophyly of Fused Clade in Both Trees
    # ====================================================================

    # Check unconstrained tree
    unconstrained_mono=$("$PYTHON_BIN" "$SRC_DIR/check_monophyly.py" \
        --tree "$UNCONSTRAINED" \
        --taxa "${FUSED_TAXA[@]}" 2>/dev/null | grep -o '"monophyletic": true' || echo "false")

    # Check constrained tree (should always be monophyletic by definition)
    constrained_mono=$("$PYTHON_BIN" "$SRC_DIR/check_monophyly.py" \
        --tree "$CONSTRAINED" \
        --taxa "${FUSED_TAXA[@]}" 2>/dev/null | grep -o '"monophyletic": true' || echo "false")

    if [[ "$au_rejected" == "1" ]]; then
        if [[ "$unconstrained_mono" != "false" ]]; then
            # AU test rejected but fused clade still monophyletic
            ((au_reject_other_reasons++))
        else
            # AU test rejected AND fused clade is NOT monophyletic
            # This is a true AU test "failure" (false positive of monophyly)
            ((au_reject_no_monophyly++))
        fi
    fi

    # Write monophyly check results
    {
        echo "Monophyly Check Results"
        echo "======================"
        echo "Fused taxa: ${FUSED_TAXA[@]}"
        echo ""
        echo "Constrained tree (should be monophyletic): $constrained_mono"
        echo "Unconstrained tree: $unconstrained_mono"
        echo ""
        echo "AU test p-value: ${au_pval:-NA}"
        echo "AU test rejected constraint: $([[ $au_rejected -eq 1 ]] && echo 'YES' || echo 'NO')"
        echo ""
        if [[ "$au_rejected" == "1" ]]; then
            if [[ "$unconstrained_mono" != "false" ]]; then
                echo "AU REJECTION REASON: Topology changed, but fused clade still monophyletic"
            else
                echo "AU REJECTION REASON: Fused clade lost monophyly (true signal of polyphyly)"
            fi
        fi
    } > "${REP_OUTPUT_DIR}/monophyly_check.txt"
done

log "Processed $rep_count replicates for AU testing"
log "AU test rejections (p < 0.05): $au_reject_count / $rep_count"
log "  - Rejections due to loss of monophyly: $au_reject_no_monophyly"
log "  - Rejections for other reasons: $au_reject_other_reasons"

# ============================================================================
# Calculate and Save Summary
# ============================================================================

SUMMARY_FILE="${OUTPUT_DIR}/au_test_summary.txt"

{
    echo "AU Test Summary with Monophyly Checking"
    echo "========================================"
    echo "Total replicates: $rep_count"
    echo ""
    echo "AU Test Results:"
    echo "  Total rejections (p < 0.05): $au_reject_count"
    echo "  Rejections with loss of monophyly: $au_reject_no_monophyly"
    echo "  Rejections for other reasons: $au_reject_other_reasons"

    if [[ $rep_count -gt 0 ]]; then
        fp_rate=$(echo "scale=4; 100 * $au_reject_count / $rep_count" | bc)
        echo ""
        echo "False Positive Rates:"
        echo "  Total FP rate (all rejections): ${fp_rate}%"

        if [[ $au_reject_no_monophyly -gt 0 ]]; then
            mono_fp_rate=$(echo "scale=4; 100 * $au_reject_no_monophyly / $rep_count" | bc)
            echo "  Monophyly-loss FP rate (true signal): ${mono_fp_rate}%"
        fi

        if [[ $au_reject_other_reasons -gt 0 ]]; then
            other_rate=$(echo "scale=4; 100 * $au_reject_other_reasons / $rep_count" | bc)
            echo "  Other topology changes: ${other_rate}%"
        fi
    fi

    echo ""
    echo "Interpretation:"
    echo "  - Expected FP rate for monophyly-loss ~5% if null model is correct"
    echo "  - Other topology changes may occur independently of monophyly"
    echo "  - Key metric: monophyly-loss FP rate, not total rejections"
} > "$SUMMARY_FILE"

log "Summary written to: $SUMMARY_FILE"
log "SUCCESS: Step 5 completed"
log "Results ready for summarization in Step 6"
