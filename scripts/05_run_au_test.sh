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
# Process Each Replicate (in parallel — replicates are independent)
# ============================================================================

JOBS="${REPLICATE_JOBS:-${HPC_THREADS:-1}}"
IQTREE_THREADS_PER_JOB="${IQTREE_THREADS_PER_JOB:-1}"
FUSED_TAXA_STR="${FUSED_TAXA[*]}"

process_one_replicate_au() {
    local rep_dir="$1"
    local rep_name
    rep_name=$(basename "$rep_dir")

    local REP_OUTPUT_DIR="${OUTPUT_DIR}/${rep_name}"
    mkdir -p "$REP_OUTPUT_DIR"
    local REP_LOG="${REP_OUTPUT_DIR}/step5_rep.log"

    # ====================================================================
    # Get alignment and tree files
    # ====================================================================

    local MASKED_ALN="${STEP3_OUTPUT_DIR}/${rep_name}_masked.aln"
    if [[ ! -f "$MASKED_ALN" ]]; then
        echo "SKIP $rep_name reason=no_masked_aln"
        return 0
    fi

    local UNCONSTRAINED="${rep_dir}/masked_unconstrained.treefile"
    local CONSTRAINED="${rep_dir}/masked_constrained.treefile"

    if [[ ! -f "$UNCONSTRAINED" ]] || [[ ! -f "$CONSTRAINED" ]]; then
        echo "SKIP $rep_name reason=missing_trees"
        return 0
    fi

    # ====================================================================
    # Run AU Test
    # ====================================================================

    local AU_OUT="${REP_OUTPUT_DIR}/au_test_result"
    local TREES_FOR_AU="${REP_OUTPUT_DIR}/trees_for_au.nwk"

    if [[ ! -f "${AU_OUT}.iqtree" || "$FORCE_RERUN" == "true" ]]; then
        # Combine both trees into one file for AU test comparison
        cat "$CONSTRAINED" "$UNCONSTRAINED" > "$TREES_FOR_AU"

        # Run AU test comparing the two topologies
        $IQTREE_BIN \
            -s "$MASKED_ALN" \
            -m "$MODEL" \
            -z "$TREES_FOR_AU" \
            -au \
            -zb "$N_BOOTSTRAP" \
            -T "$IQTREE_THREADS_PER_JOB" \
            -redo \
            --prefix "$AU_OUT" \
            > "$REP_LOG" 2>&1
    fi

    # Extract AU test p-value from IQTree output
    local au_rejected=0
    local au_pval="NA"
    if [[ -f "${AU_OUT}.iqtree" ]]; then
        # Find tree comparison table: look for "  2" followed by a negative logL value (contains "-")
        # Extract the AU p-value (second-to-last numeric value before +/- indicator)
        au_pval=$(awk '/^[ \t]+2[ \t]+-[0-9]/ {
            # Remove +/- indicators, then extract the last numeric value (AU p-value)
            for (i = NF; i >= 1; i--) {
                field = $i
                gsub(/[+-]$/, "", field)
                if (field ~ /^[0-9.]+$/) {
                    print field
                    exit
                }
            }
        }' "${AU_OUT}.iqtree" || echo "NA")

        if [[ "$au_pval" != "NA" && "$au_pval" != "" ]]; then
            # Check if p-value < 0.05 (reject constraint)
            if (( $(echo "$au_pval < 0.05" | bc -l) )); then
                au_rejected=1
            fi
        fi
    fi

    # ====================================================================
    # Check Monophyly of Fused Clade in Both Trees
    # ====================================================================

    local TAXA
    read -ra TAXA <<< "$FUSED_TAXA_STR"

    # Check unconstrained tree (true/false/NA if untestable - see check_monophyly.py)
    local unconstrained_mono_json unconstrained_mono
    unconstrained_mono_json=$("$PYTHON_BIN" "$SRC_DIR/check_monophyly.py" \
        --tree "$UNCONSTRAINED" \
        --taxa "${TAXA[@]}" 2>/dev/null)
    if grep -q '"monophyletic": true' <<< "$unconstrained_mono_json"; then
        unconstrained_mono="true"
    elif grep -q '"monophyletic": false' <<< "$unconstrained_mono_json"; then
        unconstrained_mono="false"
    else
        unconstrained_mono="NA"
    fi

    # Check constrained tree (should always be monophyletic by definition)
    local constrained_mono_json constrained_mono
    constrained_mono_json=$("$PYTHON_BIN" "$SRC_DIR/check_monophyly.py" \
        --tree "$CONSTRAINED" \
        --taxa "${TAXA[@]}" 2>/dev/null)
    if grep -q '"monophyletic": true' <<< "$constrained_mono_json"; then
        constrained_mono="true"
    elif grep -q '"monophyletic": false' <<< "$constrained_mono_json"; then
        constrained_mono="false"
    else
        constrained_mono="NA"
    fi

    local reason="none"
    if [[ "$au_rejected" == "1" ]]; then
        if [[ "$unconstrained_mono" == "false" ]]; then
            # AU test rejected AND fused clade is NOT monophyletic
            # This is a true AU test "failure" (false positive of monophyly)
            reason="no_monophyly"
        else
            # AU test rejected but fused clade still monophyletic (or untestable)
            reason="other"
        fi
    fi

    # Alignment length, taxon count, and tree length from the unconstrained IQ-TREE report
    # (used by Step 6 to summarize alignment-length/evolutionary-rate effects across genes)
    local UNCONSTRAINED_IQTREE="${rep_dir}/masked_unconstrained.iqtree"
    local ntaxa="NA" nsites="NA" treelen="NA"
    if [[ -f "$UNCONSTRAINED_IQTREE" ]]; then
        local input_data_line
        input_data_line=$(grep -E 'Input data:' "$UNCONSTRAINED_IQTREE" | head -1)
        ntaxa=$(grep -oE '[0-9]+' <<< "$input_data_line" | sed -n '1p')
        nsites=$(grep -oE '[0-9]+' <<< "$input_data_line" | sed -n '2p')
        treelen=$(grep -oE '[0-9]+\.[0-9]+' <<< "$(grep -E 'Total tree length \(sum of branch lengths\)' "$UNCONSTRAINED_IQTREE" | head -1)")
        ntaxa=${ntaxa:-NA}
        nsites=${nsites:-NA}
        treelen=${treelen:-NA}
    fi

    # Write monophyly check results
    {
        echo "Monophyly Check Results"
        echo "======================"
        echo "Fused taxa: ${TAXA[@]}"
        echo ""
        echo "Constrained tree (should be monophyletic): $constrained_mono"
        echo "Unconstrained tree: $unconstrained_mono"
        echo ""
        echo "AU test p-value: ${au_pval:-NA}"
        echo "AU test rejected constraint: $([[ $au_rejected -eq 1 ]] && echo 'YES' || echo 'NO')"
        echo ""
        if [[ "$au_rejected" == "1" ]]; then
            if [[ "$unconstrained_mono" == "false" ]]; then
                echo "AU REJECTION REASON: Fused clade lost monophyly (true signal of polyphyly)"
            else
                echo "AU REJECTION REASON: Topology changed, but fused clade still monophyletic"
            fi
        fi
    } > "${REP_OUTPUT_DIR}/monophyly_check.txt"

    echo "RESULT $rep_name rejected=$au_rejected reason=$reason mono=$unconstrained_mono ntaxa=$ntaxa nsites=$nsites treelen=$treelen pval=$au_pval"
}
export -f process_one_replicate_au
export IQTREE_BIN MODEL N_BOOTSTRAP FORCE_RERUN OUTPUT_DIR STEP3_OUTPUT_DIR PYTHON_BIN SRC_DIR IQTREE_THREADS_PER_JOB FUSED_TAXA_STR

RESULTS_FILE="${OUTPUT_DIR}/step5_results.txt"
: > "$RESULTS_FILE"

log "Running AU test for replicates (up to $JOBS concurrently)"

ls -d "$STEP4_OUTPUT_DIR"/rep_* 2>/dev/null | \
    xargs -P "$JOBS" -I{} bash -c 'process_one_replicate_au "$@"' _ {} \
    >> "$RESULTS_FILE" || true

rep_count=$(grep -c '^RESULT' "$RESULTS_FILE" || true)
au_reject_count=$(grep '^RESULT' "$RESULTS_FILE" | grep -c 'rejected=1' || true)
au_reject_no_monophyly=$(grep '^RESULT' "$RESULTS_FILE" | grep -c 'reason=no_monophyly' || true)
au_reject_other_reasons=$(grep '^RESULT' "$RESULTS_FILE" | grep -c 'reason=other' || true)
skip_count=$(grep -c '^SKIP' "$RESULTS_FILE" || true)

if [[ $skip_count -gt 0 ]]; then
    log "WARN: $skip_count replicate(s) skipped (missing alignment/trees):"
    grep '^SKIP' "$RESULTS_FILE" | while read -r line; do log "  $line"; done
fi

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
