#!/bin/bash
set -euo pipefail

# Step 4: Reconstruct gene trees for each simulated replicate
# For each of 100 replicates:
#   - Build unconstrained tree from mixed alignment
#   - Build constrained tree forcing fused sequences to be monophyletic
# Input: rep_001_mixed.aln, rep_002_mixed.aln, ... (from Step 3)
# Output: rep_*/mixed_unconstrained.treefile, rep_*/mixed_constrained.treefile

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
  -o <output_dir>       Output directory for step 4 results
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

LOG="${OUTPUT_DIR}/step4.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"
}

log "Starting Step 4: Reconstruct Gene Trees"
log "Output directory: $OUTPUT_DIR"

# Find Step 3 output directory
if [[ -z "$STEP3_OUTPUT_DIR" ]]; then
    PARENT_OUTPUT_DIR=$(dirname "$OUTPUT_DIR")
    STEP3_OUTPUT_DIR=$(ls -d "$PARENT_OUTPUT_DIR"/step3_* 2>/dev/null | head -1)

    if [[ -z "$STEP3_OUTPUT_DIR" ]]; then
        log "ERROR: Could not find Step 3 output directory"
        exit 1
    fi
fi

log "Using Step 3 output from: $STEP3_OUTPUT_DIR"

# ============================================================================
# Check Inputs
# ============================================================================

MODEL1=$(cat "${INPUT_DIR}/model_block1.txt" | tr -d '[:space:]')
MODEL2=$(cat "${INPUT_DIR}/model_block2.txt" | tr -d '[:space:]')
CONSTRAINED_TREE="${INPUT_DIR}/constrained.nwk"

if [[ ! -f "$CONSTRAINED_TREE" ]]; then
    log "ERROR: Constrained tree not found: $CONSTRAINED_TREE"
    exit 1
fi

# For mixed alignments, use a combined model or take block1's model
# (both blocks are present in the mixed alignment)
COMBINED_MODEL="${MODEL1}"
log "Using model: $COMBINED_MODEL"

# ============================================================================
# Process Each Replicate
# ============================================================================

rep_count=0
total_reps=$(ls -1 "$STEP3_OUTPUT_DIR"/rep_*_mixed.aln 2>/dev/null | wc -l)

log "Found $total_reps replicates to process"

for mixed_aln in "$STEP3_OUTPUT_DIR"/rep_*_mixed.aln; do
    if [[ ! -f "$mixed_aln" ]]; then
        continue
    fi

    rep_name=$(basename "$mixed_aln" _mixed.aln)
    ((rep_count++))

    if [[ $((rep_count % 10)) -eq 0 ]]; then
        log "Processing replicate $rep_count / $total_reps"
    fi

    REP_OUTPUT_DIR="${OUTPUT_DIR}/${rep_name}"
    mkdir -p "$REP_OUTPUT_DIR"

    # ====================================================================
    # Build Unconstrained Tree
    # ====================================================================

    UNCONSTRAINED="${REP_OUTPUT_DIR}/mixed_unconstrained"

    if [[ ! -f "${UNCONSTRAINED}.treefile" || "$FORCE_RERUN" == "true" ]]; then
        $IQTREE_BIN \
            -s "$mixed_aln" \
            -m "$COMBINED_MODEL" \
            -bb "$N_BOOTSTRAP" \
            -redo \
            -prefix "$UNCONSTRAINED" \
            >> "$LOG" 2>&1

        if [[ ! -f "${UNCONSTRAINED}.treefile" ]]; then
            log "ERROR: IQTree failed for unconstrained tree of $rep_name"
            exit 1
        fi
    fi

    # ====================================================================
    # Build Constrained Tree
    # ====================================================================

    CONSTRAINED="${REP_OUTPUT_DIR}/mixed_constrained"

    if [[ ! -f "${CONSTRAINED}.treefile" || "$FORCE_RERUN" == "true" ]]; then
        $IQTREE_BIN \
            -s "$mixed_aln" \
            -te "$CONSTRAINED_TREE" \
            -m "$COMBINED_MODEL" \
            -bb "$N_BOOTSTRAP" \
            -redo \
            -prefix "$CONSTRAINED" \
            >> "$LOG" 2>&1

        if [[ ! -f "${CONSTRAINED}.treefile" ]]; then
            log "ERROR: IQTree failed for constrained tree of $rep_name"
            exit 1
        fi
    fi
done

log "Processed $rep_count replicates"
log "SUCCESS: Step 4 completed"
log "Gene trees ready for AU testing in Step 5"
