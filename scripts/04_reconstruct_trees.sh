#!/bin/bash
set -euo pipefail

# Step 4: Reconstruct gene trees for each simulated replicate
# For each of 100 replicates:
#   - Build unconstrained tree for block 1 and block 2
#   - Build constrained tree forcing fusion sequences to be monophyletic
# Input: rep_001_block1.aln, rep_001_block2.aln, ... (from Step 3)
# Output: rep_*/block1_unconstrained.treefile, rep_*/block2_constrained.treefile, etc.

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

log "Models: Block 1: $MODEL1, Block 2: $MODEL2"

# ============================================================================
# Process Each Replicate
# ============================================================================

rep_count=0
total_reps=$(ls -1d "$STEP3_OUTPUT_DIR"/rep_* 2>/dev/null | wc -l)

log "Found $total_reps replicates to process"

for rep_dir in "$STEP3_OUTPUT_DIR"/rep_*; do
    rep_name=$(basename "$rep_dir")
    rep_num=$(echo "$rep_name" | grep -oP '(?<=rep_)\d+')

    if [[ -z "$rep_num" ]]; then
        # Handle different naming convention (might have just numbers)
        rep_num=$(ls "$rep_dir"/*_block1.aln 2>/dev/null | head -1 | grep -oP '_block1' || echo "")
        if [[ -z "$rep_num" ]]; then
            log "WARN: Could not parse replicate number from $rep_dir, skipping"
            continue
        fi
    fi

    ((rep_count++))
    if [[ $((rep_count % 10)) -eq 0 ]]; then
        log "Processing replicate $rep_count / $total_reps"
    fi

    REP_OUTPUT_DIR="${OUTPUT_DIR}/${rep_name}"
    mkdir -p "$REP_OUTPUT_DIR"

    # Locate alignment files from step 3
    ALN1=$(ls "$STEP3_OUTPUT_DIR"/*_rep_${rep_num:(-3)}_block1.aln 2>/dev/null | head -1 || \
           ls "$STEP3_OUTPUT_DIR"/rep_*${rep_num:(-3)}_block1.aln 2>/dev/null | head -1 || \
           ls "$STEP3_OUTPUT_DIR"/rep_*/block1.aln 2>/dev/null | sed -n "${rep_num}p")

    if [[ -z "$ALN1" || ! -f "$ALN1" ]]; then
        # Try alternative pattern
        ALN1=$(ls "$STEP3_OUTPUT_DIR"/rep_*_block1.aln 2>/dev/null | sed -n "${rep_num}p")
    fi

    ALN2=$(echo "$ALN1" | sed 's/_block1/_block2/')

    if [[ ! -f "$ALN1" ]] || [[ ! -f "$ALN2" ]]; then
        log "ERROR: Could not find alignment files for replicate $rep_num"
        log "  Looked for: $ALN1, $ALN2"
        exit 1
    fi

    # ========================================================================
    # Build Unconstrained Trees
    # ========================================================================

    UNCONSTRAINED1="${REP_OUTPUT_DIR}/block1_unconstrained"
    if [[ ! -f "${UNCONSTRAINED1}.treefile" || "$FORCE_RERUN" == "true" ]]; then
        $IQTREE_BIN \
            -s "$ALN1" \
            -m "$MODEL1" \
            -bb "$N_BOOTSTRAP" \
            -redo \
            -prefix "$UNCONSTRAINED1" \
            >> "$LOG" 2>&1
    fi

    UNCONSTRAINED2="${REP_OUTPUT_DIR}/block2_unconstrained"
    if [[ ! -f "${UNCONSTRAINED2}.treefile" || "$FORCE_RERUN" == "true" ]]; then
        $IQTREE_BIN \
            -s "$ALN2" \
            -m "$MODEL2" \
            -bb "$N_BOOTSTRAP" \
            -redo \
            -prefix "$UNCONSTRAINED2" \
            >> "$LOG" 2>&1
    fi

    # ========================================================================
    # Build Constrained Trees
    # ========================================================================

    CONSTRAINED1="${REP_OUTPUT_DIR}/block1_constrained"
    if [[ ! -f "${CONSTRAINED1}.treefile" || "$FORCE_RERUN" == "true" ]]; then
        $IQTREE_BIN \
            -s "$ALN1" \
            -te "$CONSTRAINED_TREE" \
            -m "$MODEL1" \
            -bb "$N_BOOTSTRAP" \
            -redo \
            -prefix "$CONSTRAINED1" \
            >> "$LOG" 2>&1
    fi

    CONSTRAINED2="${REP_OUTPUT_DIR}/block2_constrained"
    if [[ ! -f "${CONSTRAINED2}.treefile" || "$FORCE_RERUN" == "true" ]]; then
        $IQTREE_BIN \
            -s "$ALN2" \
            -te "$CONSTRAINED_TREE" \
            -m "$MODEL2" \
            -bb "$N_BOOTSTRAP" \
            -redo \
            -prefix "$CONSTRAINED2" \
            >> "$LOG" 2>&1
    fi
done

log "Processed $rep_count replicates"
log "SUCCESS: Step 4 completed"
log "Gene trees ready for AU testing in Step 5"
