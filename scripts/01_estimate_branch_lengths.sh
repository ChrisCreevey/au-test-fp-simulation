#!/bin/bash
set -euo pipefail

# Step 1: Estimate branch lengths on constrained tree for combined alignment
# Input: constrained.nwk, combined.aln, model.txt
# Output: combined_with_brlen.treefile (tree with branch lengths estimated)

# ============================================================================
# Setup
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -i <input_dir> -o <output_dir> -c <config_file> [-f]

Required:
  -i <input_dir>        Input RHG directory (e.g., inputs/RHG_001)
  -o <output_dir>       Output directory for step 1 results
  -c <config_file>      Config file

Optional:
  -f                    Force re-run even if output exists
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
    echo "Error: -i, -o, and -c arguments are required"
    usage
fi

# ============================================================================
# Load Configuration
# ============================================================================

source "$CONFIG_FILE"
mkdir -p "$OUTPUT_DIR"

LOG="${OUTPUT_DIR}/step1.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"
}

log "Starting Step 1: Estimate Branch Lengths"
log "Input directory: $INPUT_DIR"
log "Output directory: $OUTPUT_DIR"

# ============================================================================
# Check Inputs
# ============================================================================

COMBINED_ALN="${INPUT_DIR}/combined.aln"
CONSTRAINED_TREE="${INPUT_DIR}/constrained.nwk"
MODEL_FILE="${INPUT_DIR}/model.txt"

for file in "$COMBINED_ALN" "$CONSTRAINED_TREE" "$MODEL_FILE"; do
    if [[ ! -f "$file" ]]; then
        log "ERROR: Required file not found: $file"
        exit 1
    fi
done

MODEL=$(cat "$MODEL_FILE" | tr -d '[:space:]')
log "Model: $MODEL"

# ============================================================================
# Check if Tree Already Has Branch Lengths
# ============================================================================

log "Checking if constrained tree already has branch lengths..."

# Check if tree has colons followed by numbers (indicating branch lengths)
if grep -q ":[0-9]" "$CONSTRAINED_TREE"; then
    log "✓ Tree already has branch lengths"

    if [[ "$FORCE_RERUN" != "true" ]]; then
        log "Skipping Step 1 (using existing branch lengths from input tree)"
        log "Copying tree to output for Step 2"
        cp "$CONSTRAINED_TREE" "${OUTPUT_DIR}/combined_with_brlen.treefile"
        log "SUCCESS: Step 1 skipped (branch lengths already present)"
        exit 0
    else
        log "Force re-run flag set - will recalculate branch lengths anyway"
    fi
else
    log "⚠ Tree does NOT have branch lengths - will estimate them"
fi

# ============================================================================
# Estimate Branch Lengths
# ============================================================================

log "Re-estimating branch lengths on constrained tree..."

OUTPUT="${OUTPUT_DIR}/combined_with_brlen"

if [[ -f "${OUTPUT}.treefile" && "$FORCE_RERUN" != "true" ]]; then
    log "Branch length file already exists, skipping"
else
    $IQTREE_BIN \
        -s "$COMBINED_ALN" \
        -te "$CONSTRAINED_TREE" \
        -m "$MODEL" \
        -T "${IQTREE_THREADS_PER_JOB:-1}" \
        -keep-ident \
        -redo \
        -prefix "$OUTPUT" \
        >> "$LOG" 2>&1

    if [[ ! -f "${OUTPUT}.treefile" ]]; then
        log "ERROR: IQTree failed"
        exit 1
    fi
fi

log "Branch lengths estimated: ${OUTPUT}.treefile"

# ============================================================================
# Verify Output
# ============================================================================

log "Verifying output..."

if [[ -f "${OUTPUT}.treefile" ]]; then
    log "SUCCESS: Step 1 completed"
    log "Branch-length tree ready for SeqGen in Step 2"
else
    log "ERROR: Output tree file not found"
    exit 1
fi
