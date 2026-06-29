#!/bin/bash
set -euo pipefail

# Step 1: Estimate branch lengths on constrained tree for each component block
# Input: constrained.nwk, block1.aln, block2.aln, model files
# Output: block1_with_brlen.treefile, block2_with_brlen.treefile

# ============================================================================
# Setup
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -i <input_dir> -o <output_dir> -c <config_file> [-f]

Required:
  -i <dir>              Input RHG directory (e.g., inputs/RHG_001)
  -o <dir>              Output directory for step 1 results
  -c <file>             Config file

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

BLOCK1_ALN="${INPUT_DIR}/block1.aln"
BLOCK2_ALN="${INPUT_DIR}/block2.aln"
CONSTRAINED_TREE="${INPUT_DIR}/constrained.nwk"
MODEL1_FILE="${INPUT_DIR}/model_block1.txt"
MODEL2_FILE="${INPUT_DIR}/model_block2.txt"

for file in "$BLOCK1_ALN" "$BLOCK2_ALN" "$CONSTRAINED_TREE" "$MODEL1_FILE" "$MODEL2_FILE"; do
    if [[ ! -f "$file" ]]; then
        log "ERROR: Required file not found: $file"
        exit 1
    fi
done

# Read models
MODEL1=$(cat "$MODEL1_FILE" | tr -d '[:space:]')
MODEL2=$(cat "$MODEL2_FILE" | tr -d '[:space:]')

log "Model for block 1: $MODEL1"
log "Model for block 2: $MODEL2"

# ============================================================================
# Process Block 1
# ============================================================================

log "Processing Block 1..."

OUTPUT1="${OUTPUT_DIR}/block1_with_brlen"

if [[ -f "${OUTPUT1}.treefile" && "$FORCE_RERUN" != "true" ]]; then
    log "Block 1 branch length file already exists, skipping"
else
    $IQTREE_BIN \
        -s "$BLOCK1_ALN" \
        -te "$CONSTRAINED_TREE" \
        -m "$MODEL1" \
        -keep-ident \
        -redo \
        -prefix "$OUTPUT1" \
        >> "$LOG" 2>&1

    if [[ ! -f "${OUTPUT1}.treefile" ]]; then
        log "ERROR: IQTree failed for block 1"
        exit 1
    fi
fi

log "Block 1 branch lengths estimated: ${OUTPUT1}.treefile"

# ============================================================================
# Process Block 2
# ============================================================================

log "Processing Block 2..."

OUTPUT2="${OUTPUT_DIR}/block2_with_brlen"

if [[ -f "${OUTPUT2}.treefile" && "$FORCE_RERUN" != "true" ]]; then
    log "Block 2 branch length file already exists, skipping"
else
    $IQTREE_BIN \
        -s "$BLOCK2_ALN" \
        -te "$CONSTRAINED_TREE" \
        -m "$MODEL2" \
        -keep-ident \
        -redo \
        -prefix "$OUTPUT2" \
        >> "$LOG" 2>&1

    if [[ ! -f "${OUTPUT2}.treefile" ]]; then
        log "ERROR: IQTree failed for block 2"
        exit 1
    fi
fi

log "Block 2 branch lengths estimated: ${OUTPUT2}.treefile"

# ============================================================================
# Verify Output
# ============================================================================

log "Verifying output files..."

if [[ -f "${OUTPUT1}.treefile" && -f "${OUTPUT2}.treefile" ]]; then
    log "SUCCESS: Step 1 completed"
    log "Output trees ready for SeqGen in Step 2"
else
    log "ERROR: Output trees not found"
    exit 1
fi
