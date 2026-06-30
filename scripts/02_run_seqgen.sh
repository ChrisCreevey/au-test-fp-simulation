#!/bin/bash
set -euo pipefail

# Step 2: Simulate sequences using SeqGen
# Input: branch length tree from Step 1, model, original alignment
# Output: simulated_combined.phy (100 replicates, ungapped)

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
  -o <output_dir>       Output directory for step 2 results
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
# Load Configuration and Paths
# ============================================================================

source "$CONFIG_FILE"
mkdir -p "$OUTPUT_DIR"

LOG="${OUTPUT_DIR}/step2.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"
}

log "Starting Step 2: Run SeqGen"
log "Input directory: $INPUT_DIR"
log "Output directory: $OUTPUT_DIR"

# Find Step 1 output directory (look in parent's step1_* subdirectories)
PARENT_OUTPUT_DIR=$(dirname "$OUTPUT_DIR")
STEP1_OUTPUT_DIR=$(ls -d "$PARENT_OUTPUT_DIR"/step1_* 2>/dev/null | head -1)

if [[ -z "$STEP1_OUTPUT_DIR" ]]; then
    log "ERROR: Could not find Step 1 output directory"
    exit 1
fi

log "Using Step 1 output from: $STEP1_OUTPUT_DIR"

# ============================================================================
# Check Inputs
# ============================================================================

COMBINED_ALN="${INPUT_DIR}/combined.aln"
MODEL_FILE="${INPUT_DIR}/model.txt"
TREE="${STEP1_OUTPUT_DIR}/combined_with_brlen.treefile"

for file in "$COMBINED_ALN" "$MODEL_FILE" "$TREE"; do
    if [[ ! -f "$file" ]]; then
        log "ERROR: Required file not found: $file"
        exit 1
    fi
done

MODEL=$(cat "$MODEL_FILE" | tr -d '[:space:]')

# Get alignment length (ungapped)
ALN_LEN=$(grep -v "^>" "$COMBINED_ALN" | head -1 | wc -c)
ALN_LEN=$((ALN_LEN - 1))  # Remove newline character

log "Model: $MODEL"
log "Alignment length (ungapped): $ALN_LEN"

# ============================================================================
# Helper Function: Parse Model and Extract Parameters
# ============================================================================

parse_model_for_seqgen() {
    local iqtree_model=$1
    local tree_file=$2
    local iqtree_log=""

    # Remove the tree file extension and look for .iqtree log
    local tree_prefix="${tree_file%.treefile}"
    iqtree_log="${tree_prefix}.iqtree"

    # Build SeqGen model string
    local seqgen_model=""

    if [[ "$iqtree_model" =~ GTR ]]; then
        seqgen_model="-mGTR"
    elif [[ "$iqtree_model" =~ HKY ]]; then
        seqgen_model="-mHKY"
    elif [[ "$iqtree_model" =~ JC ]]; then
        seqgen_model="-mJC"
    elif [[ "$iqtree_model" =~ K80 ]]; then
        seqgen_model="-mK80"
    elif [[ "$iqtree_model" =~ TVM ]]; then
        seqgen_model="-mTVM"
    elif [[ "$iqtree_model" =~ TrN ]]; then
        seqgen_model="-mTrN"
    else
        log "WARN: Unknown model $iqtree_model, using GTR as fallback"
        seqgen_model="-mGTR"
    fi

    # Handle gamma rate variation (+G)
    if [[ "$iqtree_model" =~ \+G ]]; then
        # Extract number of categories from +G<n> or use default
        local gamma_cats=$(echo "$iqtree_model" | grep -oP '\+G\d+' | grep -oP '\d+' || echo "4")
        local alpha=""

        # Try to extract alpha from IQTree log file
        if [[ -f "$iqtree_log" ]]; then
            alpha=$(grep -oP 'Gamma shape alpha = \K[0-9.]+' "$iqtree_log" | head -1)
        fi

        if [[ -n "$alpha" ]]; then
            seqgen_model="${seqgen_model} -g${gamma_cats} -a${alpha}"
            log "Extracted alpha from IQTree: $alpha (categories: $gamma_cats)"
        else
            log "WARN: Could not extract gamma alpha from log file, using default"
            seqgen_model="${seqgen_model} -g${gamma_cats} -a1.0"
        fi
    fi

    echo "$seqgen_model"
}

# ============================================================================
# Calculate Consensus Sequence (Ancestral)
# ============================================================================

log "Calculating majority-rule consensus sequence from original alignment..."

CONSENSUS_FILE="${OUTPUT_DIR}/consensus.phy"

"$PYTHON_BIN" "$SRC_DIR/calculate_consensus.py" \
    --alignment "$COMBINED_ALN" \
    --output "$CONSENSUS_FILE"

if [[ ! -f "$CONSENSUS_FILE" ]]; then
    log "ERROR: Failed to calculate consensus"
    exit 1
fi

log "Consensus sequence calculated for ancestral state"

# ============================================================================
# Run SeqGen
# ============================================================================

SIM_MODEL=$(parse_model_for_seqgen "$MODEL" "$TREE")

log "SeqGen command: $SEQGEN_BIN $SIM_MODEL -l$ALN_LEN -n$N_REPLICATES -s"

OUTPUT_SIM="${OUTPUT_DIR}/simulated_combined.phy"

if [[ -f "$OUTPUT_SIM" && "$FORCE_RERUN" != "true" ]]; then
    log "Simulated sequences already exist, skipping"
else
    log "Simulating sequences from consensus ancestral sequence..."

    # SeqGen with -s flag reads sequences from input
    # We provide: consensus sequence (as Phylip), then the tree
    if (cat "$CONSENSUS_FILE" && cat "$TREE") | $SEQGEN_BIN $SIM_MODEL -n"$N_REPLICATES" -s > "$OUTPUT_SIM"; then
        log "Simulation completed: $OUTPUT_SIM"
    else
        log "ERROR: SeqGen failed"
        exit 1
    fi
fi

# ============================================================================
# Verify Output
# ============================================================================

log "Verifying SeqGen output..."

if [[ ! -f "$OUTPUT_SIM" ]]; then
    log "ERROR: Output file not created"
    exit 1
fi

# Check that file contains data
SIM_LINES=$(wc -l < "$OUTPUT_SIM")

log "Output file: $SIM_LINES lines"

if [[ $SIM_LINES -lt 2 ]]; then
    log "ERROR: Output file is empty or too small"
    exit 1
fi

log "SUCCESS: Step 2 completed"
log "Simulated sequences ready for gap application in Step 3"
