#!/bin/bash
set -euo pipefail

# Step 2: Simulate sequences using SeqGen
# Input: branch length trees from Step 1, model files, alignment lengths
# Output: simulated_block1.phy, simulated_block2.phy (100 replicates each)

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
STEP1_OUTPUT_DIR=""

while getopts "i:o:c:s:fh" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) CONFIG_FILE="$OPTARG" ;;
        s) STEP1_OUTPUT_DIR="$OPTARG" ;;
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
if [[ -z "$STEP1_OUTPUT_DIR" ]]; then
    PARENT_OUTPUT_DIR=$(dirname "$OUTPUT_DIR")
    STEP1_OUTPUT_DIR=$(ls -d "$PARENT_OUTPUT_DIR"/step1_* 2>/dev/null | head -1)

    if [[ -z "$STEP1_OUTPUT_DIR" ]]; then
        log "ERROR: Could not find Step 1 output directory"
        exit 1
    fi
fi

log "Using Step 1 output from: $STEP1_OUTPUT_DIR"

# ============================================================================
# Check Inputs
# ============================================================================

BLOCK1_ALN="${INPUT_DIR}/block1.aln"
BLOCK2_ALN="${INPUT_DIR}/block2.aln"
MODEL1_FILE="${INPUT_DIR}/model_block1.txt"
MODEL2_FILE="${INPUT_DIR}/model_block2.txt"
TREE1="${STEP1_OUTPUT_DIR}/block1_with_brlen.treefile"
TREE2="${STEP1_OUTPUT_DIR}/block2_with_brlen.treefile"

for file in "$BLOCK1_ALN" "$BLOCK2_ALN" "$MODEL1_FILE" "$MODEL2_FILE" "$TREE1" "$TREE2"; do
    if [[ ! -f "$file" ]]; then
        log "ERROR: Required file not found: $file"
        exit 1
    fi
done

MODEL1=$(cat "$MODEL1_FILE" | tr -d '[:space:]')
MODEL2=$(cat "$MODEL2_FILE" | tr -d '[:space:]')

# Get alignment lengths
ALN_LEN1=$(grep -v "^>" "$BLOCK1_ALN" | head -1 | wc -c)
ALN_LEN1=$((ALN_LEN1 - 1))  # Remove newline character

ALN_LEN2=$(grep -v "^>" "$BLOCK2_ALN" | head -1 | wc -c)
ALN_LEN2=$((ALN_LEN2 - 1))

log "Model for block 1: $MODEL1 (length: $ALN_LEN1)"
log "Model for block 2: $MODEL2 (length: $ALN_LEN2)"

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
    # Handle basic model (GTR, HKY, JC, etc.)
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
# Run SeqGen
# ============================================================================

SIM_MODEL1=$(parse_model_for_seqgen "$MODEL1" "$TREE1")
SIM_MODEL2=$(parse_model_for_seqgen "$MODEL2" "$TREE2")

log "SeqGen command for block 1: $SEQGEN_BIN $SIM_MODEL1 -l$ALN_LEN1 -n$N_REPLICATES"
log "SeqGen command for block 2: $SEQGEN_BIN $SIM_MODEL2 -l$ALN_LEN2 -n$N_REPLICATES"

# Block 1
OUTPUT_SIM1="${OUTPUT_DIR}/simulated_block1.phy"

if [[ -f "$OUTPUT_SIM1" && "$FORCE_RERUN" != "true" ]]; then
    log "Block 1 simulated sequences already exist, skipping"
else
    log "Simulating Block 1 sequences..."

    if cat "$TREE1" | $SEQGEN_BIN $SIM_MODEL1 -l"$ALN_LEN1" -n"$N_REPLICATES" > "$OUTPUT_SIM1"; then
        log "Block 1 simulation completed: $OUTPUT_SIM1"
    else
        log "ERROR: SeqGen failed for block 1"
        exit 1
    fi
fi

# Block 2
OUTPUT_SIM2="${OUTPUT_DIR}/simulated_block2.phy"

if [[ -f "$OUTPUT_SIM2" && "$FORCE_RERUN" != "true" ]]; then
    log "Block 2 simulated sequences already exist, skipping"
else
    log "Simulating Block 2 sequences..."

    if cat "$TREE2" | $SEQGEN_BIN $SIM_MODEL2 -l"$ALN_LEN2" -n"$N_REPLICATES" > "$OUTPUT_SIM2"; then
        log "Block 2 simulation completed: $OUTPUT_SIM2"
    else
        log "ERROR: SeqGen failed for block 2"
        exit 1
    fi
fi

# ============================================================================
# Verify Output
# ============================================================================

log "Verifying SeqGen output..."

if [[ ! -f "$OUTPUT_SIM1" ]] || [[ ! -f "$OUTPUT_SIM2" ]]; then
    log "ERROR: Output files not created"
    exit 1
fi

# Check that files contain data
SIM1_LINES=$(wc -l < "$OUTPUT_SIM1")
SIM2_LINES=$(wc -l < "$OUTPUT_SIM2")

log "Block 1 output: $SIM1_LINES lines"
log "Block 2 output: $SIM2_LINES lines"

if [[ $SIM1_LINES -lt 2 ]] || [[ $SIM2_LINES -lt 2 ]]; then
    log "ERROR: Output files are empty or too small"
    exit 1
fi

log "SUCCESS: Step 2 completed"
log "Simulated sequences ready for masking in Step 3"
