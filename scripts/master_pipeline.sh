#!/bin/bash
set -euo pipefail

# Master pipeline orchestrator for AU test false positive simulation
# Runs all 6 steps sequentially for one RHG or in parallel across multiple RHGs

# ============================================================================
# Setup
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -i <input_rhg_dir> [-c <config_file>] [-f]

Required:
  -i <dir>              Input directory for one RHG (e.g., inputs/RHG_001)

Optional:
  -c <file>             Config file (default: config/pipeline.conf)
  -f                    Force re-run (recalculate all steps, including branch lengths)
  -h                    Show this help message

Notes:
  - If constrained tree already has branch lengths, Step 1 will be skipped
  - Use -f to force re-estimation of branch lengths (overrides auto-detection)
  - This is useful if you want to recalibrate based on the alignment

Example:
  bash $0 -i inputs/RHG_001 -c config/pipeline.conf
  bash $0 -i inputs/RHG_001 -f  (recalculate branch lengths + redo all steps)

EOF
    exit 1
}

# ============================================================================
# Parse Arguments
# ============================================================================

CONFIG_FILE="${PROJECT_DIR}/config/pipeline.conf"
RHG_INPUT_DIR=""
FORCE_RERUN=false

while getopts "i:c:fh" opt; do
    case $opt in
        i) RHG_INPUT_DIR="$(cd "$OPTARG" 2>/dev/null && pwd)" || { echo "Error: input directory not found: $OPTARG"; exit 1; } ;;
        c) CONFIG_FILE="$OPTARG" ;;
        f) FORCE_RERUN=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$RHG_INPUT_DIR" ]]; then
    echo "Error: -i argument is required"
    usage
fi

# ============================================================================
# Load Configuration
# ============================================================================

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: config file not found: $CONFIG_FILE"
    exit 1
fi

source "$CONFIG_FILE"

# Set OUTPUT_DIR from config (DEFAULT_OUTPUT_DIR)
OUTPUT_DIR="${DEFAULT_OUTPUT_DIR}"

# Override with command-line if specified
FORCE_RERUN="${FORCE_RERUN:-${FORCE_RERUN:-false}}"

# ============================================================================
# Logging Setup
# ============================================================================

RHG_NAME=$(basename "$RHG_INPUT_DIR")
RHG_OUTPUT_DIR="${OUTPUT_DIR}/${RHG_NAME}"
LOG_DIR="${RHG_OUTPUT_DIR}/logs"

mkdir -p "$LOG_DIR"
MASTER_LOG="${LOG_DIR}/master_pipeline.log"

log() {
    local level=$1
    shift
    local msg="$@"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] [${level}] ${msg}" | tee -a "$MASTER_LOG"
}

# ============================================================================
# Validation
# ============================================================================

log "INFO" "Starting pipeline for RHG: $RHG_NAME"
log "INFO" "Input: $RHG_INPUT_DIR"
log "INFO" "Output: $RHG_OUTPUT_DIR"

if [[ "$VALIDATE_INPUTS" == "true" ]]; then
    log "INFO" "Validating inputs..."
    if ! "$PYTHON_BIN" "$SRC_DIR/validate_inputs.py" "$RHG_INPUT_DIR"; then
        if [[ "$VALIDATE_STRICT" == "true" ]]; then
            log "ERROR" "Input validation failed (strict mode)"
            exit 1
        else
            log "WARN" "Input validation failed (continuing anyway)"
        fi
    fi
fi

# ============================================================================
# Run Steps
# ============================================================================

run_step() {
    local step_num=$1
    local step_name=$2
    local script_path=$3
    local run_flag=$4

    if [[ "$run_flag" != "true" ]]; then
        log "INFO" "Step $step_num ($step_name) - SKIPPED"
        return 0
    fi

    local step_output="${RHG_OUTPUT_DIR}/step${step_num}_${step_name}"
    mkdir -p "$step_output"

    log "INFO" "===== STEP $step_num: $step_name ====="
    log "INFO" "Starting step script: $script_path"

    if bash "$script_path" \
        -i "$RHG_INPUT_DIR" \
        -o "$step_output" \
        -c "$CONFIG_FILE" \
        ${FORCE_RERUN:+-f}; then
        log "INFO" "Step $step_num COMPLETED successfully"
    else
        log "ERROR" "Step $step_num FAILED"
        exit 1
    fi
}

# Run the 6 steps in sequence
run_step 1 "branch_lengths" \
    "${SCRIPT_DIR}/01_estimate_branch_lengths.sh" \
    "$RUN_STEP1_BRANCH_LENGTHS"

run_step 2 "seqgen" \
    "${SCRIPT_DIR}/02_run_seqgen.sh" \
    "$RUN_STEP2_SEQGEN"

run_step 3 "masking" \
    "${SCRIPT_DIR}/03_mask_alignments.sh" \
    "$RUN_STEP3_MASKING"

run_step 4 "tree_reconstruction" \
    "${SCRIPT_DIR}/04_reconstruct_trees.sh" \
    "$RUN_STEP4_TREE_RECONSTRUCTION"

run_step 5 "au_test" \
    "${SCRIPT_DIR}/05_run_au_test.sh" \
    "$RUN_STEP5_AU_TEST"

run_step 6 "summarize" \
    "${SCRIPT_DIR}/06_summarise_results.sh" \
    "$RUN_STEP6_SUMMARIZE"

# ============================================================================
# Completion
# ============================================================================

log "INFO" "===== PIPELINE COMPLETED SUCCESSFULLY ====="
log "INFO" "Results directory: $RHG_OUTPUT_DIR"
log "INFO" "Summary: $RHG_OUTPUT_DIR/results.tsv"

echo ""
echo "Pipeline completed for RHG: $RHG_NAME"
echo "Results saved to: $RHG_OUTPUT_DIR"
