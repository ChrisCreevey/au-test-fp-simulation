#!/bin/bash
set -euo pipefail

# Step 6a: Summarize AU test results across all replicates for ONE domain
# (src/summarise_results.py) -- one row, concatenatable across every
# domain in the dataset for cross-gene analysis.
#
# Renamed from the old 06_summarise_results.sh: that version auto-discovered
# its Step 5 input via `ls $PARENT_OUTPUT_DIR/step5_*`, a naming convention
# (step5_branch_lengths-style suffixes) from the old single-tree-per-family
# layout. The current per-domain layout uses plain `step5` (no suffix)
# under <family>/<domain>/, so that glob silently matched nothing. Fixed by
# taking the Step 5 directory directly, the same convention Steps 4/5 use,
# instead of guessing a sibling path.
#
# For the cross-domain OR/AND join (once a family's domains all have Step 5
# results), see 06_summarise_family.sh instead -- this script is per-domain
# only.
#
# Input:  one domain's Step 5 output directory (step5_results.txt)
# Output: <output_dir>/results.tsv, <output_dir>/summary.txt

# ============================================================================
# Setup
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -s <step5_output_dir> -o <output_dir> -c <config_file> [-n <name>] [-f]

Required:
  -s <dir>    Step 5 output directory (contains step5_results.txt)
  -o <dir>    Output directory for step 6 results
  -c <file>   Config file

Optional:
  -n <name>   Identifier for this domain (default: derived from the path as
              <family>_<domain>, e.g. F100058_Domain1)
  -f          Force re-run
  -h          Show help

EOF
    exit 1
}

# ============================================================================
# Parse Arguments
# ============================================================================

STEP5_OUTPUT_DIR=""
OUTPUT_DIR=""
CONFIG_FILE=""
NAME=""
FORCE_RERUN=false

while getopts "s:o:c:n:fh" opt; do
    case $opt in
        s) STEP5_OUTPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) CONFIG_FILE="$OPTARG" ;;
        n) NAME="$OPTARG" ;;
        f) FORCE_RERUN=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$STEP5_OUTPUT_DIR" || -z "$OUTPUT_DIR" || -z "$CONFIG_FILE" ]]; then
    echo "Error: -s, -o, and -c are all required"
    usage
fi

# ============================================================================
# Load Configuration
# ============================================================================

CLI_FORCE_RERUN="$FORCE_RERUN"
source "$CONFIG_FILE"
if [[ "$CLI_FORCE_RERUN" == "true" ]]; then
    FORCE_RERUN=true
fi
mkdir -p "$OUTPUT_DIR"

LOG="${OUTPUT_DIR}/step6.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"
}

# Derive <family>_<domain> from .../sim_outputs/<family>/<domain>/step5
if [[ -z "$NAME" ]]; then
    DOMAIN_DIR=$(dirname "$STEP5_OUTPUT_DIR")
    FAMILY_DIR=$(dirname "$DOMAIN_DIR")
    NAME="$(basename "$FAMILY_DIR")_$(basename "$DOMAIN_DIR")"
fi

log "Starting Step 6: Summarize domain results"
log "Domain identifier: $NAME"
log "Step 5 output: $STEP5_OUTPUT_DIR"
log "Output directory: $OUTPUT_DIR"

if [[ ! -f "${STEP5_OUTPUT_DIR}/step5_results.txt" ]]; then
    log "ERROR: step5_results.txt not found in $STEP5_OUTPUT_DIR"
    exit 1
fi

# Real-data permutation comparison (scripts/real_data_delta_logl.sh), if it
# exists as a sibling of step5 -- optional, since it requires a separate,
# expensive real-data unconstrained search to have already been run.
DOMAIN_DIR=$(dirname "$STEP5_OUTPUT_DIR")
REAL_DELTA_LOGL_FILE="${DOMAIN_DIR}/real_data/real_delta_logl.txt"
if [[ -f "$REAL_DELTA_LOGL_FILE" ]]; then
    log "Found real-data comparison: $REAL_DELTA_LOGL_FILE"
else
    log "No real-data comparison found at $REAL_DELTA_LOGL_FILE (permutation columns will be NA)"
    REAL_DELTA_LOGL_FILE=""
fi

# ============================================================================
# Skip if already done (unless forced)
# ============================================================================

RESULTS_FILE="${OUTPUT_DIR}/results.tsv"

if [[ -f "$RESULTS_FILE" && "$FORCE_RERUN" != "true" ]]; then
    log "Results already exist, skipping (use -f to force)"
    log "SUCCESS: Step 6 skipped"
    exit 0
fi

# ============================================================================
# Run Python Summarization Script
# ============================================================================

log "Running Python summarization script..."

"$PYTHON_BIN" "$SRC_DIR/summarise_results.py" \
    --step5-output "$STEP5_OUTPUT_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --rhg-id "$NAME" \
    ${REAL_DELTA_LOGL_FILE:+--real-delta-logl-file "$REAL_DELTA_LOGL_FILE"} \
    --log "$LOG"

# ============================================================================
# Verify and Report Results
# ============================================================================

if [[ -f "$RESULTS_FILE" ]]; then
    log ""
    log "===== RESULTS FOR $NAME ====="
    cat "$RESULTS_FILE" | tee -a "$LOG"
    log ""
    log "SUCCESS: Step 6 completed"
    log "Results saved to: $RESULTS_FILE"
else
    log "ERROR: Results file was not created"
    exit 1
fi
