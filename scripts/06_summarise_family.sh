#!/bin/bash
set -euo pipefail

# Step 6: Join a family's per-domain Step 5 results into family-level
# OR/AND consolidation rates (src/summarise_family.py).
#
# Per the project's "report raw, decide downstream" design (see
# docs/METHODS.md §3), Step 5 never forces a single-domain-per-family
# verdict -- this step is that downstream join, done AFTER all of a
# family's domains have been through Step 5, not baked into Step 5 itself.
#
# Discovers domains automatically from the family's output directory
# structure (<family>/Domain*/step5), the same convention Steps 1-5 use --
# no need to list them manually.
#
# Input:  a family's base output directory, containing one or more
#         Domain*/step5 subdirectories (each with a step5_results.txt)
# Output: <output_dir>/family_summary.tsv (one row, concatenatable across
#         families), <output_dir>/family_replicate_detail.tsv (full
#         per-replicate transparency)

# ============================================================================
# Setup
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -i <family_output_dir> -o <output_dir> -c <config_file> [-f]

Required:
  -i <dir>    Family's base output directory (e.g. data/sim_outputs/F100058),
              containing one or more Domain*/step5 subdirectories
  -o <dir>    Output directory for family-level summary
  -c <file>   Config file

Optional:
  -f          Force re-run
  -h          Show help

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
    echo "Error: -i, -o, and -c are all required"
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

FAMILY_NAME=$(basename "$INPUT_DIR")

log "Starting Step 6: Summarize family-level results"
log "Family: $FAMILY_NAME"
log "Input directory: $INPUT_DIR"
log "Output directory: $OUTPUT_DIR"

# ============================================================================
# Discover domains (any Domain*/step5 subdirectory with a step5_results.txt)
# ============================================================================

DOMAIN_ARGS=()
for step5_dir in "$INPUT_DIR"/Domain*/step5; do
    [[ -d "$step5_dir" ]] || continue
    if [[ ! -f "${step5_dir}/step5_results.txt" ]]; then
        log "WARN: $step5_dir has no step5_results.txt, skipping"
        continue
    fi
    domain_name=$(basename "$(dirname "$step5_dir")")
    DOMAIN_ARGS+=("--domain" "${domain_name}:${step5_dir}")
    log "Found domain: $domain_name ($step5_dir)"
done

if [[ ${#DOMAIN_ARGS[@]} -eq 0 ]]; then
    log "ERROR: No domains with completed Step 5 results found under $INPUT_DIR"
    exit 1
fi

# ============================================================================
# Skip if already done (unless forced)
# ============================================================================

SUMMARY_FILE="${OUTPUT_DIR}/family_summary.tsv"

if [[ -f "$SUMMARY_FILE" && "$FORCE_RERUN" != "true" ]]; then
    log "Family summary already exists, skipping (use -f to force)"
    log "SUCCESS: Step 6 skipped"
    exit 0
fi

# ============================================================================
# Run Family-Level Join
# ============================================================================

log "Running family-level summarization..."

"$PYTHON_BIN" "$SRC_DIR/summarise_family.py" \
    --family "$FAMILY_NAME" \
    "${DOMAIN_ARGS[@]}" \
    --output-dir "$OUTPUT_DIR" \
    --log "$LOG"

# ============================================================================
# Verify and Report Results
# ============================================================================

if [[ -f "$SUMMARY_FILE" ]]; then
    log ""
    log "===== FAMILY-LEVEL RESULTS FOR $FAMILY_NAME ====="
    cat "$SUMMARY_FILE" | tee -a "$LOG"
    log ""
    log "SUCCESS: Step 6 completed"
    log "Family summary saved to: $SUMMARY_FILE"
else
    log "ERROR: Family summary file was not created"
    exit 1
fi
