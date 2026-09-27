#!/bin/bash
set -euo pipefail

# Step 3: Apply indel-gap pattern from one domain's original alignment to
# its simulated replicates.
#
# Input:  the domain's original alignment (for gap positions) + Step 2's
#         AliSim output directory
# Output: rep_001_masked.aln, rep_002_masked.aln, ... (one per replicate)

# ============================================================================
# Setup
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -a <alignment> -s <step2_output_dir> -o <output_dir> -c <config_file> [-f]

Required:
  -a <file>   Domain alignment (e.g. F100058_Domain1.trim) -- source of gap pattern
  -s <dir>    Step 2 output directory (contains simulated_domain.fa)
  -o <dir>    Output directory for step 3 results
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

ALIGNMENT=""
STEP2_OUTPUT_DIR=""
OUTPUT_DIR=""
CONFIG_FILE=""
FORCE_RERUN=false

while getopts "a:s:o:c:fh" opt; do
    case $opt in
        a) ALIGNMENT="$OPTARG" ;;
        s) STEP2_OUTPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) CONFIG_FILE="$OPTARG" ;;
        f) FORCE_RERUN=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$ALIGNMENT" || -z "$STEP2_OUTPUT_DIR" || -z "$OUTPUT_DIR" || -z "$CONFIG_FILE" ]]; then
    echo "Error: -a, -s, -o, and -c are all required"
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

LOG="${OUTPUT_DIR}/step3.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"
}

log "Starting Step 3: Apply Gap Pattern from Original Domain Alignment"
log "Alignment: $ALIGNMENT"
log "Output directory: $OUTPUT_DIR"

# ============================================================================
# Check Inputs
# ============================================================================

SIM="${STEP2_OUTPUT_DIR}/simulated_domain.fa"

for file in "$ALIGNMENT" "$SIM"; do
    if [[ ! -f "$file" ]]; then
        log "ERROR: Required file not found: $file"
        exit 1
    fi
done

log "Original alignment: $ALIGNMENT"
log "Simulated sequences: $SIM"

# ============================================================================
# Skip if already done (unless forced)
# ============================================================================

existing_count=$(ls -1 "$OUTPUT_DIR"/rep_*_masked.aln 2>/dev/null | wc -l) || true
if [[ "$existing_count" -gt 0 && "$FORCE_RERUN" != "true" ]]; then
    log "Found $existing_count existing masked alignment(s), skipping (use -f to force)"
    log "SUCCESS: Step 3 skipped"
    exit 0
fi

# ============================================================================
# Call Python Masking Script
# ============================================================================

log "Running Python masking script..."

"$PYTHON_BIN" "$SRC_DIR/apply_gap_pattern.py" \
    --original-aln "$ALIGNMENT" \
    --simulated-seqs "$SIM" \
    --output-dir "$OUTPUT_DIR" \
    --n-replicates "$N_REPLICATES" \
    --log "$LOG"

if [[ $? -ne 0 ]]; then
    log "ERROR: Masking script failed"
    exit 1
fi

# ============================================================================
# Verify Output
# ============================================================================

log "Verifying output..."

expected_files=$N_REPLICATES
actual_files=$(ls -1 "$OUTPUT_DIR"/rep_*_masked.aln 2>/dev/null | wc -l)

log "Expected files: $expected_files"
log "Actual files: $actual_files"

if [[ $actual_files -lt $expected_files ]]; then
    log "WARN: Not all expected files were created"
fi

if [[ $actual_files -gt 0 ]]; then
    log "SUCCESS: Step 3 completed"
    log "Masked alignments ready for tree reconstruction in Step 4"
else
    log "ERROR: No output files created"
    exit 1
fi
