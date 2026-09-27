#!/bin/bash
set -euo pipefail

# Step 1: fit branch lengths AND the full rate-heterogeneity model for ONE
# domain of ONE family, ready for AliSim (Step 2).
#
# History of fixes to this step:
#   1. Originally just copied the real .constr0.treefile through unchanged
#      whenever it already had branch lengths, and Step 2 (SeqGen) guessed a
#      hardcoded alpha=1.0 with +I dropped entirely. Fixed by always running
#      a fixed-topology (-te) IQTree fit here to get real fitted parameters.
#   2. That fit initially refit any +R<n> (FreeRate) model as +G<n>, because
#      SeqGen has no FreeRate support. CONFIRMED BY DIRECT TEST to be a
#      second, separate bug: fitting under a simpler substitute model
#      systematically INFLATES ML branch lengths (measured: +33% total tree
#      length, +42% on the worst branch, for a real domain) because a model
#      too simple for the true rate heterogeneity has to explain the same
#      observed divergence by assuming a uniformly higher rate. Refitting
#      under the REAL native model reproduced the original branch lengths
#      almost exactly.
#   3. Switched Step 2 from SeqGen to AliSim (IQTree's own native simulator),
#      which supports FreeRate models directly -- so the +R->+G substitution
#      is no longer needed for simulation either. This step now fits under
#      the real, unmodified model string and hands Step 2 IQTree's own
#      fully-parameterised model string (e.g. "LG+R4{0.31,0.18,...}"), read
#      straight from the "ALISIM COMMAND" line IQTree itself prints -- exact
#      fitted parameters, no re-derivation, no approximation.
#
# A residual gap versus real data remains even with these fixes (a single
# global rate-heterogeneity model, however well fit, is still a coarser
# summary than the real alignment's exact per-site pattern) -- under
# investigation separately (site-rate concentration).
#
# Input:  one domain's alignment + constrained tree + model string
# Output: <output_dir>/brlen.treefile, <output_dir>/fitted_model.txt
#         (one line: IQTree's fully-parameterised model string for AliSim)

# ============================================================================
# Setup
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -a <alignment> -t <constrained_tree> -m <model> -o <output_dir> -c <config_file> [-f]

Required:
  -a <file>   Domain alignment (e.g. F100058_Domain1.trim)
  -t <file>   Constrained tree for this domain (e.g. F100058_Domain1.constr0.treefile)
  -m <str>    IQTree model string for this domain (e.g. "LG+I+G4" or "LG+R4")
  -o <dir>    Output directory for step 1 results
  -c <file>   Config file

Optional:
  -f          Force re-fit even if output already exists
  -h          Show help

EOF
    exit 1
}

# ============================================================================
# Parse Arguments
# ============================================================================

ALIGNMENT=""
CONSTRAINED_TREE=""
MODEL=""
OUTPUT_DIR=""
CONFIG_FILE=""
FORCE_RERUN=false

while getopts "a:t:m:o:c:fh" opt; do
    case $opt in
        a) ALIGNMENT="$OPTARG" ;;
        t) CONSTRAINED_TREE="$OPTARG" ;;
        m) MODEL="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) CONFIG_FILE="$OPTARG" ;;
        f) FORCE_RERUN=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$ALIGNMENT" || -z "$CONSTRAINED_TREE" || -z "$MODEL" || -z "$OUTPUT_DIR" || -z "$CONFIG_FILE" ]]; then
    echo "Error: -a, -t, -m, -o, and -c are all required"
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

LOG="${OUTPUT_DIR}/step1.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"
}

log "Starting Step 1: Fit branch lengths + rate-heterogeneity model"
log "Alignment: $ALIGNMENT"
log "Constrained tree: $CONSTRAINED_TREE"
log "Model: $MODEL"
log "Output directory: $OUTPUT_DIR"

for file in "$ALIGNMENT" "$CONSTRAINED_TREE"; do
    if [[ ! -f "$file" ]]; then
        log "ERROR: Required file not found: $file"
        exit 1
    fi
done

# ============================================================================
# Skip only if we've already fit this exact domain
# ============================================================================

OUTPUT_TREE="${OUTPUT_DIR}/brlen.treefile"
FITTED_MODEL_FILE="${OUTPUT_DIR}/fitted_model.txt"
SITE_RATES_FILE="${OUTPUT_DIR}/site_rates.txt"

if [[ -f "$OUTPUT_TREE" && -f "$FITTED_MODEL_FILE" && -f "$SITE_RATES_FILE" && "$FORCE_RERUN" != "true" ]]; then
    log "Branch-length tree and fitted model already exist, skipping (use -f to force)"
    log "SUCCESS: Step 1 skipped"
    exit 0
fi

# ============================================================================
# Fit Branch Lengths + Rate-Heterogeneity Parameters (native model, no
# +R->+G substitution -- see header comment for why that mattered)
# ============================================================================

log "Fitting branch lengths and rate parameters ($MODEL) on constrained tree..."

PREFIX="${OUTPUT_DIR}/brlen_estimate"

# -wsr (write site rates): per-site empirical Bayes rate-category assignment,
# needed by Step 2 to concentrate simulated rate variation at the SAME real
# positions rather than spreading it i.i.d. across all sites -- a single
# global alpha/pinvar is a coarse summary that under-reproduces real data's
# actual per-site conservation pattern (confirmed: simulated pairwise
# identity far below real data's, even with correctly-fit global parameters
# and correct branch lengths).
$IQTREE_BIN \
    -s "$ALIGNMENT" \
    -te "$CONSTRAINED_TREE" \
    -m "$MODEL" \
    -T "${IQTREE_THREADS_PER_JOB:-1}" \
    -keep-ident \
    -wsr \
    -redo \
    --prefix "$PREFIX" \
    >> "$LOG" 2>&1

if [[ ! -f "${PREFIX}.treefile" ]]; then
    log "ERROR: IQTree failed to estimate branch lengths"
    exit 1
fi

cp "${PREFIX}.treefile" "$OUTPUT_TREE"

if [[ ! -f "${PREFIX}.rate" ]]; then
    log "ERROR: IQTree did not produce a .rate file (-wsr output missing)"
    exit 1
fi
cp "${PREFIX}.rate" "$SITE_RATES_FILE"

# ============================================================================
# Extract IQTree's own fully-parameterised model string (the exact command
# it suggests for AliSim), e.g. LG+R4{0.31,0.18,0.48,0.83,...} or
# LG+I{0.012}+G4{1.103} -- no re-derivation of alpha/pinvar ourselves.
# ============================================================================

# IQTree prints this in the .log, not always the .iqtree report -- check both.
FITTED_MODEL=$(grep -A2 '^ALISIM COMMAND' "${PREFIX}.iqtree" "${PREFIX}.log" 2>/dev/null | grep -oP '(?<=-m ")[^"]+' | head -1 || true)

if [[ -z "$FITTED_MODEL" ]]; then
    log "ERROR: Could not extract fitted model string from ${PREFIX}.iqtree (no ALISIM COMMAND section found)"
    exit 1
fi

echo "$FITTED_MODEL" > "$FITTED_MODEL_FILE"

log "Fitted model: ${FITTED_MODEL}"
log "SUCCESS: Step 1 completed"
log "Branch-length tree + fitted model + site rates ready for Step 2: $OUTPUT_TREE, $FITTED_MODEL_FILE, $SITE_RATES_FILE"
