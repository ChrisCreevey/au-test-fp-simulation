#!/bin/bash
set -euo pipefail

# Step 2: Simulate sequences for ONE domain of ONE family, using AliSim
# (IQTree's own native simulator) instead of SeqGen, with per-site rate
# variation CONCENTRATED at the same real positions the alignment shows it.
#
# Why AliSim, not SeqGen:
#   - SeqGen requires approximating FreeRate (+R) models as Gamma (+G), and
#     combined with an earlier Step 1 bug that also refit BRANCH LENGTHS
#     under that same substitute model, this inflated total tree length by
#     ~33% (worst branch +42%) for +R domains -- confirmed by direct
#     A/B test against the real published branch lengths.
#   - AliSim supports FreeRate natively, so Step 1 can now hand it IQTree's
#     own fully-parameterised fitted model string directly with zero
#     approximation.
#   - AliSim needs no consensus/ancestral-sequence file (calculate_consensus.py
#     is no longer used here): it draws its own root stochastically from the
#     model, which needs no real starting sequence to be representative.
#
# Why per-site concentration (not just one whole-domain AliSim call): even
# with correct branch lengths and a correctly-fitted global alpha/pinvar,
# simulated pairwise identity still came out far below real data's (e.g.
# 14.6% vs 47.2% for one real domain) -- a single global rate-heterogeneity
# summary spreads variation more evenly across sites than real alignments
# actually show. Fix: use Step 1's -wsr per-site rate-category assignment
# (src/simulate_by_site_rate.py) to simulate each rate category as its own
# uniform-rate block (via AliSim's --branch-scale) and stitch results back
# at their real positions -- reproducing the real conservation PATTERN, not
# just its average. See that script's docstring for full detail.
#
# CAUTION: do not pass -s (an alignment) to --alisim. -s does not just
# enable gap-copying -- it triggers a full, independent IQTree tree search +
# model inference on that alignment first (confirmed: 49 minutes wall-clock
# for one domain, and it crashed combined with --site-rate SAMPLING and a
# custom -t tree -- an internal IQTree bug). We don't need -s: Step 1 already
# fit the tree+model+site-rates we want, and gap masking is handled
# separately in Step 3.
#
# Input:  domain alignment (for its length only), Step 1's tree + fitted
#         model + site rates, a seed (shared across a family's domains for
#         paired replicates)
# Output: simulated_domain.fa (N_REPLICATES datasets, ungapped, one file)

# ============================================================================
# Setup
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -a <alignment> -s <step1_output_dir> -z <seed> -o <output_dir> -c <config_file> [-f]

Required:
  -a <file>   Domain alignment (e.g. F100058_Domain1.trim) -- used only for its length
  -s <dir>    Step 1 output directory (contains brlen.treefile, fitted_model.txt)
  -z <int>    Random seed for AliSim (shared across a family's domains for paired replicates)
  -o <dir>    Output directory for step 2 results
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
STEP1_OUTPUT_DIR=""
SEED=""
OUTPUT_DIR=""
CONFIG_FILE=""
FORCE_RERUN=false

while getopts "a:s:z:o:c:fh" opt; do
    case $opt in
        a) ALIGNMENT="$OPTARG" ;;
        s) STEP1_OUTPUT_DIR="$OPTARG" ;;
        z) SEED="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) CONFIG_FILE="$OPTARG" ;;
        f) FORCE_RERUN=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$ALIGNMENT" || -z "$STEP1_OUTPUT_DIR" || -z "$SEED" || -z "$OUTPUT_DIR" || -z "$CONFIG_FILE" ]]; then
    echo "Error: -a, -s, -z, -o, and -c are all required"
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

LOG="${OUTPUT_DIR}/step2.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"
}

log "Starting Step 2: Run AliSim"
log "Alignment: $ALIGNMENT"
log "Seed: $SEED"
log "Output directory: $OUTPUT_DIR"

TREE="${STEP1_OUTPUT_DIR}/brlen.treefile"
FITTED_MODEL_FILE="${STEP1_OUTPUT_DIR}/fitted_model.txt"
SITE_RATES_FILE="${STEP1_OUTPUT_DIR}/site_rates.txt"

for file in "$ALIGNMENT" "$TREE" "$FITTED_MODEL_FILE" "$SITE_RATES_FILE"; do
    if [[ ! -f "$file" ]]; then
        log "ERROR: Required file not found: $file"
        exit 1
    fi
done

FITTED_MODEL=$(cat "$FITTED_MODEL_FILE")
# Base substitution model only (e.g. "LG" from "LG+I{0.012}+G4{1.103}") --
# rate heterogeneity is applied per-category via site-rate concentration
# instead, not via the +I/+G/+R suffix here.
BASE_MODEL="${FITTED_MODEL%%+*}"
log "Fitted model: $FITTED_MODEL (base substitution model: $BASE_MODEL)"

# ============================================================================
# Run per-site-rate-concentrated simulation
# ============================================================================

OUTPUT_PREFIX="${OUTPUT_DIR}/simulated_domain"
OUTPUT_SIM="${OUTPUT_PREFIX}.fa"
WORK_DIR="${OUTPUT_DIR}/site_rate_categories"

if [[ -f "$OUTPUT_SIM" && "$FORCE_RERUN" != "true" ]]; then
    log "Simulated sequences already exist, skipping"
else
    log "Running per-site-rate-concentrated simulation..."

    if "$PYTHON_BIN" "$SRC_DIR/simulate_by_site_rate.py" \
        --tree "$TREE" \
        --base-model "$BASE_MODEL" \
        --full-model "$FITTED_MODEL" \
        --alignment "$ALIGNMENT" \
        --site-rates "$SITE_RATES_FILE" \
        --n-replicates "$N_REPLICATES" \
        --seed "$SEED" \
        --iqtree-bin "$IQTREE_BIN" \
        --work-dir "$WORK_DIR" \
        --output "$OUTPUT_SIM" \
        >> "$LOG" 2>&1; then
        log "Simulation completed: $OUTPUT_SIM"
    else
        log "ERROR: Per-site-rate simulation failed"
        exit 1
    fi
fi

# ============================================================================
# Verify Output
# ============================================================================

log "Verifying AliSim output..."

if [[ ! -f "$OUTPUT_SIM" ]]; then
    log "ERROR: Output file not created"
    exit 1
fi

N_RECORDS=$(grep -c '^>' "$OUTPUT_SIM" || true)
log "Output file: $N_RECORDS sequence records"

if [[ $N_RECORDS -lt 2 ]]; then
    log "ERROR: Output file is empty or too small"
    exit 1
fi

log "SUCCESS: Step 2 completed"
log "Simulated sequences ready for gap application in Step 3"
