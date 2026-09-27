#!/bin/bash
set -euo pipefail

# Compute the REAL data's observed delta-logL (constrained vs unconstrained)
# for one domain -- the actual test statistic to compare against the null
# delta-logL distribution built from simulated replicates (see Step 6 /
# summarise_results.py's permutation-style comparison).
#
# This is NOT part of the per-replicate null-simulation loop (Steps 1-5).
# It runs ONCE per domain, on the REAL alignment, not simulated data:
#   - constrained logL: already computed by Step 1's fixed-topology fit
#     (brlen_estimate.iqtree) -- reused here, not recomputed
#   - unconstrained logL: requires a full, real (expensive) ML tree search
#     on the real alignment -- there is no shortcut for this; it's the same
#     order of cost as one simulated replicate's unconstrained search (see
#     the thread-scaling benchmarks in data/benchmarks/thread_scaling/)
#
# Input:  domain alignment, Step 1's output dir (for the constrained logL),
#         model string
# Output: <output_dir>/real_delta_logl.txt (CONSTRAINED_LOGL=,
#         UNCONSTRAINED_LOGL=, DELTA_LOGL=)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -a <alignment> -s <step1_output_dir> -m <model> -o <output_dir> -c <config_file> [-f]

Required:
  -a <file>   Domain alignment (e.g. F100058_Domain1.trim)
  -s <dir>    Step 1 output directory (for the already-fit constrained logL)
  -m <str>    IQTree model string for this domain (e.g. "LG+I+G4")
  -o <dir>    Output directory
  -c <file>   Config file

Optional:
  -f          Force re-run
  -h          Show help
EOF
    exit 1
}

ALIGNMENT="" ; STEP1_OUTPUT_DIR="" ; MODEL="" ; OUTPUT_DIR="" ; CONFIG_FILE="" ; FORCE_RERUN=false

while getopts "a:s:m:o:c:fh" opt; do
    case $opt in
        a) ALIGNMENT="$OPTARG" ;;
        s) STEP1_OUTPUT_DIR="$OPTARG" ;;
        m) MODEL="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) CONFIG_FILE="$OPTARG" ;;
        f) FORCE_RERUN=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$ALIGNMENT" || -z "$STEP1_OUTPUT_DIR" || -z "$MODEL" || -z "$OUTPUT_DIR" || -z "$CONFIG_FILE" ]]; then
    echo "Error: -a, -s, -m, -o, and -c are all required"
    usage
fi

CLI_FORCE_RERUN="$FORCE_RERUN"
source "$CONFIG_FILE"
if [[ "$CLI_FORCE_RERUN" == "true" ]]; then
    FORCE_RERUN=true
fi
mkdir -p "$OUTPUT_DIR"

LOG="${OUTPUT_DIR}/real_delta_logl.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"; }

RESULT_FILE="${OUTPUT_DIR}/real_delta_logl.txt"

if [[ -f "$RESULT_FILE" && "$FORCE_RERUN" != "true" ]]; then
    log "Already computed, skipping (use -f to force)"
    log "SUCCESS: skipped"
    exit 0
fi

CONSTRAINED_IQTREE="${STEP1_OUTPUT_DIR}/brlen_estimate.iqtree"
if [[ ! -f "$CONSTRAINED_IQTREE" ]]; then
    log "ERROR: Step 1's brlen_estimate.iqtree not found at $CONSTRAINED_IQTREE (run Step 1 first)"
    exit 1
fi

CONSTRAINED_LOGL=$(grep -m1 "Log-likelihood of the tree" "$CONSTRAINED_IQTREE" | grep -oP '(?<=tree: )[-0-9.]+')
if [[ -z "$CONSTRAINED_LOGL" ]]; then
    log "ERROR: Could not extract constrained logL from $CONSTRAINED_IQTREE"
    exit 1
fi
log "Constrained logL (from Step 1, real alignment): $CONSTRAINED_LOGL"

log "Running full unconstrained ML search on real alignment (this is the expensive part -- one real tree search, no shortcut)..."
PREFIX="${OUTPUT_DIR}/real_unconstrained"
$IQTREE_BIN \
    -s "$ALIGNMENT" \
    -m "$MODEL" \
    -T "${IQTREE_THREADS_PER_JOB:-1}" \
    -redo \
    --prefix "$PREFIX" \
    >> "$LOG" 2>&1

if [[ ! -f "${PREFIX}.iqtree" ]]; then
    log "ERROR: IQTree failed to complete the unconstrained search"
    exit 1
fi

UNCONSTRAINED_LOGL=$(grep -m1 "Log-likelihood of the tree" "${PREFIX}.iqtree" | grep -oP '(?<=tree: )[-0-9.]+')
if [[ -z "$UNCONSTRAINED_LOGL" ]]; then
    log "ERROR: Could not extract unconstrained logL from ${PREFIX}.iqtree"
    exit 1
fi
log "Unconstrained logL (real alignment, free search): $UNCONSTRAINED_LOGL"

DELTA_LOGL=$(python3 -c "print(f'{$UNCONSTRAINED_LOGL - ($CONSTRAINED_LOGL):.4f}')")
log "Real delta-logL (unconstrained - constrained): $DELTA_LOGL"

{
    echo "CONSTRAINED_LOGL=${CONSTRAINED_LOGL}"
    echo "UNCONSTRAINED_LOGL=${UNCONSTRAINED_LOGL}"
    echo "DELTA_LOGL=${DELTA_LOGL}"
} > "$RESULT_FILE"

log "SUCCESS: written to $RESULT_FILE"
