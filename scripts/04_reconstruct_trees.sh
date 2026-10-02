#!/bin/bash
set -euo pipefail

# Step 4: Reconstruct gene trees for each simulated replicate of ONE domain.
# For each replicate:
#   - Build an unconstrained ML tree from the masked alignment
#   - Evaluate the domain's constrained topology (composite clade forced
#     monophyletic) on the masked alignment, re-fitting branch lengths (-te)
# Input:  rep_001_masked.aln, ... (from Step 3) + the domain's constrained
#         tree topology + model
# Output: rep_*/masked_unconstrained.treefile, rep_*/masked_constrained.treefile

# ============================================================================
# Setup
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -s <step3_output_dir> -t <constrained_tree> -m <model> -o <output_dir> -c <config_file> [-f]

Required:
  -s <dir>    Step 3 output directory (contains rep_*_masked.aln)
  -t <file>   Domain's constrained tree topology (e.g. F100058_Domain1.constr0.treefile)
  -m <str>    IQTree model string for this domain (e.g. "LG+I+G4")
  -o <dir>    Output directory for step 4 results
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

STEP3_OUTPUT_DIR=""
CONSTRAINED_TREE=""
MODEL=""
OUTPUT_DIR=""
CONFIG_FILE=""
FORCE_RERUN=false

while getopts "s:t:m:o:c:fh" opt; do
    case $opt in
        s) STEP3_OUTPUT_DIR="$OPTARG" ;;
        t) CONSTRAINED_TREE="$OPTARG" ;;
        m) MODEL="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) CONFIG_FILE="$OPTARG" ;;
        f) FORCE_RERUN=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$STEP3_OUTPUT_DIR" || -z "$CONSTRAINED_TREE" || -z "$MODEL" || -z "$OUTPUT_DIR" || -z "$CONFIG_FILE" ]]; then
    echo "Error: -s, -t, -m, -o, and -c are all required"
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

LOG="${OUTPUT_DIR}/step4.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"
}

log "Starting Step 4: Reconstruct Gene Trees"
log "Constrained tree: $CONSTRAINED_TREE"
log "Model: $MODEL"
log "Output directory: $OUTPUT_DIR"

if [[ ! -f "$CONSTRAINED_TREE" ]]; then
    log "ERROR: Constrained tree not found: $CONSTRAINED_TREE"
    exit 1
fi

# ============================================================================
# Process Each Replicate (in parallel — replicates are independent)
# ============================================================================
#
# Note: the unconstrained tree search intentionally omits -bb/UFBoot. Bootstrap
# support values are not consumed anywhere downstream (monophyly checking only
# looks at tree topology, and the AU test in Step 5 uses its own -zb RELL
# bootstrap), so computing UFBoot here was pure wasted compute.

JOBS="${REPLICATE_JOBS:-${HPC_THREADS:-1}}"
IQTREE_THREADS_PER_JOB="${IQTREE_THREADS_PER_JOB:-1}"

total_reps=$(ls -1 "$STEP3_OUTPUT_DIR"/rep_*_masked.aln 2>/dev/null | wc -l) || true
log "Found $total_reps replicates to process (running up to $JOBS concurrently)"

process_one_replicate() {
    local masked_aln="$1"
    local rep_name
    rep_name=$(basename "$masked_aln" _masked.aln)

    local REP_OUTPUT_DIR="${OUTPUT_DIR}/${rep_name}"
    mkdir -p "$REP_OUTPUT_DIR"
    local REP_LOG="${REP_OUTPUT_DIR}/step4_rep.log"

    # Skip-if-done must check for the FINAL .iqtree report, not .treefile --
    # IQTree updates .treefile incrementally as it finds better trees during
    # search, but only writes .iqtree on genuine completion. A run killed
    # mid-search by a time limit leaves a .treefile (whatever intermediate
    # tree it had reached) with no .iqtree; checking .treefile alone wrongly
    # treats that truncated search as done on any later resubmission,
    # silently locking in a non-converged tree (confirmed cause of ~160
    # family-domains with NA taxa/alignment-length metadata and, in the
    # worst cases, a degenerate always-zero delta-logL -- see project
    # history / docs/METHODS.md).
    local UNCONSTRAINED="${REP_OUTPUT_DIR}/masked_unconstrained"
    if [[ ! -f "${UNCONSTRAINED}.iqtree" || "$FORCE_RERUN" == "true" ]]; then
        $IQTREE_BIN \
            -s "$masked_aln" \
            -m "$MODEL" \
            -T "$IQTREE_THREADS_PER_JOB" \
            -redo \
            --prefix "$UNCONSTRAINED" \
            > "$REP_LOG" 2>&1

        if [[ ! -f "${UNCONSTRAINED}.iqtree" ]]; then
            echo "FAILED $rep_name unconstrained"
            return 1
        fi
    fi

    local CONSTRAINED="${REP_OUTPUT_DIR}/masked_constrained"
    if [[ ! -f "${CONSTRAINED}.iqtree" || "$FORCE_RERUN" == "true" ]]; then
        $IQTREE_BIN \
            -s "$masked_aln" \
            -te "$CONSTRAINED_TREE" \
            -m "$MODEL" \
            -T "$IQTREE_THREADS_PER_JOB" \
            -redo \
            --prefix "$CONSTRAINED" \
            >> "$REP_LOG" 2>&1

        if [[ ! -f "${CONSTRAINED}.iqtree" ]]; then
            echo "FAILED $rep_name constrained"
            return 1
        fi
    fi

    echo "OK $rep_name"
}
export -f process_one_replicate
export IQTREE_BIN MODEL CONSTRAINED_TREE FORCE_RERUN OUTPUT_DIR IQTREE_THREADS_PER_JOB

RESULTS_FILE="${OUTPUT_DIR}/step4_results.txt"
: > "$RESULTS_FILE"

ls "$STEP3_OUTPUT_DIR"/rep_*_masked.aln 2>/dev/null | \
    xargs -P "$JOBS" -I{} bash -c 'process_one_replicate "$@"' _ {} \
    >> "$RESULTS_FILE" || true

rep_count=$(grep -c '^OK\|^FAILED' "$RESULTS_FILE" || true)
fail_count=$(grep -c '^FAILED' "$RESULTS_FILE" || true)

if [[ $fail_count -gt 0 ]]; then
    log "ERROR: $fail_count of $rep_count replicate(s) failed tree reconstruction:"
    grep '^FAILED' "$RESULTS_FILE" | while read -r line; do log "  $line"; done
    exit 1
fi

log "Processed $rep_count replicates"
log "SUCCESS: Step 4 completed"
log "Gene trees ready for AU testing in Step 5"
