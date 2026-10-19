#!/bin/bash
set -euo pipefail

# Step 5: Run AU test on each replicate of ONE domain.
# For each replicate, run AU test comparing constrained vs unconstrained trees,
# and independently check composite-clade monophyly (rooting-invariant).
# Input:  rep_*/masked_unconstrained.treefile, rep_*/masked_constrained.treefile
#         (Step 4), rep_*_masked.aln (Step 3), the domain's -geneID.tsv
# Output: rep_*/au_test_result.* (AU test results), step5_results.txt (one
#         RESULT line per replicate: reject/no-reject, p-value, monophyly)

# ============================================================================
# Setup
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -s <step3_output_dir> -t <step4_output_dir> -m <model> -g <geneid_tsv> -o <output_dir> -c <config_file> [-f]

Required:
  -s <dir>    Step 3 output directory (contains rep_*_masked.aln)
  -t <dir>    Step 4 output directory (contains rep_*/masked_{un,}constrained.treefile)
  -m <str>    IQTree model string for this domain (e.g. "LG+I+G4")
  -g <file>   Domain's -geneID.tsv (taxon <TAB> composite|parent)
  -o <dir>    Output directory for step 5 results
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
STEP4_OUTPUT_DIR=""
MODEL=""
GENEID_FILE=""
OUTPUT_DIR=""
CONFIG_FILE=""
FORCE_RERUN=false

while getopts "s:t:m:g:o:c:fh" opt; do
    case $opt in
        s) STEP3_OUTPUT_DIR="$OPTARG" ;;
        t) STEP4_OUTPUT_DIR="$OPTARG" ;;
        m) MODEL="$OPTARG" ;;
        g) GENEID_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        c) CONFIG_FILE="$OPTARG" ;;
        f) FORCE_RERUN=true ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$STEP3_OUTPUT_DIR" || -z "$STEP4_OUTPUT_DIR" || -z "$MODEL" || -z "$GENEID_FILE" || -z "$OUTPUT_DIR" || -z "$CONFIG_FILE" ]]; then
    echo "Error: -s, -t, -m, -g, -o, and -c are all required"
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

LOG="${OUTPUT_DIR}/step5.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"
}

log "Starting Step 5: Run AU Test"
log "Model: $MODEL"
log "Output directory: $OUTPUT_DIR"

if [[ ! -f "$GENEID_FILE" ]]; then
    log "ERROR: geneID file not found: $GENEID_FILE"
    exit 1
fi

# ============================================================================
# Read geneID.tsv to get composite (fused) taxa
# ============================================================================
#
# Real data uses a 2-way composite/parent vocabulary per domain (which parent
# a "parent" taxon belongs to is implicit -- it's whichever domain this file
# is for). This replaces the old 4-way Parent_gene_1/Parent_gene_2/fused/absent
# presence_absence.txt scheme, which doesn't exist in the real per-domain data.

COMPOSITE_TAXA=()
while IFS=$'\t' read -r taxon category; do
    if [[ "$category" == "composite" ]]; then
        COMPOSITE_TAXA+=("$taxon")
    fi
done < "$GENEID_FILE"

log "Composite (fused) taxa (count: ${#COMPOSITE_TAXA[@]}): ${COMPOSITE_TAXA[@]}"

if [[ ${#COMPOSITE_TAXA[@]} -eq 0 ]]; then
    log "ERROR: No composite taxa found in $GENEID_FILE"
    exit 1
fi

# ============================================================================
# Process Each Replicate (in parallel — replicates are independent)
# ============================================================================

JOBS="${REPLICATE_JOBS:-${HPC_THREADS:-1}}"
IQTREE_THREADS_PER_JOB="${IQTREE_THREADS_PER_JOB:-1}"
COMPOSITE_TAXA_STR="${COMPOSITE_TAXA[*]}"

process_one_replicate_au() {
    local rep_dir="$1"
    local rep_name
    rep_name=$(basename "$rep_dir")

    local REP_OUTPUT_DIR="${OUTPUT_DIR}/${rep_name}"
    mkdir -p "$REP_OUTPUT_DIR"
    local REP_LOG="${REP_OUTPUT_DIR}/step5_rep.log"

    # ====================================================================
    # Get alignment and tree files
    # ====================================================================

    local MASKED_ALN="${STEP3_OUTPUT_DIR}/${rep_name}_masked.aln"
    if [[ ! -f "$MASKED_ALN" ]]; then
        echo "SKIP $rep_name reason=no_masked_aln"
        return 0
    fi

    local UNCONSTRAINED="${rep_dir}/masked_unconstrained.treefile"
    local CONSTRAINED="${rep_dir}/masked_constrained.treefile"

    # .treefile is written incrementally DURING a search, so it exists even
    # for a truncated/killed mid-search run -- .iqtree is only written on
    # genuine completion. Checking .treefile alone (as this used to do) let
    # Step 5 silently run the AU test against a non-final, possibly
    # non-converged tree -- the same bug class fixed in Step 4's skip-check
    # (see 04_reconstruct_trees.sh), which this never got updated to match.
    if [[ ! -f "$UNCONSTRAINED" ]] || [[ ! -f "$CONSTRAINED" ]] \
        || [[ ! -f "${rep_dir}/masked_unconstrained.iqtree" ]] || [[ ! -f "${rep_dir}/masked_constrained.iqtree" ]]; then
        echo "SKIP $rep_name reason=missing_trees"
        return 0
    fi

    # ====================================================================
    # Run AU Test
    # ====================================================================

    local AU_OUT="${REP_OUTPUT_DIR}/au_test_result"
    local TREES_FOR_AU="${REP_OUTPUT_DIR}/trees_for_au.nwk"

    if [[ ! -f "${AU_OUT}.iqtree" || "$FORCE_RERUN" == "true" ]]; then
        # Combine both trees into one file for AU test comparison
        cat "$CONSTRAINED" "$UNCONSTRAINED" > "$TREES_FOR_AU"

        # Run AU test comparing the two topologies.
        #
        # -n 0 is CRITICAL: without it, IQTree treats -s + -m as "search for
        # the ML tree from scratch" and only evaluates the -z trees against
        # whatever it finds AFTER a full independent NNI search completes --
        # silently redoing work as expensive as Step 4's entire unconstrained
        # search, for every replicate, before ever reaching the actual (cheap)
        # AU comparison. Confirmed via IQTree's own documented usage for
        # topology tests: "-n 0: avoid tree search and skip to parameter
        # estimation and tree evaluation". This was a pre-existing bug (not
        # introduced by the domain-split rewrite) caught by a dedicated-node
        # SLURM benchmark timing out at 45 minutes on what should be a
        # few-seconds-per-replicate test.
        $IQTREE_BIN \
            -s "$MASKED_ALN" \
            -m "$MODEL" \
            -z "$TREES_FOR_AU" \
            -au \
            -n 0 \
            -zb "$N_BOOTSTRAP" \
            -T "$IQTREE_THREADS_PER_JOB" \
            -redo \
            --prefix "$AU_OUT" \
            > "$REP_LOG" 2>&1
    fi

    # Extract AU test p-value from IQTree output.
    #
    # trees_for_au.nwk is built as `cat CONSTRAINED UNCONSTRAINED`, so Tree ID
    # 1 = constrained, Tree ID 2 = unconstrained. The question this whole test
    # exists to answer -- does AU reject the monophyly constraint? -- is
    # answered by TREE 1's p-value, not tree 2's. Tree 2 (unconstrained) is by
    # construction the best/reference tree IQTree found, so its own p-AU sits
    # at ~1 essentially always; extracting it (as this line previously did,
    # matching on "2" instead of "1") meant the pipeline's core rejection
    # count was structurally guaranteed to read ~0% regardless of the true
    # rate -- a critical, pre-existing correctness bug, not a parsing quirk.
    local au_rejected=0
    local au_pval="NA"
    local delta_logl="NA"
    local constrained_logl="NA"
    local unconstrained_logl="NA"
    if [[ -f "${AU_OUT}.iqtree" ]]; then
        # Find tree comparison table: look for "  1" followed by a negative logL value (contains "-")
        # Extract the AU p-value (last purely-numeric field on that line, ignoring +/- confidence-set markers)
        au_pval=$(awk '/^[ \t]+1[ \t]+-[0-9]/ {
            for (i = NF; i >= 1; i--) {
                field = $i
                gsub(/[+-]$/, "", field)
                if (field ~ /^[0-9.]+$/) {
                    print field
                    exit
                }
            }
        }' "${AU_OUT}.iqtree" || echo "NA")

        if [[ "$au_pval" != "NA" && "$au_pval" != "" ]]; then
            # Check if p-value < 0.05 (reject constraint)
            if (( $(echo "$au_pval < 0.05" | bc -l) )); then
                au_rejected=1
            fi
        fi

        constrained_logl=$(awk '/^[ \t]+1[ \t]+-[0-9]/ { print $2; exit }' "${AU_OUT}.iqtree" || echo "NA")
        unconstrained_logl=$(awk '/^[ \t]+2[ \t]+-[0-9]/ { print $2; exit }' "${AU_OUT}.iqtree" || echo "NA")
        constrained_logl="${constrained_logl:-NA}"
        unconstrained_logl="${unconstrained_logl:-NA}"

        # delta_logl = unconstrained - constrained, via direct subtraction
        # of the two raw log-likelihoods above -- NOT field 3 of the tree-1
        # row (IQTree's own "deltaL" column), which a previous version of
        # this script used. That column is defined relative to whichever
        # of the two compared trees has the HIGHER likelihood, and reads
        # exactly 0 for that tree by definition -- so whenever the free
        # (unconstrained) search did worse than the constrained tree (which
        # happens: ML search isn't guaranteed to converge, especially on a
        # simulated-noise replicate), the true negative ΔlogL was silently
        # replaced with 0 instead of being recorded. Confirmed on a
        # 600-family sample: ~0.65% of null replicates affected this way,
        # worst case a -539 log-likelihood difference read as 0. Raw
        # subtraction has no such floor and matches how the real-data
        # equivalent (generate_failure_report.read_real_logls) is computed.
        if [[ "$constrained_logl" != "NA" && "$unconstrained_logl" != "NA" ]]; then
            delta_logl=$(echo "$unconstrained_logl - $constrained_logl" | bc -l)
        else
            delta_logl="NA"
        fi
    fi

    # ====================================================================
    # Check Monophyly of Composite Clade in Both Trees
    # ====================================================================

    local TAXA
    read -ra TAXA <<< "$COMPOSITE_TAXA_STR"

    # Check unconstrained tree (true/false/NA if untestable - see check_monophyly.py)
    local unconstrained_mono_json unconstrained_mono
    unconstrained_mono_json=$("$PYTHON_BIN" "$SRC_DIR/check_monophyly.py" \
        --tree "$UNCONSTRAINED" \
        --taxa "${TAXA[@]}" 2>/dev/null)
    if grep -q '"monophyletic": true' <<< "$unconstrained_mono_json"; then
        unconstrained_mono="true"
    elif grep -q '"monophyletic": false' <<< "$unconstrained_mono_json"; then
        unconstrained_mono="false"
    else
        unconstrained_mono="NA"
    fi

    # Check constrained tree (should always be monophyletic by definition)
    local constrained_mono_json constrained_mono
    constrained_mono_json=$("$PYTHON_BIN" "$SRC_DIR/check_monophyly.py" \
        --tree "$CONSTRAINED" \
        --taxa "${TAXA[@]}" 2>/dev/null)
    if grep -q '"monophyletic": true' <<< "$constrained_mono_json"; then
        constrained_mono="true"
    elif grep -q '"monophyletic": false' <<< "$constrained_mono_json"; then
        constrained_mono="false"
    else
        constrained_mono="NA"
    fi

    local reason="none"
    if [[ "$au_rejected" == "1" ]]; then
        if [[ "$unconstrained_mono" == "false" ]]; then
            # AU test rejected AND composite clade is NOT monophyletic
            # This is a true AU test "failure" (false positive of monophyly)
            reason="no_monophyly"
        else
            # AU test rejected but composite clade still monophyletic (or untestable)
            reason="other"
        fi
    fi

    # Alignment length, taxon count, and tree length from the unconstrained IQ-TREE report
    # (used by Step 6 to summarize alignment-length/evolutionary-rate effects across genes)
    local UNCONSTRAINED_IQTREE="${rep_dir}/masked_unconstrained.iqtree"
    local ntaxa="NA" nsites="NA" treelen="NA"
    if [[ -f "$UNCONSTRAINED_IQTREE" ]]; then
        local input_data_line
        input_data_line=$(grep -E 'Input data:' "$UNCONSTRAINED_IQTREE" | head -1)
        ntaxa=$(grep -oE '[0-9]+' <<< "$input_data_line" | sed -n '1p')
        nsites=$(grep -oE '[0-9]+' <<< "$input_data_line" | sed -n '2p')
        treelen=$(grep -oE '[0-9]+\.[0-9]+' <<< "$(grep -E 'Total tree length \(sum of branch lengths\)' "$UNCONSTRAINED_IQTREE" | head -1)")
        ntaxa=${ntaxa:-NA}
        nsites=${nsites:-NA}
        treelen=${treelen:-NA}
    fi

    # Write monophyly check results
    {
        echo "Monophyly Check Results"
        echo "======================"
        echo "Composite taxa: ${TAXA[@]}"
        echo ""
        echo "Constrained tree (should be monophyletic): $constrained_mono"
        echo "Unconstrained tree: $unconstrained_mono"
        echo ""
        echo "AU test p-value: ${au_pval:-NA}"
        echo "AU test rejected constraint: $([[ $au_rejected -eq 1 ]] && echo 'YES' || echo 'NO')"
        echo ""
        if [[ "$au_rejected" == "1" ]]; then
            if [[ "$unconstrained_mono" == "false" ]]; then
                echo "AU REJECTION REASON: Composite clade lost monophyly (true signal of polyphyly)"
            else
                echo "AU REJECTION REASON: Topology changed, but composite clade still monophyletic"
            fi
        fi
    } > "${REP_OUTPUT_DIR}/monophyly_check.txt"

    echo "RESULT $rep_name rejected=$au_rejected reason=$reason mono=$unconstrained_mono ntaxa=$ntaxa nsites=$nsites treelen=$treelen pval=$au_pval deltalogl=$delta_logl constrained_logl=$constrained_logl unconstrained_logl=$unconstrained_logl"
}
export -f process_one_replicate_au
export IQTREE_BIN MODEL N_BOOTSTRAP FORCE_RERUN OUTPUT_DIR STEP3_OUTPUT_DIR PYTHON_BIN SRC_DIR IQTREE_THREADS_PER_JOB COMPOSITE_TAXA_STR

RESULTS_FILE="${OUTPUT_DIR}/step5_results.txt"
: > "$RESULTS_FILE"

log "Running AU test for replicates (up to $JOBS concurrently)"

ls -d "$STEP4_OUTPUT_DIR"/rep_* 2>/dev/null | \
    xargs -P "$JOBS" -I{} bash -c 'process_one_replicate_au "$@"' _ {} \
    >> "$RESULTS_FILE" || true

rep_count=$(grep -c '^RESULT' "$RESULTS_FILE" || true)
au_reject_count=$(grep '^RESULT' "$RESULTS_FILE" | grep -c 'rejected=1' || true)
au_reject_no_monophyly=$(grep '^RESULT' "$RESULTS_FILE" | grep -c 'reason=no_monophyly' || true)
au_reject_other_reasons=$(grep '^RESULT' "$RESULTS_FILE" | grep -c 'reason=other' || true)
skip_count=$(grep -c '^SKIP' "$RESULTS_FILE" || true)

if [[ $skip_count -gt 0 ]]; then
    log "WARN: $skip_count replicate(s) skipped (missing alignment/trees):"
    grep '^SKIP' "$RESULTS_FILE" | while read -r line; do log "  $line"; done
fi

log "Processed $rep_count replicates for AU testing"
log "AU test rejections (p < 0.05): $au_reject_count / $rep_count"
log "  - Rejections due to loss of monophyly: $au_reject_no_monophyly"
log "  - Rejections for other reasons: $au_reject_other_reasons"

# ============================================================================
# Calculate and Save Summary
# ============================================================================

SUMMARY_FILE="${OUTPUT_DIR}/au_test_summary.txt"

{
    echo "AU Test Summary with Monophyly Checking"
    echo "========================================"
    echo "Total replicates: $rep_count"
    echo ""
    echo "AU Test Results:"
    echo "  Total rejections (p < 0.05): $au_reject_count"
    echo "  Rejections with loss of monophyly: $au_reject_no_monophyly"
    echo "  Rejections for other reasons: $au_reject_other_reasons"

    if [[ $rep_count -gt 0 ]]; then
        fp_rate=$(echo "scale=4; 100 * $au_reject_count / $rep_count" | bc)
        echo ""
        echo "False Positive Rates:"
        echo "  Total FP rate (all rejections): ${fp_rate}%"

        if [[ $au_reject_no_monophyly -gt 0 ]]; then
            mono_fp_rate=$(echo "scale=4; 100 * $au_reject_no_monophyly / $rep_count" | bc)
            echo "  Monophyly-loss FP rate (true signal): ${mono_fp_rate}%"
        fi

        if [[ $au_reject_other_reasons -gt 0 ]]; then
            other_rate=$(echo "scale=4; 100 * $au_reject_other_reasons / $rep_count" | bc)
            echo "  Other topology changes: ${other_rate}%"
        fi
    fi

    echo ""
    echo "Interpretation:"
    echo "  - Expected FP rate for monophyly-loss ~5% if null model is correct"
    echo "  - Other topology changes may occur independently of monophyly"
    echo "  - Key metric: monophyly-loss FP rate, not total rejections"
} > "$SUMMARY_FILE"

log "Summary written to: $SUMMARY_FILE"
log "SUCCESS: Step 5 completed"
log "Results ready for summarization in Step 6"
