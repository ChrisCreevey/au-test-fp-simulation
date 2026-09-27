#!/bin/bash
set -euo pipefail

# Orchestrate the full pipeline (Steps 1-6) for ONE family across all of its
# domains, then the family-level OR/AND join. This is the entry point for
# the domain-split pipeline (see docs/METHODS.md for the architecture).
#
# Reads domain list, models, and the family's shared seed directly from the
# manifest (src/build_manifest.py output) -- no re-derivation.
#
# The real-data permutation comparison (scripts/real_data_delta_logl.sh) is
# OFF by default: it requires one additional expensive unconstrained search
# per domain, on top of everything else, and for a bulk/overview run across
# thousands of families the null-simulation reject/OR/AND rates are the
# priority. Pass -r to also run it per domain (recommended only for a
# smaller set of families you want to examine closely, not a bulk batch).
#
# Every underlying step script already skips work that's already done
# (unless -x forces a rerun), so a killed/resubmitted run_family.sh call
# resumes rather than redoing completed steps -- this is what makes it safe
# to batch many families per SLURM array task with a generous time margin.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -f <family_id> -m <manifest.json> -o <output_base_dir> -c <config_file> [-r] [-x]

Required:
  -f <id>     Family ID (e.g. F29912)
  -m <file>   manifest.json (src/build_manifest.py output)
  -o <dir>    Output base directory (per-family output goes in <dir>/<family_id>)
  -c <file>   Config file

Optional:
  -r          Also run the real-data delta-logL comparison per domain (expensive; off by default)
  -x          Force re-run of every step (default: skip steps already completed)
  -h          Show help
EOF
    exit 1
}

FAMILY="" ; MANIFEST="" ; OUTPUT_BASE="" ; CONFIG_FILE="" ; RUN_REAL_DATA=false ; FORCE_FLAG=""

while getopts "f:m:o:c:rxh" opt; do
    case $opt in
        f) FAMILY="$OPTARG" ;;
        m) MANIFEST="$OPTARG" ;;
        o) OUTPUT_BASE="$OPTARG" ;;
        c) CONFIG_FILE="$OPTARG" ;;
        r) RUN_REAL_DATA=true ;;
        x) FORCE_FLAG="-f" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$FAMILY" || -z "$MANIFEST" || -z "$OUTPUT_BASE" || -z "$CONFIG_FILE" ]]; then
    echo "Error: -f, -m, -o, and -c are all required"
    usage
fi

source "$CONFIG_FILE"

FAM_DIR="${OUTPUT_BASE}/${FAMILY}"
mkdir -p "$FAM_DIR"
LOG="${FAM_DIR}/run_family.log"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $@" | tee -a "$LOG"; }

log "Starting family: $FAMILY"

# ============================================================================
# Pull this family's domain list, models, paths, and shared seed from the
# manifest in one shot (avoids re-deriving anything build_manifest.py
# already worked out).
# ============================================================================

DOMAIN_INFO=$(python3 -c "
import json, sys
manifest = json.load(open('$MANIFEST'))
fam = manifest.get('$FAMILY')
if fam is None:
    sys.exit('FAMILY_NOT_FOUND')
print(f\"SEED={fam['seed']}\")
for dom, info in sorted(fam['domains'].items(), key=lambda kv: int(kv[0])):
    print(f\"DOMAIN\tDomain{dom}\t{info['trim']}\t{info['tree']}\t{info['model']}\t{info['geneid']}\")
")

if [[ "$DOMAIN_INFO" == "FAMILY_NOT_FOUND" ]]; then
    log "ERROR: $FAMILY not found in manifest $MANIFEST"
    exit 1
fi

SEED=$(echo "$DOMAIN_INFO" | grep '^SEED=' | cut -d= -f2)
log "Family seed: $SEED"

# ============================================================================
# Run Steps 1-6 for each domain
# ============================================================================

while IFS=$'\t' read -r tag domain trim tree model geneid; do
    [[ "$tag" == "DOMAIN" ]] || continue

    dom_dir="${FAM_DIR}/${domain}"
    mkdir -p "${dom_dir}/step1"
    log "=== ${domain} (model=${model}) ==="

    bash "$SCRIPT_DIR/01_estimate_branch_lengths.sh" \
        -a "$trim" -t "$tree" -m "$model" -o "${dom_dir}/step1" -c "$CONFIG_FILE" $FORCE_FLAG

    bash "$SCRIPT_DIR/02_run_alisim.sh" \
        -a "$trim" -s "${dom_dir}/step1" -z "$SEED" -o "${dom_dir}/step2" -c "$CONFIG_FILE" $FORCE_FLAG

    bash "$SCRIPT_DIR/03_mask_alignments.sh" \
        -a "$trim" -s "${dom_dir}/step2" -o "${dom_dir}/step3" -c "$CONFIG_FILE" $FORCE_FLAG

    bash "$SCRIPT_DIR/04_reconstruct_trees.sh" \
        -s "${dom_dir}/step3" -t "$tree" -m "$model" -o "${dom_dir}/step4" -c "$CONFIG_FILE" $FORCE_FLAG

    bash "$SCRIPT_DIR/05_run_au_test.sh" \
        -s "${dom_dir}/step3" -t "${dom_dir}/step4" -m "$model" -g "$geneid" -o "${dom_dir}/step5" -c "$CONFIG_FILE" $FORCE_FLAG

    if [[ "$RUN_REAL_DATA" == "true" ]]; then
        mkdir -p "${dom_dir}/real_data"
        bash "$SCRIPT_DIR/real_data_delta_logl.sh" \
            -a "$trim" -s "${dom_dir}/step1" -m "$model" -o "${dom_dir}/real_data" -c "$CONFIG_FILE" $FORCE_FLAG
    fi

    bash "$SCRIPT_DIR/06_summarise_domain.sh" \
        -s "${dom_dir}/step5" -o "${dom_dir}/step6" -c "$CONFIG_FILE" $FORCE_FLAG

done <<< "$DOMAIN_INFO"

# ============================================================================
# Family-level join
# ============================================================================

log "=== Family-level summary ==="
bash "$SCRIPT_DIR/06_summarise_family.sh" \
    -i "$FAM_DIR" -o "${FAM_DIR}/step6_family" -c "$CONFIG_FILE" $FORCE_FLAG

log "SUCCESS: $FAMILY completed"
