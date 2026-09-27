#!/bin/bash
set -euo pipefail

# Combine every family's step6_family/family_summary.tsv (one per
# <output_base>/<family_id>/step6_family/, produced by run_family.sh's final
# family-level join) into one dataset-wide table.
#
# Thin wrapper around src/combine_family_results.py -- the actual join needs
# real logic (a column-name union, not a plain concatenation), since families
# with different domain counts produce different columns. See that script's
# docstring for why.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

usage() {
    cat <<EOF
Usage: $(basename "$0") -o <output_base_dir> [-out <combined_file>]

Required:
  -o <dir>      Output base directory (e.g. "data/sim_outputs_bulk"),
                containing one subdirectory per family

Optional:
  -out <file>   Path for the combined TSV (default: <output_base_dir>/combined_family_results.tsv)
  -h            Show help

EOF
    exit 1
}

OUTPUT_BASE=""
COMBINED_FILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o) OUTPUT_BASE="$2"; shift 2 ;;
        -out) COMBINED_FILE="$2"; shift 2 ;;
        -h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$OUTPUT_BASE" ]]; then
    usage
fi

if [[ -z "$COMBINED_FILE" ]]; then
    COMBINED_FILE="${OUTPUT_BASE}/combined_family_results.tsv"
fi

python3 "$PROJECT_DIR/src/combine_family_results.py" \
    --output-base "$OUTPUT_BASE" \
    --out "$COMBINED_FILE"
