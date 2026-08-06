#!/bin/bash
set -euo pipefail

# Combine per-gene-family results.tsv files (one per outputs/<RHG_ID>/step6_summarize/)
# into a single master TSV for cross-gene analysis.

usage() {
    cat <<EOF
Usage: $(basename "$0") -o <outputs_dir> [-out <combined_file>]

Required:
  -o <outputs_dir>      Parent outputs directory (e.g. "outputs"), containing
                         one subdirectory per gene family run

Optional:
  -out <combined_file>  Path for the combined TSV (default: <outputs_dir>/combined_results.tsv)
  -h                    Show help

EOF
    exit 1
}

OUTPUTS_DIR=""
COMBINED_FILE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o) OUTPUTS_DIR="$2"; shift 2 ;;
        -out) COMBINED_FILE="$2"; shift 2 ;;
        -h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$OUTPUTS_DIR" ]]; then
    usage
fi

if [[ -z "$COMBINED_FILE" ]]; then
    COMBINED_FILE="${OUTPUTS_DIR}/combined_results.tsv"
fi

RESULT_FILES=()
while IFS= read -r -d '' f; do
    RESULT_FILES+=("$f")
done < <(find "$OUTPUTS_DIR" -path "*/step6_summarize/results.tsv" -print0 | sort -z)

if [[ ${#RESULT_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No results.tsv files found under $OUTPUTS_DIR" >&2
    exit 1
fi

# Header from the first file, data rows (skip header) from every file
head -1 "${RESULT_FILES[0]}" > "$COMBINED_FILE"
for f in "${RESULT_FILES[@]}"; do
    tail -n +2 "$f" >> "$COMBINED_FILE"
done

n_genes=$(( $(wc -l < "$COMBINED_FILE") - 1 ))
echo "Combined ${#RESULT_FILES[@]} results.tsv file(s) -> $n_genes gene(s) written to: $COMBINED_FILE"
