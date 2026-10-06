#!/bin/bash
set -uo pipefail

# Process one chunk of the real-unconstrained-logl target list (see
# scripts/compute_real_unconstrained_logl.sh) concurrently -- each domain
# is cheap (~1-2s), so this runs many at once via xargs rather than
# chunking by family/time-risk the way the tree-building jobs do.
#
# Usage: real_unconstrained_logl_chunk.sh <chunk_tsv> <output_base_dir> <iqtree_bin> <jobs>
#   <chunk_tsv> columns: family, domain, alignment, model, tree (no header)

CHUNK_TSV="$1"
OUTPUT_BASE="$2"
IQTREE_BIN="$3"
JOBS="${4:-16}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

run_one_line() {
    local line="$1"
    local family domain alignment model tree
    IFS=$'\t' read -r family domain alignment model tree <<< "$line"
    local output_file="${OUTPUT_BASE}/${family}/Domain${domain}/real_data/real_unconstrained_logl.txt"
    bash "$SCRIPT_DIR/compute_real_unconstrained_logl.sh" "$alignment" "$model" "$tree" "$output_file" "$IQTREE_BIN"
}
export -f run_one_line
export SCRIPT_DIR OUTPUT_BASE IQTREE_BIN

n=$(wc -l < "$CHUNK_TSV")
echo "Chunk: $CHUNK_TSV ($n domains)"

cat "$CHUNK_TSV" | xargs -P "$JOBS" -I{} bash -c 'run_one_line "$1"' _ {}

echo "Chunk complete: $CHUNK_TSV"
