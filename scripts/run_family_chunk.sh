#!/bin/bash
set -uo pipefail
# NOT -e: one family failing/timing out shouldn't stop the rest of the
# chunk from being attempted -- a partial chunk is exactly what the
# completion-check + mop-up pass is designed to clean up afterward.

# Process every family listed in one chunk file via run_family.sh,
# sequentially. This is the per-array-task worker for the bulk <=200-taxa
# batch run -- see hpc_templates/bulk_le200_array.slurm.
#
# Usage: run_family_chunk.sh <chunk_file> <manifest.json> <output_base_dir> <config_file>

CHUNK_FILE="$1"
MANIFEST="$2"
OUTPUT_BASE="$3"
CONFIG_FILE="$4"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ ! -f "$CHUNK_FILE" ]]; then
    echo "ERROR: chunk file not found: $CHUNK_FILE"
    exit 1
fi

n_families=$(wc -l < "$CHUNK_FILE")
echo "Chunk: $CHUNK_FILE ($n_families families)"

i=0
while IFS= read -r fam; do
    [[ -z "$fam" ]] && continue
    i=$((i + 1))
    echo "--- [$i/$n_families] $fam ---"
    start=$(date +%s)
    if bash "$SCRIPT_DIR/run_family.sh" -f "$fam" -m "$MANIFEST" -o "$OUTPUT_BASE" -c "$CONFIG_FILE"; then
        end=$(date +%s)
        echo "TIMING $fam $((end - start))s OK"
    else
        end=$(date +%s)
        echo "TIMING $fam $((end - start))s FAILED"
    fi
done < "$CHUNK_FILE"

echo "Chunk complete: $CHUNK_FILE"
