#!/bin/bash
set -uo pipefail

# Compute the log-likelihood of Peter's real unconstrained tree for ONE
# domain, under the real alignment and that domain's real model. This is
# a fixed-topology evaluation (-te), not a search -- cheap regardless of
# taxa count (confirmed: 0.56s at 26 taxa, 1.6s at 1955 taxa). Completely
# independent of Step 1 / the simulation pipeline's constrained-tree fit;
# the only output is the unconstrained tree's own logL. Joining this
# against a constrained logL (from Step 1 if available, or computed the
# same cheap way) to get delta-logL is a separate step, not done here.
#
# Usage: compute_real_unconstrained_logl.sh <alignment> <model> <tree> <output_file> <iqtree_bin>
#
# Skips if <output_file> already exists (resumable, same convention as
# every other step in this pipeline).

ALIGNMENT="$1"
MODEL="$2"
TREE="$3"
OUTPUT_FILE="$4"
IQTREE_BIN="$5"

if [[ -f "$OUTPUT_FILE" ]]; then
    exit 0
fi

OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

PREFIX="${OUTPUT_FILE%.txt}"

"$IQTREE_BIN" \
    -s "$ALIGNMENT" \
    -te "$TREE" \
    -m "$MODEL" \
    -T 1 \
    -redo \
    --prefix "$PREFIX" \
    > "${PREFIX}.run.log" 2>&1

if [[ ! -f "${PREFIX}.iqtree" ]]; then
    echo "FAILED (no .iqtree report) $ALIGNMENT" >&2
    exit 1
fi

LOGL=$(grep -m1 "Log-likelihood of the tree" "${PREFIX}.iqtree" | grep -oP '(?<=tree: )[-0-9.]+')
if [[ -z "$LOGL" ]]; then
    echo "FAILED (could not extract logL) $ALIGNMENT" >&2
    exit 1
fi

TREELEN=$(grep -m1 "Total tree length (sum of branch lengths)" "${PREFIX}.iqtree" | grep -oP '[-0-9.]+$')
TREELEN="${TREELEN:-NA}"

{
    echo "UNCONSTRAINED_LOGL=${LOGL}"
    echo "UNCONSTRAINED_TREELEN=${TREELEN}"
} > "$OUTPUT_FILE"
echo "OK $OUTPUT_FILE"
