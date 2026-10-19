#!/bin/bash
set -uo pipefail

# Build OUR OWN real unconstrained tree for ONE domain, for domains where
# Peter's own unconstrained tree was never fetched (no matching Newick
# file among the 10,555 trees from his Dropbox -- 879 of 10,714 manifest
# domains, see data/manifest/missing_real_tree.tsv). This is a genuine,
# single free ML tree search on the REAL alignment -- NOT a fixed-topology
# evaluation like compute_real_unconstrained_logl.sh (which needs an
# existing tree to evaluate; these domains don't have one).
#
# Deliberately minimal, matching exactly how the SIMULATED side's
# unconstrained search works (04_reconstruct_trees.sh): model string fixed
# (no ModelFinder -- the manifest's model, e.g. "LG+R6", is passed
# directly to -m, never "-m MFP"/"TEST"), no bootstrap (-bb/-alrt), single
# build (no replicates). Per explicit instruction: single build, no model
# estimation, no bootstraps.
#
# Caveat (see review_response.md, question A1): this is OUR free-search
# tree on the real alignment under the real model, not necessarily the
# same tree Peter's own AU test used -- there is no way to verify that
# from this codebase. It fills the Real_Unconstrained_LogL comparison gap
# for these domains, but with that caveat.
#
# Usage: build_real_unconstrained_tree.sh <alignment> <model> <output_prefix> <threads> <iqtree_bin>
#
# Skips if <output_prefix>.iqtree already exists (resumable, same
# convention as every other step in this pipeline).

ALIGNMENT="$1"
MODEL="$2"
OUTPUT_PREFIX="$3"
THREADS="${4:-1}"
IQTREE_BIN="${5:-iqtree2}"

if [[ -f "${OUTPUT_PREFIX}.iqtree" ]]; then
    exit 0
fi

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

"$IQTREE_BIN" \
    -s "$ALIGNMENT" \
    -m "$MODEL" \
    -T "$THREADS" \
    -redo \
    --prefix "$OUTPUT_PREFIX" \
    > "${OUTPUT_PREFIX}.run.log" 2>&1

if [[ ! -f "${OUTPUT_PREFIX}.iqtree" ]]; then
    echo "FAILED (no .iqtree report) $ALIGNMENT" >&2
    exit 1
fi

LOGL=$(grep -m1 "Log-likelihood of the tree" "${OUTPUT_PREFIX}.iqtree" | grep -oP '(?<=tree: )[-0-9.]+')
if [[ -z "$LOGL" ]]; then
    echo "FAILED (could not extract logL) $ALIGNMENT" >&2
    exit 1
fi

TREELEN=$(grep -m1 "Total tree length (sum of branch lengths)" "${OUTPUT_PREFIX}.iqtree" | grep -oP '[-0-9.]+$')
TREELEN="${TREELEN:-NA}"

{
    echo "BUILT_UNCONSTRAINED_LOGL=${LOGL}"
    echo "BUILT_UNCONSTRAINED_TREELEN=${TREELEN}"
} > "${OUTPUT_PREFIX}.txt"
echo "OK ${OUTPUT_PREFIX}"
