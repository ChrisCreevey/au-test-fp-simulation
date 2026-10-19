#!/bin/bash
#SBATCH --job-name=build_real_unc
#SBATCH --partition=k2-bioinf,k2-lowpri
#SBATCH --array=1-879
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=5-00:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/build_real_unconstrained_tree/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/build_real_unconstrained_tree/slurm_%A_%a.err

# Build our own real unconstrained tree for the 879 domains with no
# matching tree fetched from Peter's Dropbox (see
# data/manifest/missing_real_tree.tsv, data/manifest/build_real_unconstrained_targets.tsv).
# One task per domain (879 <= MaxArraySize=1001, no chunking needed --
# unlike the cheap fixed-topology case in real_unconstrained_logl_array.slurm,
# each of these is a genuine free ML search, not a ~1-2s evaluation, so
# packing many into one task via xargs would risk one task serially
# summing hours of work).
#
# k2-bioinf,k2-lowpri (whichever has a free slot first), generous uncapped
# time -- per the project's own established lesson
# (step4_surgical_lowpri_array.slurm), tiering this kind of search by
# taxa count alone isn't reliable (model complexity matters too), and
# this batch skews large anyway (median 850 taxa, 42% over 1000, max
# 3158) -- not worth re-deriving a tier split for a one-off job.
# --mem=32G matches the largest family size class already handled
# elsewhere in this project without OOM.
#
# --cpus-per-task=8, passed straight to IQTree's own -T: unlike the
# simulation steps (which run up to 10 REPLICATES concurrently and so
# keep each individual IQTree call single-threaded), this is ONE build
# per task with no sibling replicates to share cores with, so letting
# IQTree itself multithread is the only way to use the allocation.
#
# Deliberately minimal invocation (build_real_unconstrained_tree.sh):
# single build, fixed model (no ModelFinder), no bootstrap -- per
# explicit instruction.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
TARGETS="${DATA_DIR}/manifest/build_real_unconstrained_targets.tsv"

cd "$PROJECT_DIR"
source config/pipeline.conf

# +1 to skip the header line
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$TARGETS")
IFS=$'\t' read -r FAMILY DOMAIN ALIGNMENT MODEL N_TAXA <<< "$LINE"

echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, domain: ${DOMAIN}, n_taxa: ${N_TAXA}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

OUTPUT_PREFIX="${DATA_DIR}/sim_outputs_bulk/${FAMILY}/Domain${DOMAIN}/real_data/built_unconstrained"

bash scripts/build_real_unconstrained_tree.sh "$ALIGNMENT" "$MODEL" "$OUTPUT_PREFIX" "${SLURM_CPUS_PER_TASK:-8}" "$IQTREE_BIN"
