#!/bin/bash
#SBATCH --job-name=gt1000_per_fam
#SBATCH --partition=k2-hipri
#SBATCH --array=1-389
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/gt1000_per_family/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/gt1000_per_family/slurm_%A_%a.err

# Best-effort opportunistic run for the 389 families with a domain >1000
# taxa -- the tier otherwise earmarked for Peter (see
# docs/ADAPTING_TO_YOUR_CLUSTER.md, family list already sent). ONE family
# per array task (Option A from that doc), not chunked: these families have
# enormous cost variance (1001-3158 taxa), so bundling several into one
# chunk risks a slow family starving others queued behind it in the same
# task -- per-family tasks let the scheduler handle that naturally instead.
#
# This is explicitly NOT expected to complete the whole tier. Many of
# these families -- especially towards the upper end of the range -- may
# need longer than the 3h cap for even a single replicate's unconstrained
# search (see docs/METHODS.md), which no amount of retrying fixes (IQTree
# doesn't checkpoint mid-search). Submitted now opportunistically because
# k2-hipri currently has unusually high free capacity; whatever doesn't
# finish within 3h stays unfinished (not a failure, see run_family_chunk.sh
# equivalent skip-logic -- a future run with a different time limit, e.g.
# on Peter's cluster, picks up cleanly from here with nothing wasted).
#
# Same 16-core / REPLICATE_JOBS=10 / single-threaded-per-replicate profile
# as the other tranches -- NOT the "fewer replicates, more threads each"
# idea floated in docs/ADAPTING_TO_YOUR_CLUSTER.md for this size range,
# since that's untested and this run is deliberately using the one
# configuration we've already validated at scale.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/manifest/families_gt1000.txt")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family.sh -f "$FAMILY" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf
