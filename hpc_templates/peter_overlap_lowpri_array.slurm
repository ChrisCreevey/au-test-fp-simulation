#!/bin/bash
#SBATCH --job-name=peter_overlap
#SBATCH --partition=k2-lowpri
#SBATCH --array=1-334
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=5-00:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/peter_overlap_lowpri/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/peter_overlap_lowpri/slurm_%A_%a.err

# Opportunistic run of the 334 families handed to Peter (data/manifest/
# peter_handoff_v2.tsv), on our own k2-lowpri queue (no time cap, same
# node pool as hipri/medpri, just lower scheduling priority -- so this
# only eats capacity nothing higher-priority wants right now). Not a
# commitment to finish these ourselves -- just grabbing whatever we can
# in parallel with Peter's own run, since every family either of us
# finishes is one less for the other.
#
# Order is REVERSED relative to what Peter was given (his list is sorted
# ascending by family ID; this one runs descending) -- if he processes his
# copy top-to-bottom as given, we start from the opposite end, so we
# converge toward the middle instead of duplicating the same families
# first. No coordination mechanism beyond this -- if he's already done
# some of these by the time a given task here starts, run_family.sh's
# normal skip-logic just means it finishes fast with nothing to do, not
# wasted work.
#
# One family per task (not chunked, given the extreme size range here:
# up to ~2200 taxa). --mem=32G given the OOM history earlier in this
# project for large domains at 8G.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/manifest/peter_overlap_reversed.txt")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=10
export IQTREE_THREADS_PER_JOB=1

bash scripts/run_family.sh -f "$FAMILY" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf
