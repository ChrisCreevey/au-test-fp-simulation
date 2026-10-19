#!/bin/bash
#SBATCH --job-name=peter_overlap_retry
#SBATCH --partition=k2-lowpri
#SBATCH --array=1-3
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=5-00:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/peter_overlap_oom_timeout_retry/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/peter_overlap_oom_timeout_retry/slurm_%A_%a.err

# Retry for the 3 families that didn't finish under peter_overlap_lowpri_array.slurm
# (job 10126981): F1518 (OOM, 8/10 Domain1 replicates done), F1134 and
# F217 (TIMEOUT at the 5-day cap, 9/10 Domain1 replicates done).
#
# Root causes, confirmed from real sacct/log data before this submission,
# not guessed:
#   - F1518: MaxRSS=31.68GB against the original 32G request -- genuinely
#     maxed out. REPLICATE_JOBS=10 ran all 10 replicates concurrently on a
#     huge alignment (2004 taxa x 1168 sites, LG+F+R10), and 2 of those 10
#     concurrent IQTree processes got OOM-killed.
#   - F1134/F217: MaxRSS only ~19-20GB -- NOT a memory problem. These ran
#     the full 5 days with 9/10 replicates done. IQTREE_THREADS_PER_JOB=1
#     is fixed at submission time, so even once only 1 replicate remained
#     active, it stayed stuck on a single core -- 15 of 16 allocated CPUs
#     sat idle for what was likely the slowest, most important stretch.
#
# Fix for this retry: since only 1-2 replicates are actually missing per
# family now (not 10), reallocate the same 16-CPU budget into fewer,
# much faster workers instead of mostly-idle concurrency slots --
# REPLICATE_JOBS=2, IQTREE_THREADS_PER_JOB=8 (2x8=16). --mem bumped to
# 40G as cheap insurance, though concurrency dropping from 10 to ~2
# should relieve the memory pressure on its own.
#
# Run WITHOUT -x: Step 4's .iqtree-based skip-check (fixed 2026-10-02)
# auto-detects and redoes ONLY the specific missing replicates; nothing
# already complete is repeated.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
MANIFEST="${DATA_DIR}/manifest/manifest_full.json"
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${DATA_DIR}/manifest/peter_overlap_oom_timeout_retry.txt")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, family: ${FAMILY}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

export REPLICATE_JOBS=2
export IQTREE_THREADS_PER_JOB=8

bash scripts/run_family.sh -f "$FAMILY" -m "$MANIFEST" -o "$OUT_DIR" -c config/pipeline.conf
