#!/bin/bash
#SBATCH --job-name=real_unconstr_logl
#SBATCH --partition=k2-hipri
#SBATCH --array=1-20
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --output=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/real_unconstrained_logl/slurm_%A_%a.log
#SBATCH --error=/users/3053301/analyses/Mary-Peter-simulations/data/benchmarks/real_unconstrained_logl/slurm_%A_%a.err

# Compute the log-likelihood of Peter's real unconstrained tree for 9,816
# domains (every manifest domain with a matching real tree -- see
# data/manifest/real_unconstrained_targets.tsv), independent of Step 1 /
# simulation progress. Fixed-topology (-te) evaluation, not a search --
# confirmed cheap regardless of taxa count (0.56s at 26 taxa, 1.6s at
# 1955 taxa, the largest domain we have a tree for).
#
# 20 chunks of ~491 domains each, processed concurrently within each task
# (xargs -P 16) rather than one task per domain -- at ~1-2s of real work
# each, per-task SLURM overhead would otherwise dominate actual compute.
# Output: <family>/Domain<N>/real_data/real_unconstrained_logl.txt per
# domain -- independent of and does not touch anything Step 1-6 produces.

set -uo pipefail

PROJECT_DIR=/users/3053301/analyses/Mary-Peter-simulations/au-test-fp-simulation
DATA_DIR=/users/3053301/analyses/Mary-Peter-simulations/data
OUT_DIR="${DATA_DIR}/sim_outputs_bulk"
IQTREE_BIN=/mnt/scratch2/igfs-anaconda/conda-envs/iqtree/bin/iqtree2
CHUNK_TSV=$(printf "%s/manifest/real_unconstrained_chunks/chunk_%04d.tsv" "$DATA_DIR" "$SLURM_ARRAY_TASK_ID")

cd "$PROJECT_DIR"
echo "Host: $(hostname), array task: ${SLURM_ARRAY_TASK_ID}, allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown}"

bash scripts/real_unconstrained_logl_chunk.sh "$CHUNK_TSV" "$OUT_DIR" "$IQTREE_BIN" "${SLURM_CPUS_PER_TASK:-16}"
