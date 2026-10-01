# Adapting this pipeline to your own cluster

This is for running the pipeline somewhere other than the cluster it was
developed on — different scheduler, different partitions, different data
location. Nothing about the pipeline itself needs to change; what differs
is how you invoke it and how you batch the work. This doc covers both.

See [README.md](../README.md) for what the pipeline does, and
[docs/INPUTS.md](INPUTS.md) for the exact input file format.

## 1. Prerequisites

- **IQTree2 ≥ 2.2**, with AliSim support (`iqtree2 --alisim` should print
  AliSim's help text, not an "unknown option" error — some older/stripped
  builds omit it).
- **Python 3.8+**, standard library only — nothing to `pip install`.
- **bash** and **xargs** (for replicate-level parallelism within a domain).

```bash
git clone <this repo's URL>
cd au-test-fp-simulation
iqtree2 --alisim 2>&1 | head -3   # sanity check AliSim is present
python3 --version
```

## 2. Point the pipeline at your own data

The pipeline does not ship with data — it reads whatever real RHG data
directory you give it, in the format described in
[docs/INPUTS.md](INPUTS.md) (`F#####_Domain{N}.trim` /
`.constr0.treefile` / `-geneID.tsv` + a shared `CHG_model.json`).

Build a manifest against your own copy of that data:

```bash
python3 src/build_manifest.py \
    --data-dir /path/to/your/conv_fusion_data \
    --output-json /path/to/your/manifest.json \
    --output-tsv  /path/to/your/manifest.tsv
```

This is the one step that's specific to where *your* data lives — do not
try to reuse someone else's `manifest.json`, since it embeds paths to
their copy of the data, not yours. Rebuilding it against your own
`--data-dir` takes seconds.

**Filtering to a size range** (e.g. if you've been asked to take a
specific taxa-count band): `manifest.tsv` has one row per
`(family, domain)` with an `n_taxa` column. To get family IDs where the
*largest* domain exceeds some threshold:

```bash
awk -F'\t' 'NR>1 {max[$1] = ($8>max[$1] ? $8 : max[$1])} END {for (f in max) if (max[f] > 1000) print f}' \
    /path/to/your/manifest.tsv | sort > large_families.txt
```

(column 8 is `n_taxa` in the current `manifest.tsv` layout — check the
header line if that's changed.)

## 3. Point the config at your tools

Edit `config/pipeline.conf`, or just export overrides before running
(every tool path supports `${VAR:-default}` so an exported value wins):

```bash
export IQTREE_BIN=/full/path/to/your/iqtree2
export PYTHON_BIN=python3          # only if not just "python3" on your system
```

`N_REPLICATES` (default 10) and `N_BOOTSTRAP` (default 1000, IQTree's
`-zb`) are the two parameters most likely worth revisiting — see
[docs/METHODS.md](METHODS.md) §5 for the statistical reasoning behind the
current defaults.

## 4. Two ways to run — pick based on your cluster's job-count tolerance

Both call the same underlying pipeline (`scripts/run_family.sh`); the
difference is purely how work gets batched into scheduler jobs. Every step
script skips work it finds already done, so either option is safe to
resubmit after a timeout or interruption without redoing finished work.

### Option A — one job per family (simpler, recommended if your cluster allows it)

Each family becomes its own array task. This is the easier option to
reason about, and — importantly for the larger families (see §5) — it
avoids a real risk with chunking: bundling several families into one job
means one unusually slow family can starve the others queued behind it in
the same task before the time limit hits. Per-family jobs don't have that
problem; the scheduler handles load-balancing naturally.

```bash
#!/bin/bash
#SBATCH --job-name=au_test_sim
#SBATCH --array=1-<N_FAMILIES>        # adapt to your family count
#SBATCH --partition=<your_partition>   # adapt
#SBATCH --time=<your_time_limit>       # adapt -- see sizing guidance in §5
#SBATCH --cpus-per-task=<cores>        # adapt -- see §5
#SBATCH --mem=<memory>                 # adapt -- see §5

FAMILY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" large_families.txt)

bash scripts/run_family.sh \
    -f "$FAMILY" \
    -m /path/to/your/manifest.json \
    -o /path/to/your/output_dir \
    -c config/pipeline.conf
```

Only worth avoiding if your cluster enforces a low cap on simultaneously
queued/running jobs per user — some do, and a few hundred individual jobs
can trip that.

### Option B — chunked batches (fewer, larger jobs)

Group multiple families per array task, processed sequentially within the
task. This is what the ≤200-taxa bulk run on our cluster used (see
`scripts/run_family_chunk.sh`, `hpc_templates/bulk_le200_array.slurm` for
a worked example — the `.slurm` file itself won't run as-is on a different
scheduler/partition layout, but the structure is a useful template).

```bash
# one family ID per line, split into files of N_PER_CHUNK families each
split -l <N_PER_CHUNK> -d -a 4 large_families.txt chunk_ --additional-suffix=.txt
```

```bash
#!/bin/bash
#SBATCH --job-name=au_test_sim_chunked
#SBATCH --array=1-<N_CHUNKS>
#SBATCH --partition=<your_partition>
#SBATCH --time=<your_time_limit>
#SBATCH --cpus-per-task=<cores>
#SBATCH --mem=<memory>

CHUNK_FILE=$(printf "chunk_%04d.txt" "$SLURM_ARRAY_TASK_ID")

bash scripts/run_family_chunk.sh "$CHUNK_FILE" \
    /path/to/your/manifest.json /path/to/your/output_dir config/pipeline.conf
```

After a batch finishes (or times out), find what's left:

```bash
python3 scripts/check_completion.py \
    --chunks-dir <dir containing your chunk_*.txt files> \
    --output-base /path/to/your/output_dir \
    --mopup-dir <dir for mop-up chunks>
# then resubmit the array it suggests
```

**If using chunking on a size range with high cost variance** (see §5 —
this matters more for larger families than it did for our ≤200-taxa
batch), keep chunks small (a handful of families, not dozens) so one slow
family doesn't eat the whole task's time budget before reaching the rest
of its chunk.

## 5. Sizing: cores, memory, time — calibrate, don't guess

Our own chunk-sizing for the ≤200-taxa tranche was based on a real timing
calibration sample (10 families spanning that size range), not a guess —
and even then, the first production run undershot (most chunks hit the
time cap rather than finishing). Do the same on your cluster before
committing to a batch:

```bash
# time a handful of real families spanning your size range first
for fam in <5-10 representative family IDs>; do
    start=$(date +%s)
    bash scripts/run_family.sh -f "$fam" -m manifest.json -o output_dir -c config/pipeline.conf
    echo "$fam: $(( $(date +%s) - start ))s"
done
```

What we know from our own testing, to set expectations before you
calibrate:

- **Step 4 (tree reconstruction) dominates cost**, and scales with
  replicate count divided by however many replicates run concurrently
  (`REPLICATE_JOBS`, which defaults to your job's allocated core count).
  See [docs/METHODS.md](METHODS.md) §6.
- **Families above ~1000 taxa are a different cost class, not just
  "somewhat slower."** In our benchmarking, a single replicate's
  unconstrained ML search for a >1000-taxa domain alone approached a 3-hour
  wall-clock budget on a single core. Budget accordingly — likely several
  hours per family at minimum, possibly needing a higher per-task time
  limit than your cluster's default cap allows, or splitting a family's
  own replicates across more than one job.
- **Per-replicate multithreading may be worth revisiting at this scale.**
  Our ≤200-taxa tranche used `IQTREE_THREADS_PER_JOB=1` with many
  concurrent replicates (`REPLICATE_JOBS` high) — multithreading a single
  small alignment isn't worth it. For 1000+ taxa alignments, the opposite
  trade may hold: fewer concurrent replicates, each multithreaded (e.g.
  `REPLICATE_JOBS=4`, `IQTREE_THREADS_PER_JOB=8` on a 32-core node) could
  outperform many single-threaded ones. Worth A/B testing on one real
  large family before committing to a setting for the whole batch.
- **Memory**: we ran the ≤200-taxa tranche at 8GB/task without issue;
  we don't have a measured figure for 1000+ taxa alignments. Monitor your
  first few large-family jobs and adjust.

## 6. Sending results back

Once a batch of families has run, each family's output directory contains
everything needed to merge back: `<family>/step6_family/family_summary.tsv`
is the per-family rollup; `<family>/Domain*/step6/results.tsv` is the
per-domain detail. To combine many families into one table:

```bash
bash scripts/combine_results.sh -o /path/to/your/output_dir
```

produces `combined_family_results.tsv`. Send that (or the whole output
directory, if intermediate files are wanted too) back for merging with
the rest of the dataset's results.
