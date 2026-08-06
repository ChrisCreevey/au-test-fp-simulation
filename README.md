# AU Test False Positive Rate Simulation Pipeline

A computational pipeline to test whether the AU (Approximately Unbiased) test produces false positives when detecting convergently evolved remodelled (fused) gene families.

The pipeline simulates sequence data under the correct null model — single-origin gene fusion, where the fused clade is truly monophyletic — using branch lengths and gap patterns calibrated from real data. It then runs the AU test against the simulated replicates and measures how often the test incorrectly rejects monophyly, giving an empirical false positive rate to compare against observed rejection rates in real gene families.

## Dependencies

- **IQ-TREE 2.0+** — phylogenetic tree inference and the AU test
- **SeqGen** — sequence simulation
- **Python 3.8+** — data processing and summarization (standard library only, no extra packages required)
- **bash**, **xargs** — replicate parallelization (macOS/BSD and Linux/GNU both supported)
- **MAFFT** (optional) — only needed if re-aligning simulated sequences

See [docs/INSTALLATION.md](docs/INSTALLATION.md) for platform-specific install instructions (Homebrew on macOS, Conda on HPC).

## Setup

### 1. Prepare input data

For each gene family (RHG) you want to test, create a directory with these **4 required files**:

| File | Format | Description |
|------|--------|-------------|
| `combined.aln` | FASTA | Alignment of all taxa across both parent gene regions |
| `constrained.nwk` | Newick | Tree with the fused clade constrained to be monophyletic (branch lengths optional — Step 1 estimates them if missing) |
| `model.txt` | Text | One line with the substitution model, e.g. `WAG+G4` (amino acid) or `GTR+G4` (DNA) |
| `presence_absence.txt` | TSV | `species<TAB>gene_type`, where `gene_type` is one of `Parent_gene_1`, `Parent_gene_2`, `fused`, `absent` |

A worked example lives in [examples/](examples/); a self-contained copy is also in [examples/for_peter/](examples/for_peter/).

See [docs/INPUTS.md](docs/INPUTS.md) for the full format specification.

### 2. Using your own data

The fastest way to get real data into the right shape is to start from the provided template:

```bash
cp -r inputs/template inputs/<your_gene_family_name>
# then edit the 4 files inside with your real data
```

Each gene family (RHG) gets its own directory, wherever you like — `inputs/<name>/` is just the convention `master_pipeline.sh -i <dir>` and the templates use; the path is passed explicitly with `-i`, so it isn't hardcoded. For many gene families, create one directory per family, e.g. `inputs/RHG_001/`, `inputs/RHG_002/`, ...

**What each file needs to contain, precisely:**

- **`combined.aln`** — one alignment per gene family spanning *both* parent-gene regions. A taxon that only has parent gene 1 is gap-filled (`-`) across the parent-gene-2 columns, and vice versa; a fused taxon has real sequence across both regions with no internal gap block. Sequence names must match exactly between this file, `constrained.nwk`, and `presence_absence.txt`. Trim it first (e.g. trimAl `-gappyout`) — the pipeline does not trim for you.
- **`constrained.nwk`** — a Newick tree with the fused clade already constrained to be monophyletic. Branch lengths are optional: if present, Step 1 is skipped; if absent, Step 1 estimates them from `combined.aln` (or force re-estimation with `-f`).
- **`model.txt`** — a single line with the IQ-TREE model string for `combined.aln` (e.g. `Best-fit model according to BIC: GTR+G4` from an IQ-TREE `.iqtree` log → write `GTR+G4`). Case-sensitive, no rate matrix or frequency values, and amino acid vs. DNA models must match your alignment's data type.
- **`presence_absence.txt`** — a TSV (`species<TAB>gene_type`, header row included) covering every taxon in `combined.aln`, with `gene_type` one of `Parent_gene_1`, `Parent_gene_2`, `fused`, `absent`. This is what Step 5 uses to know which taxa should form the fused clade when checking monophyly.

If you're converting from a typical two-block real analysis (`block1.aln` + `block2.aln`), [docs/INPUTS.md](docs/INPUTS.md) has a ready-to-run Biopython snippet for building `combined.aln` from the two blocks, plus a one-liner for pulling the model string out of an IQ-TREE log.

**Before running the pipeline**, validate the directory:

```bash
python3 src/validate_inputs.py inputs/<your_gene_family_name>
```

This checks the required files exist, the alignment is readable, taxon names match between the alignment and tree, the tree is valid Newick, the model string is valid IQ-TREE notation, and the gap pattern actually shows a mix of fused and parent-only taxa. Fix any reported errors before proceeding — `master_pipeline.sh` also runs this validation automatically (governed by `VALIDATE_INPUTS`/`VALIDATE_STRICT` in `config/pipeline.conf`).

Then point the pipeline at your directory as normal:

```bash
bash scripts/master_pipeline.sh -i inputs/<your_gene_family_name>
```

### 3. Configure the pipeline

Edit [config/pipeline.conf](config/pipeline.conf):

- `N_REPLICATES` — number of simulated replicates per gene family (default `10`; use 100+ for a real analysis)
- `N_BOOTSTRAP` — bootstrap replicates for the AU test (`-zb`, default `1000`)
- `IQTREE_BIN` / `SEQGEN_BIN` / `PYTHON_BIN` — paths to external tools, if not on `$PATH`
- `REPLICATE_JOBS` — how many replicates to process concurrently in Steps 4/5 (defaults to `HPC_THREADS`, which is auto-detected from SLURM/SGE or falls back to `LOCAL_THREADS`)
- `IQTREE_THREADS_PER_JOB` — threads per individual IQ-TREE call (`-T`); leave at `1` unless individual alignments are large (100s of taxa), since parallelizing across replicates is far more efficient than multithreading each small run
- `DEFAULT_OUTPUT_DIR` — where results are written (default `outputs/`)

## Usage

### Run one gene family locally

```bash
bash scripts/master_pipeline.sh -i examples
```

```bash
# With a non-default config file
bash scripts/master_pipeline.sh -i examples -c config/pipeline.conf

# Force re-run of every step, including branch length re-estimation
bash scripts/master_pipeline.sh -i examples -f
```

Results land in `outputs/<input_dir_name>/`.

### Run many gene families

Run `master_pipeline.sh` once per input directory (e.g. loop over `inputs/RHG_*`, or submit as a job array — see [hpc_templates/submit_pipeline.slurm](hpc_templates/submit_pipeline.slurm) for a SLURM template).

Once every gene family has been processed, combine their per-gene `results.tsv` files into one table for cross-gene analysis:

```bash
bash scripts/combine_results.sh -o outputs
```

This writes `outputs/combined_results.tsv` (or pass `-out <path>` to choose another location) with one row per gene family — including `NA`-filled rows for gene families where the pipeline produced no usable replicates, so nothing is silently dropped.

### Quick test

```bash
# Edit config/pipeline.conf and set N_REPLICATES=1, then:
bash scripts/master_pipeline.sh -i examples
cat outputs/examples/step6_summarize/summary.txt
```

## Pipeline Overview

Each gene family is processed through 6 sequential steps:

1. **Estimate branch lengths** (`01_estimate_branch_lengths.sh`) — re-estimate on the constrained tree with IQ-TREE, unless the input tree already has them
2. **Simulate sequences** (`02_run_seqgen.sh`) — compute a majority-rule consensus ancestor and simulate `N_REPLICATES` full-length replicates with SeqGen
3. **Apply gap pattern** (`03_mask_alignments.sh`) — reinsert the original alignment's gap pattern into each simulated replicate
4. **Reconstruct trees** (`04_reconstruct_trees.sh`) — build unconstrained (ML) and constrained (monophyletic) trees per replicate with IQ-TREE, run concurrently across replicates
5. **AU test + monophyly check** (`05_run_au_test.sh`) — compare the two trees with IQ-TREE's AU test and check whether the fused clade is monophyletic in the unconstrained tree, run concurrently across replicates
6. **Summarize** (`06_summarise_results.sh`) — aggregate all replicates into a per-gene-family `results.tsv` and `summary.txt`

See [CLAUDE.md](CLAUDE.md) for implementation-level detail on each step and known bugfixes.

## Output Structure

```
outputs/<gene_family_name>/
├── logs/
│   └── master_pipeline.log
├── step1_branch_lengths/
├── step2_seqgen/
├── step3_masking/
├── step4_tree_reconstruction/
│   ├── step4_results.txt        # per-replicate OK/FAILED status
│   └── rep_*/
│       ├── masked_unconstrained.treefile
│       └── masked_constrained.treefile
├── step5_au_test/
│   ├── step5_results.txt        # per-replicate RESULT/SKIP lines (AU p-value, monophyly status, ...)
│   └── rep_*/
│       ├── au_test_result.iqtree
│       └── monophyly_check.txt
└── step6_summarize/
    ├── results.tsv               # one row per gene family, 14 columns
    └── summary.txt               # human-readable interpretation
```

`outputs/combined_results.tsv` (from `combine_results.sh`) concatenates every gene family's `results.tsv` into a single table.

## Key Output Columns (`results.tsv`)

- `N_Replicates` / `N_Skipped_Replicates` — how many replicates produced usable results
- `N_Mono_Violations` / `Pct_Mono_Violations` — replicates where the fused clade lost monophyly in the unconstrained tree
- `N_Significant_Violations` / `Pct_Significant_Violations` — of those, how many the AU test also rejected (p < 0.05) — the actual false positive signal
- `N_AU_Rejects_Total` / `Pct_AU_Rejects_Total` — all AU rejections, including ones where the clade stayed monophyletic despite a topology change
- `Mean_AU_Pvalue`, `Mean_Tree_Length`, `SD_Tree_Length` — diagnostics for correlating false positive rate with alignment/tree properties

Under a correctly calibrated AU test, `Pct_Significant_Violations` should be ≈5% (α = 0.05).

## Documentation

- [docs/INPUTS.md](docs/INPUTS.md) — input file format specification
- [docs/INSTALLATION.md](docs/INSTALLATION.md) — dependency installation
- [CLAUDE.md](CLAUDE.md) — pipeline internals, known issues, and bugfix history

## Citation

This simulation framework builds on:

- Siu-Ting et al. (2019) Inadvertent paralog inclusion drives artifactual topologies in phylogenomics. *Mol Biol Evol* 36:1344–1356
- Rambaut & Grassly (1997) Seq-Gen: an application for the Monte Carlo simulation of DNA sequence evolution along phylogenetic trees. *Comput Appl Biosci* 13:235–238
- Shimodaira (2002) An approximately unbiased test of phylogenetic tree selection. *Syst Biol* 51:492–508
