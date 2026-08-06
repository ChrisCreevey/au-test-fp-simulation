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

### 2. Configure the pipeline

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
