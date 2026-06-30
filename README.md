# AU Test False Positive Rate Simulation Pipeline

A computational pipeline to test whether the AU (Approximately Unbiased) test produces false positives in detecting convergently evolved remodelled (fused) gene families. 

This pipeline simulates single-origin gene fusions with empirically calibrated parameters and measures how often the AU test incorrectly rejects monophyly, generating a null distribution against which observed polyphyly rates can be compared.

## Quick Start

### 1. Set Up Input Data

For each RHG, create a folder in `inputs/` with these **4 required files**:

| File | Format | Description |
|------|--------|-------------|
| `combined.aln` | FASTA/Phylip | Combined alignment (first 50% = parent_gene_1, second 50% = parent_gene_2) |
| `constrained.nwk` | Newick | Monophyletic tree with fused clade constraint |
| `model.txt` | Text | Substitution model (e.g., `GTR+G4` or `HKY+G4`) |
| `presence_absence.txt` | TSV | Gene type per taxon (Parent_gene_1, Parent_gene_2, fused, absent) |

**Example structure:**
```
inputs/
├── RHG_001/
│   ├── combined.aln
│   ├── constrained.nwk
│   ├── model.txt
│   └── presence_absence.txt
├── RHG_002/
│   ├── combined.aln
│   ├── constrained.nwk
│   ├── model.txt
│   └── presence_absence.txt
```

See `docs/INPUTS.md` for detailed format specifications.

### 2. Configure Pipeline

Edit `config/pipeline.conf` to set:
- Number of replicates per RHG
- Paths to external tools (iqtree, seq-gen)
- Number of parallel jobs (for HPC)
- Output directories

### 3. Run Locally (Testing)

```bash
# Test on a single RHG
bash scripts/master_pipeline.sh -i inputs/RHG_001

# With custom config file (optional)
bash scripts/master_pipeline.sh -i inputs/RHG_001 -c config/pipeline.conf

# Force re-run (recalculate branch lengths, redo all steps)
bash scripts/master_pipeline.sh -i inputs/RHG_001 -f
```

### 4. Run on HPC (Cluster)

```bash
# Submit as SLURM job array (processes all RHGs in parallel)
sbatch hpc_templates/submit_pipeline.slurm

# Or SGE / other schedulers — see hpc_templates/
```

## Pipeline Overview

The pipeline consists of 6 sequential steps per RHG:

1. **Estimate branch lengths** — Re-estimate on constrained tree with IQTree (uses combined.aln)
2. **Calculate consensus** — Majority-rule consensus from combined.aln (ancestral sequence)
3. **Simulate sequences** — Generate 100 replicates with SeqGen from consensus ancestor
4. **Apply gap pattern** — Insert original gap pattern from combined.aln into simulated sequences
5. **Reconstruct gene trees** — IQTree unconstrained and constrained per replicate
6. **Run AU test + monophyly check** — Test constraint and verify fused clade monophyly
7. **Summarize results** — Compute false positive rate, categorize by monophyly status

See `docs/PIPELINE.md` for detailed step-by-step descriptions.

## Dependencies

- **IQTree 2.0+** — phylogenetic tree inference and AU testing
- **SeqGen** — sequence simulation
- **MAFFT** (optional) — only if re-aligning simulated sequences
- **Python 3.8+** — for data processing and summarization
- **GNU Parallel** or SLURM — for parallelization

Install with conda or system package manager:
```bash
conda install -c bioconda iqtree seq-gen
```

## Output Structure

```
outputs/
├── RHG_001/
│   ├── step1_branch_lengths/
│   │   ├── block1_with_brlen.treefile
│   │   └── block2_with_brlen.treefile
│   ├── step2_seqgen/
│   │   ├── block1_simulated.phy (100 replicates)
│   │   └── block2_simulated.phy (100 replicates)
│   ├── step3_masked/
│   │   ├── rep_001_block1.aln
│   │   ├── rep_001_block2.aln
│   │   └── ...
│   ├── step4_trees/
│   │   └── rep_*/
│   ├── step5_au_test/
│   │   └── rep_*/
│   └── results.tsv (false positive rate per RHG)
├── RHG_002/
│   └── ...
└── SUMMARY.txt (null distribution across all RHGs)
```

## Key Output

For each RHG:
- **False positive rate** — % of 100 replicates where AU test rejected monophyly
- **Bootstrap support** — mean support at fusion clade node

Across all RHGs:
- **Null distribution** — mean, median, 95th percentile of false positive rates
- **Comparison to observed 41%** — interpretation for manuscript

## Configuration and Troubleshooting

- See `docs/CONFIG.md` for all configurable parameters
- See `docs/TROUBLESHOOTING.md` for common issues
- Logs are saved to `outputs/<RHG_id>/logs/`

## Citation

This simulation framework is based on:
- Siu-Ting et al. (2019) Inadvertent paralog inclusion drives artifactual topologies in phylogenomics. *Mol Biol Evol* 36:1344–1356
- SeqGen: Rambaut & Grassly (1997) Comput Appl Biosci 13:235–238
- IQTree AU test: Shimodaira (2002) Syst Biol 51:492–508

---

For detailed documentation, see the `docs/` directory.
