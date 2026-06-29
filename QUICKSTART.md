# Quick Start Guide

Get the AU test simulation pipeline running in 5 minutes.

## 1. Prepare Input Data

Place your RHG data in the `inputs/` directory following this structure:

```
inputs/
├── RHG_001/
│   ├── block1.aln              # Alignment of component block 1 (FASTA/Phylip)
│   ├── block2.aln              # Alignment of component block 2
│   ├── constrained.nwk         # Monophyletic species tree (Newick)
│   ├── model_block1.txt        # Best-fit model (e.g., "GTR+G4")
│   ├── model_block2.txt        # Best-fit model (e.g., "HKY+G4")
│   └── fusion_taxa.txt         # List of fusion clade taxa (one per line)
├── RHG_002/
│   └── [same 6 files]
└── template/                   # [Reference template — keep as-is]
```

**See `inputs/template/` and `docs/INPUTS.md` for detailed format specifications.**

### Quick Data Prep from Real Analysis

If you already have data from the real AU test pipeline:

```bash
# Copy alignments
cp /path/to/analysis/RHG_001/block1.trimmed.aln inputs/RHG_001/block1.aln
cp /path/to/analysis/RHG_001/block2.trimmed.aln inputs/RHG_001/block2.aln

# Extract models from IQTree logs
grep "Best-fit model" /path/to/RHG_001/block1.iqtree | awk '{print $NF}' > inputs/RHG_001/model_block1.txt
grep "Best-fit model" /path/to/RHG_001/block2.iqtree | awk '{print $NF}' > inputs/RHG_001/model_block2.txt

# Use species tree with monophyletic fusion clade constraint
cp /path/to/species_tree.nwk inputs/RHG_001/constrained.nwk

# List taxa carrying the fused gene
echo "species_C_fused" > inputs/RHG_001/fusion_taxa.txt
echo "species_D_fused" >> inputs/RHG_001/fusion_taxa.txt
```

## 2. Validate Inputs (Optional but Recommended)

Before running the pipeline, check that all inputs are correct:

```bash
python3 src/validate_inputs.py inputs/RHG_001
```

You should see:
```
✓ block1.aln exists
✓ block2.aln exists
✓ constrained.nwk: N taxa
✓ model_block1.txt: GTR+G4
✓ model_block2.txt: HKY+G4
✓ fusion_taxa.txt: 2 taxa
✓ All fusion taxa found in alignments
✓ Non-fusion taxa count: 3

VALIDATION PASSED - All checks successful ✓
```

## 3. Run on a Single RHG (Test)

```bash
bash scripts/master_pipeline.sh -i inputs/RHG_001 -c config/pipeline.conf
```

This runs all 6 steps and produces:
- `outputs/RHG_001/step1_branch_lengths/` — Branch-length trees
- `outputs/RHG_001/step2_seqgen/` — Simulated sequences (100 replicates)
- `outputs/RHG_001/step3_masking/` — Masked alignments
- `outputs/RHG_001/step4_tree_reconstruction/` — IQTree results for each replicate
- `outputs/RHG_001/step5_au_test/` — AU test results
- `outputs/RHG_001/results.tsv` — **Final false positive rate summary**

Expected runtime: 10–60 minutes per RHG (depending on alignment size and compute resources)

## 4. Run All RHGs on HPC (SLURM)

Edit the job array count in the SLURM submission script:

```bash
# Check how many RHGs you have
ls -d inputs/RHG_* | wc -l
# Output: 100

# Edit submit_pipeline.slurm, change:
# #SBATCH --array=1-100
# to match your RHG count
```

Then submit:

```bash
sbatch hpc_templates/submit_pipeline.slurm
```

This submits one job per RHG, processing them in parallel across the cluster.

## 5. Collect Results

After all RHGs are complete:

```bash
python3 src/summarise_results.py \
    --step5-output outputs/RHG_001/step5_au_test \
    --output-dir outputs/RHG_001 \
    --rhg-id RHG_001
```

This produces:
- `outputs/RHG_001/results.tsv` — Per-RHG false positive rates
- `outputs/RHG_001/summary.txt` — Interpretation and statistics

For the **null distribution** across all families, you would aggregate results like:

```bash
cat outputs/*/results.tsv > all_results.tsv
# Plot distribution of false positive rates
```

## 6. Interpret Results

Your key output:
- **False positive rate per RHG** — % of 100 replicates where AU test incorrectly rejected monophyly
- **Null distribution** — Mean, median, 95th percentile across all RHGs

**Expected behavior:** If the AU test is correct, ~5% false positive rate (p < 0.05).

**Real data interpretation:**
- Observed 41% polyphyletic in real data
- If null FP rate << 5%, the 41% is likely true signal
- If null FP rate >> 5%, the AU test may be oversensitive

## Configuration

Edit `config/pipeline.conf` to adjust:
- `N_REPLICATES` — replicates per RHG (default: 100)
- `N_BOOTSTRAP` — bootstrap replicates for tree inference (default: 1000)
- Tool paths (IQTREE_BIN, SEQGEN_BIN, etc.)
- HPC settings (memory, CPUs, queue)

## Troubleshooting

**IQTree not found:**
```bash
# Check installation
which iqtree2
# or
conda list | grep iqtree
```

**SeqGen not found:**
```bash
which seq-gen
# or install:
conda install -c bioconda seq-gen
```

**Python import errors:**
```bash
# Check Python version
python3 --version  # Should be 3.8+
```

**See logs for detailed errors:**
```bash
cat outputs/RHG_001/logs/master_pipeline.log
cat outputs/RHG_001/step1_*/step1.log
```

## Full Documentation

- `README.md` — Complete overview
- `docs/INPUTS.md` — Input file format and preparation
- `docs/PIPELINE.md` — Detailed description of each step (coming soon)
- `config/pipeline.conf` — All configurable parameters
- Shell scripts (`scripts/*.sh`) — Inline comments explain each step

---

Questions? See `docs/TROUBLESHOOTING.md` or review the shell script logs.
