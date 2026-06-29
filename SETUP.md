# Complete Setup Guide

End-to-end instructions to go from cloned repo to running simulations.

## System Requirements

- **RAM:** 4+ GB per RHG (IQTree tree inference)
- **Disk:** 100 MB per RHG minimum (depends on alignment size)
- **Time:** 10–60 min per RHG locally; 2 hours wall-clock on HPC with parallelization
- **OS:** macOS, Linux, or HPC cluster with Slurm/PBS/SGE

---

## Installation: Choose Your Platform

### For macOS (Local Testing)

```bash
# 1. Install dependencies
bash install_dependencies_local.sh

# 2. Verify
python3 --version     # Should be 3.8+
iqtree2 --version     # Should be 2.0+
seq-gen -v            # Should be 1.3+
```

**Time:** ~5 minutes

### For HPC (Production Runs)

```bash
# 1. Load conda (if not already loaded)
module load anaconda3    # Or equivalent on your system

# 2. Install dependencies and create environment
bash install_dependencies_hpc.sh

# 3. Activate environment
source activate_hpc_env.sh

# 4. Verify
python3 --version
iqtree2 --version
seq-gen -v
```

**Time:** ~10 minutes (depends on network speed)

### Troubleshooting Installation

See `docs/INSTALLATION.md` for detailed troubleshooting:
- "iqtree2: command not found"
- "seq-gen not available"
- "Conda environment not found"
- And more...

---

## Data Preparation

### Step 1: Gather Input Files

From your real AU test analysis, collect for each RHG:
- **Trimmed alignments:** `block1.aln`, `block2.aln` (FASTA/Phylip)
- **Constrained tree:** `constrained.nwk` (Newick, monophyletic fusion clade)
- **Best-fit models:** `model_block1.txt`, `model_block2.txt` (e.g., "GTR+G4")
- **Fusion taxa list:** `fusion_taxa.txt` (one taxon per line)

**See `docs/INPUTS.md` for detailed format specifications.**

### Step 2: Create Directory Structure

```bash
# For one RHG
mkdir -p inputs/RHG_001
cp /path/to/analysis/block1.trimmed.aln inputs/RHG_001/block1.aln
cp /path/to/analysis/block2.trimmed.aln inputs/RHG_001/block2.aln
# ... etc

# For multiple RHGs
for rhg in RHG_001 RHG_002 RHG_003; do
    mkdir -p inputs/$rhg
    # Copy files for each
done
```

Or use the template as a starting point:

```bash
cp -r inputs/template inputs/RHG_001
# Edit each file with your actual data
```

### Step 3: Validate Inputs

Before running the full pipeline:

```bash
python3 src/validate_inputs.py inputs/RHG_001
```

Expected output:
```
✓ block1.aln exists
✓ block2.aln exists
✓ constrained.nwk: 5 taxa
✓ model_block1.txt: GTR+G4
✓ All fusion taxa found in alignments
✓ Non-fusion taxa count: 3

VALIDATION PASSED - All checks successful ✓
```

If validation fails, fix the issues (see `docs/INPUTS.md`).

---

## Configuration

Optional: Customize pipeline parameters by editing `config/pipeline.conf`:

```bash
# Number of replicates per RHG (default: 100)
N_REPLICATES=100

# Bootstrap replicates for tree inference (default: 1000)
N_BOOTSTRAP=1000

# Local thread count (default: 4)
LOCAL_THREADS=4

# HPC memory per job (default: 8G)
HPC_MEMORY_PER_JOB=8

# HPC wall time per RHG (default: 2:00:00)
HPC_TIME_PER_RHG="02:00:00"
```

For most cases, defaults are fine.

---

## Running the Pipeline

### Option A: Single RHG Locally (Testing)

```bash
bash scripts/master_pipeline.sh \
    -i inputs/RHG_001 \
    -c config/pipeline.conf
```

**Expected output:**
```
Starting pipeline for RHG: RHG_001
Input: /path/to/inputs/RHG_001
Output: /path/to/outputs/RHG_001

===== STEP 1: branch_lengths =====
...
===== STEP 6: summarize =====
...

Pipeline completed for RHG: RHG_001
Results saved to: outputs/RHG_001
```

Results in: `outputs/RHG_001/results.tsv`

### Option B: All RHGs on HPC (Production)

```bash
# 1. Count RHGs
ls -d inputs/RHG_* | wc -l
# Output: 100

# 2. Edit SLURM template to set job array
# Edit hpc_templates/submit_pipeline.slurm
# Change: #SBATCH --array=1-100

# 3. Submit job array
sbatch hpc_templates/submit_pipeline.slurm

# 4. Monitor progress
squeue -j <JOBID>
```

Each job processes one RHG in parallel.

### Option C: Process Multiple RHGs Locally (Parallel)

Use GNU Parallel to run multiple RHGs on local machine:

```bash
# Process 4 RHGs at a time (adjust -j to match CPU count)
parallel -j 4 \
    'bash scripts/master_pipeline.sh -i inputs/{} -c config/pipeline.conf' \
    ::: $(ls -d inputs/RHG_* | xargs -n1 basename)
```

---

## Monitoring Progress

### Local

Monitor logs in real-time:

```bash
tail -f outputs/RHG_001/logs/master_pipeline.log
```

### HPC

Monitor jobs:

```bash
# Check job status
squeue -u $USER

# Check specific job
squeue -j <JOBID>

# Monitor logs
tail -f outputs/RHG_001/logs/master_pipeline.log
```

---

## Results

### Per-RHG Results

After each RHG completes:

```
outputs/RHG_001/
├── step1_branch_lengths/
│   ├── block1_with_brlen.treefile
│   └── block2_with_brlen.treefile
├── step2_seqgen/
│   ├── simulated_block1.phy
│   └── simulated_block2.phy
├── step3_masking/
│   ├── rep_001_block1.aln
│   ├── rep_001_block2.aln
│   └── ... (100 replicates)
├── step4_tree_reconstruction/
│   └── rep_*/
├── step5_au_test/
│   └── rep_*/
├── results.tsv               ← **Key output**
└── logs/
```

### Key Output: results.tsv

```
RHG_ID      Block    FP_Rate_%  AU_Rejects  Total_Replicates
RHG_001     Block1   6.00       6           100
RHG_001     Block2   5.00       5           100
RHG_001     Combined 5.50       11          200
```

**False Positive Rate (FP):** % of replicates where AU test incorrectly rejected monophyly

- **Expected under null:** ~5% (α = 0.05)
- **If FP rate << 5%:** AU test is conservative (good)
- **If FP rate >> 5%:** AU test is liberal (may have false positives)

### Null Distribution

Aggregate results across all RHGs:

```bash
# Concatenate all results
cat outputs/*/results.tsv > all_results.tsv

# Summary statistics (using Python)
python3 << 'EOF'
import pandas as pd
df = pd.read_csv('all_results.tsv', sep='\t')
print("False Positive Rate Statistics")
print("=" * 50)
print(df['FP_Rate_%'].describe())
print()
print(f"95th percentile: {df['FP_Rate_%'].quantile(0.95):.2f}%")
print(f"Observed rate in real data: 41%")
print()
if df['FP_Rate_%'].quantile(0.95) < 41:
    print("→ 41% substantially exceeds null expectation")
    print("→ Strong evidence for true polyphyly signal")
EOF
```

---

## Interpreting Results

### Example Output

```
False Positive Rate Statistics
count       100.000000
mean         5.430000
std          2.100000
min          1.000000
25%          4.000000
50%          5.000000
75%          6.500000
max         14.000000

95th percentile: 10.50%
Observed rate in real data: 41%
```

### Interpretation

> **"Under a null model of single-origin fusion, the AU test produced a 
> false positive rate of 5.4% (95th percentile: 10.5%) across simulated 
> families. The observed rate of 41% polyphyletic families in the real 
> data substantially exceeds this null expectation, providing strong 
> evidence that these families represent true convergent evolution 
> rather than AU test artifacts."**

---

## Common Workflows

### Quick Test (30 minutes)

```bash
# 1. Install
bash install_dependencies_local.sh

# 2. Prepare one test RHG
cp -r inputs/template inputs/RHG_test
# Edit with your data

# 3. Validate
python3 src/validate_inputs.py inputs/RHG_test

# 4. Run
bash scripts/master_pipeline.sh -i inputs/RHG_test -c config/pipeline.conf

# 5. Check results
cat outputs/RHG_test/results.tsv
```

### Full Production Run (hours on HPC)

```bash
# 1. Install on HPC
module load anaconda3
bash install_dependencies_hpc.sh

# 2. Prepare all RHGs (100+)
# Copy all input directories

# 3. Edit job array count
vim hpc_templates/submit_pipeline.slurm
# Change: #SBATCH --array=1-200

# 4. Submit
sbatch hpc_templates/submit_pipeline.slurm

# 5. Wait for completion
watch squeue -u $USER

# 6. Aggregate results
cat outputs/*/results.tsv > all_results.tsv
```

### Resume After Interruption

```bash
# Pipeline is designed to resume from last successful step
bash scripts/master_pipeline.sh -i inputs/RHG_001 -c config/pipeline.conf

# Force re-run a step
bash scripts/master_pipeline.sh -i inputs/RHG_001 -c config/pipeline.conf -f
```

---

## Troubleshooting

### Pipeline Fails at Step 1

**Error:** "IQTree not found" or "Command failed"

**Solution:**
```bash
# Check tools are in PATH
which iqtree2
which seq-gen

# Edit config/pipeline.conf to set correct paths
IQTREE_BIN=/full/path/to/iqtree2
SEQGEN_BIN=/full/path/to/seq-gen
```

### Pipeline Fails at Step 3 (Masking)

**Error:** "Python script failed" or "Missing sequences"

**Solution:**
```bash
# Validate inputs first
python3 src/validate_inputs.py inputs/RHG_001

# Check fusion taxa match alignment
grep -c "^>" inputs/RHG_001/block1.aln
wc -l inputs/RHG_001/fusion_taxa.txt
```

### Out of Memory

**Error:** "IQTree killed" or "seq-gen: cannot allocate"

**Solution:**
```bash
# On HPC, request more memory
# Edit SLURM job:
#SBATCH --mem=16G

# Or increase alignment chunk size
# (requires splitting large alignments)
```

### Slow Performance

**Error:** Job takes > 2 hours

**Solutions:**
- Reduce N_BOOTSTRAP in config/pipeline.conf (default 1000 is conservative)
- Use smaller test dataset first
- Check alignment size (very large alignments = slow IQTree)

For more help, see `docs/TROUBLESHOOTING.md` (to be created).

---

## Next Steps

1. ✅ Install dependencies (`install_dependencies_local.sh` or `.sh_hpc`)
2. ✅ Prepare input data (see `docs/INPUTS.md`)
3. ✅ Validate inputs (`src/validate_inputs.py`)
4. ✅ Run pipeline (`scripts/master_pipeline.sh`)
5. ✅ Collect results (aggregate `.tsv` files)
6. ✅ Generate null distribution (plot FP rates)
7. ✅ Write manuscript section

---

## Getting Help

- **Installation issues?** → `docs/INSTALLATION.md`
- **Input format questions?** → `docs/INPUTS.md`
- **Pipeline not running?** → Check `outputs/RHG_*/logs/`
- **Need quick overview?** → `QUICKSTART.md`
- **Full details?** → `README.md`

---

Good luck with your simulation! 🧬
