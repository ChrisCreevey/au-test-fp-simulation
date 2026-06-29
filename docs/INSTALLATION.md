# Installation Guide

Complete guide to installing dependencies for the AU test simulation pipeline.

## Quick Start

**macOS with Homebrew:**
```bash
bash install_dependencies_local.sh
```

**HPC with Conda:**
```bash
bash install_dependencies_hpc.sh
```

---

## Option 1: macOS/Linux with Homebrew (Local)

Best for: Testing locally before running on HPC

### Prerequisites
- Homebrew installed: https://brew.sh

### Installation

```bash
bash install_dependencies_local.sh
```

This installs:
- **IQTree2** — phylogenetic inference
- **SeqGen** — DNA sequence simulation
- **MAFFT** — sequence alignment (optional)
- **GNU Parallel** — parallelization utility (optional)
- **Python 3.11**

### Verification

```bash
iqtree2 --version
seq-gen -v
python3 --version
```

### Usage

Pipeline automatically detects tools in PATH. If tools are elsewhere, edit `config/pipeline.conf`:

```bash
IQTREE_BIN=/usr/local/bin/iqtree2
SEQGEN_BIN=/usr/local/bin/seq-gen
```

---

## Option 2: HPC with Conda (Recommended for Clusters)

Best for: Slurm/PBS/SGE HPC systems, reproducible environments

### Prerequisites

Conda must be available on your HPC. Check:

```bash
module avail anaconda  # or conda, miniconda, etc.
module load anaconda3
conda --version
```

If conda is not available, contact your HPC administrator or install Miniconda locally.

### Installation

```bash
bash install_dependencies_hpc.sh
```

This script:
1. Detects conda (or loads it from modules)
2. Creates conda environment: `au_test_sim`
3. Installs all dependencies via bioconda
4. Generates activation script: `activate_hpc_env.sh`
5. Generates SLURM wrapper: `run_on_hpc.slurm`

### Verification

```bash
conda activate au_test_sim
iqtree2 --version
seq-gen -v
python3 --version
```

### Using the Environment

**In interactive HPC sessions:**
```bash
source activate_hpc_env.sh
bash scripts/master_pipeline.sh -i inputs/RHG_001 -c config/pipeline.conf
```

**In batch SLURM jobs:**
```bash
#!/bin/bash
#SBATCH --job-name=au_test_sim
#SBATCH --time=02:00:00

# Activate conda
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate au_test_sim

# Run pipeline
bash scripts/master_pipeline.sh -i inputs/RHG_001 -c config/pipeline.conf
```

Or use the provided wrapper:
```bash
sbatch run_on_hpc.slurm
```

### Sharing Environments

To export and share your environment:

```bash
# Export environment
conda env export -n au_test_sim > au_test_sim.yml

# On another system, recreate it
conda env create -f au_test_sim.yml
```

---

## Option 3: Conda from environment.yml (All Systems)

Best for: Reproducible, portable setup on any system with conda

Works on: macOS, Linux, HPC, Windows (with WSL)

### Installation

```bash
# Using provided environment file
conda env create -f environment.yml

# Or name it differently
conda env create -f environment.yml -n my_au_test_sim

# Activate
conda activate au_test_sim
```

### Updating Environment

```bash
# Update all packages
conda env update -n au_test_sim -f environment.yml --prune
```

---

## Option 4: Manual Installation (Advanced)

If your HPC doesn't have conda, install tools individually:

### IQTree2

```bash
# Homebrew
brew install iqtree2

# From source
wget https://github.com/iqtree/iqtree2/releases/download/v2.2.5/iqtree-2.2.5-Linux.tar.gz
tar xzf iqtree-2.2.5-Linux.tar.gz
export PATH="$PWD/iqtree-2.2.5-Linux/bin:$PATH"
```

### SeqGen

```bash
# Homebrew
brew install seq-gen

# From source
wget https://tree.bio.ed.ac.uk/download/SeqGen/seq-gen-1.3.4.tar.gz
tar xzf seq-gen-1.3.4.tar.gz
cd Seq-Gen-1.3.4/source
make
export PATH="$PWD:$PATH"
```

### MAFFT

```bash
# Homebrew
brew install mafft

# From source
wget https://mafft.cbrc.jp/alignment/software/mafft-7.520-without-extensions-src.tgz
tar xzf mafft-7.520-without-extensions-src.tgz
cd mafft-7.520-without-extensions/core
make
```

### Python 3

```bash
# Homebrew
brew install python@3.11

# Or use system Python (if >= 3.8)
python3 --version
```

Then edit `config/pipeline.conf` to point to your tool locations.

---

## Troubleshooting

### "iqtree2: command not found"

**On macOS:**
```bash
# Check if installed
brew list iqtree2

# If not, install
brew install iqtree2

# Add to PATH if needed
export PATH="/usr/local/bin:$PATH"
```

**On HPC:**
```bash
# Check available modules
module avail iqtree

# Load module if available
module load iqtree2

# Or activate conda environment
conda activate au_test_sim
```

### "seq-gen: command not found"

```bash
# Check installation
which seq-gen
brew list seq-gen

# If missing
brew install seq-gen

# Or in conda
conda activate au_test_sim
```

### Conda environment not found

```bash
# List all environments
conda env list

# If au_test_sim is missing, recreate it
bash install_dependencies_hpc.sh
```

### Python version too old

```bash
# Check version
python3 --version

# Need Python 3.8+
# On macOS:
brew install python@3.11

# On HPC, use conda environment (it includes Python 3.11)
bash install_dependencies_hpc.sh
```

### Module load fails on HPC

```bash
# List available modules
module avail

# Look for: anaconda, conda, miniconda, python, etc.
module load anaconda3

# Or ask administrator
# "Which modules provide conda?"
```

### Conda installation times out

On slow HPC networks, increase timeout:
```bash
conda config --set remote_connect_timeout_secs 60
bash install_dependencies_hpc.sh
```

Or manually create smaller environment:
```bash
conda create -n au_test_sim -y python=3.11 -c bioconda iqtree2 seq-gen
```

---

## Verifying Installation

Run validation script on any platform:

```bash
python3 src/validate_inputs.py inputs/template
```

Should output:
```
✓ block1.aln exists
✓ block2.aln exists
✓ constrained.nwk: 5 taxa
✓ model_block1.txt: GTR+G4
✓ model_block2.txt: HKY+G4
✓ fusion_taxa.txt: 2 taxa
✓ All fusion taxa found in alignments
✓ Non-fusion taxa count: 3

VALIDATION PASSED - All checks successful ✓
```

---

## Summary Table

| Platform | Method | Speed | Reproducibility | Notes |
|----------|--------|-------|-----------------|-------|
| macOS | Homebrew | Fast | Medium | Best for quick local testing |
| HPC | Conda | Slow (network) | High | Recommended for production |
| Any | environment.yml | Medium | Very High | Most portable option |
| Any | Manual | Medium | Low | For special HPC constraints |

---

## Recommended Setup

**For iterative development:**
1. Install locally with Homebrew (fast feedback)
2. Test on sample data: `bash scripts/master_pipeline.sh -i inputs/RHG_001`

**For production runs:**
1. Use HPC conda environment: `bash install_dependencies_hpc.sh`
2. Submit job array: `sbatch hpc_templates/submit_pipeline.slurm`
3. Monitor with: `squeue -j $JOBID`

---

For more help, see:
- `README.md` — Full overview
- `QUICKSTART.md` — Quick start guide
- `config/pipeline.conf` — Configuration parameters
