#!/bin/bash
# Install dependencies for AU test simulation pipeline on HPC systems
# Supports: conda-based and module-based systems (Slurm, PBS, SGE)

set -euo pipefail

echo "========================================================================"
echo "Installing Dependencies for AU Test Simulation Pipeline"
echo "========================================================================"
echo ""
echo "System: HPC Cluster (Slurm/PBS/SGE)"
echo ""

# ============================================================================
# Detect HPC Environment
# ============================================================================

detect_hpc_system() {
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        echo "SLURM"
    elif [[ -n "${PBS_JOBID:-}" ]]; then
        echo "PBS"
    elif [[ -n "${SGE_TASK_ID:-}" ]]; then
        echo "SGE"
    else
        echo "UNKNOWN"
    fi
}

HPC_SYSTEM=$(detect_hpc_system)
echo "Detected HPC system: $HPC_SYSTEM"
echo ""

# ============================================================================
# Check for Conda
# ============================================================================

if command -v conda &> /dev/null; then
    echo "✓ Conda found ($(conda --version))"
    CONDA_BIN=$(command -v conda)
    USE_CONDA=true
else
    echo "⚠ Conda not found in PATH"
    echo ""
    echo "Option 1: Load conda module (if available on your HPC)"
    echo "  module load anaconda3    # or similar"
    echo "  bash install_dependencies_hpc.sh"
    echo ""
    echo "Option 2: Initialize conda from scratch"
    echo "  See: https://conda.io/projects/conda/en/latest/user-guide/install/index.html"
    echo ""
    USE_CONDA=false
fi

if [[ "$USE_CONDA" == "false" ]]; then
    # Try to load conda module
    if [[ -x "$(command -v module)" ]]; then
        echo "Attempting to load conda module..."

        # Try common module names
        for module_name in "anaconda3" "anaconda" "miniconda3" "miniconda"; do
            if module load "$module_name" 2>/dev/null; then
                echo "✓ Loaded module: $module_name"
                USE_CONDA=true
                CONDA_BIN=$(command -v conda)
                break
            fi
        done
    fi

    if [[ "$USE_CONDA" == "false" ]]; then
        echo "ERROR: Could not find or load conda"
        echo ""
        echo "Please either:"
        echo "  1. Load conda/anaconda module: module load anaconda3"
        echo "  2. Initialize conda in your shell: conda init bash"
        echo "  3. Use system package manager (see Option 2 below)"
        exit 1
    fi
fi

echo ""

# ============================================================================
# Create Conda Environment
# ============================================================================

ENV_NAME="au_test_sim"

echo "Creating conda environment: $ENV_NAME"
echo ""

if conda info --envs | grep -q "^$ENV_NAME "; then
    echo "⚠ Conda environment '$ENV_NAME' already exists"
    echo ""
    read -p "Remove and recreate? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        conda env remove -n "$ENV_NAME" -y
    else
        echo "Using existing environment"
        exit 0
    fi
fi

echo "Installing packages..."
echo ""

# Create environment with bioconda packages
conda create -n "$ENV_NAME" -y \
    -c conda-forge \
    -c bioconda \
    python=3.11 \
    iqtree2 \
    seq-gen \
    mafft \
    parallel::gnu-parallel

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to create conda environment"
    exit 1
fi

echo ""
echo "✓ Conda environment created successfully"

# ============================================================================
# Activate Environment and Verify
# ============================================================================

echo ""
echo "Activating environment..."

# Source conda.sh to enable conda command
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"

conda activate "$ENV_NAME"

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to activate conda environment"
    exit 1
fi

echo "✓ Environment activated"
echo ""

# ============================================================================
# Verify Installations
# ============================================================================

echo "Verifying installations in conda environment..."
echo ""

check_tool() {
    local tool=$1
    local cmd=$2

    if command -v "$cmd" &> /dev/null; then
        echo "✓ $tool"
        "$cmd" --version 2>/dev/null || "$cmd" -v 2>/dev/null || echo "  (version check not available)"
    else
        echo "✗ $tool NOT FOUND"
        return 1
    fi
}

ALL_OK=true

check_tool "IQTree2" "iqtree2" || ALL_OK=false
check_tool "SeqGen" "seq-gen" || ALL_OK=false
check_tool "MAFFT" "mafft" || ALL_OK=false
check_tool "GNU Parallel" "parallel" || ALL_OK=false
check_tool "Python3" "python3" || ALL_OK=false

echo ""
PYTHON_BIN=$(command -v python3)
echo "Python executable: $PYTHON_BIN"
echo ""

# ============================================================================
# Create Environment Activation Script
# ============================================================================

ACTIVATE_SCRIPT="activate_hpc_env.sh"

cat > "$ACTIVATE_SCRIPT" <<'SCRIPT_EOF'
#!/bin/bash
# Activate conda environment for AU test simulation pipeline
# Run this at the start of each HPC job: source activate_hpc_env.sh

set -euo pipefail

ENV_NAME="au_test_sim"

# Source conda initialization
CONDA_BASE=$(conda info --base 2>/dev/null || echo "")

if [[ -z "$CONDA_BASE" ]]; then
    # Try to find conda
    if [[ -x "$(command -v conda)" ]]; then
        CONDA_BASE=$(conda info --base)
    else
        echo "ERROR: Could not find conda. Please load conda module:"
        echo "  module load anaconda3"
        exit 1
    fi
fi

echo "Sourcing conda from: $CONDA_BASE"
source "$CONDA_BASE/etc/profile.d/conda.sh"

echo "Activating environment: $ENV_NAME"
conda activate "$ENV_NAME"

if [[ $? -eq 0 ]]; then
    echo "✓ Environment activated"
    echo ""
    echo "Environment info:"
    conda info
else
    echo "ERROR: Failed to activate environment"
    exit 1
fi
SCRIPT_EOF

chmod +x "$ACTIVATE_SCRIPT"
echo "✓ Created activation script: $ACTIVATE_SCRIPT"

# ============================================================================
# Create SLURM Wrapper Script
# ============================================================================

SLURM_WRAPPER="run_on_hpc.slurm"

cat > "$SLURM_WRAPPER" <<'SLURM_EOF'
#!/bin/bash
#SBATCH --job-name=au_test_sim_test
#SBATCH --array=1
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# SLURM wrapper that automatically activates conda environment and runs pipeline

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "========================================================================"
echo "HPC Job Configuration"
echo "========================================================================"
echo "Project directory: $PROJECT_DIR"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}G (or ${SLURM_MEM_PER_CPU}G per CPU)"
echo ""

# Activate conda environment
echo "Activating conda environment..."
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate au_test_sim

if [[ $? -ne 0 ]]; then
    echo "ERROR: Failed to activate conda environment"
    exit 1
fi

echo "✓ Conda environment activated"
echo ""

# Get RHG to process
cd "$PROJECT_DIR"
RHG_DIRS=($(ls -d inputs/RHG_* 2>/dev/null | sort))

if [[ ${#RHG_DIRS[@]} -eq 0 ]]; then
    echo "ERROR: No RHG directories found in inputs/"
    exit 1
fi

CURRENT_RHG="${RHG_DIRS[$((SLURM_ARRAY_TASK_ID - 1))]}"

echo "Processing: $CURRENT_RHG"
echo ""

# Run pipeline
bash scripts/master_pipeline.sh \
    -i "$CURRENT_RHG" \
    -c config/pipeline.conf

echo ""
echo "========================================================================"
echo "Job completed successfully"
echo "========================================================================"
SLURM_EOF

chmod +x "$SLURM_WRAPPER"
echo "✓ Created SLURM wrapper: $SLURM_WRAPPER"

echo ""

# ============================================================================
# Summary
# ============================================================================

echo "========================================================================"

if [[ "$ALL_OK" == "true" ]]; then
    echo "SUCCESS: HPC environment set up successfully ✓"
    echo ""
    echo "Conda environment name: $ENV_NAME"
    echo "Location: $(conda info --envs | grep "$ENV_NAME" | awk '{print $2}')"
    echo ""
    echo "Next steps:"
    echo ""
    echo "1. In each HPC job, activate the environment:"
    echo "   source activate_hpc_env.sh"
    echo ""
    echo "2. Or use the provided SLURM wrapper (modifies job array):"
    echo "   sbatch run_on_hpc.slurm"
    echo ""
    echo "3. To use in your own job submission script, add:"
    echo "   # At the start of your .slurm script:"
    echo "   CONDA_BASE=\$(conda info --base)"
    echo "   source \$CONDA_BASE/etc/profile.d/conda.sh"
    echo "   conda activate au_test_sim"
    echo ""
    echo "4. Tool paths (for config/pipeline.conf):"
    echo "   IQTREE_BIN=$(command -v iqtree2)"
    echo "   SEQGEN_BIN=$(command -v seq-gen)"
    echo "   MAFFT_BIN=$(command -v mafft)"
    echo "   PYTHON_BIN=$PYTHON_BIN"
else
    echo "ERROR: Some dependencies failed to install"
    exit 1
fi

echo "========================================================================"
