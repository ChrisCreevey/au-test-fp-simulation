#!/bin/bash
# Install dependencies for AU test simulation pipeline on macOS with Homebrew

set -euo pipefail

echo "========================================================================"
echo "Installing Dependencies for AU Test Simulation Pipeline"
echo "========================================================================"
echo ""
echo "System: macOS (Homebrew)"
echo "Target: Local testing and single-node runs"
echo ""

# ============================================================================
# Check for Homebrew
# ============================================================================

if ! command -v brew &> /dev/null; then
    echo "ERROR: Homebrew not found. Install from https://brew.sh"
    exit 1
fi

echo "✓ Homebrew is installed"
echo ""

# ============================================================================
# Update Homebrew
# ============================================================================

echo "Updating Homebrew packages..."
brew update
echo "✓ Homebrew updated"
echo ""

# ============================================================================
# Core Dependencies
# ============================================================================

echo "Installing core dependencies..."

# IQTree 2 (phylogenetic inference)
if ! command -v iqtree2 &> /dev/null; then
    echo "  Installing IQTree2..."
    brew install iqtree2
    echo "  ✓ IQTree2 installed"
else
    echo "  ✓ IQTree2 already installed"
fi

echo ""

# ============================================================================
# Python and Python Dependencies
# ============================================================================

echo "Installing Python 3..."

if ! command -v python3 &> /dev/null; then
    echo "  Installing Python 3..."
    brew install python@3.11
    echo "  ✓ Python 3 installed"
else
    echo "  ✓ Python 3 already installed"
fi

PYTHON_BIN=$(command -v python3)
echo "  Using Python: $PYTHON_BIN"

# Upgrade pip
echo "  Upgrading pip..."
"$PYTHON_BIN" -m pip install --upgrade pip
echo "  ✓ pip upgraded"

echo ""

# ============================================================================
# Verify Installations
# ============================================================================

echo "Verifying installations..."
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
check_tool "Python3" "python3" || ALL_OK=false

echo ""

# ============================================================================
# Installation Summary
# ============================================================================

echo "========================================================================"

if [[ "$ALL_OK" == "true" ]]; then
    echo "SUCCESS: All dependencies installed and verified ✓"
    echo ""
    echo "Next steps:"
    echo "  1. Build a manifest against your real RHG data (see docs/INPUTS.md):"
    echo "     python3 src/build_manifest.py --data-dir <dir> \\"
    echo "         --output-json data/manifest/manifest_full.json \\"
    echo "         --output-tsv  data/manifest/manifest_full.tsv"
    echo "  2. Run one family:"
    echo "     bash scripts/run_family.sh -f <family_id> \\"
    echo "         -m data/manifest/manifest_full.json -o data/sim_outputs \\"
    echo "         -c config/pipeline.conf"
    echo ""
    echo "Tool paths (edit config/pipeline.conf if different):"
    echo "  IQTREE_BIN=$(command -v iqtree2)"
    echo "  PYTHON_BIN=$PYTHON_BIN"
else
    echo "ERROR: Some dependencies failed to install"
    exit 1
fi

echo "========================================================================"
