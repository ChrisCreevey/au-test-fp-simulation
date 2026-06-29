# Input Data Specification

## Directory Structure

Place all input data for each RHG in its own directory under `inputs/`:

```
inputs/
├── RHG_001/
│   ├── block1.aln
│   ├── block2.aln
│   ├── constrained.nwk
│   ├── model_block1.txt
│   ├── model_block2.txt
│   └── fusion_taxa.txt
├── RHG_002/
│   ├── block1.aln
│   ├── block2.aln
│   ├── constrained.nwk
│   ├── model_block1.txt
│   ├── model_block2.txt
│   └── fusion_taxa.txt
├── RHG_003/
│   └── [same structure]
└── template/
    └── [example files with descriptions]
```

## Required Files (Per RHG)

### 1. `block1.aln` and `block2.aln`

**Content:** Trimmed DNA sequence alignments (one per component homologous block)

**Format:** FASTA or Phylip (IQTree auto-detects)

**Example (FASTA):**
```
>species_A
ATGCGATCG...
>species_B
ATGCGATCG...
>species_C_fused
ATGCGATCG...
```

**Notes:**
- Must include all taxa (fusion and non-fusion species)
- Fusion clade taxa have sequences for both blocks
- Alignment should already be trimmed (using trimAl `-gappyout` as in real analysis)
- Alignment length is extracted automatically — used for SeqGen `-l` parameter
- Sequence names must match exactly between blocks and the tree/fusion_taxa file

### 2. `constrained.nwk`

**Content:** Monophyletic species-level tree in Newick format

**Format:** Newick with branch lengths optional but recommended

**Example:**
```
(species_A:0.01,species_B:0.02,(species_C_fused:0.003,species_D_fused:0.005):0.015);
```

**Notes:**
- This is the **constrained** tree that enforces monophyly of the fusion clade
- Taxon names must match exactly those in the alignments
- Branch lengths will be re-estimated in Step 1
- Can be derived from species tree with fusion clade constrained to be monophyletic

### 3. `model_block1.txt` and `model_block2.txt`

**Content:** Best-fit substitution model for each block (from ModelFinder)

**Format:** Plain text, one model string per file

**Example:**
```
GTR+G4
```

**Other valid examples:**
```
HKY+G4
JC
TVM+G4
K80+G
```

**Notes:**
- Model is extracted from IQTree log file (`.iqtree` output)
- Look for line: `Best-fit model according to BIC: GTR+G4`
- **Do not include** rate matrix details or base frequencies on separate lines
- Format must match IQTree notation exactly (case-sensitive: `GTR`, not `gtr`)
- If model includes gamma rate variation (`+G`), the alpha value will be extracted from Step 1 output
- If model includes invariant sites (`+I`), this must be preserved in the model string

### 4. `fusion_taxa.txt`

**Content:** List of taxa that carry the fused gene (one per line)

**Format:** Plain text, one taxon name per line

**Example:**
```
species_C_fused
species_D_fused
species_E_fused
```

**Notes:**
- Taxon names must match exactly those in the alignments and tree
- No comments or extra whitespace
- This defines which taxa will be **masked** in Step 3 (concatenated blocks)
- All other taxa in the alignment are treated as non-fusion (separate blocks)

## Preparing Input Data from Real Analysis

### From MAFFT and trimAl Output
```bash
# Alignments are already produced by your real pipeline
# Just copy them to the input directory with standard names:
cp my_analysis/RHG_001/block1.trimmed.aln inputs/RHG_001/block1.aln
cp my_analysis/RHG_001/block2.trimmed.aln inputs/RHG_001/block2.aln
```

### From IQTree Log Files
```bash
# Extract best-fit model from IQTree output:
grep "Best-fit model according to BIC" my_analysis/RHG_001/block1.iqtree | awk '{print $NF}' > inputs/RHG_001/model_block1.txt
grep "Best-fit model according to BIC" my_analysis/RHG_001/block2.iqtree | awk '{print $NF}' > inputs/RHG_001/model_block2.txt
```

### From Your Phylostratigraphic Output
```bash
# Constrained tree should enforce monophyly of fusion clade
# If you have a species tree and fusion clade list:
# Use a phylogenetics tool (e.g., ETE Toolkit, APE) to create constrained tree

python3 scripts/create_constrained_tree.py \
  --species_tree species_tree.nwk \
  --fusion_taxa RHG_001_fusion_taxa.txt \
  --output inputs/RHG_001/constrained.nwk
```

### From RHG Presence/Absence Matrix
```bash
# Extract taxa with the fused gene for RHG_001:
python3 scripts/extract_fusion_taxa.py \
  --rhg_presence_matrix rhg_matrix.tsv \
  --rhg_id RHG_001 \
  --output inputs/RHG_001/fusion_taxa.txt
```

## Validation

A validation script checks all inputs before starting the pipeline:

```bash
python3 src/validate_inputs.py inputs/RHG_001
```

This checks:
- ✓ All required files exist
- ✓ Alignment format is readable
- ✓ Taxon names match between files
- ✓ Tree is valid Newick format
- ✓ Model string is valid IQTree notation
- ✓ Fusion taxa exist in alignments
- ✓ At least one fusion and one non-fusion taxon

## Template

A template RHG with dummy data is provided in `inputs/template/` for reference. Copy and modify it:

```bash
cp -r inputs/template inputs/RHG_001
# Then edit the files with your actual data
```
