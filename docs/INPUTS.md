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
│   └── presence_absence.txt
├── RHG_002/
│   ├── block1.aln
│   ├── block2.aln
│   ├── constrained.nwk
│   ├── model_block1.txt
│   ├── model_block2.txt
│   └── presence_absence.txt
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

### 4. `presence_absence.txt`

**Content:** Gene presence/absence matrix for each taxon (which copy of the gene they carry)

**Format:** Tab-separated, two columns: `species` and `gene_type`

**Example:**
```
species	gene_type
species_A	Parent_gene_1
species_B	Parent_gene_2
species_C_fused	fused
species_D_fused	fused
species_E	Parent_gene_1
species_F	absent
```

**Valid gene_type values:**
- `Parent_gene_1` — Taxon carries only the first component gene
- `Parent_gene_2` — Taxon carries only the second component gene
- `fused` — Taxon carries both components concatenated (fusion)
- `absent` — Taxon is absent from this gene family

**Notes:**
- Taxon names must match exactly those in the alignments and tree
- First line is a header (species, gene_type)
- Tab-separated (not spaces)
- This determines how each taxon is represented in Step 3 (masking)
- All taxa in the tree/alignments should be represented in this file

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
# Create presence/absence file from your RHG matrix
# If you have a matrix showing which parent genes and fusions each species has:
python3 << 'EOF'
import pandas as pd

# Load your RHG presence/absence matrix
# Expected format: rows=species, columns=RHGs, values=Parent_gene_1/Parent_gene_2/fused/absent
rhg_matrix = pd.read_csv('your_rhg_matrix.tsv', sep='\t', index_col=0)

# Extract gene types for RHG_001
rhg_001 = rhg_matrix['RHG_001']

# Write to presence/absence file
with open('inputs/RHG_001/presence_absence.txt', 'w') as f:
    f.write('species\tgene_type\n')
    for species, gene_type in rhg_001.items():
        f.write(f'{species}\t{gene_type}\n')
EOF
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
