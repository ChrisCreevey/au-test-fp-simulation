# Input Data Specification

## Directory Structure

Place all input data for each RHG in its own directory under `inputs/`:

```
inputs/
├── RHG_001/
│   ├── combined.aln
│   ├── constrained.nwk
│   └── model.txt
├── RHG_002/
│   ├── combined.aln
│   ├── constrained.nwk
│   └── model.txt
├── RHG_003/
│   └── [same structure]
└── template/
    └── [example files with descriptions]
```

## Required Files (Per RHG)

### 1. `combined.aln`

**Content:** Combined DNA sequence alignment with both parent gene regions

**Format:** FASTA or Phylip (IQTree auto-detects)

**Structure:**
- Single alignment containing the full sequence space of both parent genes
- Gap characters (`-`) indicate missing regions:
  - Parent_gene_1 only: gaps where parent_gene_2 region would be
  - Parent_gene_2 only: gaps where parent_gene_1 region would be
  - Fused: sequences for both regions (no gaps between)

**Example (FASTA):**
```
>species_A
ATGCGATCG----------CGTAGGCTAGGCT
>species_B
----------CGTAGGCTAGGCTAGGCTAGGC
>species_C_fused
ATGCGATCGCGTAGGCTAGGCTAGGCTAGGC
>species_D_fused
ATGCGATCGCGTAGGCTAGGCTAGGCTAGGC
>species_E
ATGCGATCG----------CGTAGGCTAGGCT
```

**Notes:**
- All taxa must be present (gaps indicate missing regions)
- Alignment should already be trimmed (using trimAl `-gappyout` as in real analysis)
- Gap positions automatically encode which genes each taxon carries
- Sequence names must match exactly those in the constrained tree
- The alignment length includes all gap positions

### 2. `constrained.nwk`

**Content:** Monophyletic species-level tree in Newick format

**Format:** Newick, with or without branch lengths

**Examples:**
```
(species_A:0.01,species_B:0.02,(species_C_fused:0.003,species_D_fused:0.005):0.015);
(species_A,species_B,(species_C_fused,species_D_fused));
```

**Notes:**
- This is the **constrained** tree that enforces monophyly of the fused clade
- Taxon names must match exactly those in combined.aln
- **Branch lengths are optional:**
  - If tree **has** branch lengths: Step 1 will be skipped (uses them for SeqGen)
  - If tree **lacks** branch lengths: Step 1 will estimate them using IQTree
  - To force re-estimation: use `-f` flag when running pipeline
- Can be derived from species tree with fused clade constrained to be monophyletic

### 3. `model.txt`

**Content:** Best-fit substitution model for the combined alignment (from ModelFinder)

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
- Model is extracted from IQTree log file (`.iqtree` output) for the combined alignment
- Look for line: `Best-fit model according to BIC: GTR+G4`
- **Do not include** rate matrix details or base frequencies on separate lines
- Format must match IQTree notation exactly (case-sensitive: `GTR`, not `gtr`)
- If model includes gamma rate variation (`+G`), the alpha value will be extracted from Step 1 output
- If model includes invariant sites (`+I`), this must be preserved in the model string

## Preparing Input Data from Real Analysis

### Create Combined Alignment
```bash
# Your real analysis has block1.aln and block2.aln
# Combine them into a single alignment with gaps showing gene presence:

python3 << 'EOF'
from Bio import AlignIO
import sys

# Read the two component alignments
block1 = AlignIO.read('my_analysis/RHG_001/block1.trimmed.aln', 'fasta')
block2 = AlignIO.read('my_analysis/RHG_001/block2.trimmed.aln', 'fasta')

block1_len = block1.get_alignment_length()
block2_len = block2.get_alignment_length()

# Create combined alignment
combined_seqs = {}

# Add block1 sequences
for seq_record in block1:
    combined_seqs[seq_record.id] = str(seq_record.seq) + '-' * block2_len

# Add/update with block2 sequences
for seq_record in block2:
    if seq_record.id in combined_seqs:
        # Taxon in both blocks (fused)
        combined_seqs[seq_record.id] = combined_seqs[seq_record.id][:block1_len] + str(seq_record.seq)
    else:
        # Taxon only in block2
        combined_seqs[seq_record.id] = '-' * block1_len + str(seq_record.seq)

# Write combined alignment
with open('inputs/RHG_001/combined.aln', 'w') as f:
    for taxon, seq in sorted(combined_seqs.items()):
        f.write(f'>{taxon}\n{seq}\n')
EOF
```

### Extract Substitution Model
```bash
# Extract best-fit model from IQTree output:
# Use whichever block's model you prefer (or average them)
grep "Best-fit model according to BIC" my_analysis/RHG_001/block1.iqtree | awk '{print $NF}' > inputs/RHG_001/model.txt
```

### Get Constrained Tree
```bash
# Use your constrained tree from real analysis
# Should enforce monophyly of fused clade
cp my_analysis/RHG_001/constrained.nwk inputs/RHG_001/constrained.nwk
```

## Validation

A validation script checks all inputs before starting the pipeline:

```bash
python3 src/validate_inputs.py inputs/RHG_001
```

This checks:
- ✓ All required files exist (combined.aln, constrained.nwk, model.txt)
- ✓ Alignment format is readable
- ✓ Taxon names match between alignment and tree
- ✓ Tree is valid Newick format
- ✓ Model string is valid IQTree notation
- ✓ Gap patterns show mixed gene composition (at least one fused + one parent-only taxon)

## Template

A template RHG with dummy data is provided in `inputs/template/` for reference. Copy and modify it:

```bash
cp -r inputs/template inputs/RHG_001
# Then edit the combined.aln, constrained.nwk, and model.txt with your actual data
```
