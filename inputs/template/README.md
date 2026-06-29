# Input Data Template

This directory contains template files showing the required format for input data to the AU test simulation pipeline.

## Files in This Template

### `block1.aln` and `block2.aln`
- **Format**: FASTA (DNA sequences)
- **Content**: Trimmed multiple sequence alignments for each component homologous block
- **Requirements**:
  - All taxa must be present in both alignments
  - Taxon names must match exactly between blocks and tree file
  - Sequences must be aligned (same length per alignment)
  - Already trimmed (no gap-heavy columns)

### `constrained.nwk`
- **Format**: Newick tree format with branch lengths
- **Content**: Monophyletic species tree where fusion clade is constrained to be monophyletic
- **Requirements**:
  - Taxon names must exactly match alignment and presence_absence.txt
  - All taxa must be present in the tree
  - Should include branch lengths (will be re-estimated in Step 1)

### `model_block1.txt` and `model_block2.txt`
- **Format**: Plain text, one model string per file
- **Content**: Best-fit substitution model (from IQTree ModelFinder)
- **Examples**:
  - `GTR+G4` (GTR model with 4 gamma categories)
  - `HKY+G4` (HKY model with gamma variation)
  - `JC` (Jukes-Cantor, no parameters)
  - `TVM+G4` (TVM with gamma)

### `presence_absence.txt`
- **Format**: Tab-separated text with header
- **Content**: Gene type for each taxon
- **Valid gene types**:
  - `Parent_gene_1` — only the first component gene
  - `Parent_gene_2` — only the second component gene
  - `fused` — both components concatenated
  - `absent` — not present in this gene family
- **Requirements**:
  - First line: header `species	gene_type`
  - Taxon names must exactly match those in alignments
  - At least one `fused` taxon required
  - At least one `Parent_gene_1` or `Parent_gene_2` taxon required (non-fused)

## Preparing Your Data

1. **Copy this template** to start a new RHG:
   ```bash
   cp -r inputs/template inputs/RHG_001
   ```

2. **Replace each file** with your actual data:
   ```bash
   # From your real analysis:
   cp /path/to/analysis/RHG_001/block1.trimmed.aln inputs/RHG_001/block1.aln
   cp /path/to/analysis/RHG_001/block2.trimmed.aln inputs/RHG_001/block2.aln
   
   # Extract model from IQTree log:
   grep "Best-fit model" /path/to/RHG_001/block1.iqtree | awk '{print $NF}' > inputs/RHG_001/model_block1.txt
   ```

3. **Validate** your inputs:
   ```bash
   python3 src/validate_inputs.py inputs/RHG_001
   ```

## Directory Structure for Multiple RHGs

Once you have prepared data for multiple RHGs, your `inputs/` directory should look like:

```
inputs/
├── template/        (keep this as reference)
│   ├── block1.aln
│   ├── block2.aln
│   ├── constrained.nwk
│   ├── model_block1.txt
│   ├── model_block2.txt
│   └── presence_absence.txt
├── RHG_001/
│   ├── block1.aln
│   ├── block2.aln
│   ├── constrained.nwk
│   ├── model_block1.txt
│   ├── model_block2.txt
│   └── presence_absence.txt
├── RHG_002/
│   ├── (same 6 files)
├── RHG_003/
│   └── (same 6 files)
...
```

The pipeline will automatically discover all `RHG_*` directories and process them.
