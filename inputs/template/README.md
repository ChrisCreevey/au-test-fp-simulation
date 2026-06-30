# Input Data Template

This directory contains template files showing the required format for input data to the AU test simulation pipeline.

## Files in This Template

### `combined.aln`
- **Format**: FASTA (DNA sequences)
- **Content**: Combined alignment with both parent gene regions
- **Structure**:
  - Single alignment containing the full sequence space of both parent genes
  - Gap characters (`-`) indicate missing regions
  - Parent_gene_1 only: gaps in the parent_gene_2 region
  - Parent_gene_2 only: gaps in the parent_gene_1 region
  - Fused: sequences for both regions (no gaps between)

**Example:**
```
>species_A
ATGCGATCG----------CGTAGGCTAGGC  (Parent_gene_1 only)
>species_B
----------CGTAGGCTAGGCTAGGCTAGGC  (Parent_gene_2 only)
>species_C_fused
ATGCGATCGCGTAGGCTAGGCTAGGCTAGGC   (Both genes - fused)
```

- **Requirements**:
  - All taxa must be present (gaps indicate missing regions)
  - Sequence names must match exactly with the constrained tree and presence_absence.txt
  - Alignment should already be trimmed (using trimAl `-gappyout`)

### `constrained.nwk`
- **Format**: Newick tree format with branch lengths
- **Content**: Monophyletic species tree where fused clade is constrained to be monophyletic
- **Requirements**:
  - Taxon names must exactly match alignment and presence_absence.txt
  - All taxa must be present in the tree
  - Should include branch lengths (will be re-estimated in Step 1 if missing)

### `model.txt`
- **Format**: Plain text, one model string
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
  - Taxon names must exactly match those in combined.aln and constrained.nwk
  - At least one `fused` taxon required
  - At least one `Parent_gene_1` or `Parent_gene_2` taxon required (non-fused)
- **Used for**:
  - Identifying fused taxa in AU test results
  - Checking monophyly of fused clade in reconstructed trees

## Preparing Your Data

1. **Copy this template** to start a new RHG:
   ```bash
   cp -r inputs/template inputs/RHG_001
   ```

2. **Create combined alignment** from your block1 and block2:
   ```bash
   python3 << 'EOF'
   from Bio import AlignIO
   
   # Read component alignments
   block1 = AlignIO.read('my_analysis/block1.aln', 'fasta')
   block2 = AlignIO.read('my_analysis/block2.aln', 'fasta')
   
   block1_len = block1.get_alignment_length()
   block2_len = block2.get_alignment_length()
   
   # Create combined
   combined_seqs = {}
   for seq in block1:
       combined_seqs[seq.id] = str(seq.seq) + '-' * block2_len
   for seq in block2:
       if seq.id in combined_seqs:
           combined_seqs[seq.id] = combined_seqs[seq.id][:block1_len] + str(seq.seq)
       else:
           combined_seqs[seq.id] = '-' * block1_len + str(seq.seq)
   
   # Write
   with open('inputs/RHG_001/combined.aln', 'w') as f:
       for taxon, seq in sorted(combined_seqs.items()):
           f.write(f'>{taxon}\n{seq}\n')
   EOF
   ```

3. **Extract model**:
   ```bash
   grep "Best-fit model" my_analysis/block1.iqtree | awk '{print $NF}' > inputs/RHG_001/model.txt
   ```

4. **Copy constrained tree**:
   ```bash
   cp my_analysis/constrained.nwk inputs/RHG_001/constrained.nwk
   ```

5. **Validate** your inputs:
   ```bash
   python3 src/validate_inputs.py inputs/RHG_001
   ```

## Directory Structure for Multiple RHGs

Once you have prepared data for multiple RHGs, your `inputs/` directory should look like:

```
inputs/
├── template/        (keep this as reference)
│   ├── combined.aln
│   ├── constrained.nwk
│   ├── model.txt
│   └── presence_absence.txt
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
├── RHG_003/
│   └── [same structure]
...
```

The pipeline will automatically discover all `RHG_*` directories and process them.
