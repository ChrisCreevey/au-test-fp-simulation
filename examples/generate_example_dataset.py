#!/usr/bin/env python3
"""
Generate example dataset for testing AU test simulation pipeline.

Creates a clean example with:
- 10 taxa: 3 Parent_gene_1, 3 Parent_gene_2, 4 Fused
- Combined alignment: first 50% = parent1, second 50% = parent2
- Balanced tree with monophyletic fused clade
- Appropriate gap patterns
"""

import random
import math
from pathlib import Path

# Set seed for reproducibility
random.seed(42)

# ============================================================================
# Configuration
# ============================================================================

OUTPUT_DIR = Path(__file__).parent
TAXA = {
    'P1_G1': 'Parent_gene_1',
    'P1_G2': 'Parent_gene_1',
    'P1_G3': 'Parent_gene_1',
    'P2_G1': 'Parent_gene_2',
    'P2_G2': 'Parent_gene_2',
    'P2_G3': 'Parent_gene_2',
    'FG_G1': 'fused',
    'FG_G2': 'fused',
    'FG_G3': 'fused',
    'FG_G4': 'fused',
}

PARENT1_LEN = 500  # Length of parent_gene_1 region
PARENT2_LEN = 500  # Length of parent_gene_2 region
TOTAL_LEN = PARENT1_LEN + PARENT2_LEN

# ============================================================================
# Generate Random Sequences (Amino Acids)
# ============================================================================

def generate_sequence(length: int, seed_offset: int = 0) -> str:
    """Generate a random amino acid sequence."""
    random.seed(seed_offset)
    # Standard 20 amino acids
    amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
    return ''.join(random.choice(amino_acids) for _ in range(length))


def evolve_sequence(ancestor: str, branch_length: float, mutation_rate: float = 0.02) -> str:
    """
    Evolve a sequence down a branch.
    branch_length controls number of mutations (approx).
    mutation_rate is per-site per-unit-branch-length (for amino acids).
    """
    seq = list(ancestor)
    num_mutations = int(len(seq) * branch_length * mutation_rate)
    amino_acids = list('ACDEFGHIKLMNPQRSTVWY')

    for _ in range(num_mutations):
        pos = random.randint(0, len(seq) - 1)
        available_aas = [aa for aa in amino_acids if aa != seq[pos]]
        seq[pos] = random.choice(available_aas)

    return ''.join(seq)


# ============================================================================
# Create Sequences
# ============================================================================

print("Generating example dataset (amino acids)...")
print(f"  Parent_gene_1 region: {PARENT1_LEN} aa")
print(f"  Parent_gene_2 region: {PARENT2_LEN} aa")
print(f"  Total alignment length: {TOTAL_LEN} aa")
print()

# Create ancestral sequences
ancestor_p1 = generate_sequence(PARENT1_LEN, seed_offset=1)
ancestor_p2 = generate_sequence(PARENT2_LEN, seed_offset=2)

# Evolve sequences down the tree
sequences = {}

# P1 clade (P1_G1, (P1_G2, P1_G3))
# Branch to P1_G1: length 0.05
sequences['P1_G1'] = evolve_sequence(ancestor_p1, 0.05)

# P1 internal node
p1_internal = evolve_sequence(ancestor_p1, 0.02)
# P1_G2 and P1_G3: branch length 0.03
sequences['P1_G2'] = evolve_sequence(p1_internal, 0.03)
sequences['P1_G3'] = evolve_sequence(p1_internal, 0.03)

# P2 clade (P2_G1, (P2_G2, P2_G3))
p2_internal = evolve_sequence(ancestor_p2, 0.02)
sequences['P2_G1'] = evolve_sequence(p2_internal, 0.05)
# P2_G2 and P2_G3
p2_internal2 = evolve_sequence(p2_internal, 0.02)
sequences['P2_G2'] = evolve_sequence(p2_internal2, 0.03)
sequences['P2_G3'] = evolve_sequence(p2_internal2, 0.03)

# Fused clade ((FG_G1, FG_G2), (FG_G3, FG_G4))
fg_ancestor_p1 = evolve_sequence(ancestor_p1, 0.08)
fg_ancestor_p2 = evolve_sequence(ancestor_p2, 0.08)

# FG internal
fg_internal1_p1 = evolve_sequence(fg_ancestor_p1, 0.02)
fg_internal1_p2 = evolve_sequence(fg_ancestor_p2, 0.02)

fg_internal2_p1 = evolve_sequence(fg_ancestor_p1, 0.02)
fg_internal2_p2 = evolve_sequence(fg_ancestor_p2, 0.02)

# FG leaves
sequences['FG_G1'] = evolve_sequence(fg_internal1_p1, 0.03) + evolve_sequence(fg_internal1_p2, 0.03)
sequences['FG_G2'] = evolve_sequence(fg_internal1_p1, 0.03) + evolve_sequence(fg_internal1_p2, 0.03)
sequences['FG_G3'] = evolve_sequence(fg_internal2_p1, 0.03) + evolve_sequence(fg_internal2_p2, 0.03)
sequences['FG_G4'] = evolve_sequence(fg_internal2_p1, 0.03) + evolve_sequence(fg_internal2_p2, 0.03)

# ============================================================================
# Create Combined Alignment with Gap Pattern
# ============================================================================

print("Creating combined alignment with gap pattern...")

combined_alignment = {}

for taxon, gene_type in TAXA.items():
    if gene_type == 'Parent_gene_1':
        # Only parent1: sequence in first half, gaps in second half
        combined_alignment[taxon] = sequences[taxon][:PARENT1_LEN] + '-' * PARENT2_LEN
    elif gene_type == 'Parent_gene_2':
        # Only parent2: gaps in first half, sequence in second half
        combined_alignment[taxon] = '-' * PARENT1_LEN + sequences[taxon][:PARENT2_LEN]
    else:  # fused
        # Both: full sequence (first half parent1, second half parent2)
        combined_alignment[taxon] = sequences[taxon][:PARENT1_LEN] + sequences[taxon][PARENT1_LEN:]

# Write combined alignment
combined_file = OUTPUT_DIR / 'combined.aln'
with open(combined_file, 'w') as f:
    for taxon in sorted(combined_alignment.keys()):
        f.write(f">{taxon}\n{combined_alignment[taxon]}\n")

print(f"  Written: {combined_file}")

# ============================================================================
# Create Tree with Branch Lengths
# ============================================================================

print("Creating constrained tree...")

# Tree structure from user with balanced branch lengths
tree_str = "((P1_G1:0.05,(P1_G2:0.03,P1_G3:0.03):0.02):0.1,((P2_G1:0.05,(P2_G2:0.03,P2_G3:0.03):0.02):0.1,((FG_G1:0.03,FG_G2:0.03):0.02,(FG_G3:0.03,FG_G4:0.03):0.02):0.08):0.05):0.0;"

tree_file = OUTPUT_DIR / 'constrained.nwk'
with open(tree_file, 'w') as f:
    f.write(tree_str)

print(f"  Written: {tree_file}")

# ============================================================================
# Create Model File
# ============================================================================

print("Creating model file...")

model_file = OUTPUT_DIR / 'model.txt'
with open(model_file, 'w') as f:
    f.write("GTR+G4\n")

print(f"  Written: {model_file}")

# ============================================================================
# Create Presence/Absence File
# ============================================================================

print("Creating presence/absence file...")

pa_file = OUTPUT_DIR / 'presence_absence.txt'
with open(pa_file, 'w') as f:
    f.write("species\tgene_type\n")
    for taxon in sorted(TAXA.keys()):
        f.write(f"{taxon}\t{TAXA[taxon]}\n")

print(f"  Written: {pa_file}")

# ============================================================================
# Summary
# ============================================================================

print()
print("=" * 70)
print("Example Dataset Summary")
print("=" * 70)
print()
print("Tree structure:")
print("  ((P1_G1,(P1_G2,P1_G3)),((P2_G1,(P2_G2,P2_G3)),((FG_G1,FG_G2),(FG_G3,FG_G4))))")
print()
print("Taxa:")
for gene_type in ['Parent_gene_1', 'Parent_gene_2', 'fused']:
    taxa = [t for t, gt in TAXA.items() if gt == gene_type]
    print(f"  {gene_type}: {', '.join(taxa)}")
print()
print("Alignment structure (amino acids):")
print(f"  Total length: {TOTAL_LEN} aa")
print(f"  Parent_gene_1 region: positions 1-{PARENT1_LEN}")
print(f"  Parent_gene_2 region: positions {PARENT1_LEN+1}-{TOTAL_LEN}")
print()
print("Gap pattern:")
print(f"  P1_G* taxa: sequences in region 1, gaps in region 2")
print(f"  P2_G* taxa: gaps in region 1, sequences in region 2")
print(f"  FG_G* taxa: sequences in both regions (concatenated)")
print()
print("Files created:")
print(f"  ✓ {combined_file}")
print(f"  ✓ {tree_file}")
print(f"  ✓ {model_file}")
print(f"  ✓ {pa_file}")
print()
print("To test the pipeline:")
print(f"  bash scripts/master_pipeline.sh -i {OUTPUT_DIR}")
print()
