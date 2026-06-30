#!/usr/bin/env python3
"""
Validate input data for AU test simulation pipeline.
Checks:
  - All required files exist
  - Alignment format is readable
  - Taxon names match between files
  - Tree is valid Newick format
  - Model string is valid IQTree notation
  - Gap patterns show mixed gene composition
"""

import sys
import os
import argparse
from pathlib import Path


def read_fasta(filename):
    """Read FASTA file and return dict of {taxon: sequence}."""
    seqs = {}
    current_taxon = None
    current_seq = []

    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_taxon:
                        seqs[current_taxon] = ''.join(current_seq)
                    current_taxon = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_taxon:
                seqs[current_taxon] = ''.join(current_seq)
    except Exception as e:
        print(f"ERROR reading FASTA {filename}: {e}", file=sys.stderr)
        return None

    return seqs


def read_newick(filename):
    """Read Newick format tree."""
    try:
        with open(filename, 'r') as f:
            content = f.read().strip()
            if not content.endswith(';'):
                content += ';'
            return content
    except Exception as e:
        print(f"ERROR reading Newick {filename}: {e}", file=sys.stderr)
        return None


def extract_taxa_from_newick(newick):
    """Extract taxon names from Newick tree."""
    taxa = set()
    current_taxon = []

    for char in newick:
        if char in '():;,':
            if current_taxon:
                taxon = ''.join(current_taxon).strip()
                taxon = taxon.split(':')[0]
                if taxon and taxon not in ['', ' ']:
                    taxa.add(taxon)
                current_taxon = []
        elif char != ' ':
            current_taxon.append(char)
        else:
            if current_taxon and current_taxon[-1] not in ['(', ',', ')']:
                current_taxon.append(char)

    return taxa


def validate_model(model_str):
    """Check if model string is valid IQTree notation."""
    valid_models = [
        'JC', 'F81', 'K80', 'HKY', 'TrN', 'TPM1', 'TPM2', 'TPM3',
        'TVM', 'GTR', 'SYM'
    ]

    base_model = model_str.split('+')[0]

    if base_model not in valid_models:
        print(f"WARNING: Unknown base model '{base_model}' in '{model_str}'", file=sys.stderr)
        return False

    return True


def analyze_gap_patterns(seqs):
    """
    Analyze gap patterns to infer gene composition.
    Returns dict: {taxon: 'parent_1', 'parent_2', or 'both'}
    """
    if not seqs:
        return {}

    # Find regions by looking for gap patterns across all sequences
    aln_length = len(list(seqs.values())[0])
    gap_counts = [0] * aln_length

    for seq in seqs.values():
        for i, char in enumerate(seq):
            if char == '-':
                gap_counts[i] += 1

    # Find likely boundaries (regions where almost all taxa have gaps)
    # This is a heuristic approach
    taxa_composition = {}

    for taxon, seq in seqs.items():
        gap_positions = [i for i, c in enumerate(seq) if c == '-']

        if not gap_positions:
            # No gaps = has both regions
            taxa_composition[taxon] = 'both'
        elif len(gap_positions) == aln_length:
            # All gaps = absent (shouldn't happen)
            taxa_composition[taxon] = 'absent'
        else:
            # Some gaps - need to check if they're at start or end
            gap_ranges = []
            in_gap = False
            gap_start = None

            for i in range(aln_length):
                if seq[i] == '-':
                    if not in_gap:
                        gap_start = i
                        in_gap = True
                else:
                    if in_gap:
                        gap_ranges.append((gap_start, i - 1))
                        in_gap = False

            if in_gap:
                gap_ranges.append((gap_start, aln_length - 1))

            # If gaps are only at the start, it's parent_2
            # If gaps are only at the end, it's parent_1
            # Otherwise, complex pattern
            if len(gap_ranges) == 1:
                gap_start, gap_end = gap_ranges[0]
                if gap_start == 0:
                    taxa_composition[taxon] = 'parent_2'
                elif gap_end == aln_length - 1:
                    taxa_composition[taxon] = 'parent_1'
                else:
                    taxa_composition[taxon] = 'mixed'
            else:
                taxa_composition[taxon] = 'mixed'

    return taxa_composition


def main():
    parser = argparse.ArgumentParser(
        description='Validate input data for AU test simulation pipeline'
    )
    parser.add_argument('input_dir', help='Input RHG directory')
    parser.add_argument('--strict', action='store_true', help='Exit on any warning')

    args = parser.parse_args()
    input_dir = Path(args.input_dir)

    if not input_dir.exists():
        print(f"ERROR: Input directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)

    errors = []
    warnings = []

    print(f"Validating inputs in: {input_dir}")
    print("-" * 60)

    # ========================================================================
    # Check Required Files
    # ========================================================================

    required_files = [
        'combined.aln',
        'constrained.nwk',
        'model.txt'
    ]

    for fname in required_files:
        fpath = input_dir / fname
        if not fpath.exists():
            errors.append(f"Missing required file: {fname}")
        else:
            print(f"✓ {fname} exists")

    if errors:
        print("\nFATAL ERRORS:", file=sys.stderr)
        for err in errors:
            print(f"  ✗ {err}", file=sys.stderr)
        sys.exit(1)

    # ========================================================================
    # Read and Validate Alignment
    # ========================================================================

    print("\nValidating combined alignment...")

    aln = read_fasta(str(input_dir / 'combined.aln'))

    if not aln:
        errors.append("Failed to read combined.aln (invalid FASTA format)")
    else:
        aln_length = len(list(aln.values())[0])
        n_taxa = len(aln)
        print(f"✓ combined.aln: {n_taxa} sequences, {aln_length} bp (with gaps)")

    if errors:
        print("\nFATAL ERRORS:", file=sys.stderr)
        for err in errors:
            print(f"  ✗ {err}", file=sys.stderr)
        sys.exit(1)

    # ========================================================================
    # Validate Tree
    # ========================================================================

    print("\nValidating tree...")

    tree_content = read_newick(str(input_dir / 'constrained.nwk'))
    if not tree_content:
        errors.append("Failed to read constrained.nwk")
    else:
        tree_taxa = extract_taxa_from_newick(tree_content)
        print(f"✓ constrained.nwk: {len(tree_taxa)} taxa")
        print(f"  Tree: {tree_content[:80]}...")

    if errors:
        print("\nFATAL ERRORS:", file=sys.stderr)
        for err in errors:
            print(f"  ✗ {err}", file=sys.stderr)
        sys.exit(1)

    # ========================================================================
    # Validate Model
    # ========================================================================

    print("\nValidating model...")

    with open(input_dir / 'model.txt', 'r') as f:
        model = f.read().strip()

    if not validate_model(model):
        warnings.append(f"Unusual model: {model}")
    else:
        print(f"✓ model.txt: {model}")

    # ========================================================================
    # Analyze Gap Patterns
    # ========================================================================

    print("\nAnalyzing gap patterns...")

    composition = analyze_gap_patterns(aln)

    parent_1_count = sum(1 for c in composition.values() if c == 'parent_1')
    parent_2_count = sum(1 for c in composition.values() if c == 'parent_2')
    both_count = sum(1 for c in composition.values() if c == 'both')
    mixed_count = sum(1 for c in composition.values() if c == 'mixed')

    print(f"✓ Gene composition:")
    print(f"    Parent_gene_1 only: {parent_1_count}")
    print(f"    Parent_gene_2 only: {parent_2_count}")
    print(f"    Both (fused): {both_count}")
    print(f"    Mixed/complex: {mixed_count}")

    # Check for realistic composition
    if both_count == 0:
        errors.append("No taxa with both genes found (need fused genes)")

    if parent_1_count + parent_2_count == 0:
        errors.append("No parent-only taxa found (need non-fused genes)")

    if errors:
        print("\nFATAL ERRORS:", file=sys.stderr)
        for err in errors:
            print(f"  ✗ {err}", file=sys.stderr)
        sys.exit(1)

    # ========================================================================
    # Cross-Validation: Taxon Name Consistency
    # ========================================================================

    print("\nValidating taxon name consistency...")

    aln_taxa = set(aln.keys())

    if tree_taxa and tree_taxa != aln_taxa:
        warnings.append("Tree taxa don't match alignment taxa")
        print(f"  tree only: {tree_taxa - aln_taxa}")
        print(f"  aln only: {aln_taxa - tree_taxa}")
    else:
        print("✓ All taxon names match between alignment and tree")

    # ========================================================================
    # Report Results
    # ========================================================================

    print("\n" + "=" * 60)

    if errors:
        print("VALIDATION FAILED - FATAL ERRORS:", file=sys.stderr)
        for err in errors:
            print(f"  ✗ {err}", file=sys.stderr)
        sys.exit(1)

    if warnings:
        print("VALIDATION PASSED with warnings:")
        for warn in warnings:
            print(f"  ⚠ {warn}")
        if args.strict:
            sys.exit(1)
    else:
        print("VALIDATION PASSED - All checks successful ✓")

    sys.exit(0)


if __name__ == '__main__':
    main()
