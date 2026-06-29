#!/usr/bin/env python3
"""
Validate input data for AU test simulation pipeline.
Checks:
  - All required files exist
  - Alignment format is readable
  - Taxon names match between files
  - Tree is valid Newick format
  - Model string is valid IQTree notation
  - Fusion taxa exist in alignments
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
                    current_taxon = line[1:].split()[0]  # Use first word as taxon name
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_taxon:
                seqs[current_taxon] = ''.join(current_seq)
    except Exception as e:
        print(f"ERROR reading FASTA {filename}: {e}", file=sys.stderr)
        return None

    return seqs


def read_phylip(filename):
    """Read Phylip format alignment."""
    seqs = {}

    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            if len(lines) < 2:
                print(f"ERROR: {filename} too short for Phylip format", file=sys.stderr)
                return None

            # Skip header line
            for line in lines[1:]:
                line = line.strip()
                if not line:
                    continue
                parts = line.split(None, 1)
                if len(parts) == 2:
                    taxon, seq = parts
                    seqs[taxon] = seq

    except Exception as e:
        print(f"ERROR reading Phylip {filename}: {e}", file=sys.stderr)
        return None

    return seqs if seqs else None


def read_alignment(filename):
    """Try to read alignment file in FASTA or Phylip format."""
    seqs = read_fasta(filename)
    if not seqs or len(seqs) == 0:
        seqs = read_phylip(filename)

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
                # Remove branch length notation
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
        print(f"  Valid models: {', '.join(valid_models)}", file=sys.stderr)
        return False

    return True


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
        'block1.aln',
        'block2.aln',
        'constrained.nwk',
        'model_block1.txt',
        'model_block2.txt',
        'fusion_taxa.txt'
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
    # Read and Validate Alignments
    # ========================================================================

    print("\nValidating alignments...")

    aln1 = read_alignment(str(input_dir / 'block1.aln'))
    aln2 = read_alignment(str(input_dir / 'block2.aln'))

    if not aln1:
        errors.append("Failed to read block1.aln (invalid FASTA/Phylip format)")
    else:
        print(f"✓ block1.aln: {len(aln1)} sequences, {len(list(aln1.values())[0])} bp")

    if not aln2:
        errors.append("Failed to read block2.aln (invalid FASTA/Phylip format)")
    else:
        print(f"✓ block2.aln: {len(aln2)} sequences, {len(list(aln2.values())[0])} bp")

    if errors:
        print("\nFATAL ERRORS:", file=sys.stderr)
        for err in errors:
            print(f"  ✗ {err}", file=sys.stderr)
        sys.exit(1)

    # ========================================================================
    # Validate Trees
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
    # Validate Models
    # ========================================================================

    print("\nValidating models...")

    with open(input_dir / 'model_block1.txt', 'r') as f:
        model1 = f.read().strip()
    with open(input_dir / 'model_block2.txt', 'r') as f:
        model2 = f.read().strip()

    if not validate_model(model1):
        warnings.append(f"Unusual model for block 1: {model1}")
    else:
        print(f"✓ model_block1.txt: {model1}")

    if not validate_model(model2):
        warnings.append(f"Unusual model for block 2: {model2}")
    else:
        print(f"✓ model_block2.txt: {model2}")

    # ========================================================================
    # Validate Fusion Taxa
    # ========================================================================

    print("\nValidating fusion taxa...")

    with open(input_dir / 'fusion_taxa.txt', 'r') as f:
        fusion_taxa = set(
            line.strip()
            for line in f
            if line.strip() and not line.startswith('#')
        )

    print(f"✓ fusion_taxa.txt: {len(fusion_taxa)} taxa")
    for taxon in sorted(fusion_taxa):
        print(f"    - {taxon}")

    # ========================================================================
    # Cross-Validation: Taxon Name Consistency
    # ========================================================================

    print("\nValidating taxon name consistency...")

    aln1_taxa = set(aln1.keys())
    aln2_taxa = set(aln2.keys())

    if aln1_taxa != aln2_taxa:
        warnings.append("Taxon sets differ between block1 and block2")
        print(f"  block1 only: {aln1_taxa - aln2_taxa}")
        print(f"  block2 only: {aln2_taxa - aln1_taxa}")

    if tree_taxa and tree_taxa != aln1_taxa:
        warnings.append("Tree taxa don't match alignment taxa")
        print(f"  tree only: {tree_taxa - aln1_taxa}")
        print(f"  aln only: {aln1_taxa - tree_taxa}")

    # Check fusion taxa exist in alignments
    missing_fusion_taxa = fusion_taxa - aln1_taxa
    if missing_fusion_taxa:
        errors.append(f"Fusion taxa not in alignments: {missing_fusion_taxa}")
    else:
        print("✓ All fusion taxa found in alignments")

    # Check that there are both fusion and non-fusion taxa
    non_fusion = aln1_taxa - fusion_taxa
    if not non_fusion:
        errors.append("No non-fusion taxa found (need at least one for comparison)")
    else:
        print(f"✓ Non-fusion taxa count: {len(non_fusion)}")

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
