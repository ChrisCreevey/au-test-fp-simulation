#!/usr/bin/env python3
"""
Calculate majority-rule consensus sequence from alignment.

For each column in the alignment:
  - Count all bases (excluding gaps)
  - Take the most common one
  - Result: consensus sequence with no gaps, same length as alignment

This consensus will be used as the ancestral sequence for SeqGen simulations.
"""

import sys
import argparse
from pathlib import Path
from collections import Counter


def read_fasta(filename: str) -> tuple:
    """
    Read FASTA alignment.
    Returns: (sequences_dict, alignment_length)
    """
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
        return None, None

    if not seqs:
        return None, None

    aln_length = len(list(seqs.values())[0])
    return seqs, aln_length


def calculate_consensus(seqs: dict) -> str:
    """
    Calculate majority-rule consensus sequence.

    For each position, take the most common base (excluding gaps).
    If all gaps at a position, use 'N' (shouldn't happen with trimmed alignment).
    """
    aln_length = len(list(seqs.values())[0])
    consensus = []

    for pos in range(aln_length):
        # Collect all bases at this position (excluding gaps)
        bases = []
        for taxon, seq in seqs.items():
            if pos < len(seq):
                base = seq[pos]
                if base != '-':
                    bases.append(base)

        if not bases:
            # All gaps at this position (shouldn't happen with trimmed alignment)
            consensus.append('N')
        else:
            # Most common base
            counter = Counter(bases)
            most_common = counter.most_common(1)[0][0]
            consensus.append(most_common)

    return ''.join(consensus)


def main():
    parser = argparse.ArgumentParser(
        description='Calculate majority-rule consensus sequence from alignment'
    )
    parser.add_argument('--alignment', '-a', required=True,
                        help='Input alignment (FASTA format)')
    parser.add_argument('--output', '-o', required=True,
                        help='Output consensus sequence file')

    args = parser.parse_args()

    # Read alignment
    print(f"Reading alignment: {args.alignment}")
    seqs, aln_length = read_fasta(args.alignment)

    if not seqs:
        print("ERROR: Failed to read alignment", file=sys.stderr)
        sys.exit(1)

    print(f"  {len(seqs)} sequences, {aln_length} positions")

    # Calculate consensus
    print("Calculating majority-rule consensus...")
    consensus = calculate_consensus(seqs)

    print(f"  Consensus length: {len(consensus)}")
    print(f"  Expected length: {aln_length}")

    if len(consensus) != aln_length:
        print("ERROR: Consensus length doesn't match alignment", file=sys.stderr)
        sys.exit(1)

    # Write consensus in SeqGen format (Phylip all taxa with consensus sequence)
    # SeqGen's -k flag (use taxon as ancestral) works better with all taxa present
    # All taxa get the consensus sequence; SeqGen will evolve from the first one
    with open(args.output, 'w') as f:
        f.write(f"{len(seqs)} {len(consensus)}\n")
        for taxon in sorted(seqs.keys()):
            f.write(f"{taxon} {consensus}\n")

    print(f"Consensus written to: {args.output}")
    print(f"  Format: SeqGen Phylip ({len(seqs)} sequences, all identical consensus)")
    print(f"  Sequence: {consensus[:50]}...")

    sys.exit(0)


if __name__ == '__main__':
    main()
