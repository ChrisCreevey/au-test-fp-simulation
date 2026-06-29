#!/usr/bin/env python3
"""
Apply single-origin fusion masking pattern to simulated sequences.

For each replicate:
  - Fusion clade taxa: concatenate block1 + block2 sequences into one "fused" sequence
  - Non-fusion taxa: keep block1 and block2 as separate sequences

This creates the null model where the fusion has a single origin.
"""

import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, Set, List, Tuple


def read_phylip_seqgen(filename: str) -> Dict[int, Dict[str, str]]:
    """
    Read SeqGen output (Phylip format with multiple datasets).
    Returns: {replicate_number: {taxon: sequence}}
    """
    replicates = {}
    current_rep = None
    current_seqs = {}
    current_taxon = None

    try:
        with open(filename, 'r') as f:
            lines = f.readlines()

        i = 0
        rep_num = 1

        while i < len(lines):
            line = lines[i].strip()

            # Header line: <ntaxa> <nchars>
            if line and not line[0].isspace() and ' ' in line:
                parts = line.split()
                if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                    if current_seqs:
                        replicates[rep_num - 1] = current_seqs
                    rep_num += 1
                    current_seqs = {}
                    current_taxon = None
                    i += 1
                    continue

            # Sequence line
            if line:
                if not line[0].isspace():
                    # New taxon
                    parts = line.split(None, 1)
                    if len(parts) == 2:
                        current_taxon = parts[0]
                        current_seqs[current_taxon] = parts[1]
                else:
                    # Continuation of previous taxon
                    if current_taxon:
                        current_seqs[current_taxon] += line.strip()

            i += 1

        # Add last replicate
        if current_seqs:
            replicates[rep_num - 1] = current_seqs

    except Exception as e:
        logging.error(f"Error reading SeqGen output {filename}: {e}")
        return None

    return replicates if replicates else None


def write_fasta(sequences: Dict[str, str], output_file: str) -> None:
    """Write sequences in FASTA format."""
    with open(output_file, 'w') as f:
        for taxon, seq in sorted(sequences.items()):
            f.write(f">{taxon}\n")
            # Write sequence in 80-char lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def apply_fusion_mask(
    block1_seqs: Dict[str, str],
    block2_seqs: Dict[str, str],
    fusion_taxa: Set[str],
    rep_num: int,
    output_dir: Path
) -> Tuple[bool, str]:
    """
    Apply masking for one replicate.
    Returns (success, message)
    """

    try:
        output_block1 = output_dir / f"rep_{rep_num:03d}_block1.aln"
        output_block2 = output_dir / f"rep_{rep_num:03d}_block2.aln"

        # For fusion clade: create concatenated "fused" sequence
        # For non-fusion: keep separate
        masked_block1 = {}
        masked_block2 = {}

        all_taxa = set(block1_seqs.keys())

        for taxon in all_taxa:
            seq1 = block1_seqs.get(taxon, '')
            seq2 = block2_seqs.get(taxon, '')

            if taxon in fusion_taxa:
                # This taxon carries the fusion
                # Concatenate into block1, leave block2 empty
                masked_block1[taxon] = seq1 + seq2
                # Don't include in block2 (or use placeholder)
            else:
                # Non-fusion taxon: keep blocks separate
                masked_block1[taxon] = seq1
                masked_block2[taxon] = seq2

        # Only write block1 if there are sequences
        if masked_block1:
            write_fasta(masked_block1, str(output_block1))

        if masked_block2:
            write_fasta(masked_block2, str(output_block2))

        return True, f"rep_{rep_num:03d} masked successfully"

    except Exception as e:
        return False, f"rep_{rep_num:03d} failed: {str(e)}"


def main():
    parser = argparse.ArgumentParser(
        description='Apply single-origin fusion masking to simulated sequences'
    )
    parser.add_argument('--sim-block1', required=True,
                        help='SeqGen simulated block 1 sequences')
    parser.add_argument('--sim-block2', required=True,
                        help='SeqGen simulated block 2 sequences')
    parser.add_argument('--fusion-taxa-file', required=True,
                        help='File with fusion clade taxa (one per line)')
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for masked alignments')
    parser.add_argument('--n-replicates', type=int, default=100,
                        help='Number of replicates')
    parser.add_argument('--log', default=None,
                        help='Log file path')

    args = parser.parse_args()

    # Setup logging
    if args.log:
        logging.basicConfig(
            level=logging.INFO,
            format='%(message)s',
            handlers=[
                logging.FileHandler(args.log, mode='a'),
                logging.StreamHandler()
            ]
        )
    else:
        logging.basicConfig(level=logging.INFO, format='%(message)s')

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ========================================================================
    # Load Fusion Taxa
    # ========================================================================

    logging.info("Reading fusion taxa...")

    fusion_taxa = set()
    try:
        with open(args.fusion_taxa_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    fusion_taxa.add(line)
    except Exception as e:
        logging.error(f"Failed to read fusion taxa: {e}")
        sys.exit(1)

    logging.info(f"Fusion taxa (n={len(fusion_taxa)}): {', '.join(sorted(fusion_taxa)[:5])}...")

    # ========================================================================
    # Load SeqGen Output
    # ========================================================================

    logging.info(f"Reading SeqGen output...")

    block1_reps = read_phylip_seqgen(args.sim_block1)
    block2_reps = read_phylip_seqgen(args.sim_block2)

    if not block1_reps or not block2_reps:
        logging.error("Failed to read SeqGen output files")
        sys.exit(1)

    logging.info(f"Block 1: {len(block1_reps)} replicates")
    logging.info(f"Block 2: {len(block2_reps)} replicates")

    n_reps = min(len(block1_reps), len(block2_reps), args.n_replicates)
    logging.info(f"Processing {n_reps} replicates")

    # ========================================================================
    # Apply Masking
    # ========================================================================

    success_count = 0
    fail_count = 0

    for rep_idx in range(n_reps):
        if rep_idx % 10 == 0:
            logging.info(f"  Processing replicate {rep_idx + 1} / {n_reps}")

        # SeqGen replicates are 0-indexed in dict
        block1_seqs = block1_reps.get(rep_idx, {})
        block2_seqs = block2_reps.get(rep_idx, {})

        if not block1_seqs or not block2_seqs:
            logging.warning(f"Replicate {rep_idx} missing sequence data, skipping")
            continue

        success, message = apply_fusion_mask(
            block1_seqs,
            block2_seqs,
            fusion_taxa,
            rep_idx + 1,  # Use 1-indexed numbering for output files
            output_dir
        )

        if success:
            success_count += 1
        else:
            logging.error(message)
            fail_count += 1

    # ========================================================================
    # Report Results
    # ========================================================================

    logging.info("")
    logging.info("=" * 60)
    logging.info(f"Masking Complete")
    logging.info(f"  Success: {success_count} / {n_reps}")
    logging.info(f"  Failed: {fail_count} / {n_reps}")

    if fail_count > 0:
        sys.exit(1)

    sys.exit(0)


if __name__ == '__main__':
    main()
