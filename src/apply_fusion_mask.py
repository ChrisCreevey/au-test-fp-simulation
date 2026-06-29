#!/usr/bin/env python3
"""
Apply single-origin fusion masking pattern to simulated sequences.

Creates a SINGLE mixed alignment per replicate where:
  - Fused taxa: concatenated block1 + block2 sequences
  - Parent_gene_1 taxa: only block1 sequence
  - Parent_gene_2 taxa: only block2 sequence
  - Absent taxa: not included

This reflects the real data structure where different species carry
different versions of the gene(s).
"""

import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, Set


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


def read_presence_absence(filename: str) -> Dict[str, str]:
    """
    Read presence/absence matrix.
    Returns: {taxon: gene_type} where gene_type is one of:
      - Parent_gene_1, Parent_gene_2, fused, absent
    """
    pa_dict = {}

    try:
        with open(filename, 'r') as f:
            lines = f.readlines()

        # Skip header
        for line in lines[1:]:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) != 2:
                logging.warning(f"Skipping malformed line: {line}")
                continue

            taxon, gene_type = parts[0].strip(), parts[1].strip()
            pa_dict[taxon] = gene_type

    except Exception as e:
        logging.error(f"Error reading presence/absence file {filename}: {e}")
        return None

    return pa_dict


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
    presence_absence: Dict[str, str],
    rep_num: int,
    output_dir: Path
) -> tuple:
    """
    Apply masking for one replicate.
    Creates a single mixed alignment based on presence/absence.
    Returns (success, message)
    """

    try:
        output_mixed = output_dir / f"rep_{rep_num:03d}_mixed.aln"

        # Create mixed alignment
        mixed_alignment = {}

        for taxon, gene_type in presence_absence.items():
            if gene_type == 'absent':
                # Don't include absent taxa
                continue
            elif gene_type == 'fused':
                # Concatenate both blocks
                seq1 = block1_seqs.get(taxon, '')
                seq2 = block2_seqs.get(taxon, '')
                if seq1 and seq2:
                    mixed_alignment[taxon] = seq1 + seq2
            elif gene_type == 'Parent_gene_1':
                # Use only block1
                seq1 = block1_seqs.get(taxon, '')
                if seq1:
                    mixed_alignment[taxon] = seq1
            elif gene_type == 'Parent_gene_2':
                # Use only block2
                seq2 = block2_seqs.get(taxon, '')
                if seq2:
                    mixed_alignment[taxon] = seq2

        # Write mixed alignment
        if mixed_alignment:
            write_fasta(mixed_alignment, str(output_mixed))
            return True, f"rep_{rep_num:03d} masked successfully ({len(mixed_alignment)} taxa)"
        else:
            return False, f"rep_{rep_num:03d} produced empty alignment"

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
    parser.add_argument('--presence-absence-file', required=True,
                        help='File with gene presence/absence per taxon')
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
    # Load Presence/Absence Matrix
    # ========================================================================

    logging.info("Reading presence/absence matrix...")

    presence_absence = read_presence_absence(args.presence_absence_file)

    if not presence_absence:
        logging.error("Failed to read presence/absence file")
        sys.exit(1)

    fused_count = sum(1 for gt in presence_absence.values() if gt == 'fused')
    parent1_count = sum(1 for gt in presence_absence.values() if gt == 'Parent_gene_1')
    parent2_count = sum(1 for gt in presence_absence.values() if gt == 'Parent_gene_2')

    logging.info(f"Presence/absence matrix loaded:")
    logging.info(f"  Fused taxa: {fused_count}")
    logging.info(f"  Parent_gene_1 only: {parent1_count}")
    logging.info(f"  Parent_gene_2 only: {parent2_count}")

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
            presence_absence,
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
    logging.info(f"  Output files: rep_NNN_mixed.aln")

    if fail_count > 0:
        sys.exit(1)

    sys.exit(0)


if __name__ == '__main__':
    main()
