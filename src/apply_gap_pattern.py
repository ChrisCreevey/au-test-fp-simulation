#!/usr/bin/env python3
"""
Apply gap pattern from original alignment to simulated sequences.

The original alignment contains gaps that indicate which regions each taxon has.
This script:
1. Reads the original alignment to extract gap patterns per taxon
2. Reads SeqGen output (ungapped simulated sequences)
3. Inserts gaps at the same positions as in the original
4. Outputs masked alignments that mimic the original structure

This preserves the real gene distribution (which taxa have parent_gene_1,
parent_gene_2, or both) while using simulated sequences.
"""

import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple


def read_fasta(filename: str) -> Dict[str, str]:
    """Read FASTA alignment and return {taxon: sequence}."""
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
        logging.error(f"Error reading FASTA {filename}: {e}")
        return None

    return seqs


def extract_gap_pattern(sequence: str) -> List[int]:
    """
    Extract gap positions from a sequence.
    Returns list of indices where gaps occur.
    """
    return [i for i, char in enumerate(sequence) if char == '-']


def remove_gaps(sequence: str) -> str:
    """Remove all gaps from sequence."""
    return sequence.replace('-', '')


def apply_gap_pattern(ungapped_seq: str, gap_positions: List[int]) -> str:
    """
    Insert gaps into ungapped sequence at specified positions.

    Args:
        ungapped_seq: Sequence without gaps
        gap_positions: List of indices where gaps should be inserted

    Returns:
        Sequence with gaps inserted
    """
    if not gap_positions:
        return ungapped_seq

    # Build result by iterating through original positions
    result = []
    ungapped_idx = 0
    all_positions = set(range(len(ungapped_seq) + len(gap_positions)))
    gap_pos_set = set(gap_positions)

    for pos in sorted(all_positions):
        if pos in gap_pos_set:
            result.append('-')
        else:
            if ungapped_idx < len(ungapped_seq):
                result.append(ungapped_seq[ungapped_idx])
                ungapped_idx += 1

    return ''.join(result)


def read_phylip_seqgen(filename: str) -> Dict[int, Dict[str, str]]:
    """
    Read SeqGen output (Phylip format with multiple datasets).
    Returns: {replicate_number: {taxon: sequence}}
    """
    replicates = {}
    current_seqs = {}
    current_taxon = None
    rep_num = 1

    try:
        with open(filename, 'r') as f:
            lines = f.readlines()

        i = 0
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

    return replicates


def write_fasta(sequences: Dict[str, str], output_file: str) -> None:
    """Write sequences in FASTA format."""
    with open(output_file, 'w') as f:
        for taxon, seq in sorted(sequences.items()):
            f.write(f">{taxon}\n")
            # Write sequence in 80-char lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Apply gap pattern from original alignment to simulated sequences'
    )
    parser.add_argument('--original-aln', required=True,
                        help='Original alignment with gap patterns')
    parser.add_argument('--simulated-seqs', required=True,
                        help='SeqGen simulated sequences (ungapped)')
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
    # Load Original Alignment and Extract Gap Patterns
    # ========================================================================

    logging.info("Reading original alignment and extracting gap patterns...")

    original = read_fasta(args.original_aln)

    if not original:
        logging.error("Failed to read original alignment")
        sys.exit(1)

    # Extract gap patterns per taxon
    gap_patterns = {}
    ungapped_seqs = {}

    for taxon, seq in original.items():
        gap_patterns[taxon] = extract_gap_pattern(seq)
        ungapped_seqs[taxon] = remove_gaps(seq)

    logging.info(f"Loaded {len(original)} taxa from original alignment")
    logging.info(f"Alignment length (with gaps): {len(list(original.values())[0])}")
    logging.info(f"Alignment length (ungapped): {len(list(ungapped_seqs.values())[0])}")

    # ========================================================================
    # Load SeqGen Output
    # ========================================================================

    logging.info(f"Reading SeqGen output...")

    sim_reps = read_phylip_seqgen(args.simulated_seqs)

    if not sim_reps:
        logging.error("Failed to read SeqGen output")
        sys.exit(1)

    logging.info(f"Loaded {len(sim_reps)} replicates from SeqGen")

    n_reps = min(len(sim_reps), args.n_replicates)
    logging.info(f"Processing {n_reps} replicates")

    # ========================================================================
    # Apply Gap Pattern to Each Replicate
    # ========================================================================

    success_count = 0
    fail_count = 0

    for rep_idx in range(n_reps):
        if rep_idx % 10 == 0:
            logging.info(f"  Processing replicate {rep_idx + 1} / {n_reps}")

        rep_num = rep_idx + 1  # Use 1-indexed numbering
        sim_seqs = sim_reps.get(rep_idx, {})

        if not sim_seqs:
            logging.warning(f"Replicate {rep_idx} has no sequences, skipping")
            continue

        try:
            # Apply gap patterns
            masked_seqs = {}

            for taxon in original.keys():
                if taxon not in sim_seqs:
                    logging.warning(f"Taxon {taxon} not in SeqGen output for rep {rep_num}")
                    continue

                sim_seq = sim_seqs[taxon]
                gap_pos = gap_patterns[taxon]

                # Apply original gap pattern to simulated sequence
                masked_seq = apply_gap_pattern(sim_seq, gap_pos)
                masked_seqs[taxon] = masked_seq

            # Write masked alignment
            output_file = output_dir / f"rep_{rep_num:03d}_masked.aln"
            write_fasta(masked_seqs, str(output_file))

            success_count += 1

        except Exception as e:
            logging.error(f"Failed to process replicate {rep_num}: {e}")
            fail_count += 1

    # ========================================================================
    # Report Results
    # ========================================================================

    logging.info("")
    logging.info("=" * 60)
    logging.info(f"Gap Pattern Application Complete")
    logging.info(f"  Success: {success_count} / {n_reps}")
    logging.info(f"  Failed: {fail_count} / {n_reps}")
    logging.info(f"  Output files: rep_NNN_masked.aln")

    if fail_count > 0:
        sys.exit(1)

    sys.exit(0)


if __name__ == '__main__':
    main()
