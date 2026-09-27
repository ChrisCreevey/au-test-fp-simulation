#!/usr/bin/env python3
"""
Apply gap pattern from one domain's original alignment to its simulated
replicates.

The original per-domain alignment (e.g. F100058_Domain1.trim) contains
ordinary indel gaps -- the normal product of aligning real, variable-length
sequences. SeqGen/AliSim simulate full alignment-length, ungapped-but-
POSITIONALLY-SYNCHRONIZED sequences for every taxon (column i in every
taxon's output descends from the same ancestral column i), so this script:
1. Reads the original domain alignment to extract each taxon's gap pattern
2. Reads the simulator's FASTA output (N_REPLICATES datasets, one
   contiguous block of records per replicate)
3. Overlays each taxon's real gap pattern onto its simulated sequence, IN
   PLACE at the same column positions (see apply_gap_pattern()'s docstring
   for a critical bug this replaced: an earlier version extracted a
   taxon's first N simulated characters instead of the characters actually
   at its real positions, silently misaligning almost every comparison
   between taxa with different gap counts)
4. Writes one masked alignment per replicate
"""

import sys
import argparse
import logging
import time
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


def apply_gap_pattern(sim_seq: str, gap_positions: List[int], aln_length: int) -> str:
    """
    Overlay a taxon's real gap pattern onto its simulated sequence.

    CRITICAL: SeqGen/AliSim output is already full alignment-length and
    POSITIONALLY SYNCHRONIZED across every taxon -- position i in every
    taxon's simulated sequence descends from the same ancestral position i,
    evolved via that position's own site-specific rate. This is NOT a
    variable-length "ungapped" sequence needing gaps reinserted; it is
    already the right length and already in the right column order.

    An earlier version of this function treated sim_seq as if it needed
    gaps re-inserted -- extracting sim_seq's FIRST N characters (N = this
    taxon's real non-gap count) and placing THOSE at the taxon's real
    non-gap positions. That silently misaligns almost every comparison
    between two taxa with different gap counts (comparing arbitrary,
    evolutionarily-unrelated column pairs instead of the same real
    position), which for a sparse/mosaic alignment (most taxa cover
    different, only partially-overlapping regions) corrupts nearly every
    pairwise comparison. Confirmed directly: at a branch-scale near zero
    (sequences should be ~100% identical to the ancestor and each other),
    the old logic gave ~4% identity between two taxa with different gap
    counts; positionally-correct masking gives the expected ~98%.

    Args:
        sim_seq: Simulated sequence, same length as the original alignment,
            already positionally aligned to it (index i in sim_seq
            corresponds to alignment column i)
        gap_positions: Positions in the original alignment that are gaps
            for this taxon
        aln_length: Length of the original alignment

    Returns:
        sim_seq with '-' at this taxon's real gap positions, unchanged
        elsewhere -- a direct positional overlay, nothing extracted or
        reordered.
    """
    gap_pos_set = set(gap_positions)
    return ''.join('-' if pos in gap_pos_set else sim_seq[pos] for pos in range(aln_length))


def read_fasta_seqgen_records(filename: str) -> List[Tuple[str, str]]:
    """
    Read SeqGen's FASTA output as an ORDERED list of (taxon, sequence)
    records, preserving duplicates. SeqGen (this build's default output
    format, confirmed by direct testing -- not the Phylip-with-headers
    format an earlier version of this script assumed) writes N_REPLICATES
    datasets back to back, each a contiguous block of one record per taxon,
    always in the same fixed tree-traversal order.
    """
    records = []
    current_taxon = None
    current_seq = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_taxon is not None:
                    records.append((current_taxon, ''.join(current_seq)))
                current_taxon = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_taxon is not None:
            records.append((current_taxon, ''.join(current_seq)))

    return records


def group_into_replicates(records: List[Tuple[str, str]], n_taxa: int,
                           expected_taxa: set) -> Dict[int, Dict[str, str]]:
    """
    Chunk an ordered (taxon, seq) record list into per-replicate dicts,
    n_taxa records per replicate. Validates each block's taxon set matches
    the original alignment's -- if SeqGen ever changed its per-replicate
    ordering this would silently mix taxa across replicates, so we fail
    loudly instead.
    Returns: {replicate_index (0-based): {taxon: sequence}}
    """
    if n_taxa == 0:
        raise ValueError("n_taxa is 0, cannot chunk into replicates")
    if len(records) % n_taxa != 0:
        raise ValueError(
            f"SeqGen output has {len(records)} sequences, not a multiple of "
            f"n_taxa={n_taxa} -- cannot cleanly chunk into replicate blocks"
        )

    n_reps = len(records) // n_taxa
    replicates = {}
    for rep_idx in range(n_reps):
        block = dict(records[rep_idx * n_taxa:(rep_idx + 1) * n_taxa])
        if set(block.keys()) != expected_taxa:
            raise ValueError(
                f"Replicate {rep_idx} taxon set doesn't match the original "
                f"alignment's -- SeqGen output isn't in the fixed "
                f"per-replicate blocks this script assumes"
            )
        replicates[rep_idx] = block

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
                        help='SeqGen simulated sequences (ungapped, FASTA, N_REPLICATES blocks)')
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

    start_time = time.time()

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

    for taxon, seq in original.items():
        gap_patterns[taxon] = extract_gap_pattern(seq)

    logging.info(f"Loaded {len(original)} taxa from original alignment")
    logging.info(f"Alignment length: {len(list(original.values())[0])}")

    # ========================================================================
    # Load SeqGen Output
    # ========================================================================

    logging.info(f"Reading SeqGen output...")

    records = read_fasta_seqgen_records(args.simulated_seqs)

    if not records:
        logging.error("Failed to read SeqGen output")
        sys.exit(1)

    try:
        sim_reps = group_into_replicates(records, len(original), set(original.keys()))
    except ValueError as e:
        logging.error(f"Failed to group SeqGen output into replicates: {e}")
        sys.exit(1)

    logging.info(f"Loaded {len(sim_reps)} replicates from SeqGen")

    if len(sim_reps) != args.n_replicates:
        logging.warning(
            f"Expected {args.n_replicates} replicates, SeqGen produced {len(sim_reps)}"
        )

    n_reps = len(sim_reps)
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
                aln_len = len(list(original.values())[0])  # All sequences same length

                # Apply original gap pattern to simulated sequence
                masked_seq = apply_gap_pattern(sim_seq, gap_pos, aln_len)
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

    elapsed = time.time() - start_time

    logging.info("")
    logging.info("=" * 60)
    logging.info(f"Gap Pattern Application Complete")
    logging.info(f"  Success: {success_count} / {n_reps}")
    logging.info(f"  Failed: {fail_count} / {n_reps}")
    logging.info(f"  Output files: rep_NNN_masked.aln")
    logging.info(f"Time taken: {elapsed:.6f} seconds")

    if fail_count > 0:
        sys.exit(1)

    sys.exit(0)


if __name__ == '__main__':
    main()
