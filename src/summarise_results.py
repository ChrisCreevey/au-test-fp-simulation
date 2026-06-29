#!/usr/bin/env python3
"""
Summarize AU test results across all replicates.
Computes false positive rates and generates null distribution.
"""

import sys
import argparse
import logging
import re
from pathlib import Path
from typing import Dict, List
from collections import defaultdict


def extract_au_pvalue(au_text_file: str) -> float:
    """
    Extract AU test p-value from IQTree output.
    Returns p-value or None if not found.
    """
    try:
        with open(au_text_file, 'r') as f:
            content = f.read()

        # Look for AU test p-value line
        # Pattern: "AU      0.05  0.001  0.50  0.999"
        match = re.search(r'AU\s+([\d.]+)', content)
        if match:
            return float(match.group(1))

        return None
    except Exception as e:
        logging.debug(f"Could not extract AU p-value from {au_text_file}: {e}")
        return None


def summarize_per_rhg(
    step5_output_dir: Path,
    rhg_id: str,
    output_dir: Path
) -> Dict:
    """
    Summarize results for one RHG across all replicates.
    Returns dict with false positive rates and bootstrap support.
    """

    results = {
        'rhg_id': rhg_id,
        'block1_au_rejects': 0,
        'block2_au_rejects': 0,
        'block1_total': 0,
        'block2_total': 0,
        'block1_bootstrap_support': [],
        'block2_bootstrap_support': [],
    }

    # Find all replicate directories
    rep_dirs = sorted([d for d in step5_output_dir.iterdir() if d.is_dir() and d.name.startswith('rep_')])

    logging.info(f"Found {len(rep_dirs)} replicates for {rhg_id}")

    for rep_dir in rep_dirs:
        rep_name = rep_dir.name

        # ====================================================================
        # Block 1 AU Test
        # ====================================================================

        au_file1 = rep_dir / 'au_test_block1.txt'
        if au_file1.exists():
            au_pval1 = extract_au_pvalue(str(au_file1))

            if au_pval1 is not None:
                results['block1_total'] += 1

                # AU test rejects monophyly if p < 0.05
                if au_pval1 < 0.05:
                    results['block1_au_rejects'] += 1
                    logging.debug(f"{rep_name} block1: AU rejected (p={au_pval1:.4f})")
                else:
                    logging.debug(f"{rep_name} block1: AU accepted (p={au_pval1:.4f})")

        # ====================================================================
        # Block 2 AU Test
        # ====================================================================

        au_file2 = rep_dir / 'au_test_block2.txt'
        if au_file2.exists():
            au_pval2 = extract_au_pvalue(str(au_file2))

            if au_pval2 is not None:
                results['block2_total'] += 1

                if au_pval2 < 0.05:
                    results['block2_au_rejects'] += 1
                    logging.debug(f"{rep_name} block2: AU rejected (p={au_pval2:.4f})")
                else:
                    logging.debug(f"{rep_name} block2: AU accepted (p={au_pval2:.4f})")

    # ========================================================================
    # Calculate False Positive Rates
    # ========================================================================

    results['block1_fp_rate'] = (
        100.0 * results['block1_au_rejects'] / results['block1_total']
        if results['block1_total'] > 0 else None
    )

    results['block2_fp_rate'] = (
        100.0 * results['block2_au_rejects'] / results['block2_total']
        if results['block2_total'] > 0 else None
    )

    return results


def write_results_tsv(results: Dict, output_file: Path) -> None:
    """Write results to TSV file."""

    with open(output_file, 'w') as f:
        # Header
        f.write("RHG_ID\tBlock\tFP_Rate_%\tAU_Rejects\tTotal_Replicates\n")

        # Block 1 results
        if results['block1_total'] > 0:
            fp_rate = results['block1_fp_rate'] or 0
            f.write(
                f"{results['rhg_id']}\tBlock1\t{fp_rate:.2f}\t"
                f"{results['block1_au_rejects']}\t{results['block1_total']}\n"
            )

        # Block 2 results
        if results['block2_total'] > 0:
            fp_rate = results['block2_fp_rate'] or 0
            f.write(
                f"{results['rhg_id']}\tBlock2\t{fp_rate:.2f}\t"
                f"{results['block2_au_rejects']}\t{results['block2_total']}\n"
            )

        # Combined (if both blocks present)
        if results['block1_total'] > 0 and results['block2_total'] > 0:
            total_rejects = results['block1_au_rejects'] + results['block2_au_rejects']
            total_reps = results['block1_total'] + results['block2_total']
            combined_rate = 100.0 * total_rejects / total_reps
            f.write(
                f"{results['rhg_id']}\tCombined\t{combined_rate:.2f}\t"
                f"{total_rejects}\t{total_reps}\n"
            )


def write_summary(results: Dict, output_file: Path) -> None:
    """Write human-readable summary."""

    with open(output_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write(f"AU Test False Positive Rate Summary: {results['rhg_id']}\n")
        f.write("=" * 70 + "\n\n")

        # Block 1
        if results['block1_total'] > 0:
            fp_rate = results['block1_fp_rate'] or 0
            f.write(f"Block 1:\n")
            f.write(f"  Total replicates: {results['block1_total']}\n")
            f.write(f"  AU test rejections (p < 0.05): {results['block1_au_rejects']}\n")
            f.write(f"  False positive rate: {fp_rate:.2f}%\n\n")

        # Block 2
        if results['block2_total'] > 0:
            fp_rate = results['block2_fp_rate'] or 0
            f.write(f"Block 2:\n")
            f.write(f"  Total replicates: {results['block2_total']}\n")
            f.write(f"  AU test rejections (p < 0.05): {results['block2_au_rejects']}\n")
            f.write(f"  False positive rate: {fp_rate:.2f}%\n\n")

        # Combined
        if results['block1_total'] > 0 and results['block2_total'] > 0:
            total_rejects = results['block1_au_rejects'] + results['block2_au_rejects']
            total_reps = results['block1_total'] + results['block2_total']
            combined_rate = 100.0 * total_rejects / total_reps
            f.write(f"Combined (Block 1 + Block 2):\n")
            f.write(f"  Total replicates: {total_reps}\n")
            f.write(f"  Total AU rejections: {total_rejects}\n")
            f.write(f"  Combined false positive rate: {combined_rate:.2f}%\n\n")

        f.write("=" * 70 + "\n")
        f.write("Interpretation:\n")
        f.write("  Under a null model of single-origin gene fusion, the AU test\n")
        f.write("  should reject monophyly approximately 5% of the time due to\n")
        f.write("  stochastic variation. Rates substantially higher than 5%\n")
        f.write("  suggest true signal of polyphyly in the real data.\n")
        f.write("=" * 70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Summarize AU test results and generate null distribution'
    )
    parser.add_argument('--step5-output', required=True,
                        help='Step 5 output directory (with AU test results)')
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for summary files')
    parser.add_argument('--rhg-id', required=True,
                        help='RHG identifier')
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

    step5_dir = Path(args.step5_output)
    if not step5_dir.exists():
        logging.error(f"Step 5 output directory not found: {step5_dir}")
        sys.exit(1)

    # ========================================================================
    # Summarize Results
    # ========================================================================

    logging.info(f"Summarizing results for {args.rhg_id}")

    results = summarize_per_rhg(step5_dir, args.rhg_id, output_dir)

    # ========================================================================
    # Write Output Files
    # ========================================================================

    tsv_file = output_dir / 'results.tsv'
    summary_file = output_dir / 'summary.txt'

    write_results_tsv(results, tsv_file)
    write_summary(results, summary_file)

    logging.info(f"Results written to:")
    logging.info(f"  - {tsv_file}")
    logging.info(f"  - {summary_file}")

    logging.info("")
    if results['block1_total'] > 0:
        logging.info(f"Block 1 FP rate: {results['block1_fp_rate']:.2f}%")
    if results['block2_total'] > 0:
        logging.info(f"Block 2 FP rate: {results['block2_fp_rate']:.2f}%")

    sys.exit(0)


if __name__ == '__main__':
    main()
