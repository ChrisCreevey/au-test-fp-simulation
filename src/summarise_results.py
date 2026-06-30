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
) -> dict:
    """
    Summarize results for one RHG across all replicates.
    Returns dict with false positive rates.
    """

    results = {
        'rhg_id': rhg_id,
        'au_rejects': 0,
        'total': 0,
    }

    # Find all replicate directories
    rep_dirs = sorted([d for d in step5_output_dir.iterdir() if d.is_dir() and d.name.startswith('rep_')])

    logging.info(f"Found {len(rep_dirs)} replicates for {rhg_id}")

    for rep_dir in rep_dirs:
        rep_name = rep_dir.name

        # ====================================================================
        # Mixed Alignment AU Test
        # ====================================================================

        # Try .txt first (legacy), then .iqtree (IQTree2 output)
        au_file = rep_dir / 'au_test_result.txt'
        if not au_file.exists():
            au_file = rep_dir / 'au_test_result.iqtree'

        if au_file.exists():
            au_pval = extract_au_pvalue(str(au_file))

            if au_pval is not None:
                results['total'] += 1

                # AU test rejects monophyly if p < 0.05
                if au_pval < 0.05:
                    results['au_rejects'] += 1
                    logging.debug(f"{rep_name}: AU rejected (p={au_pval:.4f})")
                else:
                    logging.debug(f"{rep_name}: AU accepted (p={au_pval:.4f})")

    # ========================================================================
    # Calculate False Positive Rate
    # ========================================================================

    results['fp_rate'] = (
        100.0 * results['au_rejects'] / results['total']
        if results['total'] > 0 else None
    )

    return results


def write_results_tsv(results: dict, output_file: Path) -> None:
    """Write results to TSV file."""

    with open(output_file, 'w') as f:
        # Header
        f.write("RHG_ID\tFP_Rate_%\tAU_Rejects\tTotal_Replicates\n")

        # Results
        if results['total'] > 0:
            fp_rate = results['fp_rate'] or 0
            f.write(
                f"{results['rhg_id']}\t{fp_rate:.2f}\t"
                f"{results['au_rejects']}\t{results['total']}\n"
            )


def write_summary(results: dict, output_file: Path) -> None:
    """Write human-readable summary."""

    with open(output_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write(f"AU Test False Positive Rate Summary: {results['rhg_id']}\n")
        f.write("=" * 70 + "\n\n")

        if results['total'] > 0:
            fp_rate = results['fp_rate'] or 0
            f.write(f"Total replicates: {results['total']}\n")
            f.write(f"AU test rejections (p < 0.05): {results['au_rejects']}\n")
            f.write(f"False positive rate: {fp_rate:.2f}%\n\n")

            f.write("=" * 70 + "\n")
            f.write("Interpretation:\n")
            f.write("  Under null model of single-origin fusion, the AU test\n")
            f.write("  should reject monophyly approximately 5% of the time (α=0.05).\n")
            f.write("  \n")
            if fp_rate < 2:
                f.write("  → Conservative: AU test rarely rejects monophyly\n")
            elif fp_rate < 7:
                f.write("  → Well-calibrated: AU test behaves as expected\n")
            else:
                f.write("  → Liberal: AU test may reject monophyly too often\n")
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
    if results['total'] > 0:
        logging.info(f"False positive rate: {results['fp_rate']:.2f}%")

    sys.exit(0)


if __name__ == '__main__':
    main()
