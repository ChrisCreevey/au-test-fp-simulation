#!/usr/bin/env python3
"""
Summarize AU test results across all replicates for one gene family (RHG).
Produces a single-row, wide TSV designed to be concatenated across
thousands of gene families for cross-gene analysis.
"""

import sys
import argparse
import logging
import statistics
from pathlib import Path
from typing import Optional

NA = "NA"

# Columns written to results.tsv, in order. Kept as a module constant so the
# header and the row-writer can never drift apart.
COLUMNS = [
    "RHG_ID",
    "N_Taxa",
    "Alignment_Length",
    "N_Replicates",
    "N_Skipped_Replicates",
    "N_Mono_Violations",
    "Pct_Mono_Violations",
    "N_Significant_Violations",
    "Pct_Significant_Violations",
    "N_AU_Rejects_Total",
    "Pct_AU_Rejects_Total",
    "Mean_AU_Pvalue",
    "Mean_Tree_Length",
    "SD_Tree_Length",
]


def parse_result_line(line: str) -> Optional[dict]:
    """Parse one 'RESULT rep_xxx key=value key=value ...' line from step5_results.txt."""
    parts = line.split()
    if len(parts) < 2 or parts[0] != "RESULT":
        return None

    fields = {"rep_name": parts[1]}
    for token in parts[2:]:
        if "=" in token:
            key, value = token.split("=", 1)
            fields[key] = value
    return fields


def to_float(value: str) -> Optional[float]:
    if value is None or value == NA:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def summarize_per_rhg(step5_output_dir: Path, rhg_id: str) -> dict:
    """
    Summarize results for one RHG across all replicates by reading the
    per-replicate RESULT/SKIP lines Step 5 already wrote to step5_results.txt
    (rejected=, reason=, mono=, ntaxa=, nsites=, treelen=, pval=).
    """

    results_file = step5_output_dir / "step5_results.txt"

    n_replicates = 0
    n_skipped = 0
    n_mono_violations = 0
    n_significant_violations = 0
    n_au_rejects = 0
    ntaxa_values = []
    nsites_values = []
    pvals = []
    treelens = []

    if results_file.exists():
        for line in results_file.read_text().splitlines():
            if line.startswith("SKIP"):
                n_skipped += 1
                continue

            fields = parse_result_line(line)
            if fields is None:
                continue

            n_replicates += 1

            rejected = fields.get("rejected") == "1"
            mono = fields.get("mono")

            if rejected:
                n_au_rejects += 1

            if mono == "false":
                n_mono_violations += 1
                if rejected:
                    n_significant_violations += 1

            ntaxa = to_float(fields.get("ntaxa"))
            if ntaxa is not None:
                ntaxa_values.append(ntaxa)

            nsites = to_float(fields.get("nsites"))
            if nsites is not None:
                nsites_values.append(nsites)

            pval = to_float(fields.get("pval"))
            if pval is not None:
                pvals.append(pval)

            treelen = to_float(fields.get("treelen"))
            if treelen is not None:
                treelens.append(treelen)
    else:
        logging.warning(f"step5_results.txt not found in {step5_output_dir}")

    logging.info(f"Found {n_replicates} replicates for {rhg_id} ({n_skipped} skipped)")

    def pct(n):
        return 100.0 * n / n_replicates if n_replicates > 0 else None

    return {
        "RHG_ID": rhg_id,
        "N_Taxa": int(ntaxa_values[0]) if ntaxa_values else None,
        "Alignment_Length": int(nsites_values[0]) if nsites_values else None,
        "N_Replicates": n_replicates,
        "N_Skipped_Replicates": n_skipped,
        "N_Mono_Violations": n_mono_violations,
        "Pct_Mono_Violations": pct(n_mono_violations),
        "N_Significant_Violations": n_significant_violations,
        "Pct_Significant_Violations": pct(n_significant_violations),
        "N_AU_Rejects_Total": n_au_rejects,
        "Pct_AU_Rejects_Total": pct(n_au_rejects),
        "Mean_AU_Pvalue": statistics.mean(pvals) if pvals else None,
        "Mean_Tree_Length": statistics.mean(treelens) if treelens else None,
        "SD_Tree_Length": statistics.stdev(treelens) if len(treelens) > 1 else (0.0 if treelens else None),
    }


def format_value(value) -> str:
    if value is None:
        return NA
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def write_results_tsv(results: dict, output_file: Path) -> None:
    """Write a single-row TSV. Always emits a row (NA-filled if untestable) so
    concatenating results.tsv across many gene families yields one row per gene,
    including genes where the pipeline produced no usable replicates."""

    with open(output_file, "w") as f:
        f.write("\t".join(COLUMNS) + "\n")
        f.write("\t".join(format_value(results[col]) for col in COLUMNS) + "\n")


def interpret_fp_rate(results: dict) -> str:
    if results["N_Replicates"] == 0:
        return "untestable (no usable replicates)"
    fp_rate = results["Pct_Significant_Violations"] or 0
    if fp_rate < 2:
        return "conservative (rarely rejects monophyly)"
    elif fp_rate < 7:
        return "well-calibrated (behaves as expected under the null, ~5%)"
    else:
        return "liberal (rejects monophyly more often than expected)"


def write_summary(results: dict, output_file: Path) -> None:
    """Write a single-row (header + one data row) higher-level summary for this gene family."""

    with open(output_file, "w") as f:
        f.write("\t".join(COLUMNS) + "\tInterpretation\n")
        row = [format_value(results[col]) for col in COLUMNS] + [interpret_fp_rate(results)]
        f.write("\t".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Summarize AU test results and generate a per-gene-family summary row"
    )
    parser.add_argument("--step5-output", required=True,
                        help="Step 5 output directory (with step5_results.txt)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for summary files")
    parser.add_argument("--rhg-id", required=True,
                        help="RHG identifier")
    parser.add_argument("--log", default=None,
                        help="Log file path")

    args = parser.parse_args()

    if args.log:
        logging.basicConfig(
            level=logging.INFO,
            format="%(message)s",
            handlers=[
                logging.FileHandler(args.log, mode="a"),
                logging.StreamHandler()
            ]
        )
    else:
        logging.basicConfig(level=logging.INFO, format="%(message)s")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    step5_dir = Path(args.step5_output)
    if not step5_dir.exists():
        logging.error(f"Step 5 output directory not found: {step5_dir}")
        sys.exit(1)

    logging.info(f"Summarizing results for {args.rhg_id}")

    results = summarize_per_rhg(step5_dir, args.rhg_id)

    tsv_file = output_dir / "results.tsv"
    summary_file = output_dir / "summary.txt"

    write_results_tsv(results, tsv_file)
    write_summary(results, summary_file)

    logging.info("Results written to:")
    logging.info(f"  - {tsv_file}")
    logging.info(f"  - {summary_file}")

    logging.info("")
    if results["N_Replicates"] > 0:
        logging.info(f"Significant violation rate: {format_value(results['Pct_Significant_Violations'])}%")

    sys.exit(0)


if __name__ == "__main__":
    main()
