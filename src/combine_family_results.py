#!/usr/bin/env python3
"""
Combine every family's step6_family/family_summary.tsv (src/summarise_family.py
output) into one dataset-wide table.

Why not a plain `cat`/header-then-concat: summarise_family.py's header is
dynamic -- it has one Pct_Significant_<Domain>/Real_DeltaLogL_<Domain>/
Permutation_P_Upper_<Domain> column per domain the family actually has, and
most families have 2 domains but some (treated as ordinary representatives
of the family, see docs/INPUTS.md) have 3+. A naive concatenation would silently
misalign columns between a 2-domain and a 3-domain family's rows. This does a
column-name union instead: every row gets every column that ANY family
produced, "NA" where that family didn't have that domain.

Usage:
  python3 src/combine_family_results.py --output-base data/sim_outputs_bulk \
      --out data/sim_outputs_bulk/combined_family_results.tsv
"""
import argparse
import csv
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output-base", required=True, help="Directory containing one subdirectory per family (e.g. data/sim_outputs_bulk)")
    ap.add_argument("--out", required=True, help="Path for the combined TSV")
    args = ap.parse_args()

    output_base = Path(args.output_base)
    summary_paths = sorted(output_base.glob("*/step6_family/family_summary.tsv"))

    if not summary_paths:
        raise SystemExit(f"ERROR: no family_summary.tsv files found under {output_base}/*/step6_family/")

    rows = []
    all_cols = []
    seen_cols = set()

    for path in summary_paths:
        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                rows.append(row)
                for col in reader.fieldnames:
                    if col not in seen_cols:
                        seen_cols.add(col)
                        all_cols.append(col)

    # Fixed leading columns first (present in every file), then whatever
    # per-domain columns showed up, in first-seen order.
    leading = [c for c in ("Family_ID", "N_Domains", "N_Replicates_Paired") if c in seen_cols]
    trailing = [c for c in all_cols if c not in leading]
    header = leading + trailing

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t", restval="NA")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Combined {len(summary_paths)} family_summary.tsv file(s) -> "
          f"{len(rows)} families, {len(header)} columns written to: {args.out}")


if __name__ == "__main__":
    main()
