#!/usr/bin/env python3
"""
Simulate a domain's replicates with per-site rate variation concentrated at
the SAME real positions the alignment actually shows it, instead of spread
i.i.d. across all sites under a single smoothed gamma+invariant summary.

Background: simulating a domain's real alignment
length under its correctly-fitted branch lengths and a correctly-fitted
global alpha/pinvar still reproduces far less pairwise sequence identity
than the real data shows (e.g. 14.6% vs 47.2% for one real domain). Ruled
out: SeqGen/AliSim implementation bugs, guessed vs fitted rate parameters,
FreeRate-vs-Gamma branch-length inflation, and model distortion from the
monophyly constraint. What's left: a single global rate-heterogeneity
summary is too coarse. Real alignments typically have a MINORITY of sites
carrying most of the real conservation signal (near-invariant), with the
rest genuinely free to vary -- a global alpha spreads variation more evenly
across all sites than that, diluting exactly the signal a tree-reconstruction
method needs.

Fix: use IQTree's own per-site empirical-Bayes rate-category assignment
(-wsr, run once as part of Step 1's already-cheap fixed-topology fit) to
learn which REAL columns are in which rate category, then simulate each
category's sites SEPARATELY (as a block of uniform-rate sites, via AliSim's
--branch-scale) and stitch the results back together at their original
positions. This reproduces the real conservation PATTERN, not just its
average.

Usage:
  python3 src/simulate_by_site_rate.py \
      --tree brlen.treefile --base-model LG \
      --site-rates site_rates.txt --n-replicates 10 --seed 12345 \
      --iqtree-bin /path/to/iqtree2 --work-dir /tmp/work \
      --output simulated_domain.fa
"""

import argparse
import subprocess
import sys
from pathlib import Path
from collections import defaultdict


def read_site_rates(path):
    """Parse an IQTree .rate file (from -wsr).
    Returns ({category: [1-indexed site numbers]}, {category: rate})."""
    cats = defaultdict(list)
    cat_rates = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('Site'):
                continue
            parts = line.strip().split('\t')
            if len(parts) != 4:
                continue
            site, rate, cat, crate = parts
            cats[cat].append(int(site))
            cat_rates[cat] = float(crate)
    return cats, cat_rates


def read_fasta_records(path):
    """Ordered list of (taxon, seq), preserving duplicates (multi-replicate files)."""
    records = []
    cur_taxon = None
    cur_seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if cur_taxon is not None:
                    records.append((cur_taxon, ''.join(cur_seq)))
                cur_taxon = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_taxon is not None:
            records.append((cur_taxon, ''.join(cur_seq)))
    return records


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tree', required=True)
    ap.add_argument('--base-model', required=True, help='Substitution model with no rate-heterogeneity suffix, e.g. LG')
    ap.add_argument('--site-rates', required=True, help="IQTree's .rate file (from -wsr)")
    ap.add_argument('--n-replicates', type=int, required=True)
    ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--iqtree-bin', required=True)
    ap.add_argument('--work-dir', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    cats, cat_rates = read_site_rates(args.site_rates)
    if not cats:
        print("ERROR: no sites read from site-rates file", file=sys.stderr)
        sys.exit(1)

    total_sites = sum(len(v) for v in cats.values())
    print(f"Loaded {total_sites} sites across {len(cats)} rate categories: "
          f"{ {c: len(s) for c, s in cats.items()} }", file=sys.stderr)

    # ------------------------------------------------------------------
    # Simulate each category's sites as one block, all N_REPLICATES at once
    # ------------------------------------------------------------------
    cat_records = {}
    for cat, sites in sorted(cats.items()):
        n_sites = len(sites)
        rate = cat_rates[cat]
        # AliSim rejects a literal 0 branch-scale ("<SCALE> must be
        # positive!"); a tiny positive value is functionally invariant
        # (confirmed: produces identical sequences across taxa at that
        # scale, including across a taxon with an otherwise very long
        # branch) without triggering the error.
        scale = rate if rate > 0 else 1e-7
        prefix = work_dir / f"cat{cat}"
        cmd = [
            args.iqtree_bin, "--alisim", str(prefix),
            "-t", args.tree,
            "-m", args.base_model,
            "--length", str(n_sites),
            "--branch-scale", str(scale),
            "--seed", str(args.seed),
            "--num-alignments", str(args.n_replicates),
            "--single-output",
            "-af", "fasta",
        ]
        print(f"Category {cat}: {n_sites} sites, rate={rate}", file=sys.stderr)
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if result.returncode != 0:
            print(f"ERROR: AliSim failed for category {cat}:\n{result.stdout}", file=sys.stderr)
            sys.exit(1)

        out_file = f"{prefix}.fa"
        if not Path(out_file).exists():
            print(f"ERROR: expected output file not created: {out_file}", file=sys.stderr)
            sys.exit(1)
        cat_records[cat] = read_fasta_records(out_file)

    # ------------------------------------------------------------------
    # Sanity check: every category's block must have the same taxon set,
    # in the same order (same tree -> same fixed traversal order every
    # call), and the same number of records (n_taxa * n_replicates).
    # ------------------------------------------------------------------
    cat_list = sorted(cat_records.keys())
    first_cat = cat_list[0]
    n_total_records = len(cat_records[first_cat])
    if n_total_records % args.n_replicates != 0:
        print(f"ERROR: category {first_cat} produced {n_total_records} records, "
              f"not a multiple of n_replicates={args.n_replicates}", file=sys.stderr)
        sys.exit(1)
    n_taxa = n_total_records // args.n_replicates

    reference_taxa_order = [t for t, _ in cat_records[first_cat][:n_taxa]]
    for cat in cat_list[1:]:
        if len(cat_records[cat]) != n_total_records:
            print(f"ERROR: category {cat} has {len(cat_records[cat])} records, "
                  f"expected {n_total_records} (mismatched with category {first_cat})", file=sys.stderr)
            sys.exit(1)
        this_taxa_order = [t for t, _ in cat_records[cat][:n_taxa]]
        if this_taxa_order != reference_taxa_order:
            print(f"ERROR: category {cat}'s taxon order doesn't match category {first_cat}'s -- "
                  f"AliSim's per-call traversal order assumption doesn't hold", file=sys.stderr)
            sys.exit(1)

    # ------------------------------------------------------------------
    # Map each real site position -> (category, index within that category's block)
    # ------------------------------------------------------------------
    site_to_catpos = {}
    for cat, sites in cats.items():
        for pos, site in enumerate(sites):
            site_to_catpos[site] = (cat, pos)
    max_site = max(site_to_catpos.keys())
    if max_site != total_sites:
        print(f"WARN: max site number ({max_site}) != total sites counted ({total_sites}) "
              f"-- site numbering may have gaps", file=sys.stderr)

    # ------------------------------------------------------------------
    # Stitch: for each replicate, for each taxon, rebuild the full-length
    # sequence by placing each category's simulated characters back at
    # their original real positions.
    # ------------------------------------------------------------------
    with open(args.output, 'w') as out:
        for rep_idx in range(args.n_replicates):
            rep_blocks = {}
            for cat in cat_list:
                block = cat_records[cat][rep_idx * n_taxa:(rep_idx + 1) * n_taxa]
                rep_blocks[cat] = dict(block)

            for taxon in reference_taxa_order:
                full_seq = [''] * max_site
                for site in range(1, max_site + 1):
                    cat, pos = site_to_catpos[site]
                    full_seq[site - 1] = rep_blocks[cat][taxon][pos]
                out.write(f">{taxon}\n{''.join(full_seq)}\n")

    print(f"Wrote {args.n_replicates} stitched replicates ({n_taxa} taxa, {max_site} sites) to {args.output}",
          file=sys.stderr)


if __name__ == '__main__':
    main()
