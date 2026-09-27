#!/usr/bin/env python3
"""
Build a manifest of (family, domain) simulation units from Peter's real data.

Real data layout (per family F#####, per domain N):
  F#####_Domain{N}.trim              trimmed FASTA alignment for that domain
  F#####_Domain{N}.constr0.treefile  constrained (composite-clade-monophyletic)
                                      tree with branch lengths already estimated
  F#####_Domain{N}-geneID.tsv        per-taxon label: "composite" or "parent"
  CHG_model.json                     {"F#####_Domain{N}": "<IQTree model string>", ...}
                                      one file, shared across all families

A family can have any number of domains (in the current data drop it is always
2, but nothing here assumes that). This script discovers whatever domains are
actually present as files, ignoring any extra domain keys in CHG_model.json
that have no corresponding .trim/.treefile/.tsv triple -- these are
orphaned metadata from an earlier upstream step.

Output:
  --output-json  nested manifest: family -> {seed, domains: {n: {...}}, ...}
  --output-tsv   flat one-row-per-(family,domain) table for shell/xargs use
"""

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

DOMAIN_RE = re.compile(r'^(F\d+)_Domain(\d+)\.trim$')


def family_seed(family: str) -> int:
    """Deterministic per-family seed so every domain of a family can share
    the same SeqGen draw (replicate i means the same underlying simulation
    scenario across all of a family's domains)."""
    h = hashlib.sha256(family.encode()).hexdigest()
    return int(h[:8], 16) % (2**31 - 1) or 1  # avoid 0, some tools treat 0 as "unset"


def read_fasta_taxa_and_length(path: Path):
    taxa = []
    length = None
    cur_len = 0
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if length is None and taxa:
                    length = cur_len
                taxa.append(line[1:].split()[0])
                cur_len = 0
            else:
                cur_len += len(line)
        if length is None:
            length = cur_len
    return taxa, length


def read_geneid_categories(path: Path):
    """Returns {taxon: category}. Real data uses 'composite' / 'parent'."""
    cats = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) != 2:
                continue
            cats[parts[0]] = parts[1]
    return cats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--data-dir', required=True, help='Directory containing *_Domain*.trim etc. (conv_fusion_data)')
    ap.add_argument('--model-json', default=None, help='Path to CHG_model.json (default: <data-dir>/CHG_model.json)')
    ap.add_argument('--output-json', required=True)
    ap.add_argument('--output-tsv', required=True)
    ap.add_argument('--family', default=None, help='Restrict to one family ID (e.g. F100058), for quick testing')
    args = ap.parse_args()

    data_dir = Path(args.data_dir)
    model_json_path = Path(args.model_json) if args.model_json else data_dir / 'CHG_model.json'

    if not data_dir.is_dir():
        print(f"ERROR: data dir not found: {data_dir}", file=sys.stderr)
        sys.exit(1)

    with open(model_json_path) as f:
        models = json.load(f)

    # Discover families/domains from actual .trim files on disk (source of truth,
    # NOT the model json's key set -- see docstring).
    families = {}  # family -> {domain_num: trim_path}
    for p in sorted(data_dir.glob('*_Domain*.trim')):
        m = DOMAIN_RE.match(p.name)
        if not m:
            continue
        fam, dom = m.group(1), int(m.group(2))
        if args.family and fam != args.family:
            continue
        families.setdefault(fam, {})[dom] = p

    manifest = {}
    tsv_rows = []
    n_missing_companion = 0
    n_missing_model = 0
    n_composite_inconsistent = 0

    for fam in sorted(families.keys(), key=lambda s: int(s[1:])):
        domains_found = families[fam]
        seed = family_seed(fam)
        domain_entries = {}
        composite_sets = {}

        for dom_num, trim_path in sorted(domains_found.items()):
            stem = f"{fam}_Domain{dom_num}"
            tree_path = data_dir / f"{stem}.constr0.treefile"
            geneid_path = data_dir / f"{stem}-geneID.tsv"

            if not tree_path.exists() or not geneid_path.exists():
                n_missing_companion += 1
                print(f"WARN: {stem} missing companion file(s), skipping this domain", file=sys.stderr)
                continue

            model = models.get(stem)
            if model is None:
                n_missing_model += 1
                print(f"WARN: {stem} has no entry in CHG_model.json, skipping this domain", file=sys.stderr)
                continue

            taxa, aln_len = read_fasta_taxa_and_length(trim_path)
            cats = read_geneid_categories(geneid_path)
            composite = {t for t, c in cats.items() if c == 'composite'}
            parent = {t for t, c in cats.items() if c == 'parent'}
            composite_sets[dom_num] = composite

            domain_entries[str(dom_num)] = {
                "trim": str(trim_path),
                "tree": str(tree_path),
                "geneid": str(geneid_path),
                "model": model,
                "n_taxa": len(taxa),
                "aln_length": aln_len,
                "n_composite": len(composite),
                "n_parent": len(parent),
            }

            tsv_rows.append([
                fam, str(dom_num), str(trim_path), str(tree_path), str(geneid_path),
                model, str(seed), str(len(taxa)), str(aln_len), str(len(composite)), str(len(parent)),
            ])

        if not domain_entries:
            continue

        composite_union = set()
        for s in composite_sets.values():
            composite_union |= s
        consistent = len(composite_sets) <= 1 or all(s == composite_union for s in composite_sets.values())
        if not consistent:
            n_composite_inconsistent += 1

        manifest[fam] = {
            "seed": seed,
            "n_domains": len(domain_entries),
            "composite_taxa_consistent": consistent,
            "composite_taxa": sorted(composite_union),
            "domains": domain_entries,
        }

    Path(args.output_json).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_json, 'w') as f:
        json.dump(manifest, f, indent=1)

    Path(args.output_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output_tsv, 'w') as f:
        f.write('\t'.join([
            'family', 'domain', 'trim', 'tree', 'geneid', 'model', 'seed',
            'n_taxa', 'aln_length', 'n_composite', 'n_parent',
        ]) + '\n')
        for row in tsv_rows:
            f.write('\t'.join(row) + '\n')

    print(f"Families in manifest: {len(manifest)}", file=sys.stderr)
    print(f"Domain entries: {len(tsv_rows)}", file=sys.stderr)
    print(f"Domains skipped (missing companion file): {n_missing_companion}", file=sys.stderr)
    print(f"Domains skipped (missing model): {n_missing_model}", file=sys.stderr)
    print(f"Families with inconsistent composite-taxon sets across domains: {n_composite_inconsistent}", file=sys.stderr)


if __name__ == '__main__':
    main()
