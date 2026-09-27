#!/usr/bin/env python3
"""
Join a family's per-domain Step 5 results (paired by replicate, via the
shared per-family seed) into family-level OR/AND consolidation rates.

Per the project's "report raw, decide downstream" design, Step 5 never
forces a single-domain-per-family verdict. This script is that downstream
step: for each replicate, a domain counts as a "significant violation" only
if the AU test rejected AND that rejection reflects real loss of monophyly
(rejected=1 AND mono=false) -- matching how Peter's real analysis presumably
calls a rejection (a break in monophyly), not any AU rejection regardless of
cause. See summarise_results.py's identical Pct_Significant_Violations
definition, reused here per-domain before joining.

OR-rate: fraction of replicates where AT LEAST ONE domain shows a
significant violation.
AND-rate: fraction of replicates where ALL domains show a significant
violation.

Usage:
  python3 src/summarise_family.py --family F100058 \
      --domain Domain1:/path/to/Domain1/step5 \
      --domain Domain2:/path/to/Domain2/step5 \
      --output-dir /path/to/output
"""

import argparse
import logging
import sys
from pathlib import Path


def parse_result_line(line: str):
    parts = line.split()
    if len(parts) < 2 or parts[0] != "RESULT":
        return None
    fields = {"rep_name": parts[1]}
    for token in parts[2:]:
        if "=" in token:
            key, value = token.split("=", 1)
            fields[key] = value
    return fields


def read_domain_replicates(step5_dir: Path):
    """Returns {rep_name: {'rejected': bool, 'mono': str, 'reason': str, 'pval': str, 'deltalogl': float|None}}."""
    results_file = step5_dir / "step5_results.txt"
    if not results_file.exists():
        logging.error(f"step5_results.txt not found in {step5_dir}")
        return {}

    reps = {}
    for line in results_file.read_text().splitlines():
        fields = parse_result_line(line)
        if fields is None:
            continue
        deltalogl = fields.get("deltalogl")
        try:
            deltalogl = float(deltalogl) if deltalogl not in (None, "NA") else None
        except ValueError:
            deltalogl = None
        reps[fields["rep_name"]] = {
            "rejected": fields.get("rejected") == "1",
            "mono": fields.get("mono"),
            "reason": fields.get("reason"),
            "pval": fields.get("pval"),
            "deltalogl": deltalogl,
        }
    return reps


def read_real_delta_logl(step5_dir: Path):
    """Real-data comparison lives at <domain_dir>/real_data/real_delta_logl.txt,
    a sibling of step5 (see scripts/real_data_delta_logl.sh). Optional --
    returns None if the (expensive, separate) real-data search hasn't been run."""
    path = step5_dir.parent / "real_data" / "real_delta_logl.txt"
    if not path.exists():
        return None
    for line in path.read_text().splitlines():
        if line.startswith("DELTA_LOGL="):
            try:
                return float(line.split("=", 1)[1])
            except ValueError:
                return None
    return None


def significant_violation(rep_data: dict) -> bool:
    """A domain's replicate counts as a real (monophyly-breaking) violation
    only if AU rejected AND the composite clade actually isn't monophyletic
    -- not any AU rejection regardless of cause (see module docstring)."""
    return rep_data["rejected"] and rep_data["mono"] == "false"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family", required=True)
    ap.add_argument("--domain", action="append", required=True,
                     help="NAME:STEP5_OUTPUT_DIR, repeatable, one per domain")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--log", default=None)
    args = ap.parse_args()

    if args.log:
        logging.basicConfig(level=logging.INFO, format="%(message)s",
                             handlers=[logging.FileHandler(args.log, mode="a"), logging.StreamHandler()])
    else:
        logging.basicConfig(level=logging.INFO, format="%(message)s")

    domains = {}
    domain_paths = {}
    for spec in args.domain:
        if ":" not in spec:
            logging.error(f"--domain must be NAME:PATH, got: {spec}")
            sys.exit(1)
        name, path = spec.split(":", 1)
        domain_paths[name] = Path(path)
        domains[name] = read_domain_replicates(Path(path))

    if not domains:
        logging.error("No domains provided")
        sys.exit(1)

    domain_names = sorted(domains.keys())

    # Only replicates present in EVERY domain can be jointly evaluated --
    # this should be all of them, given the shared-seed batch design, but
    # don't silently assume it.
    rep_sets = [set(d.keys()) for d in domains.values()]
    common_reps = sorted(set.intersection(*rep_sets)) if rep_sets else []
    all_reps = sorted(set.union(*rep_sets)) if rep_sets else []
    if len(common_reps) != len(all_reps):
        missing = set(all_reps) - set(common_reps)
        logging.warning(f"{len(missing)} replicate(s) not present in all domains, excluded from paired analysis: {sorted(missing)}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Per-replicate detail table
    # ------------------------------------------------------------------
    detail_file = output_dir / "family_replicate_detail.tsv"
    n_or = 0
    n_and = 0
    with open(detail_file, "w") as f:
        header = ["family", "replicate"] + [f"{d}_significant_violation" for d in domain_names] + ["OR", "AND"]
        f.write("\t".join(header) + "\n")
        for rep in common_reps:
            sig = {d: significant_violation(domains[d][rep]) for d in domain_names}
            or_result = any(sig.values())
            and_result = all(sig.values())
            n_or += or_result
            n_and += and_result
            row = [args.family, rep] + [str(int(sig[d])) for d in domain_names] + [str(int(or_result)), str(int(and_result))]
            f.write("\t".join(row) + "\n")

    n_reps = len(common_reps)

    # ------------------------------------------------------------------
    # Family-level summary row
    # ------------------------------------------------------------------
    def pct(n):
        return round(100.0 * n / n_reps, 4) if n_reps > 0 else "NA"

    per_domain_pct = {}
    for d in domain_names:
        n_sig = sum(1 for rep in common_reps if significant_violation(domains[d][rep]))
        per_domain_pct[d] = pct(n_sig)

    # Per-domain permutation comparison (real ΔlogL vs this domain's null
    # ΔlogL distribution) -- optional per domain, since it requires the
    # separate, expensive real-data unconstrained search to have been run
    # (scripts/real_data_delta_logl.sh). NA for any domain missing it.
    perm_p_by_domain = {}
    real_delta_by_domain = {}
    for d in domain_names:
        real_delta = read_real_delta_logl(domain_paths[d])
        real_delta_by_domain[d] = real_delta
        null_deltas = [domains[d][rep]["deltalogl"] for rep in common_reps if domains[d][rep]["deltalogl"] is not None]
        if real_delta is not None and null_deltas:
            n_extreme = sum(1 for v in null_deltas if v >= real_delta)
            perm_p_by_domain[d] = round((n_extreme + 1) / (len(null_deltas) + 1), 4)
        else:
            perm_p_by_domain[d] = "NA"
            real_delta_by_domain[d] = "NA" if real_delta is None else real_delta

    summary_file = output_dir / "family_summary.tsv"
    with open(summary_file, "w") as f:
        cols = ["Family_ID", "N_Domains", "N_Replicates_Paired"]
        cols += [f"Pct_Significant_{d}" for d in domain_names]
        cols += ["N_OR_Violations", "Pct_OR_Violations", "N_AND_Violations", "Pct_AND_Violations"]
        cols += [f"Real_DeltaLogL_{d}" for d in domain_names]
        cols += [f"Permutation_P_Upper_{d}" for d in domain_names]
        f.write("\t".join(cols) + "\n")
        row = [args.family, str(len(domain_names)), str(n_reps)]
        row += [str(per_domain_pct[d]) for d in domain_names]
        row += [str(n_or), str(pct(n_or)), str(n_and), str(pct(n_and))]
        row += [str(real_delta_by_domain[d]) for d in domain_names]
        row += [str(perm_p_by_domain[d]) for d in domain_names]
        f.write("\t".join(row) + "\n")

    logging.info(f"Family: {args.family}")
    logging.info(f"Domains: {domain_names}")
    logging.info(f"Paired replicates: {n_reps}")
    for d in domain_names:
        logging.info(f"  {d} significant-violation rate: {per_domain_pct[d]}%")
    logging.info(f"OR rate (either domain): {pct(n_or)}%")
    logging.info(f"AND rate (all domains): {pct(n_and)}%")
    logging.info(f"Written: {detail_file}, {summary_file}")


if __name__ == "__main__":
    main()
