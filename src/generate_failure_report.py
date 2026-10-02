#!/usr/bin/env python3
"""
Build a running summary of AU-test failure rates across every family
processed so far in a bulk run, and optionally email it.

For each family, scans however many domains currently have a completed
Step 6 domain summary (data/sim_outputs_bulk/<family>/Domain<N>/step6/results.tsv)
-- this is live/partial-safe: a family with only Domain1 done shows Domain1's
numbers and NA for Domain2, rather than waiting for the whole family to
finish. One row per family, one block of columns per domain (a column-name
union across families, since some have more domains than others -- see
src/combine_family_results.py for the same pattern).

"Failed the AU test" is reported two ways per domain, since it's ambiguous
which one a reader wants:
  - AU_Reject_* : raw AU rejections (p < 0.05), regardless of cause
  - Significant_* : AU rejections that also reflect genuine loss of
    monophyly (rejected AND mono=false) -- the definition used everywhere
    else in this pipeline (see docs/METHODS.md), matching how Peter's real
    analysis presumably calls a rejection

Also includes, per domain: every replicate's own delta-logL (not just the
mean/min/max), for downstream comparison against the real delta-logL once
available; and a real-vs-simulated %ID sanity check (skipped by default --
pass --with-pctid to include it; it reads every alignment file and is by
far the slowest part of this script once the dataset reaches thousands of
families).

Usage:
  python3 src/generate_failure_report.py --output-base data/sim_outputs_bulk \
      --manifest data/manifest/manifest_full.json \
      --out-tsv data/reports/failure_report.tsv \
      [--email you@example.com] [--with-pctid]
"""
import argparse
import csv
import glob
import json
import re
import smtplib
import sys
from collections import Counter
from datetime import datetime
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path

DOMAIN_RESULTS_RE = re.compile(r'^(F\d+)/Domain(\d+)/step6/results\.tsv$')


def read_fasta_seqs(path):
    """Ordered list of sequences (ignoring headers) from a FASTA alignment."""
    seqs = []
    cur = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if cur:
                    seqs.append(''.join(cur))
                cur = []
            else:
                cur.append(line)
        if cur:
            seqs.append(''.join(cur))
    return seqs


def pct_identity(seqs):
    """Pooled pairwise identity across an alignment: for every column, count
    matching vs comparable non-gap residue pairs, sum across columns, divide.
    Equivalent to averaging over all taxon pairs but O(taxa x sites) instead
    of O(taxa^2 x sites) -- matters once families have 100+ taxa. Returns
    None if there's nothing comparable (e.g. <2 sequences)."""
    if len(seqs) < 2:
        return None
    total_match = 0
    total_compare = 0
    for col in zip(*seqs):
        counts = Counter(c for c in col if c != '-')
        n = sum(counts.values())
        total_compare += n * (n - 1) // 2
        for cnt in counts.values():
            total_match += cnt * (cnt - 1) // 2
    if total_compare == 0:
        return None
    return 100.0 * total_match / total_compare


def compute_real_pct_id(trim_path):
    if not trim_path or not Path(trim_path).is_file():
        return None
    return pct_identity(read_fasta_seqs(trim_path))


def compute_mean_sim_pct_id(step3_dir):
    """Mean %ID across every rep_*_masked.aln found for this domain (however
    many replicates have made it through Step 3 so far)."""
    rep_files = sorted(glob.glob(str(Path(step3_dir) / "rep_*_masked.aln")))
    if not rep_files:
        return None
    vals = []
    for rf in rep_files:
        pid = pct_identity(read_fasta_seqs(rf))
        if pid is not None:
            vals.append(pid)
    if not vals:
        return None
    return sum(vals) / len(vals)


def read_replicate_deltalogl(step5_dir):
    """Per-replicate delta-logL values (unconstrained - constrained) from
    step5_results.txt, ordered by replicate number (rep_001, rep_002, ...).
    Returns a list, possibly shorter than N_REPLICATES if some replicates
    were skipped -- missing ones are left as None, not silently dropped,
    so position in the list always corresponds to that replicate number."""
    path = Path(step5_dir) / "step5_results.txt"
    if not path.is_file():
        return []
    values = {}
    max_rep = 0
    for line in open(path):
        m = re.match(r'^RESULT rep_(\d+) .*?\bdeltalogl=([^\s]+)', line)
        if not m:
            continue
        rep_num = int(m.group(1))
        max_rep = max(max_rep, rep_num)
        try:
            values[rep_num] = float(m.group(2))
        except ValueError:
            values[rep_num] = None
    return [values.get(i) for i in range(1, max_rep + 1)]


def family_sort_key(fam):
    m = re.match(r'^F(\d+)$', fam)
    return int(m.group(1)) if m else fam


def load_domain_results(output_base):
    """Returns {family: {domain_num: {column: value}}}."""
    families = {}
    for path in output_base.glob("*/Domain*/step6/results.tsv"):
        rel = path.relative_to(output_base)
        m = DOMAIN_RESULTS_RE.match(str(rel))
        if not m:
            continue
        fam, dom = m.group(1), m.group(2)
        with open(path, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            row = next(reader, None)
        if row is None:
            continue
        families.setdefault(fam, {})[dom] = row
    return families


def load_manifest_info(manifest_path):
    """Returns ({family: n_domains}, {(family, domain): trim_path})."""
    if not manifest_path or not Path(manifest_path).is_file():
        return {}, {}
    manifest = json.load(open(manifest_path))
    counts = {fam: len(info["domains"]) for fam, info in manifest.items()}
    trim_paths = {
        (fam, dom): dom_info["trim"]
        for fam, info in manifest.items()
        for dom, dom_info in info["domains"].items()
    }
    return counts, trim_paths


def build_table(families, expected_counts, trim_paths, output_base, skip_pctid=True):
    all_domain_nums = sorted({d for doms in families.values() for d in doms.keys()}, key=int)

    # Max replicate count, from the N_Replicates column already loaded in
    # memory (load_domain_results) -- no extra file reads needed. This is
    # a safe upper bound on the actual replicate NUMBER ever seen in a
    # step5_results.txt (N_REPLICATES controls how many get simulated in
    # the first place), so columns never come up short even though
    # N_Replicates itself counts usable (non-skipped) replicates, not the
    # max rep number.
    max_reps = 10
    for doms in families.values():
        for r in doms.values():
            try:
                max_reps = max(max_reps, int(r["N_Replicates"]))
            except (KeyError, ValueError):
                pass

    header = ["Family_ID", "N_Domains_Expected", "N_Domains_Completed"]
    for d in all_domain_nums:
        header += [
            f"Domain{d}_N_Taxa", f"Domain{d}_Alignment_Length", f"Domain{d}_N_Replicates",
            f"Domain{d}_N_AU_Reject_Total", f"Domain{d}_Pct_AU_Reject_Total",
            f"Domain{d}_N_Significant_Violations", f"Domain{d}_Pct_Significant_Violations",
            f"Domain{d}_Mean_AU_Pvalue",
            f"Domain{d}_Mean_Null_DeltaLogL", f"Domain{d}_Real_DeltaLogL",
            f"Domain{d}_Real_PctID", f"Domain{d}_Mean_Sim_PctID",
        ]
        # Every replicate's own ΔlogL (unconstrained - constrained), not just
        # the mean -- lets downstream analysis compare the real ΔlogL against
        # this null DISTRIBUTION directly (e.g. a permutation p-value or any
        # other summary stat), rather than only the mean/min/max already in
        # Domain{d}_Mean_Null_DeltaLogL. NA for any replicate that was
        # skipped or hasn't run yet.
        header += [f"Domain{d}_DeltaLogL_Rep{i}" for i in range(1, max_reps + 1)]

    rows = []
    for fam in sorted(families.keys(), key=family_sort_key):
        doms = families[fam]
        n_expected = expected_counts.get(fam, "NA")
        row = {
            "Family_ID": fam,
            "N_Domains_Expected": n_expected,
            "N_Domains_Completed": len(doms),
        }
        for d in all_domain_nums:
            r = doms.get(d)
            if r is None:
                for col in ("N_Taxa", "Alignment_Length", "N_Replicates",
                            "N_AU_Reject_Total", "Pct_AU_Reject_Total",
                            "N_Significant_Violations", "Pct_Significant_Violations",
                            "Mean_AU_Pvalue", "Mean_Null_DeltaLogL", "Real_DeltaLogL",
                            "Real_PctID", "Mean_Sim_PctID"):
                    row[f"Domain{d}_{col}"] = "NA"
                for i in range(1, max_reps + 1):
                    row[f"Domain{d}_DeltaLogL_Rep{i}"] = "NA"
            else:
                row[f"Domain{d}_N_Taxa"] = r["N_Taxa"]
                row[f"Domain{d}_Alignment_Length"] = r["Alignment_Length"]
                row[f"Domain{d}_N_Replicates"] = r["N_Replicates"]
                row[f"Domain{d}_N_AU_Reject_Total"] = r["N_AU_Rejects_Total"]
                row[f"Domain{d}_Pct_AU_Reject_Total"] = r["Pct_AU_Rejects_Total"]
                row[f"Domain{d}_N_Significant_Violations"] = r["N_Significant_Violations"]
                row[f"Domain{d}_Pct_Significant_Violations"] = r["Pct_Significant_Violations"]
                row[f"Domain{d}_Mean_AU_Pvalue"] = r["Mean_AU_Pvalue"]
                # Mean of ALL replicates' delta-logL (unconstrained - constrained),
                # recorded unconditionally in Step 5 regardless of AU reject/accept
                # -- see docs/METHODS.md. Real_DeltaLogL is NA until a real-data
                # comparison (real_data_delta_logl.sh, or Peter's real trees) has
                # been run for this domain -- present now so it populates
                # automatically once that's wired in, no report-format change needed.
                row[f"Domain{d}_Mean_Null_DeltaLogL"] = r["Null_DeltaLogL_Mean"]
                row[f"Domain{d}_Real_DeltaLogL"] = r["Real_DeltaLogL"]

                # Pairwise %ID: real alignment vs. mean across whatever
                # simulated replicates have reached Step 3 so far -- a
                # sanity check that the simulation is actually mimicking the
                # real data's divergence, not just running without error.
                # Reads full alignment files (real + up to N_REPLICATES
                # simulated) per domain -- by far the slowest part of this
                # script once the dataset reaches thousands of families;
                # skipped by default, pass --with-pctid to include it.
                if skip_pctid:
                    row[f"Domain{d}_Real_PctID"] = "NA"
                    row[f"Domain{d}_Mean_Sim_PctID"] = "NA"
                else:
                    trim_path = trim_paths.get((fam, d))
                    real_pid = compute_real_pct_id(trim_path)
                    step3_dir = output_base / fam / f"Domain{d}" / "step3"
                    sim_pid = compute_mean_sim_pct_id(step3_dir)
                    row[f"Domain{d}_Real_PctID"] = round(real_pid, 2) if real_pid is not None else "NA"
                    row[f"Domain{d}_Mean_Sim_PctID"] = round(sim_pid, 2) if sim_pid is not None else "NA"

                step5_dir = output_base / fam / f"Domain{d}" / "step5"
                rep_vals = read_replicate_deltalogl(step5_dir)
                for i in range(1, max_reps + 1):
                    v = rep_vals[i - 1] if i <= len(rep_vals) else None
                    row[f"Domain{d}_DeltaLogL_Rep{i}"] = v if v is not None else "NA"
        rows.append(row)

    return header, rows, all_domain_nums


def pooled_stats(families, all_domain_nums):
    """Pooled (population-level) significant-violation rate per domain slot,
    across every replicate of every family -- the actual number that matters
    for comparing against the real 41% figure, not a mean-of-percentages."""
    stats = {}
    for d in all_domain_nums:
        n_sig, n_reps, n_fams = 0, 0, 0
        for doms in families.values():
            r = doms.get(d)
            if r is None:
                continue
            n_sig += int(r["N_Significant_Violations"])
            n_reps += int(r["N_Replicates"])
            n_fams += 1
        stats[d] = {
            "n_families": n_fams,
            "n_replicates": n_reps,
            "n_significant": n_sig,
            "pct_significant": round(100.0 * n_sig / n_reps, 2) if n_reps else None,
        }
    return stats


def send_email(to_addr, subject, body, attachment_path):
    msg = MIMEMultipart()
    msg["Subject"] = subject
    msg["From"] = f"au-test-sim@{Path('/etc/hostname').read_text().strip() if Path('/etc/hostname').is_file() else 'localhost'}"
    msg["To"] = to_addr
    msg.attach(MIMEText(body, "plain"))

    with open(attachment_path, "rb") as f:
        part = MIMEApplication(f.read(), Name=Path(attachment_path).name)
    part["Content-Disposition"] = f'attachment; filename="{Path(attachment_path).name}"'
    msg.attach(part)

    with smtplib.SMTP("localhost") as smtp:
        smtp.send_message(msg)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--output-base", required=True)
    ap.add_argument("--manifest", default=None)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--email", default=None, help="If given, email the report to this address")
    ap.add_argument("--with-pctid", action="store_true",
                     help="Include the real-vs-simulated %%ID check (reads every alignment file; "
                          "by far the slowest part once the dataset reaches thousands of families -- "
                          "skipped by default)")
    args = ap.parse_args()

    output_base = Path(args.output_base)
    families = load_domain_results(output_base)
    expected_counts, trim_paths = load_manifest_info(args.manifest)

    if not families:
        print(f"No domain results found yet under {output_base}", file=sys.stderr)
        sys.exit(1)

    header, rows, all_domain_nums = build_table(families, expected_counts, trim_paths, output_base,
                                                 skip_pctid=not args.with_pctid)
    stats = pooled_stats(families, all_domain_nums)

    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t", restval="NA")
        writer.writeheader()
        writer.writerows(rows)

    n_complete_families = sum(
        1 for fam, doms in families.items()
        if expected_counts.get(fam) is not None and len(doms) >= expected_counts.get(fam, 0)
    )

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    summary_lines = [
        f"AU test simulation -- running failure-rate report, generated {timestamp}",
        "",
        f"Families with at least one domain completed: {len(families)}",
        f"Families with ALL expected domains completed: {n_complete_families}",
        "",
        "Pooled (population-level) significant-violation rate per domain slot",
        "(AU rejected AND composite clade genuinely lost monophyly -- the",
        "headline null false-positive rate, pooled across every replicate of",
        "every family with that domain completed so far, not a mean of",
        "per-family percentages):",
        "",
    ]
    for d in all_domain_nums:
        s = stats[d]
        pct = f"{s['pct_significant']}%" if s['pct_significant'] is not None else "NA"
        summary_lines.append(
            f"  Domain{d}: {s['n_significant']} / {s['n_replicates']} replicates "
            f"({pct}) across {s['n_families']} families"
        )

    summary_lines += [
        "",
        "Mean pairwise %ID, real data vs. simulated replicates (sanity check",
        "that the simulation is mimicking the real data's divergence -- big",
        "gaps here would mean the null distribution isn't trustworthy,",
        "regardless of what the AU-reject rate says):",
        "",
    ]
    for d in all_domain_nums:
        real_vals = [row[f"Domain{d}_Real_PctID"] for row in rows if row[f"Domain{d}_Real_PctID"] != "NA"]
        sim_vals = [row[f"Domain{d}_Mean_Sim_PctID"] for row in rows if row[f"Domain{d}_Mean_Sim_PctID"] != "NA"]
        if real_vals and sim_vals:
            mean_real = sum(real_vals) / len(real_vals)
            mean_sim = sum(sim_vals) / len(sim_vals)
            summary_lines.append(
                f"  Domain{d}: real={mean_real:.2f}%  simulated={mean_sim:.2f}%  "
                f"(diff={mean_real - mean_sim:+.2f} points, n={len(real_vals)} families)"
            )
        else:
            summary_lines.append(f"  Domain{d}: NA (no data yet)")

    summary_lines += [
        "",
        f"Full per-family, per-domain table attached: {out_path.name}",
        f"(also written to: {out_path})",
    ]
    body = "\n".join(summary_lines)

    print(body)

    if args.email:
        subject = f"AU test sim report: {len(families)} families ({n_complete_families} fully complete), {timestamp}"
        send_email(args.email, subject, body, out_path)
        print(f"\nEmailed to {args.email}", file=sys.stderr)


if __name__ == "__main__":
    main()
